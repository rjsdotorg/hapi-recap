#!/usr/bin/env python
"""
HAPI-RECAP: Reconstruct Parent Genomes from Sibling IBD and HAPI2 Phasing.

This module combines phased haplotypes from HAPI2 with IBD segments to
reconstruct parental DNA when parent genotypes are unavailable or incomplete.
It identifies which phased haplotypes belong to which parent using IBD overlap
patterns with relatives, and infers parent sexes using sex-specific genetic maps
and detected crossover events.

References:
    Qiao Y, et al. "Reconstructing parent genomes using siblings and other
    relatives." bioRxiv (2024). https://doi.org/10.1101/2024.05.10.593578
"""

import argparse
import itertools
import json
import logging
import sys
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import piso
from pysam import VariantFile, VariantHeader


REQUIRED_HAPI_KEYS = {"marker", "physpos", "chr", "chrstr"}


def _require_columns(df, required_columns, label):
    """Validate that a DataFrame contains required columns."""
    missing = [column for column in required_columns if column not in df.columns]
    if missing:
        raise KeyError(f"{label} is missing required columns: {', '.join(missing)}")


def _get_relative_ibd_groups(parent_ibd, chrom):
    """Return IBD segments on one chromosome grouped by relative id."""
    rel_ibd_chrom = parent_ibd[parent_ibd.chromosome == chrom]
    if len(rel_ibd_chrom.index) == 0:
        return {}
    return {
        rel_id: rel_df.reset_index(drop=True)
        for rel_id, rel_df in rel_ibd_chrom.groupby("id2", sort=False)
    }


def _annotate_overlap_coordinates(overlap_df, p_seg_start, p_seg_end, chrom, genetic_map):
    """Attach clipped overlap coordinates and genetic positions."""
    overlap_rel_ibd = overlap_df.copy()
    overlap_rel_ibd["relative_start"] = overlap_rel_ibd["start"].clip(lower=p_seg_start)
    overlap_rel_ibd["relative_end"] = overlap_rel_ibd["end"].clip(upper=p_seg_end)
    overlap_rel_ibd["rel_start_cm"] = [genetic_map[chrom][index] for index in overlap_rel_ibd["relative_start"]]
    overlap_rel_ibd["rel_end_cm"] = [genetic_map[chrom][index] for index in overlap_rel_ibd["relative_end"]]
    return overlap_rel_ibd


def parse_args():
    """
    Parse command-line arguments for HAPI-RECAP.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with the following attributes:
        - inf_hapi_json (str): Path to JSON file output by HAPI2
        - ibd_feather (str): Path pattern to IBD Feather files
        - bim (str): Path to PLINK BIM file
        - sex_avg_map (str): Path to sex-averaged genetic map
        - male_map (str): Path to male-specific genetic map
        - female_map (str): Path to female-specific genetic map
        - co_dir (str): Directory containing HAPI2 crossover files
        - out (str): Output directory for VCF files
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-inf_hapi_json",
        help="parent json file output by HAPI2 run on the inferred data",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-ibd_feather",
        help="feather file containing the phased IBD calls",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-bim",
        help="bim file used for running HAPI2",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-sex_avg_map",
        help="sex averaged genetic map file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-co_dir",
        help="directory containing HAPI2-formatted crossover event files",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-female_map",
        help="female genetic map file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-male_map",
        help="male genetic map file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-out",
        help="file to print results to",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    return args


def merge_ibd(segments, genetic_map):
    """
    Merge overlapping IBD segments for each (id1, id2, chromosome) group.

    Uses interval union operations to combine overlapping IBD segments within
    each group, updating both SNP-based and genetic map-based coordinates.
    Adapted from 23andMe's RACMACHINE GitHub repository.

    Parameters
    ----------
    segments : pandas.DataFrame
        IBD segments with columns: id1, id2, chromosome, start, end, start_cm, end_cm
        where start/end are SNP indices and start_cm/end_cm are genetic positions.
    genetic_map : dict
        Dictionary mapping chromosome IDs to lists of genetic positions (cM).

    Returns
    -------
    pandas.DataFrame
        Merged IBD segments with same structure as input, grouped by
        (id1, id2, chromosome) and sorted accordingly.

    Notes
    -----
    Handles corner cases where SNP-based and genetic map-based coordinate
    counts may differ due to rounding.
    """
    # code adapted from https://github.com/23andme-private/racmachine/blob/348a856f4baa638cefdbdee9bacc08f506644aa6/racmachine/phasedibd/ibd.py#L227-L271
    def _merge_groupby_segments(groupby_df, genetic_map, chrom):
        interval_arr_snps = pd.arrays.IntervalArray.from_arrays(
            groupby_df["start"], groupby_df["end"], closed="left"
        )
        interval_arr_cms = pd.arrays.IntervalArray.from_arrays(
            groupby_df["start_cm"], groupby_df["end_cm"], closed="left"
        )
        merged_interval_arr_snps = piso.union(interval_arr_snps)
        merged_interval_arr_cms = piso.union(interval_arr_cms)
        start = [intv.left for intv in merged_interval_arr_snps]
        end = [intv.right for intv in merged_interval_arr_snps]
        start_cm = [intv.left for intv in merged_interval_arr_cms]
        end_cm = [intv.right for intv in merged_interval_arr_cms]
        if len(start) != len(start_cm):
            # corner case: because of rounding, can get fewer start/end cMs than SNPs
            start_cm = [genetic_map[chrom][s] for s in start]
            end_cm = [genetic_map[chrom][e] for e in end]
        df = pd.DataFrame({"start": start, "end": end, "start_cm": start_cm, "end_cm": end_cm})
        return df

    groupby_vars = ["id1", "id2", "chromosome"]
    merged_dfs = []
    for groupby_values, groupby_df in segments.groupby(groupby_vars, sort=False):
        df = _merge_groupby_segments(groupby_df, genetic_map, groupby_values[2])  # 2 is index of "chromosome" in groupby element
        for var, value in zip(groupby_vars, groupby_values):
            df[var] = value
        merged_dfs.append(df)
    merged_ibd = pd.concat(merged_dfs, ignore_index=True)
    new_cols = groupby_vars + [col for col in merged_ibd.columns if col not in groupby_vars]
    merged_ibd = merged_ibd[new_cols]
    return merged_ibd


def filter_overlap_ibd(segments, genetic_map):
    """
    Filter IBD segments to remove overlaps between different parents for same relative.

    For relatives (id2) with IBD to multiple parent IDs (id1), removes overlapping
    regions so that segments are non-overlapping. Keeps non-overlapping segments intact.
    Adapted from 23andMe's RACMACHINE GitHub repository.

    Parameters
    ----------
    segments : pandas.DataFrame
        IBD segments with columns: id1, id2, chromosome, start, end, start_cm, end_cm.
    genetic_map : dict
        Dictionary mapping chromosome IDs to lists of genetic positions (cM).

    Returns
    -------
    pandas.DataFrame
        Filtered IBD segments with overlaps removed, grouped by (id2, chromosome).

    Notes
    -----
    Uses interval difference operations to remove overlapping regions.
    Handles corner cases with rounding differences between coordinates.
    """
    # code adapted from https://github.com/23andme-private/racmachine/blob/348a856f4baa638cefdbdee9bacc08f506644aa6/racmachine/phasedibd/ibd.py#L227-L271
    def _filter_groupby_segments(groupby_df, genetic_map, chrom):
        interval_arrs_snps = []
        interval_arrs_cms = []
        parent_ids = groupby_df.id1.unique().tolist()
        for parent_id in parent_ids:
            interval_arrs_snps.append(
                pd.arrays.IntervalArray.from_arrays(
                    groupby_df[groupby_df.id1 == parent_id]["start"],
                    groupby_df[groupby_df.id1 == parent_id]["end"],
                    closed="left"  # actually it should be "both", but that causes a bug
                )
            )
            interval_arrs_cms.append(
                pd.arrays.IntervalArray.from_arrays(
                    groupby_df[groupby_df.id1 == parent_id]["start_cm"],
                    groupby_df[groupby_df.id1 == parent_id]["end_cm"],
                    closed="left"
                )
            )
        assert len(interval_arrs_snps) == 2
        dfs = []
        intersect_interval_arr_snps = piso.intersection(interval_arrs_snps[0], interval_arrs_snps[1])
        intersect_interval_arr_cms = piso.intersection(interval_arrs_cms[0], interval_arrs_cms[1])
        for parent_id, interval_arr_snps, interval_arr_cms in zip(parent_ids, interval_arrs_snps, interval_arrs_cms):
            subtracted_interval_arr_snps = piso.difference(interval_arr_snps, intersect_interval_arr_snps)
            subtracted_interval_arr_cms = piso.difference(interval_arr_cms, intersect_interval_arr_cms)
            start = [intv.left for intv in subtracted_interval_arr_snps]
            end = [intv.right for intv in subtracted_interval_arr_snps]
            start_cm = [intv.left for intv in subtracted_interval_arr_cms]
            end_cm = [intv.right for intv in subtracted_interval_arr_cms]
            if len(start) != len(start_cm):
                # corner case: because of rounding, can get fewer start/end cMs than SNPs
                start_cm = [genetic_map[chrom][s] for s in start]
                end_cm = [genetic_map[chrom][e] for e in end]
            df = pd.DataFrame({"start": start, "end": end, "start_cm": start_cm, "end_cm": end_cm})
            df["id1"] = parent_id
            dfs.append(df)
        return pd.concat(dfs)

    groupby_vars = ["id2", "chromosome"]
    filtered_dfs = []
    for groupby_values, groupby_df in segments.groupby(groupby_vars, sort=False):
        if len(groupby_df.id1.unique()) > 1:
            df = _filter_groupby_segments(groupby_df, genetic_map, groupby_values[1])  # 1 is index of "chromosome" groupby element
            if len(df.index) == 0:
                continue  # completely filtered
            for var, value in zip(groupby_vars, groupby_values):
                df[var] = value
            filtered_dfs.append(df)
        else:
            filtered_dfs.append(groupby_df)
    filtered_ibd = pd.concat(filtered_dfs, ignore_index=True)
    #new_cols = groupby_vars + [col for col in merged_ibd.columns if col not in groupby_vars]
    #merged_ibd = merged_ibd[new_cols]
    return filtered_ibd


def extract_parent_segments(parents, inf_hapi_data):
    """
    Extract genomic segments where HAPI2 phased two distinct parents.

    Identifies regions between P (phased parent) and other code sites where
    HAPI2 determined distinct parental haplotypes. Segments must span at least
    100 markers to be included.

    Parameters
    ----------
    parents : str
        Parent identifiers in format 'id1-id2'.
    inf_hapi_data : dict
        HAPI2 output dictionary containing 'codes', 'chr', 'chrstr', and other metadata.

    Returns
    -------
    dict
        Dictionary keyed by chromosome, with values as lists of [start, end] positions
        for each parent segment in marker coordinates (relative to chromosome start).

    Notes
    -----
    P sites indicate phase boundaries. Segments occur between P sites or at
    chromosome boundaries. Only segments with >= 100 markers are retained.
    """
    the_parent_segments = defaultdict(list)

    prev_chrom = None
    most_recent_P_site = None
    for marker_idx, (inf_codes, cur_chrom) in enumerate(zip(inf_hapi_data[parents]["codes"], inf_hapi_data["chr"])):
        # Two ways parent segment can end:
        # (1) prev chromosome ended
        if cur_chrom != prev_chrom:
            if prev_chrom is not None:
                assert most_recent_P_site is not None
                # marker indexes are relative to the first site on a given chromosome:
                first_marker = most_recent_P_site - inf_hapi_data["chrstr"][prev_chrom]
                last_marker = inf_hapi_data["chrstr"][cur_chrom]-1 - inf_hapi_data["chrstr"][prev_chrom]
                if last_marker - first_marker >= 100:
                    the_parent_segments[prev_chrom].append(
                        [first_marker, last_marker]
                    )
            most_recent_P_site = None
            prev_chrom = cur_chrom
        # (2) encountering a P site
        if (
                inf_codes is not None and
                (inf_codes[0] == 'P' or inf_codes[0] == 'PC' or inf_codes[0] == 'PA')
            ):
            if most_recent_P_site is not None:
                # marker indexes are relative to the first site on a given chromosome:
                assert prev_chrom == cur_chrom
                first_marker = most_recent_P_site - inf_hapi_data["chrstr"][cur_chrom]
                last_marker = marker_idx-1 - inf_hapi_data["chrstr"][cur_chrom]
                if last_marker - first_marker >= 100:
                    the_parent_segments[prev_chrom].append(
                        [first_marker, last_marker]
                    )
            most_recent_P_site = marker_idx

    # finished analyzing all chromosomes: put last segment into the_parent_segments
    assert most_recent_P_site is not None
    # marker indexes are relative to the first site on a given chromosome:
    first_marker = most_recent_P_site - inf_hapi_data["chrstr"][prev_chrom]
    last_marker = marker_idx-1 - inf_hapi_data["chrstr"][prev_chrom]
    if last_marker - first_marker >= 100:
        the_parent_segments[prev_chrom].append(
            [first_marker, last_marker]
        )

    return the_parent_segments


def find_purple(close_rels, parent_ibd, parent_segments, genetic_map):
    """
    Identify "purple" relatives with IBD to both parents.

    Purple relatives show evidence of being related to both parents, which may
    indicate parental relationships or misphasing. Relatives are considered purple
    if more than 10% of their IBD overlaps with parent segments from both parents.

    Parameters
    ----------
    close_rels : pandas.DataFrame
        Relatives with total IBD > 600 cM, indexed by relative ID.
    parent_ibd : pandas.DataFrame
        IBD segments for the parent pair with relatives.
    parent_segments : dict
        Parent segments keyed by chromosome with [start, end] positions.
    genetic_map : dict
        Dictionary mapping chromosome IDs to genetic positions (cM).

    Returns
    -------
    set
        Set of relative IDs identified as purple (related to both parents).

    Notes
    -----
    Filters overlaps by requiring > 6 cM of genetic distance to be counted.
    Logs information for each relative about purple overlap counts.
    """
    purple_rels = set()
    for rel_id in close_rels.index:
        num_purple_overlap = 0
        num_overlap = 0
        for chrom, chr_parent_segments in parent_segments.items():
            rel_ibd_chrom = _get_relative_ibd_groups(parent_ibd, chrom).get(rel_id)
            if rel_ibd_chrom is None:
                continue
            for p_seg in chr_parent_segments:
                p_seg_start = p_seg[0]
                p_seg_end = p_seg[1]

                overlap_rel_ibd = rel_ibd_chrom[
                    (rel_ibd_chrom.start < p_seg_end) &
                    (rel_ibd_chrom.end > p_seg_start)
                ].reset_index(drop=True)
                if len(overlap_rel_ibd.index) == 0:
                    continue

                overlap_rel_ibd = _annotate_overlap_coordinates(overlap_rel_ibd, p_seg_start, p_seg_end, chrom, genetic_map)
                # TODO: could tune this -- what fraction of the parent segments does it overlap, e.g.?
                overlap_rel_ibd = overlap_rel_ibd[
                    overlap_rel_ibd.rel_end_cm - overlap_rel_ibd.rel_start_cm > 6
                ]
                if len(overlap_rel_ibd.index) > 0:
                    # have >= 1 overlapping segment
                    num_overlap += 1
                    if len(overlap_rel_ibd.id1.unique()) > 1:
                        # purple overlap
                        num_purple_overlap += 1

        logging.info(f"  {rel_id} has {num_purple_overlap} purple overlaps of {num_overlap} overlapping IBD segments")
        if num_overlap == 0:
            continue
        if num_purple_overlap / num_overlap > .10:
            purple_rels.add(rel_id)

    return purple_rels


def calc_co_probs(chrom, the_cos, p_seg_start, p_seg_end, ss_genetic_map):
    """
    Calculate crossover log-LOD scores for parent sex inference.

    Computes Poisson log-likelihood (base 10) for observing observed crossovers
    under male vs. female-specific genetic maps. Used to infer parent sexes.

    Parameters
    ----------
    chrom : str or int
        Chromosome identifier.
    the_cos : list of dict
        Two dictionaries containing crossover positions for each parent.
        Each dict maps chromosome to list of (marker_up, marker_down) tuples.
    p_seg_start : int
        Start position of parent segment in marker coordinates.
    p_seg_end : int
        End position of parent segment in marker coordinates.
    ss_genetic_map : list of dict
        Two genetic maps (male, female) mapping chromosome to genetic positions.

    Returns
    -------
    list of float
        Two-element list of LOD scores:
        - Index 0: Log-likelihood if parent 0 is male and parent 1 is female
        - Index 1: Log-likelihood if parent 0 is female and parent 1 is male

    Notes
    -----
    Only crossovers spanning <= 45 markers and falling within the segment are counted.
    Uses Poisson distribution with expected rate based on genetic map distance.
    """
    # probabilities that:
    # index 0: parent 0 is male and 1 female (default)
    # index 1: parent 0 is female and 1 male (swapped)
    probs = [ 0.0, 0.0 ]
    for parent_idx, parent_cos in enumerate(the_cos):
        for (marker_up, marker_down) in parent_cos[chrom]:
            if marker_down <= p_seg_start or marker_up > p_seg_end:
                continue  # the CO does not fall in this parent segment
            if marker_down - marker_up > 45:
                continue

            # the following inappropriately shortens the CO bounds: even though it
            # may go outside the segment, we shouldn't imply that the CO is
            # more tightly localized than what the raw data says
            #marker_up = max(p_seg_start, marker_up)
            #marker_down = min(p_seg_end, marker_down)

            for sex, sex_map in enumerate(ss_genetic_map):
                genet_up = sex_map[chrom][marker_up]
                genet_down = sex_map[chrom][marker_down]
                genet_length = genet_down - genet_up
                if genet_length == 0.0:
                    EPSILON = 1e-6
                    genet_length = EPSILON
                assert(genet_length > 0.0)
                # Poisson log (base 10) probability of 1 event in
                # genet_length region:
                # (using parent_idx ^ sex gives 0 if parent_idx == sex
                #  and 1 otherwise)
                probs[parent_idx ^ sex] += (np.log(genet_length) - genet_length) / np.log(10)

    return probs


def collect_update_segment_overlaps(the_parent_segments, rels_to_analyze, parent_ibd, genetic_map, father):
    """
    Link parent segments to relatives and count co-inheritance patterns.

    For each parent segment, determines which relatives share IBD with it and
    counts how many segments relative pairs share with the same parent vs.
    opposite parents. Also calculates total overlap length per relative.

    Parameters
    ----------
    the_parent_segments : dict
        Parent segments keyed by chromosome with [start, end] positions.
        Modified in-place to append relative_parent_linkage information.
    rels_to_analyze : set
        Set of relative IDs to analyze (excludes purple relatives).
    parent_ibd : pandas.DataFrame
        IBD segments for parent pair with relatives.
    genetic_map : dict
        Dictionary mapping chromosome IDs to genetic positions (cM).
    father : int
        Integer ID of the father.

    Returns
    -------
    tuple
        - rel_pairs_side_counts (dict): Maps frozenset of relative pairs to
          [count_same_parent, count_opposite_parent]
        - rel_overlap_length (dict): Maps relative ID to total overlap length in cM

    Notes
    -----
    Updates the_parent_segments in-place by appending relative_parent_linkage
    dictionaries (or None for purple segments).
    Segments marked purple (overlapping both parents) are skipped.
    """
    # indexed on a frozenset of two relative ids, each element is a 2-length list of counts
    # of how many parent segments the two relatives share IBD with the same parent vs opposite
    # parents
    rel_pairs_side_counts = dict()
    # per-sample total length of overlapping parent segments
    rel_overlap_length = defaultdict(lambda: 0)

    for chrom, chr_parent_segments in the_parent_segments.items():
        rel_ibd_by_id = _get_relative_ibd_groups(parent_ibd, chrom)
        # Below, we will append a dict to this_parent_segments containing the parent each relative has an IBD segment to
        for this_parent_seg in chr_parent_segments:
            [p_seg_start, p_seg_end] = this_parent_seg
            relative_parent_linkage = dict()

            parent_rels = [set(), set()]  # relatives connected to parent 0/1
            segment_is_purple = False

            for rel_id in rels_to_analyze:
                rel_ibd_chrom = rel_ibd_by_id.get(rel_id)
                if rel_ibd_chrom is None:
                    continue
                overlap_rel_ibd = rel_ibd_chrom[
                    (rel_ibd_chrom.start < p_seg_end) &
                    (rel_ibd_chrom.end > p_seg_start)
                ].reset_index(drop=True)
                if len(overlap_rel_ibd.index) == 0:
                    # No IBD segment overlap
                    continue

                overlap_rel_ibd = _annotate_overlap_coordinates(overlap_rel_ibd, p_seg_start, p_seg_end, chrom, genetic_map)
                # TODO: could tune this -- what fraction of the parent segments does it overlap, e.g.?
                overlap_rel_ibd = overlap_rel_ibd[
                    overlap_rel_ibd.rel_end_cm - overlap_rel_ibd.rel_start_cm > 6
                ]
                if len(overlap_rel_ibd.index) == 0:
                    continue
                if len(overlap_rel_ibd.id1.unique()) > 1:
                    # purple overlap: skip
                    segment_is_purple = True
                    break

                overlap_parent = overlap_rel_ibd.id1.iloc[0]
                parent_idx = 0 if overlap_parent == father else 1
                parent_rels[parent_idx].add(rel_id)
                relative_parent_linkage[rel_id] = parent_idx

            if segment_is_purple:
                this_parent_seg.append(None)
                continue
            this_parent_seg.append(relative_parent_linkage)

            # pairs with IBD to the same parent:
            for parent_idx in range(2):
                for pair in itertools.combinations(parent_rels[parent_idx], 2):
                    if frozenset(pair) not in rel_pairs_side_counts:
                        rel_pairs_side_counts[frozenset(pair)] = [0, 0]
                    rel_pairs_side_counts[frozenset(pair)][0] += 1
            # pairs with IBD to opposite parents:
            for pair in itertools.product(parent_rels[0], parent_rels[1]):
                if frozenset(pair) not in rel_pairs_side_counts:
                    rel_pairs_side_counts[frozenset(pair)] = [0, 0]
                rel_pairs_side_counts[frozenset(pair)][1] += 1

            # add to overlap length:
            length_cm = genetic_map[chrom][p_seg_end] - genetic_map[chrom][p_seg_start]
            for rel in parent_rels[0] | parent_rels[1]:
                rel_overlap_length[rel] += length_cm

    return rel_pairs_side_counts, rel_overlap_length


def determine_overlap_rel_parent_orient(rels_to_analyze, rel_overlap_length, rel_pairs_side_counts):
    """
    Determine parent orientation based on relative IBD co-inheritance patterns.

    Identifies which relatives are informative for parent assignment by selecting
    the relative with maximum overlap length and finding best secondary relative
    based on co-inheritance scores.

    Parameters
    ----------
    rels_to_analyze : set
        Set of relative IDs to analyze.
    rel_overlap_length : dict
        Maps relative ID to total overlap length in cM.
    rel_pairs_side_counts : dict
        Maps frozenset of relative pairs to [count_same, count_opposite].

    Returns
    -------
    dict
        Maps relative ID to orientation (0 or 1) indicating which parent they
        are primarily linked to.

    Notes
    -----
    Selects maximum overlap relative as anchor (always orientation 0).
    Secondary relative selected based on co-inheritance score:
    score = max(same, opposite) - 3*min(same, opposite)
    Skips pairs with min(same, opposite) > 2 (too ambiguous).
    """
    overlap_rel_parent_orient = dict()
    if len(rels_to_analyze) > 0 and len(rel_overlap_length) > 0:
        max_overlap_rel = max(rel_overlap_length, key=rel_overlap_length.get)
        overlap_rel_parent_orient[max_overlap_rel] = 0

        best_other_rel = None
        best_other_score = None
        best_other_val = None
        for other_rel in rel_overlap_length.keys():
            key = frozenset((max_overlap_rel, other_rel))
            if key not in rel_pairs_side_counts:
                continue

            val = rel_pairs_side_counts[key]
            if min(val) > 2:
                continue  # too ambiguous
            score = max(val) - 3 * min(val)
            if best_other_score is None or score > best_other_score:
                best_other_rel = other_rel
                best_other_score = score
                best_other_val = val
        if best_other_rel is not None:
            orient = 0 if best_other_val[0] > best_other_val[1] else 1
            overlap_rel_parent_orient[best_other_rel] = orient

    return overlap_rel_parent_orient


def link_segs_parents(the_parent_segments, overlap_rel_parent_orient, the_cos, ss_genetic_map):
    """
    Link parent segments to specific parents based on IBD and crossover evidence.

    For each parent segment, determines which parent (0 or 1) it belongs to by
    combining IBD linkage information with crossover probability scores. Segments
    without IBD linkage are marked with None parent_idx but still assigned CO scores.

    Parameters
    ----------
    the_parent_segments : dict
        Modified parent segments containing relative_parent_linkage info.
        Keyed by chromosome with values [start, end, relative_parent_linkage].
    overlap_rel_parent_orient : dict
        Maps relative ID to parent orientation (0 or 1).
    the_cos : list of dict
        Two dictionaries with crossover positions per parent per chromosome.
    ss_genetic_map : list of dict
        Two sexual-specific genetic maps (male index 0, female index 1).

    Returns
    -------
    dict
        Dictionary keyed by chromosome with values as lists of tuples:
        (p_seg_start, p_seg_end, parent_idx, co_probs)
        where parent_idx is 0/1 if IBD-linked or None if not.

    Notes
    -----
    parent_idx values: 0 indicates one parent, 1 the other, None if unlinked.
    All segments include crossover probability scores for later sex inference.
    """
    # dictionary keyed on chromosome that stores a list of tupes;
    # each tuple stores information about the parent segments as:
    #   (p_seg_start, p_seg_end, parent_idx, co_probs)
    # see below for how parent_idx is defined
    parent_segs_linkage = defaultdict(list)

    for chrom, chr_parent_segments in the_parent_segments.items():
        for this_parent_seg in chr_parent_segments:
            [p_seg_start, p_seg_end, relative_parent_linkage] = this_parent_seg

            # probabilities that:
            # index 0: parent 0 is male and 1 female (default)
            # index 1: parent 0 is female and 1 male (swapped)
            co_probs = calc_co_probs(chrom, the_cos, p_seg_start, p_seg_end, ss_genetic_map)

            if relative_parent_linkage is None:
                # No IBD segment overlap, but we may be able to use the CO probabilities to add this
                parent_segs_linkage[chrom].append(
                    (p_seg_start, p_seg_end, None, co_probs)
                )
                continue

            overlap_rel_this_seg = relative_parent_linkage.keys() & overlap_rel_parent_orient.keys()
            if len(overlap_rel_this_seg) == 0:
                # No IBD segment overlap, but we may be able to use the CO probabilities to add this
                parent_segs_linkage[chrom].append(
                    (p_seg_start, p_seg_end, None, co_probs)
                )
                continue

            parent_idx = None
            for rel_id in overlap_rel_this_seg:
                # relative_parent_linkage is 0 if rel_id is IBD to the (locally labeled) father at this segment and 1 otherwise
                # overlap_rel_parent_orient represents which relatives are linked to the same parents 0 or 1
                this_parent_idx = relative_parent_linkage[rel_id] ^ overlap_rel_parent_orient[rel_id]
                if parent_idx is None:
                    parent_idx = this_parent_idx
                elif parent_idx != this_parent_idx:
                    # suggests purple segment: skip
                    continue

            parent_segs_linkage[chrom].append(
                (p_seg_start, p_seg_end, parent_idx, co_probs)
            )

    return parent_segs_linkage


def print_segment(inf_hapi_data, parents, snp_to_alleles, this_parent_orient,
                  chrom, p_seg_start, p_seg_end, num_data_sites, vcf_out):
    """
    Write reconstructed parent genotypes for a chromosomal segment to VCF.

    Extracts inferred parent haplotypes from HAPI2 output for a segment,
    orients them according to determined parent sexes, and writes VCF records.
    Handles missing and half-missing sites appropriately.

    Parameters
    ----------
    inf_hapi_data : dict
        HAPI2 output JSON data with phased haplotypes and metadata.
    parents : str
        Parent identifier string (format 'id1-id2').
    snp_to_alleles : dict
        Maps SNP IDs to (ref, alt) allele tuples.
    this_parent_orient : int
        0 if parent 0 is father (parent 1 is mother), 1 if reversed.
    chrom : str or int
        Chromosome identifier.
    p_seg_start : int
        Segment start position (marker index relative to chromosome).
    p_seg_end : int
        Segment end position (marker index relative to chromosome).
    num_data_sites : list of list
        Modified in-place: counts of [full, half, missing] sites for each parent.
    vcf_out : pysam.VariantFile
        Open VCF file for writing records.

    Notes
    -----
    HAPI2 error codes (E, R, ?) result in missing parent genotypes.
    Half-missing sites (one allele is 0) are coded as homozygous for known allele.
    this_parent_orient=0: parent 0 is father, parent 1 is mother
    this_parent_orient=1: parent 0 is mother, parent 1 is father (reversed)
    """
    inf_parent_haps = inf_hapi_data[parents]["parhaps"]
    inf_codes = inf_hapi_data[parents]["codes"]
    chr_start_marker_idx = inf_hapi_data["chrstr"][chrom]
    for marker_idx in range(chr_start_marker_idx + p_seg_start, chr_start_marker_idx + p_seg_end + 1):
        assert inf_hapi_data["chr"][marker_idx] == chrom

        set_missing = False
        if inf_codes[marker_idx] is not None:
            if inf_codes[marker_idx][0] in {'E', 'R', '?'}:
                # Mendelian error, recombination-detected error, or ambiguous
                # site; parents should be missing
                set_missing = True
            else:
                # pretty sure the above codes always come first -- confirm:
                assert set(inf_codes[marker_idx]) & {'E', 'R', '?'} == set()

        rec = vcf_out.new_record()
        rec.chrom = chrom
        rec.pos = inf_hapi_data["physpos"][marker_idx]
        rec.id = inf_hapi_data["marker"][marker_idx]
        rec.alleles = snp_to_alleles[rec.id]
        # inf_parent_haps[0] are the alleles assigned to parent 0
        # inf_parent_haps[1] are the alleles assigned to parent 1
        # so inf_parent_haps[0^this_parent_orient] give the inferred paternal haplotypes
        # and index 1^this_parent_orient give the inferred maternal haplotypes
        for parent_idx, sample in enumerate(rec.samples.values()):

            parent_alleles = (inf_parent_haps[parent_idx^this_parent_orient][0][marker_idx],
                    inf_parent_haps[parent_idx^this_parent_orient][1][marker_idx])
            if set_missing or parent_alleles == ('0', '0'):  # fully missing
                num_data_sites[parent_idx][0] += 1  # site with 0 alleles reconstructed
                # assigning .alleles here gives an error, so we use allele_indices
                sample.allele_indices = (None, None)
                sample.phased = True
                continue
            if parent_alleles[0] == '0' or parent_alleles[1] == '0':  # half-missing
                num_data_sites[parent_idx][1] += 1  # site with 1 allele reconstructed
                # make homozygous for the allele we do have:
                if parent_alleles[0] != '0':
                    parent_alleles = (parent_alleles[0], parent_alleles[0])
                else:
                    parent_alleles = (parent_alleles[1], parent_alleles[1])
            else:
                num_data_sites[parent_idx][2] += 1  # site with 2 alleles reconstructed

            sample.alleles = parent_alleles
            sample.phased = True
        vcf_out.write(rec)


def reconstruct(parents, inf_hapi_data, snp_to_alleles, ibd, genetic_map, ss_genetic_map,
                chrom_names, args):
    """
    Main reconstruction pipeline for a parent pair.

    Orchestrates the full pipeline: filters IBD segments, identifies close relatives,
    detects purple relatives, determines parent orientations using IBD and crossovers,
    links segments to parents, infers parent sexes, and outputs VCF with reconstructed
    genotypes and reconstruction statistics.

    Parameters
    ----------
    parents : str
        Parent identifier string (format 'father_id-mother_id').
    inf_hapi_data : dict
        HAPI2 output JSON containing phased haplotypes and metadata.
    snp_to_alleles : dict
        Maps SNP IDs to (ref, alt) allele tuples.
    ibd : pandas.DataFrame
        IBD segments for all relative pairs.
    genetic_map : dict
        Sex-averaged genetic map (chromosome -> genetic positions).
    ss_genetic_map : list of dict
        Sex-specific genetic maps [male_map, female_map].
    chrom_names : list
        List of chromosome identifiers in order.
    args : argparse.Namespace
        Command-line arguments containing output directory and crossover directory.

    Returns
    -------
    None
        Outputs VCF file to args.out/[parents].vcf and logs statistics.

    Notes
    -----
    Parent segments are classified by IBD linkage evidence vs. crossover-only evidence.
    Reconstructed segments output with estimated fraction of parent genomes recovered.
    Statistics distinguish fully reconstructed (both alleles) vs. half-reconstructed sites.
    """
    logging.info(f"Analayzing parents {parents}")

    # TODO: document assumption re: '-' in parent ids
    father, mother = [int(p) for p in parents.split('-')]

    # because the parent haplotype assignments can be swapped, IBD to
    # either parent is what we're interested in
    parent_ibd = ibd[(ibd.id1 == father) | (ibd.id1 == mother)]
    # remove segments between the parents
    parent_ibd = parent_ibd[(parent_ibd.id2 != father) & (parent_ibd.id2 != mother)]
    # TODO: combine merge_ibd() with filter_overlap_ibd()
    parent_ibd = merge_ibd(parent_ibd, genetic_map)
    # filter segments that overlap both parents
    parent_ibd = filter_overlap_ibd(parent_ibd, genetic_map)
    parent_ibd = parent_ibd.astype({"start": int, "end": int})
    parent_ibd["length_cm"] = parent_ibd["end_cm"] - parent_ibd["start_cm"]
    parent_ibd = parent_ibd[parent_ibd.length_cm > 9]
    close_rels = (
        parent_ibd
        .groupby(["id2"], sort=False)[["length_cm"]]
        .sum()
        .sort_values(by=["length_cm"], ascending=False)
    )
    close_rels = close_rels[close_rels.length_cm > 600]

    # returns a dictionary keyed on chromosome that stores a list of lists;
    # each inner list stores [first_marker, last_marker] for the corresponding
    # segment
    the_parent_segments = extract_parent_segments(parents, inf_hapi_data)

    # find individuals that may be related to both parents. We term these
    # purple relatives: if we code one parent's relatives in red and the other
    # in blue, the analogy is that purple relatives are related to both parents
    purple_rels = set()
    if len(close_rels.index) > 0:
        # find purple relatives
        purple_rels = find_purple(close_rels, parent_ibd, the_parent_segments, genetic_map)

    # decide which relatives are connected to the same/different parents
    rels_to_analyze = set(close_rels.index) - purple_rels
    rel_pairs_side_counts, rel_overlap_length = collect_update_segment_overlaps(the_parent_segments, rels_to_analyze, parent_ibd, genetic_map, father)

    overlap_rel_parent_orient = determine_overlap_rel_parent_orient(rels_to_analyze, rel_overlap_length, rel_pairs_side_counts)

    # read in the crossover locations, one dictionary per parent
    the_cos = [ defaultdict(list), defaultdict(list) ]
    for chrom in range(1, 23):
        with open(f"{args.co_dir}/co-{parents}.{chrom}", "r") as fin:
            for line in fin:
                if line[0] == "#":
                    continue

                fields = line.strip().split()
                parent_idx = int(fields[2])
                marker_up = int(fields[6])
                marker_down = int(fields[8])
                the_cos[parent_idx][str(chrom)].append((marker_up, marker_down))

    # link the parent segments to the parents and get CO probs
    parent_segs_linkage = link_segs_parents(the_parent_segments, overlap_rel_parent_orient, the_cos, ss_genetic_map)

    # use CO probs to determine which parent corresponds to the parent_orient = 0 value
    overall_co_probs = 0.0  # summed over all IBD-linked segments:
    for chrom, chr_parent_seg_link in parent_segs_linkage.items():
        for (p_seg_start, p_seg_end, parent_orient, co_probs) in chr_parent_seg_link:
            # Note: the paper uses a formula that's -1* this
            co_prob_diff = co_probs[0] - co_probs[1]

            if parent_orient is not None:  # have IBD-linkage to this segment
                overall_co_probs += co_prob_diff if parent_orient == 0 else -co_prob_diff

    # if overall_co_probs > 0, crossovers suggest that parent_orient 0 segments have the dad assigned as parent 0
    # and mom as parent 1. The reverse is true if overall_co_probs < 0
    overall_parent_orient = 0 if overall_co_probs > 0 else 1
    # TODO: treat sexes as ambiguous if abs(overall_co_probs) < 3?

    # Assign segments to parents and quantify reconstruction
    num_markers = len(inf_hapi_data["chr"])
    # numbers for IBD-based reconstruction
    # top-level list corresponds to parents 0 and 1 (father and mother) and
    # second-level list counts the number of fully missing, half-missing, and
    #   full data sites _within_ reconstructed segments
    num_ibd_data_sites = [[0, 0, 0], [0, 0, 0]]
    num_segments = 0
    # code doesn't currently do SS map-only reconstruction (i.e., we use IBD if
    # we have it)
    # numbers for additional placed segments using crossovers (COs) and SS maps:
    num_co_data_sites = [[0, 0, 0], [0, 0, 0]]
    num_co_segments = 0

    # setup header
    header = VariantHeader()
    for chrom in chrom_names:
        header.contigs.add(chrom)
    header.formats.add("GT", 1, "String", "Genotype")
    for parent_id in (father, mother):
        header.add_sample(str(parent_id))
    # print VCF:
    with VariantFile(f"{args.out}/{parents}.vcf", mode="w", header=header) as vcf_out:
        for chrom, chr_parent_seg_link in parent_segs_linkage.items():
            for (p_seg_start, p_seg_end, parent_orient, co_probs) in chr_parent_seg_link:
                co_prob_diff = co_probs[0] - co_probs[1]

                if parent_orient is not None:  # have IBD-linkage to this segment
                    num_segments += 1
                    this_parent_orient = overall_parent_orient ^ parent_orient
                    logging.info(f"  Linked segment {chrom} {p_seg_start}-{p_seg_end} with orientation {this_parent_orient}")
                    # flip sign for overall opposite parent orientation
                    if overall_parent_orient == 1:
                        co_prob_diff *= -1
                    if co_prob_diff <= -3:
                        # Crossovers in this segment suggest inconsistent sex relative to the other
                        # IBD-linked segments
                        logging.info(f"  NOTE: segment {chrom} {p_seg_start}-{p_seg_end} with orientation {this_parent_orient} has LOD {co_prob_diff:.1f}")

                    #if abs(co_prob_diff) >= 3:
                    #    num_map_segments += 1  # would place this segment using SS map only

                    print_segment(inf_hapi_data, parents, snp_to_alleles, this_parent_orient, chrom, p_seg_start, p_seg_end, num_ibd_data_sites, vcf_out)

                elif abs(co_prob_diff) >= 3:  # no IBD-linkage, but |LOD| >= 3
                    # TODO: add option to only incorporate these segments if they cover a specified fraction of the genome
                    num_co_segments += 1
                    logging.info(f"  Assigning segment {chrom} {p_seg_start}-{p_seg_end} to parents only using crossover LOD {abs(co_prob_diff):.1f}")

                    this_parent_orient = 0 if co_prob_diff > 0 else 1
                    print_segment(inf_hapi_data, parents, snp_to_alleles, this_parent_orient, chrom, p_seg_start, p_seg_end, num_co_data_sites, vcf_out)

    # quantify reconstruction
    num_ibd_reconstructed = [0, 0]
    num_co_reconstructed = [0, 0]
    frac = [0, 0]
    for parent_idx in range(2):
        num_ibd_reconstructed[parent_idx] = num_ibd_data_sites[parent_idx][1] +\
                                        2 * num_ibd_data_sites[parent_idx][2]
        num_co_reconstructed[parent_idx] = num_co_data_sites[parent_idx][1] +\
                                        2 * num_co_data_sites[parent_idx][2]
        frac[parent_idx] = (num_ibd_reconstructed[parent_idx] + num_co_reconstructed[parent_idx]) / (2 * num_markers)
    logging.info(f"Reconstructed {frac[0]*100:.1f}% of father's and {frac[1]*100:.1f}% of mother's variants")
    logging.info(f"from {num_segments} segments using IBD and {num_co_segments} segments using crossovers (CO):")
    logging.info("\tFather Full\tFather Half\tMother Full\tMother Half")
    ibd_strings = []
    for parent_idx in range(2):
        ibd_strings.append(f"\t{num_ibd_data_sites[parent_idx][2]/num_markers*100:.1f}%\t\t{num_ibd_data_sites[parent_idx][1]/num_markers*100:.1f}%\t")
    logging.info(f"IBD{ibd_strings[0]}{ibd_strings[1]}")
    co_strings = []
    for parent_idx in range(2):
        co_strings.append(f"\t{num_co_data_sites[parent_idx][2]/num_markers*100:.1f}%\t\t{num_co_data_sites[parent_idx][1]/num_markers*100:.1f}%\t")
    logging.info(f"CO{co_strings[0]}{co_strings[1]}")


def main():
    """
    Main entry point for HAPI-RECAP.

    Loads all required input files (HAPI2 JSON, IBD segments, genetic maps),
    processes each parent pair, and outputs reconstructed VCF files.

    Returns
    -------
    None
        Exits with code 2 on error creating output directory, otherwise exits 0.

    Raises
    ------
    SystemExit
        If output directory cannot be created.
    """
    args = parse_args()

    # if necessary, create output directory
    try:
        Path(args.out).mkdir(parents=True, exist_ok=True)
    except Exception:
        print("Error creating output directory")
        print(traceback.format_exc(), end="")
        sys.exit(2)

    logging.info("Loading inferred json")
    with open(args.inf_hapi_json) as fin:
        inf_hapi_data = json.load(fin)
    missing_hapi_keys = REQUIRED_HAPI_KEYS - set(inf_hapi_data)
    if missing_hapi_keys:
        raise KeyError(f"HAPI2 JSON is missing required keys: {', '.join(sorted(missing_hapi_keys))}")

    logging.info("Loading alleles from bim file")
    snp_to_alleles = dict()
    with open(args.bim) as fin:
        for line in fin:
            fields = line.strip().split()
            # order of alleles here is alternate, reference, but for VCF, we
            # want the reverse:
            snp_to_alleles[fields[1]] = (fields[5], fields[4])

    # get genetic positions
    logging.info("Loading genetic map")
    genetic_map = defaultdict(list)
    with open(args.sex_avg_map, "r") as fin:
        for line in fin:
            fields = line.strip().split()
            genetic_map[fields[0]].append(float(fields[2]))

    # get genetic positions
    logging.info("Loading sex-specific genetic maps")
    ss_genetic_map = [ defaultdict(list), defaultdict(list) ]
    for sex_idx, sex_map in enumerate((args.male_map, args.female_map)):
        with open(sex_map, "r") as fin:
            for line in fin:
                fields = line.strip().split()
                ss_genetic_map[sex_idx][fields[0]].append(float(fields[2]))

    logging.info("Loading IBD segments")
    ibd_by_chrom = []
    for chrom in range(1, 23):
        filename = args.ibd_feather.replace("chr1", f"chr{chrom}")
        if not Path(filename).exists():
            raise FileNotFoundError(f"Missing IBD feather file: {filename}")
        chrom_ibd = pd.read_feather(filename)
        _require_columns(
            chrom_ibd,
            ["id1", "id2", "chromosome", "start", "end", "start_cm", "end_cm"],
            f"IBD feather file {filename}",
        )
        ibd_by_chrom.append(chrom_ibd)
    if not ibd_by_chrom:
        raise ValueError("No IBD feather files were loaded")
    ibd = pd.concat(ibd_by_chrom, ignore_index=True)

    chrom_names = [k for k, v in sorted(inf_hapi_data["chrstr"].items(), key=lambda item: item[1])]

    # analyze each family in turn
    for key in inf_hapi_data.keys():
        if key == "marker" or key == "physpos" or key == "chr" or key == "chrstr":
            continue

        parents = key
        reconstruct(parents, inf_hapi_data, snp_to_alleles, ibd, genetic_map, ss_genetic_map,
                    chrom_names, args)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
    main()
