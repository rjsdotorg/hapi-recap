#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze Reconstructed Parents Sex using Cousin IBD

Given two reconstructed parents and known paternal cousins, determines which
parent is paternal by comparing IBD segment totals to the cousins.

Logic:
  - Paternal parent should have more/longer/stronger IBD to paternal cousins
  - Autosomal IBD (chr1-22) is primary signal
  - X-chromosome IBD (chr23) is secondary confirmation
  - Outputs sexing report with confidence metrics
"""
import argparse
import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


def _merge_intervals(intervals):
    """Merge overlapping [start_cm, end_cm] intervals and return merged list."""
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [list(sorted_intervals[0])]
    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def _pick_column(columns, aliases):
    normalized = {str(col).strip().lower(): col for col in columns}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def _load_min_map(path):
    raw = pd.read_csv(path, sep='\t', header=0)
    chrom_col = _pick_column(raw.columns, ['chromosome', 'chrom', 'chr'])
    pos_col = _pick_column(raw.columns, ['position', 'pos'])
    cm_col = _pick_column(raw.columns, ['cm', 'length (cm)'])
    if chrom_col is None or pos_col is None or cm_col is None:
        raise KeyError('min_map is missing required columns: Chromosome, Position, cM')

    raw = raw[[chrom_col, pos_col, cm_col]].copy()
    raw.columns = ['chrom', 'position', 'cm']
    raw['chrom'] = pd.to_numeric(raw['chrom'], errors='coerce')
    raw['position'] = pd.to_numeric(raw['position'], errors='coerce')
    raw['cm'] = pd.to_numeric(raw['cm'], errors='coerce')
    raw = raw.dropna()
    raw['chrom'] = raw['chrom'].astype(int)
    raw = raw[raw['chrom'].between(1, 23)]

    by_pos = {}
    by_cm = {}
    for chrom, dfc in raw.groupby('chrom', sort=False):
        dfc = dfc.sort_values('position')
        by_pos[str(chrom)] = dfc['position'].to_numpy(dtype=np.float64)
        by_cm[str(chrom)] = dfc['cm'].to_numpy(dtype=np.float64)
    return by_pos, by_cm


def _load_bim_positions(bim_path):
    by_chrom = defaultdict(list)
    with open(bim_path, 'r', encoding='utf-8', errors='ignore') as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom = str(parts[0]).strip().upper().replace('CHR', '')
            if chrom == 'X':
                chrom = '23'
            if not chrom.isdigit():
                continue
            c = int(chrom)
            if not (1 <= c <= 23):
                continue
            try:
                bp = int(parts[3])
            except ValueError:
                continue
            by_chrom[str(c)].append(bp)
    return {chrom: np.asarray(vals, dtype=np.int64) for chrom, vals in by_chrom.items()}


def _segment_cm_from_idx(chrom, start_idx, end_idx, bim_positions, dmap_positions, dmap_cms):
    positions = bim_positions.get(chrom)
    map_pos = dmap_positions.get(chrom)
    map_cm = dmap_cms.get(chrom)
    if positions is None or map_pos is None or map_cm is None:
        return None
    if len(positions) == 0 or len(map_pos) == 0 or len(map_cm) == 0:
        return None
    s = int(min(max(start_idx, 0), len(positions) - 1))
    e = int(min(max(end_idx - 1, 0), len(positions) - 1))
    if e < s:
        s, e = e, s
    start_bp = float(positions[s])
    end_bp = float(positions[e])
    start_cm = float(np.interp(start_bp, map_pos, map_cm))
    end_cm = float(np.interp(end_bp, map_pos, map_cm))
    return max(0.0, end_cm - start_cm)


def load_feather_files(ibd_dir, chromosomes):
    """Load all IBD feather files for given chromosomes."""
    feathers = {}
    for chrom in chromosomes:
        feather_path = os.path.join(ibd_dir, f"chr{chrom}.feather")
        if not os.path.exists(feather_path):
            print(f"WARNING: {feather_path} not found, skipping chr{chrom}", flush=True)
            continue
        try:
            df = pd.read_feather(feather_path)
            feathers[chrom] = df
            print(f"Loaded chr{chrom}: {len(df)} segments", flush=True)
        except Exception as e:
            print(f"ERROR loading {feather_path}: {e}", flush=True)
    return feathers


def aggregate_ibd_by_parent(feathers, parent_ids, cousins):
    """
    Aggregate IBD segments by parent and cousin pairs.

    Returns dict: {parent_id: {cousin_id: {'count': int, 'total_cm': float, 'mean_cm': float}}}
    """
    aggregates = {pid: {} for pid in parent_ids}

    for chrom, df in feathers.items():
        if df.empty:
            continue

        for _, row in df.iterrows():
            id1 = str(row['id1'])
            id2 = str(row['id2'])
            segment_cm = float(row['end_cm']) - float(row['start_cm'])

            # Check if this pair involves a parent and a cousin
            parent_match = None
            cousin_match = None

            if id1 in parent_ids and id2 in cousins:
                parent_match = id1
                cousin_match = id2
            elif id2 in parent_ids and id1 in cousins:
                parent_match = id2
                cousin_match = id1

            if parent_match and cousin_match:
                if cousin_match not in aggregates[parent_match]:
                    aggregates[parent_match][cousin_match] = {
                        'count': 0,
                        'total_cm_raw': 0.0,
                        'segments': [],
                        'intervals_by_chrom': defaultdict(list),
                    }

                aggregates[parent_match][cousin_match]['count'] += 1
                aggregates[parent_match][cousin_match]['total_cm_raw'] += segment_cm
                aggregates[parent_match][cousin_match]['segments'].append(segment_cm)

                # Prefer merged non-overlapping cM spans to avoid double-counting.
                start_cm = float(row['start_cm'])
                end_cm = float(row['end_cm'])
                if np.isfinite(start_cm) and np.isfinite(end_cm):
                    lo, hi = sorted((start_cm, end_cm))
                    if hi > lo:
                        chrom_key = str(row['chromosome'])
                        aggregates[parent_match][cousin_match]['intervals_by_chrom'][chrom_key].append((lo, hi))

    # Calculate mean cM for each parent-cousin pair
    for parent_id in aggregates:
        for cousin_id in aggregates[parent_id]:
            seg_list = aggregates[parent_id][cousin_id]['segments']
            if seg_list:
                aggregates[parent_id][cousin_id]['mean_cm'] = np.mean(seg_list)
                aggregates[parent_id][cousin_id]['median_cm'] = np.median(seg_list)
                aggregates[parent_id][cousin_id]['max_cm'] = np.max(seg_list)
            else:
                aggregates[parent_id][cousin_id]['mean_cm'] = 0.0
                aggregates[parent_id][cousin_id]['median_cm'] = 0.0
                aggregates[parent_id][cousin_id]['max_cm'] = 0.0

            merged_total = 0.0
            merged_count = 0
            for _, intervals in aggregates[parent_id][cousin_id]['intervals_by_chrom'].items():
                merged_intervals = _merge_intervals(intervals)
                merged_count += len(merged_intervals)
                merged_total += sum((e - s) for s, e in merged_intervals)
            aggregates[parent_id][cousin_id]['merged_interval_count'] = merged_count
            aggregates[parent_id][cousin_id]['total_cm_merged'] = merged_total

    return aggregates


def aggregate_sibling_cousin_ibd(feathers, siblings, cousins, bim_positions, dmap_positions, dmap_cms):
    """Aggregate sibling-cousin evidence from pairwise IBD feathers.

    Uses BIM+min_map to recompute cM lengths from marker-index intervals so this
    mode remains valid even if feather start_cm/end_cm are all zero.
    """
    aggregates = {sid: {} for sid in siblings}

    for chrom, df in feathers.items():
        if df.empty:
            continue
        chrom_key = str(chrom)
        for _, row in df.iterrows():
            id1 = str(row['id1'])
            id2 = str(row['id2'])
            sib = id1 if id1 in siblings else (id2 if id2 in siblings else None)
            cous = id2 if id2 in cousins else (id1 if id1 in cousins else None)
            if sib is None or cous is None:
                continue

            cm_len = _segment_cm_from_idx(
                chrom_key,
                int(row['start']),
                int(row['end']),
                bim_positions,
                dmap_positions,
                dmap_cms,
            )
            if cm_len is None:
                continue

            if cous not in aggregates[sib]:
                aggregates[sib][cous] = {
                    'count': 0,
                    'total_cm_raw': 0.0,
                    'segments': [],
                    'intervals_by_chrom': defaultdict(list),
                }

            aggregates[sib][cous]['count'] += 1
            aggregates[sib][cous]['total_cm_raw'] += cm_len
            aggregates[sib][cous]['segments'].append(cm_len)

            # Intervals in cM space for overlap-merged totals
            positions = bim_positions.get(chrom_key)
            if positions is not None and len(positions) > 0:
                s = int(min(max(int(row['start']), 0), len(positions) - 1))
                e = int(min(max(int(row['end']) - 1, 0), len(positions) - 1))
                if e < s:
                    s, e = e, s
                start_bp = float(positions[s])
                end_bp = float(positions[e])
                start_cm = float(np.interp(start_bp, dmap_positions[chrom_key], dmap_cms[chrom_key]))
                end_cm = float(np.interp(end_bp, dmap_positions[chrom_key], dmap_cms[chrom_key]))
                lo, hi = sorted((start_cm, end_cm))
                if hi > lo:
                    aggregates[sib][cous]['intervals_by_chrom'][chrom_key].append((lo, hi))

    # finalize stats
    for sid in aggregates:
        for cous in aggregates[sid]:
            seg_list = aggregates[sid][cous]['segments']
            if seg_list:
                aggregates[sid][cous]['mean_cm'] = float(np.mean(seg_list))
                aggregates[sid][cous]['median_cm'] = float(np.median(seg_list))
                aggregates[sid][cous]['max_cm'] = float(np.max(seg_list))
            else:
                aggregates[sid][cous]['mean_cm'] = 0.0
                aggregates[sid][cous]['median_cm'] = 0.0
                aggregates[sid][cous]['max_cm'] = 0.0

            merged_total = 0.0
            merged_count = 0
            for _, intervals in aggregates[sid][cous]['intervals_by_chrom'].items():
                merged_intervals = _merge_intervals(intervals)
                merged_count += len(merged_intervals)
                merged_total += sum((e - s) for s, e in merged_intervals)
            aggregates[sid][cous]['merged_interval_count'] = merged_count
            aggregates[sid][cous]['total_cm_merged'] = merged_total

    return aggregates


def aggregate_x_chromosome(feathers, parent_ids, cousins):
    """
    Analyze X chromosome IBD separately (chr 23).
    """
    x_agg = {pid: {} for pid in parent_ids}

    if 23 not in feathers:
        return x_agg

    df = feathers[23]
    if df.empty:
        return x_agg

    for _, row in df.iterrows():
        id1 = str(row['id1'])
        id2 = str(row['id2'])
        segment_cm = float(row['end_cm']) - float(row['start_cm'])

        parent_match = None
        cousin_match = None

        if id1 in parent_ids and id2 in cousins:
            parent_match = id1
            cousin_match = id2
        elif id2 in parent_ids and id1 in cousins:
            parent_match = id2
            cousin_match = id1

        if parent_match and cousin_match:
            if cousin_match not in x_agg[parent_match]:
                x_agg[parent_match][cousin_match] = {
                    'count': 0,
                    'total_cm': 0.0,
                    'segments': []
                }

            x_agg[parent_match][cousin_match]['count'] += 1
            x_agg[parent_match][cousin_match]['total_cm'] += segment_cm
            x_agg[parent_match][cousin_match]['segments'].append(segment_cm)

    # Calculate stats
    for parent_id in x_agg:
        for cousin_id in x_agg[parent_id]:
            seg_list = x_agg[parent_id][cousin_id]['segments']
            if seg_list:
                x_agg[parent_id][cousin_id]['mean_cm'] = np.mean(seg_list)
                x_agg[parent_id][cousin_id]['median_cm'] = np.median(seg_list)
            else:
                x_agg[parent_id][cousin_id]['mean_cm'] = 0.0
                x_agg[parent_id][cousin_id]['median_cm'] = 0.0

    return x_agg


def build_chromosome_totals(aggregates, entity_ids, cousins):
    """Summarize merged cM totals per chromosome for each entity."""
    chrom_totals = {entity_id: defaultdict(float) for entity_id in entity_ids}

    for entity_id in entity_ids:
        entity_data = aggregates.get(entity_id, {})
        for cousin in cousins:
            cousin_data = entity_data.get(cousin)
            if not cousin_data:
                continue
            for chrom_key, intervals in cousin_data.get('intervals_by_chrom', {}).items():
                merged = _merge_intervals(intervals)
                chrom_totals[entity_id][str(chrom_key)] += sum((end - start) for start, end in merged)

    return {
        entity_id: dict(sorted(chrom_totals[entity_id].items(), key=lambda item: int(item[0])))
        for entity_id in entity_ids
    }


def compute_entity_scores(aggregates, entity_ids, cousins, primary_cousins, secondary_weight):
    """Rank any set of entities by weighted cousin IBD evidence."""
    primary_set = set(primary_cousins)
    chrom_totals = build_chromosome_totals(aggregates, entity_ids, cousins)
    scores = []

    for entity_id in entity_ids:
        total_cm = 0.0
        weighted_cm = 0.0
        segment_count = 0
        cousin_details = {}

        for cousin in cousins:
            data = aggregates.get(entity_id, {}).get(cousin)
            if not data:
                continue
            cousin_cm = float(data.get('total_cm_merged', data.get('total_cm_raw', 0.0)))
            total_cm += cousin_cm
            segment_count += int(data.get('merged_interval_count', data.get('count', 0)))
            weight = 1.0 if cousin in primary_set else float(secondary_weight)
            weighted_cm += cousin_cm * weight
            cousin_details[cousin] = data

        scores.append(
            {
                'entity_id': entity_id,
                'total_cm': total_cm,
                'weighted_cm': weighted_cm,
                'segment_count': segment_count,
                'cousin_details': cousin_details,
                'chromosome_totals': chrom_totals.get(entity_id, {}),
            }
        )

    scores.sort(key=lambda item: (item['weighted_cm'], item['total_cm']), reverse=True)
    return scores


def compute_sibling_anchor_summary(aggregates, siblings, cousins, primary_cousins, secondary_weight):
    """Summarize sibling-level cousin evidence and select the strongest anchor."""
    ranking = compute_entity_scores(aggregates, siblings, cousins, primary_cousins, secondary_weight)
    summary = {
        'ranking': ranking,
        'paternal_anchor_sibling': 'UNKNOWN',
        'difference_cm': 0.0,
        'difference_pct': 0.0,
        'confidence': 'none',
    }

    if not ranking:
        return summary

    summary['paternal_anchor_sibling'] = ranking[0]['entity_id']
    if len(ranking) == 1:
        summary['difference_cm'] = ranking[0]['weighted_cm']
        summary['difference_pct'] = 100.0 if ranking[0]['weighted_cm'] > 0 else 0.0
        summary['confidence'] = 'high' if ranking[0]['weighted_cm'] > 0 else 'none'
        return summary

    top = ranking[0]
    second = ranking[1]
    diff = float(top['weighted_cm']) - float(second['weighted_cm'])
    summary['difference_cm'] = diff
    if second['weighted_cm'] > 0:
        summary['difference_pct'] = (diff / float(second['weighted_cm'])) * 100.0
    elif top['weighted_cm'] > 0:
        summary['difference_pct'] = 100.0

    if summary['difference_pct'] > 50:
        summary['confidence'] = 'high'
    elif summary['difference_pct'] > 25:
        summary['confidence'] = 'medium'
    elif top['weighted_cm'] > 0:
        summary['confidence'] = 'low'

    return summary


def diagnose_parent_symmetry(aggregates, parent1, parent2, cousins):
    """Explain when both reconstructed parents carry near-identical cousin signal."""
    chrom_totals = build_chromosome_totals(aggregates, [parent1, parent2], cousins)
    chrom_keys = sorted(
        set(chrom_totals.get(parent1, {}).keys()) | set(chrom_totals.get(parent2, {}).keys()),
        key=lambda value: int(value),
    )

    per_chrom = []
    balanced = 0
    winner_counts = {parent1: 0, parent2: 0, 'TIE': 0}
    total1 = 0.0
    total2 = 0.0

    for chrom_key in chrom_keys:
        total_p1 = float(chrom_totals.get(parent1, {}).get(chrom_key, 0.0))
        total_p2 = float(chrom_totals.get(parent2, {}).get(chrom_key, 0.0))
        total1 += total_p1
        total2 += total_p2

        max_total = max(total_p1, total_p2)
        ratio = 1.0 if max_total <= 0 else min(total_p1, total_p2) / max_total
        is_balanced = max_total > 0 and ratio >= 0.90
        if is_balanced:
            balanced += 1

        diff = abs(total_p1 - total_p2)
        if diff <= 5.0:
            winner = 'TIE'
        elif total_p1 > total_p2:
            winner = parent1
        else:
            winner = parent2
        winner_counts[winner] += 1

        per_chrom.append(
            {
                'chromosome': chrom_key,
                parent1: total_p1,
                parent2: total_p2,
                'ratio_smaller_to_larger': ratio,
                'winner': winner,
                'balanced': is_balanced,
            }
        )

    overall_ratio = 1.0 if max(total1, total2) <= 0 else min(total1, total2) / max(total1, total2)
    balanced_fraction = (balanced / len(chrom_keys)) if chrom_keys else 0.0
    is_symmetric = overall_ratio >= 0.95 or balanced_fraction >= 0.70

    return {
        'chromosome_totals': chrom_totals,
        'per_chromosome': per_chrom,
        'balanced_chromosomes': balanced,
        'total_chromosomes': len(chrom_keys),
        'balanced_fraction': balanced_fraction,
        'overall_ratio': overall_ratio,
        'winner_counts': winner_counts,
        'is_symmetric': is_symmetric,
    }


def compute_parent_comparison(autosomal_agg, parent1, parent2, cousins, primary_cousins, secondary_weight):
    """
    Compare two parents' IBD totals to cousins.
    Returns dict with comparison metrics.
    """
    ranked = compute_entity_scores(
        autosomal_agg,
        [parent1, parent2],
        cousins,
        primary_cousins,
        secondary_weight,
    )
    by_id = {entry['entity_id']: entry for entry in ranked}
    parent1_score = by_id.get(parent1, {'total_cm': 0.0, 'weighted_cm': 0.0, 'segment_count': 0, 'cousin_details': {}, 'chromosome_totals': {}})
    parent2_score = by_id.get(parent2, {'total_cm': 0.0, 'weighted_cm': 0.0, 'segment_count': 0, 'cousin_details': {}, 'chromosome_totals': {}})

    results = {
        'parent1': parent1,
        'parent2': parent2,
        'parent1_total_cm': float(parent1_score['total_cm']),
        'parent2_total_cm': float(parent2_score['total_cm']),
        'parent1_weighted_cm': float(parent1_score['weighted_cm']),
        'parent2_weighted_cm': float(parent2_score['weighted_cm']),
        'parent1_segment_count': int(parent1_score['segment_count']),
        'parent2_segment_count': int(parent2_score['segment_count']),
        'parent1_cousin_details': parent1_score['cousin_details'],
        'parent2_cousin_details': parent2_score['cousin_details'],
        'parent1_chromosome_totals': parent1_score['chromosome_totals'],
        'parent2_chromosome_totals': parent2_score['chromosome_totals'],
        'difference_cm': 0.0,
        'difference_pct': 0.0,
        'paternal_parent': None,
        'maternal_parent': None,
        'confidence': 'low'
    }

    # Determine paternal parent and confidence
    if results['parent1_weighted_cm'] == 0 and results['parent2_weighted_cm'] == 0:
        results['paternal_parent'] = 'UNKNOWN'
        results['maternal_parent'] = 'UNKNOWN'
        results['difference_cm'] = 0.0
        results['difference_pct'] = 0.0
        results['confidence'] = 'none'
        return results

    if results['parent1_weighted_cm'] > results['parent2_weighted_cm']:
        results['paternal_parent'] = parent1
        results['maternal_parent'] = parent2
        difference = results['parent1_weighted_cm'] - results['parent2_weighted_cm']
        results['difference_cm'] = difference
        if results['parent2_weighted_cm'] > 0:
            results['difference_pct'] = (difference / results['parent2_weighted_cm']) * 100
        else:
            results['difference_pct'] = 100.0
    elif results['parent2_weighted_cm'] > results['parent1_weighted_cm']:
        results['paternal_parent'] = parent2
        results['maternal_parent'] = parent1
        difference = results['parent2_weighted_cm'] - results['parent1_weighted_cm']
        results['difference_cm'] = difference
        if results['parent1_weighted_cm'] > 0:
            results['difference_pct'] = (difference / results['parent1_weighted_cm']) * 100
        else:
            results['difference_pct'] = 100.0
    else:
        # Tie case with non-zero data: report as inconclusive.
        results['paternal_parent'] = 'UNKNOWN'
        results['maternal_parent'] = 'UNKNOWN'
        results['difference_cm'] = 0.0
        results['difference_pct'] = 0.0
        results['confidence'] = 'low'
        return results

    # Assign confidence based on difference
    if results['difference_pct'] > 50:
        results['confidence'] = 'high'
    elif results['difference_pct'] > 25:
        results['confidence'] = 'medium'
    else:
        results['confidence'] = 'low'

    return results


def determine_propagation_status(comparison, sibling_summary, symmetry):
    """Convert combined evidence into an upstream-friendly suggested role map."""
    suggested_paternal = comparison.get('paternal_parent')
    suggested_maternal = comparison.get('maternal_parent')

    payload = {
        'status': 'unresolved',
        'suggested_paternal_parent': suggested_paternal,
        'suggested_maternal_parent': suggested_maternal,
        'confidence': comparison.get('confidence', 'none'),
        'reason': '',
    }

    if suggested_paternal in (None, 'UNKNOWN') or suggested_maternal in (None, 'UNKNOWN'):
        payload['reason'] = 'Parent-cousin analysis did not produce a usable paternal candidate.'
        return payload

    if comparison['confidence'] in ('high', 'medium') and not symmetry.get('is_symmetric', False):
        payload['status'] = 'confirmed' if comparison['confidence'] == 'high' else 'provisional'
        payload['reason'] = 'Parent-cousin separation is not symmetric and the reconstructed-parent score is decisive enough to use upstream.'
        return payload

    if sibling_summary and sibling_summary.get('confidence') == 'high':
        payload['status'] = 'low-confidence'
        payload['reason'] = (
            'Sibling-cousin evidence is strong, but reconstructed parents remain symmetric across cousin IBD. '
            'Use the suggested parent labels only as provisional aliases.'
        )
        return payload

    payload['reason'] = 'Both reconstructed parents show near-identical cousin IBD totals, so father/mother propagation is unresolved.'
    return payload


def compute_chromosome_local_labels(
    aggregates,
    parent1,
    parent2,
    cousins,
    min_diff_cm=10.0,
    min_diff_pct=5.0,
):
    """Assign a provisional paternal parent independently per chromosome."""
    chrom_totals = build_chromosome_totals(aggregates, [parent1, parent2], cousins)
    chrom_keys = sorted(
        set(chrom_totals.get(parent1, {}).keys()) | set(chrom_totals.get(parent2, {}).keys()),
        key=lambda value: int(value),
    )

    rows = []
    confident_counts = {parent1: 0, parent2: 0}
    weak_counts = {parent1: 0, parent2: 0}
    ambiguous = 0
    no_data = 0

    for chrom_key in chrom_keys:
        total1 = float(chrom_totals.get(parent1, {}).get(chrom_key, 0.0))
        total2 = float(chrom_totals.get(parent2, {}).get(chrom_key, 0.0))

        if total1 <= 0 and total2 <= 0:
            row = {
                'chromosome': chrom_key,
                parent1: total1,
                parent2: total2,
                'winner': 'NO_DATA',
                'difference_cm': 0.0,
                'difference_pct': 0.0,
                'confidence': 'none',
                'suggested_paternal_parent': 'UNKNOWN',
            }
            rows.append(row)
            no_data += 1
            continue

        diff = abs(total1 - total2)
        bigger = max(total1, total2)
        smaller = min(total1, total2)
        if smaller > 0:
            diff_pct = (diff / smaller) * 100.0
        else:
            diff_pct = 100.0 if diff > 0 else 0.0

        if diff <= 0:
            winner = 'TIE'
            confidence = 'ambiguous'
            suggested = 'UNKNOWN'
            ambiguous += 1
        else:
            winner = parent1 if total1 > total2 else parent2
            if diff >= float(min_diff_cm) and diff_pct >= float(min_diff_pct):
                confidence = 'confident'
                suggested = winner
                confident_counts[winner] += 1
            else:
                confidence = 'weak'
                suggested = winner
                weak_counts[winner] += 1

        rows.append(
            {
                'chromosome': chrom_key,
                parent1: total1,
                parent2: total2,
                'winner': winner,
                'difference_cm': diff,
                'difference_pct': diff_pct,
                'confidence': confidence,
                'suggested_paternal_parent': suggested,
            }
        )

    if confident_counts[parent1] > 0 and confident_counts[parent2] == 0:
        consistent_with = parent1
    elif confident_counts[parent2] > 0 and confident_counts[parent1] == 0:
        consistent_with = parent2
    elif confident_counts[parent1] == 0 and confident_counts[parent2] == 0:
        consistent_with = 'UNKNOWN'
    else:
        consistent_with = 'MIXED'

    return {
        'cousins_used': list(cousins),
        'min_diff_cm': float(min_diff_cm),
        'min_diff_pct': float(min_diff_pct),
        'per_chromosome': rows,
        'confident_counts': confident_counts,
        'weak_counts': weak_counts,
        'ambiguous_count': ambiguous,
        'no_data_count': no_data,
        'total_chromosomes': len(chrom_keys),
        'consistent_paternal_parent': consistent_with,
    }


def validate_x_chromosome(x_agg, paternal_parent, female_cousins):
    """
    Validate paternal assignment using X-chromosome.
    Female cousins can share X with both parents; strong X match to paternal
    parent from Diane (female sibling) confirms paternal assignment.
    """
    validation = {
        'x_paternal_match_cm': 0.0,
        'x_maternal_match_cm': 0.0,
        'x_confidence': 'unknown',
        'notes': []
    }

    if paternal_parent == 'UNKNOWN':
        validation['notes'].append("X-chromosome validation skipped because paternal parent is unknown")
        return validation

    if not x_agg or paternal_parent not in x_agg:
        validation['notes'].append("No X-chromosome data available for validation")
        return validation

    paternal_match_cm = 0.0

    for cousin in female_cousins:
        if cousin in x_agg[paternal_parent]:
            paternal_match_cm += x_agg[paternal_parent][cousin]['total_cm']

    validation['x_paternal_match_cm'] = paternal_match_cm

    if paternal_match_cm > 10:
        validation['x_confidence'] = 'high'
        validation['notes'].append(f"Strong X-chromosome match to paternal parent ({paternal_match_cm:.1f} cM)")
    elif paternal_match_cm > 5:
        validation['x_confidence'] = 'medium'
        validation['notes'].append(f"Moderate X-chromosome match to paternal parent ({paternal_match_cm:.1f} cM)")
    else:
        validation['x_confidence'] = 'low'
        validation['notes'].append(f"Weak X-chromosome match to paternal parent ({paternal_match_cm:.1f} cM)")

    return validation


def generate_report(
    comparison,
    x_validation,
    output_file,
    cousins,
    sibling_summary=None,
    symmetry=None,
    propagation=None,
    chromosome_labels=None,
):
    """Generate human-readable sexing report."""
    report_lines = []

    report_lines.append("=" * 80)
    report_lines.append("RECONSTRUCTED PARENT SEXING ANALYSIS REPORT")
    report_lines.append("=" * 80)
    report_lines.append("")

    # Summary
    report_lines.append("SUMMARY")
    report_lines.append("-" * 80)
    report_lines.append(f"Paternal Parent (FATHER):    {comparison['paternal_parent']}")
    report_lines.append(f"Maternal Parent (MOTHER):    {comparison['maternal_parent']}")
    report_lines.append(f"Confidence Level:            {comparison['confidence'].upper()}")
    report_lines.append("")

    if sibling_summary is not None:
        report_lines.append("SIBLING ANCHOR EVIDENCE")
        report_lines.append("-" * 80)
        report_lines.append(f"Paternal Anchor Sibling:     {sibling_summary['paternal_anchor_sibling']}")
        report_lines.append(f"Sibling Anchor Confidence:   {sibling_summary['confidence'].upper()}")
        if sibling_summary.get('override'):
            report_lines.append("Sibling Anchor Source:       EXPLICIT OVERRIDE")
        report_lines.append(f"Sibling Anchor Margin:       {sibling_summary['difference_cm']:.2f} cM ({sibling_summary['difference_pct']:.1f}%)")
        if sibling_summary.get('ranking'):
            report_lines.append("Sibling Weighted Ranking:")
            for entry in sibling_summary['ranking']:
                report_lines.append(
                    f"  {entry['entity_id']:15s} Weighted: {entry['weighted_cm']:8.2f} cM  Total: {entry['total_cm']:8.2f} cM  Segments: {entry['segment_count']:4d}"
                )
        report_lines.append("")

    # Autosomal Analysis
    report_lines.append("AUTOSOMAL IBD ANALYSIS (Chr 1-22)")
    report_lines.append("-" * 80)
    report_lines.append(f"{comparison['parent1']:20s} Total cM (merged): {comparison['parent1_total_cm']:10.2f}  Segments: {comparison['parent1_segment_count']:4d}")
    report_lines.append(f"{comparison['parent2']:20s} Total cM (merged): {comparison['parent2_total_cm']:10.2f}  Segments: {comparison['parent2_segment_count']:4d}")
    report_lines.append(f"{comparison['parent1']:20s} Weighted cM:       {comparison['parent1_weighted_cm']:10.2f}")
    report_lines.append(f"{comparison['parent2']:20s} Weighted cM:       {comparison['parent2_weighted_cm']:10.2f}")
    report_lines.append(f"Weighted Difference:         {comparison['difference_cm']:10.2f} cM ({comparison['difference_pct']:.1f}%)")
    report_lines.append("")

    if symmetry is not None:
        report_lines.append("PARENT SYMMETRY DIAGNOSTICS")
        report_lines.append("-" * 80)
        report_lines.append(f"Balanced Chromosomes:        {symmetry['balanced_chromosomes']}/{symmetry['total_chromosomes']} ({symmetry['balanced_fraction'] * 100.0:.1f}%)")
        report_lines.append(f"Overall Smaller/Larger Ratio:{symmetry['overall_ratio']:.3f}")
        report_lines.append(f"Chromosome Winner Counts:    {comparison['parent1']}={symmetry['winner_counts'][comparison['parent1']]}  {comparison['parent2']}={symmetry['winner_counts'][comparison['parent2']]}  TIE={symmetry['winner_counts']['TIE']}")
        report_lines.append(f"Symmetric Parent Signal:     {'YES' if symmetry['is_symmetric'] else 'NO'}")
        report_lines.append("Chromosome Winners:")
        for row in symmetry['per_chromosome']:
            report_lines.append(
                f"  Chr{row['chromosome']:>2s}: {comparison['parent1']}={row[comparison['parent1']]:7.2f}  {comparison['parent2']}={row[comparison['parent2']]:7.2f}  winner={row['winner']}  ratio={row['ratio_smaller_to_larger']:.3f}"
            )
        report_lines.append("")

    if chromosome_labels is not None:
        report_lines.append("CHROMOSOME-LOCAL PARENT LABELING")
        report_lines.append("-" * 80)
        report_lines.append(f"Cousins Used:                {', '.join(chromosome_labels['cousins_used'])}")
        report_lines.append(f"Thresholds:                  diff >= {chromosome_labels['min_diff_cm']:.1f} cM and diff >= {chromosome_labels['min_diff_pct']:.1f}%")
        report_lines.append(
            f"Confident Winners:           {comparison['parent1']}={chromosome_labels['confident_counts'][comparison['parent1']]}  "
            f"{comparison['parent2']}={chromosome_labels['confident_counts'][comparison['parent2']]}"
        )
        report_lines.append(
            f"Weak Winners:                {comparison['parent1']}={chromosome_labels['weak_counts'][comparison['parent1']]}  "
            f"{comparison['parent2']}={chromosome_labels['weak_counts'][comparison['parent2']]}"
        )
        report_lines.append(
            f"Ambiguous/No-data:           {chromosome_labels['ambiguous_count']} / {chromosome_labels['no_data_count']}"
        )
        report_lines.append(f"Consistent Paternal Parent:  {chromosome_labels['consistent_paternal_parent']}")
        report_lines.append("Per-chromosome calls:")
        for row in chromosome_labels['per_chromosome']:
            report_lines.append(
                f"  Chr{row['chromosome']:>2s}: "
                f"{comparison['parent1']}={row[comparison['parent1']]:7.2f} "
                f"{comparison['parent2']}={row[comparison['parent2']]:7.2f} "
                f"winner={row['winner']} conf={row['confidence']} "
                f"diff={row['difference_cm']:.2f}cM ({row['difference_pct']:.1f}%)"
            )
        report_lines.append("")

    # Per-cousin breakdown for paternal parent
    report_lines.append(f"IBD Details for PATERNAL PARENT ({comparison['paternal_parent']}):")
    report_lines.append("-" * 80)
    if comparison['parent1'] == comparison['paternal_parent']:
        paternal_details = comparison['parent1_cousin_details']
    elif comparison['parent2'] == comparison['paternal_parent']:
        paternal_details = comparison['parent2_cousin_details']
    else:
        paternal_details = {}

    for cousin in cousins:
        if cousin in paternal_details:
            data = paternal_details[cousin]
            merged_cm = float(data.get('total_cm_merged', data.get('total_cm_raw', 0.0)))
            merged_count = int(data.get('merged_interval_count', data.get('count', 0)))
            report_lines.append(f"  {cousin:15s} Count: {merged_count:4d}  Total: {merged_cm:8.2f} cM  Mean(raw): {data['mean_cm']:6.2f} cM")
        else:
            report_lines.append(f"  {cousin:15s} Count: {0:4d}  Total: {0.0:8.2f} cM  Mean: {0.0:6.2f} cM")
    report_lines.append("")

    # X-chromosome validation
    report_lines.append("X-CHROMOSOME VALIDATION (Chr 23)")
    report_lines.append("-" * 80)
    report_lines.append(f"X-chr Match to Paternal Parent: {x_validation['x_paternal_match_cm']:.2f} cM")
    report_lines.append(f"X-chr Confidence:               {x_validation['x_confidence'].upper()}")
    for note in x_validation['notes']:
        report_lines.append(f"  - {note}")
    report_lines.append("")

    # Recommendation
    report_lines.append("RECOMMENDATION")
    report_lines.append("-" * 80)
    if propagation is not None:
        report_lines.append(f"Propagation Status:          {propagation['status'].upper()}")
        report_lines.append(f"Propagation Reason:          {propagation['reason']}")
        if propagation.get('suggested_paternal_parent') not in (None, 'UNKNOWN'):
            report_lines.append(f"Suggested Paternal Parent:   {propagation['suggested_paternal_parent']}")
            report_lines.append(f"Suggested Maternal Parent:   {propagation['suggested_maternal_parent']}")
        report_lines.append("")

    if comparison['confidence'] == 'none':
        report_lines.append("[WARN] NO MATCHING PARENT-COUSIN SEGMENTS FOUND")
        report_lines.append("[WARN] Verify that --parent1/--parent2 IDs appear in id1/id2 columns of the feather files.")
    elif comparison['confidence'] == 'high':
        report_lines.append(f"[OK] HIGH CONFIDENCE: {comparison['paternal_parent']} is the PATERNAL PARENT")
    elif comparison['confidence'] == 'medium':
        report_lines.append(f"[WARN] MEDIUM CONFIDENCE: {comparison['paternal_parent']} is likely the PATERNAL PARENT")
    else:
        report_lines.append(f"[WARN] LOW CONFIDENCE: Results are inconclusive. Manual review recommended.")
    report_lines.append("")

    report_lines.append("=" * 80)

    report_text = "\n".join(report_lines)

    # Print to console
    print(report_text, flush=True)

    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report_text)

    print(f"\nReport saved to: {output_file}", flush=True)


def write_role_map(output_file, comparison, sibling_summary, symmetry, propagation, chromosome_labels=None):
    """Write a machine-readable role suggestion payload for upstream use."""
    payload = {
        'reconstructed_parent_comparison': {
            'parent1': comparison['parent1'],
            'parent2': comparison['parent2'],
            'paternal_parent': comparison['paternal_parent'],
            'maternal_parent': comparison['maternal_parent'],
            'confidence': comparison['confidence'],
            'difference_cm': comparison['difference_cm'],
            'difference_pct': comparison['difference_pct'],
        },
        'sibling_anchor': sibling_summary,
        'symmetry_diagnostics': symmetry,
        'chromosome_local_labeling': chromosome_labels,
        'propagation': propagation,
        'roles': [
            {
                'entity_id': comparison['parent1'],
                'suggested_role': 'father' if propagation.get('suggested_paternal_parent') == comparison['parent1'] else 'mother' if propagation.get('suggested_maternal_parent') == comparison['parent1'] else 'unknown',
                'status': propagation.get('status', 'unresolved'),
            },
            {
                'entity_id': comparison['parent2'],
                'suggested_role': 'father' if propagation.get('suggested_paternal_parent') == comparison['parent2'] else 'mother' if propagation.get('suggested_maternal_parent') == comparison['parent2'] else 'unknown',
                'status': propagation.get('status', 'unresolved'),
            },
        ],
    }

    with open(output_file, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=2)

    print(f"Role map saved to: {output_file}", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze reconstructed parent sex using cousin IBD segments"
    )
    parser.add_argument(
        '--ibd-dir',
        required=True,
        help='Directory containing feather IBD files (chr1.feather, chr2.feather, ...)'
    )
    parser.add_argument(
        '--parent1',
        required=True,
        help='First reconstructed parent ID'
    )
    parser.add_argument(
        '--parent2',
        required=True,
        help='Second reconstructed parent ID'
    )
    parser.add_argument(
        '--cousins',
        required=True,
        help='Comma-separated list of paternal cousin IDs (e.g., "Sue,Wendy,LS")'
    )
    parser.add_argument(
        '--female-cousins',
        default='',
        help='Comma-separated list of female cousins for X-validation (e.g., "Sue,Wendy,LS")'
    )
    parser.add_argument(
        '--primary-cousins',
        default='',
        help='Comma-separated cousins to weight as primary evidence (e.g., "Sue,Wendy"); others get --secondary-weight'
    )
    parser.add_argument(
        '--secondary-weight',
        type=float,
        default=0.5,
        help='Weight for non-primary cousins in autosomal scoring (default: 0.5)'
    )
    parser.add_argument(
        '--out-report',
        default='parent_sexing_report.txt',
        help='Output report file name'
    )
    parser.add_argument(
        '--analysis-mode',
        choices=['parent-cousin', 'sibling-cousin'],
        default='parent-cousin',
        help='Analyze either reconstructed parents vs cousins or sibling vs cousin pairwise evidence'
    )
    parser.add_argument(
        '--siblings',
        default='',
        help='Comma-separated sibling IDs for --analysis-mode sibling-cousin (e.g., "Ray,Tom,Diane")'
    )
    parser.add_argument(
        '--min-map',
        default='',
        help='Path to min_map.txt (required for --analysis-mode sibling-cousin)'
    )
    parser.add_argument(
        '--bim',
        default='',
        help='Path to BIM file used for HAPI/HAPI-RECAP marker indexing (required for --analysis-mode sibling-cousin)'
    )
    parser.add_argument(
        '--sibling-ibd-dir',
        default='',
        help='Optional sibling-cousin feather directory used to anchor parent labeling and diagnose parent symmetry'
    )
    parser.add_argument(
        '--write-role-map',
        default='',
        help='Optional JSON output path for upstream role propagation metadata'
    )
    parser.add_argument(
        '--paternal-anchor-sibling',
        default='',
        help='Optional trusted sibling anchor override (for example, "Ray") when sibling-cousin evidence was already validated separately'
    )
    parser.add_argument(
        '--anchor-confidence',
        choices=['none', 'low', 'medium', 'high'],
        default='high',
        help='Confidence attached to --paternal-anchor-sibling (default: high)'
    )
    parser.add_argument(
        '--chromosome-cousins',
        default='',
        help='Optional comma-separated cousin IDs to use for chromosome-local labeling (default: --primary-cousins)'
    )
    parser.add_argument(
        '--chromosome-min-diff-cm',
        type=float,
        default=10.0,
        help='Minimum cM difference to call a chromosome-level winner as confident (default: 10.0)'
    )
    parser.add_argument(
        '--chromosome-min-diff-pct',
        type=float,
        default=5.0,
        help='Minimum percent difference to call a chromosome-level winner as confident (default: 5.0)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.isdir(args.ibd_dir):
        print(f"ERROR: IBD directory not found: {args.ibd_dir}", flush=True)
        sys.exit(1)

    parent_ids = [args.parent1, args.parent2]
    cousins = [c.strip() for c in args.cousins.split(',')]
    female_cousins = [c.strip() for c in args.female_cousins.split(',')] if args.female_cousins else cousins
    primary_cousins = [c.strip() for c in args.primary_cousins.split(',') if c.strip()] if args.primary_cousins else cousins

    print(f"Loading IBD feather files from: {args.ibd_dir}", flush=True)
    print(f"Parents: {args.parent1}, {args.parent2}", flush=True)
    print(f"Paternal Cousins: {', '.join(cousins)}", flush=True)
    print(f"Female Cousins: {', '.join(female_cousins)}", flush=True)
    print(f"Primary Cousins (weight 1.0): {', '.join(primary_cousins)}", flush=True)
    print(f"Secondary Cousin Weight: {args.secondary_weight}", flush=True)
    print()

    # Load feathers: chromosomes 1-23 (23 = X)
    feathers = load_feather_files(args.ibd_dir, list(range(1, 24)))

    if not feathers:
        print("ERROR: No feather files loaded!", flush=True)
        sys.exit(1)

    # Separate autosomal (1-22) from X (23)
    autosomal_feathers = {c: feathers[c] for c in feathers if c != 23}

    print()
    print("Aggregating IBD by parent and cousin...", flush=True)

    # Aggregate autosomal IBD
    if args.analysis_mode == 'sibling-cousin':
        if not args.siblings or not args.min_map or not args.bim:
            print('ERROR: --analysis-mode sibling-cousin requires --siblings, --min-map, and --bim', flush=True)
            sys.exit(1)
        siblings = [s.strip() for s in args.siblings.split(',') if s.strip()]
        dmap_positions, dmap_cms = _load_min_map(args.min_map)
        bim_positions = _load_bim_positions(args.bim)
        autosomal_agg = aggregate_sibling_cousin_ibd(
            autosomal_feathers,
            siblings,
            cousins,
            bim_positions,
            dmap_positions,
            dmap_cms,
        )

        # In sibling-cousin mode, compare top sibling evidence as proxy and map to report slots.
        sibling_scores = []
        for sid in siblings:
            total = 0.0
            for cous in cousins:
                if cous in autosomal_agg.get(sid, {}):
                    total += float(autosomal_agg[sid][cous].get('total_cm_merged', 0.0))
            sibling_scores.append((sid, total))
        sibling_scores.sort(key=lambda x: x[1], reverse=True)
        if len(sibling_scores) >= 2:
            args.parent1, args.parent2 = sibling_scores[0][0], sibling_scores[1][0]
        elif len(sibling_scores) == 1:
            args.parent1, args.parent2 = sibling_scores[0][0], 'UNKNOWN'
        else:
            args.parent1, args.parent2 = 'UNKNOWN', 'UNKNOWN'
        parent_ids = [args.parent1, args.parent2]
    else:
        autosomal_agg = aggregate_ibd_by_parent(autosomal_feathers, parent_ids, cousins)

    sibling_summary = None
    if args.sibling_ibd_dir:
        if not os.path.isdir(args.sibling_ibd_dir):
            print(f"ERROR: sibling IBD directory not found: {args.sibling_ibd_dir}", flush=True)
            sys.exit(1)
        if not args.siblings or not args.min_map or not args.bim:
            print('ERROR: --sibling-ibd-dir requires --siblings, --min-map, and --bim', flush=True)
            sys.exit(1)

        siblings = [s.strip() for s in args.siblings.split(',') if s.strip()]
        dmap_positions, dmap_cms = _load_min_map(args.min_map)
        bim_positions = _load_bim_positions(args.bim)
        sibling_feathers = load_feather_files(args.sibling_ibd_dir, list(range(1, 24)))
        sibling_autosomal = {c: sibling_feathers[c] for c in sibling_feathers if c != 23}
        sibling_agg = aggregate_sibling_cousin_ibd(
            sibling_autosomal,
            siblings,
            cousins,
            bim_positions,
            dmap_positions,
            dmap_cms,
        )
        sibling_summary = compute_sibling_anchor_summary(
            sibling_agg,
            siblings,
            cousins,
            primary_cousins,
            args.secondary_weight,
        )

    if args.paternal_anchor_sibling:
        if sibling_summary is None:
            sibling_summary = {
                'ranking': [],
                'paternal_anchor_sibling': args.paternal_anchor_sibling,
                'difference_cm': 0.0,
                'difference_pct': 0.0,
                'confidence': args.anchor_confidence,
                'override': True,
            }
        else:
            sibling_summary['paternal_anchor_sibling'] = args.paternal_anchor_sibling
            sibling_summary['confidence'] = args.anchor_confidence
            sibling_summary['override'] = True

    chromosome_cousins = [
        c.strip() for c in args.chromosome_cousins.split(',') if c.strip()
    ] if args.chromosome_cousins else primary_cousins

    # Helpful diagnostic: show which IDs are present in the loaded feather files.
    observed_id1 = set()
    observed_id2 = set()
    for _, chrom_df in autosomal_feathers.items():
        if chrom_df.empty:
            continue
        observed_id1.update(chrom_df['id1'].astype(str).tolist())
        observed_id2.update(chrom_df['id2'].astype(str).tolist())
    observed_ids = sorted(observed_id1.union(observed_id2))
    print(f"Observed IDs in feather files (sample up to 20): {observed_ids[:20]}", flush=True)

    # Compare parents
    comparison = compute_parent_comparison(
        autosomal_agg,
        args.parent1,
        args.parent2,
        cousins,
        primary_cousins,
        args.secondary_weight,
    )

    # Aggregate and validate X-chromosome
    x_agg = aggregate_x_chromosome(feathers, parent_ids, female_cousins)
    x_validation = validate_x_chromosome(x_agg, comparison['paternal_parent'], female_cousins)
    symmetry = diagnose_parent_symmetry(autosomal_agg, args.parent1, args.parent2, cousins)
    chromosome_labels = compute_chromosome_local_labels(
        autosomal_agg,
        args.parent1,
        args.parent2,
        chromosome_cousins,
        min_diff_cm=args.chromosome_min_diff_cm,
        min_diff_pct=args.chromosome_min_diff_pct,
    )
    propagation = determine_propagation_status(comparison, sibling_summary, symmetry)

    # Generate and save report
    print()
    generate_report(
        comparison,
        x_validation,
        args.out_report,
        cousins,
        sibling_summary=sibling_summary,
        symmetry=symmetry,
        propagation=propagation,
        chromosome_labels=chromosome_labels,
    )

    if args.write_role_map:
        write_role_map(args.write_role_map, comparison, sibling_summary, symmetry, propagation, chromosome_labels=chromosome_labels)

    # Return exit code based on confidence
    if comparison['confidence'] == 'high':
        sys.exit(0)
    elif comparison['confidence'] == 'medium':
        sys.exit(0)  # Still success, but user should review
    else:
        sys.exit(1)  # Low confidence or no evidence


if __name__ == '__main__':
    main()
