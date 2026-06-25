#!/usr/bin/env python3
import argparse

from analyze_reconstructed_parents_sex import (
    _load_bim_positions,
    _load_min_map,
    aggregate_sibling_cousin_ibd,
    compute_entity_scores,
    load_feather_files,
)


def confidence_from_scores(ranking):
    if not ranking:
        return ("none", 0.0, 0.0, "UNKNOWN", "UNKNOWN")
    if len(ranking) == 1:
        top = ranking[0]
        return ("high" if top["weighted_cm"] > 0 else "none", float(top["weighted_cm"]), 100.0 if top["weighted_cm"] > 0 else 0.0, top["entity_id"], "UNKNOWN")

    top = ranking[0]
    second = ranking[1]
    diff = float(top["weighted_cm"] - second["weighted_cm"])
    if second["weighted_cm"] > 0:
        diff_pct = (diff / float(second["weighted_cm"])) * 100.0
    elif top["weighted_cm"] > 0:
        diff_pct = 100.0
    else:
        diff_pct = 0.0

    if diff_pct > 50:
        conf = "high"
    elif diff_pct > 25:
        conf = "medium"
    elif top["weighted_cm"] > 0:
        conf = "low"
    else:
        conf = "none"
    return (conf, diff, diff_pct, top["entity_id"], second["entity_id"])


def main():
    parser = argparse.ArgumentParser(description="Per-chromosome sibling+cousin confidence report")
    parser.add_argument("--ibd-dir", required=True)
    parser.add_argument("--bim", required=True)
    parser.add_argument("--min-map", required=True)
    parser.add_argument("--siblings", required=True)
    parser.add_argument("--cousins", required=True)
    parser.add_argument("--primary-cousins", default="")
    parser.add_argument("--secondary-weight", type=float, default=0.5)
    args = parser.parse_args()

    siblings = [s.strip() for s in args.siblings.split(",") if s.strip()]
    cousins = [c.strip() for c in args.cousins.split(",") if c.strip()]
    primary = [c.strip() for c in args.primary_cousins.split(",") if c.strip()] if args.primary_cousins else cousins

    dmap_positions, dmap_cms = _load_min_map(args.min_map)
    bim_positions = _load_bim_positions(args.bim)

    print("Chr  Winner  Confidence  Diff(cM)  Diff(%)  TopWeighted  RunnerUp  Ranking")
    print("---  ------  ----------  --------  -------  -----------  --------  -------")

    conf_counts = {"high": 0, "medium": 0, "low": 0, "none": 0}
    winner_counts = {sid: 0 for sid in siblings}

    for chrom in range(1, 23):
        feathers = load_feather_files(args.ibd_dir, [chrom])
        if chrom not in feathers:
            print(f"{chrom:>3}  MISSING  none       0.00      0.0      0.00        UNKNOWN   -")
            conf_counts["none"] += 1
            continue

        agg = aggregate_sibling_cousin_ibd(
            {chrom: feathers[chrom]},
            siblings,
            cousins,
            bim_positions,
            dmap_positions,
            dmap_cms,
        )
        ranking = compute_entity_scores(agg, siblings, cousins, primary, args.secondary_weight)
        conf, diff, diff_pct, top_id, second_id = confidence_from_scores(ranking)
        conf_counts[conf] = conf_counts.get(conf, 0) + 1
        if top_id in winner_counts:
            winner_counts[top_id] += 1

        top_weighted = float(ranking[0]["weighted_cm"]) if ranking else 0.0
        rank_txt = ", ".join(f"{r['entity_id']}:{r['weighted_cm']:.1f}" for r in ranking)
        print(f"{chrom:>3}  {top_id:<6}  {conf:<10}  {diff:>8.2f}  {diff_pct:>7.1f}  {top_weighted:>11.2f}  {second_id:<8}  {rank_txt}")

    print("\nSummary")
    print("-------")
    print("Confidence counts:", ", ".join(f"{k}={v}" for k, v in conf_counts.items()))
    print("Winner counts:", ", ".join(f"{k}={v}" for k, v in winner_counts.items()))


if __name__ == "__main__":
    main()
