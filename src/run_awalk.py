"""A script to run an adaptive walk to find an mRNA coding sequence (CDS) for a given protein."""
import pickle
import protein
import awalk
import objective_functions as objectives
import argparse
import vienna


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=1000,
                    help='Number of steps in the adaptive walk')
    ap.add_argument('--aa_seq', type=str, required=True,
                    help='Amino acid sequence to find a CDS for. A string using the standard one-letter code.')
    ap.add_argument('--stability', type=str, choices=['aup', 'efe', 'none'], default='aup',
                    help='Stability objective to use. Set to none to only use CAI')
    ap.add_argument('--freq_table_path', type=str,
                    default='../codon_tables/homosapiens.txt', help='Path to codon frequency table')
    ap.add_argument('--cai_threshold', type=float, default=0.8,
                    help='Objective function forces CAI to be at least this')
    ap.add_argument('--cai_exp_scale', type=float, default=1.0,
                    help='Scaling factor for CAI. Increase to make CAI more important')
    ap.add_argument("--save_path", type=str, default=None,
                    help="The path to save the result. Saved in pickle format.")
    ap.add_argument("--verbose", action="store_true", help="Log all progress")
    ap.add_argument("--load_path", type=str, default=None,
                    help="Loads the initial CDS from the given pickle file. If not specified, the initial CDS is randomly generated.")
    args = ap.parse_args()

    # Load frequency table
    freq_table = protein.CodonFrequencyTable(args.freq_table_path)

    obj_config = objectives.CAIThresholdObjectiveConfig(
        freq_table, args.cai_threshold, args.cai_exp_scale, verbose=args.verbose)

    # Get obj function
    if args.stability == 'aup':
        obj = objectives.make_cai_and_aup_obj(obj_config)
    elif args.stability == 'efe':
        obj = objectives.make_cai_and_efe_obj(obj_config)
    elif args.stability == 'none':
        obj = objectives.make_cai_threshold_obj(obj_config)

    # Load initial CDS if specified
    init_cds = None
    if args.load_path is not None:
        # Read previous result as pickle
        with open(args.load_path, "rb") as f:
            init_cds = pickle.load(f).cds

    # Create walk config
    walk_config = awalk.WalkConfig(
        args.aa_seq, freq_table, obj, args.steps, init_cds=init_cds, verbose=args.verbose)

    # Run walker
    res = awalk.adapative_random_walk(walk_config)

    # Output results
    cai = freq_table.codon_adaptation_index(res.cds)
    aup, efe = vienna.aup_and_efe(vienna.cds_to_rna(res.cds))
    print(
        f"Adaptive walk result for {args.aa_seq} with {args.stability} objective and {args.steps} steps:")
    print(
        f"Result CDS: {res.cds} \n\t Fitness: {res.fitness}, CAI: {cai}, AUP: {aup}, EFE: {efe}")

    # Save the result as pickle
    if args.save_path is not None:
        with open(args.save_path, "wb") as f:
            pickle.dump(res, f)


if __name__ == '__main__':
    main()
