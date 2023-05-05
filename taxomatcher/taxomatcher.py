import argparse
import sys 

def main():    
    ap = argparse.ArgumentParser(description="Taxomatcher")
    subparsers = ap.add_subparsers(title="mode", dest="mode")
    gbif_parser = subparsers.add_parser("gbif")
    gbif_parser.add_argument(
        "-i",
        "--csv",
        dest="input_csv",
        required=True,
        help="CSV where first column is a list of target species names to look up.",
    )
    gbif_parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        help="Path to output.",
    )
    gbif_parser.add_argument("-t", "--threads", dest="threads", required=False, default=4)
    ncbi_parser = subparsers.add_parser("ncbi")
    ncbi_parser.add_argument(
        "-i",
        "--csv",
        dest="input_csv",
        required=True,
        help="CSV where first column is a list of target species names to look up.",
    )
    ncbi_parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        help="Path to output.",
    )
    ncbi_parser.add_argument("-e", "--email", dest="email", required=True)

    trait_parser = subparsers.add_parser("trait")
    trait_parser.add_argument(
        "-t",
        "--traitfile",
        dest="traitfile",
        required=True,
        help="CSV of trait values, first column must be species names.",
    )
    trait_parser.add_argument(
        "-s",
        "--speciesfile",
        dest="speciesfile",
        required=True,
        help="CSV of species synonyms output by the gbif or ncbi modules.",
    )
    trait_parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        help="Path to output.",
    )
    ua = ap.parse_args()
    
    if ua.mode == "gbif":
        from . import gbif_matcher      
        gbif_matcher.main(ua.input_csv, ua.outfile, ua.threads)
    elif ua.mode == "ncbi":
        from . import entrez_matcher
        entrez_matcher.main(ua.input_csv, ua.outfile, ua.email)
    elif ua.mode == "trait":
        from . import trait_matcher
        trait_matcher.main(ua.traitfile, ua.speciesfile, ua.outfile)

if __name__=="__main__":
    main()