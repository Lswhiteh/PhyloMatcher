import argparse
import sys 

def main():    
    ap = argparse.ArgumentParser(description="Taxomatcher")
    ap.add_argument(dest="mode", choices=["gbif", "ncbi"], help="Which database to query.")
    ap.add_argument(
        "-c",
        "--csv",
        dest="input_csv",
        required=True,
        help="CSV where first column is a list of target species names to look up.",
    )
    ap.add_argument("-t", "--threads", dest="threads", required=False, default=4)
    ap.add_argument("-e", "--email", dest="email", required="ncbi" in sys.argv)

    ua = ap.parse_args()
    
    if ua.mode == "gbif":
        from . import gbif_matcher      
        gbif_matcher.main(ua.input_csv, ua.threads)
    elif ua.mode == "ncbi":
        from . import entrez_matcher
        entrez_matcher.main(ua)


if __name__=="__main__":
    main()