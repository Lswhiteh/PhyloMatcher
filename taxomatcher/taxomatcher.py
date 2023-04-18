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

    return ap.parse_args()

if __name__=="__main__":
    main()