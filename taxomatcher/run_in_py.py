from gbif_matcher import main as gbif_main

# Call the gbif_main function to fetch GBIF information for the input CSV file
gbif_main(input_csv = '../data/taxolist_first_7.csv',
          threads = 4
           )