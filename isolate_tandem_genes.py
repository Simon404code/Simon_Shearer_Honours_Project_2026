#####
##### This file will take in a .tsv file and isolate tandemly duplicated arrays
##### It will do that by looking for genes that are in a common orthogroup and are located on the 
##### same chromosome.
##### After finding one it will make a tandem array and in it store the following:    chr	start	end	id	ofID	pepLen	ord	genome	arrayID	isArrayRep	globOG  globHOG	noAnchor	og
##### It will split the tandem arrays based on genome
##### And then output a .tsv file for each genome with the tandem arrays in it and one with no tandem arrays in it
##### Each genome inputted will have its own output .tsv files, 1 for tandem arrays and 1 for no tandem arrays.


if __name__ == "__main__":
    print("This script is not meant to be run directly. Please run the main script instead.")

