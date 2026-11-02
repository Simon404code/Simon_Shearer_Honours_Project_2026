import pandas as pd


if __name__ == "__main__":
    combBed_v5_df = pd.read_csv("C:/Users/Simon_Asus_Laptop/Studies/Fourth_year/project/Data/GENESPACE_output/results-20260723T095822Z-1-001/results/combBed_v5_complete-2.txt", sep="\t")
    table2_df = pd.read_csv("C:/Users/Simon_Asus_Laptop/Studies/Fourth_year/project/Data/GENESPACE_output/results-20260723T095822Z-1-001/results/out_table2.5.csv")
    table3_df = pd.read_csv("C:/Users/Simon_Asus_Laptop/Studies/Fourth_year/project/Data/GENESPACE_output/results-20260723T095822Z-1-001/results/out_table3.5.csv")
    #print(combBed_v5_df.head())
    table2_df.set_index('HOG_ID',inplace=True)
    table3_df.set_index('HOG_ID(array)',inplace=True)
    sorted = table3_df.sort_values('max_gene_length', ascending=False)

    #print("table2_df")
    #print(table2_df.head())
    #print(table3_df.head())
    #print(sorted.head())
    
    ##locate specific row in table2_df
    table2_df.loc['N0.HOG0000001 OG0000000']
    output = table2_df.loc['N0.HOG0000001 OG0000000']
    #print(output)

    #print(table2_df.avg_CNV)
    #spefic row on loc
    #print(table2_df.loc['2.0'])
    #print(table2_df['avg_CNV']  > 2.0) #true/ flase
    #print(table2_df.loc[(table2_df['avg_CNV']  > 2.0) & (table2_df['number_of_Genomes'] > 5)])   # use True/False as indexs

    majority = ["core", "soft-core"]
    #print(table2_df.loc[table2_df['classifier'].isin(majority)])  ~ is negating
    #print(table2_df.loc[~table2_df['classifier'].isin(majority)])

    total_tandem_array_count = len(table2_df)
    core_tandem_array_count = len(table2_df.loc[table2_df["classifier"] == "core"])
    softcore_tandem_array_count = len(table2_df.loc[table2_df["classifier"] == "soft-core"])
    core_softcore_tandem_array_count = len(table2_df.loc[table2_df['classifier'].isin(majority)])
    singleton_tandem_array_count = len(table2_df.loc[table2_df["number_of_Genomes"] == 1])
        

    print("Total tandem arrays")
    print(total_tandem_array_count)
    print("core tandem arrays")
    print(core_tandem_array_count)
    print("soft-core tandem arrays")
    print(softcore_tandem_array_count)
    print("core & soft-core tandem arrays")
    print(core_softcore_tandem_array_count)
    print("Singleton tandem array count")
    print(singleton_tandem_array_count)

    # array with each value, can pass to R for pie chart.
    tandem_array_distribution = [total_tandem_array_count, core_tandem_array_count, softcore_tandem_array_count, core_softcore_tandem_array_count]


