

library(tidyverse)
str(out_table2)




##box_plot of number of tandem arrays per HOG_ID

ggplot(data = out_table2.5,
       aes(
         x= classifier,
         y= tandem_array_rate_percent
       ))+
  geom_boxplot()
## explanations per graph

## flow diagram -> all genes -> tandem array genes
## size of tandem arrays vs core and per 
## per functional grouping 
## all genes, all duplicates, all tandems (HOG)
## devided tandems by classifiaction?
## bar graph/ pie-chart with number per classification
## rate of tandem genes <50% wh
## Histogram. number of tandem arrays per individual. vairability of tandem array count.
## Histogram number of tandem CNV. per CNV eg 2,3 ..47. number of tandems
## box and whisker sd stable copy number
## CNV range per tandem array vs genome number. exclude stable(0 variation)
## in Methods how the genes were filtered.
## pangenome era. data, tandem data across individuals.

##violin of number of tandem arrays per HOG_ID

ggplot(data = out_table2.5,
       aes(
         x = classifier,
         y = tandem_array_rate_percent
       ))+
  geom_violin()
## remove singleton



## dot plot of genomes vs rate of tandems (number of tandems/number of single copies)

ggplot(data = out_table2.5,
       aes(
         x= number_of_Genomes,
         y = tandem_array_rate_percent
         
       ))+
  geom_point(size = 3,
             alpha = 0.05)




## dot plot of genomes vs rate of tandems (number of tandems/number of single copies)

ggplot(data = out_table2.5,
       aes(
         x= number_of_Genomes,
         y = tandem_array_rate_percent
         
       ))+
  geom_point(size = 3,
             alpha = 0.05)





##dot plot genomes vs genes

ggplot(data = out_table2.5,
       aes(
         x = number_of_Genomes,
         y = number_of_Genes
       ))+
  geom_point(size = 3,
             alpha = 0.005) 



## avg_CNV vs genomes

ggplot(data =  out_table2.5,
       aes(x = number_of_Genomes,
           y = avg_CNV)) +
  geom_point(
    size = 3,
    alpha = 0.005
  )

## avg_CNV vs genomes

ggplot(data =  out_table2.5,
       aes(x = number_of_Genomes,
           y = avg_CNV)) +
  geom_point(
    size = 3,
    alpha = 0.5
  )
## take a few examples and make bar-graphs of how the CNV is shown across
## the 50 genomes. extreme examples. 
## per array show the CNV for each individual. For each tandem array. 
## HOG annotations. from gffs? Give annotation of extreme examples
## add graph that shows the 
## graph CNV per array, add avg line
## every figure should answer a question.
## remember the story you are trying to tell.
## take provenince ID an array where the CNV is high between province
## assign each genome with provenice. Region.
## how many arrays only occour in 1 haplotype/region/provenince 
## differ in sequence IDentity. 
## think about including ID plotting
## chromosome map R2R3-MYB blast ID. 

## frequency vs gene length

ggplot(data = out_table2.5,
       aes(
         x = number_of_Genomes/50,
         y = avg_gene_length
       ))+
  geom_point(
    size = 3,
    alpha = 0.005
  )



