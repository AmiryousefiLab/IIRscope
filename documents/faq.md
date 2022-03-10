#### Q1: Why the gene names are not displayed properly?

**A1**: Check your input file(s), IRscope is based on the gene names specified in them.

</br>

#### Q2: Why some of the genes are extended to the large portion of the tracks?

**A2**: This is due to the bad annotation implemented in the input file(s). Often this is related to the genes with introns or trans-spliced genes like rps12. Check and fix your input files.

</br>

#### Q3: I have a well annotated genome but not getting the plot?

**A3**: IIRscope can handle the redundant base pairs but if the sequence is of a very poor quality the program may fail in detecting the inverted regions and hence not produce any output.

</br>


#### Q4: I see genes overplotting on each other, how can I fix this?

**A4**: This is mostly because at least one of your species differs significantly from the rest so that it creates too large radius for finding the genes in the vicinity of the junction site and consequently causes the populated genes in these respective areas. Try reducing the species sampling to a smaller group of more closely related species.

</br>

#### Q5: There is a warning about the spacing around the junction not being in scale. Why?

**A5**: If at least one of your species differs significantly from the rest so that it creates too large radius for finding the genes in the vicinity of the junction site, the radius in that junction will adapt to each of the species best, instead of using the same radius. Try reducing the species sampling to a smaller group of more closely related species.

</br>

#### Q6: Can I modify and run the codes on my own computer

**A6**: Yes, please download the files from [GitHub](https://github.com/AmiryousefiLab/IIRscope), run the whole script in your R session. This will set up a pseudo web interface on your computer which resemles and works exactly like the online version. You can further tune the functions as you wish.

</br>

#### Q7: What should be the format of the annotation files in the 'Manual Files' section?

**A7**: This should be a plain text format (.txt) of the four tab separated columns, as start of the gene, end of the gene, its name, and the direction of it (+ or -). Note that the file needs to be without header.

</br>

#### Q8: My species names are not plotted correctly, why?

**A8**: The genome names in the 'GB Files' are read from the 'Organism' line of the file while the first space separated text of the genomes first line in the manual part is a determinant of species name in this section. In the 'Manual Files' section, the names are read from the first line after the '>' of the fasta file uploaded. Please edit them if they do not follow your expectations and rerun.

</br>

#### Q9: My analysis took excessively long time and even after that the plot is of an overall low quality, why is that?

**A9**: The IIRscope is primarily optimized for the angiosperms and it often turns reliable results for other seed plants as well. However, further departure from the Embryophyta will extend this program to its limits. In such cases, you may want to consider use of the manual section with providing the the IR coordinates.

</br>

#### Q10: How can I cite the program

**A10**: The paper describing this web app in detail with the title *IRscope+: The updated program to visualise the junction sites of chloroplast genomes* is submitted to [** journal](https://doi.org/10.1093/bioinformatics/bty220). Please, cite the paper when you use the program.