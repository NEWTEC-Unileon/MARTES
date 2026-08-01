hkraken={}
hkrakeng={}
hgenomadv={}
hgenomadp={}
hplasmidfinder={}
hplaton={}	
hintegron={}
hmefinder={}
hfam={}

tt=File.open("gene_list.txt").each_line do |line|
line.chomp!
#Beta-lactam     blaOXA-161_1_GQ202693
col=line.split("\t")
hfam[col[1].split("\_")[0]]=col[0]
end
tt.close

# to try to add all the taxonomy levels... but it is complicate and we normally never use it (UNFINISHED CODE)
=begin
aa=File.open("ARGcontigs_kraken2_PlusPF_20251015_report.txt").each_line do |line|
line.chomp!
if line =~ /(\w\_\_\w+\s+\d+)/ or  line =~ /(\w\_\_\w+\s+\w+\s+\d+)/ or  line =~ /(\w\_\_\w+\s+sp\..*\s+\d+)/
puts $1
end
end
aa.close
=end

genera=[]
jj=File.open("genera_names.txt").each_line do |line|
line.chomp!
genera<<line
end
jj.close


aa=File.open("contigs_filtered.tsv").each_line do |line|
line.chomp!
#C       E01_k141_404296_length_1235     Bacillota (taxid 1239)  1235
#C       E01_k141_250372_length_3906     Enterococcus cecorum (taxid 44008)      3906
#C       E01_k141_27587_length_12420     Chryseobacterium sp. POL2 (taxid 2713414)       12420
col=line.split("\t")
hkraken[col[1]]=col[2].split("\s\(")[0]
	if genera.include?(col[2].split("\s")[0])
	hkrakeng[col[1]]=col[2].split("\s")[0]
	else hkrakeng[col[1]]="unclassified"
	end
end
aa.close
=begin
bb=File.open("blastn_filtered8080.txt").each_line do |line|
line.chomp!
#E01     E01_k141_488314_length_4216     lsa(E)_1_JX560992       97.374  1485    39      0       1767    3251    1485    1       0.0     2527    1485
col=line.split("\t")
samples << col[0]
contigs << col[1]
#hgenome[$2]=$1
genes << col[2].split("_")[0]
#icontig << "#{$1}_#{$2}"
end
bb.close
=end

#icontig_tmp=''
cc=File.open("00.platon/AMR_contigs/AMR_contigs.tsv").each_line do |line|
line.chomp!
#ID      Length  Coverage        # ORFs  RDS     Circular        Inc Type(s)     # Replication   # Mobilization  # OriT  # Conjugation   # AMRs  # rRNAs # Plasmid Hits
#E36_k141_113748_length_70571    70571   NA      81      -0.9    no      0       0       2       0       0       1       0       0 
col=line.split("\t")
hplaton[col[0]]="plasmid detected\t#{col[5]}"	#col[5] --> Circular?
end
cc.close

=begin
`grep -P "\t1\t" 03.integron_finder/Results_Integron_Finder_AMR_contigs/AMR_contigs.summary > integron_finder_detected.txt`
dd=File.open("integron_finder_detected.txt").each_line do |line|
line.chomp!
col=line.split("\t")
hintegron[col[0]]="integron found"
end
1dd.close
=end

n=0
dd=File.open("03.integron_finder/Results_Integron_Finder_AMR_contigs/AMR_contigs.integrons").each_line do |line|
  line.chomp!
  n+=1
  cols=line.split("\t")
  if n!=1
    contig = cols[1]        # ID_replicon
    start  = cols[3].to_i   # pos_beg
    stop   = cols[4].to_i   # pos_end
    strand  = cols[5].to_i  # strand
      if hintegron[contig].nil?
      hintegron[contig] = ["integron found", start, stop, strand]
      else
      # actualizamos rango mínimo y máximo
      hintegron[contig][1] = [hintegron[contig][1], start].min
      hintegron[contig][2] = [hintegron[contig][2], stop].max
      end
#puts hintegron[contig][3]
   end
end
dd.close


ee=File.open("01.plasmidfinder/results_tab.tsv").each_line do |line|
line.chomp!
#Database        Plasmid Identity        Query / Template length Contig  Position in contig      Note    Accession number
#enterobacteriaceae      Col3M   98.09   157 / 157       E76_k141_355649_length_3008     2643..2799              JX514065
col=line.split("\t")
hplasmidfinder[col[4]]="plasmid detected\t#{col[1]}"  #plasmid identity
end
ee.close

ff=File.open("04.genomad_results/AMR_contigs_summary/AMR_contigs_plasmid_summary.tsv").each_line do |line|
line.chomp!
col=line.split("\t")
hgenomadp[col[0]]="plasmid detected\t#{col[2]}\t#{col[9]}"
end
ff.close

gg=File.open("04.genomad_results/AMR_contigs_summary/AMR_contigs_virus_summary.tsv").each_line do |line|
line.chomp!
col=line.split("\t")
hgenomadv[col[0]]="virus detected\t#{col[2]}\t#{col[10]}"
end
gg.close

hh=File.open("02.mefinder/03.mefinder.csv").each_line do |line|
line.chomp!
col=line.split("\,")
hmefinder[col[12]]="virus detected\t#{col[1]}\t#{col[2]}\t#{col[4]}\t#{col[13]}\t#{col[14]}"
end
hh.close



out=File.new("megatabla_ARGcontigs.txt","w")
out2=File.new("temp.txt","w")
out3=File.new("temp2.txt","w")

#out.puts "Sample\tContig\tARG\tintegron\tPlasmid_platon\tLength\tCoverage\tORFs\tRDS\tCircular\tInc_Type(s)\tReplication\tMobilization\tOriT\tConjugation\tAMRs\trRNAs\tPlasmid Hits"
out.puts "Sample\tContig\tARG\tstart\tend\tstrand\tAntibiotic_family\tkraken_taxonomy\tkraken_genus\tintegron\tintegron_start\tintegron_end\tintegron_strand\tintegron_in\tplaton_Plasmid\tplaton_Circular\tplasmidfinder_Plasmid\tplasmidfinder_identity\tgenomad_Plasmid\tgenomad\tgenomad_circular\tgenomad_conjugation_gene\tgenomad_Virus\tgenomad_Virus_topology\tgenomad_Virus_taxonomy\tmefinder_mge\tmefinder_name\tmefinder_synonym\tmefinder_type\tmefinder_start\tmefinder_end\tmefinder_in"

bb=File.open("AMR_contigs_blastn_fixed.txt").each_line do |line|
line.chomp!
#E01     E01_k141_488314_length_4216     lsa(E)_1_JX560992       97.374  1485    39      0       1767    3251    1485    1       0.0     2527    1485
col=line.split("\t")
	if hkraken[col[1]]==nil
        hkraken[col[1]]= "-"
        end
	if hkrakeng[col[1]]==nil
        hkrakeng[col[1]]= "-"
        end
	if hintegron[col[1]]==nil
        #hintegron[col[1]]= "not found\t-\t-\t-"
	hintegron[col[1]] = ["not found", "-", "-", "-"]
	end
	if hplaton[col[1]]==nil
	#hplaton[icontig_temp]="no detected\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-"
	hplaton[col[1]]="no detected\t-"
	end
 	if hplasmidfinder[col[1]]==nil
        hplasmidfinder[col[1]]= "no detected\t-"
        end
	if hgenomadp[col[1]]==nil
        hgenomadp[col[1]]= "no detected\t-\t-"
        end
        if hgenomadv[col[1]]==nil
        hgenomadv[col[1]]= "no detected\t-\t-"
        end
	if hmefinder[col[1]]==nil
        hmefinder[col[1]]= "no detected\t-\t-\t-\t-\t-"
        end
  if col[9].to_i < col[10].to_i
  out2.puts "#{col[0]}\t#{col[1]}\t#{col[2].split("_")[0]}\t#{col[7]}\t#{col[8]}\t1\t#{hfam[col[2].split("_")[0]]}\t#{hkraken[col[1]]}\t#{hkrakeng[col[1]]}\t#{hintegron[col[1]].join("\t")}\t#{hplaton[col[1]]}\t#{hplasmidfinder[col[1]]}\t#{hgenomadp[col[1]]}\t#{hgenomadv[col[1]]}\t#{hmefinder[col[1]]}"
  else
  out2.puts "#{col[0]}\t#{col[1]}\t#{col[2].split("_")[0]}\t#{col[7]}\t#{col[8]}\t0\t#{hfam[col[2].split("_")[0]]}\t#{hkraken[col[1]]}\t#{hkrakeng[col[1]]}\t#{hintegron[col[1]].join("\t")}\t#{hplaton[col[1]]}\t#{hplasmidfinder[col[1]]}\t#{hgenomadp[col[1]]}\t#{hgenomadv[col[1]]}\t#{hmefinder[col[1]]}"
  end


#out.puts "#{col[0]}\t#{col[1]}\t#{col[2].split("_")[0]}\t#{hfam[col[2].split("_")[0]]}\t#{hkraken[col[1]]}\t#{hkrakeng[col[1]]}\t#{hintegron[col[1]]}\t#{hplaton[col[1]]}\t#{hplasmidfinder[col[1]]}\t#{hgenomadp[col[1]]}\t#{hgenomadv[col[1]]}\t#{hmefinder[col[1]]}"
#out.puts "#{col[0]}\t#{col[1]}\t#{col[2].split("_")[0]}\t#{hfam[col[2].split("_")[0]]}\t#{hkraken[col[0]]}\t#{hkrakeng[col[0]]}\t#{hintegron[col[0]]}\t#{hplaton[col[0]]}\t#{hplasmidfinder[col[0]]}\t#{hgenomadp[col[0]]}\t#{hgenomadv[col[0]]}\t#{hmefinder[col[0]]}"
end
bb.close

out2.close

cc=File.open("temp.txt").each_line do |line|
line.chomp!
col=line.split("\t")
#puts col[1..4].join("\t")
if col[9]=="integron found"
	if col[3].to_i>col[10].to_i-1 and  col[4].to_i<col[11].to_i+1
	puts col[9]
	puts "#{col[3]}\t#{col[4]}\t#{col[10]}\t#{col[11]}\t"
	out3.puts "#{col[0..12].join("\t")}\tin\t#{col[13..28].join("\t")}"
	else out3.puts "#{col[0..12].join("\t")}\tout\t#{col[13..28].join("\t")}"
	end
else out3.puts "#{col[0..12].join("\t")}\t-\t#{col[13..28].join("\t")}" 
end 
end 
cc.close
out3.close

dd=File.open("temp2.txt").each_line do |line|
line.chomp!
col=line.split("\t")
#puts col[1..4].join("\t")
if col[27] != "-"
        if col[3].to_i>col[28].to_i-1 and  col[4].to_i<col[29].to_i+1
        puts col[27]
        puts "#{col[3]}\t#{col[4]}\t#{col[28]}\t#{col[29]}\t"
        out.puts "#{line}\tin"
        else out.puts "#{line}\tout"
        end
else out.puts "#{line}\t-"
end
end
dd.close


