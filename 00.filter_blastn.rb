
=begin
out1=File.new("blastn_filtered8080.txt","w")
out2=File.new("list_temp.txt","w")

aa=File.open("blastn_allcontigs_resfinder20240402.txt").each_line do |line|
line.chomp!
#E01_k141_488314_length_4216     lsa(E)_1_JX560992       97.374  1485    39      0       1767    3251    1485    1       0.0     2527    1485
col=line.split("\t")
if col[3].to_f*100/col[12].to_f >= 80.0
puts line
out1.puts "#{col[0].split("\_")[0]}\t#{line}"
out2.puts col[0]
end
end
aa.close

out1.close
out2.close
=end
`cat list_temp.txt | sort | uniq > list_ARGcontigs.txt`
`seqtk subseq all_contigs.fna list_ARGcontigs.txt > ARG_contigs.fna`
#`rm list_temp.txt`
