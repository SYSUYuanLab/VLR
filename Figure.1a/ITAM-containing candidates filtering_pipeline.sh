#1.Searching for ITAM motif
less 3iso_all_pep_1line.fasta |egrep "Y\w\w[IL]\w{6,12}Y\w\w[IL]" -B1 >1st_ITAM_pep.fasta
sed -i 's/--//g;/^$/d' 1st_ITAM_pep.fasta


#2.Transmembrane prediction
#2.1 TMHMM 2.0  
https://services.healthtech.dtu.dk/services/TMHMM-2.0/ 
#online Output format :Extensive,no grphics

#2.2 Phobius
https://phobius.sbc.su.se/ 
#online Output format: Short

#2.3 Selecting protein sequences with TM values of 1 or 2 for downstream analysis
awk -F' ' 'NF==4' Phobius_prediction.html |awk '$2=="1" || $2=="2"'|awk '{print $1}' >ID1  #Phobius
less TMHMM_prediction.html |egrep "Number of predicted TMHs:  [12]$" |awk '{print $2}' >ID2  #TMHMM 2.1

cat ID1 ID2 >ID
sort ID|uniq >ID.txt
mv ID.txt 2nd_Phobiu+TMHMM_1or2TM_ID.txt
grep -F -f 2nd_Phobius+TMHMM_1or2TM_ID.txt 3iso_all_pep_1line.fasta -A1 > 2nd_Phobius+TMHMM_1or2TM_pep.fasta


#3.Signal peptide prediction
#SignalP4/5/6 input:2nd_Phobius+TMHMM_1or2TM_pep.fasta
#downloading processed mature fasta withought signal peptides。
#Since the SignalP4's Processed.fasta file will add record.description, special processing is required for SignalP4. All - SignalP-noTM and SignalP-TM output (no graphics).
python clean_fasta_id.py -i 3_SignalP4_output_mature.fasta -o 3_SignalP4_output_mature
mv 3_SignalP4_output_mature 3_SignalP4_output_mature.fasta

#recorrecting to the right ID.
sed -i 's/__/::/g' 3rd_SignalP4_processed.fasta
sed -i 's/__/::/g' 3rd_SignalP5_processed.fasta
sed -i 's/__/::/g' 3rd_SignalP6_processed.fasta

sed -i 's/F01_transcript_/F01_transcript\//g' 3rd_SignalP4_processed.fasta
sed -i 's/F01_transcript_/F01_transcript\//g' 3rd_SignalP5_processed.fasta
sed -i 's/F01_transcript_/F01_transcript\//g' 3rd_SignalP6_processed.fasta

#merge and extracting unique values
python merge_dedup_id_seq.py -i 3rd_SignalP4_processed.fasta 3rd_SignalP5_processed.fasta 3rd_SignalP6_processed.fasta -o 3rd_SignalP456_processed.fasta

#4.Predicting transmembrane again
#4.1 TMHMM 2.0
https://services.healthtech.dtu.dk/services/TMHMM-2.0/ 
#online Output format :Extensive,no grphics

less 4th_TMHMM_prediction.html |egrep "Number of predicted TMHs:  [1]$" |awk '{print $2}' >ID1

#4.2 Phobius
https://phobius.sbc.su.se/ 
#online Output format: Short

awk -F' ' 'NF==4' 4th_Phobius_prediction.html |awk '$2=="1"'|awk '{print $1}' >ID2

#4.3 merge
cat ID1 ID2>ID
sort ID|uniq >ID.txt
mv ID.txt 4th_Phobius+TMHMM+DeepTMHMM_1TM_ID.txt
grep -F -f 4th_Phobius+TMHMM_1TM_ID.txt 3iso_all_pep_1line.fasta -A1 > 4th_Phobius+TMHMM_1TM_pep.fasta

#5.Identify ITAM regions
#Re-predict transmembrane regions using full-length sequences with TMHMM, Phobius (Long without Graphics), and DeepTMHMM.
#Check if the amino terminus is extracellular, the carboxyl terminus is intracellular, and if the intracellular region contains ITAM motifs.

#5.1 TMHMM2
python ITAMoverTM_TMHMM.py -f 4th_Phobius+TMHMM_1TM_pep.fasta -t 5th_87pep_TMHMM\ result.html > 5th_TMHMM.txt

#5.2 DeepTMHMM1
#ITAM
sed -i 's/__/::/g' 5th_pep_TMRs.gff3
sed -i 's/F01_transcript_/F01_transcript\//g' 5th_pep_TMRs.gff3
python ITAMoverTM_DeepTMHMM.py -f 4th_Phobius+TMHMM_1TM_87pep.fasta -t 5th_pep_TMRs.gff3  > 5th_DeepTMHMM.txt

#5.3 Phobius (long,no graphic)
python ITAMoverTM_Phobius.py -f 4th_Phobius+TMHMM_1TM_87pep.fasta -t 5th_87pep_Phobius_prediction.html > 5th_Phobius.txt

#5.4 
python parse_TM_motif.py -i 5th_TMHMM.txt -o 5th_TMHMM_ITAM.txt
python parse_TM_motif.py -i 5th_DeepTMHMM.txt -o 5th_DeepTMHMM_ITAM.txt
python parse_TM_motif.py -i 5th_Phobius.txt -o 5th_Phobius_ITAM.txt
#Summarizing results where the ITAM motif is intracellular."
