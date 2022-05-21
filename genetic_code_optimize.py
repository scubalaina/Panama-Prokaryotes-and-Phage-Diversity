import os, re, sys, shlex, subprocess, argparse
from Bio import SeqIO

args_parser = argparse.ArgumentParser(description="Script for predicting genes with prodigal using the code with optimal density", epilog="Virginia Tech Department of Biological Sciences")
args_parser.add_argument('-i','--inputdir', required=True, help='directory with sequences as individual nucleotide files')
args_parser.add_argument('-o','--outputdir', required=True, help='directory that gene predictions will be deposited')
args_parser.add_argument('-f','--file_tbl', required=True, help='file with table of each genome and its optimal genetic code, coding density, and genome length')

args_parser = args_parser.parse_args()

inputdir = args_parser.inputdir
outputdir = args_parser.outputdir
file_tbl = args_parser.file_tbl


codes = ['1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23', '24', '25']

if inputdir.endswith("/"):
	pass
else:
	inputdir = inputdir + "/"

if outputdir.endswith("/"):
	pass
else:
	outputdir = outputdir + "/"

if os.path.isdir(outputdir):
	pass
else:
	os.mkdir(outputdir)

if os.path.isfile(file_tbl):
	os.remove(file_tbl)
	outfile = open(file_tbl,'a')
else:
	outfile = open(file_tbl,'a')



headerlist = ["genome","code","density","genome_length"]
headers = "\t".join(headerlist)
outfile.write(headers + "\n")

keep_list = []
for i in os.listdir(inputdir):
	if i.endswith(".fna") or i.endswith(".fasta"):
		fasta = os.path.join(inputdir, i)
		nucl_length = float(0)
		genome = re.sub(".fna","",i)

		for j in SeqIO.parse(fasta, "fasta"):
			nucl_length += len(j.seq)

		code2density = {}

		for code in codes:

			#protein = os.path.join(outputdir, re.sub(".fna", "."+str(code)+".faa", i))
			protein = outputdir + genome + "." + code + ".prot.faa"
			orfs = outputdir + genome + "." + code + ".orfs.fna"
			#print(protein)

			if nucl_length < 20000:
				cmd = "prodigal -p meta -i "+ fasta +" -a "+ protein +" -g "+ code + " -d " + orfs

			else:
				cmd = "prodigal -i "+ fasta +" -a "+ protein +" -g "+ code  + " -d " + orfs

			cmd2 = shlex.split(cmd)
			#print(cmd)
			subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))

			prot_length = float(0)
			for j in SeqIO.parse(protein, "fasta"):
				prot_length += len(j.seq)

			#print(prot_length)
			ratio = (prot_length * 3) / nucl_length
			#sys.stdout.write(i +"\t"+ str(code) +"\t"+ str(ratio) +"\t"+ str(nucl_length) +"\n")
			code2density[code] = float(ratio)

		standard = float(code2density['11'])
		#print(standard)
		if standard < float(0.8):
			max_key = max(code2density, key=code2density. get)
			print(genome +"\t"+ str(max_key) +"\t"+ str(code2density[max_key]) +"\t"+ str(nucl_length))
			infolist = [genome,str(max_key),str(code2density[max_key]),str(nucl_length)]
			info = "\t".join(infolist)
			outfile.write(info + "\n")
			new_protein = outputdir + genome + "." + max_key + ".prots.faa"
			new_orfs = outputdir + genome + "." + max_key + ".orfs.fna"
			keep_list.append(new_protein)
			keep_list.append(new_orfs)
			#print(os.path.basename(protein))
			
		else:
			max_key = '11'
			print(genome +"\t"+ str(max_key) +"\t"+ str(code2density[max_key]) +"\t"+ str(nucl_length))
			infolist = [genome,str(max_key),str(code2density[max_key]),str(nucl_length)]
			info = "\t".join(infolist)
			outfile.write(info + "\n")
			new_protein = outputdir + genome + "." + max_key + ".prots.faa"
			new_orfs = outputdir + genome + "." + max_key + ".orfs.fna"
			keep_list.append(new_orfs)
			keep_list.append(new_protein)
			#keep_list.append(os.path.basename(protein))
			#keep_list.append(outprotfile)

		
		for protfiles in os.listdir(outputdir):
			if protfiles.endswith(".faa") or protfiles.endswith("orfs.fna"):
				outprot = outputdir + protfiles
				if outprot in keep_list:
					pass
					a = "b"
				else:
					os.remove(outprot)
					c = "d"

