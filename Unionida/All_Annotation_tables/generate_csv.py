import sys,os, subprocess
import argparse 

parser = argparse.ArgumentParser(prog="GenerateCSV")
parser = argparse.ArgumentParser(description='With naming a folder where the annotation tables are, this scripts creates a single csv with all information about each file: species names, order of the genes, their length and lists how much bp the genes distances from eachother')
parser.add_argument('-dir','--Directory',help = "Choose folder with annotation tables")
parser.add_argument('-outdir','--OutDirectory',help = "(optional) Name the directory you wish to output the pretended file")
args = parser.parse_args()

if args.OutDirectory == None:
	args.OutDirectory = os.getcwd()

final_file = open(os.path.join(str(args.OutDirectory),'Annotations.csv'),'w')

for filename in os.listdir(args.Directory):
	current_file = open(os.path.join(str(args.Directory),filename),'r')
	genes = []
	positions = []
	size = []
	cont = 0
	for lines in current_file:
		if '---' not in lines or '' not in lines:
			molecule, type_,gene,start,stop,strand = lines.split('\t')
			genes.append(gene)
			positions.append([start,stop])
			size.append(str(int(positions[cont][1])-int(positions[cont][0])+1))
			cont += 1
	dif = []
	i = 0
	for values in positions:
		if i != len(positions)-1:
			dif.append(str(int(positions[i+1][0])-int(positions[i][1])))
			i += 1
		else:
			continue
	final_file.write(molecule + ',')
	final_file.write(','.join([x for x in genes]))
	final_file.write('\n'+ ' ,'+','.join([x for x in size]))
	final_file.write('\n'+ ' ,'+','.join([x for x in dif]))
	final_file.write('\n\n')
	current_file.close()

print('\n\n..Done!')
print('Output file written to: ',str(args.OutDirectory))
final_file.close())