import sys, os, subprocess
import argparse
from ete3 import EvolTree

def parse_args():
	"""Parsing arguments via terminal"""
	parser = argparse.ArgumentParser(prog="selection_analysis")
	parser = argparse.ArgumentParser(description='With an imput of MSA and a phylogeny, plus an additional txt detailing the phylogenetic tree structure, this script runs models of codeml (whole, site, site-branch), can do LRT tests, can download images, and can process the output folders of what you want')
	parser.add_argument('-tree', '--Tree', help='Give this command the phylogenetic tree in newick format')
	parser.add_argument('-msa', '--MSA', help = 'Give this command a MSA file that produced the phylogeny')
	parser.add_argument('-tree_structure', '--TreeStruct', help = 'Give this command a phylogenetic tree structure file (example is in the examples directory), this is only necessary when running any models that use branches')
	parser.add_argument('-binpath', "--BinPath", help = '(Optional) For some reason, sometimes codeml is installed in another unexpected place, if getting an error message that there is no codeml in a directory, find it (probably in ete3_apps) and input the full directory here')
	parser.add_argument('-models', "--Models", help = 'You can choose any number of models, however be aware that you have to specify with - the models that will be compared')
	parser.add_argument("-compare", "--Compare", help = "Choose the models you want to compare, seperating them by commas first and then by -")
	parser.add_argument("-load", "--LoadedModels", help = "In case you already ran the models you want to analyse, you can load them with this command")
	args = parser.parse_args()
	return args

def main(args):
	if args.BinPath:
		tree = EvolTree(args.Tree, binpath = args.BinPath)
	else: 
		tree = EvolTree(args.Tree)
	if args.MSA[:-3] == ".phy":
		tree.link_to_alignment(args.MSA, format = "phylip")
	elif args.MSA: 
		tree.link_to_alignment(args.MSA)
	print(tree)
	tree.workdir = os.getcwd()
	if args.LoadedModels:
		load_model(args.LoadedModels, tree)
		compare_models(models = args.LoadedModels, tree = tree, args = args)
	if args.Models:
		run_models(args.models, tree)
	if args.Compare:
		if args.TreeStruct:
			tree_structure = parse_structure_file(args.TreeStruct)
			compare_models(models = args.Compare, tree = tree, tree_structure = tree_structure)
		else: compare_models(models = args.Compare, tree = tree, args = args)
	
def run_models(models, tree):
	models = models.split(",")
	if len(models) == 1:
		model_1 = models
		print("\nNow running model {}\n" .format(model_1))
		tree.run_model(model_1)
	elif len(models) == 2:
		model_1 = models[0]
		print("\nNow running model {}\n" .format(model_1))
		tree.run_model(model_1)
		model_2 = models[1]
		print("\nNow running model {}\n" .format(model_2))
		tree.run_model(model_2)
	elif len(models) == 3:
		model_1 = models[0]
		print("\nNow running model {}\n" .format(model_1))
		tree.run_model(model_1)
		model_2 = models[1]
		print("\nNow running model {}\n" .format(model_2))
		tree.run_model(model_2)
		model_3 = models[2]
		print("\nNow running model {}\n" .format(model_3))
		tree.run_model(model_3)
	elif len(models) == 4:
		model_1 = models[0]
		print("\nNow running model {}\n" .format(model_1))
		tree.run_model(model_1)
		model_2 = models[1]
		print("\nNow running model {}\n" .format(model_2))
		tree.run_model(model_2)
		model_3 = models[2]
		print("\nNow running model {}\n" .format(model_3))
		tree.run_model(model_3)
		model_4 = models[3]
		print("\nNow running model {}\n" .format(model_4))
		tree.run_model(model_4)
	else:
		sys.exit("You gave none or more than 4 models, exiting.")
	
def compare_models(models = "", tree = "", tree_structure = "", args = ""):
	model_A, model_B = models.split(",")
	if args.LoadedModels:
		print('Comparing models..')
		pvalue = tree.get_most_likely(model_B, model_A)
		print(str(pvalue))
		test_significance(model_B, pvalue, tree)
	elif model_A.startswith("b") or model_B.startswith("b"):
		for branch in tree_structure:
			name = branch[2]
			print("\nWorking on branch: " + name)
			mark_phylogeny(branch, tree)	
			print("\nNow running model {}\n" .format(model_A))
			tree.run_model(model_A + '.' + name[:-1])
			print("\nNow running model {}\n" .format(model_B))
			tree.run_model(model_B + '.' + name[:-1])
			print('Comparing models..')
			pvalue = tree.get_most_likely(model_B + '.' + name[:-1] , model_A + '.' + name[:-1])
			print(str(pvalue))
			unmark_phylogeny(tree)
			test_significance(model_B + '.' + name[:-1], pvalue, tree)

def test_significance(selected_model, pvalue, tree):
	model = tree.get_evol_model(selected_model)
	print(model)
	print(tree.write())
	if pvalue <= 0.05 and float(model.classes["foreground w"][2]) > 1:
		print('we have positive selection on sites on this branch')
		#tree.show(histfaces=[model])

def mark_phylogeny(branch, tree):
	# tenho de por apenas um #1 em cada branch, limpar e repetir para todos, bota la programar!!!
	ancestor = tree.get_common_ancestor(branch[0], branch[1])
	for leaf in ancestor:
		tree.mark_tree([leaf.node_id], marks = ["#1"])

def unmark_phylogeny(tree):
	tree.mark_tree(map(lambda x: x.node_id, tree.get_descendants()),
						marks = [""] * len(tree.get_descendants()))

def parse_structure_file(file):
	tree_structure = []
	with open(file, "r") as f:
		for line in f:
			leaf_1, leaf_2, name = line.split(";")
			tree_structure.append((leaf_1, leaf_2, name))
	return tree_structure

def load_model(model_dir, tree):
	if len(model_dir.split(",")) == 2:
		modeldir_1, modeldir_2 = model_dir.split(",")
		tree.link_to_evol_model(modeldir_1 + "/out", modeldir_1) # ver o segundo argumento
		tree.link_to_evol_model(modeldir_2 + "/out", modeldir_2)
	elif len(model_dir.split(",")) == 3:
		modeldir_1, modeldir_2, modeldir_3 = model_dir.split(",")
		tree.link_to_evol_model(modeldir_1 + "/out", modeldir_1)
		tree.link_to_evol_model(modeldir_2 + "/out", modeldir_2)
		tree.link_to_evol_model(modeldir_3 + "/out", modeldir_3)
	elif len(model_dir.split(",")) == 4:
		modeldir_1, modeldir_2, modeldir_3, modeldir_4 = model_dir.split(",")
		tree.link_to_evol_model(modeldir_1 + "/out", modeldir_1)
		tree.link_to_evol_model(modeldir_2 + "/out", modeldir_2)
		tree.link_to_evol_model(modeldir_3 + "/out", modeldir_3)
		tree.link_to_evol_model(modeldir_4 + "/out", modeldir_4)


'''
tree = EvolTree("tree.nw", binpath = "/home/edu/miniconda3/envs/ete3/bin/ete3_apps/bin")
tree.link_to_alignment("infile.phy", alg_format = "phylip")
tree.workdir = os.getcwd()

print(tree)

print('running model M0, for comparison with branch-site models...')

tree.run_model('M0', keep = True)
#tree.link_to_evol_model("/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Phylogenetic_Tree","M0")
chimaeriformes =tree.get_common_ancestor("HM147138.1","HM147135.1")
#chimaeriformes =tree.get_common_ancestor("Human_ECP","Goril_ECP")

for leaf in chimaeriformes:
    tree.mark_tree([leaf.node_id], marks = ["#1"])
#tree.run_model("bsA." + chimaeriformes)
#tree.mark_tree([leaf.node_id], marks = ["#1"])
print("Running")
print(tree.write())
tree.run_model('bsA.Chimaeriformes')
tree.run_model("bsA1.Chimaeriformes")

print ('p-value of positive selection for sites on this branch is: ')
ps = tree.get_most_likely('bsA.Chimaeriformes' , 'bsA1.Chimaeriformes')
print(str(ps))
rx = tree.get_most_likely('bsA1.Chimaeriformes', 'M0')
print(str(rx))
model = tree.get_evol_model("bsA.Chimaeriformes")
if ps < 0.05 and float(model.classes['foreground w'][2]) > 1:
    print ('we have positive selection on sites on this branch')
    tree.show(histfaces=['bsA1.Chimaeriformes'])
elif rx<0.05 and ps>=0.05:
    print ('we have relaxation on sites on this branch')
else:
    print ('no signal detected on this branch, best fit for M0')
#tree.show(histfaces=['bsA1.'])

for models in tree._models:
    print(tree.get_evol_model(models))

from _pickle import dump

#out = open('my_tree.pik', 'w')
#dump(tree, out)
#out.close()

#tree.mark_tree(map(lambda x: x.node_id, tree.get_descendants()),
#                    marks=[''] * len(tree.get_descendants()), verbose=True)

print("The End.")
'''

if __name__ == "__main__":
	args = parse_args()
	main(args)
