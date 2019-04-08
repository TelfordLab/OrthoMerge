import argparse, math, os, sys
from Bio import SeqIO

work_dir = os.getcwd()+"/"

parser = argparse.ArgumentParser(description='Merges the orthology predictions of OrthoMCL, Orthofinder and Orthoinspector')
parser.add_argument("-p", "--phylum", required=True, help="phylum of interest")
parser.add_argument("--phyla", default="phyla.txt", help="path to file defining phyla and superphyla (default: ./phyla.txt)")
parser.add_argument("--orthomcl", required=True, help="path to OrthoMCL groups file")
parser.add_argument("--orthofinder", required=True, help="path to OrthoFinder groups file")
parser.add_argument("--orthoinspector", required=True, help="path to OrthoInspector groups file")
parser.add_argument("-fd","--fastadir", help="path to sequence file directory (REQUIRED unless --no-groups)")
parser.add_argument("--only-identical", action="store_true", help="only produce FASTA files for orthologous groups identical in all methods")
parser.add_argument("--include-ancestral", action="store_true", help="produces additional FASTA files for all ancestral orthologous groups")
parser.add_argument("-o","--output", default=work_dir+"OrthologousGroups_identical+split/", help="output folder for orthologous groups (default: ./OrthologousGroups_identical+split)")
parser.add_argument("--no-output", action="store_true", help="do not produce FASTA files for merged orthologous groups")
args = parser.parse_args()

if not args.no_output and not args.fastadir:
    print(parser.format_help())
    print("Option '--fastadir'/'-fd' needed for producing output!")
    sys.exit(0)

if not args.no_output and args.output and os.path.exists(args.output):
    print("ABORT: directory {} already exists!".format(args.output))
    print("Specify a different directory with '--output'/'-o'")
    sys.exit(0)

groups_dir = args.output + ("" if args.output.endswith("/") else "/")
log_file = groups_dir[:-1] + ".log"

# read phyla information
all_species = set()
phyla = dict()
with open(args.phyla, "r") as f:
    for line in f:
        if len(line.strip(" \n")) < 1:
            continue # ignore empty lines
        phylum = line.strip(" \n").split(":")[0]
        if "|" in line: # definition of superphylum
            phyla[phylum] = []
            subphyla = [{p} for p in line.strip(" \n").split(":")[1].split("|")]
            for subphylum in subphyla:
                unresolved = True
                while unresolved:
                    unresolved = False
                    add = set()
                    remove = set()
                    for species in subphylum:
                        if species in phyla:
                            unresolved = True
                            remove.add(species)
                            if isinstance(phyla[species],list): # 2nd level phylum
                                add = {s for sphylum in phyla[species] for s in sphylum} # flatten subphylum
                            else:
                                for s in phyla[species]:
                                    add.add(s)
                    subphylum.update(add)
                    subphylum = subphylum - remove
                phyla[phylum].append(subphylum)
        else: # definition of phylum
            phyla[phylum] = set(line.strip(" \n").split(":")[1].split(","))
            all_species.update(phyla[phylum])

# create set of non-phylum organisms
phylum_of_interest = args.phylum
phyla["non-"+phylum_of_interest] = all_species - {s for sphylum in phyla[phylum_of_interest] for s in sphylum}

# set up dictionaries for group results
groups = dict()
sequences2groups = dict()

groups["OrthoMCL"] = dict()
groups["OrthoMCL"]["file"] = args.orthomcl

groups["OrthoFinder"] = dict()
groups["OrthoFinder"]["file"] = args.orthofinder

groups["OrthoInspector"] = dict()
groups["OrthoInspector"]["file"] = args.orthoinspector

for method in groups:
    sequences2groups[method] = dict()

for method in groups:
    
    groups[method]["groups"] = dict()
    groups[method]["specific"] = set()
    groups[method]["clade-specific"] = set()
    groups[method]["ancestral"] = set()
    groups[method]["nonspecific"] = set()
    
    # read all inferred orthogroups
    with open(groups[method]["file"], "r") as f:
        for line in f:
            if len(line) > 2 and len(line.replace("\n", "").split(" ")) > 2:
                groups[method]["groups"][line.split(": ")[0]] = line.strip("\n").split(" ")[1:]
                for sequence in line.strip("\n").split(" ")[1:]:
                    sequences2groups[method][sequence] = line.split(": ")[0]
    
    # identify clade-specific, specific, ancestral and non-specific groups wrt chosen phylum
    for group_name in groups[method]["groups"]:
        specific = True
        for gene in groups[method]["groups"][group_name]:
            if gene.split("|")[0] in phyla["non-"+phylum_of_interest]:
                specific = False
                break
        if specific:
            clade_specific = True
            represented = []
            for i in range(0, len(phyla[phylum_of_interest])):
                represented.append(False)
            for gene in groups[method]["groups"][group_name]:
                for i in range(0, len(phyla[phylum_of_interest])):
                    if gene.split("|")[0] in phyla[phylum_of_interest][i]:
                        represented[i] = True
            clade_specific = True
            for i in range(0, len(phyla[phylum_of_interest])):
                clade_specific = clade_specific and represented[i]
            
            if clade_specific:
                groups[method]["clade-specific"].add(frozenset(groups[method]["groups"][group_name]))
            groups[method]["specific"].add(frozenset(groups[method]["groups"][group_name]))
        else:
            represented = False
            for gene in groups[method]["groups"][group_name]:
                for i in range(0, len(phyla[phylum_of_interest])):
                    if gene.split("|")[0] in phyla[phylum_of_interest][i]:
                        represented = True
                        break
                if represented:
                    break
            if represented:
                groups[method]["ancestral"].add(frozenset(groups[method]["groups"][group_name]))
            groups[method]["nonspecific"].add(frozenset(groups[method]["groups"][group_name]))
            

# identify groups which have been identified identically and groups that are in split between methods
identical_groups = dict()
split_groups = dict()
for set_name in ["clade-specific", "ancestral"]:
    identical_groups[set_name] = set()
    first_method = True
    for method in groups:
        if first_method:
            identical_groups[set_name] = groups[method][set_name]
            first_method = False
        else:
            identical_groups[set_name] = identical_groups[set_name] & groups[method][set_name]
    
    split_groups[set_name] = identical_groups[set_name] | set()
    
    for method1 in groups:
        for group in groups[method1][set_name]:
            if group in split_groups[set_name]:
                continue
            split_group = True
            method_groups = dict()
            method_singletons = dict()
            for method2 in sequences2groups:
                method_groups[method2] = set()
                method_singletons[method2] = set()
                for sequence in group:
                    if sequence not in sequences2groups[method2]:
                        method_singletons[method2].add(sequence)
                    else:
                        group_name = sequences2groups[method2][sequence]
                        for seq in groups[method2]["groups"][group_name]:
                            method_groups[method2].add(seq)
            for method2 in method_groups:
                if method_groups[method2] != method_groups[method1] - method_singletons[method2]:
                    split_group = False
                    break
            
            if split_group:
                split_groups[set_name].add(group)

# write numerical information about inferred orthogroups
out = ""
for method in groups:
    out += method + "_total:"
    for group in groups[method]["groups"]:
        if len(groups[method]["groups"][group]) > 1:
            out += " " + str(len(groups[method]["groups"][group]))
    out += "\n"

for method in groups:
    out += method + "_" + phylum_of_interest + "_specific:"
    for group in groups[method]["clade-specific"]:
        if len(group) > 1:
            out += " " + str(len(group))
    out += "\n"

with open("total.stats", "w") as f:
    f.write(out)

# write log file
out = "Total number of orthologous groups:\n"
for method in groups:
    num_groups = 0
    num_sequences = 0
    for group in groups[method]["groups"]:
        if len(groups[method]["groups"][group]) > 1:
            num_groups += 1
            num_sequences += len(groups[method]["groups"][group])
        
    out += "   {}\t{} (mean: {})\n".format(method, len(groups[method]["groups"]), (num_sequences/num_groups))

out += "Number of orthologous groups not specific to '{}':\n".format(phylum_of_interest)
for method in groups:
    out += "   {}\t{}\n".format(method, len(groups[method]["nonspecific"]))

out += "Number of orthologous groups specific to '{}':\n".format(phylum_of_interest)
for method in groups:
    out += "   {}\t{}\n".format(method, len(groups[method]["specific"]))

out += "Number of orthologous groups present in and specific to the ancestor of '{}' :\n".format(phylum_of_interest)
for method in groups:
    num_groups = 0
    num_sequences = 0
    for group in groups[method]["clade-specific"]:
        if len(group) > 1:
            num_groups += 1
            num_sequences += len(group)
    
    out += "   {}\t{} (mean: {})\n".format(method, len(groups[method]["clade-specific"]), (num_sequences/num_groups))
out += " {} identical groups found by all methods\n".format(len(identical_groups["clade-specific"]))
out += " {} identical or split groups found by all methods\n".format(len(split_groups["clade-specific"]))

out += "Number of orthologous groups present in the ancestor of '{}', but not specific :\n".format(phylum_of_interest)
for method in groups:
    out += "   {}\t{}\n".format(method, len(groups[method]["ancestral"]))
out += " {} identical groups found by all methods\n".format(len(identical_groups["ancestral"]))
out += " {} identical or split groups found by all methods\n".format(len(split_groups["ancestral"]))

with open(log_file, "w") as f:
    f.write(out)

# exit if identified groups should _not_ be written to FASTA files
if args.no_output:
    sys.exit(0)

# read sequences into dictionary
fasta_dir = args.fastadir + "/"
proteomes = dict()
for fasta_file in os.listdir(fasta_dir):
    if not (fasta_file.endswith(".fa") or fasta_file.endswith(".fasta")):
        print("SKIPPING: {} is not a FASTA file.".format(fasta_file))
        continue
    sequences = list(SeqIO.parse(fasta_dir + fasta_file, "fasta"))
    for sequence in sequences:
        species = sequence.id.split("|")[0]
        if species not in proteomes:
            proteomes[species] = dict()
        proteomes[species].update(SeqIO.to_dict([sequence]))

# export sequences according to merged groups

if args.only_identical:
    if args.include_ancestral:
        orthologous_groups = identical_groups["clade-specific"] | identical_groups["ancestral"]
    else:
        orthologous_groups = identical_groups["clade-specific"]
    if not args.output:
        if args.include_ancestral:
            groups_dir = work_dir + "Ancestral_Orthologous_Groups_identical/"
        else:
            groups_dir = work_dir + "Orthologous_Groups_identical/"
else:
    if args.include_ancestral:
        orthologous_groups = split_groups["clade-specific"] | split_groups["ancestral"]
    else:
        orthologous_groups = split_groups["clade-specific"]
    if not args.output:
        if args.include_ancestral:
            groups_dir = work_dir + "Ancestral_Orthologous_Groups_identical+split/"
        else:
            groups_dir = work_dir + "Orthologous_Groups_identical+split/"

if not os.path.exists(groups_dir):
    os.mkdir(groups_dir)
    print("Created new directory " + groups_dir)

log10 = int(math.log10(len(orthologous_groups)))
counter = 1
for group in orthologous_groups:
        sequences = []
        for sequence in group:
            species = sequence.split("|")[0]
            sequences.append(proteomes[species][sequence])
        
        file_name = groups_dir + "G"+"".join([ "0" for i in range(log10 - int(math.log10(counter)))]) + str(counter) + ".fasta"
        SeqIO.write(sequences, file_name, "fasta")
        counter += 1


