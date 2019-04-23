import argparse, os, subprocess, concurrent.futures, multiprocessing, sys
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from dendropy import Tree


work_dir = os.getcwd() + "/"

parser = argparse.ArgumentParser(description='Tests orthologous group for potential orthology to closest potential homologues')
parser.add_argument("-od", "--ogroups-dir", required=True, help="folder containing orthologous groups (FASTA format)")
parser.add_argument("-hd", "--hits-dir", required=True, help="folder containing blast hits from monophyly check")
parser.add_argument("-m", "--monophyly", required=True, help="monophyly test results file")
parser.add_argument("-dbd", "--database-dir", required=True, help="directory containing species BLAST databases")
parser.add_argument("-n", "--num-threads", default=1, help="Number of threads to run in parallel (default: 1)")
args = parser.parse_args()

groups_dir = args.ogroups_dir + "/"
hits_dir = args.hits_dir + "/"

rBBH_hits_dir = groups_dir.rstrip("/") + "_rBBHits/"
if not os.path.exists(rBBH_hits_dir):
    os.mkdir(rBBH_hits_dir)
    print("Created new directory " + rBBH_hits_dir)
rBBH_results = groups_dir.rstrip("/") + "_rBBH.txt"

monophyletic_groups = []
with open(args.monophyly, "r") as f:
    for line in f:
        if line.strip("\n").split("\t")[2] == "True":
            monophyletic_groups.append(line.split("\t")[0])

def find_node_for_monophyly(tree, node, sequences):
    
    if node == None:
        node = tree.seed_node
    
    if len(node.child_nodes()) == 0:
        return None
    
    leaf_nodes = set(l.taxon.label for l in node.leaf_iter())
    if leaf_nodes == set(sequences):
        return node
    
    all_leaves = set(l.taxon.label for l in tree.seed_node.leaf_iter())
    if all_leaves - leaf_nodes == set(sequences):
        return node
    
    mono_node = None
    for child_node in node.child_nodes():
        mono_node = find_node_for_monophyly(tree, child_node, sequences)
        if mono_node != None:
            break
    return mono_node 

def find_closest_homologues(tree, sequences):
    
    for edge in tree.preorder_edge_iter():
        if edge.length > 0:
            edge.length = 1
    
    pdm = tree.phylogenetic_distance_matrix()
    min_dist = sys.maxsize
    closest_homologues = set()
    for leaf1 in tree.leaf_node_iter():
        if leaf1.taxon.label not in sequences:
            continue
        for leaf2 in tree.leaf_node_iter():
            if leaf2.taxon.label not in sequences:
                if pdm(leaf1.taxon, leaf2.taxon) < min_dist:
                    min_dist = pdm(leaf1.taxon, leaf2.taxon)
                    closest_homologues = { leaf2.taxon.label.replace(" ", "_") }
                elif pdm(leaf1.taxon, leaf2.taxon) == min_dist:
                    closest_homologues.add(leaf2.taxon.label.replace(" ", "_"))
    
    return closest_homologues

def rBBH(group):
    
    orthologues = set()
    with open(groups_dir + group + ".fasta", "r") as f:
        for line in f:
            if line.startswith(">"):
                orthologues.add(line.strip("> \n").replace("_", " "))
    tree = Tree.get(path=hits_dir + group + ".tree", schema="newick")
    if len(list(tree.leaf_node_iter())) < 4:
        closest_homologues = set([node.taxon.label for node in tree.leaf_node_iter()]) - orthologues
        print(closest_homologues)
    else :
        closest_homologues = find_closest_homologues(tree, orthologues)
    
    hit_sequences = SeqIO.to_dict(SeqIO.parse(hits_dir + group + ".fasta", "fasta"))
    closest_homolog_sequences = []
    for seq in hit_sequences:
        if hit_sequences[seq].id in closest_homologues:
            closest_homolog_sequences.append(hit_sequences[seq])
    SeqIO.write(closest_homolog_sequences, hits_dir + group + "_closest_homologues.fasta", "fasta")
    
    rBBH_found = False
    species_set = set([o.split("|")[0] for o in orthologues])
    for species in species_set:
        db = args.database_dir + "/" + species
        #db = "blast_dbs/" + species
        max_hits = 1
        blast_out = rBBH_hits_dir + group + "_against_" + species + ".blast.out"
        if os.path.exists(blast_out):
            print("SKIPPING search: " + blast_out + " already exists!")
        else:
            blastp = NcbiblastpCommandline(query=hits_dir + group + "_closest_homologues.fasta", db=db, max_target_seqs=max_hits, outfmt=6, out=blast_out)
            stdout, stderr = blastp()
        with open(blast_out, "r") as f:
            for line in f:
                if line.split("\t")[1].replace("_", " ") in orthologues:
                    rBBH_found = True
                    break
        if rBBH_found:
            break
    print(group + " done!")
    return "{}\t{}\n".format(group, rBBH_found)

with concurrent.futures.ProcessPoolExecutor(max_workers=int(args.num_threads)) as executor:
    futures = [executor.submit(rBBH, group) for group in monophyletic_groups]

results = concurrent.futures.wait(futures)
out = "Group\trBBH found?\n"
for result in results.done:
    out += result.result()
    
with open(rBBH_results, "w") as f:
    f.write(out)

