import argparse, os, subprocess, concurrent.futures, multiprocessing, sys
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
from dendropy import Tree

work_dir = os.getcwd() + "/"

parser = argparse.ArgumentParser(description='Tests orthologous group for monophyly after searching for potential homologues')
parser.add_argument("-od", "--ogroups-dir", required=True, help="folder containing orthologous groups (FASTA format)")
parser.add_argument("-db", "--database", required=True, help="database location to query for potential homologues")
parser.add_argument("-e", "--e-value", default=10.0, help="e-value used to search for potential homologues")
parser.add_argument("-o", "--output", default=work_dir+"output.log", help="output log file (default: output.log)")
parser.add_argument("-n", "--num-threads", default=1, help="Number of threads to run in parallel (default: 1)")
args = parser.parse_args()

groups_dir = args.ogroups_dir + "/"

blast_dir = groups_dir.rstrip("/") + "_blast_results/"
hits_dir = groups_dir.rstrip("/") + "_hits/"
no_hits_dir = groups_dir.rstrip("/") + "_no_hits/"
db = args.database

if not os.path.exists(hits_dir):
    os.mkdir(hits_dir)
    print("Created new directory " + hits_dir)

if not os.path.exists(no_hits_dir):
    os.mkdir(no_hits_dir)
    print("Created new directory " + no_hits_dir)

if not os.path.exists(blast_dir):
    os.mkdir(blast_dir)
    print("Created new directory " + blast_dir)

def homFind(group_file):
    if not group_file.endswith(".fasta"):
        return
    blast_file = blast_dir + group_file.replace(".fasta", ".blast.out")
    if os.path.exists(blast_file):
        print("SKIPPING reblast for " + group_file)
    else:
        with open(groups_dir + group_file, "r") as f:
            content = f.read()
            max_hits = int(100/int(content.count(">")) + 1)
        blastp = NcbiblastpCommandline(query=groups_dir+group_file, db=db, max_target_seqs=max_hits, outfmt=5, out=blast_file, evalue=float(args.e_value))
        stdout, stderr = blastp()

    if  os.path.exists(hits_dir + group_file):
        print("SKIPPING collecting results for " + group_file)
    else:
        found_hits = dict()
        blast_results = list(NCBIXML.parse(open(blast_file)))
        
        for result in blast_results:
            for alignment in result.alignments:
                name = alignment.title.split(" ", 1)[1]
                if name not in found_hits:
                    found_hits[name] = alignment.hsps[0].sbjct.replace("-", "").replace("\n", "")
        genes = list(SeqIO.parse(groups_dir + group_file, "fasta"))
        for sequence in found_hits:
            # replace non-standard characters (J, U, O) and only use RefSeq ID as sequence name
            genes.append(SeqRecord(Seq(found_hits[sequence].replace("J", "X").replace("U", "X").replace("O", "X")), id=sequence.split(" ")[0], description=""))
        if len(found_hits) < 1:
            print("No potential homologues found for " + group_file + "!")
            SeqIO.write(genes, no_hits_dir + group_file, "fasta")
            return
        else:
            SeqIO.write(genes, hits_dir + group_file, "fasta")
    
    msa = hits_dir + group_file.replace(".fasta", ".msa")
    if os.path.exists(msa):
        print("SKIPPING multiple sequence alignment for " + group_file)
    else:
        clustalomega_cline = ClustalOmegaCommandline(infile=hits_dir+group_file, outfile=msa, auto=True)
        stdout,stderr = clustalomega_cline()
    
    sequences = []
    with open(hits_dir + group_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                sequences.append(line.strip("> \n").split(" ")[0])
    
    tree_file = msa.replace(".msa", ".tree")
    tree_file_name = tree_file.rsplit("/", 1)[1]

    if len(sequences) < 4:
        print("SKIPPING tree building for " + group_file + ": not enough sequences - " + str(len(sequences)))
        with open(tree_file, "w") as f:
            f.write("(" + ":1,".join(sequences) +":1);")

    if os.path.exists(tree_file):
        print("SKIPPING tree building for " + group_file)
    else:
        p = subprocess.Popen(["raxml-ng", "--msa", msa, "--model", "PROTGTR+G4", "--data-type", "AA", "--threads", "1"])
        print(tree_file, tree_file_name, p.args)
        p.wait()
        os.rename(msa + ".raxml.bestTree", tree_file)
        os.rename(msa + ".raxml.startTree", tree_file.replace(".tree", ".raxml.startTree"))
        os.rename(msa + ".raxml.bestModel", tree_file.replace(".tree", ".raxml.bestModel"))
        os.rename(msa + ".raxml.log", tree_file.replace(".tree", ".raxml.log"))

    print(group_file + " done!")

with concurrent.futures.ProcessPoolExecutor(max_workers=int(args.num_threads)) as executor:
    futures = [executor.submit(homFind, fasta) for fasta in os.listdir(groups_dir) if fasta.endswith(".fasta")]

results = concurrent.futures.wait(futures)
for rd in results.done:
    print(rd.result())

def check_for_monophyly(tree, node, sequences):
    if node == None:
        node = tree.seed_node
    if len(node.child_nodes()) == 0:
        return False
    leaf_nodes = set(l.taxon.label for l in node.leaf_iter())
    if leaf_nodes == set(sequences):
        return True
    all_leaves = set(l.taxon.label for l in tree.seed_node.leaf_iter())
    if all_leaves - leaf_nodes == set(sequences):
        return True
    monophyletic = False
    for child_node in node.child_nodes():
        monophyletic = monophyletic or check_for_monophyly(tree, child_node, sequences)
        if monophyletic:
            break
    return monophyletic

out = "Group\t#Orthologs\tMonophyletic?\n"
for group_file in os.listdir(hits_dir):
    if not group_file.endswith(".fasta") or group_file.endswith("_closest_homologues.fasta"):
        continue

    orthologs = []
    with open(groups_dir + group_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                orthologs.append(line.strip("> \n").replace("_", " "))
    tree = Tree.get(path=hits_dir + group_file.replace(".fasta", ".tree"), schema="newick")
    if len(list(tree.leaf_node_iter())) < 4:
        monophyletic = True
    else:
        monophyletic = check_for_monophyly(tree, None, orthologs)
    out += "{}\t{}\t{}\n".format(group_file.replace(".fasta", ""), len(orthologs), monophyletic)
    if not monophyletic:
        for node in tree.mrca(taxon_labels = orthologs).leaf_nodes():
            if node.taxon.label not in orthologs:
                node.taxon.label = "BREAKING:"+node.taxon.label
        tree.write(path=hits_dir + group_file.replace(".fasta", ".tree.tagged"), schema="newick")
    print(group_file)

with open(groups_dir.rstrip("/")+"_monophyly.txt", "w") as f:
    f.write(out)
