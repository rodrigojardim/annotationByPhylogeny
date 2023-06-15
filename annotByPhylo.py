from Bio import Phylo
from Bio import SeqIO
import argparse
import logging

logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.INFO)
logging.info('Started')

ap = argparse.ArgumentParser()

# Add the arguments to the parser
ap.add_argument("-t", "--tree", required=True, help="input tree file newick format")
ap.add_argument("-f", "--fasta", required=True, help="fasta file")
ap.add_argument("-a", "--terms", required=True, help="terms for annotation")
ap.add_argument("-o", "--output", required=True, help="output fasta")
ap.add_argument("-u", "--unknown", required=True, help="sequences with no annotation")

args = vars(ap.parse_args())

#treeFile = "../data/uniprot_mammal_2022-12-27.fasta.filtered.rpsblast.domains.cdhit.renamed.mafft.fasttree"
#fastaFile = "../data/uniprot_mammal_2022-12-27.fasta.filtered.rpsblast.domains.cdhit.renamed"
#termsFile = "../data/headers_itol_annotation.txt"

treeFile = args['tree']
fastaFile = args['fasta']
termsFile = args['terms']


# Lista de termos para a anotacao
terms = open(termsFile).readlines()

# Lista de sequencias sem anotacao
unknowns = open(args['unknown']).readlines()

def getHeaders():
    ret = {}
    for sequence in SeqIO.parse(fastaFile,'fasta'):
        description = sequence.description
        ret[sequence.id] = {}
        ret[sequence.id]['description'] = description
        ret[sequence.id]['seq'] = sequence.seq
    return ret

def getParent(tree, clade):
    # Path ate o clado
    st = tree.get_path(clade)
    #print(st)
    
    # Ancestral
    ancestor = st[len(st)-2]
    #print(tree.get_path(ancestor))
   
    # Clado ancestral
    return next(tree.find_clades(ancestor))

def reannotated(cladeAncestor, seqId): 
    #ret = ''
    # Encontrando a anotacao
    for elements in cladeAncestor.clades:
        if (elements.name and elements.name != seqId):
            #print("{}\t{}\t{}".format(seqId, elements.name, headers[elements.name]['description']))
            ## Verifica se o header consta na lista de termos
            ## Se estiver, muda a anotacao
            for term in terms:
                h = term.split(',')[0]
                if (h in headers[elements.name]['description']):
                    return headers[elements.name]['description']
    return reannotated(getParent(tree,cladeAncestor), seqId)

###############################################

headers = getHeaders()

print('Search for unknowns')
for line in unknowns:
    seqId = line.split("\t")[0]
    logging.debug(seqId)

    # Abrir arvore
    tree = Phylo.read(treeFile, 'newick')
    logging.debug(tree)
    
    # Buscar o clado
    clade = next(tree.find_clades(name=seqId))
    logging.debug(tree.get_path(clade))
    
    # Clado ancestral
    cladeAncestor = getParent(tree, clade)
   
    # Reanotando
    try:
        newAnnotation = reannotated(cladeAncestor, seqId)
    except:
        newAnnotation = ''
    if (newAnnotation != ''):
        temp = ' '.join(newAnnotation.split(' ')[1:len(newAnnotation)-1])
        newAnnotation = temp + "\told_annotation\t" + headers[seqId]['description']
        logging.debug(newAnnotation)
        headers[seqId]['description'] = newAnnotation 
    else:
        logging.debug('No')
        pass
    #print(headers[seqId]['description'])
                    

## Reescrevendo o fasta
logging.info('Write reannotated fasta')
output = open(args['output'],"w")
for k in headers.keys():
    newHeader = ' '.join(headers[k]['description'].split(' ')[1:len(headers[k]['description'])])
    output.write(">{}\t{}\n{}\n".format(k, newHeader,headers[k]['seq'] ))    
output.close()
