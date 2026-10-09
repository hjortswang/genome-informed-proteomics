from biocommons.seqrepo import SeqRepo
from biocommons.seqrepo.fastadir.fabgz import FabgzReader
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils.CheckSum import seguid

from pathlib import Path
import subprocess
import os
import ast
import re
import pandas as pd
import numpy as np
import cdhit_reader
from statistics import mean,stdev
from scipy.stats import norm,kurtosis,skew

from scipy.signal import find_peaks
from scipy.stats import gaussian_kde

def load_SeqRepo(repo_path:str):

    basedir =Path(repo_path) 
    basedir.parent.mkdir(exist_ok=True)

    if not basedir.is_dir():
        proc = subprocess.Popen(["seqrepo", f"--root-directory", str(basedir.parent), "init"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
    sr = SeqRepo(str(basedir), writeable=True)
    return sr

def make_SeqRepo(fastapath:str,sr:SeqRepo,exclude=[],genome_ids=[]):
    """
     Here, one can either specify which genomes to use, or which genomes not to use. 
     genome_ids: what genomes to use
     exclude: what genomes to exclude
     if all files in fastapath should be used, one leaves exclude and genome_ids empty.
    """

    if not genome_ids:
        fastafiles=os.listdir(fastapath)
        genome_ids=[fastafile.replace('.fasta','') for fastafile in fastafiles]
    
    if exclude:
        for id in exclude:
            genome_ids.drop(exclude)
    
    for id in genome_ids:
        for record in SeqIO.parse(f'{fastapath}/{id}.fasta',"fasta"):
            description = record.description.replace(' ','__',1) # put a double underline in front of protein description, such that it is easily extracted downstream
            sr.store(nsaliases=[{"namespace": id, "alias": description}], seq=str(record.seq))
            sr.commit()

def get_repofiles(basedir):
    dir=Path(basedir)
    filelist=list(dir.glob('**/*.fa.bgz'))
    return filelist

def get_clean_seq(filelist:list,sr:SeqRepo):
    clean_sequences=[]
    for ix in range(0,len(filelist)):
        temp_fabgz = FabgzReader(filelist[ix]) # read file
        seq_id=temp_fabgz.keys() # get seguid
        seq=temp_fabgz.fetch(seq_id[0]) # translate sequid into sequence


        if not re.search('[XJOUZ]',seq):

            clean_sequences.append(filelist[ix]) #add sequence to list of "clean" sequences
    return clean_sequences  

def write_fasta(
    filelist: list,
    sr: SeqRepo,
    outfile: str,
    threshold=0,
    protein_id=False,
    exclude=None
):

    if exclude is None:
        exclude = []

    exclude = set(exclude)

    clean_alias = []
    clean_genomes = []
    clean_genomes_count = []
    clean_seguid = []
    clean_sequences = []
    clean_repo_files = []
    records = []

    for ix in range(len(filelist)):

        temp_fabgz = FabgzReader(filelist[ix])

        # Get sequence ID
        seq_id = temp_fabgz.keys()

        # Get sequence
        seq = temp_fabgz.fetch(seq_id[0])

        if len(seq) >= threshold:

            if not re.search("X", seq):

                # Get aliases for sequence
                alias = list(
                    sr.translate_alias(
                        seguid(seq)
                    )
                )

                # Keep PATRIC aliases
                peg = [
                    a for a in alias
                    if ".peg." in a
                ]

                # ------------------------------------------------
                # Exclude selected genomes
                # ------------------------------------------------

                if exclude:

                    peg = [
                        alias
                        for alias in peg
                        if not any(
                            genome in alias
                            for genome in exclude
                        )
                    ]

                # If the sequence only occurred in excluded
                # genomes, do not include the sequence
                if len(peg) == 0:
                    continue

                # ------------------------------------------------
                # Header
                # ------------------------------------------------

                if protein_id is False:

                    header = re.split(
                        "__",
                        peg[0]
                    )[1]

                else:

                    header = peg[0]

                # ------------------------------------------------
                # Genome IDs
                # ------------------------------------------------

                genome_id = [
                    a.split(":")[0]
                    for a in peg
                ]

                res = list(
                    set(genome_id)
                )

                # ------------------------------------------------
                # Save sequence information
                # ------------------------------------------------

                clean_sequences.append(
                    seq
                )

                clean_genomes.append(
                    res
                )

                clean_genomes_count.append(
                    len(res)
                )

                clean_alias.append(
                    peg
                )

                clean_seguid.append(
                    seguid(seq)
                )

                record = SeqRecord(
                    Seq(seq),
                    id=f"sequid:{seguid(seq)}",
                    description=header
                )

                records.append(
                    record
                )

                clean_repo_files.append(
                    filelist[ix]
                )

    # ============================================================
    # Repository dataframe
    # ============================================================

    d = {
        "seguid": clean_seguid,
        "sequence": clean_sequences,
        "aliases": clean_alias,
        "repository": clean_repo_files,
        "genome_id": clean_genomes,
        "nbr_genomes": clean_genomes_count
    }

    repo_frame = pd.DataFrame(
        d
    )

    # ============================================================
    # Save
    # ============================================================

    frame_name = (
        f"{outfile[:-6]}_frame.csv"
    )

    repo_frame.to_csv(
        frame_name
    )

    SeqIO.write(
        records,
        outfile,
        "fasta"
    )

    return repo_frame




def run_cdhit(inputfile:str,cutoff:str, sr:SeqRepo):#,edit_header=True):

    """
    make new cdhit repo, if it does not already exist.

    make new name for the edited fasta-file 
    Run cdhit
    load clusters from cdhit
    for every one of the clusters:
    
    find sequence corresponding to the 
    """
    newdir = f'../data/cdhit_{cutoff}'
    try:
        os.mkdir(newdir)
    except FileExistsError:
        pass
        
    
    outputfile= inputfile.replace('.fasta',f'_{cutoff}.fasta')
    outputfile_edited= inputfile.replace('.fasta',f'_{cutoff}_edited.fasta')
    
    os.system(f'cdhit -i {inputfile} -o {newdir}/{outputfile} -c 0.{cutoff} -d 100 > {newdir}/log')

  
    clusters = cdhit_reader.read_cdhit(f'{newdir}/{outputfile}.clstr').read_items()
    
    records=[]
    for cluster in clusters:
        seq=sr.fetch(cluster.refname)
        cluster.name= cluster.name.replace(' ','_')
        record=SeqRecord(Seq(seq),description=cluster.name, id=cluster.name)
        records.append(record)
    
    SeqIO.write(records, f'{newdir}/{outputfile_edited}', "fasta")


def mk_ClusterTable(
        cutoff: str,
        clusterfile: str,
        repofile: str,
        sr: SeqRepo):

    clusters = cdhit_reader.read_cdhit(
        clusterfile
    ).read_items()

    repo_frame = pd.read_csv(repofile)

    output_file = clusterfile.replace(
        '.fasta.clstr',
        '_cluster_frame.csv'
    )

    try:
        repo_frame['genome_id'] = (
            repo_frame['genome_id'].apply(eval)
        )
    except:
        pass

    ref_frame= pd.DataFrame(columns=['cluster', 'cluster size', 'names','nbr unique names', 'nbr_genomes','genomes',
            'duplicates', 'nbr_duplicates', 'seguid reference',  'seq reference',
           'cluster seguids', 'cluster sequences',  'occurences','seguid weight','seq length', 'all length','mean length', 'std length',
           'length interval', 'length interval width'
           ],index=range(0, len(clusters)))

    count =0

    for cluster in clusters:
        ref_frame.loc[count,'cluster']=cluster.name
        ref_frame.loc[count,'cluster size']=len(cluster)

        # ====================================================
        # Get SEGUIDs
        # ====================================================

        if cluster.refname.startswith('seguid:') or cluster.refname.startswith('sequid:'):
            cluster.refname=re.split(':',cluster.refname)[1]
            seguids= [member.name[7:] for member in cluster.sequences]

        else:
            seguids= [member.name for member in cluster.sequences]

        # ====================================================
        # Sequences
        # ====================================================

        seqs = [sr.fetch(id) for id in seguids]
        length= [member.length for member in cluster.sequences]
        ref_frame.loc[count,'seguid reference']=cluster.refname
        ref_frame.loc[count,'seq reference']=sr.fetch(cluster.refname)
        ref_frame.at[count,'seq length']=length


        # ====================================================
        # Get ALL names associated with the unique sequences
        # in the cluster
        # ====================================================

        cluster_names = []

        for seguid in seguids:

            matches = repo_frame.loc[
                repo_frame["seguid"] == seguid,
                "aliases"
            ]

            if len(matches) == 0:
                continue

            aliases = matches.iloc[0]

            # Convert string representation to list
            try:
                aliases = ast.literal_eval(
                    aliases
                )
            except (ValueError, SyntaxError, TypeError):
                continue

            names = [
                alias.split("__", 1)[1]
                for alias in aliases
                if "__" in alias
            ]

            cluster_names.extend(
                names
            )

        # Remove duplicates while preserving order
        cluster_names = list(
            dict.fromkeys(cluster_names)
        )

        ref_frame.loc[
            count,
            "names"
        ] = "; ".join(
            cluster_names
        )

        ref_frame.loc[
            count,
            "nbr unique names"
        ] = len(
            cluster_names
        )

        seq_length=[]
        for i in range(len(seguids)):
            genomes = repo_frame.loc[repo_frame['seguid']==seguids[i],'genome_id'].values[0]
            j=0
            while j< len(genomes):
                seq_length.append(length[i])
                j+=1
        ref_frame.at[count,'all length']=seq_length
        ref_frame.at[count,'cluster seguids']=seguids
        ref_frame.at[count,'cluster sequences']=seqs   
        
        
        
        ref_frame.loc[count,'mean length']= mean(seq_length)
        length_width= max(seq_length)-min(seq_length)
        ref_frame.loc[count,'length interval width']= length_width
        ref_frame.loc[count,'length interval']= [min(seq_length), max(seq_length)]

        # ====================================================
        # Optional statistics
        # ====================================================
        if False:
            #standard deviation
            if (len(seq_length)>=2) & (length_width>2):
                ref_frame.loc[count,'std length']= stdev(seq_length)
            elif len(seq_length)<2: 
                ref_frame.loc[count,'std length']= 0
    
            # Skewness
            if len(seq_length) >= 3 and np.ptp(seq_length) > 0:
                ref_frame.loc[count, 'skewness length'] = skew(seq_length)
            else:
                ref_frame.loc[count, 'skewness length'] = 0
    
            # Kurtosis
            if len(seq_length) >= 4 and np.ptp(seq_length) > 0:
                ref_frame.loc[count, 'kurtosis length'] = kurtosis(seq_length)
            else:
                ref_frame.loc[count, 'kurtosis length'] = 0
    
            if (len(seq_length)>10) & (length_width>2):
                
                try:
                    # KDE estimation
                    kde = gaussian_kde(seq_length)
                    x = np.linspace(min(seq_length), max(seq_length), 2000)
                    kde_values = kde(x)
    
                    # Detect peaks
                    peaks, _ = find_peaks(kde_values)
    
                    ref_frame.loc[count,'nbr length peaks']= len(peaks)
                    ref_frame.loc[count,'peak points']= peaks
                except Exception as e: 
                    print(f'Error in cluster {count}: {e}, length: {seq_length}')
                    ref_frame.to_csv(clusterfile.replace('.fasta.clstr','cluster_frame_broken.csv'))
                    break


        # ====================================================
        # Genome occurrence information
        # ====================================================

        duplicates = 0
        duplicate_ids=[]
        genome_ids=[]
        nbr_ids=[]

        for id in seguids:
            genomes = repo_frame.loc[repo_frame['seguid']==id,'genome_id'].values[0]
            nbr_ids.append(len(genomes)) # number of ids connected to seguid
            for genome in genomes:
                if genome not in genome_ids:
                    genome_ids.append(genome)
                elif genome in genome_ids:
                    duplicate_ids.append(genome)
                    duplicates+=1
        norm_nbr_ids = [float(i)/sum(nbr_ids) for i in nbr_ids]
        ref_frame.at[count,'genomes']=genome_ids
        ref_frame.at[count,'occurences']=nbr_ids
        ref_frame.at[count,'seguid weight']=norm_nbr_ids
        ref_frame.loc[count,'nbr_genomes']=int(len(genome_ids))
        ref_frame.at[count,'duplicates']=duplicate_ids
        ref_frame.loc[count,'nbr_duplicates']=duplicates
        count+=1

    # ========================================================
    # Save
    # ========================================================

    output_file=clusterfile.replace('.fasta.clstr','_cluster_frame.csv')
    ref_frame.to_csv(f'{output_file}')