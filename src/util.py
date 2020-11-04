import os
# import sys
import numpy as np
import struct
import matplotlib
from scipy.sparse import csr_matrix
from scipy.special import logsumexp
from scipy.stats import multivariate_normal
import pandas as pd
import struct
import pickle

def iter_loadtxt(filename, delimiter='\t', skiprows=0, dtype=float):
    """ Read from text line by line.
    
    Args:
        skiprows (int): Skip the first $skiprows rows.
        dtype: Data type of loaded data. Default: float.
    
    Returns:
        data: Numpy array of loaded data.
        
    """
    def iter_func():
        with open(filename, 'r') as infile:
            for _ in range(skiprows):
                next(infile)
            for line in infile:
                line = line.rstrip().split(delimiter)
                for item in line:
                    yield dtype(item)
        iter_loadtxt.rowlength = len(line)

    data = np.fromiter(iter_func(), dtype=dtype)
    data = data.reshape((-1, iter_loadtxt.rowlength))
    return data


def read_file(file):
    """ Load data with numpy.genfromtxt. """
    data = np.genfromtxt(file, dtype=None, encoding=None)
    return data


def readBedGraph(file):
    """ Read bedgraph format data. """
    data = np.genfromtxt(file, dtype=None, encoding=None)
    return data

def readData(file):
    """ Read input data (no missing data). """
    data = np.loadtxt(file)
    return data

def readHiC(file):
    """ Read hi-c/interaction input. """
    # data = iter_loadtxt((file)
    # hic_data = np.loadtxt(args.hic)
    # hic_data = util.iter_loadtxt(args.hic)
    data = pd.read_csv(file, delimiter="\t").values

    return data
 
def create_hic_matrix(hic_input, n):
    """ Create hi-c input matrix, consider upper-left/lower-right case. 
        
    Args:
        hic_input (str): File name of hi-c input file.
        n (int) : Number of nodes/bins.
    Returns:
        hic_data_merge: Numpy array of hi-c interaction data.
        
    """
    hic_data = readHiC(hic_input)
    #hic_data = hic_data[hic_data[:,0]<hic_data[:,1]]
    hic_data_swap = hic_data.copy()
    hic_data_swap[:,[0, 1]] = hic_data_swap[:,[1, 0]]
    hic_data_merge = np.concatenate((hic_data, hic_data_swap), axis=0)

    return hic_data_merge

def log_transform(data, pseudocount=1e-100):
    """ log transform of data. 
    
    Args:
        data : Numpy array.
        pseudocount : pseudocount to avoid divide by 0.
    Returns:
        Numpy array of log transformed data.
        
    """
    return np.log(data + pseudocount)

def sparse_logsumexp(m):
    """ Calculate the log of the sum of exponentials. """
    return logsumexp(m.toarray())

def sparse_logsumexp_row(m):
    """ Calculate the log of the sum of exponentials for row. """
    return logsumexp(m.toarray(),axis=0)

def get_multi_normal_pdf(mean, cov, data, log=True):
    """ Get Multivariate normal distribution pdf. 
    
    Args:
        mean : Numpy array of mean.
        cov : Numpy array of covariance.
        data : Numpy array of data points.
        log : Return pdf in log space. Default: True.
        
    Returns:
        Multivariate normal distribution pdf.
    """
    multi_normal = multivariate_normal(mean=mean, cov=cov, allow_singular=True)
    if log:
        return multi_normal.logpdf(data)
    else:
        return multi_normal.pdf(data)

def get_hic_chr(file):
    """ Get chromosome list from .hic file. Modified from straw: https://github.com/theaidenlab/straw. """
    # print(file)

    hic_file = open(file, 'rb')
    magic_string = struct.unpack('<3s', hic_file.read(3))[0]
    hic_file.read(1)
    version = struct.unpack('<i', hic_file.read(4))[0]
    print('HiC version:')
    print('  {0}'.format(str(version)))

    """ Genome version """
    masterindex = struct.unpack('<q', hic_file.read(8))[0]
    genome = ""
    c = hic_file.read(1)

    while (c != b'\0'):
        genome += c.decode('ISO-8859-1')
        c = hic_file.read(1)
    print('Genome ID:')
    print('  {0}'.format(str(genome)))

    """ read and throw away attribute dictionary (stats+graphs) """
    # print('Attribute dictionary:')
    nattributes = struct.unpack('<i', hic_file.read(4))[0]
    for x in range(0, nattributes):
        key = readcstr(hic_file)
        value = readcstr(hic_file)
    nChrs = struct.unpack('<i', hic_file.read(4))[0]
    print("Chromosomes: ")
    
    """ get chromosome list """
    chr_list = np.array([])

    for x in range(0, nChrs):
        name = readcstr(hic_file)
        length = struct.unpack('<i', hic_file.read(4))[0]
        # print('  {0}  {1}'.format(name, length))
        # print(name)
        name = name[2:name.__len__() - 1]
        print(name, length)
        chr_list = np.append(chr_list, name)

    chr_list_new = np.delete(chr_list, 0)
    for i in range(chr_list_new.size):
        chr_list_new[i] = "chr" + chr_list_new[i]

    return chr_list_new

def save_variable(variable, output_name):
    """ Save model/variable to pickle file.
    
    Args:
        variable : Name of variable to save.
        output_name : Name of saved file.
        
    """
    outfile = open(output_name, "wb")
    pickle.dump(variable, outfile)
    outfile.close()


def print_log(text, file_name):
    """ Print text to log file. """
    outfile = open(file_name, "a+")
    outfile.write(text + "\n")
    outfile.close()
    

def juicer_dump(juicer, type, norm_method, hic_file, chr1, chr2, resolution, output):
    """ Extract data from .hic files. 
        
    Args:
        juicer : Path to juicer tools.
        type : Type of data to extract. Should be "observed" or "oe".
        norm_method: Normalization method to use.
        hic_file : Path to .hic file.
        chr1 : Name of chrmosome 1.
        chr2 : Name of chrmosome 2.
        resolution : Hi-C resolution to use.
        output :  Path to output file.
    
    """
    print("Running juicer tools dump:")
    command = "java -jar " + juicer + " dump " + type + " " + norm_method + " " + hic_file + " " \
              + chr1 + " " + chr2 + " BP " + str(resolution) + " " + output

    print(command)
    os.system(command)

    return


def hic_add_bin_num(hic_file, chr1, chr2, genomic_bin_file, bin_size, output):
    """ Add bin number to hic interactions.
            
    Args:
        hic_file : Path to Hi-C interaction file.
        chr1 : Name of chrmosome 1.
        chr2 : Name of chrmosome 2.
        genomic_bin_file : Path to genomic bin file. 4th column should be bin number.
        bin_size : Bin size to use.
        output :  Path to output file.
    
    """ 
    bins = readGenomicBin(genomic_bin_file)
    hic = np.genfromtxt(hic_file, dtype=None, encoding=None)
    # colnames_bins = bins.dtype.names
    # colnames_hic = hic.dtype.names
    (n_total,) = hic.shape
    # print(n_total)

    bin_dict = {(chr, start): index for (chr, start, end, index) in bins}

    # print(bin_dict)

    if os.path.exists(output):
        os.remove(output)
    else:
        pass

    f_handle = open(output, 'wb')

    """ Write to file. """ 
    n = 0
    for (start1, start2, interaction) in hic:
        if (chr1, start1) in bin_dict and (chr2, start2) in bin_dict and (not np.isnan(interaction)):

            tmp = np.array([[str(bin_dict[(chr1, start1)]), str(bin_dict[(chr2, start2)]), str(interaction)]])
            # hic_index = np.vstack((hic_index, tmp))
            # print(tmp)
            np.savetxt(f_handle, tmp, delimiter='\t', fmt='%s')
            n = n + 1
            if n % 100000 == 0:
                print(str(n) + "/" + str(n_total) + " processed.")

    # np.savetxt(output, hic_index, delimiter='\t', fmt='%s')
    f_handle.close()

    return

def plot_oe_matrix(hic_file, chr1, chr2, start1, end1, start2, end2, output):
    """ Save hic interactions in requested region. 
    
    Args:
        hic_file : Path to Hi-C interaction file.
        chr1 : Name of chrmosome for region 1.
        chr2 : Name of chrmosome for region 2.
        start1 : Start of region 1.
        end1 : End of region 1.
        start2 : Start of region 2.
        end2 : End of region 2.
        output :  Path to output file.
    
    """ 
    hic = np.genfromtxt(hic_file, dtype=None, encoding=None)
    matrix = np.zeros((end1 - start1 + 1, end2 - start2 + 1))

    for (bin1, bin2, interaction) in hic:
        if bin1 >= start1 and bin1 <= end1:
            if bin2 >= start2 and bin2 <= end2:
                matrix[bin1 - start1][bin2 - start2] = interaction
        if bin2 >= start1 and bin2 <= end1:
            if bin1 >= start2 and bin1 <= end2:
                matrix[bin2 - start1][bin1 - start2] = interaction


    np.savetxt(output + ".tmp", matrix, delimiter='\t', fmt='%s')


def dump_hic_all(juicer, hic_file, chr_list, genomic_bin_file, resolution, output):
    """ Dump all hic interactions. 
    
    Args:
        juicer : Path to juicer tools.
        hic_file : Path to .hic file.
        chr_list : List of chrmosomes to dump.
        genomic_bin_file : Path to genomic bin file. 
        resolution : Hi-C resolution to use.
        output :  Path to output file.
    
    """ 
    if not os.path.isdir(output):
        os.makedirs(output)
    else:
        pass

    hic_chr_list = get_hic_chr(hic_file)
    # print(hic_chr_list)
    for index1 in range(hic_chr_list.size):
        for index2 in range(index1, hic_chr_list.size):
            if (hic_chr_list[index1] in chr_list) and (hic_chr_list[index2] in chr_list):
                print(hic_chr_list[index1], hic_chr_list[index2])
                type = "observed"
                norm_method = "KR"
                chr1 = hic_chr_list[index1]
                chr2 = hic_chr_list[index2]
                output_target = output + "/" + chr1 + "_" + chr2 + ".observed"

                juicer_dump(juicer, type, norm_method, hic_file,
                            chr1, chr2, resolution, output_target)
                hic_add_bin_num(output_target, chr1, chr2,
                                genomic_bin_file, resolution, output_target + ".txt")
                os.remove(output_target)

                if chr1 == chr2:
                    type = "oe"
                    output_target = output + "/" + chr1 + "_" + chr2 + ".oe"

                    juicer_dump(juicer, type, norm_method, hic_file,
                                chr1, chr2, resolution, output_target)
                    hic_add_bin_num(output_target,
                                    chr1, chr2, genomic_bin_file, resolution, output_target + ".txt")
                    os.remove(output_target)


