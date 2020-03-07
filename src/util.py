import os
# import sys
import numpy as np
import struct
import matplotlib
from scipy.sparse import csr_matrix
from scipy.special import logsumexp
import pandas as pd
import struct
import pickle


def iter_loadtxt(filename, delimiter='\t', skiprows=0, dtype=float):
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
    data = np.genfromtxt(file, dtype=None, encoding=None)
    return data


# Read bedgraph format data.
def readBedGraph(file):
    data = np.genfromtxt(file, dtype=None, encoding=None)
    return data


def readData(file):
    data = np.loadtxt(file)
    return data


def readHiC(file):
    # data = iter_loadtxt((file)
    # hic_data = np.loadtxt(args.hic)
    # hic_data = util.iter_loadtxt(args.hic)
    data = pd.read_csv(file, delimiter="\t").values

    return data

def create_hic_matrix(hic_input, n, weight=False):
    hic_data = readHiC(hic_input)
    #hic_data = hic_data[hic_data[:,0]<hic_data[:,1]]
    hic_data_swap = hic_data.copy()
    hic_data_swap[:,[0, 1]] = hic_data_swap[:,[1, 0]]
    hic_data_merge = np.concatenate((hic_data, hic_data_swap), axis=0)
    #if weight:
    #    hic_matrix = csr_matrix((hic_data[:, 2], (hic_data[:, 0], hic_data[:, 1])), shape=(n, n))
    #else:
    #    hic_matrix = csr_matrix((np.ones(np.size(hic_data,0)), (hic_data[:, 0], hic_data[:, 1])), shape=(n, n))
    #print(edges.toarray()[1:20, 1:20])

    return hic_data_merge

def log_transform(data, pseudocount=1e-100):
    #pseudocount = 1e-10
    return np.log(data + pseudocount)

def sparse_logsumexp(m):
    return logsumexp(m.toarray())

def sparse_logsumexp_row(m):
    return logsumexp(m.toarray(),axis=0)

# Modified from straw: https://github.com/theaidenlab/straw
def get_hic_chr(file):
    # print(file)

    hic_file = open(file, 'rb')
    magic_string = struct.unpack('<3s', hic_file.read(3))[0]
    hic_file.read(1)
    version = struct.unpack('<i', hic_file.read(4))[0]
    print('HiC version:')
    print('  {0}'.format(str(version)))

    masterindex = struct.unpack('<q', hic_file.read(8))[0]
    genome = ""
    c = hic_file.read(1)

    while (c != b'\0'):
        genome += c.decode('ISO-8859-1')
        c = hic_file.read(1)
    print('Genome ID:')
    print('  {0}'.format(str(genome)))

    # read and throw away attribute dictionary (stats+graphs)
    # print('Attribute dictionary:')
    nattributes = struct.unpack('<i', hic_file.read(4))[0]
    for x in range(0, nattributes):
        key = readcstr(hic_file)
        value = readcstr(hic_file)
    nChrs = struct.unpack('<i', hic_file.read(4))[0]
    print("Chromosomes: ")

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
    outfile = open(output_name, "wb")
    pickle.dump(variable, outfile)
    outfile.close()


def print_log(text, file_name):
    outfile = open(file_name, "a+")
    outfile.write(text + "\n")
    outfile.close()



def juicer_dump(juicer, type, norm_method, hic_file, chr1, chr2, resolution, output):
    print("Running juicer tools dump:")
    command = "java -jar " + juicer + " dump " + type + " " + norm_method + " " + hic_file + " " \
              + chr1 + " " + chr2 + " BP " + str(resolution) + " " + output

    print(command)
    os.system(command)

    return


# Add bin number to hic interactions
def hic_add_bin_num(hic_file, chr1, chr2, genomic_bin_file, bin_size, output):
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


# Dump all hic interactions.
def dump_hic_all(juicer, hic_file, chr_list, genomic_bin_file, resolution, output):
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


