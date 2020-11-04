# SPIN - An PGM framework to predict nuclear organization
# By Yuchuan Wang
# yuchuanw@andrew.cmu.edu

import sys
import os
import numpy as np
import time
import argparse
import util
import mrf

import pandas as pd
from scipy.sparse import csr_matrix


# Parse arguments
def parse_arg():
    p = argparse.ArgumentParser()
    p.add_argument('-i', help="1D genomic measurements of nuclear organization.")
    p.add_argument('--hic', help="Hi-C interactions.")
    p.add_argument('-w', help="Resolution.", type=int, default=100000)
    p.add_argument('-n', help="Number of states.", type=int, default=5)
    #p.add_argument('-m', help="Choose mode. Supported: full, hic", type=string, default="full")
    p.add_argument('-o', help="Output dir.")
    p.add_argument('-p', help="Number of processes.", type=int)
    p.add_argument('-g', help="Genome bin file.")
    p.add_argument("--prev", help="reload existing model.")
    p.add_argument("--save", help="save model.", action='store_true')
    # p.add_argument('-j', help="Juicer tools dir.")

    # missing arguments
    if len(sys.argv) < 5:
        p.print_help()
        exit(1)
    return p.parse_args()

def run(args):
    
    # create output folder
    if not os.path.isdir(args.o):
        os.makedirs(args.o)
        print("Create output dir: %s" % args.o)
    else:
        sys.exit("Error: Output dir %s already exits." % args.o)

    util.print_log(time.ctime() + " Start running.", args.o + "/log.txt")

    # Print args().
    util.print_log("Input 1D signal file: %s" % (args.i), args.o + "/log.txt")
    util.print_log("Input Hi-C file: %s" % (args.hic), args.o + "/log.txt")
    util.print_log("Input bin size file: %s" % (args.g), args.o + "/log.txt")
    util.print_log("Resolution: %s nt" % (args.w), args.o + "/log.txt")
    util.print_log("Number of states: %s " % (args.n), args.o + "/log.txt")
    util.print_log("Number of cores to use: %s" % (args.p), args.o + "/log.txt")

    # Read data
    bin_data = util.readBedGraph(args.g)
    input_data = util.readData(args.i)
    util.print_log(time.ctime() + " Finished reading input.", args.o + "/log.txt")

    # Create Hi-C Matrix
    (n, d) = input_data.shape
    edges = util.create_hic_matrix(args.hic, n)
    util.print_log(time.ctime() + " Finished creating edges.", args.o + "/log.txt")

    # Create graph object
    hmrf = mrf.MarkovRandomField(n=n, edges=edges, obs=input_data,args=args)
    
    # initialization
    hmrf.init_gmm()
    print(hmrf.label)
    util.print_log(time.ctime() + " Init GMM.", args.o + "/log.txt")
    hmrf.init_trans()
    util.print_log(time.ctime() + " Init trans matrix.", args.o + "/log.txt")
    print(hmrf.edge_potential)
    
    # inference of states
    hmrf.solve()
    
    # Save file
    np.savetxt(args.o + "/state_" + str(hmrf.n_state), hmrf.label, delimiter='\n', fmt='%d')
    
    if args.save:
        util.print_log(time.ctime() + " Save model.", args.o + "/log.txt")
        util.save_variable(hmrf, args.o + "/model.pkl")

# Main function
if __name__ == '__main__':
    
    start = time.time()
    
    args = parse_arg()

    run(args)

    # End.
    end = time.time()
    util.print_log("Total running time: %.2f s." % (end - start), args.o + "/log.txt")
    