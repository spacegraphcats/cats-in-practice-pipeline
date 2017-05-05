#! /usr/bin/env python
import sys
import argparse
import os.path
import sourmash_lib
from sourmash_lib import signature
import csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('benchmark_sig')
    p.add_argument('signature_location')
    p.add_argument('overhead_values', nargs='+', type=float)
    p.add_argument('-o', '--output', default=None)
    args = p.parse_args()

    # load reference sig
    benchmark_sig = signature.load_one_signature(args.benchmark_sig)
    print('loaded benchmark {}'.format(args.benchmark_sig), file=sys.stderr)
    benchmark = benchmark_sig.minhash


    # load frontier sigs
    mh_list = []
    for overhead in args.overhead_values:
        overhead = float(overhead)
        p = int(overhead * 100)
        filename = 'frontier.p{:02d}.sig'.format(p)
        filename = os.path.join(args.signature_location, filename)
        with open(filename) as fp:
            sig = signature.load_one_signature(fp)
            mh_list.append((overhead, sig.minhash))
            print('loaded {}'.format(filename), file=sys.stderr)

    print(len(mh_list))

    mh_list.sort()

    fieldnames = ['overhead', 'containment', 'similarity']
    outfp = sys.stdout
    if args.output:
        outfp = open(args.output, 'w')
    w = csv.DictWriter(outfp, fieldnames)
    w.writeheader()
                           
    for (overhead, mh) in mh_list:
        containment = benchmark.contained_by(mh)
        similarity = benchmark.similarity(mh)
        w.writerow(dict(overhead=overhead, containment=containment,
                        similarity=similarity))


if __name__ == '__main__':
    main()
    
