#!/usr/bin/env python

import argparse
import re

def escape_string(s):
    return re.sub(r'"', r'\"', s)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType(), help='KEGG brite file')

    args = parser.parse_args()
    print_indent = ' '*4
    level = 1
    print '---'
    for line in args.infile:
        line = line.rstrip()
        match = re.search(r'\s+(K\d+)\s+(.*);\s+(.*)', line)
        if match:
            genes = match.group(2).split(',')
            genes = [x.lstrip() for x in genes]
            product = match.group(3)
            kid = match.group(1)
            ec_match = re.match(r'(.*)\s+\[EC:(.*)\]$', product)
            #print "    %s- %s" % (print_indent, kid)
            if ec_match:
                ecs = ec_match.group(2).split()
                print "%s:\n%sgene: %s\n%sproduct: \"%s\"\n%sEC_number: %s" % (kid, print_indent, str(genes), print_indent, escape_string(ec_match.group(1)), print_indent, str(ecs))
            else:
                print "%s:\n%sgene: %s\n%sproduct: \"%s\"" % (kid, print_indent, str(genes), print_indent, escape_string(product))
