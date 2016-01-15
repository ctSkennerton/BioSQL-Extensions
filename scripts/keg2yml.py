#!/usr/bin/env python

import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType(), help='KEGG brite file')

    args = parser.parse_args()
    indent = 4
    level = 0
    print '---'
    for line in args.infile:
        line = line.rstrip()
        if line[0] in '!#+':
            continue

        if line[0] == 'A':
            name = re.sub(r'<.?b>', '', line[1:])
            if name == 'Signature module':
                break

            print '"%s":' % name
        elif line[0] == 'B':
            level = 1
            fields = line.split(None, 1)
            if len(fields) == 1:
                continue
            else:
                print_indent = ' '*level*indent
                name = re.sub(r'<.?b>', '', fields[1])
                print '%s"%s":' % (print_indent, name)
        elif line[0] == 'C':
            level = 2
            print_indent = ' '*level*indent
            fields = line.split(None, 1)
            if len(fields) == 1:
                continue
            else:
                match = re.match(r'(\d+)\s+(.*)\s+\[.*\]$', fields[1])
                if match:
                    print '%sko%s:\n    %sdefinition: "%s"\n%s    orthologs:' %(print_indent, match.group(1), print_indent, match.group(2), print_indent)
                else:
                    print '%s"%s":' % (print_indent, fields[1])
        elif line[0] == 'D':
            level = 3
            fields = line.split(None, 1)
            print_indent = ' '*level*indent
            if len(fields) == 1:
                continue
            else:
                match = re.match(r'(K\d+)\s+(.*);\s+(.*)', fields[1])
                if match:
                    genes = match.group(2).split(',')
                    genes = [x.lstrip() for x in genes]
                    product = match.group(3)
                    kid = match.group(1)
                    ec_match = re.match(r'(.*)\s+\[EC:(.*)\]$', product)
                    print "    %s- %s" % (print_indent, kid)
                    if ec_match:
                        ecs = ec_match.group(2).split()
                        #print "%s- %s:\n%sgene: %s\n%sproduct: \"%s\"\n%sEC_number: %s" % (print_indent, kid, print_indent, str(genes), print_indent, ec_match.group(1), print_indent, str(ecs))
                    else:
                        pass
                        #print "%s%s:\n%sgene: %s\n%sproduct: \"%s\"" % (print_indent, kid, print_indent, str(genes), print_indent, product)
                else:
                    match = re.match(r'(M\d+)\s+(.*)', fields[1])
                    if match:
                        print '%s%s:\n%s    definition: "%s"\n%s    orthologs:' % (print_indent, match.group(1), print_indent, match.group(2), print_indent)
                    else:
                        raise RuntimeError('Line does not fit the form:\n%s' % fields[1])

        elif line[0] == 'E':
            level = 4
            fields = line.split(None, 1)
            if len(fields) == 1:
                continue
            else:
                match = re.match(r'(K\d+)\s+(.*);\s+(.*)', fields[1])
                if match:
                    genes = match.group(2).split(',')
                    genes = [x.lstrip() for x in genes]
                    product = match.group(3)
                    kid = match.group(1)
                    ec_match = re.match(r'(.*)\s+\[EC:(.*)\]$', product)
                    print_indent = ' '*level*indent
                    print "%s    - %s" % (print_indent, kid)
                    #if ec_match:
                    #    ecs = ec_match.group(2).split()
                    #    print "%s%s:\n%sgene: %s\n%sproduct: \"%s\"\n%sEC_number: %s" % (print_indent, kid, print_indent, str(genes), print_indent, ec_match.group(1), print_indent, str(ecs))
                    #else:
                    #    print "%s%s:\n%sgene: %s\n%sproduct: \"%s\"" % (print_indent, kid, print_indent, str(genes), print_indent, product)
                else:
                    raise RuntimeError('Line does not fit the form:\n%s\n%s' % (fields[1], line))
        else:
            raise RuntimeError('Line does not fit the form:\n%s', line)
