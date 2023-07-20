#!/usr/bin/python3

import argparse, re


parser = argparse.ArgumentParser()
parser.add_argument("filename")

args = parser.parse_args()
filename = args.filename

post = open(args.filename, "r").read()

post_out = open(args.filename + ".out", "w")

for match in re.findall("\$[^$]*?\$", post):
    print(match)

#print(post)
print(args)
