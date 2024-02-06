#!/bin/env python

import argparse
import requests
from bs4 import BeautifulSoup
import re
import sys
import os

def get_browsershot(server_url, session_file, position_string, ofname):
	if "genometest.gs.washington.edu" in server_url:
		hgsid = "?hgsid=100000"
	else: 
		hgsid = ""
	url = "".join([server_url, "/cgi-bin/hgTracks", hgsid, "?hgS_doLoadUrl=submit&hgS_loadUrlName=", session_file, "&hgt.psOutput=on&pix=2000", position_string])
	page = requests.get(url)
	print(url)
	if page.status_code != requests.codes.ok:
		print ("Invalid page URL: %s\n" % url)
		print ("Make sure your session file is globally readable and in a web-accessible directory\n")
		sys.exit(1)

	soup = BeautifulSoup(page.text, "html.parser")
	relative_url = None
	for entry in soup.find_all(href=re.compile("pdf")):
		if entry.parent.find(text=re.compile("the current browser graphic in PDF")) is not None:
			relative_url = entry.get("href")
			break

	if relative_url is None:
		print ("Could not find browsershot pdf at %s\n" % url)
		print ("Make sure your session file is globally readable and in a web-accessible directory\n")
		sys.exit(1)

	pdf_url = server_url + relative_url.replace("../", "/")
	print (ofname)
	with open(ofname, "wb") as outfile:
		r = requests.get(pdf_url)
		if r.status_code == requests.codes.ok:
			outfile.write(r.content)
		else:
			print (r.headers)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
#	group = parser.add_mutually_exclusive_group(required=True)
	parser.add_argument("--regions_file", type = argparse.FileType("r"), help="Gets browsershots for multiple regions. Must be a path to tab-delimited file with list of region, outfile name, and session name (chr,start,end,ofname,sessionname).")
	parser.add_argument("--region", default = None, nargs=5, metavar=("chr", "start", "end", "outfile", "session_file"), help="Get browsershot for single region.")
	parser.add_argument("--region_highlight", default = None, nargs=3, metavar=("chr", "start", "end"), help="Highlight a region in the browsershot.")
	parser.add_argument("--tempSessionFile", default = None, help="path to a tempSessionFile.")
	parser.add_argument("--server_url", default="https://genome.gs.washington.edu", help="Server URL (May require username and password)")

	args = parser.parse_args()

	if args.regions_file is not None:
		regions = []
		for line in args.regions_file:
			regions.append(line.rstrip().split())

		for region in regions:
			position_string = "&position=%s:%s-%s" % (region[0], region[1], region[2])
			outfile = region[3]
			session_file = os.path.abspath(region[4])
			get_browsershot(args.server_url, session_file, position_string, outfile)
	else:
		position = "%s:%s-%s" % tuple(args.region[0:3])
		outfile = args.region[3]
		session_file = os.path.abspath(args.region[4])
		position_string = "&position=" + position

		if args.region_highlight is not None and args.tempSessionFile is not None:
			region_highlight = "%s:%s-%s" % tuple(args.region_highlight) + "#ffeda0"
			
			with open(session_file) as fin , open(args.tempSessionFile, "w") as fout:
				for line in fin:
					if line.startswith("db"):
						db = line.strip().split()[1]
						region_highlight = db + "." + region_highlight
						fout.write(line)
					elif line.startswith("highlight"):
						fout.write("highlight " + region_highlight + "\n")
					else:
						fout.write(line)
			
			session_file = args.tempSessionFile
		
		get_browsershot(args.server_url, session_file, position_string, outfile)
		
