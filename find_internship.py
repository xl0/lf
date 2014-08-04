#!/usr/bin/env python2.7

import os
import sys
import string
from Bio import Entrez
Entrez.email = 'alexey.zaytsev@gmail.com'
import urllib2

from collections import defaultdict

import json
import unicodedata
import re
import datetime
import time
import gc
import cgi
import cgitb; cgitb.enable()

article_half_life = 365 * 2

global chunk_size
chunk_size = 5000
global search_chunk_size
search_chunk_size = 100000
global pmid_cache_file
pmid_cache_file = 'pmid_cache.json'

cachedir = "./cache/"
def mkfilename(pmid):
	return cachedir + pmid + ".json"

rec = re.compile('[\w\.]+@[\w\.]$')

def write(string):
	sys.stdout.write(string)
#	sys.stdout.write('')
	sys.stdout.flush()

# Some lab entries start with mis-parsed garbage, remove it.
# Also, remove accents/umlauts/etc.
# Also, remove email, if found in the end.
def sanitize_lab_name(lab):
	# First - de-accent it.
	lab = unicodedata.normalize('NFD', lab).encode('ascii', 'ignore')

	orig = lab

	# Strip non-alphas before the name. We loose a few bits, but get
	# rid of a lot of crap.
	m = 0
	while m < len(lab) and not str.isalpha(lab[m]):
		m += 1
	lab = lab[m:]

	while True:
		done = True
		# Remove crap from the end as well.
		n = len(lab) - 1
		while n >= 0 and not (str.isalpha(lab[n])):
			n -= 1

		lab = lab[:n+1]
		res = re.search('[\w\.\-]+@[\w\.\-]+$', lab)
		if not res == None:
			lab = lab[:res.start()]
		else:
			# Break if no email, otherwise, look for an other.
			break

#	print "orig:", orig
#	print "xxxx:", lab

	return lab

def score_date(ymd_tuple):
	date = datetime.date(ymd_tuple[0], ymd_tuple[1], ymd_tuple[2])
	local_time = time.localtime()
	local_date = datetime.date(local_time.tm_year, local_time.tm_mon, local_time.tm_mday)

	delta_d = local_date - date

	score = float(article_half_life) / (article_half_life + delta_d.days)

	return score

def load_pubmed_data(pmid):
#	print "Meditating on ", pmid
	filename = mkfilename(pmid)
	assert(os.path.exists(filename))

	handle = open(filename)
	entry = json.load(handle)
	handle.close

	return entry

def get_lab_list(entry):
	lab_list = []
	if entry.has_key('MedlineCitation'):
		citation = entry['MedlineCitation']
#		citation_date = get_article_date(citation)
		if citation.has_key('Article'):
			article = citation['Article']
#			print pmid, entry['MedlineCitation']['Article']['ArticleDate']
			if article.has_key('AuthorList'):
				author_list = article['AuthorList']

				for author in author_list:
					if author.has_key('Affiliation'):
						lab = author['Affiliation']
#						lab = sanitize_lab_name(author['Affiliation'])
						lab_list.append(lab)

#				print "========="
#				print json.dumps(author_list, sort_keys=True, indent=4)

				pass
#				print "Looks good, ", len(entry['MedlineCitation']['Article']['AuthorList'])
		else:
			pass
#			print "Wtf? Medline, but not Article?"
#		print 'An article!'
#	elif entry.has_key('BookDocument'):
#		print "A book!"

	return lab_list

def month_to_int(month):
	return{
		'Jan' : 1,
		'Feb' : 2,
		'Mar' : 3,
		'Apr' : 4,
		'May' : 5,
		'Jun' : 6,
		'Jul' : 7,
		'Aug' : 8,
		'Sep' : 9, 
		'Oct' : 10,
		'Nov' : 11,
		'Dec' : 12
	}[month]

def get_citation_date(entry):
	citation = entry['MedlineCitation']
	day = int(citation['DateCreated']['Day'])
	try:
		month = int(citation['DateCreated']['Month'])
	except ValueError:
		month = month_to_int(citation['DateCreated']['Month'])

	year = int(citation['DateCreated']['Year'])

	return (year, month, day)


def score_citation(entry):
	# For now, just score based on publication date

	date = get_citation_date(entry)
	score = score_date(date)

	return score

def fetch_me_data(record_list):
	need_to_fetch = []
	for pmid in record_list:
		filename = cachedir + "/" + pmid + ".json"
		if not os.path.isfile(filename):
			need_to_fetch.append(pmid)

	if len(need_to_fetch) > 0:
		num_records = len(need_to_fetch)
		write("I need to consult the literature (%d of %d articles new to me)...\n" % (num_records, len(record_list)))
		n = 0
		while True:
			try:
				handle = Entrez.epost(db='pubmed', id=','.join(need_to_fetch))
			except urllib2.HTTPError:
				n += 1
				if n > 3:
					write('NCBI is not happy with us\n')
					False
				write('Ugh\n')
				continue
			break

		result = Entrez.read(handle)
		handle.close()
		del handle
		webenv = result['WebEnv']
		query_key = result['QueryKey']

		for start in range(0, num_records, chunk_size):
			end = min(start+chunk_size, num_records)
			write('Getting abstracts %d-%d ... ' % (start, end))
			n = 0
			while True:
				try:
					handle = Entrez.efetch(db="pubmed", rettype="summary", version="2.0",
							retstart=start, retmax=chunk_size, webenv=webenv, query_key=query_key)
				except urllib2.HTTPError:
					n += 1
					if n > 3:
						write('NCBI is not happy with us\n')
						False
					write('Ugh\n')
					continue
				break

			data_list = Entrez.read(handle)
			handle.close()
			del handle
			write('Done. Memorizing them...')

			for entry in data_list:
				if entry.has_key('MedlineCitation'):
					pmid = entry['MedlineCitation']['PMID']
				elif entry.has_key('BookDocument'):
					pmid = entry['BookDocument']['PMID']
				else:
					raise KeyError

				filename = cachedir + "/" + pmid + ".json"
				assert(not os.path.exists(filename))
				handle = open(filename, 'w+')
				json_string = json.dumps(entry, sort_keys=True, indent=4)
				handle.write(json_string)
				del json_string
				handle.close()
				del handle
			del data_list
			write("Done!\n")
	else:
		write('And I\'ve already read all of them!\n')
		

replacement_dict = [
	('[De]ept\.?', 'Department'),
	('[Ll]ab\.?\s', 'Laboratory '),
	('\s&\s', ' and ')
]

def heuristic_fix(string, dictionary):
	for entry in dictionary:
		n_string = re.sub(entry[0], entry[1], string)
#		if not n_string == string:
#			print string
#			print n_string
		string = n_string
	return string

def fact():
	return {'Variants' : [],
		'Score' : 0
		}

def sum_pub_scores(pub_list):
	score = 0
	for pub in pub_list:
		score += pub['Score']

	return score

def merge_labs(lab_dict):
	merged_dict = defaultdict(list)
	merged_dict.default_factory = fact

	for lab_key in lab_dict.keys():
		lab_name = sanitize_lab_name(lab_key)
		# Determine the separator they use
		comas = lab_name.count(',')
		semicol = lab_name.count(';')
#		print lab_key
#		print comas, semicol
		if comas >= semicol:
			tmp_lab_name_array = lab_name.split(',')
		else:
			tmp_lab_name_array = lab_name.split(';')

		lab_name_array = []
		# Now fix up the names - clean up and do replacements
		for tmp in tmp_lab_name_array:
			tmp = tmp.strip()
			tmp = heuristic_fix(tmp, replacement_dict)
			lab_name_array.append(tmp)


		new_entry = { 'NameVariant' : lab_name_array[1:],
				'RawName' : lab_key,
				'Publications' : lab_dict[lab_key],
				'Score' : sum_pub_scores(lab_dict[lab_key])
				}
		merged_dict[lab_name_array[0]]['Variants'].append(new_entry)
		merged_dict[lab_name_array[0]]['Score'] += new_entry['Score']

	return merged_dict

#	for lab_key in merged_dict.keys():
#		print lab_key
#	print json.dumps(merged_dict, sort_keys=True, indent=4)

#		print lab_key, lab_dict[lab_key]

def main():
	form = cgi.FieldStorage() # instantiate only once!
	request = form.getfirst('q', None)

	# Avoid script injection escaping the user input
	request = cgi.escape(request)

	pmid_cache = {}
	if os.path.exists(pmid_cache_file):
		h = open(pmid_cache_file)
		pmid_cache = json.load(h)
		h.close()
		del h

	write("Content-Type: text/html\n\n")

	h = open('index.html')
	index = h.read()
	index = re.sub("name=q value=\"\"", "name=q value=\"%s\"" % request, index)

	write(index)


	write("<html><body>\n")
	write("<p>Meditating on the request...</p>\n")
	write("<pre>\n")


	if pmid_cache.has_key(request):
		write('Query found in cache.\n')
		pmid_list = pmid_cache[request]
	else:		
		write('Let me ask NCBI...')

		# Do the request
		n = 0
		while True:
			try:
				handle = Entrez.esearch(db="pubmed", term=request, retmax=0, usehistory='y')
			except urllib2.HTTPError:
				n += 1
				if (n > 3):
					write('NCBI is not happy, giving up')
					return 1
					
				write('Oops\n')
				continue
			break

			
		records = Entrez.read(handle)
		num_records = int(records['Count'])
		webenv = records['WebEnv']
		query_key = records['QueryKey']

		write("Done. %d articles to consider.\n" % num_records)
		write('webenv=%s, query_key=%s\n' % (webenv, str(query_key)))
		write("Getting the PMIDs...\n")
		n = 0
		pmid_list = []
		for start in range(0, num_records, search_chunk_size):
			end = min(start+search_chunk_size, num_records)
			while True:
				write('Getting %d-%d\n' % (start, end))
				try:
					handle = Entrez.esearch(db="pubmed", retstart=start, retmax=search_chunk_size,
						 webenv=webenv, query_key=query_key)
					records = Entrez.read(handle)
					new_pmid_list = records['IdList']
					
					pmid_list += new_pmid_list
				except urllib2.HTTPError:
					n += 1
					if (n > 3):
						write('NCBI is not happy, giving up')
						return 1
					write('Oops\n')
					continue
				break
			pmid_cache[request] = pmid_list
			h = open(pmid_cache_file + '.tmp', 'wr+')
			json.dump(pmid_cache, h, sort_keys=True, indent=4)
			h.close()
			os.rename(pmid_cache_file + '.tmp', pmid_cache_file)

	data = fetch_me_data(pmid_list)


#	pmid_list = pmid_list_get()


	labs_dict = defaultdict(list)

	write('1\n')

	n = 0
	for pmid in pmid_list:
		n += 1
		if not n % 10000:
			write(str(n) + '\n')
		pubmed_data = load_pubmed_data(pmid)
		labs = get_lab_list(pubmed_data)
		if len(labs) > 0:
			pubdate = get_citation_date(pubmed_data)
			score = score_citation(pubmed_data)
		for lab in labs:
			labs_dict[lab].append({ 'PMID' : pmid,
						'PibDate' : pubdate,
						'Score' : score
						})


	write('2\n')
	merged_labs_dict = merge_labs(labs_dict)


	write('3\n')
	#print json.dumps(merged_labs_dict, sort_keys=True, indent=4)

	#print labs_dict
	sorted_lab_list = sorted(merged_labs_dict, key=lambda entry: merged_labs_dict[entry]['Score'])

	write('4\n')
	for lab in sorted_lab_list:
		print merged_labs_dict[lab]['Score'], lab
	#	print json.dumps(labs_dict[lab], sort_keys=True, indent=4)

	print "</pre>"
	print "<p>Done!</p>"
	print "</body></html>"

if __name__ == "__main__":
    main()


#def pmid_list_get():
#	h = open('files.txt')
#	for line in h.readlines():
#		pmid = line.split('.')
#		print pmid[0]
#		yield pmid[0]



#print sorted_lab_list

#print labs_dict

#handle = Entrez.efetch(db="pubmed", id="15491125", rettype="summary")
#print handle.read()
#record = Entrez.read(handle)
#print json.dumps(record, sort_keys=True, indent=4)
#handle.close()

#print record


#handle = Entrez.efetch(db="pubmed", id="25063984", rettype="summary", version="2.0")
#entry = Entrez.read(handle)

#affiliations = get_affiliations(entry)

#record = Entrez.read(handle)

