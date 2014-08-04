#!/usr/bin/env python2.7

import os
import string
from Bio import Entrez
Entrez.email = 'alexey.zaytsev@gmail.com'

from collections import defaultdict

import json
import unicodedata
import re
import datetime
import time


article_half_life = 365 * 2
cachedir = "./cache/"
def mkfilename(pmid):
	return cachedir + pmid + ".json"

rec = re.compile('[\w\.]+@[\w\.]$')

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
	while not str.isalpha(lab[m]) and m < len(lab):
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
			"Wtf? Medline, but not Article?"
#		print 'An article!'
#	elif entry.has_key('BookDocument'):
#		print "A book!"

	return lab_list

def get_citation_date(entry):
	citation = entry['MedlineCitation']
	day = int(citation['DateCreated']['Day'])
	month = int(citation['DateCreated']['Month'])
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
		print "To fetch: ", need_to_fetch
		handle = Entrez.efetch(db="pubmed", id=need_to_fetch, rettype="summary", version="2.0")
		data_list = Entrez.read(handle)
		handle.close()


		for entry in data_list:
			if entry.has_key('MedlineCitation'):
				pmid = entry['MedlineCitation']['PMID']
			elif entry.has_key('BookDocument'):
				pmid = entry['BookDocument']['PMID']
			else:
				raise KeyError

			print "Saving entry ", pmid
			filename = cachedir + "/" + pmid + ".json"
			assert(not os.path.exists(filename))
			handle = open(filename, 'w+')

			handle.write(json.dumps(entry, sort_keys=True, indent=4))
			handle.close()



replacement_dict = [
	('[De]ept\.?', 'Department'),
	('[Ll]ab\.?\s', 'Laboratory '),
	('\s&\s', ' and ')
]

def heuristic_fix(string, dictionary):
	for entry in dictionary:
		n_string = re.sub(entry[0], entry[1], string)
		if not n_string == string:
			print string
			print n_string
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




def pmid_list_get():
	h = open('files.txt')
	for line in h.readlines():
#		print line
		pmid = line.split('.')
#		print pmid[0]
		yield pmid[0]

#handle = Entrez.esearch(db="pubmed", term="\"synthetic biology\"", retmax=100000)
#records = Entrez.read(handle)
#pmid_list = records['IdList']
#data = fetch_me_data(pmid_list)


pmid_list = pmid_list_get()


labs_dict = defaultdict(list)


for pmid in pmid_list:
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


merged_labs_dict = merge_labs(labs_dict)


#print json.dumps(merged_labs_dict, sort_keys=True, indent=4)

#print labs_dict
sorted_lab_list = sorted(merged_labs_dict, key=lambda entry: merged_labs_dict[entry]['Score'])

for lab in sorted_lab_list:
	print merged_labs_dict[lab]['Score'], lab
#	print json.dumps(labs_dict[lab], sort_keys=True, indent=4)


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

