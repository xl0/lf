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

def sum_pub_scores(pub_list):
	score = 0
	for pub in pub_list:
		score += pub['Score']

	return score

def fact():
	return {
		'PubList' : [],
		'Score' : 0
	}

def merge_labs(pmid_list, medline_extract):
	cached_extract = medline_extract.distinct('_id')
	need_to_extract = list(set(pmid_list) - set(cached_extract))

	assert(len(need_to_extract) == 0)

	write('Merging labs from %d publications...\n' % len(pmid_list))
	num_labs = 0
	merged_lab_dict = defaultdict(dict)
	merged_lab_dict.default_factory = fact
	for extract in medline_extract.find({'_id' : {'$in' : pmid_list}}):
		# Not interested in books
		if not extract['article']:
			continue
	
		for author_lab in extract['AuthorLabList']:
			num_labs += 1
			# Use the first chunk of address as key
			key = author_lab['ProcessedName'][0]
			pub_score = score_date(extract['PubDate'])
			merged_entry = {'LabName' : tuple(author_lab['ProcessedName']),
				'PMID' : extract['_id'],
				'Score' : pub_score,
				'Title' : extract['Title'],
				'Author' : author_lab['Author'],
				'Authors' : author_lab['Authors']}

			merged_lab_dict[key]['PubList'].append(merged_entry)
			merged_lab_dict[key]['Score'] += pub_score
	write('Found %d labs\n' % num_labs)
	return merged_lab_dict	


def show_labs(merged_lab_dict, sorted_list):
	write('<ol>\n')
	for lab_cluster in sorted_list:
		split_lab_dict = defaultdict(list)
		pub_list = merged_lab_dict[lab_cluster]['PubList']

		for pub in pub_list:
			entry = {
				'PMID' : pub['PMID'],
				'Title' : pub['Title'],
				'Author' : pub['Author'],
				'Authors' : pub['Authors']
			}
			split_lab_dict[pub['LabName']].append(entry)
		
		write('<li>' + lab_cluster + '\n')
		write('<ul>')
		for lab_name in split_lab_dict.keys():
			write('<li>%s<ul>\n' % ', '.join(lab_name))
			for publication in split_lab_dict[lab_name]:
				pub_string = ''
				for author in publication['Authors']:
					if author == publication['Author']:
						pub_string += "<b>" + author + " </b>"
					else:
						pub_string += author + ' '
				pub_string += '<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/%s\">' % publication['PMID']
				pub_string += publication['Title'] + '</a>\n'
				write('<li>' + pub_string.encode('utf-8') + '</li>')
			write('</ul></li>')
		write('</ul></li>')

	write('</ol>')

#			print ', '.join(lab_name)
#			print split_lab_dict[lab_name]

	

#def split_labs(merged_lab_dict):
#	split_lab_dict = defaultdict(dict)
#	for key in merged_lab_dict.keys():
#		labs = merged_lab_dict[key]
#		write('key ' + str(key) + str(labs) + '\n')
#		variants = defaultdict(list)
#		for lab in labs:
#			new_entry = {
#				'pubmed' : lab['pubmed']
#			}



#			print lab
#			split_la_entry = {'pubmed' : lab}
#			split_lab_dict[key].append()
		


def fetch_me_data(record_list, raw_medline_col):
	
	cached_medline = raw_medline_col.distinct('_id')
	need_to_fetch = list(set(record_list) - set(cached_medline))

	if len(need_to_fetch) > 0:
		num_records = len(need_to_fetch)
		write("I need to consult the literature (%d of %d articles new to me)...\n" % (num_records, len(record_list)))
		string_pmid_list = [ str(x) for x in need_to_fetch]

		n = 0
		while True:
			try:
				handle = Entrez.epost(db='pubmed', id=','.join(string_pmid_list))
			except urllib2.HTTPError:
				n += 1
				if n > 3:
					write('NCBI is not happy with us\n')
					return False
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
						return False
					write('Ugh\n')
					continue
				break

			data_list = Entrez.read(handle)
			handle.close()
			del handle
			write('Done. Memorizing them...')

			for entry in data_list:
				if entry.has_key('MedlineCitation'):
					pmid = int(entry['MedlineCitation']['PMID'])
					is_article = True
				elif entry.has_key('BookDocument'):
					pmid = int(entry['BookDocument']['PMID'])
					is_article = False
				else:
					raise KeyError

				raw_medline_col.save({'_id' : pmid, 'entry' : entry, 'is_article': is_article})

			del data_list
			write("Done!\n")
	else:
		write('Nothing to fetch, all in cache\n')
	return True

def get_medline_author_list(entry):
	author_list = []
	if entry.has_key('MedlineCitation'):
		citation = entry['MedlineCitation']
		if citation.has_key('Article'):
			article = citation['Article']
			if article.has_key('AuthorList'):
				authors = article['AuthorList']
				for author in authors:
					if author.has_key('Affiliation'):
						lab = author['Affiliation']
						forename = author.get('ForeName', '')
						last_name = author.get('LastName', '')
#						write('author = ' + forename + last_name + '\n')
			
						author_list.append((forename + ' ' + last_name, lab))
	return author_list

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

def get_medline_date(entry):
	citation = entry['MedlineCitation']
	day = int(citation['DateCreated']['Day'])
	try:
		month = int(citation['DateCreated']['Month'])
	except ValueError:
		month = month_to_int(citation['DateCreated']['Month'])

	year = int(citation['DateCreated']['Year'])

	return (year, month, day)

def get_medline_title(entry):
	citation = entry['MedlineCitation']
	article = article = citation['Article']
	title = article['ArticleTitle']
	return title

def score_citation(entry):
	# For now, just score based on publication date

	date = get_citation_date(entry)
	score = score_date(date)

	return score


replacement_dict = [
	('[De]ept\.?', 'Department'),
	('[Ll]ab\.?\s', 'Laboratory '),
	('\s&\s', ' and '),
	('Electronic address', '')
]

def heuristic_fix(string, dictionary):
	string = string.strip()
	for entry in dictionary:
		n_string = re.sub(entry[0], entry[1], string)
#		if not n_string == string:
#			print string
#			print n_string
		string = n_string.strip()
	return string

def extract_data(pmid_list, raw_medline, medline_extract):
	cached_extract = medline_extract.distinct('_id')
	need_to_extract = list(set(pmid_list) - set(cached_extract))

	if len(need_to_extract) > 0:
		write('Need to analyze %d of %d articles...\n' % (len(need_to_extract), len(pmid_list)))
		fetch_me_data(need_to_extract, raw_medline)

		n = 0
		write('need to extract: ' +  str(need_to_extract) + '\n')
		for data in raw_medline.find({'_id' : {'$in' : need_to_extract}}):
			if not n % 1000:
				write(str(n) + '\n')
			n += 1
			entry = data['entry']
			# If it's a book, just save the id to not look at it again.
			if not data['is_article']:
				entry = {'_id' : data['_id'],
					'article' : False}
				medline_extract.save(entry)
				continue

			entry_list = []
			author_list = get_medline_author_list(entry)
			# Full author list is saved for each lab entry's publications
			full_author_list = []
			for author in author_list:
				full_author_list.append(author[0])

			for author in author_list:
				lab_name = sanitize_lab_name(author[1])
				# Determine the separator they use
				comas = lab_name.count(',')
				semicol = lab_name.count(';')
				if comas >= semicol:
					tmp_lab_name_array = lab_name.split(',')
				else:
					tmp_lab_name_array = lab_name.split(';')

				lab_name_array = []
				# Now fix up the names - clean up and do replacements
				for tmp in tmp_lab_name_array:
					tmp = tmp.strip()
					tmp = heuristic_fix(tmp, replacement_dict)
					tmp = tmp.strip()
					lab_name_array.append(tmp)


				new_entry = { 'ProcessedName' : lab_name_array,
						'RawName' : author[1],
						'Author' : author[0],
						'Authors' : full_author_list
				}
				entry_list.append(new_entry)
			pubdate = get_medline_date(entry)
			title = get_medline_title(entry)
			extract_entry = {'_id' : data['_id'],
					'PubDate' : pubdate,
					'Title' : title,
					'AuthorLabList' : entry_list,
					'article' : True}
			medline_extract.save(extract_entry)
		write('Done with premilinary analysis\n')
	else:
		write('All %d abstracts analyzed before\n' % len(pmid_list))

	return True

def run_search(request, request_cache):

	pmid_list = []
	res = request_cache.find_one(request)
	if res:
		pmid_list = res['pmid_list']
		num_records = len(pmid_list)
		write('Request found in cache.\n')
		write('%d PMIDs to consider\n' % num_records)
	else:
		write('Let me ask NCBI...\n')

		# Do the request
		n = 0
		while True:
			try:
				handle = Entrez.esearch(db="pubmed", term=request, retmax=0, usehistory='y')
			except urllib2.HTTPError:
				n += 1
				if (n > 3):
					write('\nNCBI is not happy, giving up.\n')
					return None

				write('Request to NCBI failed, trying again...\n')
				continue
			break


		records = Entrez.read(handle)
		write(str(records) + '\n')
		num_records = int(records['Count'])
		webenv = records['WebEnv']
		query_key = records['QueryKey']

		write("Done. %d articles to consider.\n" % num_records)
		write("Getting the PMIDs...\n")
		n = 0
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
						return None
					write('Oops\n')
					continue
				break
		pmid_list = [ int(pmid) for pmid in pmid_list]
		request_cache.save({'_id' : request, 'pmid_list' : pmid_list})

	return pmid_list

def header(request):
	write("Content-Type: text/html\n\n")

	h = open('index.html')
	index = h.read()
	index = re.sub("name=q value=\"\"", "name=q value=\"%s\"" % request, index)

	write(index)

	write("<html><body>\n")
	write("<p>Meditating on the request...</p>\n")
	write("<pre>\n")
	
def footer_good():
	write("<p>Done!</p>\n")
	write("</body></html>\n")
	
def footer_bad():
	write("</pre>\n")
	write("<p>Giving up. Try again later.</p>\n")
	write("</body></html>\n")

def main():
	from pymongo import MongoClient

	mongo = MongoClient()
	db = mongo['labfinder']

	request_cache = db['request_cache']
	raw_medline = db['raw_medline']
	medline_extract = db['medline_extract']

	form = cgi.FieldStorage() # instantiate only once!
	request = form.getfirst('q', None)

	# Avoid script injection escaping the user input
	request = cgi.escape(request)

	request_escaped_quotes = cgi.escape(request, quote=True)
	header(request_escaped_quotes)

	pmid_list = run_search(request, request_cache)
	if not pmid_list:
		footer_bad()
		return
	del request_cache

	res = extract_data(pmid_list, raw_medline, medline_extract)
	if not res:
		footer_bad()
		return

	merged_lab_dict = merge_labs(pmid_list, medline_extract)
#	for i in merged_lab_dict.keys():
#		write(i + json.dumps(merged_lab_dict[i], sort_keys=True, indent=2) + '\n')


	sorted_lab_list = sorted(merged_lab_dict, key=lambda lab_key: merged_lab_dict[lab_key]['Score'], reverse=True)

	write('</pre>')
	show_labs(merged_lab_dict, sorted_lab_list)

#	for i in sorted_lab_list:
#		print merged_lab_dict[i]['Score'], i
#	order_list = 


	del pmid_list
	del medline_extract

#	split_lab_dict = split_labs(merged_lab_dict)
	footer_good()



#	data = fetch_me_data(pmid_list, raw_medline)


#	pmid_list = pmid_list_get()

'''
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
'''
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

