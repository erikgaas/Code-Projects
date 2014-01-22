import re
import math
import itertools

def read_data_file(txt_file):
	"""Given txt file with data, read them into appropriate
	datatypes"""
	f = open(txt_file).read()
	f = f.split("\n")
	return [x for x in f]



def print_list_to_file(ls):
	"""Given a list, output list as carriage returns in txt file"""
	f = open("output.txt", "w")
	for i in range(len(ls)):
		if not i == len(ls) - 1:
			f.write(str(ls[i]) + " ")
		else:
			f.write(str(ls[i]))
	return


def rev(s):
	fin = ""
	for i in range(len(s)):
		fin = fin + s[len(s) - 1 - i]
	return fin

def comp(s):
	fin = ""
	for i in s:
		if i == 'A':
			fin = fin + 'T'
		elif i == 'T':
			fin = fin + 'A'
		elif i == 'C':
			fin = fin + 'G'
		elif i == 'G':
			fin = fin + 'C'
	return fin

def rev_comp(s):
	"""Find the reverse-complement of a DNA sequence"""
	return rev(comp(s))

def sub_sequence_indicies(sequence, sub_seq):
	"""Find all indicies of subseq in sequence"""
	sub_len = len(sub_seq)
	return [i for i in range(len(sequence) - sub_len + 1) if sequence[i:i+sub_len] == sub_seq]

def kmer_frequencies(sequence, k):
	"""Given DNA sequence and subseq length, return
	dict of all subseq:frequency pairs"""
	all_kmer = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
	d = {sub_seq: all_kmer.count(sub_seq) for sub_seq in all_kmer}
	return d


def kmer(sequence, k):
	"""Given a DNA sequence and a subsequence length,
	return all of the most common kmer sequences"""
	d = kmer_frequencies(sequence, k)
	max_freq = d[max(d, key=d.get)]
	return [x for x in d if d.get(x) == max_freq]


def find_all_indicies(sequence, k, L):
	"""Given a sequence and L, give all of the indicies of kmers
	from sequence[0] to sequence[2L]"""
	d = {}
	if len(sequence) > 2*L:
		frame_size = 2*L
	else:
		frame_size = len(sequence)

	for i in range(frame_size - k + 1):
		current = sequence[i:i+k]
		if current not in d:
			d[current] = sub_sequence_indicies(sequence[0:frame_size], current)
	return d



def kmer_clumps(sequence, k, L, t): #better_kmer_clumps is way better. Use that one
	"""Given a DNA sequence, find all distinct kmers that
	form (L,t)-clumps in sequence."""
	found = []
	for step in range(int(math.ceil(float(len(sequence))/float(L))) - 1):
		print step * L
		indicies = find_all_indicies(sequence[step*L:], k, L)
		for i in indicies:

			if len(indicies[i]) >= t and i not in found:
				for j in range(len(indicies[i]) - t + 1):
					if indicies[i][j+t - 1] - indicies[i][j] < L - k + 1:
						found.append(i)
						break


	return found

def better_kmer_clumps(sequence, k, L, t):
	"""Given a DNA sequence, find all distinct kmers that
	form (L,t)-clumps in sequence."""
	num_done = 0
	finished = []
	d = {}
	for i in range(len(sequence) - k + 1):
		if i%10000 == 0:
			print i
		if i - L + k > 0:
			d[sequence[i-L+k-1:i-L+k+k-1]].remove(i-L+k-1) #remove elements out of L frame

		if sequence[i:i+k] in d:
			d[sequence[i:i+k]].append(i) #add position to index
		else:
			d[sequence[i:i+k]] = [i]

		if len(d[sequence[i:i+k]]) >= t:
			finished.append(sequence[i:i+k])

	f = open("output.txt", "w")
	f.write(str(finished))
	return len(set(finished))


def skew(sequence):
	"""Given a DNA sequence, return a skew list of the number of Gs
	minus the number of Cs"""
	current_num = 0
	ls = [0]
	for i in sequence:
		if i == "C":
			current_num = current_num - 1
		elif i == "G":
			current_num = current_num + 1
		ls.append(current_num)
	return ls

def find_mins(answer):
	"""Given a skew list, return the positions of the minimum skew position"""
	current_min = 10000
	mins = []
	for i in range(len(answer)):
		if answer[i] < current_min:
			mins = [i]
			current_min = answer[i]
		elif answer[i] == current_min:
			mins.append(i)
	return mins

def approx_match(s1,s2,d):
	"""Given two strings, return whether they are d changes appart from each other"""
	number_of_mistakes = 0
	for i in range(len(s1)):
		if s1[i] != s2[i]:
			number_of_mistakes = number_of_mistakes + 1
		if number_of_mistakes > d:
			return False
	return True

def sequence_mismatch(str, sequence, d):
	"""Give index of all positions where str is at most d changes appart from a 
	subsequence of sequence"""
	ls = []
	for i in range(0, len(sequence) - len(str) + 1):
		if approx_match(str, sequence[i:i+len(str)], d):
			ls.append(i)
	return ls

def create_graph(d_freq, mismathes):
	"""Given a dictionary of frequencies of kmers, construct a graph where
	all "similar" sequences are paired with each other"""
	done = []
	graph = {}
	for i in d_freq:
		graph[i] = [d_freq[i], []]
		for j in d_freq:
			if j != i and approx_match(i,j, mismathes):
				graph[i][1].append(j)
	return graph

def all_mismatch(sequence, k, d):
	"""Given a sequence, k, and d mismatches, return the dictionary list, pairing every occurance
	with a list of every matching string"""
	log = {}
	for i in range(len(sequence) - k + 1):
		if sequence[i:i+k] not in log:
			log[sequence[i:i+k]] = map(lambda x: sequence[x:x+k] , sequence_mismatch(sequence[i:i+k], sequence, d))
	return log

def max_mismatch(sequence, k, d):
	"""Find the string within the sequence that mismatch-pairs with the greatest number
	of subsequences within the string"""
	elements = all_mismatch(sequence, k, d)
	ls = []
	current_max = 0
	for i in elements:
		if len(elements[i]) > current_max:
			current_max = len(elements[i])
			ls = [i]
		elif len(elements[i]) == current_max:
			ls.append(i)
	return ls

def consensus(ls):
	"""Find the consensus sequence of a list of sequences"""
	char = ""
	for j in range(len(ls[0])):
		
		a = 0
		g = 0
		c = 0
		t = 0
		for i in ls:
			if i[j] == "A":
				a = a + 1
			elif i[j] == "G":
				g = g + 1
			elif i[j] == "C":
				c = c + 1
			elif i[j] == "T":
				t = t + 1
		char = char + "("
		if a == max(a,g,c,t):
			char = char + "A"
		if g == max(a,g,c,t):
			char = char + "G"
		if c == max(a,g,c,t):
			char = char + "C"
		if t == max(a,g,c,t):
			char = char + "T"
		char = char + ")"
	return char

def mutations(word, hamming_distance, charset='ATCG'):
	"""function was retrieved from stack overflow"""
	for indices in itertools.combinations(range(len(word)), hamming_distance):
		for replacements in itertools.product(charset, repeat=hamming_distance):
			mutation = list(word)
			for index, replacement in zip(indices, replacements):
				mutation[index] = replacement
				yield "".join(mutation)


def find_kmer_mismatch(seed, sequence, d):
	current_max = len(sequence_mismatch(seed, sequence, d))
	ls = []
	for i in set([x for x in mutations(seed, d)]):
		test = sequence_mismatch(i, sequence, d)
		if len(test) > current_max:
			return find_kmer_mismatch(i, sequence, d)
		elif len(test) == current_max:
			ls.append(i)
	return ls

def all_max_mismatch(ls, sequence, d):
	answer = []
	current_max = 0
	for i in ls:
		temp = sequence_mismatch(i, sequence, d)
		if len(temp) > current_max:
			answer = [i]
			current_max = len(temp)
		elif len(temp) == current_max:
			answer.append(i)
	return answer

def final_mismatch(sequence, k, d):
	ls = max_mismatch(sequence, k, d)
	rest = []
	for i in ls:
		rest = rest + find_kmer_mismatch(i, sequence, d)
	return all_max_mismatch(set(rest), sequence, d)


data = read_data_file("data.txt")
e_coli = open("E-coli.txt").read()
#print_list_to_file(sub_sequence_indicies(data[0], data[1]))
#print sub_sequence_indicies("GATATATGCATATACTT", "ATAT")
# test = create_graph(kmer_frequencies("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4), 2)
# to_delete = []
# for i in test:
# 	if test[i][1] == []:
# 		to_delete.append(i)

# for i in to_delete:
# 	del test[i]

# print test

# for i in mutations("TTAAA", 2):
# 	if i == "AAAAA":
# 		print "hello"



test = final_mismatch("CGCCTACGCTGCTGCTTACCGGGCGCCGCTGCTGCTGCTGCTCGCCTAATACCGCCCGCCCGCCCGCCGCTCGGGCGCCGCTTACGCTTACCGGGCGCCTACGCTTACGCTGCTGCTTACTAATAAGCTTACTAACGCCCGCCTAAGCTTACTAACGCCTACCGGGTAATACTACTACCGCCTAACGCCCGCCGCTTAACGGGGCTCGCCGCTGCTTACTACGCTCGGGGCTCGGGCGCCTACTACTAACGCCCGCCCGGGTAACGCCCGGGTAACGCCGCTCGGGGCTTACGCTTACGCTGCTCGGGTACTAAGCTGCTCGCCTACTACTAATACCGGGTAATAATAATACCGGGTAATAATAA", 9, 3)
print_list_to_file(test)
#print max_mismatch("CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC", 10, 2)
#print find_kmer_mismatch("GCGCACACAC", "CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC", 2, 0)