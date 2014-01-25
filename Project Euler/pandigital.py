import itertools

perms = itertools.permutations("0123456789", 10)

print [x for x in perms if x[0] != "0"]