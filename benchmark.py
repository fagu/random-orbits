#!/usr/bin/python3

import os
import sys
import time

ps = [[2**10, int(pow(2, 0.25*k))] for k in range(3*4,10*4)]
ps += [[2**min(10,23-k), 2**k] for k in range(10,19)]

for [a, b] in ps:
	print(b, file=sys.stderr)
	# print(a, b)
	start = time.time()
	os.system(f"/usr/bin/python3 -c 'import sys; sys.set_int_max_str_digits(0); print(2**{b})' | ./random-cubic {a} 3 - > /dev/null")
	end = time.time()
	# print(end - start)
	print(b, (end - start) / a)
	sys.stdout.flush()
