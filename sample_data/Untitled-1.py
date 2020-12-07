def factorize(number):
	if number/5 == 1 or number/2 ==1:
		return True
	elif number % 5 != 0 and number %2 != 0:
		return False
	elif number % 5 == 0:
		return factorize(number/5)
	else: 
		return factorize(number/2)

print(factorize(9))