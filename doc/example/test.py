#!/usr/bin/python

f = open("tmp1","r")
t = f.read()
f.close()

sum1 = 0
sum2 = 0
sum3 = 0

for i in t.split("]"):
	tmp =[j for j in i.split("\n")[1:] if "(" in j]
	for r in tmp:
		a =r.split("(")[1].replace(")","").split(",")
		print a
		sum1+=float(a[0])
		sum2+=float(a[1])
	 	sum3+=float(a[2])

print sum1, " ", sum1/3.0
print sum2, " ", sum2/3.0
print sum3, " ", sum3/3.0






