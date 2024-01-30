m = [1, 2, 3, 4, 5]
cp = [11, 12, 13, 14, 15]
R = 8.31446261815324 

num = 0
den = 0

for i in range (5):
    num = num + m[i] * cp[i]
    den = den + m[i]

cpmix = num/den

gamma = (cpmix) / (cpmix - R)