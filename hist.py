import bisect
import csv

import matplotlib.pyplot as plt


def pixel_to_value(y):
    if y == 0:
        return 0
    return (1290 - y) / 75.2


def accumu(l):
    total = 0
    for x in l:
        total += x
        yield total

def estimate_percentile(value, vals_cum):
    i = bisect.bisect_right(val, value)
    ratio = (value - val[i - 1]) / (val[i] - val[i - 1])
    return vals_cum[i - 1] + ratio * (vals_cum[i] - vals_cum[i - 1])

val = []
t1d = []
t2d = []
with open("data/dist.csv") as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        val.append(float(row[0]))
        t2d.append(pixel_to_value(float(row[1])))
        t1d.append(pixel_to_value(float(row[2])))

t1d_cum = [x * 100 /sum(t1d) for x in list(accumu(t1d))]
t2d_cum = [x * 100 /sum(t2d) for x in list(accumu(t2d))]

cutoff = 0.231
print("T1D percentile for %.3f estimated %.0f" % (cutoff, estimate_percentile(cutoff, t1d_cum)))
print("T2D percentile for %.3f estimated %.0f" % (cutoff, estimate_percentile(cutoff, t2d_cum)))

individual = 0.266
print("T1D percentile for %.3f estimated %.0f" % (individual, estimate_percentile(individual, t1d_cum)))
print("T2D percentile for %.3f estimated %.0f" % (individual, estimate_percentile(individual, t2d_cum)))

bins = val[:]
bins.append(bins[-1] + 0.005)

plt.hist(val, bins, weights=t2d, alpha=1.0, label='Type 2 diabetes', edgecolor='black', color='#9DD0D9')
plt.hist(val, bins, weights=t1d, alpha=0.7, label='Type 1 diabetes', edgecolor='black', color='#F0B298')
plt.axvline(x=cutoff, color='red', linestyle='--', label='Cutoff')
plt.axvline(x=individual, color='blue', linestyle='--', label='Individual')
plt.xlabel('Type 1 diabetes genetic risk score')
plt.ylabel('Population density')
plt.yticks([0, 5, 10, 15])
plt.legend(loc='upper left')
plt.savefig('t1d-grs.png')
