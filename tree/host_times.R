# Gamma for Scotti times
# JS 15/05/17
# generate random start host times using same params as in demo_*.R

# accordint to the R documentation, mean is a*s, that is shape*scale.
# try shape=1.45 scale=1/0.3 so mean about 5 (years), as used for sampling dist.

set.seed(2310)

vegard_offsets = rgamma(13,shape=1.45,scale=1/0.3)
roetzer_offsets = rgamma(86,shape=1.45,scale=1/0.3)

print(vegard_offsets)
print(roetzer_offsets)

print(mean(vegard_offsets))
print(mean(roetzer_offsets))
