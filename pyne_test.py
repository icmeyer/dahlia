from pyne import data, nucname
import numpy as np

print(data.decay_const('U-235'))
print(data.decay_const('922350'))
print(data.decay_const('922350000'))
print(data.branch_ratio('932390000','942390000', use_metastable=False))
print(data.decay_children('932390000'))

print(data.decay_const('942420000'))
print(data.decay_children('942420000'))
print(np.log(2)/data.decay_const('922340000')/3.15e7)

print('-------Np-239 to Pu-239 test--------')
print(data.decay_const('932390'))
print(data.decay_children('932390'))

print(data.decay_const('932390'))
print(data.decay_children('932390'))
print(data.branch_ratio(932390,942390))

print('-------U-240 decay test--------')
print(np.log(2)/data.decay_const('922400')/3600)
print(data.branch_ratio('922400','932400', use_metastable=False))

print('-----U234 Capture Test-----')
print(float('922350')-float('922340') == 10)

print('-----Mass Test-----')
print(nucname.anum('922350'))
print('-----Name Test-----')
print(nucname.serpent('922350'))

print('-------U-236 test--------')
print(data.decay_const('922360'))
print(data.half_life('922360')/(3.15e7))

print('-------Np-234 to U-234 test--------')
print(data.decay_const('932340'))
print(data.decay_children('932340'))
print(data.branch_ratio(932340,922340))
