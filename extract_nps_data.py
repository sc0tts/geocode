import numpy as np

data = np.fromfile('../../Desktop/nt_20180101_f18_nrt_n.bin',
                   dtype=np.uint8)[300:].reshape((448, 304))

extract = data[177:185+1, 44:53+1]
print(extract)
print(' ')
print(repr(extract))
#data[177:185+1, 44:53+1].tofile('nps_data')
