import sys
import numpy as np
n=120;k=2*np.pi*np.fft.fftfreq(n);w=2.6*abs(np.sin(k/2));sgn=np.sign(k);sgn[60]=0
j=np.arange(n)-18
carrier=.9594572754 if "--printed" in sys.argv else 2*np.arcsin(6/13)
q=.002*np.exp(-j*j/98)*np.cos(carrier*j)
v=np.fft.ifft(-1j*sgn*w*np.fft.fft(q)).real
print(' '.join(format(float(z),'.17g') for z in np.r_[q,v]))
