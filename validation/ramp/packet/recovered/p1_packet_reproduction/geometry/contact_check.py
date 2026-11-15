"""Independent pair-force and energy/flux test for Paper One's six-segment geometry."""
import numpy as np

gen=np.random.default_rng(1866)
q=gen.normal(size=7); v=gen.normal(size=7)
u=gen.normal(size=(6,4,3)); w=gen.normal(size=(6,4,3))
G=(.1,.2,.2,.2,.2,.1)
corners=[(sy,sz) for sy in (-1,1) for sz in (-1,1)]
force_q=np.zeros(7); force_u=np.zeros_like(u); printed_q=np.zeros(7);printed_u=np.zeros_like(u)
H=H_reduced=J=dH=0.
ex=np.array([1.,0,0])
for i,g in enumerate(G):
    for c,(sy,sz) in enumerate(corners):
        tc=np.array([0.,sy,sz])
        A=u[i,c,0]-(q[i]+q[i+1])/2
        B=tc@u[i,c]+(q[i+1]-q[i])/2
        H_reduced+=g/3*(A*A+B*B)
        printed_q[i]+=g/3*(A+B); printed_q[i+1]+=g/3*(A-B)
        printed_u[i,c]=-2*g/3*(A*ex+B*tc)
        for end in range(2):
            n=(tc+(1-2*end)*ex)/np.sqrt(3); e=n@(u[i,c]-q[i+end]*ex)
            H+=g*e*e/2; actual=g*e*n
            force_u[i,c]-=actual;force_q[i+end]+=actual@ex
            J-=(w[i,c]+v[i+end]*ex)@actual/2
            dH+=g*e*(w[i,c]-v[i+end]*ex)@n
Pb=np.sum(w*force_u);Pq=v@force_q
print('source_pairs',6*4*2,'energy_agreement',abs(H-H_reduced),'scalar_force_error',np.max(abs(printed_q-force_q)),'bulk_force_error',np.max(abs(printed_u-force_u)))
print('source_work_plus_storage',Pq+Pb+dH,'bulk_assigned_current_error',abs(J-(Pb+dH/2)))
K=1.69; m=1; k=2*np.arcsin(6/13); omega=np.sqrt(2*K*(1-np.cos(k))/m); phase=np.exp(1j*k)
f=K*(phase-1);vi=-1j*omega;vj=vi*phase
jmean=-.25*np.real((vi+vj)*np.conj(f))
print('mean_flux',jmean,'dispersion',omega,'gamma_flux',omega*(2*K*np.sin(k))/4,'group_speed',K*np.sin(k)/(m*omega))
