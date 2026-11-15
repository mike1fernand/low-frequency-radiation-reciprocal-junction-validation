from pathlib import Path
import json, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image

T=Path(__file__).resolve().parent
load=lambda p:json.loads((T/p).read_text())
chain=load('tolerance/VALIDATED_PROOF_CHAIN.json')
assert chain['status']=='PASS_RIGOROUS_INFINITE_LATTICE_CANDIDATE_RESPONSE_CERTIFICATES'
point=chain['nominal_infinite_lattice_interval'];box=chain['whole_six_motif_plus_minus5_percent_infinite_lattice_floor']
fresh=load('fresh_kernel/SUMMARY.json');spec=load('packet_spectrum/PACKET_SPECTRUM.json')
base=load('packet/CANDIDATE_PACKET_BOX80_RESULT.json')
half=load('packet/CONTROL_DT002_BOX80_RESULT.json')
large=load('packet/CONTROL_DT004_BOX96_RESULT.json')
for x in [base,half,large]:assert x['status']=='PASS'
primary=next(x for x in spec['families'] if x['label']=='surface_N160')
periodic_spectral=primary['new_ramp_full_energy_in_band_contribution']['qT']
spectral=periodic_spectral*spec['periodic_to_initial_open_total_energy_factor']
mean=lambda d:float((d['mean_low_high_total'][0]+d['mean_low_high_total'][1])/2)
qb,qh,ql=map(mean,[base,half,large])
assert base['longest']==[57.4,64.2]
assert half['selection_mode']!='coefficient-blind-primary-selector' or half['longest']==[57.4,64.2]
plot=T/'packet_validation.png'
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans'})
fig,axes=plt.subplots(1,2,figsize=(7.15,2.6),layout='constrained')
for d,label,color,style in [(base,'80 box; dt=0.04','#17628C','-'),(half,'80 box; dt=0.02','#AA6C17','--'),(large,'96 box; dt=0.04','#26765A','-.')]:
    b=np.asarray(d['branch_samples']);axes[0].plot(b[:,0],100*(b[:,1]+b[:,2])/2,style,color=color,label=label,lw=1.2)
axes[0].axhline(100*spectral,color='#6B4C86',linestyle=':',lw=1.3,label='Weighted stationary prediction')
axes[0].set(xlabel='Structure time τ',ylabel='Stored qT energy / initial energy (%)',title='Independent packet corroboration')
axes[0].legend(fontsize=6.3,loc='best')
r=base['rows'];t=np.array([x['t'] for x in r]);criteria=[('link',.0025,'Contact energy'),('J',.0005,'Current'),('shell',.002,'Kinetic shell'),('periodic_difference',.0025,'Periodic/open'),('drift',1e-5,'Energy drift')]
for key,threshold,label in criteria:axes[1].plot(t,np.maximum(1e-6,np.abs([x[key] for x in r])/threshold),label=label,lw=1.0)
axes[1].axhline(1,color='#444444',linestyle='--',lw=1)
axes[1].axvspan(57.4,64.2,color='#B6CBD7',alpha=.3)
axes[1].set(xlabel='Structure time τ',ylabel='Diagnostic / original threshold',yscale='log',ylim=(1e-4,100),title='Coefficient-blind window selection')
axes[1].legend(fontsize=6.1,ncol=2,loc='upper right')
fig.savefig(plot,dpi=220);fig.savefig(T/'packet_validation.svg');plt.close(fig)

pdfmetrics.registerFont(TTFont('DS','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DSBold','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='BodyDS',fontName='DS',fontSize=9.4,leading=13.1,spaceAfter=6,textColor=colors.HexColor('#263747')))
styles.add(ParagraphStyle(name='TitleDS',fontName='DSBold',fontSize=20,leading=25,spaceAfter=12,textColor=colors.HexColor('#153B57')))
styles.add(ParagraphStyle(name='HeadDS',fontName='DSBold',fontSize=12.5,leading=16,spaceBefore=8,spaceAfter=7,textColor=colors.HexColor('#153B57')))
styles.add(ParagraphStyle(name='SmallDS',fontName='DS',fontSize=8,leading=11,spaceAfter=7,textColor=colors.HexColor('#344758')))
story=[]
def clean(s):
    s=re.sub(r'(?<=[a-z])(?=\d)',' ',str(s))
    s=re.sub(r'(?<=\d)(?=[a-zA-Z])',' ',s)
    s=re.sub(r'(?<=[⁰¹²³⁴⁵⁶⁷⁸⁹%])(?=[A-Za-z])',' ',s)
    return s
def p(s,style='BodyDS'):story.append(Paragraph(clean(s),styles[style]))
def tab(rows,widths):
    ts=Table([[Paragraph(clean(c),styles['SmallDS']) for c in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
    ts.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E7EFF4')),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#7594AA')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F6F8FA')]),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('VALIGN',(0,0),(-1,-1),'TOP')]))
    story.append(ts);story.append(Spacer(1,9))
def footer(c,d):
    c.setStrokeColor(colors.HexColor('#CBD7E0'));c.line(47,36,565,36);c.setFont('DS',8);c.setFillColor(colors.HexColor('#546B7B'));c.drawString(47,23,'Paper One | Candidate-specific validation | 8 November 2026');c.drawRightString(565,23,str(d.page))

p('Strong validation of the 73.6% candidate','TitleDS')
p('<b>Yes: the new candidate now has rigorous response certificates and independent numerical corroboration of the same kinds used for the 12.6% benchmark.</b> The work applies to Paper One’s stated passive harmonic model. The manuscript and prior accepted results were not edited.')
tab([['Evidence','New result and scope'],['Nominal infinite-lattice certificate',"<b>73.2166% to 73.9652%</b>, displayed outward from the exact rational endpoints in the proof files."],['Continuous ±5% stiffness certificate',f"<b>qT conversion greater than 69%</b> throughout the six independent complete-motif factor intervals. The certified lower bound exceeds {100*box['lower_display']:.8f}%."],['Independent stationary calculations','Fresh surface/PV and a freshly executed tetrahedron integration support approximately73.6%. Fine-grid independent tetrahedron checks agree.'],['Independent real-time packet','Direct runs give71.6267%–71.6456%, compared with a correctly spectrum-weighted stationary prediction71.6306%. The packet is broadband, so its value differs from the carrier coefficient.']],[160,358])
p('Exact physical operating point','HeadDS')
p('The actual per-link stiffnesses are <b>g = (1,2,3,4,5,6)/7</b>, increasing along the left-incident propagation direction. All six complete eight-spring motifs remain. Frequency is omega = 6/5 in structure time; line stiffness is169/100. Lattice spacing, both masses and both bulk stiffnesses are1. Each of the48springs is positive; total physical contact stiffness is24.')
p('What “proof of73.6%” means','HeadDS')
p('The best numerical coefficient is approximately73.6%. The deterministic certificate rigorously encloses the exact model coefficient; it does <b>not</b> prove that the exact value rounds to73.6% at one decimal place. The guaranteed nominal statement is stronger than73.2%, and the guaranteed ±5%box statement is stronger than69%. These are mathematical bounds, not95%statistical confidence estimates.')
p('The old12.6%benchmark likewise combined a central numerical estimate, a robust lower-bound certificate and a separate packet calculation. Its rigorous infinite-lattice ±5%whole-motif floor was10.37%, rather than a theorem that the coefficient equals exactly12.6%. The same distinction is maintained here.')
story.append(PageBreak())

p('The mathematical proof chain','TitleDS')
p('1. Authenticate the reusable bare-bulk bounds','HeadDS')
p('The complete original R4 saved-certificate release was recovered from the authoritative R8 archive. Its unmodified verifier replayed all222saved-job states, reproduced the complete ledger exactly, and passed fresh compact physical/restart/corruption controls and the original exact-lattice bridge. The449,698integration leaves were not reintegrated. The actual G and qT-metric targets match the new candidate’s cached12×12arrays exactly.')
p('The saved global error bounds are approximately0.00390345 for eG and0.00196691 for eQ; their exact rational forms are retained in the proof chain. The new proof safely uses the larger exact rational caps <b>eG ≤1/250 and eQ ≤1/500</b>. They apply to1.69 times the Green-matrix and qT-metric differences. The old reference-contact floor is not reused as the new floor.')
p('2. Recenter the response at the actual ramp','HeadDS')
p('For the full19×19coupled system, write A = A_c + D_g ΔG P_eta, R_c = A_c⁻¹, H = P_eta R_c D_g and J = F_c R_c D_g. Subtracting the systems gives eta = (I + H ΔG)⁻¹eta_c and F_c y = F_c y_c − J ΔG eta. With verified bounds h,j,z on the corresponding norms, Z = z/(1−h eG) and d = j eG Z bound the source and output change.')
p('The true power lies between (a_low−d)²−2eQ Z² and (a_high+d)²+2eQ Z². Exact Fraction inequalities certify all matrix caps, radical enclosures and inverse-residual corrections. NumPy supplies only an inverse proposal; it is not trusted as exact. A rational LDL decomposition verifies the cached metric’s positivity. No rounded matrix square root is trusted.')
tab([['Verified conservative quantity','Bound'],['||H||₂ ≤ ||H||F','947/1000'],['||J||₂ ≤ ||J||F','3/5'],['||eta_c||','569/1000'],['Cached output amplitude a_c','4289/5000 < a_c <8579/10000'],['Exact nominal coefficient enclosure','0.7321660760687072… ≤ T_qT ≤0.7396512062343026…']],[280,238])
p('3. Certify the entire tolerance box','HeadDS')
p('For g_j = (j/7)(1+delta_j), each delta_j independently spans[−1/20,1/20]. An energy-first-order remainder bound and uniform source/output resolvent bound cover the continuous box. The64vertices bound convex auxiliary norms and an absolute linear-energy term; transmission extrema are not assumed to occur at vertices. Exact arithmetic proves the infinite-lattice floor42767052725496289/61936276900000000 >0.69. Independent source review and isolated exact replays passed.')
p('This covers six independent motif factors, with each motif’s eight links varying together. It does not cover48independently varied springs or a frequency band. Full derivations, exact endpoints and the authenticated file bindings are in tolerance/ and bound_review/.','SmallDS')
story.append(PageBreak())

p('Independent numerical corroboration','TitleDS')
p('The new proof is supported by genuinely fresh quadrature, rather than repeated solves using one cached matrix. The same actual ramp, source basis, power normalization and operating frequency were retained.')
tab([['Calculation','Incident-to-qT conversion'],['Original selected N192 / M256 kernel',f"{100*fresh['candidate']['selected_cached_qT']:.5f}%"],['Fresh N128 surface + fresh M256 PV, two-seed mean',f"{100*fresh['fresh_marching_and_PV']['new_N128_M256_seedmean_qT']:.5f}%"],['Fresh independent N48 tetrahedron integration',f"{100*fresh['fresh_independent_tetra']['qT']:.5f}%"],['Recovered independent N192 tetrahedron kernel',f"{100*fresh['independent_tetra_reused_fine_kernels']['largest_grid_qT']:.5f}%"]],[360,158])
p('Executed and preserved','HeadDS')
p('Three new marching-surface grids(N48,N80,N128); four new Sobol principal-value histograms(two independent seeds at2¹⁸ and2²⁰samples); PV orders192,256,320; one freshly executed independent N48tetrahedron integration; and a recovered independent tetrahedron ladder N96–N192 including a half-shift control. The recovered reference12.6%rows reproduced exactly before the new ramp was evaluated.')
p('All27fresh surface/PV response variants lie between73.5791% and73.6227%. The independent fine tetrahedron values lie around73.57%–73.58%. These comparisons support the approximate73.6%central estimate; they are sensitivity evidence, not the rigorous error enclosure. The exact enclosure comes from the separate certificate on page2.')
p('Separate packet-spectrum calculation','HeadDS')
p('The original120-site Gaussian line packet was retained. Its1%–99%energy interval is approximately0.932–1.473, so the correct comparison is an energy-weighted spectrum, not the transmission at omega1.2. All14distinct retained DFT frequencies were evaluated directly using five bound kernel families; both signed conjugates remain included. The numerator excludes nonpropagating bins, while the denominator retains the full initial spectral energy.')
tab([['Packet spectral method','qT / full initial spectral energy'],['Surface N128','71.63305%'],['Surface N160(primary)','71.63065%'],['Tetrahedron N96','71.61939%'],['Tetrahedron N128','71.62027%'],['Tetrahedron N160','71.62270%']],[315,203])
p('The70new-rampresponses and70reference regressions retain their internal vectors; maximum reference discrepancy is1.00×10⁻¹⁵. The retained band covers99.9683493%of initial spectral energy. The remaining passivity allowance is0.0316507percentage points. Its interval bounds omitted spectral weight only, not quadrature, finite-box, timestep or finite-time estimator errors. Initial contact energy and the small periodic/open normalization difference were retained.')
story.append(PageBreak())

p('Independent packet calculation and controls','TitleDS')
p('An autonomous open-box calculation used the actual48contact springs and initially unexcited bulk. The original56transverse box did not provide a long enough accepted interval for the stronger contact. This failed numerical control was preserved. A larger80box passed the unchanged five-predicate selector; the interval <b>57.4 ≤ tau ≤64.2</b> was frozen before subsequent controls.')
story.append(Image(str(plot),width=518,height=180));story.append(Spacer(1,7))
p('Left: stored qT energy relative to initial energy in the frozen window. The stationary overlay is added after selection. Right: the five original diagnostic thresholds; all are satisfied throughout the shaded accepted interval. These global stored-mode diagnostics are not a remote local branch flux.','SmallDS')
tab([['Direct run on the same frozen window','qT / initial energy'],['80 box; dt0.04(primary)',f'{100*qb:.8f}%'],['80 box; dt0.02',f'{100*qh:.8f}%'],['96 box; dt0.04',f'{100*ql:.8f}%'],['Weighted stationary prediction(surfaceN160)',f'{100*spectral:.8f}%']],[355,163])
p('The half-step changes the result by0.0003255%relative; enlarging the transverse box changes it by0.02647%relative. The larger-box direct value differs from the weighted stationary prediction by0.00397percentage points(0.00554%relative). These are observed comparisons, not a combined error bar. The larger-box value is slightly below the tail-only stationary interval’s lower endpoint; that interval never included these numerical sensitivities.')
p('All five reporting predicates pass in both controls. The full six-face node-energy shell fraction falls from0.18743%to0.03758%when the box grows. The current/assigned-energy residual improves from0.03062%to0.00765%when the timestep halves. Independent physical-spring one-step controls, modal-energy checks and exact checkpoint restart controls passed. A secondary cadence field omitted the final endpoint; it was corrected from saved data without changing the primary mean or selected window.')
p('Limits and deliverables','HeadDS')
p('The nominal and tolerance certificates concern the stated harmonic infinite-lattice model at one frequency. The packet evidence is separate finite-time corroboration. There is no exact one-decimal rounding proof, continuous-band guarantee,48independent-link certificate, global optimum or experimental validation. The report and core repository contain the new proofs, source-bound inputs, code and compact results. Two source parts reconstruct the unchanged complete baseline certificate for a full independent replay. No new foundation problem or physics-model change is introduced.','SmallDS')

out=T/'Paper1_736_Strong_Validation_2026-11-08.pdf'
doc=SimpleDocTemplate(str(out),pagesize=(612,792),leftMargin=47,rightMargin=47,topMargin=47,bottomMargin=47)
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(json.dumps({'report':str(out),'packet_primary':qb,'packet_halfdt':qh,'packet_largerbox':ql,'weighted_prediction':spectral},indent=2))
