cnumpy.core.multiarray
_reconstruct
p0
(cnumpy
ndarray
p1
(I0
tp2
S'b'
p3
tp4
Rp5
(I1
(I51
tp6
cnumpy
dtype
p7
(S'O8'
p8
I0
I1
tp9
Rp10
(I3
S'|'
p11
NNNI-1
I-1
I63
tp12
bI00
(lp13
(lp14
I1
aS'H   + e- ->  H-  + hv'
p15
aS'[H]'
p16
aS'[e-]'
p17
aS''
p18
aS'[H-]'
p19
aS'[hv]'
p20
ag18
aS'k[0](y[8])'
p21
aS'lambda T_m: 1.4*(10**-18)*(T_m**0.928)*exp(-T_m/16200.0)'
p22
aS'Galli D & Palla 1998'
p23
aS'0-8000'
p24
aa(lp25
I2
aS'H-  + H  ->  H2  + e-'
p26
ag19
ag16
ag18
aS'[H2]'
p27
ag17
ag18
aS'k[1](y[8])'
p28
aS'lambda T_m: 4.0*(10**-9)*(T_m**-0.17) if T_m > 300 else (1.5*10**-9)'
p29
ag23
ag24
aa(lp30
I3
aS'H+  + H  ->  H2+ + hv'
p31
aS'[H+]'
p32
ag16
ag18
aS'[H2+]'
p33
ag20
ag18
aS'k[2](y[8])'
p34
aS'lambda T_m: 10**(-19.38-(1.523*log10(T_m))+(1.118*((log10(T_m))**2))-(0.1269*((log10(T_m))**3)))'
p35
aS'Glover & Abel 2008'
p36
ag24
aa(lp37
I4
aS'H2+ + H  ->  H2  + H+'
p38
ag33
ag16
ag18
ag27
ag32
ag18
aS'k[3](y[8])'
p39
aS'lambda T_m: 6.4*(10**-10)'
p40
ag23
ag24
aa(lp41
I6
aS'H2  + H+ ->  H2+ + H '
p42
ag27
ag32
ag18
ag33
ag16
ag18
aS'k[4](y[8])'
p43
aS'lambda T_m: 3.0*(10**-10)*exp(-21050.0/T_m)'
p44
ag23
ag24
aa(lp45
I7
aS'H2  + e- ->  2H  + e-'
p46
ag27
ag17
ag18
ag16
ag16
ag17
aS'k[5](y[8])'
p47
aS'lambda T_m: 1.91*(10**-9)*(T_m**0.136)*exp(-53407.1/T_m)'
p48
aS'Trevian, C.S & Tennyson J.2002'
p49
ag24
aa(lp50
I9
aS'H-  + e- ->  H   + 2e-'
p51
ag19
ag17
ag18
ag16
ag17
ag17
aS'k[6](y[8]*kbev)'
p52
aS'lambda T_m: exp((-18.01849334)+(2.3608522*log(T_m))-(0.28274430*((log(T_m))**2))+(1.62331664*(10**-2)*(log(T_m))**3)-(3.36501203*(10**-2)*((log(T_m))**4))+(1.17832978*(10**-2)*((log(T_m))**5))-(1.65619470*(10**-3)*((log(T_m))**6))+(1.06827520*(10**-4)*((log(T_m))**7))-(2.63128581*(10**-6)*((log(T_m))**8)))'
p53
ag36
ag24
aa(lp54
I10
aS'H-  + H  ->  2H  + e-'
p55
ag19
ag16
ag18
ag16
ag16
ag17
aS'k[7](y[8]*kbev)'
p56
aS'lambda T_m: 2.5634*(10**-9)*(T_m**1.78186) if T_m < 0.1 else exp((-20.37260896) + (1.13944933*log(T_m)) - (0.14210135*((log(T_m))**2)) + (8.4644554*(10**-3)*(log(T_m))**3) - (1.4327641*(10**-3)*((log(T_m))**4)) + (2.0122503*(10**-4)*((log(T_m))**5)) + (8.6639632*(10**-5)*((log(T_m))**6)) - (2.5850097*(10**-5)*((log(T_m))**7)) + (2.4555012*(10**-6)*((log(T_m))**8)) - (8.0683825*(10**-8)*((log(T_m))**9))) '
p57
ag36
ag24
aa(lp58
I11
aS'H-  + H+ ->  H   + H '
p59
ag19
ag32
ag18
ag16
ag16
ag18
aS'k[8](y[8])'
p60
aS'lambda T_m: 1.40*(10**-7)*((T_m/300.0)**-0.487)*exp(T_m/29300.0)'
p61
aS'Lepp,S. Stancil P.C., & Dalgarno 2002'
p62
ag24
aa(lp63
I12
aS'H-  + H+ ->  H2+ + e-'
p64
ag19
ag32
ag18
ag33
ag17
ag18
aS'k[9](y[8])'
p65
aS'lambda T_m: 6.9*(10**-9)*(T_m**-0.35) if T_m < 8000 else 9.6*(10**-7)*(T_m**-0.9)'
p66
ag36
ag24
aa(lp67
I13
aS'H2+ + e- ->  H   + H '
p68
ag33
ag17
ag18
ag16
ag16
ag18
aS'k[10](y[8])'
p69
aS'lambda T_m: 1.0*(10**-8) if T_m < 617 else 1.32*(10**-6)*(T_m**-0.76)'
p70
ag36
ag24
aa(lp71
I14
aS'H2+ + H- ->  H   + H2'
p72
ag33
ag19
ag18
ag16
ag27
ag18
aS'k[11](y[8])'
p73
aS'lambda T_m: 1.4*(10**-7)*((T_m/300.0)**-0.5)'
p74
ag36
ag24
aa(lp75
I15
aS'H2  + e- ->  H   + H-'
p76
ag27
ag17
ag18
ag16
ag19
ag18
aS'k[12](y[8])'
p77
aS'lambda T_m: 36.7*(T_m**-2.28)*exp(-47172.0/T_m)'
p78
aS'Capitelli, M., Coppola, C. M., Diomede, P., & Longo, S. 2007'
p79
ag24
aa(lp80
I16
aS'H-  + hv ->  H   + e-'
p81
ag19
ag20
ag18
ag16
ag17
ag18
aS'k[13](y[7])'
p82
aS'lambda T_r: 0.144*(T_r**2.13)*exp(-8650.0/T_r)'
p83
aS'Abell et al 1997'
p84
ag24
aa(lp85
I17
aS'H2+ + hv ->  H   + H+'
p86
ag33
ag20
ag18
ag16
ag32
ag18
aS'k[14](y[7])'
p87
aS'lambda T_r: 2.0*10*(T_r**1.59)*exp(-82000.0/T_r)'
p88
aS'Stancil et al 1994'
p89
ag24
aa(lp90
I18
aS'H2  + hv ->  H2+ + e-'
p91
ag27
ag20
ag18
ag33
ag17
ag18
aS'k[15](y[7])'
p92
aS'lambda T_r: 2.9*(10**2)*(T_r**1.56)*exp(-178500.0/T_r)'
p93
ag23
ag24
aa(lp94
I20
aS'H2  + hv ->  H   + H '
p95
ag27
ag20
ag18
ag16
ag16
ag18
aS'k[16](y[7])'
p96
aS'lambda T_r: 1.13*(10**6)*(T_r**0.369)*exp(-140000.0/T_r)'
p97
aS'Glover, S.C.0 & Jappsen 2007'
p98
ag24
aa(lp99
I21
aS'D-  + hv ->  D   + e-'
p100
aS'[D-]'
p101
ag20
ag18
aS'[D]'
p102
ag17
ag18
aS'k[17](y[7])'
p103
aS'lambda T_r: 1.1*(10**-1)*(T_r**2.13)*exp(-8823.0/T_r)'
p104
ag23
ag24
aa(lp105
I25
aS'HD  + hv ->  HD+ + e-'
p106
aS'[HD]'
p107
ag20
ag18
aS'[HD+]'
p108
ag17
ag18
aS'k[18](y[7])'
p109
aS'lambda T_r: 2.9*100*(T_r**1.56)*exp(-178500.0/T_r)'
p110
ag23
ag24
aa(lp111
I26
aS'H+  +  e- -> H  +  hv'
p112
aS'[D+]'
p113
ag17
ag18
ag102
ag20
ag18
aS'k[19](Z)'
p114
aS'lambda Z: (8.76*(10**-11))*((1+Z)**-0.58)'
p115
aS'Galli & Palla 1998'
p116
ag24
aa(lp117
I27
aS'D   + H+ ->  D+  + H'
p118
ag102
ag32
ag18
ag113
ag16
ag18
aS'k[20](y[8])'
p119
aS'lambda T_m: 2.0*(10**-10)*(T_m**0.402)*exp(-37.1/T_m) - 3.31*(10**-17)*(T_m**1.48)'
p120
aS'Savin 2002'
p121
ag24
aa(lp122
I28
aS'D+  + H  ->  D   + H+'
p123
ag113
ag16
ag18
ag102
ag32
ag18
aS'k[21](y[8])'
p124
aS'lambda T_m: 2.06*(10**-10)*(T_m**0.396)*exp(-33.0/T_m) + 2.03*(10**-9)*(T_m**-0.332)'
p125
ag121
ag24
aa(lp126
I29
aS'D   + H  ->  HD  + hv'
p127
ag102
ag16
ag18
ag107
ag20
ag18
aS'k[22](y[8])'
p128
aS'lambda T_m: (10**-25)*(2.80202-(6.63697*log(T_m))+(4.75619*(log(T_m)**2))-(1.39325*(log(T_m)**3))+(0.178259*(log(T_m)**4))-(0.00817097*(log(T_m)**5)))'
p129
ag36
ag24
aa(lp130
I30
aS'D   + H2 ->  HD  + H'
p131
ag102
ag27
ag18
ag107
ag16
ag18
aS'k[23](y[8])'
p132
aS'lambda T_m: 9.0*(10**-11)*exp(-3876.0/T_m) if T_m < 200 else 1.69*(10**-10)*exp((-4680/T_m) +(198800.0/(T_m**2)))'
p133
ag116
ag24
aa(lp134
I31
aS'HD+ + H  ->  HD  + H+'
p135
ag108
ag16
ag18
ag107
ag32
ag18
aS'k[24](y[8])'
p136
ag40
aS'SLD98'
p137
ag24
aa(lp138
I32
aS'D+  + H2 ->  HD  + H+'
p139
ag113
ag27
ag18
ag107
ag32
ag18
aS'k[25](y[8])'
p140
aS'lambda T_m: 1.0*(10**-9)*(0.417+0.846*log10(T_m)-0.137*(log10(T_m)**2))'
p141
aS'GP02'
p142
ag24
aa(lp143
I33
aS'HD  + H  ->  D   + H2'
p144
ag107
ag16
ag18
ag102
ag27
ag18
aS'k[26](y[8])'
p145
aS'lambda T_m: 3.2*(10**-11)*exp(-3624.0/T_m) if T_m < 200 else 5.25*(10**-11)*exp((-4430.0/T_m)+(173900.0/(T_m**2)))'
p146
ag116
ag24
aa(lp147
I34
aS'HD  + H+ ->  D+  + H2'
p148
ag107
ag32
ag18
ag113
ag27
ag18
aS'k[27](y[8])'
p149
aS'lambda T_m: 1.1*(10**-9)*exp(-488.0/T_m)'
p150
aS'Galli & Palla 2002'
p151
ag24
aa(lp152
I35
aS'D   + H+ ->  HD+ + hv'
p153
ag102
ag32
ag18
ag108
ag20
ag18
aS'k[28](y[8])'
p154
aS'lambda T_m: 10**(-19.38-1.523*log(T_m)+1.118*(log(T_m)**2)-0.1269*(log(T_m)**3))'
p155
ag116
ag24
aa(lp156
I36
aS'D+  + H  ->  HD+ + hv'
p157
ag113
ag16
ag18
ag108
ag20
ag18
aS'k[29](y[8])'
p158
ag155
ag116
ag24
aa(lp159
I38
aS'D   + e- ->  D-  + hv'
p160
ag102
ag17
ag18
ag101
ag20
ag18
aS'k[30](y[8])'
p161
aS'lambda T_m: 3.0*(10**-16)*(T_m/300.0)**0.95*exp(-T_m/9320.0)'
p162
ag137
ag24
aa(lp163
I41
aS'H-  + D  ->  H   + D-'
p164
ag19
ag102
ag18
ag16
ag101
ag18
aS'k[31](y[8])'
p165
aS'lambda T_m: 6.4*(10**-9)*(T_m/300.0)**0.41'
p166
ag137
ag24
aa(lp167
I42
aS'D-  + H  ->  D   + H-'
p168
ag101
ag16
ag18
ag102
ag19
ag18
aS'k[32](y[8])'
p169
ag166
ag137
ag24
aa(lp170
I43
aS'D-  + H  ->  HD  + e-'
p171
ag101
ag16
ag18
ag107
ag17
ag18
aS'k[33](y[8])'
p172
aS'lambda T_m: 1.5*(10**-9)*(T_m/300.0)**-0.1'
p173
ag137
ag24
aa(lp174
I44
aS'D   + H- ->  HD  + e-'
p175
ag102
ag19
ag18
ag107
ag17
ag18
aS'k[34](y[8])'
p176
ag166
ag137
ag24
aa(lp177
I51
ag112
ag32
ag17
ag18
ag16
ag20
ag18
aS'k[35](Z)'
p178
aS'lambda Z: (8.76*(10**-11))*((1.0+Z)**-0.58)'
p179
ag116
ag24
aa(lp180
I52
ag112
ag16
ag20
ag18
ag32
ag17
ag18
aS'k[36](y[7],Z)'
p181
aS'lambda T_r,Z: ((2.41*(10**+15)*(T_r**1.5)*exp(-179472/T_r))*((8.76*(10**-11))*((1.0+Z)**-0.58)))'
p182
aS'Galli & Palla 19980-8000'
p183
aa(lp184
I53
aS'H   + e- ->  H+  + 2e-'
p185
ag16
ag17
ag18
ag32
ag17
ag17
aS'k[37](y[8])'
p186
aS'lambda T_m: 5.85*(10**-11)*(T_m**0.5)*exp(-157809.1/T_m)'
p187
aS'Janev et al 1987'
p188
ag24
aa(lp189
I94
ag112
ag102
ag20
ag18
ag113
ag17
ag18
aS'k[38](y[7],Z)'
p190
aS'lambda T_r, Z: 2.41*(10**15)*(T_r**1.5)*exp(-179472.0/T_r)*(8.76*(10**-11))*((1+Z)**-0.58)'
p191
ag116
ag24
aa(lp192
I100
aS'H   + hv -> H+   + e-'
p193
ag16
ag18
ag18
ag32
ag17
ag18
aS'flux[0]'
p194
aa(lp195
I101
aS'He+ + hv -> He++ + e-'
p196
aS'[He+]'
p197
ag18
ag18
aS'[He++]'
p198
ag17
ag18
aS'flux[1]'
p199
aa(lp200
I102
aS'He +  hv -> He+  + e-'
p201
aS'[He]'
p202
ag18
ag18
ag197
ag17
ag18
aS'flux[2]'
p203
aa(lp204
I103
aS'H- +  hv -> H    + e-'
p205
ag19
ag18
ag18
ag16
ag17
ag18
aS'flux[3]'
p206
aa(lp207
I104
aS'D +  hv -> D+  + e-'
p208
ag102
ag18
ag18
ag113
ag17
ag18
aS'flux[4]'
p209
aa(lp210
I105
aS'H2 +  hv -> H2+  + e-'
p211
ag27
ag18
ag18
ag33
ag17
ag18
aS'flux[5]'
p212
aa(lp213
I106
ag211
ag27
ag18
ag18
ag33
ag17
ag18
aS'flux[6]'
p214
aa(lp215
I107
ag211
ag27
ag18
ag18
ag33
ag17
ag18
aS'flux[7]'
p216
aa(lp217
I108
aS'H2+ + hv -> H    + H+'
p218
ag33
ag18
ag18
ag16
ag32
ag18
aS'flux[8]'
p219
aa(lp220
I109
aS'H2+ + hv -> 2H+  + e-'
p221
ag33
ag18
ag18
ag32
ag32
ag17
aS'flux[9]'
p222
aa(lp223
I110
aS'H2  + hv -> H2* -> H + H'
p224
ag27
ag18
ag18
ag16
ag16
ag18
aS'flux[10]'
p225
aa(lp226
I111
aS'H2  + hv -> H + H'
p227
ag27
ag18
ag18
ag16
ag16
ag18
aS'flux[11]'
p228
aatp229
b.