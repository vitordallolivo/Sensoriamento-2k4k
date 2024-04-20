from scipy import integrate
import numpy as np
import pandas as pd
from Filtros_passa_baixa_passa_alta import *

beta = 0.5 # Ganho entre os dados
FrequenciaDados = 2 # 100 hz

invSampleFreq = 1/FrequenciaDados # 1/dt

# Ler dados do excel

df_dados=pd.read_excel('Dados.xlsx')

ax  = df_dados['ax'].to_list()
ay  = df_dados['ay'].to_list()
az  = df_dados['az'].to_list()
gx  = df_dados['gx'].to_list()
gy  = df_dados['gy'].to_list()
gz  = df_dados['gz'].to_list()
mx  = df_dados['mx'].to_list()
my  = df_dados['my'].to_list()
mz  = df_dados['mz'].to_list()
 
# mx = np.zeros(len(ax))
# my = np.zeros(len(ax))
# mz = np.zeros(len(ax))


q0 = 1
q1 = 0
q2 = 0
q3 = 0

roll = []
pitch = []
yaw = []


An =[]
Al = []
Ab = []


time =[]

flag = 0

for i in range(0,len(ax)-1):

    # print(ax[i])

    # Convertando para rad/s
    gx[i]=gx[i]*np.pi/180
    gy[i]=gy[i]*np.pi/180
    gz[i]=gz[i]*np.pi/180

    # Derivada dos quartenions
    qDot1 = 0.5 * (-q1 * gx[i] - q2 * gy[i] - q3 * gz[i])
    qDot2 = 0.5 * (q0 * gx[i] + q2 * gz[i] - q3 * gy[i])
    qDot3 = 0.5 * (q0 * gy[i] - q1 * gz[i] + q3 * gx[i])
    qDot4 = 0.5 * (q0 * gz[i] + q1 * gy[i] - q2 * gx[i])


    recipNorm= 1/(np.sqrt(ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i]))

    ax[i] = ax[i]*recipNorm
    ay[i] = ay[i]*recipNorm
    az[i] = az[i]*recipNorm

    if((mx[i] == 0) and (my[i] == 0 ) and (mz[i] ==0)):
        flag = 1
        

    if (flag == 0):

        if(not(((ax[i]==0) and (ay[i]==0) and (az[i]==0)))):

            recipNorm= 1/(np.sqrt(mx[i] * mx[i] + my[i] * my[i] + mz[i] * mz[i]))
            mx[i] = mx[i]*recipNorm
            my[i] = my[i]*recipNorm
            mz[i] = mz[i]*recipNorm


            _2q0mx = 2.0 * q0 * mx[i]
            _2q0my = 2.0 * q0 * my[i]
            _2q0mz = 2.0 * q0 * mz[i]
            _2q1mx = 2.0 * q1 * mx[i]
            _2q0 = 2.0 * q0
            _2q1 = 2.0 * q1
            _2q2 = 2.0 * q2
            _2q3 = 2.0 * q3
            _2q0q2 = 2.0 * q0 * q2
            _2q2q3 = 2.0 * q2 * q3
            q0q0 = q0 * q0
            q0q1 = q0 * q1
            q0q2 = q0 * q2
            q0q3 = q0 * q3
            q1q1 = q1 * q1
            q1q2 = q1 * q2
            q1q3 = q1 * q3
            q2q2 = q2 * q2
            q2q3 = q2 * q3
            q3q3 = q3 * q3

            #  Correção do referêncial
            hx = mx[i] * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx[i] * q1q1 + _2q1 * my[i] * q2 + _2q1 * mz[i] * q3 - mx[i] * q2q2 - mx[i] * q3q3
            hy = _2q0mx * q3 + my[i] * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my[i] * q1q1 + my[i] * q2q2 + _2q2 * mz[i] * q3 - my[i] * q3q3
            _2bx = np.sqrt(hx * hx + hy * hy)
            _2bz = -_2q0mx * q2 + _2q0my * q1 + mz[i] * q0q0 + _2q1mx * q3 - mz[i] * q1q1 + _2q2 * my[i] * q3 - mz[i] * q2q2 + mz[i] * q3q3
            _4bx = 2.0 * _2bx
            _4bz = 2.0 * _2bz

            # Gradiente de F
            s0 = -_2q2 * (2.0 * q1q3 - _2q0q2 - ax[i]) + _2q1 * (2.0 * q0q1 + _2q2q3 - ay[i]) - _2bz * q2 * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx[i]) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my[i]) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz[i])
            s1 = _2q3 * (2.0 * q1q3 - _2q0q2 - ax[i]) + _2q0 * (2.0* q0q1 + _2q2q3 - ay[i]) - 4.0 * q1 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az[i]) + _2bz * q3 * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx[i]) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my[i]) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz[i])
            s2 = -_2q0 * (2.0 * q1q3 - _2q0q2 - ax[i]) + _2q3 * (2.0 * q0q1 + _2q2q3 - ay[i]) - 4.0 * q2 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az[i]) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx[i]) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my[i]) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz[i])
            s3 = _2q1 * (2.0 * q1q3 - _2q0q2 - ax[i]) + _2q2 * (2.0 * q0q1 + _2q2q3 - ay[i]) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx[i]) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my[i]) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz[i])
            recipNorm = 1/(np.sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3)) 
            s0 = s0*recipNorm
            s1 = s1*recipNorm
            s2 = s2*recipNorm
            s3 = s3*recipNorm

            #  Feedback

            qDot1 = qDot1 - beta * s0
            qDot2 = qDot2 - beta * s1
            qDot3 = qDot3 - beta * s2
            qDot4 = qDot4 - beta * s3
        
        q0=q0 + qDot1 * invSampleFreq
        q1=q1 + qDot2 * invSampleFreq
        q2=q2 + qDot3 * invSampleFreq
        q3=q3 + qDot4 * invSampleFreq
   
        recipNorm = 1/(np.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3))
        q0=q0 * recipNorm
        q1=q1 * recipNorm
        q2=q2 * recipNorm
        q3=q3 * recipNorm

    else:
        if(not((ax[i]==0) and (ay[i]==0) and (az[i]==0))):
            
            # Variaveis Auxiliares
            
            _2q0 = 2.0 * q0
            _2q1 = 2.0 * q1
            _2q2 = 2.0 * q2
            _2q3 = 2.0 * q3
            _4q0 = 4.0 * q0
            _4q1 = 4.0 * q1
            _4q2 = 4.0 * q2
            _8q1 = 8.0 * q1
            _8q2 = 8.0 * q2
            q0q0 = q0 * q0
            q1q1 = q1 * q1
            q2q2 = q2 * q2
            q3q3 = q3 * q3

            s0 = _4q0 * q2q2 + _2q2 * ax[i] + _4q0 * q1q1 - _2q1 * ay[i]
            s1 = _4q1 * q3q3 - _2q3 * ax[i] + 4.0 * q0q0 * q1 - _2q0 * ay[i] - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az[i]
            s2 = 4.0 * q0q0 * q2 + _2q0 * ax[i] + _4q2 * q3q3 - _2q3 * ay[i] - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az[i]
            s3 = 4.0 * q1q1 * q3 - _2q1 * ax[i] + 4.0 * q2q2 * q3 - _2q2 * ay[i]
            recipNorm = 1/(np.sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3)) 
            s0 = s0*recipNorm
            s1 = s1*recipNorm
            s2 = s2*recipNorm
            s3 = s3*recipNorm

            qDot1 = qDot1 - beta * s0
            qDot2 = qDot2 - beta * s1
            qDot3 = qDot3 - beta * s2
            qDot4 = qDot4 - beta * s3
        
        q0=q0 + qDot1 * invSampleFreq
        q1=q1 + qDot2 * invSampleFreq
        q2=q2 + qDot3 * invSampleFreq
        q3=q3 + qDot4 * invSampleFreq
   
        recipNorm = 1/(np.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3))
        q0=q0 * recipNorm
        q1=q1 * recipNorm
        q2=q2 * recipNorm
        q3=q3 * recipNorm


        flag = 0
 


    An.append((ax[i]*q0*q0)+(ax[i]*q1*q1)-(ax[i]*q2*q2)+(2*az[i]*q0*q2)+(2*ay[i]*q1*q2)-(ax[i]*q3*q3)-(2*ay[i]*q0*q3)+(2*az[i]*q1*q3))
    Al.append((ay[i]*q0*q0)-(ay[i]*q1*q1)-(2*az[i]*q0*q1)+(ay[i]*q2*q2)+(2*ax[i]*q1*q2)-(ay[i]*q3*q3)+(2*ax[i]*q0*q3)+(2*az[i]*q2*q3))
    Ab.append((az[i]*q0*q0)-(az[i]*q1*q1)+(2*ay[i]*q0*q1)-(az[i]*q2*q2)-(2*ax[i]*q0*q2)+(az[i]*q3*q3)+(2*ax[i]*q1*q3)+(2*ay[i]*q2*q3))

    roll.append (np.arctan2(q0*q1 + q2*q3, 0.5 - q1*q1 - q2*q2))
    pitch.append(np.arcsin(-2.0 * (q1*q3 - q0*q2)))
    yaw.append(np.arctan2(q1*q2 + q0*q3, 0.5 - q2*q2 - q3*q3))

    time.append(invSampleFreq*i)


for i in range(0,len(An)):
    An[i] = An[i]*9.81
    Al[i] = Al[i]*9.81
    Ab[i] = Ab[i]*9.81 

    roll[i] = roll[i]*180/np.pi
    pitch[i] = pitch[i]*180/np.pi
    yaw[i] = yaw[i]*180/np.pi

Vn=[]
Vl=[]
Vb=[]



Vn = (integrate.cumtrapz(highpass(An,invSampleFreq,1))).tolist()
Vl = (integrate.cumtrapz(highpass(Al,invSampleFreq,1))).tolist()
Vb = (integrate.cumtrapz(highpass(Ab,invSampleFreq,1))).tolist()

Vn.append(0)
Vl.append(0)
Vb.append(0)

vn_=np.array(Vn)
vl_=np.array(Vl)
vb_=np.array(Vb)

V_modulo = (np.sqrt(vn_**2+vl_**2+vb_**2)).tolist()

data_frame = {
    'An[m/s²]' : An,
    'Al[m/s²]' : Al,
    'Ab[m/s²]' : Ab,
    'Roll[º]':roll,
    'Pitch[º]':pitch,
    'Yaw[º]' : yaw,
    'Vn[m/s]' : Vn,
    'Vl[m/s]' : Vl,
    'Vb[m/s]' : Vb,
    'V_mod[m/s]': V_modulo,
    'time[s]' : time
}

out_name='output_file.xlsx'

df_out = pd.DataFrame(data_frame)
df_out.to_excel(out_name,index=False)