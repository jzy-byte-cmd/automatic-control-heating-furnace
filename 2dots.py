"""
2dots.py
@brief 利用两点法（多次平均）计算一阶时滞系统的K,时滞系数Td>0和惯性时间系数Ts
@author  23010341  
"""
import math
#时间关系 不考虑程序的规范性和可移植性etc.....


#已知vlote 为const=3.5v
vlote=3.5
time=[]     #时间
tem=[]  #温度



with open("./data.csv","r",encoding="utf-8") as f:
    data=f.read()
    data=data.split("\n")
    for line in data[1:]:       #跳过第1行的列名
        if line !='':
            line=line.split(",")
            #print(line)
            time.append(float(line[0]))
            tem.append(float(line[1]))

print((len(time)==len(tem)))        #检查数据是否对齐  对齐一般无遗漏

#由参考资料知 K=y(inf)/(U=3.5)
#  t=10400+时 系统稳态  取10400-end 的平均温度 作为y(inf)
y_inf_average=0

for i in range(int(10400*2),len(time)):
    y_inf_average+=tem[i]

y_inf_average/=(len(time)-10400*2)


print("y_inf_average=%.4f"%y_inf_average)
K=y_inf_average/3.5
print("K=%.4f"%K)


#由于时间跨度大,总数据量较大,因此我们计算1000对点 的两点法结果后取平均
#要确保两点法的正确性  t1和t2 间隔要尽可能大
#t1从 t=10, t2从t=6000开始        t=10000+时 系统稳态

index1,index2=20,12000
#print(time[index1],time[index2])       #验证

T_sum,tao_sum=0,0


for i in range(1000):
    i1=index1+i
    i2=index2+i
    t1,t2=time[i1],time[i2]     #t1 t2
    y1,y2=tem[i1],tem[i2]

    M1,M2=math.log(1-y1/y_inf_average),math.log(1-y2/y_inf_average)
    T=(t2-t1)/(M1-M2)
    tao=(t2*M1-t1*M2)/(M1-M2)
    if tao<=0: tao=0        #tao<0 则指定为0
    T_sum+=T
    tao_sum+=tao

T_ave,tao_ave=T_sum/1000,tao_sum/1000


print("tao_ave=%.4f"%tao_ave)
print("T_ave=%.4f"%T_ave)