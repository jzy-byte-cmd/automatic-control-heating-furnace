#programed by 23010341 金哲宇
import numpy as np
import matplotlib.pyplot as plt
import os,sys,time
from scipy.interpolate import make_interp_spline    #原spline函数

x,y=[],[]
x_name,y_name,x_unit,y_unit,hd="","","","",""

def read_file(file_name):   #读取 第1,2列数据,并进行数据处理
    global x_name,y_name,x,y
    tmp_1,tmp_2=[],[]
    count=0
    with open(file_name,"r") as file:   #读取
        for i in file:      
            tmp_1.append(i.strip().split(',')[0])
            tmp_2.append(i.strip().split(',')[1])
    if len(tmp_1)!=len(tmp_2):
        print("输入的数据无法一一对应.")
        os.system("pause")
        os._exit(-2)
    #print("in")
    
    for i in range(1,len(tmp_1)):       #转成小数
        tmp_1[i]=float(tmp_1[i])
    for i in range(1,len(tmp_2)):
        tmp_2[i]=float(tmp_2[i])



    x_name,y_name=tmp_1[0],tmp_2[0]
    tmp_1.remove(x_name)        #存储标签,移除标签
    tmp_2.remove(y_name)

    x=np.array(tmp_1)
    y=np.array(tmp_2)
    return 

def print_data():
    global x,y,x_name,y_name
    print("读取完成,数据:")
    print("|{:>10}|{:>10}|".format(x_name,y_name))
    for i in range(len(x)):
        print("|{:>10.6f}|{:>10.6f}|".format(x[i],y[i]))
    
    print("\n注:出来的图横纵坐标标签和csv第一行的2列一致(作为数据名)")
    return



def linar_pic():
    global x,y,x_name,y_name


    plt.rcParams['font.sans-serif'] = [u'SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    k,b=np.polyfit(x,y,1)

    y_p=0.364   #y实际值
    x_p=(y_p-b)/k     #x的预测值

    #print(k,b)
    plt.scatter(x,y,color="red",s=10,label="原始数据")
    plt.title(hd)
    plt.xlabel(x_name+"/"+x_unit)
    plt.ylabel(y_name+"/"+y_unit)
    #plt.xlim(75,40,5)     #反转
    plt.plot(x,x*k+b,label="拟合曲线\ny=%.6lfx%+.6lf"%(k,b))
    #plt.plot(np.linspace(x.min(),x_p,20),[y_p]*20,'-.',label="A=0.364",color="orange")
    #plt.plot([x_p]*20,np.linspace(y.min(),y_p,20),"-.",label="V(Fe)=%.6lf"%x_p,color="orange")
    #plt.scatter(x_p,y_p,color="green",s=15)
    plt.legend()
    plt.show()


def spline():
    plt.rcParams['font.sans-serif'] = [u'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    global x,y,x_name,y_name,x_unit,y_unit


    x_new,y_new=[],[]
    for i in range(0,len(x)-1):
        x_new.extend(np.linspace(x[i],x[i+1],10))
    y_new=make_interp_spline(x,y,3)(x_new) 

    plt.scatter(x,y,color="red",s=10,label="原始数据")
    plt.title(hd)
    plt.xlabel(x_name+"/"+x_unit)
    plt.ylabel(y_name+"/"+y_unit)
    plt.plot(x_new,y_new,label="拟合曲线\n")
    #plt.plot(np.linspace(x.min(),x_p,20),[y_p]*20,'-.',label="A=0.364",color="orange")
    #plt.plot([x_p]*20,np.linspace(y.min(),y_p,20),"-.",label="V(Fe)=%.6lf"%x_p,color="orange")
    #plt.scatter(x_p,y_p,color="green",s=15)
    plt.legend()
    plt.show()

    return

def gather_inform():
    global x_unit,y_unit,hd
    x_unit=input("横坐标x数据的单位:")
    y_unit=input("纵坐标y数据的单位:")
    hd=input("图像的标题:")
    return

#try:
print("曲线拟合软件v0.0.1 programed by 23010341金哲宇")
file_name=input("把csv文件拖到这里(然后回车):")
read_file(file_name)
gather_inform()
print_data()
tmp=int(input("1.拟合直线\n2.拟合曲线\n3.退出\n请选择(输入序号){ENTER}结束:"))

if tmp==1:  #拟合直线
    linar_pic()
elif tmp==2:    #拟合曲线
    spline()
elif tmp==3:
    print("退出中...")
    time.sleep(2)
    os._exit(0)
"""
except:
    print("Unexpected error:", sys.exc_info()[0])       #打印错误信息
    print("exiting...")
    os.system("pause")
    os._exit(-1)
    """


