/*
Sod算精确解
*/
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>

//初始指定的常量
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //激波管初始状态  假设两边速度=0(后边计算中速度区分方向!!)
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 2, rou2 = 0.125, p2 = 0.1;
double u1 = 0.698, rou1 = 0.445, p1 = 3.528, u2 = 0, rou2 = 0.5, p2 = 0.571; //https://blog.csdn.net/lusongno1/article/details/90636623

double gama = 1.4;
double Percise = 1e-5; //函数Central_pressure_pstar 二分法迭代精度
const int m = 200, n = 600; //网格数 空间网格点数2*m+1  时间网格点数n+1
double Sod_length = 1, Time = 1;//指定激波管长度和演化时间(隔膜位置x=0)
//待求物理量
double c1, c2; //c*c = gama * p /rou
double ustar, pstar, roustarL, roustarR;//中心区的 速度 压力 左侧密度 右侧密度
double Z1, Z2, Zhead1, Ztail1, Zhead2, Ztail2;//激波、稀疏波的速度
//网格物理量
Eigen::MatrixXd desnity = Eigen::MatrixXd::Zero(2*m + 1, n + 1); //密度精确解（位置，时层）
Eigen::MatrixXd pressure = Eigen::MatrixXd::Zero(2*m + 1, n + 1); //压力精确解（位置，时层）
Eigen::MatrixXd velocity = Eigen::MatrixXd::Zero(2*m + 1, n + 1); //速度精确解（位置，时层）

//激波(激波关系式)、稀疏波前后(等熵关系式) 速度u-压力p 依赖关系   
//ustar = u1 - f(pstar,p1,rou1)   ustar = u2 + f(pstar,p2,rou2)  
//联立两个方程,f1+f2 = u1-u2 再牛顿迭代出pstar
double f(double pstar, int ind) {
//输入pstar为未知数，ind为下标1或2
    double pi, roui, ci;
    if (ind == 1) pi = p1, roui = rou1, ci = c1;
    else pi = p2, roui = rou2, ci = c2;
	if (pstar > pi) //激波
		return (pstar - pi) / (roui * ci * sqrt((gama + 1) * pstar / (2 * gama * pi) + (gama - 1) / (2 * gama)));
	else //稀疏波
		return 2 * ci * (pow((pstar / pi), (gama - 1) / (2 * gama)) - 1) / (gama - 1);
}
double F(double pstar) {return f(pstar, 1) + f(pstar, 2);} //迭代函数F，单调，性质很好
double Central_pressure_pstar() {
/*
二分法求方程f1+f2 = u1-u2的解  返回pstar
左边界点x,右边界点y
*/
    double x, y; x = std::min(p1, p2); y = std::max(p1, p2);
    double middle = (x + y) / 2;
    double percise = Percise; // 迭代精度
    int i = 0;
    while (fabs(F(middle)) > percise) {
        F(middle)* F(x) > 0 ? x = middle : y = middle;
        middle = (x + y) / 2;
        i++;
        //printf("迭代次数为:%d,左边界点%f,右边界点%f\n", i, x, y);
    }
    //printf("最终结果为%f\n", middle);
    return middle;
}
double Central_pressure_ustar() {return (u1 + u2 + f(pstar, 2) - f(pstar, 1))/2;}

/*
//中心区接触间断左侧的物理量
double shock_wave_Z(double ui, double roui, double pi, double ci) { //需要指定波前状态 波前速度 密度 压力 声速 
    //激波的速度Z
    double Ai = roui * ci * sqrt(((gama + 1) * pstar) / (2 * gama * pi) + (gama - 1) / (2 * gama));
    return ui - Ai / roui;
}
double shock_wave_roustar(double ui, double roui, double pi, double ci) {
    //激波的波后密度roustar L or R
    double Ai = roui * ci * sqrt(((gama + 1) * pstar) / (2 * gama * pi) + (gama - 1) / (2 * gama));
    return roui * Ai / (Ai - roui * (ui - ustar));
}
double rarefaction_wave_Zhead() {
    //稀疏波的波头速度   稀疏波以声速传播

}
double rarefaction_wave_Ztail() {

}
double rarefaction_wave_roustar() {

}
//中心区接触间断右侧的物理量
*/

//case1 左侧激波   右侧激波
void Sodcase1() {
    std::cout << "1111";
}
//case2 左侧稀疏波 右侧激波
void Sodcase2() {
    //中心区接触间断左侧的物理量
    double cstar1 = c1 + (gama - 1) * (u1 - ustar) / 2;
    roustarL = gama * pstar / (cstar1 * cstar1);
    Zhead1 = u1 - c1; Ztail1 = ustar - cstar1;

    //中心区接触间断右侧的物理量
    double A2 = rou2 * c2 * sqrt(((gama + 1) * pstar) / (2 * gama * p2) + (gama - 1) / (2 * gama));
    roustarR = (rou2 * A2) / (A2 + rou2 * (u2 - ustar));
    Z2 = u2 + A2 / rou2;

    //网格划分 时间推进
    int i, k;
    for (k = 0; k <= n; ++k) {//时间层推进
        //std::cout << k << std::endl;
        double t = (double)k * Time / n; //确定时间
        //if (k == 50) std::cout << Zhead1 * t << " " << Ztail1 * t << " " << ustar * t <<" " << Z2 * t << std::endl;
        //-0.591608 -0.0351364 0.463726 -0.876078(Z2是正数！！！！)
        for (i = -m; i <= m; ++i) {
            double x = (double)i * Sod_length / m; //确定空间位置
            i += m;
            
            if (x < Zhead1 * t) {  //左侧稀疏波前
                desnity(i, k) = rou1;
                pressure(i, k) = p1;
                velocity(i, k) = u1;
            }
            else if (x >= Zhead1 * t && x <= Ztail1 * t) { //左侧稀疏波区域
                double c = (gama - 1) * (u1 - x / t) / (gama + 1) + 2 * c1 / (gama + 1); //中心区x，t位置处的声速c
                velocity(i, k) = x / t + c;
                pressure(i, k) = p1 * pow(c / c1, (2 * gama) / (gama - 1));
                desnity(i, k) = gama * pressure(i, k) / (c * c);
            }
            else if (x > Ztail1 * t && x <= ustar * t) { //接触间断左边 稀疏波波后
                desnity(i, k) = roustarL;
                pressure(i, k) = pstar;
                velocity(i, k) = ustar;
            }
            else if (x > ustar * t && x <= Z2 * t) { //接触间断右边 激波波后
                desnity(i, k) = roustarR;
                pressure(i, k) = pstar;
                velocity(i, k) = ustar;
            }
            else { //激波波前
                desnity(i, k) = rou2;
                pressure(i, k) = p2;
                velocity(i, k) = u2;
            }
            i -= m;
        }
    }
}
//case3 左侧激波   右侧稀疏波
void Sodcase3() {
    std::cout << "3333";
}
//case4 左侧稀疏波 右侧稀疏波
void Sodcase4() {
    std::cout << "4444";
}
//case5 左侧稀疏波 右侧稀疏波 中间真空
void Sodcase5() {
    std::cout << "5555";
}

void Sod()
{
    //需要计算的常量
    c1 = sqrt(gama * p1 / rou1);
    c2 = sqrt(gama * p2 / rou2);

    //求解中心区的压力和速度
    pstar = Central_pressure_pstar();
    ustar = Central_pressure_ustar();

    //分5种情况
    double F0 = F(0);
    double Fp1 = F(p1);
    double Fp2 = F(p2);
    double u12 = u1 - u2;
    if (p2 >= p1) {
        if (u12 > Fp2) Sodcase1();
        else if (u12 > Fp1 && u12 <= Fp2) Sodcase3();
        else if (u12 > F0 && u12 <= Fp1) Sodcase4();
        else Sodcase5();
    }
    else {
        if (u12 > Fp1) Sodcase1();
        else if (u12 > Fp2 && u12 <= Fp1) Sodcase2();
        else if (u12 > F0 && u12 <= Fp2) Sodcase4();
        else Sodcase5();
    }
    //结果保存
    std::ofstream outfile;
    outfile.open("desnity.dat");
    outfile << desnity << std::endl;
    outfile.close();
    outfile.open("pressure.dat");
    outfile << pressure << std::endl;
    outfile.close();
    outfile.open("velocity.dat");
    outfile << velocity << std::endl;
    outfile.close();
}

