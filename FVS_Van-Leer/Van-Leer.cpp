/*
基于马赫数的分裂
*/
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::Array;
using std::cout;
using std::endl;
using std::max;
//初始常量
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //test1 Sod
double u1 = -0.5, rou1 = 1, p1 = 0.4, u2 = 0.5, rou2 = 1, p2 = 0.4; //test2     两侧逃离，中间真空
//double u1 = 0, rou1 = 1, p1 = 1000, u2 = 0, rou2 = 1, p2 = 0.01; //test3  左侧高压
//double u1 = 0, rou1 = 1, p1 = 0.01, u2 = 0, rou2 = 1, p2 = 100; //test4   右侧高压
//double u1 = 19.5975, rou1 = 5.99924, p1 = 460.894, u2 = -6.19633, rou2 = 5.99242, p2 = 46.0950; //test5
//double u1 = 0.698, rou1 = 0.445, p1 = 3.528, u2 = 0, rou2 = 0.5, p2 = 0.571; //https://blog.csdn.net/lusongno1/article/details/90636623

double CFL = 0.5;
double gama = 1.4;
#define m 500
#define n 2000//无所谓
//const int m = 100, n = 100; //网格数 空间网格点数m+1  时间网格点数n+1
double Sod_length = 1, Time = 1;//指定激波管长度和演化时间(隔膜位置x=Sod_length/2)
double tau, h;//时间步长与空间步长
double Now_Time = 0;
//网格物理量
MatrixXd desnity = MatrixXd::Zero(m + 1, n + 1); //密度数值解（位置，时层）
MatrixXd pressure = MatrixXd::Zero(m + 1, n + 1); //压力数值解（位置，时层）
MatrixXd velocity = MatrixXd::Zero(m + 1, n + 1); //速度数值解（位置，时层）
// U中的三项
MatrixXd U1_rou = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U2_rou_u = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U3_e = MatrixXd::Zero(m + 1, n + 1);

//中间变量 每个点不一样
VectorXd c = VectorXd::Zero(m + 1);//c*c = gama * p /rou//c*c = gama * p /rou
VectorXd Ma = VectorXd::Zero(m + 1);
MatrixXd ffd = MatrixXd::Zero(m + 1, 3);//总通量导数

//分裂前的通量f
VectorXd Funf(int i, int k) { //计算每个点的通量 f
	VectorXd result = VectorXd::Zero(3);
	result(0) = desnity(i, k) * c(i) * Ma(i);
	result(1) = desnity(i, k) * c(i) * c(i) * (Ma(i) * Ma(i) + 1 / gama);
	result(2) = desnity(i, k) * c(i) * c(i) * c(i) * Ma(i) * (Ma(i) * Ma(i) / 2 + 1 / (gama - 1));
	return result;
}
//分裂后的通量f+  f-
VectorXd Fun_f(int i, int k, bool flag) //flag==1表示计算f+  否则计算f-
{
	VectorXd result = VectorXd::Zero(3);
	if (flag) {
		result(0) = desnity(i, k) * c(i) * (Ma(i) + 1) * (Ma(i) + 1) / 4;
		result(1) = result(0) * ((gama - 1) * velocity(i, k) + 2 * c(i)) / gama;
		result(2) = result(0) * ((gama - 1) * velocity(i, k) + 2 * c(i)) * ((gama - 1) * velocity(i, k) + 2 * c(i)) / (2 * gama * gama - 2);
	}
	else {
		result(0) = -desnity(i, k) * c(i) * (Ma(i) - 1) * (Ma(i) - 1) / 4;
		result(1) = result(0) * ((gama - 1) * velocity(i, k) - 2 * c(i)) / gama;
		result(2) = result(0) * ((gama - 1) * velocity(i, k) - 2 * c(i)) * ((gama - 1) * velocity(i, k) - 2 * c(i)) / (2 * gama * gama - 2);
	}
	return result;
}

//第k个时间层的通量导数ffd的计算
void layerk(int layerk)//k=0~n
{
	//迎风格式 计算ffd
	for (int i = 0; i <= m; ++i) {
		double temp = gama * pressure(i, layerk) / desnity(i, layerk);
		c(i) = sqrt(temp);
		Ma(i) = velocity(i, layerk) / c(i);
	}
	for (int i = 1; i < m; ++i) {
		
		//基于Ma的分裂
		if ( Ma(i) >= 1) { //三波全部是右传，用后差
			//cout << "layer" << layerk << " " << i << endl;
			ffd(i, Eigen::seq(0, 2)) = (Funf(i, layerk) - Funf(i - 1, layerk)) / (h);
		}
		else if (Ma(i) <= -1) {//三波全部左传，用前差
			ffd(i, Eigen::seq(0, 2)) = (Funf(i + 1, layerk) - Funf(i, layerk)) / (h);
		}
		else {//亚声速流，正常分裂
			ffd(i, Eigen::seq(0, 2)) = (Fun_f(i, layerk, 1) - Fun_f(i - 1, layerk, 1) + Fun_f(i + 1, layerk, 0) - Fun_f(i, layerk, 0)) / h;
		}
		////基于Ma的分裂  二阶迎风
		//if (Ma(i) >= 1 || i==1) { //三波全部是右传，用后差
		//	ffd(i, Eigen::seq(0, 2)) = (Funf(i, layerk) - Funf(i - 1, layerk)) / (h);
		//}
		//else if (Ma(i) <= -1 || i==m-1) {//三波全部左传，用前差
		//	ffd(i, Eigen::seq(0, 2)) = (Funf(i + 1, layerk) - Funf(i, layerk)) / (h);
		//}
		//else {//亚声速流，正常分裂
		//	ffd(i, Eigen::seq(0, 2)) = (3 * Fun_f(i, layerk, 1) - 4 * Fun_f(i - 1, layerk, 1) + Fun_f(i - 2, layerk, 1)
		//		+ (-1)*(3 * Fun_f(i, layerk, 0) - 4 * Fun_f(i + 1, layerk, 0) + Fun_f(i + 2, layerk, 0))) / (2 * h);
		//}
			
	}

	//重新计算时间步长
	double max_a = 0;//计算该层最大特征速度
	double max_c = 0;
	for (int i = 0; i <= m; ++i) {
		max_a = std::max(max_a, abs(velocity(i, layerk)));
		max_a = std::max(max_a, abs(velocity(i, layerk) - c(i)));
		max_a = std::max(max_a, abs(velocity(i, layerk) + c(i)));

		max_c = max(max_c, c(i));
	}
	//   max_a*tau<=h
	tau = CFL * h / max_a;
	//保证cfl数<=1
	std::cout <<"max_c"<<max_c<<" "<< tau * layerk << std::endl;
}

//时间离散采用显式Euler，推进到下一层
void next_layer_Euler(int layerk)
{
	//把U1、U2、U3推进到下一个时间层
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = U2_rou_u(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = U3_e(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 2);

	//计算k+1时间层的密度压力速度  注意边界
	desnity(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk + 1);
	for (int i = 1; i < m; ++i)
		velocity(i, layerk + 1) = U2_rou_u(i, layerk + 1) / desnity(i, layerk + 1);
	for (int i = 1; i < m; ++i) {
		pressure(i, layerk + 1) = (U3_e(i, layerk + 1) - desnity(i, layerk + 1) * velocity(i, layerk + 1) * velocity(i, layerk + 1) / 2) * (gama - 1);
	}
}

//时间离散采用三阶R-K，推进到下一层
void next_layer_R_K(int layerk)
{
	//把U1、U2、U3推进到下一个时间层
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = U2_rou_u(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = U3_e(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 2);
	//修正
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = (3.0 / 4) * U1_rou(Eigen::seq(1, m - 1), layerk) + (1.0 / 4) * U1_rou(Eigen::seq(1, m - 1), layerk + 1) - (1.0 / 4) * tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = (3.0 / 4) * U2_rou_u(Eigen::seq(1, m - 1), layerk) + (1.0 / 4) * U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) - (1.0 / 4) * tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = (3.0 / 4) * U3_e(Eigen::seq(1, m - 1), layerk) + (1.0 / 4) * U3_e(Eigen::seq(1, m - 1), layerk + 1) - (1.0 / 4) * tau * ffd(Eigen::seq(1, m - 1), 2);
	//再修正
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = (1.0 / 3) * U1_rou(Eigen::seq(1, m - 1), layerk) + (2.0 / 3) * U1_rou(Eigen::seq(1, m - 1), layerk + 1) - (2.0 / 3) * tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = (1.0 / 3) * U2_rou_u(Eigen::seq(1, m - 1), layerk) + (2.0 / 3) * U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) - (2.0 / 3) * tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = (1.0 / 3) * U3_e(Eigen::seq(1, m - 1), layerk) + (2.0 / 3) * U3_e(Eigen::seq(1, m - 1), layerk + 1) - (2.0 / 3) * tau * ffd(Eigen::seq(1, m - 1), 2);




	//计算k+1时间层的密度压力速度  注意边界
	desnity(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk + 1);
	for (int i = 1; i < m; ++i)
		velocity(i, layerk + 1) = U2_rou_u(i, layerk + 1) / desnity(i, layerk + 1);
	for (int i = 1; i < m; ++i)
		pressure(i, layerk + 1) = (U3_e(i, layerk + 1) - desnity(i, layerk + 1) * velocity(i, layerk + 1) * velocity(i, layerk + 1) / 2) * (gama - 1);
}


void solve()
{
	tau = Time / 100; h = Sod_length / m;
	//初始条件，初始化第0层
	for (int i = 0; i <= m; ++i) {
		if (i * Sod_length / m < Sod_length / 2) {
			desnity(i, 0) = rou1;
			pressure(i, 0) = p1;
			velocity(i, 0) = u1;
		}
		else {
			desnity(i, 0) = rou2;
			pressure(i, 0) = p2;
			velocity(i, 0) = u2;
		}
		U1_rou(i, 0) = desnity(i, 0);
		U2_rou_u(i, 0) = desnity(i, 0) * velocity(i, 0);
		U3_e(i, 0) = pressure(i, 0) / (gama - 1) + (desnity(i, 0) * velocity(i, 0) * velocity(i, 0) / 2);
	}
	//边界条件
	for (int k = 1; k <= n; ++k) {
		desnity(0, k) = rou1;
		pressure(0, k) = p1;
		velocity(0, k) = u1;
		desnity(m, k) = rou2;
		pressure(m, k) = p2;
		velocity(m, k) = u2;
	}

	for (int k = 0; k <= n; ++k) {
		std::cout << k << std::endl;
		layerk(k);
		if (k != n) {
			next_layer_Euler(k);
		}
	}

	////结果保存(选择最后1层)
	//Eigen::MatrixXd desnity1 = Eigen::MatrixXd::Zero(m + 1, 1); //密度数值解（位置，时层）
	//Eigen::MatrixXd pressure1 = Eigen::MatrixXd::Zero(m + 1, 1); //压力数值解（位置，时层）
	//Eigen::MatrixXd velocity1 = Eigen::MatrixXd::Zero(m + 1, 1); //速度数值解（位置，时层）
	//desnity1(seq(0, m), 0) = desnity(seq(0, m), n);
	//pressure1(seq(0, m), 0) = pressure(seq(0, m), n);
	//velocity1(seq(0, m), 0) = velocity(seq(0, m), n);
	//std::ofstream outfile;
	//outfile.open("desnity.dat");
	//outfile << desnity << std::endl;
	//outfile.close();
	//outfile.open("pressure.dat");
	//outfile << pressure << std::endl;
	//outfile.close();
	//outfile.open("velocity.dat");
	//outfile << velocity << std::endl;
	//outfile.close();

}
/*
u取负数算出来NaN 压力算出了负数
*/
