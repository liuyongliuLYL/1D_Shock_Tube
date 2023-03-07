/*
Lax-Friedrichs格式可直接推广到Euler方程
线性的稳定条线 步长比小于1
啊啊啊，不能直接推广吧 不是这么推广的？不能用守恒形式？
*/

#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>

//初始指定的常量
double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //激波管初始状态  假设两边速度=0
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 2, rou2 = 0.125, p2 = 0.1; //双稀疏波
double gama = 1.4;
#define m 100
#define n 200
//const int m = 100, n = 100; //网格数 空间网格点数m+1  时间网格点数n+1
double Sod_length = 1, Time = 1;//指定激波管长度和演化时间(隔膜位置x=Sod_length/2)
double tau, h;//时间步长与空间步长
//网格物理量
Eigen::MatrixXd desnity = Eigen::MatrixXd::Zero(m + 1, n + 1); //密度数值解（位置，时层）
Eigen::MatrixXd pressure = Eigen::MatrixXd::Zero(m + 1, n + 1); //压力数值解（位置，时层）
Eigen::MatrixXd velocity = Eigen::MatrixXd::Zero(m + 1, n + 1); //速度数值解（位置，时层）
// U中的三项
Eigen::MatrixXd U1_rou = Eigen::MatrixXd::Zero(m + 1, n + 1);
Eigen::MatrixXd U2_rou_u = Eigen::MatrixXd::Zero(m + 1, n + 1);
Eigen::MatrixXd U3_e = Eigen::MatrixXd::Zero(m + 1, n + 1);

//中间变量 每个点不一样
Eigen::VectorXd ff = Eigen::VectorXd::Zero(m + 1, 3);

void solve()
{
	tau = Time / n; h = Sod_length / m;
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
	//指定边界条件  假设波不会传到那里
	for (int k = 1; k <= n; ++k) {
		desnity(0, k) = rou1;
		pressure(0, k) = p1;
		velocity(0, k) = u1;
		//边界上的U也要指定
		U1_rou(0, k) = desnity(0, k);
		U2_rou_u(0, k) = desnity(0, k) * velocity(0, k);
		U3_e(0, k) = pressure(0, k) / (gama - 1) + desnity(0, k) * velocity(0, k) * velocity(0, k) / 2;

		desnity(m, k) = rou2;
		pressure(m, k) = p2;
		velocity(m, k) = u2;
		//边界上的U也要指定
		U1_rou(m, k) = desnity(m, k);
		U2_rou_u(m, k) = desnity(m, k) * velocity(m, k);
		U3_e(m, k) = pressure(m, k) / (gama - 1) + desnity(m, k) * velocity(m, k) * velocity(m, k) / 2;
	}
	

	for (int k = 0; k < n; ++k) {
		//计算第k层的通量 , ff(i,0~2) 表示第i个点的通量
		for (int i = 0; i <= m; ++i) {
			ff(i, 0) = desnity(i, k) * velocity(i, k);
			ff(i, 1) = desnity(i, k) * velocity(i, k) * velocity(i, k) + pressure(i, k);
			ff(i, 2) = velocity(i, k) * (U3_e(i, k) + pressure(i, k));
		}

		//显式向前推进至k+1层
		for (int i = 1; i < m; ++i) {
			U1_rou(i, k + 1) = (ff(i + 1, 0) - ff(i - 1, 0)) * tau / (-2 * h) + (U1_rou(i - 1, k) + U1_rou(i + 1, k)) / 2;
			U2_rou_u(i, k + 1) = (ff(i + 1, 1) - ff(i - 1, 1)) * tau / (-2 * h) + (U2_rou_u(i - 1, k) + U2_rou_u(i + 1, k)) / 2;
			U3_e(i, k + 1) = (ff(i + 1, 2) - ff(i - 1, 2)) * tau / (-2 * h) + (U3_e(i - 1, k) + U3_e(i + 1, k)) / 2;
		}
		
		//计算k+1层物理量
		//计算k+1时间层的密度压力速度  注意边界
		desnity(Eigen::seq(1, m-1), k + 1) = U1_rou(Eigen::seq(1, m-1), k + 1);
		for (int i = 1; i < m; ++i)
			velocity(i, k + 1) = U2_rou_u(i, k + 1) / desnity(i, k + 1);
		for (int i = 1; i < m; ++i)
			pressure(i, k + 1) = (U3_e(i, k + 1) - desnity(i, k + 1) * velocity(i, k + 1) * velocity(i, k + 1) / 2) * (gama - 1);
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

	//std::cout << "done" << std::endl;

}
/*

*/
