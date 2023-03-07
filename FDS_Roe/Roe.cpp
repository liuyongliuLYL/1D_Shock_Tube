/*
Roe近似
*/
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::Array;
using std::cout;
using std::endl;
using std::max;

//初始常量
double u1 = 0.75, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1;
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //test1 Sod
//double u1 = -2, rou1 = 1, p1 = 0.4, u2 = 2, rou2 = 1, p2 = 0.4; //test2     两侧逃离，中间真空
//double u1 = 0, rou1 = 1, p1 = 1000, u2 = 0, rou2 = 1, p2 = 0.01; //test3  左侧高压
//double u1 = 0, rou1 = 1, p1 = 0.01, u2 = 0, rou2 = 1, p2 = 100; //test4   右侧高压
//double u1 = 19.5975, rou1 = 5.99924, p1 = 460.894, u2 = -6.19633, rou2 = 5.99242, p2 = 46.0950; //test5
//double u1 = 0.698, rou1 = 0.445, p1 = 3.528, u2 = 0, rou2 = 0.5, p2 = 0.571; //https://blog.csdn.net/lusongno1/article/details/90636623

double CFL = 0.4;
double gama = 1.4;
#define m 100
#define n 130
double Sod_length = 1;
double tau, h;//时间步长与空间步长
double Now_Time = 0;

//网格物理量ρ p u c
MatrixXd desnity = MatrixXd::Zero(m + 1, n + 1); //密度数值解（位置，时层）
MatrixXd pressure = MatrixXd::Zero(m + 1, n + 1); //压力数值解（位置，时层）
MatrixXd velocity = MatrixXd::Zero(m + 1, n + 1); //速度数值解（位置，时层）
MatrixXd c = MatrixXd::Zero(m + 1, n + 1); //声速数值解（位置，时层） c*c = gama * p /rou
// U中的三项  焓H
MatrixXd U1_rou = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U2_rou_u = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U3_e = MatrixXd::Zero(m + 1, n + 1);
MatrixXd H = MatrixXd::Zero(m + 1, n + 1);

//每个网格中心的左值和右值
MatrixXd rou_L = MatrixXd::Zero(m, n + 1);  MatrixXd rou_R = MatrixXd::Zero(m, n + 1);
MatrixXd u_L = MatrixXd::Zero(m, n + 1);    MatrixXd u_R = MatrixXd::Zero(m, n + 1);
MatrixXd H_L = MatrixXd::Zero(m, n + 1);    MatrixXd H_R = MatrixXd::Zero(m, n + 1);
//每个网格中心的Roe平均rou u H p c   m+1个点对应m个网格  下标i表示 i~i+1 点之间的网格
MatrixXd ave_rou = MatrixXd::Zero(m, n + 1);
MatrixXd ave_u = MatrixXd::Zero(m, n + 1);
MatrixXd ave_H = MatrixXd::Zero(m, n + 1);
MatrixXd ave_p = MatrixXd::Zero(m, n + 1);
MatrixXd ave_c = MatrixXd::Zero(m, n + 1);
//每个网格中心的数值通量ff(i) 表示 f (i+1/2)
MatrixXd ff1 = MatrixXd::Zero(m, n + 1);  //向量三维
MatrixXd ff2 = MatrixXd::Zero(m, n + 1);
MatrixXd ff3 = MatrixXd::Zero(m, n + 1);
//每个网格节点处的通量导数ffd
MatrixXd ffd1 = MatrixXd::Zero(m + 1, n + 1);  //向量三维
MatrixXd ffd2 = MatrixXd::Zero(m + 1, n + 1);
MatrixXd ffd3 = MatrixXd::Zero(m + 1, n + 1);

//计算初始层的 rou p u U c
void InI()
{
	h = Sod_length / m;
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
		c(i, 0) = sqrt(gama * pressure(i, 0) / desnity(i, 0));

		U1_rou(i, 0) = desnity(i, 0);
		U2_rou_u(i, 0) = desnity(i, 0) * velocity(i, 0);
		U3_e(i, 0) = pressure(i, 0) / (gama - 1) + (desnity(i, 0) * velocity(i, 0) * velocity(i, 0) / 2);
		H(i, 0) = (U3_e(i, 0) + pressure(i, 0)) / desnity(i, 0);
	}
	////边界条件 要给 因为时间推进时边界无法推进
	//for (int k = 1; k <= n; ++k) {
	//	desnity(0, k) = rou1;
	//	pressure(0, k) = p1;
	//	velocity(0, k) = u1;

	//	desnity(m, k) = rou2;
	//	pressure(m, k) = p2;
	//	velocity(m, k) = u2;
	//	U1_rou(m, k) = desnity(m, k);
	//	U2_rou_u(m, k) = desnity(m, k) * velocity(m, k);
	//	U3_e(m, k) = pressure(m, k) / (gama - 1) + (desnity(m, k) * velocity(m, k) * velocity(m, k) / 2);
	//	if (desnity(m, k) == 0) c(m, k) = 0;
	//	else c(m, k) = sqrt(gama * pressure(m, k) / desnity(m, k));
	//}
}

//计算当前层的所有网格中心处的左值 右值
void cal_L_R(int k)
{
	//方案1：一阶迎风（直接取边上的一个节点的值）
	for (int i = 0; i < m; ++i) {
		//计算i+1/2处的左值L和右值R
		rou_L(i, k) = desnity(i, k); rou_R(i, k) = desnity(i + 1, k);
		u_L(i, k) = velocity(i, k);  u_R(i, k) = velocity(i + 1, k);
		H_L(i, k) = H(i, k);		 H_R(i, k) = H(i + 1, k);
	}
	//// 方案2：二阶迎风
	//for (int i = 0; i < m; ++i) {
	//	//计算i+1/2处的左值L和右值R
	//	if (i == 0 || i == m - 1) {
	//		rou_L(i, k) = desnity(i, k); rou_R(i, k) = desnity(i + 1, k);
	//		u_L(i, k) = velocity(i, k);  u_R(i, k) = velocity(i + 1, k);
	//		H_L(i, k) = H(i, k);		 H_R(i, k) = H(i + 1, k);
	//	}
	//	else {
	//		rou_L(i, k) = abs(3 * desnity(i, k) - desnity(i - 1, k)) / 2;		rou_R(i, k) = abs(3 * desnity(i + 1, k) - desnity(i + 2, k)) / 2;
	//		u_L(i, k) = (3 * velocity(i, k) - velocity(i - 1, k)) / 2;		u_R(i, k) = (3 * velocity(i + 1, k) - velocity(i + 2, k)) / 2;
	//		H_L(i, k) = abs(3 * H(i, k) - H(i - 1, k)) / 2;					H_R(i, k) = abs(3 * H(i + 1, k) - H(i + 2, k)) / 2;
	//	}
	//}
	////方案3：简单三点平均
	//for (int i = 0; i < m; ++i) {
	//	//计算i+1/2处的左值L和右值R
	//	if (i == 0 || i == m - 1) {
	//		rou_L(i, k) = desnity(i, k); rou_R(i, k) = desnity(i + 1, k);
	//		u_L(i, k) = velocity(i, k);  u_R(i, k) = velocity(i + 1, k);
	//		H_L(i, k) = H(i, k);		 H_R(i, k) = H(i + 1, k);
	//	}
	//	else {
	//		rou_L(i, k) = (desnity(i + 1, k) + desnity(i, k) + desnity(i - 1, k)) / 3;		rou_R(i, k) = (desnity(i, k) + desnity(i + 1, k) + desnity(i + 2, k)) / 3;
	//		u_L(i, k) = (velocity(i + 1, k) + velocity(i, k) + velocity(i - 1, k)) / 3;		u_R(i, k) = (velocity(i, k) + velocity(i + 1, k) + velocity(i + 2, k)) / 3;
	//		H_L(i, k) = (H(i + 1, k) + H(i, k) + H(i - 1, k)) / 3;					        H_R(i, k) = (H(i, k) + H(i + 1, k) + H(i + 2, k)) / 3;
	//	}
	//}
}
//计算当前层Roe平均，每个网格点中心的平均rou u H p c
void Roe_ave(int k)
{
	cal_L_R(k);
	for (int i = 0; i < m; ++i) {
		//Roe平均
		ave_rou(i, k) = ((sqrt(rou_L(i, k)) + sqrt(rou_R(i, k))) / 2) * ((sqrt(rou_L(i, k)) + sqrt(rou_R(i, k))) / 2);
		ave_u(i, k) = (sqrt(rou_L(i, k)) * u_L(i, k) + sqrt(rou_R(i, k)) * u_R(i, k)) / (sqrt(rou_L(i, k)) + sqrt(rou_R(i, k)));
		ave_H(i, k) = (sqrt(rou_L(i, k)) * H_L(i, k) + sqrt(rou_R(i, k)) * H_R(i, k)) / (sqrt(rou_L(i, k)) + sqrt(rou_R(i, k)));
		
		ave_p(i, k) = (ave_rou(i, k) * ave_H(i, k) - ave_rou(i, k) * ave_u(i, k) * ave_u(i, k) / 2) * (gama - 1) / gama;
		ave_c(i, k) = sqrt((gama - 1) * (ave_H(i, k) - ave_u(i, k) * ave_u(i, k) / 2));
	}
}

//物理通量f计算  传入每个节点坐标 前提：当前层所有物理量已经算出
VectorXd Fun_f(int i,int k) {
	VectorXd result = VectorXd::Zero(3);
	result(0) = desnity(i, k) * velocity(i, k);
	result(1) = desnity(i, k) * velocity(i, k) * velocity(i, k) + pressure(i, k);
	result(2) = velocity(i, k) * (U3_e(i, k) + pressure(i, k));
	return result;
}
//每个网格内部单独计算数值通量（Riemann解）（Euler方程已经近似为三个单波方程，算三个单波方程的Riemann问题）
void Riemann(int k)
{
	//计算f(i+1/2)
	for (int i = 0; i < m; ++i) {
		//第一步 计算网格i~i+1 处的平均矩阵A   A=S1 diag S
		MatrixXd A = MatrixXd::Zero(3, 3);
		MatrixXd S1 = MatrixXd::Zero(3, 3);//逆矩阵
		MatrixXd S = MatrixXd::Zero(3, 3);
		MatrixXd diag = MatrixXd::Zero(3, 3);//对角    
		{
			S(0, 0) = ave_rou(i, k) * ave_rou(i, k) / 2 - ave_c(i, k) * ave_c(i, k) / (gama - 1);
			S(0, 1) = -ave_u(i, k);
			S(0, 2) = 1;
			S(1, 0) = -ave_u(i, k) - (gama - 1) * ave_u(i, k) * ave_u(i, k) / (2 * ave_c(i, k));
			S(1, 1) = 1 + (gama - 1) * ave_u(i, k) / ave_c(i, k);
			S(1, 2) = (1 - gama) / ave_c(i, k);
			S(2, 0) = -ave_u(i, k) + (gama - 1) * ave_u(i, k) * ave_u(i, k) / (2 * ave_c(i, k));
			S(2, 1) = 1 - (gama - 1) * ave_u(i, k) / ave_c(i, k);
			S(2, 2) = (gama - 1) / ave_c(i, k);
			S1(0, 0) = (1 - gama) / (ave_c(i, k) * ave_c(i, k));
			S1(0, 1) = -1 / (2 * ave_c(i, k));
			S1(0, 2) = 1 / (2 * ave_c(i, k));
			S1(1, 0) = (1 - gama) * ave_u(i, k) / (ave_c(i, k) * ave_c(i, k));
			S1(1, 1) = -1 * (ave_u(i, k) - ave_c(i, k)) / (2 * ave_c(i, k));
			S1(1, 2) = (ave_u(i, k) + ave_c(i, k)) / (2 * ave_c(i, k));
			S1(2, 0) = (1 - gama) * ave_u(i, k) * ave_u(i, k) / (2 * ave_c(i, k) * ave_c(i, k));
			S1(2, 1) = (ave_u(i, k) * ave_u(i, k) / 2 + ave_c(i, k) * ave_c(i, k) / (gama - 1) - ave_u(i, k) * ave_c(i, k)) / (-2 * ave_c(i, k));
			S1(2, 2) = (ave_u(i, k) * ave_u(i, k) / 2 + ave_c(i, k) * ave_c(i, k) / (gama - 1) + ave_u(i, k) * ave_c(i, k)) / (2 * ave_c(i, k));
			diag(0, 0) = abs(ave_u(i, k));
			diag(1, 1) = abs(ave_u(i, k) - ave_c(i, k));
			diag(2, 2) = abs(ave_u(i, k) + ave_c(i, k));
		}
		A = S1 * diag * S;//平均增长率矩阵 绝对值

		VectorXd tempL(3); 
		tempL(0) = U1_rou(i, k); tempL(1) = U2_rou_u(i, k); tempL(2) = U3_e(i, k);
		VectorXd tempR(3);
		tempR(0) = U1_rou(i+1, k); tempR(1) = U2_rou_u(i+1, k); tempR(2) = U3_e(i+1, k);

		VectorXd temp1(3);
		temp1 = (Fun_f(i, k) + Fun_f(i + 1, k) - A * (tempR - tempL)) / 2;
		ff1(i, k) = temp1(0); ff2(i, k) = temp1(1); ff3(i, k) = temp1(2);
	}
}

//计算通量的空间导数
void ffd(int k)
{
	for (int i = 1; i < m; ++i) {
		ffd1(i, k) = (ff1(i, k) - ff1(i - 1, k)) / h;
		ffd2(i, k) = (ff2(i, k) - ff2(i - 1, k)) / h;
		ffd3(i, k) = (ff3(i, k) - ff3(i - 1, k)) / h;
	}
}

void next_layer_Euler(int k)
{
	//计算时间步长
	double max_a = 0;//计算该层最大特征速度
	for (int i = 0; i <= m; ++i) {
		max_a = std::max(max_a, abs(velocity(i, k)));
		max_a = std::max(max_a, abs(velocity(i, k) - c(i,k)));
		max_a = std::max(max_a, abs(velocity(i, k) + c(i,k)));
	}
	//   max_a*tau<=h
	tau = CFL * h / max_a;
	Now_Time += tau;

	//推进到下一个时层,计算下一个时层需要的所有物理量 边界怎么推进
	for (int i = 1; i < m; ++i) {
		U1_rou(i, k + 1) = U1_rou(i, k) - tau * ffd1(i, k);
		U2_rou_u(i, k + 1) = U2_rou_u(i, k) - tau * ffd2(i, k);
		U3_e(i, k + 1) = U3_e(i, k) - tau * ffd3(i, k);
		//根据U计算k+1层所有物理量
		desnity(i, k + 1) = U1_rou(i, k + 1);
		velocity(i, k + 1) = U2_rou_u(i, k + 1) / desnity(i, k + 1);
		pressure(i, k + 1) = (U3_e(i, k + 1) - desnity(i, k + 1) * velocity(i, k + 1) * velocity(i, k + 1) / 2) * (gama - 1);
		c(i, k + 1) = sqrt(gama * pressure(i, k + 1) / desnity(i, k + 1));
		H(i, k + 1) = (U3_e(i, k + 1) + pressure(i, k + 1)) / desnity(i, k + 1);
	}
	//边界推进 假设边界上物理量不改变
	{
		desnity(0, k + 1) = desnity(0, k); desnity(m, k + 1) = desnity(m, k);
		velocity(0, k + 1) = velocity(0, k); velocity(m, k + 1) = velocity(m, k);
		pressure(0, k + 1) = pressure(0, k); pressure(m, k + 1) = pressure(m, k);
		c(0, k + 1) = c(0, k); c(m, k + 1) = c(m, k);
		U1_rou(0, k + 1) = U1_rou(0, k); U1_rou(m, k + 1) = U1_rou(m, k);
		U2_rou_u(0, k + 1) = U2_rou_u(0, k); U2_rou_u(m, k + 1) = U2_rou_u(m, k);
		U3_e(0, k + 1) = U3_e(0, k); U3_e(m, k + 1) = U3_e(m, k);
		H(0, k + 1) = H(0, k); H(m, k + 1) = H(m, k);
	}

}

void main()
{
	InI();
	for (int k = 0; k <= n; ++k) {
		Roe_ave(k);
		Riemann(k);
		ffd(k);
		if(k!=n) next_layer_Euler(k);
		cout <<k<<" " << Now_Time << endl;
	}
	//结果保存
	std::ofstream outfile;
	outfile.open("desnity.dat");
	outfile << desnity << std::endl; cout << "done\n";
	outfile.close();
	outfile.open("pressure.dat");
	outfile << pressure << std::endl; cout << "done\n";
	outfile.close();
	outfile.open("velocity.dat");
	outfile << velocity << std::endl; cout << "done\n";
	outfile.close();
	outfile.open("c.dat");
	outfile << c << std::endl; cout << "done\n";
	outfile.close();

	//outfile.open("ave_rou.dat");
	//outfile << ave_rou << std::endl; cout << "done\n";
	//outfile.close();
	//outfile.open("ave_u.dat");
	//outfile << ave_u << std::endl; cout << "done\n";
	//outfile.close();
	//outfile.open("ave_p.dat");
	//outfile << ave_p << std::endl; cout << "done\n";
	//outfile.close();
	//outfile.open("ave_c.dat");
	//outfile << ave_c << std::endl; cout << "done\n";
	//outfile.close();
}
