/*
采用Steger-Warming通量分裂技术，离散Euler方程
一阶迎风格式 
五种初值计算全部正确
*/
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::Array;
using std::cout;
using std::endl;
using std::max;

//初始常量
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //test1 Sod
double u1 = -2, rou1 = 1, p1 = 0.4, u2 = 2, rou2 = 1, p2 = 0.4; //test2     两侧逃离，中间真空
//double u1 = 0, rou1 = 1, p1 = 1000, u2 = 0, rou2 = 1, p2 = 0.01; //test3  左侧高压
//double u1 = 0, rou1 = 1, p1 = 0.01, u2 = 0, rou2 = 1, p2 = 100; //test4   右侧高压
//double u1 = 19.5975, rou1 = 5.99924, p1 = 460.894, u2 = -6.19633, rou2 = 5.99242, p2 = 46.0950; //test5

double CFL = 0.5;
double gama = 1.4;
#define m 500
#define n 20
//const int m = 100, n = 100; //网格数 空间网格点数m+1  时间网格点数n+1
double Sod_length = 1, Time = 1;//指定激波管长度和演化时间(隔膜位置x=Sod_length/2)
double tau, h;//时间步长与空间步长
double epsilon = 0.0;// 特征分裂参数
//网格物理量
MatrixXd desnity = MatrixXd::Zero(m + 1, n + 1); //密度数值解（位置，时层）
MatrixXd pressure = MatrixXd::Zero(m + 1, n + 1); //压力数值解（位置，时层）
MatrixXd velocity = MatrixXd::Zero(m + 1, n + 1); //速度数值解（位置，时层）
// U中的三项
MatrixXd U1_rou = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U2_rou_u = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U3_e = MatrixXd::Zero(m + 1, n + 1);

//中间变量 每层的数据
VectorXd c = VectorXd::Zero(m + 1);//c*c = gama * p /rou
MatrixXd ff1 = MatrixXd::Zero(m + 1, 3);//+
MatrixXd ff2 = MatrixXd::Zero(m + 1, 3);//-
MatrixXd ffd = MatrixXd::Zero(m + 1, 3);//总通量导数

MatrixXd lambda = MatrixXd::Zero(m+1,3);
MatrixXd lambda1 = MatrixXd::Zero(m+1,3); // +
MatrixXd lambda2 = MatrixXd::Zero(m+1,3); // -
VectorXd Ma = VectorXd::Zero(m + 1);

//分裂前的通量f
VectorXd Funf(int i, int k) { //计算每个点的通量 f
	VectorXd result = VectorXd::Zero(3);
	result(0) = U2_rou_u(i, k);
	result(1) = (gama - 1) * U3_e(i, k) + (3 - gama) * U2_rou_u(i, k) * U2_rou_u(i, k) / (2 * U1_rou(i, k));
	result(2) = gama * U2_rou_u(i, k) * U3_e(i, k) / U1_rou(i, k) + (gama - 1) * U2_rou_u(i, k) * U2_rou_u(i, k) * U2_rou_u(i, k) / (2 * U1_rou(i, k) * U1_rou(i, k));
	return result;
}
//分裂后
VectorXd Fun_f(int i,int k,int flag) { //计算每个点的通量 f 或者 f+ 或者 f-
	Eigen::VectorXd result = Eigen::VectorXd::Zero(3);
	Eigen::VectorXd Lambda = Eigen::VectorXd::Zero(3);
	if (flag) { //算正通量
		Lambda(seq(0, 2)) = lambda1(i, seq(0, 2));
	}
	else {
		Lambda(seq(0, 2)) = lambda2(i, seq(0, 2));
	}
	result(0) = 
		2 * (gama - 1) * Lambda(0)
		+ Lambda(1)
		+ Lambda(2);
	result(1) = 
		2 * (gama - 1) * Lambda(0) * velocity(i, k)
		+ Lambda(1) * (velocity(i, k) - c(i))
		+ Lambda(2) * (velocity(i, k) + c(i));
	result(2) = 
		(gama - 1) * Lambda(0) * velocity(i, k) * velocity(i, k)
		+ Lambda(1) * (velocity(i, k) - c(i)) * (velocity(i, k) - c(i)) / 2
		+ Lambda(2) * (velocity(i, k) + c(i)) * (velocity(i, k) + c(i)) / 2
		+ (3 - gama) * (Lambda(1) + Lambda(2)) * c(i) * c(i) / (2 * (gama - 1));
	result(0) *= (desnity(i, k) / (2 * gama));
	result(1) *= (desnity(i, k) / (2 * gama));
	result(2) *= (desnity(i, k) / (2 * gama));
	return result;
}


//第k个时间层的通量导数ffd的计算
void layer_k(int layerk)//k=0~n
{
	//先计算每层的参数
	for (int i = 0; i <= m; ++i) {
		c(i) = sqrt(gama * pressure(i, layerk) / desnity(i, layerk));
		Ma(i) = velocity(i, layerk) / c(i);

		lambda(i,0) = velocity(i, layerk);
		lambda(i,1) = velocity(i, layerk) - c(i);
		lambda(i,2) = velocity(i, layerk) + c(i);
		//特征值分裂
		for (int k = 0; k < 3; ++k) {
			lambda1(i,k) = (lambda(i,k) + sqrt(lambda(i,k) * lambda(i,k) + epsilon * epsilon)) / 2;
			lambda2(i,k) = (lambda(i,k) - sqrt(lambda(i,k) * lambda(i,k) + epsilon * epsilon)) / 2;
		}
	}
	//迎风格式 计算ffd  这个循环一定要单独拿出来，不能并入上一个循环
	for (int i = 1; i < m; ++i) {
		//对正分裂通量f+，相当于来流往右，用后差
		//对正分裂通量f-，相当于来流往左，用前差
		if (Ma(i) >= 1) {
			ffd(i, seq(0, 2)) = (Funf(i, layerk) - Funf(i - 1, layerk)) / (h);
		}
		else if (Ma(i) <= -1) {
			ffd(i, seq(0, 2)) = (Funf(i + 1, layerk) - Funf(i, layerk)) / (h);
		}
		else {
			ffd(i, seq(0, 2)) = (Fun_f(i, layerk, 1) - Fun_f(i - 1, layerk, 1)
				+ (Fun_f(i + 1, layerk, 0) - Fun_f(i, layerk, 0))) / (h);
		}
	}
}

//时间离散采用显式Euler，推进到下一层
void next_layer_Euler(int layerk)
{
	//重新计算时间步长
	double max_a = 0;//计算该层最大特征速度
	for (int i = 0; i <= m; ++i) {
		max_a = std::max(max_a, abs(velocity(i, layerk)));
		max_a = std::max(max_a, abs(velocity(i, layerk) - c(i)));
		max_a = std::max(max_a, abs(velocity(i, layerk) + c(i)));
	}
	//   max_a*tau<=h
	tau = CFL * h / max_a;
	//保证cfl数<=1
	std::cout << "ceng:" << layerk << " " << tau * layerk << std::endl;


	//把U1、U2、U3推进到下一个时间层
	U1_rou(Eigen::seq(1, m-1 ), layerk + 1) = U1_rou(Eigen::seq(1, m-1 ), layerk) - tau * ffd(Eigen::seq(1, m-1 ), 0);
	U2_rou_u(Eigen::seq(1, m-1 ), layerk + 1) = U2_rou_u(Eigen::seq(1, m-1), layerk) - tau * ffd(Eigen::seq(1, m-1), 1);
	U3_e(Eigen::seq(1, m-1 ), layerk + 1) = U3_e(Eigen::seq(1, m-1 ), layerk) - tau * ffd(Eigen::seq(1, m-1 ), 2);
	
	//计算k+1时间层的密度压力速度  注意边界
	desnity(Eigen::seq(1, m-1), layerk + 1) = U1_rou(Eigen::seq(1, m-1), layerk + 1);
	for (int i = 1; i < m; ++i) 
		velocity(i, layerk + 1) = U2_rou_u(i, layerk + 1) / desnity(i, layerk + 1);
	for (int i = 1; i < m; ++i)
		pressure(i, layerk + 1) = (U3_e(i, layerk + 1) - desnity(i, layerk + 1)*velocity(i, layerk + 1) * velocity(i, layerk + 1) / 2)* (gama - 1);
}


void main()
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
	//边界条件
	for (int k = 1; k <= n; ++k) {
		desnity(0, k) = rou1;
		pressure(0, k) = p1;
		velocity(0, k) = u1;
		U1_rou(0, k) = desnity(0, k);
		U2_rou_u(0, k) = desnity(0, k) * velocity(0, k);
		U3_e(0, k) = pressure(0, k) / (gama - 1) + (desnity(0, k) * velocity(0, k) * velocity(0, k) / 2);

		desnity(m, k) = rou2;
		pressure(m, k) = p2;
		velocity(m, k) = u2;
		U1_rou(m, k) = desnity(m, k);
		U2_rou_u(m, k) = desnity(m, k) * velocity(m, k);
		U3_e(m, k) = pressure(m, k) / (gama - 1) + (desnity(m, k) * velocity(m, k) * velocity(m, k) / 2);

	}

	for (int k = 0; k <= n; ++k) {
		layer_k(k);
		if(k!=n)
			next_layer_Euler(k);
	}
	
	////结果保存(选择最后1层)
	//Eigen::MatrixXd desnity1 = Eigen::MatrixXd::Zero(m + 1, 1); //密度数值解（位置，时层）
	//Eigen::MatrixXd pressure1 = Eigen::MatrixXd::Zero(m + 1, 1); //压力数值解（位置，时层）
	//Eigen::MatrixXd velocity1 = Eigen::MatrixXd::Zero(m + 1, 1); //速度数值解（位置，时层）
	//desnity1(seq(0, m), 0) = desnity(seq(0, m), n);
	//pressure1(seq(0, m), 0) = pressure(seq(0, m), n);
	//velocity1(seq(0, m), 0) = velocity(seq(0, m), n);
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
* nan可能原因：除零 根号负数
压力算出了负数 pressure(49,1)=-0.6
原因：产生了超声速流，怎么处理？？？ 处理了
隔膜打开。下一个时间层隔膜位置产生超声速激波，ffd需要重新计算
啊啊啊 超声速处理没错啊，为什么压力会算出负数，导致声速计算开根号负数了
bug：lambda 没有 /2
结果：现在压力正常了 但后面还是nan 出现了负数
bug2：计算f函数写错了 多了个括号
结果：修改了
应该是不稳定，修改了时间步长，算出来勉强离谱
return 不是0
会不会是数太小，发生除零错误
应该不是nan的问题，全部打印出来没看到啊

啊啊啊啊啊撒啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊 k=n下一步不能推进了啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊啊

最终错误：lambda写反了

大错：每点声速不一样c(i) ！！！！
最终bug：循环写错了，计算第i个点处的通量导数的时候用到了第i+1个点的通量，但它还没算出来
*/
