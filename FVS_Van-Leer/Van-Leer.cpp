/*
����������ķ���
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
//��ʼ����
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //test1 Sod
double u1 = -0.5, rou1 = 1, p1 = 0.4, u2 = 0.5, rou2 = 1, p2 = 0.4; //test2     �������룬�м����
//double u1 = 0, rou1 = 1, p1 = 1000, u2 = 0, rou2 = 1, p2 = 0.01; //test3  ����ѹ
//double u1 = 0, rou1 = 1, p1 = 0.01, u2 = 0, rou2 = 1, p2 = 100; //test4   �Ҳ��ѹ
//double u1 = 19.5975, rou1 = 5.99924, p1 = 460.894, u2 = -6.19633, rou2 = 5.99242, p2 = 46.0950; //test5
//double u1 = 0.698, rou1 = 0.445, p1 = 3.528, u2 = 0, rou2 = 0.5, p2 = 0.571; //https://blog.csdn.net/lusongno1/article/details/90636623

double CFL = 0.5;
double gama = 1.4;
#define m 500
#define n 2000//����ν
//const int m = 100, n = 100; //������ �ռ��������m+1  ʱ���������n+1
double Sod_length = 1, Time = 1;//ָ�������ܳ��Ⱥ��ݻ�ʱ��(��Ĥλ��x=Sod_length/2)
double tau, h;//ʱ�䲽����ռ䲽��
double Now_Time = 0;
//����������
MatrixXd desnity = MatrixXd::Zero(m + 1, n + 1); //�ܶ���ֵ�⣨λ�ã�ʱ�㣩
MatrixXd pressure = MatrixXd::Zero(m + 1, n + 1); //ѹ����ֵ�⣨λ�ã�ʱ�㣩
MatrixXd velocity = MatrixXd::Zero(m + 1, n + 1); //�ٶ���ֵ�⣨λ�ã�ʱ�㣩
// U�е�����
MatrixXd U1_rou = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U2_rou_u = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U3_e = MatrixXd::Zero(m + 1, n + 1);

//�м���� ÿ���㲻һ��
VectorXd c = VectorXd::Zero(m + 1);//c*c = gama * p /rou//c*c = gama * p /rou
VectorXd Ma = VectorXd::Zero(m + 1);
MatrixXd ffd = MatrixXd::Zero(m + 1, 3);//��ͨ������

//����ǰ��ͨ��f
VectorXd Funf(int i, int k) { //����ÿ�����ͨ�� f
	VectorXd result = VectorXd::Zero(3);
	result(0) = desnity(i, k) * c(i) * Ma(i);
	result(1) = desnity(i, k) * c(i) * c(i) * (Ma(i) * Ma(i) + 1 / gama);
	result(2) = desnity(i, k) * c(i) * c(i) * c(i) * Ma(i) * (Ma(i) * Ma(i) / 2 + 1 / (gama - 1));
	return result;
}
//���Ѻ��ͨ��f+  f-
VectorXd Fun_f(int i, int k, bool flag) //flag==1��ʾ����f+  �������f-
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

//��k��ʱ����ͨ������ffd�ļ���
void layerk(int layerk)//k=0~n
{
	//ӭ���ʽ ����ffd
	for (int i = 0; i <= m; ++i) {
		double temp = gama * pressure(i, layerk) / desnity(i, layerk);
		c(i) = sqrt(temp);
		Ma(i) = velocity(i, layerk) / c(i);
	}
	for (int i = 1; i < m; ++i) {
		
		//����Ma�ķ���
		if ( Ma(i) >= 1) { //����ȫ�����Ҵ����ú��
			//cout << "layer" << layerk << " " << i << endl;
			ffd(i, Eigen::seq(0, 2)) = (Funf(i, layerk) - Funf(i - 1, layerk)) / (h);
		}
		else if (Ma(i) <= -1) {//����ȫ���󴫣���ǰ��
			ffd(i, Eigen::seq(0, 2)) = (Funf(i + 1, layerk) - Funf(i, layerk)) / (h);
		}
		else {//������������������
			ffd(i, Eigen::seq(0, 2)) = (Fun_f(i, layerk, 1) - Fun_f(i - 1, layerk, 1) + Fun_f(i + 1, layerk, 0) - Fun_f(i, layerk, 0)) / h;
		}
		////����Ma�ķ���  ����ӭ��
		//if (Ma(i) >= 1 || i==1) { //����ȫ�����Ҵ����ú��
		//	ffd(i, Eigen::seq(0, 2)) = (Funf(i, layerk) - Funf(i - 1, layerk)) / (h);
		//}
		//else if (Ma(i) <= -1 || i==m-1) {//����ȫ���󴫣���ǰ��
		//	ffd(i, Eigen::seq(0, 2)) = (Funf(i + 1, layerk) - Funf(i, layerk)) / (h);
		//}
		//else {//������������������
		//	ffd(i, Eigen::seq(0, 2)) = (3 * Fun_f(i, layerk, 1) - 4 * Fun_f(i - 1, layerk, 1) + Fun_f(i - 2, layerk, 1)
		//		+ (-1)*(3 * Fun_f(i, layerk, 0) - 4 * Fun_f(i + 1, layerk, 0) + Fun_f(i + 2, layerk, 0))) / (2 * h);
		//}
			
	}

	//���¼���ʱ�䲽��
	double max_a = 0;//����ò���������ٶ�
	double max_c = 0;
	for (int i = 0; i <= m; ++i) {
		max_a = std::max(max_a, abs(velocity(i, layerk)));
		max_a = std::max(max_a, abs(velocity(i, layerk) - c(i)));
		max_a = std::max(max_a, abs(velocity(i, layerk) + c(i)));

		max_c = max(max_c, c(i));
	}
	//   max_a*tau<=h
	tau = CFL * h / max_a;
	//��֤cfl��<=1
	std::cout <<"max_c"<<max_c<<" "<< tau * layerk << std::endl;
}

//ʱ����ɢ������ʽEuler���ƽ�����һ��
void next_layer_Euler(int layerk)
{
	//��U1��U2��U3�ƽ�����һ��ʱ���
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = U2_rou_u(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = U3_e(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 2);

	//����k+1ʱ�����ܶ�ѹ���ٶ�  ע��߽�
	desnity(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk + 1);
	for (int i = 1; i < m; ++i)
		velocity(i, layerk + 1) = U2_rou_u(i, layerk + 1) / desnity(i, layerk + 1);
	for (int i = 1; i < m; ++i) {
		pressure(i, layerk + 1) = (U3_e(i, layerk + 1) - desnity(i, layerk + 1) * velocity(i, layerk + 1) * velocity(i, layerk + 1) / 2) * (gama - 1);
	}
}

//ʱ����ɢ��������R-K���ƽ�����һ��
void next_layer_R_K(int layerk)
{
	//��U1��U2��U3�ƽ�����һ��ʱ���
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = U2_rou_u(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = U3_e(Eigen::seq(1, m - 1), layerk) - tau * ffd(Eigen::seq(1, m - 1), 2);
	//����
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = (3.0 / 4) * U1_rou(Eigen::seq(1, m - 1), layerk) + (1.0 / 4) * U1_rou(Eigen::seq(1, m - 1), layerk + 1) - (1.0 / 4) * tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = (3.0 / 4) * U2_rou_u(Eigen::seq(1, m - 1), layerk) + (1.0 / 4) * U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) - (1.0 / 4) * tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = (3.0 / 4) * U3_e(Eigen::seq(1, m - 1), layerk) + (1.0 / 4) * U3_e(Eigen::seq(1, m - 1), layerk + 1) - (1.0 / 4) * tau * ffd(Eigen::seq(1, m - 1), 2);
	//������
	U1_rou(Eigen::seq(1, m - 1), layerk + 1) = (1.0 / 3) * U1_rou(Eigen::seq(1, m - 1), layerk) + (2.0 / 3) * U1_rou(Eigen::seq(1, m - 1), layerk + 1) - (2.0 / 3) * tau * ffd(Eigen::seq(1, m - 1), 0);
	U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) = (1.0 / 3) * U2_rou_u(Eigen::seq(1, m - 1), layerk) + (2.0 / 3) * U2_rou_u(Eigen::seq(1, m - 1), layerk + 1) - (2.0 / 3) * tau * ffd(Eigen::seq(1, m - 1), 1);
	U3_e(Eigen::seq(1, m - 1), layerk + 1) = (1.0 / 3) * U3_e(Eigen::seq(1, m - 1), layerk) + (2.0 / 3) * U3_e(Eigen::seq(1, m - 1), layerk + 1) - (2.0 / 3) * tau * ffd(Eigen::seq(1, m - 1), 2);




	//����k+1ʱ�����ܶ�ѹ���ٶ�  ע��߽�
	desnity(Eigen::seq(1, m - 1), layerk + 1) = U1_rou(Eigen::seq(1, m - 1), layerk + 1);
	for (int i = 1; i < m; ++i)
		velocity(i, layerk + 1) = U2_rou_u(i, layerk + 1) / desnity(i, layerk + 1);
	for (int i = 1; i < m; ++i)
		pressure(i, layerk + 1) = (U3_e(i, layerk + 1) - desnity(i, layerk + 1) * velocity(i, layerk + 1) * velocity(i, layerk + 1) / 2) * (gama - 1);
}


void solve()
{
	tau = Time / 100; h = Sod_length / m;
	//��ʼ��������ʼ����0��
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
	//�߽�����
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

	////�������(ѡ�����1��)
	//Eigen::MatrixXd desnity1 = Eigen::MatrixXd::Zero(m + 1, 1); //�ܶ���ֵ�⣨λ�ã�ʱ�㣩
	//Eigen::MatrixXd pressure1 = Eigen::MatrixXd::Zero(m + 1, 1); //ѹ����ֵ�⣨λ�ã�ʱ�㣩
	//Eigen::MatrixXd velocity1 = Eigen::MatrixXd::Zero(m + 1, 1); //�ٶ���ֵ�⣨λ�ã�ʱ�㣩
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
uȡ���������NaN ѹ������˸���
*/
