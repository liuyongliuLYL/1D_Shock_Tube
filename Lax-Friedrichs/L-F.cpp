/*
Lax-Friedrichs��ʽ��ֱ���ƹ㵽Euler����
���Ե��ȶ����� ������С��1
������������ֱ���ƹ�� ������ô�ƹ�ģ��������غ���ʽ��
*/

#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>

//��ʼָ���ĳ���
double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //�����ܳ�ʼ״̬  ���������ٶ�=0
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 2, rou2 = 0.125, p2 = 0.1; //˫ϡ�貨
double gama = 1.4;
#define m 100
#define n 200
//const int m = 100, n = 100; //������ �ռ��������m+1  ʱ���������n+1
double Sod_length = 1, Time = 1;//ָ�������ܳ��Ⱥ��ݻ�ʱ��(��Ĥλ��x=Sod_length/2)
double tau, h;//ʱ�䲽����ռ䲽��
//����������
Eigen::MatrixXd desnity = Eigen::MatrixXd::Zero(m + 1, n + 1); //�ܶ���ֵ�⣨λ�ã�ʱ�㣩
Eigen::MatrixXd pressure = Eigen::MatrixXd::Zero(m + 1, n + 1); //ѹ����ֵ�⣨λ�ã�ʱ�㣩
Eigen::MatrixXd velocity = Eigen::MatrixXd::Zero(m + 1, n + 1); //�ٶ���ֵ�⣨λ�ã�ʱ�㣩
// U�е�����
Eigen::MatrixXd U1_rou = Eigen::MatrixXd::Zero(m + 1, n + 1);
Eigen::MatrixXd U2_rou_u = Eigen::MatrixXd::Zero(m + 1, n + 1);
Eigen::MatrixXd U3_e = Eigen::MatrixXd::Zero(m + 1, n + 1);

//�м���� ÿ���㲻һ��
Eigen::VectorXd ff = Eigen::VectorXd::Zero(m + 1, 3);

void solve()
{
	tau = Time / n; h = Sod_length / m;
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
	//ָ���߽�����  ���貨���ᴫ������
	for (int k = 1; k <= n; ++k) {
		desnity(0, k) = rou1;
		pressure(0, k) = p1;
		velocity(0, k) = u1;
		//�߽��ϵ�UҲҪָ��
		U1_rou(0, k) = desnity(0, k);
		U2_rou_u(0, k) = desnity(0, k) * velocity(0, k);
		U3_e(0, k) = pressure(0, k) / (gama - 1) + desnity(0, k) * velocity(0, k) * velocity(0, k) / 2;

		desnity(m, k) = rou2;
		pressure(m, k) = p2;
		velocity(m, k) = u2;
		//�߽��ϵ�UҲҪָ��
		U1_rou(m, k) = desnity(m, k);
		U2_rou_u(m, k) = desnity(m, k) * velocity(m, k);
		U3_e(m, k) = pressure(m, k) / (gama - 1) + desnity(m, k) * velocity(m, k) * velocity(m, k) / 2;
	}
	

	for (int k = 0; k < n; ++k) {
		//�����k���ͨ�� , ff(i,0~2) ��ʾ��i�����ͨ��
		for (int i = 0; i <= m; ++i) {
			ff(i, 0) = desnity(i, k) * velocity(i, k);
			ff(i, 1) = desnity(i, k) * velocity(i, k) * velocity(i, k) + pressure(i, k);
			ff(i, 2) = velocity(i, k) * (U3_e(i, k) + pressure(i, k));
		}

		//��ʽ��ǰ�ƽ���k+1��
		for (int i = 1; i < m; ++i) {
			U1_rou(i, k + 1) = (ff(i + 1, 0) - ff(i - 1, 0)) * tau / (-2 * h) + (U1_rou(i - 1, k) + U1_rou(i + 1, k)) / 2;
			U2_rou_u(i, k + 1) = (ff(i + 1, 1) - ff(i - 1, 1)) * tau / (-2 * h) + (U2_rou_u(i - 1, k) + U2_rou_u(i + 1, k)) / 2;
			U3_e(i, k + 1) = (ff(i + 1, 2) - ff(i - 1, 2)) * tau / (-2 * h) + (U3_e(i - 1, k) + U3_e(i + 1, k)) / 2;
		}
		
		//����k+1��������
		//����k+1ʱ�����ܶ�ѹ���ٶ�  ע��߽�
		desnity(Eigen::seq(1, m-1), k + 1) = U1_rou(Eigen::seq(1, m-1), k + 1);
		for (int i = 1; i < m; ++i)
			velocity(i, k + 1) = U2_rou_u(i, k + 1) / desnity(i, k + 1);
		for (int i = 1; i < m; ++i)
			pressure(i, k + 1) = (U3_e(i, k + 1) - desnity(i, k + 1) * velocity(i, k + 1) * velocity(i, k + 1) / 2) * (gama - 1);
	}




	//�������
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
