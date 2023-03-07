/*
����Steger-Warmingͨ�����Ѽ�������ɢEuler����
һ��ӭ���ʽ 
���ֳ�ֵ����ȫ����ȷ
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

//��ʼ����
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //test1 Sod
double u1 = -2, rou1 = 1, p1 = 0.4, u2 = 2, rou2 = 1, p2 = 0.4; //test2     �������룬�м����
//double u1 = 0, rou1 = 1, p1 = 1000, u2 = 0, rou2 = 1, p2 = 0.01; //test3  ����ѹ
//double u1 = 0, rou1 = 1, p1 = 0.01, u2 = 0, rou2 = 1, p2 = 100; //test4   �Ҳ��ѹ
//double u1 = 19.5975, rou1 = 5.99924, p1 = 460.894, u2 = -6.19633, rou2 = 5.99242, p2 = 46.0950; //test5

double CFL = 0.5;
double gama = 1.4;
#define m 500
#define n 20
//const int m = 100, n = 100; //������ �ռ��������m+1  ʱ���������n+1
double Sod_length = 1, Time = 1;//ָ�������ܳ��Ⱥ��ݻ�ʱ��(��Ĥλ��x=Sod_length/2)
double tau, h;//ʱ�䲽����ռ䲽��
double epsilon = 0.0;// �������Ѳ���
//����������
MatrixXd desnity = MatrixXd::Zero(m + 1, n + 1); //�ܶ���ֵ�⣨λ�ã�ʱ�㣩
MatrixXd pressure = MatrixXd::Zero(m + 1, n + 1); //ѹ����ֵ�⣨λ�ã�ʱ�㣩
MatrixXd velocity = MatrixXd::Zero(m + 1, n + 1); //�ٶ���ֵ�⣨λ�ã�ʱ�㣩
// U�е�����
MatrixXd U1_rou = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U2_rou_u = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U3_e = MatrixXd::Zero(m + 1, n + 1);

//�м���� ÿ�������
VectorXd c = VectorXd::Zero(m + 1);//c*c = gama * p /rou
MatrixXd ff1 = MatrixXd::Zero(m + 1, 3);//+
MatrixXd ff2 = MatrixXd::Zero(m + 1, 3);//-
MatrixXd ffd = MatrixXd::Zero(m + 1, 3);//��ͨ������

MatrixXd lambda = MatrixXd::Zero(m+1,3);
MatrixXd lambda1 = MatrixXd::Zero(m+1,3); // +
MatrixXd lambda2 = MatrixXd::Zero(m+1,3); // -
VectorXd Ma = VectorXd::Zero(m + 1);

//����ǰ��ͨ��f
VectorXd Funf(int i, int k) { //����ÿ�����ͨ�� f
	VectorXd result = VectorXd::Zero(3);
	result(0) = U2_rou_u(i, k);
	result(1) = (gama - 1) * U3_e(i, k) + (3 - gama) * U2_rou_u(i, k) * U2_rou_u(i, k) / (2 * U1_rou(i, k));
	result(2) = gama * U2_rou_u(i, k) * U3_e(i, k) / U1_rou(i, k) + (gama - 1) * U2_rou_u(i, k) * U2_rou_u(i, k) * U2_rou_u(i, k) / (2 * U1_rou(i, k) * U1_rou(i, k));
	return result;
}
//���Ѻ�
VectorXd Fun_f(int i,int k,int flag) { //����ÿ�����ͨ�� f ���� f+ ���� f-
	Eigen::VectorXd result = Eigen::VectorXd::Zero(3);
	Eigen::VectorXd Lambda = Eigen::VectorXd::Zero(3);
	if (flag) { //����ͨ��
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


//��k��ʱ����ͨ������ffd�ļ���
void layer_k(int layerk)//k=0~n
{
	//�ȼ���ÿ��Ĳ���
	for (int i = 0; i <= m; ++i) {
		c(i) = sqrt(gama * pressure(i, layerk) / desnity(i, layerk));
		Ma(i) = velocity(i, layerk) / c(i);

		lambda(i,0) = velocity(i, layerk);
		lambda(i,1) = velocity(i, layerk) - c(i);
		lambda(i,2) = velocity(i, layerk) + c(i);
		//����ֵ����
		for (int k = 0; k < 3; ++k) {
			lambda1(i,k) = (lambda(i,k) + sqrt(lambda(i,k) * lambda(i,k) + epsilon * epsilon)) / 2;
			lambda2(i,k) = (lambda(i,k) - sqrt(lambda(i,k) * lambda(i,k) + epsilon * epsilon)) / 2;
		}
	}
	//ӭ���ʽ ����ffd  ���ѭ��һ��Ҫ�����ó��������ܲ�����һ��ѭ��
	for (int i = 1; i < m; ++i) {
		//��������ͨ��f+���൱���������ң��ú��
		//��������ͨ��f-���൱������������ǰ��
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

//ʱ����ɢ������ʽEuler���ƽ�����һ��
void next_layer_Euler(int layerk)
{
	//���¼���ʱ�䲽��
	double max_a = 0;//����ò���������ٶ�
	for (int i = 0; i <= m; ++i) {
		max_a = std::max(max_a, abs(velocity(i, layerk)));
		max_a = std::max(max_a, abs(velocity(i, layerk) - c(i)));
		max_a = std::max(max_a, abs(velocity(i, layerk) + c(i)));
	}
	//   max_a*tau<=h
	tau = CFL * h / max_a;
	//��֤cfl��<=1
	std::cout << "ceng:" << layerk << " " << tau * layerk << std::endl;


	//��U1��U2��U3�ƽ�����һ��ʱ���
	U1_rou(Eigen::seq(1, m-1 ), layerk + 1) = U1_rou(Eigen::seq(1, m-1 ), layerk) - tau * ffd(Eigen::seq(1, m-1 ), 0);
	U2_rou_u(Eigen::seq(1, m-1 ), layerk + 1) = U2_rou_u(Eigen::seq(1, m-1), layerk) - tau * ffd(Eigen::seq(1, m-1), 1);
	U3_e(Eigen::seq(1, m-1 ), layerk + 1) = U3_e(Eigen::seq(1, m-1 ), layerk) - tau * ffd(Eigen::seq(1, m-1 ), 2);
	
	//����k+1ʱ�����ܶ�ѹ���ٶ�  ע��߽�
	desnity(Eigen::seq(1, m-1), layerk + 1) = U1_rou(Eigen::seq(1, m-1), layerk + 1);
	for (int i = 1; i < m; ++i) 
		velocity(i, layerk + 1) = U2_rou_u(i, layerk + 1) / desnity(i, layerk + 1);
	for (int i = 1; i < m; ++i)
		pressure(i, layerk + 1) = (U3_e(i, layerk + 1) - desnity(i, layerk + 1)*velocity(i, layerk + 1) * velocity(i, layerk + 1) / 2)* (gama - 1);
}


void main()
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
	//�߽�����
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
	
	////�������(ѡ�����1��)
	//Eigen::MatrixXd desnity1 = Eigen::MatrixXd::Zero(m + 1, 1); //�ܶ���ֵ�⣨λ�ã�ʱ�㣩
	//Eigen::MatrixXd pressure1 = Eigen::MatrixXd::Zero(m + 1, 1); //ѹ����ֵ�⣨λ�ã�ʱ�㣩
	//Eigen::MatrixXd velocity1 = Eigen::MatrixXd::Zero(m + 1, 1); //�ٶ���ֵ�⣨λ�ã�ʱ�㣩
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
* nan����ԭ�򣺳��� ���Ÿ���
ѹ������˸��� pressure(49,1)=-0.6
ԭ�򣺲����˳�����������ô�������� ������
��Ĥ�򿪡���һ��ʱ����Ĥλ�ò��������ټ�����ffd��Ҫ���¼���
������ �����ٴ���û����Ϊʲôѹ��������������������ټ��㿪���Ÿ�����
bug��lambda û�� /2
���������ѹ�������� �����滹��nan �����˸���
bug2������f����д���� ���˸�����
������޸���
Ӧ���ǲ��ȶ����޸���ʱ�䲽�����������ǿ����
return ����0
�᲻������̫С�������������
Ӧ�ò���nan�����⣬ȫ����ӡ����û������

���������������������������������������������������������� k=n��һ�������ƽ��˰���������������������������������

���մ���lambdaд����

���ÿ�����ٲ�һ��c(i) ��������
����bug��ѭ��д���ˣ������i���㴦��ͨ��������ʱ���õ��˵�i+1�����ͨ����������û�����
*/
