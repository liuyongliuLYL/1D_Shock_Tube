//������Lax-Wendroff ������ʽ

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

//��ʼ����
double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //test1 Sod
//double u1 = -0.5, rou1 = 1, p1 = 0.4, u2 = 0.5, rou2 = 1, p2 = 0.4; //test2     �������룬�м����
//double u1 = 0, rou1 = 1, p1 = 1000, u2 = 0, rou2 = 1, p2 = 0.01; //test3  ����ѹ
//double u1 = 0, rou1 = 1, p1 = 0.01, u2 = 0, rou2 = 1, p2 = 100; //test4   �Ҳ��ѹ
//double u1 = 19.5975, rou1 = 5.99924, p1 = 460.894, u2 = -6.19633, rou2 = 5.99242, p2 = 46.0950; //test5
//double u1 = 0.698, rou1 = 0.445, p1 = 3.528, u2 = 0, rou2 = 0.5, p2 = 0.571; //https://blog.csdn.net/lusongno1/article/details/90636623

double gama = 1.4;
#define m 400
#define n 2000
//const int m = 100, n = 100; //������ �ռ��������m+1  ʱ���������n+1
double Sod_length = 1, Time = 1;//ָ�������ܳ��Ⱥ��ݻ�ʱ��(��Ĥλ��x=Sod_length/2)
double tau, h;//ʱ�䲽����ռ䲽��
double CFL = 0.3;
double Now_Time = 0;

VectorXd c = VectorXd::Zero(m + 1);
MatrixXd U1_rou = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U2_rou_u = MatrixXd::Zero(m + 1, n + 1);
MatrixXd U3_e = MatrixXd::Zero(m + 1, n + 1);
//��ǰ���ͨ��������
MatrixXd B = MatrixXd::Zero(m + 1, 3);//3��ʾ����

MatrixXd desnity = MatrixXd::Zero(m + 1, n + 1); //�ܶ���ֵ�⣨λ�ã�ʱ�㣩
MatrixXd pressure = MatrixXd::Zero(m + 1, n + 1); //ѹ����ֵ�⣨λ�ã�ʱ�㣩
MatrixXd velocity = MatrixXd::Zero(m + 1, n + 1); //�ٶ���ֵ�⣨λ�ã�ʱ�㣩

VectorXd fun_f(VectorXd U) {
	VectorXd result = VectorXd::Zero(3);
	result(0) = U(1);
	double temp_u = U(1) / U(0);
	double temp_p = (U(2) - U(0) * temp_u * temp_u / 2) * (gama - 1);
	result(1) = U(0) * temp_u * temp_u + temp_p;
	result(2) = temp_u * (U(2) + temp_p);
	return result;
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

	//�м䲽,�±�i��ʾi+1/2  ʱ�����k+1/2   i=0~m-1
	VectorXd temp_U1_rou = VectorXd::Zero(m + 1);
	VectorXd temp_U2_rou_u = VectorXd::Zero(m + 1);
	VectorXd temp_U3_e = VectorXd::Zero(m + 1);
	//�м䲽,�±�i��ʾi-1/2  ʱ�����k+1/2   i=1~m
	VectorXd Temp_U1_rou = VectorXd::Zero(m + 1);
	VectorXd Temp_U2_rou_u = VectorXd::Zero(m + 1);
	VectorXd Temp_U3_e = VectorXd::Zero(m + 1);

	for (int k = 0; k < n; ++k) {
		//���㵱������
		for (int i = 0; i <= m; ++i) {
			double temp = gama * pressure(i, k) / desnity(i, k);
			c(i) = sqrt(temp);
		}
		//���¼���ʱ�䲽��
		double max_a = 0;//����ò���������ٶ�
		for (int i = 0; i <= m; ++i) {
			max_a = std::max(max_a, abs(velocity(i, k)));
			max_a = std::max(max_a, abs(velocity(i, k) - c(i)));
			max_a = std::max(max_a, abs(velocity(i, k) + c(i)));
		}
		//   max_a*tau<=h
		tau = CFL * h / max_a; Now_Time += tau;
		//��֤cfl��<=1
		std::cout << "max_a:" << max_a <<" "<<k<<" " << Now_Time << std::endl;


		//����ò��м䲽temp_U i+1/2
		for (int i = 0; i < m; ++i) {
			VectorXd temp1 = VectorXd::Zero(3);  VectorXd temp2 = VectorXd::Zero(3); //�ֱ��Ӧf(u i+1 n) - f(u i n)
			temp1(0) = U1_rou(i + 1, k);
			temp1(1) = U2_rou_u(i + 1, k);
			temp1(2) = U3_e(i + 1, k);
			temp2(0) = U1_rou(i, k);
			temp2(1) = U2_rou_u(i, k);
			temp2(2) = U3_e(i, k);
			VectorXd temp3 = fun_f(temp1);
			VectorXd temp4 = fun_f(temp2);
			temp_U1_rou(i) = (U1_rou(i + 1, k) + U1_rou(i, k)) / 2 - (temp3(0) - temp4(0)) * tau / (2 * h);
			temp_U2_rou_u(i) = (U2_rou_u(i + 1, k) + U2_rou_u(i, k)) / 2 - (temp3(1) - temp4(1)) * tau / (2 * h);
			temp_U3_e(i) = (U3_e(i + 1, k) + U3_e(i, k)) / 2 - (temp3(2) - temp4(2)) * tau / (2 * h);
		}
		//����ò��м䲽Temp_U i-1/2
		for (int i = 1; i <= m; ++i) {
			VectorXd temp1 = VectorXd::Zero(3);  VectorXd temp2 = VectorXd::Zero(3); //�ֱ��Ӧf(u i-1 n) - f(u i n)
			temp1(0) = U1_rou(i - 1, k);
			temp1(1) = U2_rou_u(i - 1, k);
			temp1(2) = U3_e(i - 1, k);
			temp2(0) = U1_rou(i, k);
			temp2(1) = U2_rou_u(i, k);
			temp2(2) = U3_e(i, k);
			VectorXd temp3 = fun_f(temp1);
			VectorXd temp4 = fun_f(temp2);
			Temp_U1_rou(i) = (U1_rou(i - 1, k) + U1_rou(i, k)) / 2 + (temp3(0) - temp4(0)) * tau / (2 * h);
			Temp_U2_rou_u(i) = (U2_rou_u(i - 1, k) + U2_rou_u(i, k)) / 2 + (temp3(1) - temp4(1)) * tau / (2 * h);
			Temp_U3_e(i) = (U3_e(i - 1, k) + U3_e(i, k)) / 2 + (temp3(2) - temp4(2)) * tau / (2 * h);
		}
		//����ڶ��� ������B��f(U)� B=f(temp_U) - f(Temp_U)  /h
		for (int i = 1; i < m; ++i) {
			VectorXd temp = VectorXd::Zero(3);  VectorXd Temp = VectorXd::Zero(3);  //������f�Ĳ���
			temp(0) = temp_U1_rou(i); temp(1) = temp_U2_rou_u(i); temp(2) = temp_U3_e(i);
			Temp(0) = Temp_U1_rou(i); Temp(1) = Temp_U2_rou_u(i); Temp(2) = Temp_U3_e(i);
			B(i, seq(0, 2)) = (fun_f(temp) - fun_f(Temp)) / h;
		}
		//��ǰ�ƽ�
		for (int i = 1; i < m; ++i) {
			U1_rou(i, k + 1) = U1_rou(i, k) - B(i, 0) * tau;
			U2_rou_u(i, k + 1) = U2_rou_u(i, k) - B(i, 1) * tau;
			U3_e(i, k + 1) = U3_e(i, k) - B(i, 2) * tau;
		}
		//������һ����ܶ�ѹ���ٶ�
		for (int i = 1; i < m; ++i) {
			desnity(i, k + 1) = U1_rou(i, k + 1);
			velocity(i, k + 1) = U2_rou_u(i, k + 1) / U1_rou(i, k + 1);
			pressure(i, k + 1) = (U3_e(i, k + 1) - desnity(i, k + 1) * velocity(i, k + 1) * velocity(i, k + 1) / 2) * (gama - 1);
			if (pressure(i, k + 1) < 0) pressure(i, k + 1) = 0;
		}
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
}



