/*
Sod�㾫ȷ��
*/
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>

//��ʼָ���ĳ���
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 0, rou2 = 0.125, p2 = 0.1; //�����ܳ�ʼ״̬  ���������ٶ�=0(��߼������ٶ����ַ���!!)
//double u1 = 0, rou1 = 1, p1 = 1, u2 = 2, rou2 = 0.125, p2 = 0.1;
double u1 = 0.698, rou1 = 0.445, p1 = 3.528, u2 = 0, rou2 = 0.5, p2 = 0.571; //https://blog.csdn.net/lusongno1/article/details/90636623

double gama = 1.4;
double Percise = 1e-5; //����Central_pressure_pstar ���ַ���������
const int m = 200, n = 600; //������ �ռ��������2*m+1  ʱ���������n+1
double Sod_length = 1, Time = 1;//ָ�������ܳ��Ⱥ��ݻ�ʱ��(��Ĥλ��x=0)
//����������
double c1, c2; //c*c = gama * p /rou
double ustar, pstar, roustarL, roustarR;//�������� �ٶ� ѹ�� ����ܶ� �Ҳ��ܶ�
double Z1, Z2, Zhead1, Ztail1, Zhead2, Ztail2;//������ϡ�貨���ٶ�
//����������
Eigen::MatrixXd desnity = Eigen::MatrixXd::Zero(2*m + 1, n + 1); //�ܶȾ�ȷ�⣨λ�ã�ʱ�㣩
Eigen::MatrixXd pressure = Eigen::MatrixXd::Zero(2*m + 1, n + 1); //ѹ����ȷ�⣨λ�ã�ʱ�㣩
Eigen::MatrixXd velocity = Eigen::MatrixXd::Zero(2*m + 1, n + 1); //�ٶȾ�ȷ�⣨λ�ã�ʱ�㣩

//����(������ϵʽ)��ϡ�貨ǰ��(���ع�ϵʽ) �ٶ�u-ѹ��p ������ϵ   
//ustar = u1 - f(pstar,p1,rou1)   ustar = u2 + f(pstar,p2,rou2)  
//������������,f1+f2 = u1-u2 ��ţ�ٵ�����pstar
double f(double pstar, int ind) {
//����pstarΪδ֪����indΪ�±�1��2
    double pi, roui, ci;
    if (ind == 1) pi = p1, roui = rou1, ci = c1;
    else pi = p2, roui = rou2, ci = c2;
	if (pstar > pi) //����
		return (pstar - pi) / (roui * ci * sqrt((gama + 1) * pstar / (2 * gama * pi) + (gama - 1) / (2 * gama)));
	else //ϡ�貨
		return 2 * ci * (pow((pstar / pi), (gama - 1) / (2 * gama)) - 1) / (gama - 1);
}
double F(double pstar) {return f(pstar, 1) + f(pstar, 2);} //��������F�����������ʺܺ�
double Central_pressure_pstar() {
/*
���ַ��󷽳�f1+f2 = u1-u2�Ľ�  ����pstar
��߽��x,�ұ߽��y
*/
    double x, y; x = std::min(p1, p2); y = std::max(p1, p2);
    double middle = (x + y) / 2;
    double percise = Percise; // ��������
    int i = 0;
    while (fabs(F(middle)) > percise) {
        F(middle)* F(x) > 0 ? x = middle : y = middle;
        middle = (x + y) / 2;
        i++;
        //printf("��������Ϊ:%d,��߽��%f,�ұ߽��%f\n", i, x, y);
    }
    //printf("���ս��Ϊ%f\n", middle);
    return middle;
}
double Central_pressure_ustar() {return (u1 + u2 + f(pstar, 2) - f(pstar, 1))/2;}

/*
//�������Ӵ��������������
double shock_wave_Z(double ui, double roui, double pi, double ci) { //��Ҫָ����ǰ״̬ ��ǰ�ٶ� �ܶ� ѹ�� ���� 
    //�������ٶ�Z
    double Ai = roui * ci * sqrt(((gama + 1) * pstar) / (2 * gama * pi) + (gama - 1) / (2 * gama));
    return ui - Ai / roui;
}
double shock_wave_roustar(double ui, double roui, double pi, double ci) {
    //�����Ĳ����ܶ�roustar L or R
    double Ai = roui * ci * sqrt(((gama + 1) * pstar) / (2 * gama * pi) + (gama - 1) / (2 * gama));
    return roui * Ai / (Ai - roui * (ui - ustar));
}
double rarefaction_wave_Zhead() {
    //ϡ�貨�Ĳ�ͷ�ٶ�   ϡ�貨�����ٴ���

}
double rarefaction_wave_Ztail() {

}
double rarefaction_wave_roustar() {

}
//�������Ӵ�����Ҳ��������
*/

//case1 ��༤��   �Ҳ༤��
void Sodcase1() {
    std::cout << "1111";
}
//case2 ���ϡ�貨 �Ҳ༤��
void Sodcase2() {
    //�������Ӵ��������������
    double cstar1 = c1 + (gama - 1) * (u1 - ustar) / 2;
    roustarL = gama * pstar / (cstar1 * cstar1);
    Zhead1 = u1 - c1; Ztail1 = ustar - cstar1;

    //�������Ӵ�����Ҳ��������
    double A2 = rou2 * c2 * sqrt(((gama + 1) * pstar) / (2 * gama * p2) + (gama - 1) / (2 * gama));
    roustarR = (rou2 * A2) / (A2 + rou2 * (u2 - ustar));
    Z2 = u2 + A2 / rou2;

    //���񻮷� ʱ���ƽ�
    int i, k;
    for (k = 0; k <= n; ++k) {//ʱ����ƽ�
        //std::cout << k << std::endl;
        double t = (double)k * Time / n; //ȷ��ʱ��
        //if (k == 50) std::cout << Zhead1 * t << " " << Ztail1 * t << " " << ustar * t <<" " << Z2 * t << std::endl;
        //-0.591608 -0.0351364 0.463726 -0.876078(Z2��������������)
        for (i = -m; i <= m; ++i) {
            double x = (double)i * Sod_length / m; //ȷ���ռ�λ��
            i += m;
            
            if (x < Zhead1 * t) {  //���ϡ�貨ǰ
                desnity(i, k) = rou1;
                pressure(i, k) = p1;
                velocity(i, k) = u1;
            }
            else if (x >= Zhead1 * t && x <= Ztail1 * t) { //���ϡ�貨����
                double c = (gama - 1) * (u1 - x / t) / (gama + 1) + 2 * c1 / (gama + 1); //������x��tλ�ô�������c
                velocity(i, k) = x / t + c;
                pressure(i, k) = p1 * pow(c / c1, (2 * gama) / (gama - 1));
                desnity(i, k) = gama * pressure(i, k) / (c * c);
            }
            else if (x > Ztail1 * t && x <= ustar * t) { //�Ӵ������� ϡ�貨����
                desnity(i, k) = roustarL;
                pressure(i, k) = pstar;
                velocity(i, k) = ustar;
            }
            else if (x > ustar * t && x <= Z2 * t) { //�Ӵ�����ұ� ��������
                desnity(i, k) = roustarR;
                pressure(i, k) = pstar;
                velocity(i, k) = ustar;
            }
            else { //������ǰ
                desnity(i, k) = rou2;
                pressure(i, k) = p2;
                velocity(i, k) = u2;
            }
            i -= m;
        }
    }
}
//case3 ��༤��   �Ҳ�ϡ�貨
void Sodcase3() {
    std::cout << "3333";
}
//case4 ���ϡ�貨 �Ҳ�ϡ�貨
void Sodcase4() {
    std::cout << "4444";
}
//case5 ���ϡ�貨 �Ҳ�ϡ�貨 �м����
void Sodcase5() {
    std::cout << "5555";
}

void Sod()
{
    //��Ҫ����ĳ���
    c1 = sqrt(gama * p1 / rou1);
    c2 = sqrt(gama * p2 / rou2);

    //�����������ѹ�����ٶ�
    pstar = Central_pressure_pstar();
    ustar = Central_pressure_ustar();

    //��5�����
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

