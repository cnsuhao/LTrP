#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <opencv2/opencv.hpp>
#include <time.h>

#define N 16

bool maxV(float a, float b){
	return a > b;
}

bool maxAbsV(float a, float b){
	return abs(a) > abs(b);
}

void selectKMax(float m[], int maxLoc[], int choose, int n){
	float *maxValue = new float[n];
	for (int i = 0; i < n; i++){
		maxLoc[i] = 0;
		maxValue[i] = m[0];
	}

	for (int i = 0; i < 8; i++){
		int j = 0;
		for (; j < n; j++){
			if (choose == 0){
				if (!maxAbsV(m[i], maxValue[j])){
					break;
				}
			}
			else if (choose == 1){
				if (!maxV(m[i], maxValue[j])){
					break;
				}
			}

		}

		//printf("%d\n", j);

		if (j > 0){
			for (int k = j - 1; k > 0; k--){
				maxValue[k - 1] = maxValue[k];
				maxLoc[k - 1] = maxLoc[k];
			}
			maxValue[j - 1] = m[i];
			maxLoc[j - 1] = i;
		}
	}
	delete[]maxValue;
}

void computeM(uchar a, uchar b, uchar c, float &m){
	m = sqrt((b - a) * (b - a) * 1.0f + (c - a) * (c - a) * 1.0f);
}

void computeDir(uchar h1, uchar v1, uchar h2, uchar v2, uchar &d){
	char difV = v1 - v2;
	char difH = h1 - h2;
	if (difV >= 0 && difH >= 0){
		d = 0;
	}
	else if (difV >= 0 && difH < 0){
		d = 1;
	}
	else if (difV < 0 && difH < 0){
		d = 2;
	}
	else{
		d = 3;
	}
	//printf("(%d	%d	%d	%d)\n", v1, h1, v2, h2);
	//printf("(%d	%d)\n", difV, difH);
}

template <typename _Tp1> static
void LTrP_(cv::InputArray _src, uchar &code, int i, int j){
	int x[9] = { 1, 2, 4, 8, 16, 32, 64, 128, 256 };
	int dir[8][2] = { { -1, 0 }, { -1, -1 }, { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 } };
	int dir2[2][2] = { { 0, 1 }, { -1, 0 } };
	float m[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	uchar codeBit[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	code = 0;

	// get matrices  
	cv::Mat src = _src.getMat();
	//_Tp center = src.at<_Tp>(i, j);
	//cout<<"center"<<(int)center<<"  ";  
	float gc, gp;
	computeM(src.ptr<_Tp1>(i)[j], src.ptr<_Tp1>(dir2[0][0] + i)[dir2[0][1] + j], src.ptr<_Tp1>(dir2[1][0] + i)[dir2[1][1] + j], gc);

	for (int d = 0; d < 8; d++){
		computeM(src.ptr<_Tp1>(dir[d][0] + i)[dir[d][1] + j], src.ptr<_Tp1>(dir2[0][0] + i)[dir2[0][1] + j], src.ptr<_Tp1>(dir2[1][0] + i)[dir2[1][1] + j], gp);
		if (gp >= gc){
			code += x[d];
		}
	}
}

template <typename _Tp1> static
void CLTrP_(cv::InputArray _src, uchar &code, int i, int j){
	int x[5] = { 1, 4, 16, 64, 256 };
	//int dir[8][2] = { { -1, 0 }, { -1, -1 }, { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 } };
	int loc[4][2] {{-1, 1}, { 0, 1 }, { 1, 1 }, { 1, 0 }};
	int dir2[2][2] = { { 0, 1 }, { -1, 0 } };
	float m[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	uchar codeBit[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	code = 0;

	// get matrices  
	cv::Mat src = _src.getMat();

	int p1_col, p1_row, p2_col, p2_row;
	cv::Point2i p1h, p1v, p2h, p2v;

	for (int k = 0; k < 4; k++){
		uchar d = 0;
		p1_row = i + loc[k][0];
		p1_col = j + loc[k][1];
		p2_row = i - loc[k][0];
		p2_col = j - loc[k][1];

		p1h.y = p1_row + dir2[0][0];
		p1h.x = p1_col + dir2[0][1];
		p1v.y = p1_row + dir2[1][0];
		p1v.x = p1_col + dir2[1][1];

		p2h.y = p2_row + dir2[0][0];
		p2h.x = p2_col + dir2[0][1];
		p2v.y = p2_row + dir2[1][0];
		p2v.x = p2_col + dir2[1][1];

		//printf("(%d	%d	%d	%d)	(%d	%d	%d	%d)\n", p1h.y, p1h.x, p1v.y, p1v.x, p2h.y, p2h.x, p2v.y, p2v.x);

		computeDir(src.ptr<uchar>(p1h.y)[p1h.x],
			src.ptr<uchar>(p1v.y)[p1v.x],
			src.ptr<uchar>(p2h.y)[p2h.x],
			src.ptr<uchar>(p2v.y)[p2v.x],
			d);
		code += d*x[k];
	}
}

template <typename _Tp1> static
void DLTrP_(cv::InputArray _src, uchar &code, int i, int j){
	int x[9] = { 1, 2, 4, 8, 16, 32, 64, 128, 256 };
	int x2[5] = { 1, 4, 16, 64, 256 };
	int dir[8][2] = { { -1, 0 }, { -1, -1 }, { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 } };
	int dir2[2][2] = { { 0, 1 }, { -1, 0 } };
	float m[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	uchar codeBit[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	code = 0;

	// get matrices  
	cv::Mat src = _src.getMat();
	//_Tp center = src.at<_Tp>(i, j);
	//cout<<"center"<<(int)center<<"  ";  
	float gc, gp;
	computeM(src.ptr<_Tp1>(i)[j], src.ptr<_Tp1>(dir2[0][0] + i)[dir2[0][1] + j], src.ptr<_Tp1>(dir2[1][0] + i)[dir2[1][1] + j], gc);

	for (int d = 0; d < 8; d++){
		computeM(src.ptr<_Tp1>(dir[d][0] + i)[dir[d][1] + j], src.ptr<_Tp1>(dir2[0][0] + i)[dir2[0][1] + j], src.ptr<_Tp1>(dir2[1][0] + i)[dir2[1][1] + j], gp);
		/*if (gp >= gc){
			codeBit[d] = 1;
		}*/
		m[d] = gp;
	}

	for (int d = 0; d < 7; d++){
		if (m[d] >= m[ (d + 1) % 8]){
			code += x[d];
		}
	}
	//code += codeBit[7] * x[7];

}

template <typename _Tp1> static
void OLTrP_(cv::InputArray _src, cv::OutputArray _dst, int choose, int n){
	// allocate memory for result  
	cv::Mat src = _src.getMat();
	_dst.create(src.rows, src.cols, CV_8UC1);
	cv::Mat dst = _dst.getMat();
	// zero the result matrix  
	dst.setTo(0);
	uchar code = 0;
#pragma omp parallel for
	for (int i = 2; i < src.rows - 1; i++) {
		//std::cout << i << std::endl;
		for (int j = 2; j < src.cols - 1; j++) {
			if (choose == 0){
				CLTrP_<uchar>(src, code, i, j);
				//printf("code: %d\n", code);
			}
			else if (choose == 1){
				LTrP_<uchar>(src, code, i, j);
			}
			else if (choose == 2){
				DLTrP_<uchar>(src, code, i, j);
			}
			dst.at<unsigned char>(i, j) = code;
			//cout<<(int)code<<" ";  
			//cout<<(int)code<<endl;  
		}
	}
	std::cout << "rows " << src.rows << " cols " << src.cols << std::endl;
	std::cout << "channels " << src.channels();
	//getchar();
	// calculate patterns  
}

int main(){
	std::string OUTPATH = ".//result//";
	std::string imname[N] = { "im", "lbp", "im_cslbp", "im_csltp_t1",
		"im_csltp_t5", "im_csltp_t10", "im_csltp_2t1", "im_csltp_4t1", "im_csltp_6t1", "Baby1", "Baby2", "Baby3"
		, "Bowling1", "Cloth1", "Flowerpots", "Rocks1" };

	std::string impath[N][2] = { { "im2.ppm", "im6.ppm" }, { ".//lbp//lbpl.png", ".//lbp//lbpr.png" },
	{ ".//cslbp//iml_cslbp.png", ".//cslbp//imr_cslbp.png" }, { ".//csltp/iml_csltp_t1.png", ".//csltp/imr_csltp_t1.png" },
	{ ".//csltp/iml_csltp_t5.png", ".//csltp/imr_csltp_t5.png" }, { ".//csltp/iml_csltp_t10.png", ".//csltp/imr_csltp_t10.png" },
	{ ".//csltp/iml_csltp_2t1.png", ".//csltp/imr_csltp_2t1.png" }, { ".//csltp/iml_csltp_4t1.png", ".//csltp/imr_csltp_4t1.png" },
	{ ".//csltp/iml_csltp_6t5.png", ".//csltp/imr_csltp_6t5.png" },
	{ "C://Users//dream//Desktop//ALL-2views//Baby1//view1.png", "C://Users//dream//Desktop//ALL-2views//Baby1//view5.png" },
	{ "C://Users//dream//Desktop//ALL-2views//Baby2//view1.png", "C://Users//dream//Desktop//ALL-2views//Baby2//view5.png" },
	{ "C://Users//dream//Desktop//ALL-2views//Baby3//view1.png", "C://Users//dream//Desktop//ALL-2views//Baby3//view5.png" },

	{ "C://Users//dream//Desktop//ALL-2views//Bowling1//view1.png", "C://Users//dream//Desktop//ALL-2views//Bowling1//view5.png" },

	{ "C://Users//dream//Desktop//ALL-2views//Cloth1//view1.png", "C://Users//dream//Desktop//ALL-2views//Cloth1//view5.png" },

	{ "C://Users//dream//Desktop//ALL-2views//Flowerpots//view1.png", "C://Users//dream//Desktop//ALL-2views//Flowerpots//view5.png" },

	{ "C://Users//dream//Desktop//ALL-2views//Rocks1//view1.png", "C://Users//dream//Desktop//ALL-2views//Rocks1//view5.png" },
	};


	for (int i = 0; i < N; i++){
		int at = 0;
		printf("N : %d imput at:", N);
		scanf("%d", &at);
		printf("\n");
		impath[at][0];

		cv::Mat ldp_facel, ldp_facer;

		//2 dltrp, 1 ltrp, 0 cltrp
		int choose = 0;
		//int choose = 1;
		//int choose = 2;

		cv::Mat imageL = cv::imread(impath[at][0].c_str(), 0);
		cv::Mat imageR = cv::imread(impath[at][1].c_str(), 0);

		if (imageL.empty()) printf("read error\n");
		if (imageR.empty()) printf("read error\n");

		char strbuffer[256];
		std::string sl, sr;
		//printf("input step and radius\n");
		//sprintf_s(strbuffer, "_ltrp_L");
		//sprintf_s(strbuffer, "_cltrp_L");
		sprintf_s(strbuffer, "_dltrp_L");
		sl = strbuffer;

		//sprintf_s(strbuffer, "_ltrp_R");
		//sprintf_s(strbuffer, "_cltrp_R");
		sprintf_s(strbuffer, "_dltrp_R");
		sr = strbuffer;
		//0表示cslbp，1表示lbp
		
		clock_t startc, stopc;
		startc = clock();
		OLTrP_<uchar>(imageL, ldp_facel, choose, 3);
		OLTrP_<uchar>(imageR, ldp_facer, choose, 3);
		stopc = clock();
		printf("\n%fs\n", (stopc - startc)*1.0 / CLOCKS_PER_SEC);

		sl = OUTPATH + imname[at] + sl + ".png";
		sr = OUTPATH + imname[at] + sr + ".png";
		cv::Mat resultl = ldp_facel.clone();
		cv::Mat resultr = ldp_facer.clone();
		cv::imwrite(sl.c_str(), resultl);
		cv::imwrite(sr.c_str(), resultr);
		cv::imshow("result", resultl);
		cv::imshow("result", resultr);
	}
	cv::waitKey(0);
	return 0;

}
