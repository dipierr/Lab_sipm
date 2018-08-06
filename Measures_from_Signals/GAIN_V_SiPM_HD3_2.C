/******************************************************************************\
 * GAIN_V_SiPM_HD3_2.cxx
 *
 * GAIN values obtained by Ana_Traces_SiPM.cxx (version of 06/08/2018)
 *
 * KEY POINTS:
 *  > Ana1(...)
 *  > dleddt = 5
 *  > NO trace smoothing
 *  > thr (parameter of Ana1(...)) set depending on the situation
 *
\******************************************************************************/




#define n_GAIN_1 6
#define n_GAIN_2 6
#define n_GAIN_3 6


#define h 600
#define w 800

void GAIN_V_SiPM_HD3_2();
int find_index(double *vect, double value);

void GAIN_V_SiPM_HD3_2(){

    // HV:
    // SiPM1:
    double HV_1[n_GAIN_1], errHV_1[n_GAIN_1];
    double GAIN_1[n_GAIN_1], errGAIN_1[n_GAIN_1];
    // SiPM2:
    double HV_2[n_GAIN_2], errHV_2[n_GAIN_2];
    double GAIN_2[n_GAIN_2], errGAIN_2[n_GAIN_2];
    // SiPM3:
    double HV_3[n_GAIN_3], errHV_3[n_GAIN_3];
    double GAIN_3[n_GAIN_3], errGAIN_3[n_GAIN_3];

    double HV = 0.;

    // index:
    int index = 0;

    // Initialization
    for(int i=0; i<n_GAIN_1; i++){
        HV_1[i] = errHV_1[i] = GAIN_1[i] = errGAIN_1[i] = 0.;
    }
    for(int i=0; i<n_GAIN_2; i++){
        HV_2[i] = errHV_2[i] = GAIN_2[i] = errGAIN_2[i] = 0.;
    }
    for(int i=0; i<n_GAIN_3; i++){
        HV_3[i] = errHV_3[i] = GAIN_3[i] = errGAIN_3[i] = 0.;
    }

    // HV
    HV_1[0]    = 32.00;
    errHV_1[0] =  0.01;
    for(int i=1; i<n_GAIN_1; i++){
        HV_1[i]    = HV_1[i-1]+1.;
        errHV_1[i] = errHV_1[0];
    }

    HV_2[0]    = 32.00;
    errHV_2[0] =  0.01;
    for(int i=1; i<n_GAIN_2; i++){
        HV_2[i]    = HV_2[i-1]+1.;
        errHV_2[i] = errHV_2[0];
    }

    HV_3[0]    = 32.00;
    errHV_3[0] =  0.01;
    for(int i=1; i<n_GAIN_3; i++){
        HV_3[i]    = HV_3[i-1]+1.;
        errHV_3[i] = errHV_3[0];
    }


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM1
    ///////////////////////////////////////////////////////////////////////////
    HV = 32;
    GAIN_1[find_index(GAIN_1, HV)] = 10;


    ///////////////////////////////////////////////////////////////////////////
    //      SiPM2
    ///////////////////////////////////////////////////////////////////////////



    ///////////////////////////////////////////////////////////////////////////
    //      SiPM3
    ///////////////////////////////////////////////////////////////////////////


    //------------------------------

    TGraphErrors *gV_GAIN_1  = new TGraphErrors(n_GAIN_1, HV_1, GAIN_1, errHV_1, errGAIN_1);
    TGraphErrors *gV_GAIN_2  = new TGraphErrors(n_GAIN_2, HV_2, GAIN_2, errHV_2, errGAIN_2);
    TGraphErrors *gV_GAIN_3  = new TGraphErrors(n_GAIN_3, HV_3, GAIN_3, errHV_3, errGAIN_3);


    //------------------------------

    gV_GAIN_1->SetMarkerStyle(20);
    gV_GAIN_1->SetMarkerColor(kOrange+2);
    gV_GAIN_1->SetTitle();
    gV_GAIN_1->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_1->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_2->SetMarkerStyle(20);
    gV_GAIN_2->SetMarkerColor(kOrange+2);
    gV_GAIN_2->SetTitle();
    gV_GAIN_2->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_2->GetYaxis()->SetTitle("GAIN (mV)");

    gV_GAIN_3->SetMarkerStyle(20);
    gV_GAIN_3->SetMarkerColor(kOrange+2);
    gV_GAIN_3->SetTitle();
    gV_GAIN_3->GetXaxis()->SetTitle("Bias Voltage (V)");
    gV_GAIN_3->GetYaxis()->SetTitle("GAIN (mV)");

    //------------------------------

    TCanvas *cV_GAIN_1 = new TCanvas("cV_GAIN_1", "cV_GAIN_1",w,h);
    cV_GAIN_1->SetGrid();
    gV_GAIN_1->Draw("AP");

    TCanvas *cV_GAIN_2 = new TCanvas("cV_GAIN_2", "cV_GAIN_2",w,h);
    cV_GAIN_2->SetGrid();
    gV_GAIN_2->Draw("AP");

    TCanvas *cV_GAIN_3 = new TCanvas("cV_GAIN_3", "cV_GAIN_3",w,h);
    cV_GAIN_3->SetGrid();
    gV_GAIN_3->Draw("AP");




}


int find_index(double *vect, double value){
    int index, n;
    n = sizeof(vect)/sizeof(vect[0]);

    double epsilon = 0.00001;

    for(int i=0; i<n; i++){
        if((vect[i]>value-epsilon) && (vect[i]<value+epsilon)){
            index = i;
            break;
        }
    }

    return index;

}
