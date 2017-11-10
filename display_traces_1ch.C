//void display_traces(string file_path,string file_num)

void display_traces_1ch(string file_1){
gROOT->Reset();

ifstream fileIn1 (file_1.c_str());

char dummy[10];
int n_samples1=0;
int board_id1, channel1, event_n1;
unsigned long time1;
double dato1;

string quit_tasto="q";
string prossimo="n";
string altro_eve="a";
int eve_trova;
string aaa="z"; 

int ymin = 80;
int ymax = 120;
//int ymin = 870;
//int ymax = 910;

TCanvas *c1=new TCanvas ("EVENTS"," ",0,0,1000,600);

while(aaa.compare(quit_tasto.c_str())){
	aaa="z";
	cout<<"Evento: "<<endl;
	cin>>eve_trova;

	while (aaa.compare(altro_eve.c_str())){ // loop lettura file 
		
		fileIn1>>dummy>>dummy>>n_samples1; 	 
		fileIn1>>dummy>>board_id1;
		fileIn1>>dummy>>channel1;
		fileIn1>>dummy>>dummy>>event_n1;
		fileIn1>>dummy>>dummy;
		fileIn1>>dummy>>dummy>>dummy>>time1;
		fileIn1>>dummy>>dummy>>dummy>>dummy;
   
		cout<<" event = "<<event_n1<<"  time - "<<time1<<endl;
		
		char title_hist1[20];
		sprintf(title_hist1,"channel %d - event %d ",channel1,event_n1); 
		
		TH1F* h1 = new TH1F("h1",title_hist1,n_samples1,0,n_samples1-1);

		for (int j=0; j<n_samples1; j++){   
         fileIn1>>dato1;				 		
			h1->Fill(j,dato1);	
   	}
		
		if (fileIn1.eof()) break;
   		 
 		if (eve_trova == event_n1){
			h1->SetLineColor(1);
			h1->GetYaxis()->SetRangeUser(ymin,ymax);
			h1->Draw();			
			c1->Update();
 			cout<<"n = next event ; a = other specific event  ;  q = quit"<<endl;
			cin>>aaa;
		}

		if(!aaa.compare(quit_tasto.c_str())) break; //exit the inner while loop
		
		if(!aaa.compare(prossimo.c_str())) eve_trova=eve_trova+1;
		
		//if aaa is "a", exit this while loop and goes to selected event in the outer loop
		
		delete h1;

	} // end loop file

} // 

}//end main
