//void display_traces(string file_path,string file_num)

void display_traces_1ch(string file_1){
gROOT->Reset();

ifstream fileIn1 (file_1.c_str());

char dummy[20];
int n_samples1=0;
int event_n1, eve_trova;
double data;

string quit_tasto="q";
string prossimo="n";
string altro_eve="a";

string aaa="z"; 

int ymin = 80;
int ymax = 120;
//int ymin = 870;
//int ymax = 910;
//int ymin = 50;
//int ymax = 300;

TCanvas *c1=new TCanvas ("EVENTS"," ",0,0,1000,600);

while(aaa.compare(quit_tasto.c_str())){
	aaa="z";
	cout<<"Evento: "<<endl;
	cin>>eve_trova;

	while (aaa.compare(altro_eve.c_str())){ // loop lettura file 
		
		fileIn1>>dummy>>dummy>>n_samples1; 	 
		fileIn1>>dummy>>dummy;
		fileIn1>>dummy>>dummy;
		fileIn1>>dummy>>dummy>>event_n1;
		fileIn1>>dummy>>dummy;
		fileIn1>>dummy>>dummy>>dummy>>dummy;
		fileIn1>>dummy>>dummy>>dummy>>dummy;
   
		cout<<" event: "<<event_n1<<" , eve_trova: "<<eve_trova<<endl;
				
		if (fileIn1.eof()){ 
			aaa="q";
			break;
		}
		
		char title_hist1[20];
		sprintf(title_hist1,"event %d ",event_n1); 
   	TH1F* h1 = new TH1F("h1",title_hist1,n_samples1,0,n_samples1-1);		
		for (int j=0; j<n_samples1; j++){   
         fileIn1>>data;				 		
			h1->Fill(j,data);	
   	}	 
		
 		if (eve_trova == event_n1){
			h1->SetLineColor(1);
			h1->GetYaxis()->SetRangeUser(ymin,ymax);
			h1->Draw();			
			c1->Update();
 			cout<<"n = next event ; a = other specific event  ;  q = quit"<<endl;
			cin>>aaa;
		}
		else{
			if(eve_trova < event_n1) eve_trova = event_n1+1;
		}
	
		if(!aaa.compare(quit_tasto.c_str())) break; //exit the inner while loop
		
		if(!aaa.compare(prossimo.c_str()) && eve_trova <= event_n1) eve_trova++;
		
		delete h1;

	} // end loop file

} // 

}//end main
