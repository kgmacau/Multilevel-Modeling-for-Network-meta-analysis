/*****************************************************************************************/
/*Part 1: Multilevel network meta-analysis of diabetes data with hba1c and fpg********/
/*****************************************************************************************/

proc import out=diabetes
  datafile="U:\scans\diabetes.csv"
  dbms=csv replace;
run;


data fpghba1c;
  length drug outvar $100 mean sd n 8;
  set diabetes;
  drug=drug1; mean=hba1c_mean1; sd=hba1c_sd1; n=hba1c_n1; outvar="hbA1c"; output;
  drug=drug2; mean=hba1c_mean2; sd=hba1c_sd2; n=hba1c_n2; outvar="hbA1c"; output;
  drug=drug1; mean=fpg_mean1;   sd=fpg_sd1;   n=fpg_n1;   outvar="FPG"; output;
  drug=drug2; mean=fpg_mean2;   sd=fpg_sd2;   n=fpg_n2;   outvar="FPG"; output;
  rename drug1=drug1c drug2=drug2c;
run;


%macro data_prepare(indata=, outdata=);

data &outdata._1 ;
  set &indata.;
  if nmiss(mean, sd, n)=3  then delete;
run;

proc sort data=&outdata._1 out=&outdata._2 nodupkey; by author outvar year  drug mean sd n; run;


data design;
  set &outdata._2;
  length design $200;
  retain design;
  by author outvar year  drug;
  if first.year then design=drug;
  else design=catx("-", design, drug);

  if last.year;
  keep author outvar year design;
run;

data &outdata._3;
  merge &outdata._2(in=in1) design;
  by author outvar year;
  if in1;
run;
proc sort data=&outdata._3; by design author year; run;


data &outdata.;
  set &outdata._3;
  length trial outcome es se var weight drug1 - drug11 8;
  retain trial 0;
  by design author year;
  if first.year then trial+1;
  else trial=trial;

  outcome=_n_;
  if ~missing(mean) then es=mean;
  if ~missing(sd) then se=sd/sqrt(n);
  if ~missing(se) then do; var=se**2; weight=1/var; end;

  array drugarm(11) drug1 - drug11;
  do i=1 to 11;
     drugarm(i)=0;
  end;
  if drug="pioglitazone" then drug1=1;
  else if drug="rosiglitazone" then drug2=1;
  else if drug="glimepiride" then drug3=1;
  else if drug= "acarbose" then drug4=1;
  else if drug= "vildagliptin" then drug5=1;
  else if drug="metformin" then drug6=1;
  else if drug= "gliclazide" then drug7=1;
  else if drug="sitagliptin" then drug8=1;
  else if drug="glyburide" then drug9=1;
  else if drug="rapaglinide" then drug10=1;
  else if drug="nateglinide" then drug11=1;
  keep author year drug drug1 drug2 drug3 drug4 drug5 drug6 drug7 drug8 drug9 drug10 drug11 mean sd n var es weight design trial outcome outvar;
run;
%mend ;

%data_prepare(indata=fpghba1c, outdata=analysis);

proc export data=analysis outfile="U:\DATA\fpg-hba1c.csv"
     dbms=csv replace;
run;

proc mixed data=analysis method=reml;
  class design trial outcome;
  model es= drug1 drug2 drug3 drug4 drug5 drug6 drug7 drug8 drug9 drug10 drug11/ solution ddfm=satterthwaite;
  weight weight;
  random intercept /sub=design;
  random intercept /sub=trial;
  random intercept /sub=outcome;
  parms 1 1 1 1  /hold=4;
run;
 

/*****************************************************************************************/
/*Part 2: Multilevel network meta-analysis of diabetes data with FPG              ********/
/*****************************************************************************************/

data fpg;
  length drug outvar $100 mean sd n 8;
  set diabetes;
  drug=drug1; mean=fpg_mean1;   sd=fpg_sd1;   n=fpg_n1;   outvar="FPG"; output;
  drug=drug2; mean=fpg_mean2;   sd=fpg_sd2;   n=fpg_n2;   outvar="FPG"; output;
  rename drug1=drug1c drug2=drug2c;
run;


%data_prepare(indata=fpg, outdata=analysis);

proc mixed data=analysis method=reml;
  class design trial;
  model es= drug1 drug2 drug3 drug4 drug5 drug6 drug7 drug8 drug9 drug10 drug11/ solution ddfm=satterthwaite;
  weight weight;
  random intercept /sub=design;
  random intercept /sub=trial;
  parms 1 1 1  /hold=3;
run;


proc export data=analysis outfile="U:\DATA\fpg-3level.csv"
     dbms=csv replace;
run;


/*****************************************************************************************/
/*Part 3: Multilevel network meta-analysis of diabetes data with hba1c            ********/
/*****************************************************************************************/

data hba1c;
  length drug outvar $100 mean sd n 8;
  set diabetes;
  drug=drug1; mean=hba1c_mean1; sd=hba1c_sd1; n=hba1c_n1; outvar="hba1c"; output;
  drug=drug2; mean=hba1c_mean2; sd=hba1c_sd2; n=hba1c_n2; outvar="hba1c"; output;
  rename drug1=drug1c drug2=drug2c;
run;


%data_prepare(indata=hba1c, outdata=analysis);


proc mixed data=analysis method=reml;
  class design trial;
  model es= drug1 drug2 drug3 drug4 drug5 drug6 drug7 drug8 drug9 drug10 drug11/ solution ddfm=satterthwaite;
  weight weight;
  random intercept /sub=design;
  random intercept /sub=trial;
  parms 1 1 1  /hold=3;
run;


proc export data=analysis outfile="U:\DATA\hba1c-3level.csv"
     dbms=csv replace;
run;


/*****************************************************************************************/
/*Part 4: Lumley LME network meta-analysis of diabetes data with FPG              ********/
/*****************************************************************************************/

data fpglme;
  set diabetes;
  rename drug1=drug1c drug2=drug2c;
run;

data fpglme1;
  set fpglme;
  length placebo drug1 - drug11 es se var weight 8 trtpair $40;
  if nmiss(fpg_mean1, fpg_sd1, fpg_n1, fpg_mean2, fpg_sd2, fpg_n2)=6 then delete;

  array drugarm(12) placebo drug1 - drug11;
  do i=1 to 12;
      if missing(drugarm(i)) then drugarm(i)=0;
  end;
  
  se=sqrt(fpg_sd1**2/fpg_n1 + fpg_sd2**2/fpg_n2);
  var=se**2;
  weight=1/var;

  if drug1c<drug2c then do;
     trtpair=catx("-", drug1c , drug2c);
	 es=(fpg_mean1-fpg_mean2);
	 if drug1c="placebo" then placebo=1;
     else if drug1c="pioglitazone" then drug1=1;
     else if drug1c="rosiglitazone" then drug2=1;
     else if drug1c="glimepiride" then drug3=1;
     else if drug1c= "acarbose" then drug4=1;
     else if drug1c= "vildagliptin" then drug5=1;
     else if drug1c="metformin" then drug6=1;
     else if drug1c= "gliclazide" then drug7=1;
     else if drug1c="sitagliptin" then drug8=1;
     else if drug1c="glyburide" then drug9=1;
     else if drug1c="rapaglinide" then drug10=1;
     else if drug1c="nateglinide" then drug11=1;

     if drug2c="placebo" then placebo=-1;
     else if drug2c="pioglitazone" then drug1=-1;
     else if drug2c="rosiglitazone" then drug2=-1;
     else if drug2c="glimepiride" then drug3=-1;
     else if drug2c= "acarbose" then drug4=-1;
     else if drug2c= "vildagliptin" then drug5=-1;
     else if drug2c="metformin" then drug6=-1;
     else if drug2c= "gliclazide" then drug7=-1;
     else if drug2c="sitagliptin" then drug8=-1;
     else if drug2c="glyburide" then drug9=-1;
     else if drug2c="rapaglinide" then drug10=-1;
     else if drug2c="nateglinide" then drug11=-1;
  end;

  else if drug1c>drug2c then do;
     trtpair=catx("-", drug2c , drug1c);
	 es=(fpg_mean2-fpg_mean1);
	 if drug1c="placebo" then placebo=-1;
     else if drug1c="pioglitazone" then drug1=-1;
     else if drug1c="rosiglitazone" then drug2=-1;
     else if drug1c="glimepiride" then drug3=-1;
     else if drug1c= "acarbose" then drug4=-1;
     else if drug1c= "vildagliptin" then drug5=-1;
     else if drug1c="metformin" then drug6=-1;
     else if drug1c= "gliclazide" then drug7=-1;
     else if drug1c="sitagliptin" then drug8=-1;
     else if drug1c="glyburide" then drug9=-1;
     else if drug1c="rapaglinide" then drug10=-1;
     else if drug1c="nateglinide" then drug11=-1;

     if drug2c="placebo" then placebo=1;
     else if drug2c="pioglitazone" then drug1=1;
     else if drug2c="rosiglitazone" then drug2=1;
     else if drug2c="glimepiride" then drug3=1;
     else if drug2c= "acarbose" then drug4=1;
     else if drug2c= "vildagliptin" then drug5=1;
     else if drug2c="metformin" then drug6=1;
     else if drug2c= "gliclazide" then drug7=1;
     else if drug2c="sitagliptin" then drug8=1;
     else if drug2c="glyburide" then drug9=1;
     else if drug2c="rapaglinide" then drug10=1;
     else if drug2c="nateglinide" then drug11=1;
  end;

  drop i;
run;

proc mixed data=fpglme1 method=reml;
  class trtpair;
  model es= drug1 drug2 drug3 drug4 drug5 drug6 drug7 drug8 drug9 drug10 drug11 / solution noint;
  weight weight;
  random intercept/sub=trtpair;
  parms 1 1 /hold=2;
run;


proc export data=fpglme1 outfile="U:\DATA\fpg-lme.csv"
     dbms=csv replace;
run;




/*****************************************************************************************/
/*Part 5: Lumley LME network meta-analysis of diabetes data with hba1c            ********/
/*****************************************************************************************/

data hba1clme;
  set diabetes;
  rename drug1=drug1c drug2=drug2c;
run;

data hba1clme1;
  set hba1clme;
  length placebo drug1 - drug11 es se var weight 8 trtpair $40;
  if nmiss(hba1c_mean1, hba1c_sd1, hba1c_n1, hba1c_mean2, hba1c_sd2, hba1c_n2)=6 then delete;

  array drugarm(12) placebo drug1 - drug11;
  do i=1 to 12;
      if missing(drugarm(i)) then drugarm(i)=0;
  end;

  se=sqrt(hba1c_sd1**2/hba1c_n1 + hba1c_sd2**2/hba1c_n2);
  var=se**2;
  weight=1/var;
  
    if drug1c<drug2c then do;
     trtpair=catx("-", drug1c , drug2c);
	 es=(hba1c_mean1-hba1c_mean2);
	 if drug1c="placebo" then placebo=1;
     else if drug1c="pioglitazone" then drug1=1;
     else if drug1c="rosiglitazone" then drug2=1;
     else if drug1c="glimepiride" then drug3=1;
     else if drug1c= "acarbose" then drug4=1;
     else if drug1c= "vildagliptin" then drug5=1;
     else if drug1c="metformin" then drug6=1;
     else if drug1c= "gliclazide" then drug7=1;
     else if drug1c="sitagliptin" then drug8=1;
     else if drug1c="glyburide" then drug9=1;
     else if drug1c="rapaglinide" then drug10=1;
     else if drug1c="nateglinide" then drug11=1;

     if drug2c="placebo" then placebo=-1;
     else if drug2c="pioglitazone" then drug1=-1;
     else if drug2c="rosiglitazone" then drug2=-1;
     else if drug2c="glimepiride" then drug3=-1;
     else if drug2c= "acarbose" then drug4=-1;
     else if drug2c= "vildagliptin" then drug5=-1;
     else if drug2c="metformin" then drug6=-1;
     else if drug2c= "gliclazide" then drug7=-1;
     else if drug2c="sitagliptin" then drug8=-1;
     else if drug2c="glyburide" then drug9=-1;
     else if drug2c="rapaglinide" then drug10=-1;
     else if drug2c="nateglinide" then drug11=-1;
  end;

  else if drug1c>drug2c then do;
     trtpair=catx("-", drug2c , drug1c);
	 es=(hba1c_mean2-hba1c_mean1);
	 if drug1c="placebo" then placebo=-1;
     else if drug1c="pioglitazone" then drug1=-1;
     else if drug1c="rosiglitazone" then drug2=-1;
     else if drug1c="glimepiride" then drug3=-1;
     else if drug1c= "acarbose" then drug4=-1;
     else if drug1c= "vildagliptin" then drug5=-1;
     else if drug1c="metformin" then drug6=-1;
     else if drug1c= "gliclazide" then drug7=-1;
     else if drug1c="sitagliptin" then drug8=-1;
     else if drug1c="glyburide" then drug9=-1;
     else if drug1c="rapaglinide" then drug10=-1;
     else if drug1c="nateglinide" then drug11=-1;

     if drug2c="placebo" then placebo=1;
     else if drug2c="pioglitazone" then drug1=1;
     else if drug2c="rosiglitazone" then drug2=1;
     else if drug2c="glimepiride" then drug3=1;
     else if drug2c= "acarbose" then drug4=1;
     else if drug2c= "vildagliptin" then drug5=1;
     else if drug2c="metformin" then drug6=1;
     else if drug2c= "gliclazide" then drug7=1;
     else if drug2c="sitagliptin" then drug8=1;
     else if drug2c="glyburide" then drug9=1;
     else if drug2c="rapaglinide" then drug10=1;
     else if drug2c="nateglinide" then drug11=1;
  end;


  drop i;
run;

proc mixed data=hba1clme1 method=reml;
  class trtpair;
  model es= drug1 drug2 drug3 drug4 drug5 drug6 drug7 drug8 drug9 drug10 drug11 / solution noint;
  weight weight;
  random intercept/sub=trtpair;
  parms 1 1 /hold=2;
run;


proc export data=hba1clme1 outfile="U:\DATA\hba1c-lme.csv"
     dbms=csv replace;
run;
