libname Socmob4 odbc noprompt="dsn=CCCC;database=CCCC" schema=CLEAN_2022;
%macro connect;
    connect to odbc(noprompt="dsn=CCCC;database=CCCC");
%mend;
%macro age(date,birth);
	floor ((intck('month',&birth,&date) - (day(&date) < day(&birth))) / 12) 
%mend age;
%macro cts_age(date,birth);
	(((&date) - (&birth)) / 365.25) 
%mend cts_age;
* macros for SQL Server;
%macro sql_generate(start,finish,generate=generate,step=1);
&generate AS (
    SELECT &start AS &generate
    UNION ALL
    SELECT &generate+&step FROM &generate WHERE &generate+&step<=&finish)
%mend;
%macro sql_age(date,birth);
	(DATEDIFF(year, &birth, &date) - CASE WHEN (dateadd(year, datediff (year, &birth, &date), &birth) > &date) THEN 1 ELSE 0 END)
%mend;
%macro least(a,b,c,d,e);
%if "&c"="" %then
	%do; (case when (&a) is null then (&b) when (&b) is null then (&a) when (&a)<=(&b) then (&a) else (&b) end) %end;
%else %do; %least(&a,%least(&b,&c,&d,&e)) %end;
%mend;
%macro greatest(a,b,c,d,e);
%if "&c"="" %then
	%do; (case when (&a) is null then (&b) when (&b) is null then (&a) when (&a)>=(&b) then (&a) else (&b) end) %end;
%else %do; %greatest(&a,%greatest(&b,&c,&d,&e)) %end;
%mend;
%macro sql_cts_age(date,birth);
	(DATEDIFF(day, &birth, &date)/365.25)
%mend;
%macro trailing(expression); %if "&expression"~="" %then &expression,; %mend;
%macro sql_lexis(data,entry,exit,fail,start=0,finish=100,step=10,lhs=lhs,rhs=rhs,breaks=breaks,out=lexis,keep=,pt=pt);
&breaks as (select &start as &lhs, &start+&step as &rhs union all select &lhs+(&step),&lhs+2*(&step) from &breaks where &lhs+2*(&step)<=&finish),
&out as (select %trailing(&keep) &lhs as &lhs, &rhs as &rhs, %greatest(&lhs,&entry) as &entry, %least(&rhs,&exit) as &exit,
	%least(&rhs,&exit) - %greatest(&lhs,&entry) as &pt,
	case when &lhs<=&exit and &exit<=&rhs then &fail else 0 end as &fail
from &data, &breaks
where %greatest(&lhs,&entry)<%least(&rhs,&exit))
%mend;
%let bs = +(-1);

proc sql;
    %connect;

%macro with_statement(eof=datefromparts(2022,12,31),sof=datefromparts(1991,1,1));
	with 
		mid(yr,mid) as (select year(&sof) as yr, datefromparts(year(&sof),7,1) as mid
			union all
			select yr+1, datefromparts(yr+1,7,1) from mid WHERE yr+1<=year(&eof)),
		deaths as (select lopnr, cast(case 
			when substring(death_date,5,4)='0000' then concat(substring(death_date,1,4), '0701') 
			when substring(death_date,7,4)='00' then concat(substring(death_date,1,6), '15') 
			when len(death_date)=6 then concat(death_date, '15') 
			when death_date='19610229' then '19610228'
			else death_date end as date) as dod 
		from clean_2022.death),		
		t1 as (select lopnr, dereg_reas, cast(dereg_date as char(8)) as s from clean_2022.nkc_pop_2022 where dereg_reas is not null), 
		t2 as (select lopnr, dereg_reas, case when substring(s,5,4)='0000' then concat(substring(s,1,4), '0701') 
			when substring(s,7,2)='00' then concat(substring(s,1,6),'15') else s end as s2 from t1),
		dereg as (select lopnr, dereg_reas, try_cast(s2 as date) as ddereg from t2),
		t3(lopnr,s) as (select lopnr, cast(dereg_date as char(8)) as s from clean_2022.nkc_pop_2022 where dereg_reas='UV'), 
		t4(lopnr,s2) as (select lopnr, case when substring(s,5,4) in ('0000','0015') then concat(substring(s,1,4), '0701') 
			when substring(s,7,2)='00' then concat(substring(s,1,6),'15') else s end as s2 from t3),
		old_emigration(lopnr,doem) as (select lopnr, try_cast(s2 as date) as doem from t4),
		emi_imm as (select lopnr, case when emi_imm_date in ('0--') then null 
			when substring(emi_imm_date,6,5)='00-00' then cast(concat(substring(emi_imm_date,1,5), '07-01') as date) 
			when substring(emi_imm_date,9,2)='00' then cast(concat(substring(emi_imm_date,1,8), '15') as date)
			when emi_imm_date='1966-02-29' then cast('1966-02-28' as date)
			when emi_imm_date='1963-11-31' then cast('1963-11-30' as date)
			else cast(emi_imm_date as date) end as dt, emi_imm_code
		from clean_2022.emi_imm),
		old_births as (select lopnr, min(case when birth_date is not null then try_cast(birth_date as date) when birth_yr is not null then datefromparts(birth_yr,7,1) else null end) as dob
			from clean_2022.nkc_person_2022 where lopnr<>4048225 and valid_pnr=1
			group by lopnr),
		births as (select lopnr, cast(concat(birth_date,'15') as date) as dob from clean_2022.birth where sex=2),
		invalid_lopnr(lopnr) as (select lopnr from births where dob is null 
			union select lopnr from deaths where dod is null),
		cytology as (select lopnr, x_sample_date, x_sample_yr, snomed_severity, 
			case when snomed_severity is null then null when snomed_severity in (6,7) then 1 else null end as lg,
			case when snomed_severity is null then null when snomed_severity between 8 and 15 then 1 else null end as hg
		from clean_2022.Nkc_trans_cell_2022 
		where x_sample_date is not null),
		histology as (select lopnr, x_sample_date, x_sample_yr,
			case when topo3 = 'T83' and snomed_translated in ('M80703','M81401','M81403','M85601','M85603','M80413') then 1 else null end as cancer,
			case when topo3 = 'T83' and snomed_translated in ('M74007','M80772','M80702','M80762','M81402','M85602','M80703','M81401','M81403','M85601','M85603','M80413') 
				then 1 else null end as hsil
		from clean_2022.nkc_pad_translated_2022 
		where x_sample_date is not null),
		first_cancer as (select lopnr, min(cast(x_sample_date as date)) as dia_date, min(x_sample_yr) as dia_yr from histology where cancer=1 group by lopnr),
		person as (select b.lopnr, b.dob, d.dod, ddereg, dereg_reas, e.doem, i.doim, &eof as eof, %greatest(dob,doim) as entry, %least(ddereg,dod,doem,&eof) as ext,
			case when dod=%least(ddereg,dod,doem,&eof) then 1 else 0 end as died,
			c.dia_date
		from births as b 
			left join deaths as d on b.lopnr=d.lopnr 
			left join dereg on b.lopnr=dereg.lopnr
			left join (select lopnr, min(dt) as doem from emi_imm where emi_imm_code='Utv' and dt is not null group by lopnr) as e on b.lopnr=e.lopnr
			left join (select lopnr, min(dt) as doim from emi_imm where emi_imm_code='Inv' and dt is not null group by lopnr) as i on b.lopnr=i.lopnr
				and (e.lopnr is null or doim<doem)
			left join first_cancer as c on b.lopnr=c.lopnr
		where b.lopnr not in (select lopnr from invalid_lopnr)),
		events as ((select lopnr, x_sample_date as x_date, null as cancer, null as hsil, 'cyto' as event from cytology) 
			union all (select lopnr, x_sample_date as x_date, cancer, hsil, 'hist' as event from histology))
%mend;

proc sql;
    %connect;
	create table events as select * from connection to odbc (
	%with_statement select p.*, x_date, cancer, hsil, event from person as p left join events as e on p.lopnr=e.lopnr order by p.lopnr, e.x_date
	);
select year(input(dob,e8601da.)) as year, count(*) from events where cancer group by year;
quit;

%macro as_date(var);
if &var~='' then &var.1=input(&var,e8601da.);
format &var.1 e8601da.;
rename &var.1 = &var;
drop &var;
%mend;

data events2;
	set events;
	%as_date(dob);
	%as_date(dod);
	%as_date(ddereg);
	%as_date(doem);
	%as_date(doim);
	%as_date(eof);
	%as_date(entry);
	%as_date(ext);
	%as_date(dia_date);
	if x_date~=. then date=datepart(x_date);
	format date e8601da.;	
	drop x_date;
	if cancer=. then cancer=0;
	if hsil=. then hsil=0;
	cohort=year(dob1);
run;
proc sql; create index cohort on events2; quit;

*   1. \text{No cancer detected and no HSIL detected} && W(t)+X(t)+Y(t) \\
    2. \text{First HSIL detected at $t_j$ with no subsequent detection} && W^*(t|t_j)+X^*(t|t_j)+Y^*(t|t_j)\\
    3. \text{Second HSIL detected with first HSIL at $t_j$} && X^*(t|t_j)(1-\beta_2) \\
    4. \text{Screen-detected cancer with no previous HSIL} && Y(t)(1-\beta_3)\\
    5. \text{Screen-detected cancer with one previous HSIL detected} && Y^*(t|t_j)(1-\beta_3)\\
    6. \text{Interval cancer with no HSIL detected} && I(t) \\
    7. \text{Interval cancer with HSIL detected at time $t_j$} && I^*(t|t_j);

* Combine: persons(lopnr, entry, dob, date=ext, event=eof)  + cyto(lopnr,event=cyto,date) + histo(lopnr, event=histo, date, I(HSIL), I(cancer)) + diagnoses(lopnr, date, event=cancer);

proc sql;
* select sum(cancer) from events2 where %cts_age(entry,dob)<15 and %cts_age(ext,dob)>=15 and cohort=1950; * n=672;
* Repeated tests on the same day;
create table events3 as
select lopnr, date, max(cohort) as cohort, max(entry) as entry, max(dob) as dob, max(ext) as ext, min(event) as event, max(hsil) as hsil, max(event='cyto') as cyto,
	max(event='hist') as hist, max(cancer) as cancer
from events2
where 1960<=cohort<=1960
group by lopnr, date
order by lopnr, date;

* collapse for histo and cyto at close times;
data events4;
set events3;
by lopnr;
array lags[5];
array old[5];
array curr[5] date hsil cancer cyto hist;
do k=1 to 5; lags[k]=lag(curr[k]); end;
retain completed;
if first.lopnr then do; completed=0; end;
else if completed then do; * nothing more; end;
else do;
	* Case: cancer histology within 30 days of a cytology => interval cancer (do not merge);
	if cancer=1 and 0<=(date-lags[1])<=30 then do;
		do k=1 to 5; old[k]=curr[k]; curr[k]=lags[k]; end;
		output;
		do k=1 to 5; curr[k]=old[k]; end;
		completed=1;
		output;
	end;
	* Case: merge cyto and histo within 90 days;
	else if 0<=(date - lags[1])<=90 then do; * merge records;
		* old_date = date;
		do k=1 to 5; curr[k]=max(curr[k],lags[k]); end; * NB: => max(date)!;
		* date = min(lags[1],old_date);
	end; else do; * output *previous* record;
		do k=1 to 5; old[k]=curr[k]; curr[k]=lags[k]; end;
		output;
		do k=1 to 5; curr[k]=old[k]; end;
		if cancer=1 then do;
			output;
			completed=1;
		end;
	end;
end;
if last.lopnr and completed=0 then output;
keep lopnr cohort dob entry ext date cyto hist hsil cancer;
run;

* proc sql;
* create table checks as
select cancers.lopnr, freq3.freq as freq3, freq4.freq as freq4 from
(select distinct lopnr from events3 where cancer=1) as cancers
inner join (select lopnr, count(*) as freq from events3 group by lopnr) as freq3 on cancers.lopnr=freq3.lopnr
inner join (select lopnr, count(*) as freq from events4 group by lopnr) as freq4 on freq3.lopnr=freq4.lopnr and freq3.freq>freq4.freq;

filename temp '\\freja/homedir$/marcle/Downloads/temp-20231217.R';
data temp;
set events4 end=last;
file temp;
if _n_=1 then do;
	put 'screening=list(' @;
end;
* ext = min(ext, dia_date);
* age_entry = %cts_age(entry,dob); * from age 15 years;
age_ext = %cts_age(ext,dob);
age_event = %cts_age(date,dob);
* if dia_date~=. then age_dx = %cts_age(dia_date,dob);
by lopnr;
array ti[200];
array tc [200];
retain ti1-ti200 tc1-tc200 i tj j tj2 j2 nhsil type reported;
if first.lopnr then do;
	* clear arrays;
	do k=1 to max(i,1); ti[k]=.; tc[k]=.; end;
	i=1;           * counter for the output;
	j=.;           * index for the first HSIL;
	tj=.;          * age for the first HSIL;
	j2=.;          * index for the second HSIL;
	tj2=.;         * age for the second HSIL;
	type=1;        * final record type - starts with no HSIL and no cancer;
	reported=0;
	nhsil=0;       * number of previous HSILs;	
end; 
* Case: no history;
if age_event=. and reported=0 then do;
	reported=1;
	put 'list(t=' age_ext &bs ',' cohort= ', type=1, nhsil=0, ti=Inf, j=NA, tj=NA, j2=NA, tj2=NA)'@;
	if last then put ')'; else put ',';
	output;
end;
if .<age_event<=age_ext and type<4 then do;
	ti[i]=age_event;
	tc[i]=cyto;
	if hist=1 and nhsil=0 and cancer=1 then do; type=4; end; * cancer before HSIL (defaults to interval cancer);
	else if hist=1 and nhsil>0 and cancer=1 then do; type=5; end; * cancer after previous HSIL (defaults to interval cancer);
	else if hist=1 and nhsil=0 and hsil=1 then do; nhsil+1; type=2; tj=age_event; j=i; end; * prev HSIL with no other history: ;
	else if hist=1 and nhsil=1 and hsil=1 then do; nhsil+1; type=3; tj2=age_event; j2=i; end; * prev HSIL with no other history: ;
	else if hist=1 and nhsil>1 and hsil=1 then do; nhsil+1; type=3; end; * subsequent HSIL;
	* screen-detected cancer?;
	if hist=1 and cancer=1 and cyto=1 then do;
		if nhsil=0 then type=6; else type=7;
	end;
	else if hist=1 and cancer=1 and i>1 then 
		do k=1 to i-1;
			if tc[k]=1 and 1/12<=(age_event-ti[k])<=6/12 then do; if nhsil=0 then type=6; else type=7; end;
		end;
    i+1; * NB: incremented here!;
end; * else ignore the event;
if reported=0 and (last.lopnr or type>3) then do;
	reported=1;
	put 'list('@;
	if age_event=. or age_ext<age_event then put 't=' age_ext &bs ',' @;
	else if type<=3 and age_ext>age_event then put 't=' age_ext &bs ',' @; else put 't=' age_event &bs ','@;
	put 'cohort= ' cohort ','@;
	put 'type=' type ', nhsil=' nhsil ','@;
	put 'ti=c( '@; do k=1 to i-1; put ti[k]@; if k<i-1 then put ','@; end; put '),'@;
	/*put 'te=c( '@; do k=1 to i-1; put "'" te[k] &bs "'" @; if k<i-1 then put ','@; end; put '),'@;*/
	if j=. then put 'j=NA, tj=NA,'@; else put 'j=' j &bs ', tj=' tj ','@;
	if j2=. then put 'j2=NA, tj2=NA'@; else put 'j2=' j2 &bs ', tj2=' tj2@;
	put ')'@;
	if last then put ')'; else put ',';
	output;
end;
where %cts_age(entry,dob)<15 and %cts_age(ext,dob)>=15;
run;
