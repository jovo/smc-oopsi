function report=errorsSMC(P,est,tgt)
%report=errorsSMC(P,estimators,targets)
%Calculates mean error in spike positions and mean fractions of false
%negative and false positive spikes from collection of "targets" and
%their corresponding "estimator" spike trains. Both "targets" and
%"estimators" are cell-arrays containing spike trains. For example, 
%targets=cell-array of 100 logical arrays specifying spike trains;
%estimators=cell-array of 100 logical arrays specifying estimated spike 
%trains corresponding to targets. P is a structure of parameters containing 
%fields P.mismatch_penalty, specifying distance between spikes that are 
%considered mismatched, and P.expected_slack, specifying an upper expected 
%bound on the fraction of mismatched spikes.


%spikes discrepancy penalty 
DM=P.mismatch_penalty^2;
slack_f=P.expected_slack;
v=1; %verbosity

errs=zeros(3,length(est));
for i=1:length(est)
  %prepare spike trains and match costs
  X=find(tgt{i}(:));     %true spike train
  Y=find(est{i}(:));  %estimator spike train
  A=-2*repmat(X,1,length(Y)).*repmat(Y',length(X),1);
  A=A+repmat(X,1,length(Y)).^2;
  A=A+repmat(Y',length(X),1).^2;  %distances matrix for matching

  %add slack-matches: allow each X to be matched out
  slack_n=ceil(slack_f*length(X));  
  A=[A,repmat(DM,length(X),slack_n)];
  slack_n=ceil(slack_f*length(Y));  
  A=[A;repmat(DM,slack_n,size(A,2))];
  
  %prepare optimal assignments
  [a,b]=assignmentoptimal(A);     %Hungarion alg
  a=a(1:length(X)); a(a>length(Y))=0;
  x=(X(a>0)-Y(a(a>0))).^2;       %mismatches
  dn_lost=length(X)-sum(a>0);
  dn_false=length(Y)-sum(a>0);
  errs(1,i)=sqrt(sum(x(x<DM))/sum(x<DM));   %mean spike inaccuracy
  errs(2,i)=(dn_lost)/length(X);  %spikes lost
  errs(3,i)=(dn_false)/length(X); %false spikes
  if(v && mod(i,10)==0) fprintf('\b\b\b%3i',i); end
end
if(v) fprintf('\n'); end

report=[];
report.time_inaccuracy=mean(errs(1,:));
report.false_negatives=mean(errs(2,:));
report.false_positives=mean(errs(3,:));
report.full_data=errs;
