GeneontObj = geneont('LIVE', true); %prezemanje na tekovna GO

array_of_terms = GeneontObj.terms; %kolona od termini od tekovna GO, potrebna za manipulacija so konkretni atributi na termini

ontologies = get(array_of_terms,'ontology'); %kolona od vrednosti za atributot 'ontology' za sekoj termin

mask_CC=strcmp(ontologies,'cellular component'); %kolona od logicki vrednosti: 0-praviloto ne e ispolneto, 1-praviloto e ispolneto
mask_MF=strcmp(ontologies,'molecular function');
mask_BP=strcmp(ontologies,'biological process');

GO_CC_terms=GeneontObj.terms(mask_CC);%niza od termini
GO_CC=GeneontObj(GO_CC_terms);%nova ontologija samo za cellular component
GO_CC_obs=get(GO_CC_terms,'obsolete'); %kolona od vrednosti za atributot obsolete, za da mozeme da gi isfrlime onie so vrednost 1
GO_CC_obs=not(cell2mat(GO_CC_obs)); %cell2mat za da dobieme logicki vrednosti i not zatoa sto treba obratna logika
GO_CC_terms=GO_CC.terms(GO_CC_obs); %niza od termini koi ne se obsolete
GO_CC=GO_CC(GO_CC_terms); %nova ontologija za non-obsolete cellular component


GO_MF_terms=GeneontObj.terms(mask_MF);%niza od termini
GO_MF=GeneontObj(GO_MF_terms);%nova ontologija samo za molecular function
GO_MF_obs=get(GO_MF_terms,'obsolete'); %kolona od vrednosti za atributot obsolete, za da mozeme da gi isfrlime onie so vrednost 1
GO_MF_obs=not(cell2mat(GO_MF_obs)); %cell2mat za da dobieme logicki vrednosti i not zatoa sto treba obratna logika
GO_MF_terms=GO_MF.terms(GO_MF_obs); %niza od termini koi ne se obsolete
GO_MF=GO_MF(GO_MF_terms); %nova ontologija za non-obsolete molecular function

GO_BP_terms=GeneontObj.terms(mask_BP);%niza od termini
GO_BP=GeneontObj(GO_BP_terms);%nova ontologija samo za biological process
GO_BP_obs=get(GO_BP_terms,'obsolete'); %kolona od vrednosti za atributot obsolete, za da mozeme da gi isfrlime onie so vrednost 1
GO_BP_obs=not(cell2mat(GO_BP_obs)); %cell2mat za da dobieme logicki vrednosti i not zatoa sto treba obratna logika
GO_BP_terms=GO_BP.terms(GO_BP_obs); %niza od termini koi ne se obsolete
GO_BP=GO_BP(GO_BP_terms); %nova ontologija za non-obsolete cellular component


[CCmatrix,CCid,CCrel]=getmatrix(GO_CC);
[CCi,CCj,CCval]=find(CCmatrix);
CCgo_id=num2goid(CCid);
CCtermi=CCgo_id(CCi);
CCtermj=CCgo_id(CCj);
CCnameval=CCrel(CCval);
CCnameval=CCnameval';
CCfull=[CCtermi,CCtermj,CCnameval];
fidCC1=fopen('CCGOterms.txt','wt');
fidCC2=fopen('CCGOfull.txt','wt');
fprintf(fidCC1,'%s\n',CCgo_id{:});
len=size(CCtermi);
for counter=1:len(1)
  fprintf(fidCC2,'%s %s %s\n',CCfull{counter,2},CCfull{counter,3},CCfull{counter,1});
end
fclose(fidCC1);
fclose(fidCC2);

[MFmatrix,MFid,MFrel]=getmatrix(GO_MF);
[MFi,MFj,MFval]=find(MFmatrix);
MFgo_id=num2goid(MFid);
MFtermi=MFgo_id(MFi);
MFtermj=MFgo_id(MFj);
MFnameval=MFrel(MFval);
MFnameval=MFnameval';
MFfull=[MFtermi,MFtermj,MFnameval];
fidMF1=fopen('MFGOterms.txt','wt');
fidMF2=fopen('MFGOfull.txt','wt');
fprintf(fidMF1,'%s\n',MFgo_id{:});
len=size(MFtermi);
for counter=1:len(1)
  fprintf(fidMF2,'%s %s %s\n',MFfull{counter,2},MFfull{counter,3},MFfull{counter,1});
end
fclose(fidMF1);
fclose(fidMF2);

[BPmatrix,BPid,BPrel]=getmatrix(GO_BP);
[BPi,BPj,BPval]=find(BPmatrix);
BPgo_id=num2goid(BPid);
BPtermi=BPgo_id(BPi);
BPtermj=BPgo_id(BPj);
BPnameval=BPrel(BPval);
BPnameval=BPnameval';
BPfull=[BPtermi,BPtermj,BPnameval];
fidBP1=fopen('BPGOterms.txt','wt');
fidBP2=fopen('BPGOfull.txt','wt');
fprintf(fidBP1,'%s\n',BPgo_id{:});
len=size(BPtermi);
for counter=1:len(1)
  fprintf(fidBP2,'%s %s %s\n',BPfull{counter,2},BPfull{counter,3},BPfull{counter,1});
end
fclose(fidBP1);
fclose(fidBP2);