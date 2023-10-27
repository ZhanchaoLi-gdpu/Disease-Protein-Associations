% select negative samples with SVDD 
% clear
% clc

load HumanProteinDiseaseAssociationNetwork.mat
clear DiseSimeLin HumanProteinDiseaseAssociationMatrix
clear WeiAdjMat

% negative sample
for i = 1 : length(ProtDiseOMIMMeSHAsso)
    PostProtDisePair{i,1} = [ProtDiseOMIMMeSHAsso{i,1} '-' ProtDiseOMIMMeSHAsso{i,5}(6:end)];
end
dise = MeSHName(IndicesDiseProt,:);
prot = HumanProteID(IndicesProtDise,:);
Dise = repmat(dise,length(prot),1);
Prot = repmat(prot,length(dise),1);
for i = 1 : length(Dise)
    NegProtDisePair{i,1} = [Prot{i,1} '-' Dise{i,1}];
end
[NegProtDisePair,ia] = setdiff(NegProtDisePair,PostProtDisePair);
clear i dise prot Dise Prot ia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load DeepWalk_embedding128.mat

% positive sample feature
PosFeat = zeros(length(ProtDiseOMIMMeSHAsso),2*size(DeepWalk_embedding,2));
for i = 1 : length(ProtDiseOMIMMeSHAsso)
    Pline = find(strcmp(HumanProteID,ProtDiseOMIMMeSHAsso{i,1}));
    Dline = find(strcmp(MeSHName,ProtDiseOMIMMeSHAsso{i,5}(6:end)));
    PosFeat(i,:) = [DeepWalk_embedding(Pline,:) DeepWalk_embedding(length(HumanProteID)+Dline,:)];
end
clear i Pline Dline 
PosAimout = ones(size(PosFeat,1),1);
% [PosFeat,ps] = mapminmax(PosFeat',0,1);
% PosFeat = PosFeat';

% Negative sample feature
NegFeat = zeros(length(NegProtDisePair),2*size(DeepWalk_embedding,2));
for i = 1 : length(NegProtDisePair)
    temp = regexp(NegProtDisePair{i,1},'-','split');
    Pline = find(strcmp(HumanProteID,temp{1,1}));
    Dline = find(strcmp(MeSHName,temp{1,2}));
    NegFeat(i,:) = [DeepWalk_embedding(Pline,:) DeepWalk_embedding(length(HumanProteID)+Dline,:)];
end
clear i Pline Dline temp
NegAimout = -1*ones(size(NegFeat,1),1);

C = 2.^(-5 : 2 : 15);
G = 2.^(3 : -2 : -15);
Indices = crossvalind('Kfold',PosAimout,10);
for i = 1 : length(C)
    for j = 1 : length(G)
        % parameter setting
         cost = C(i);
         kernel = Kernel('type', 'gaussian', 'gamma', G(j));
         svddParameter = struct('cost', cost,'kernelFunc', kernel);

        for KK = 1 : 10
            test = (Indices == KK);
            train = ~test;

            Dataxl = PosFeat(train,:);
            Aimoutxl = PosAimout(train,:);

            Datayc = PosFeat(test,:);
            Aimoutyc = PosAimout(test,:);
            
            % creat an SVDD object
            svdd = BaseSVDD(svddParameter);
            % train SVDD model
            svdd.train(Dataxl, Aimoutxl);
            % test SVDD model
            results = svdd.test(Datayc, Aimoutyc);
            Y_hat(test,:) = results.predictedLabel;
        end
        ConMat = confusionmat(PosAimout,Y_hat,'order',[1,-1]);
        TP = ConMat(1,1);
        TN = ConMat(2,2);
        FN = ConMat(1,2);
        FP = ConMat(2,1);
        Acc(i,j) = (TP + TN)/(TP+FN+TN+FP);
        Sen(i,j) = TP/(TP+FN);
        Spe(i,j) = TN/(TN+FP);
        Pre(i,j) = TP/(TP + FP);
        Mcc(i,j) = (TP*TN-FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP));
    end

end

[i,j] = find(max(max(Sen)) == Sen);
% parameter setting
cost = C(i(1));
kernel = Kernel('type', 'gaussian', 'gamma', G(j(1)));
svddParameter = struct('cost', cost,'kernelFunc', kernel);
% creat an SVDD object
svdd = BaseSVDD(svddParameter);
% train SVDD model
svdd.train(PosFeat, PosAimout);
% test SVDD model
for i = 1 : length(NegFeat)
    results = svdd.test(NegFeat(i,:), NegAimout(i,1));
    Y_Hat(i,1) = results.predictedLabel;
    Y_Distance(i,1) = results.distance;
end

Label = (Y_Distance >= 1.0*svdd.radius);
Label = find(Label);

[Test1,NegLabel] = crossvalind('LeaveMOut',length(Label),length(PostProtDisePair)); 
NegLabel = find(NegLabel);

SelNegFeat = zeros(length(NegLabel),2*size(DeepWalk_embedding,2));
for i = 1 : length(NegLabel)
    temp = regexp(NegProtDisePair{Label(NegLabel(i)),1},'-','split');
    Pline = find(strcmp(HumanProteID,temp{1,1}));
    Dline = find(strcmp(MeSHName,temp{1,2}));
    SelNegFeat(i,:) = [DeepWalk_embedding(Pline,:) DeepWalk_embedding(length(HumanProteID)+Dline,:)];
end
clear i Pline Dline temp

Data = [PosFeat;SelNegFeat];
Aimout = [ones(length(ProtDiseOMIMMeSHAsso),1);-1*ones(length(NegLabel),1)];
Indices = crossvalind('Kfold',Aimout,10);

Ntree = 500; %: 100 : 1000;
Fnum = ceil(sqrt(size(Data,2))); % [2.^(1:7)  ceil(sqrt(size(Feature,2)))];%ceil(sqrt(size(Feature,2))); %ceil(size(NetTopFea,2)/3);   %RF鐗瑰緛鏁扮洰榛樿鍊兼槸鎬荤壒寰佹暟鐩殑寮?鏍瑰彿
for i = 1 : length(Ntree)
    for j = 1 : length(Fnum)
        for KK = 1 : 10
            test = (Indices == KK);
            train = ~test;

            Dataxl = Data(train,:);
            Aimoutxl = Aimout(train,:);

            Datayc = Data(test,:);
            Aimoutyc = Aimout(test,:);

            model = TreeBagger(Ntree(i),Dataxl,Aimoutxl,'Method','classification','NumPredictorsToSample',Fnum(j));
            [RFPred,RFvotes] = predict(model,Datayc);
            RFY_hat(test,:) = RFPred;
            RFVotes(test,:) = RFvotes;
        end
        RFY_Hat = str2double(RFY_hat);
        ConMat = confusionmat(Aimout,RFY_Hat,'order',[1,-1]);
        TP = ConMat(1,1);
        TN = ConMat(2,2);
        FN = ConMat(1,2);
        FP = ConMat(2,1);
        RFAcc(i,j) = (TP + TN)/(TP+FN+TN+FP);
        RFSen(i,j) = TP/(TP+FN);
        RFSpe(i,j) = TN/(TN+FP);
        RFPre(i,j) = TP/(TP + FP);
        RFMcc(i,j) = (TP*TN-FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP));
        [RFTXR,RFTYR,RFTTR,RFTAUCR] = perfcurve(Aimout,RFVotes(:,2),1);
        [RFTXP,RFTYP,RFTTP,RFTAUCP] = perfcurve(Aimout,RFVotes(:,2),1,'xcrit','reca','ycrit','prec');
        RFAUCR(i,j) = RFTAUCR;
        RFAUCP(i,j) = RFTAUCP;
    end
end

