% construct human disease-protein association network
clear
clc

load HIPPIEMaxConnComp.mat
HumanProteID = ProtID;
HumanProteSeq = ProtSeq;
clear ProtID ProtSeq

load DiseSimeLinMaxConSubGrap.mat

load ProtDiseMeSHAssoResult.mat
DiseProtID = ProtID;
clear ProtAC ProtID

% delete protein-disease associaion in which protein unincluded in
% human PPI and diseae unincluded in Mesh similar netwrok
count = 0;
for i = 1 : length(ProtDiseOMIMMeSHAsso)
    Pline = find(strcmp(HumanProteID,ProtDiseOMIMMeSHAsso{i,1}));
    Dline = find(strcmp(MeSHName,ProtDiseOMIMMeSHAsso{i,5}(6:end)));
    if  ~isempty(Pline) && ~isempty(Dline)
        count = count + 1;
        Indices(count,1) = i;
    else

    end
end    
Temp(:,1)  = ProtDiseOMIMMeSHAsso(Indices,1);
Temp(:,2)  = ProtDiseOMIMMeSHAsso(Indices,2);
Temp(:,3)  = ProtDiseOMIMMeSHAsso(Indices,3);
Temp(:,4)  = ProtDiseOMIMMeSHAsso(Indices,4);
Temp(:,5)  = ProtDiseOMIMMeSHAsso(Indices,5);
ProtDiseOMIMMeSHAsso = Temp;
clear count clear Dline i Indices Pline Temp
clear DiseLabel DiseProtID MeSHDiseID OMIMDiseID
clear ProtDiseMeSHAssoMatr UnipDiseID

% construct protein-disease associatin network
ProtDiseAssoMatr = zeros(length(HumanProteID),length(MeSHName));
for i = 1 : length(ProtDiseOMIMMeSHAsso)
    Pline = find(strcmp(HumanProteID,ProtDiseOMIMMeSHAsso{i,1}));
    Dline = find(strcmp(MeSHName,ProtDiseOMIMMeSHAsso{i,5}(6:end)));
    ProtDiseAssoMatr(Pline,Dline) = 1;
end
clear i Pline Dline

% find known protein related to disease
KnowProtDise = unique(ProtDiseOMIMMeSHAsso(:,1));
[temp,ia,ib] = intersect(HumanProteID,KnowProtDise);
IndicesProtDise = ia;
clear temp ia ib KnowProtDise
% find known disease with associated protein
for i = 1 : length(ProtDiseOMIMMeSHAsso(:,5))
    KnowDiseProt{i,1} = ProtDiseOMIMMeSHAsso{i,5}(6:end);
end
KnowDiseProt = unique(KnowDiseProt);
[temp,ia,ib] = intersect(MeSHName,KnowDiseProt);
IndicesDiseProt = ia;
clear i temp ia ib KnowDiseProt

clear AdjMat C PPI ProtDiseAssoMatr HumanProteSeq
DiseSimeLin = sparse(DiseSimeLin);

save HumanProteinDiseaseAssociationNetwork

TriuAdjMat = triu(HumanProteinDiseaseAssociationMatrix);
HumanProteinDiseaseAssociationMatrix(find(HumanProteinDiseaseAssociationMatrix)) = 1;
fid = fopen('HumanPDNetwork.txt','w');
for i = 1 : length(TriuAdjMat) 
    Line = find(TriuAdjMat(i,:));
    for j = 1 : length(Line)
        fprintf(fid,'%d %d %.2f\n',[i-1 Line(j)-1 full(HumanProteinDiseaseAssociationMatrix(i,Line(j)))]);
    end
end
fclose(fid);
clear fid TriuAdjMat i j Line

