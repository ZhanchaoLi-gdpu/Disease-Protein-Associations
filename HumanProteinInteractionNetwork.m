% Construction of human protein-protein interaction netwrok based on the 
% data in the HIPPIE database
clear
clc

fid = fopen('HIPPIE-current.mitab.txt'); 
C = textscan(fid,'%*s %*s %s %s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %s %*s %*s %*s','delimiter','\t','HeaderLines',1);
clear fid

% delete interactions without an ID in teh UniprotKB database 
Line = find(~cellfun(@isempty,regexp(C{1,1},'uniprotkb:\w{1,}_HUMAN','match')));
D{1,1} = C{1,1}(Line,1);
D{1,2} = C{1,2}(Line,1);
D{1,3} = C{1,3}(Line,1);
Line = find(~cellfun(@isempty,regexp(D{1,2},'uniprotkb:\w{1,}_HUMAN','match')));
E{1,1} = D{1,1}(Line,1);
E{1,2} = D{1,2}(Line,1);
E{1,3} = D{1,3}(Line,1);
C = E;
clear D E Line

% remove interactions with score = 0
score = str2double(C{1,3});
Indices = (score ~= 0);
C{1,1} = C{1,1}(Indices,1);
C{1,2} = C{1,2}(Indices,1);
C{1,3} = C{1,3}(Indices,1);
clear Indices score

% abolish self-interactions 
Line = strcmp(C{1,1},C{1,2});
Line = ~Line;
C{1,1} = C{1,1}(Line,1);
C{1,2} = C{1,2}(Line,1);
C{1,3} = C{1,3}(Line,1);
clear Line

% delete duplicate interactions 
for i = 1 : length(C{1,1})
    D{1,1} = C{1,1}{i,1};
    D{1,2} = C{1,2}{i,1}; 
    D = sort(D);
    PPI{i,1} = [D{1,1} ' ' D{1,2} ' ' C{1,3}{i,1}];
end
PPI = unique(PPI);
clear i D C

% reorganization of collected protein-protein interaction data
PP = regexp(PPI,' ','split');
for i = 1 : length(PPI)
    C{i,1} = PP{i,1}{1,1};
    C{i,2} = PP{i,1}{1,2};
    C{i,3} = PP{i,1}{1,3};
end
clear i PP

for i = 1 : length(C)
    C{i,1} = C{i,1}(11:end);
    C{i,2} = C{i,2}(11:end);
end

% Collect sequence information of the proteins in the network, 
% if there is no sequence information, then further remove 
% the protein in the constructed network.
[H, S] = fastaread('uniprot_sprot.fasta');

ID = unique([C(:,1);C(:,2)]);
count = 0;
delt = 0;
for i = 1 : length(ID)
    Line = find(~cellfun(@isempty,strfind(H,ID{i,1}))); 
    if  ~isempty(Line)
        count = count + 1;
        ProtID{count,1} = ID{i,1};
        ProtSeq{count,1} = S{1,Line(1)};
    else
        delt = delt + 1;
        DeltID{delt,1} = ID{i,1};
    end
end
clear i count H S Line delt ID

for i = 1 : length(C)
    PPI{i,1} = [C{i,1} ' ' C{i,2} ' ' C{i,3}];
end
clear i 

% remove protein-protein interaction pairs, where a protein has no 
% sequence information
Line = [];
for i = 1 : length(DeltID)
    Lin = find(~cellfun(@isempty,regexp(PPI,['\<' DeltID{i,1} '\>'],'match'))); 
    Line = [Lin;Line];
end
Line = sort(Line);
clear Lin i DeltID
Indices = ones(length(PPI),1);
Indices(Line,1) = 0;
PPI = PPI(logical(Indices),:);
clear Indices Line C

% reorganization of collected protein-protein interaction data
PP = regexp(PPI,' ','split');
for i = 1 : length(PPI)
    C{i,1} = PP{i,1}{1,1};
    C{i,2} = PP{i,1}{1,2};
    C{i,3} = PP{i,1}{1,3};
end
clear i PP
ID = unique([C(:,1);C(:,2)]);
for i = 1 : length(ID)
    Line = find(strcmp(ProtID,ID{i,1}));
    Seq{i,1} = ProtSeq{Line,1};
end
clear i Line
ProtID = ID;
ProtSeq = Seq;
clear ID Seq

% generate an adjacency matrix
AdjMat = zeros(length(ProtID));
for i = 1 : length(C)
    LineA = find(strcmp(ProtID,C{i,1}));
    LineB = find(strcmp(ProtID,C{i,2}));
    AdjMat(LineA,LineB) = 1;
    AdjMat(LineB,LineA) = 1;
end
clear i LineA LineB
AdjMat = sparse(AdjMat);

% generate an adjacency matrix with edge weighted (i.e. interaction confidence score)
WeiAdjMat = sparse(zeros(length(ProtID)));
for i = 1 : length(C)
    LineA = find(strcmp(ProtID,C{i,1}));
    LineB = find(strcmp(ProtID,C{i,2}));
    WeiAdjMat(LineA,LineB) = str2num(C{i,3});
    WeiAdjMat(LineB,LineA) = str2num(C{i,3});
end
clear i LineA LineB

save HumanPPI.mat













