% retrieve protein-disease association
clear
clc

Obj = BioIndexedFile('FLAT','uniprot_sprot_human.dat','EntryDelimiter','//');

Count = 0;
for i = 1 : Obj.NumEntries
    Entries = getEntryByIndex(Obj,i);
    Entries = strread(Entries,'%s','delimiter','\n','whitespace','');

    Ln = 1;
    [tempid,Ln] = extractfield(Entries,Ln,'ID   ');
    tempid = regexp(tempid,'\w{1,}_HUMAN','match');
    
    [temp,Ln] = extractfield(Entries,Ln,'AC   ');
    temp = regexp(cellstr(temp),'\w{6}','match');
    count = 0;
    tempac = [];
    for j = 1 : length(temp)
        for k = 1 : length(temp{j,:})
            count = count + 1;
            tempac{1,count} = temp{j,1}{1,k};
        end
    end
        
    LineDise = find(strncmp(Entries,'CC   -!- DISEASE:',17));
    tempmim = [];
    if  ~isempty(LineDise)
        count = 0;
        for j = 1 : length(LineDise)
            Ln = LineDise(j) + 1;     
            while matchstart(Entries{Ln,1},'CC       ')
                  Ln = Ln + 1;
            end
            Ln = Ln - 1;
            
            Temp = regexp(Entries(LineDise(j) : Ln,1),'\[MIM\:\d{6}\]','match');
            line = find(~cellfun(@isempty,Temp));
            if  ~isempty(line)
                for k = 1 : length(line)
                    count = count + 1;
                    tempmim{1,count} = Temp{line(k),1}{1,1}(2:end-1);
                end
            else
                
            end
            
        end
         
    else
        
    end
    
    tempmim = unique(tempmim);
    if  ~isempty(tempmim)
        Count = Count + 1;
        ProtID{Count,1} = tempid{1,1};
        ProtAC{Count,1} = tempac;
        ProtDiseOMIM{Count,1} = tempmim;
    else
        
    end
        
end
clear count Count Entries i j k line LineDise Ln Obj
clear temp Temp tempac tempid tempmim

count = 0;
for i = 1 : length(ProtDiseOMIM)
    for j = 1 : length(ProtDiseOMIM{i,1})
        count = count + 1;
        OMIM{count,1} = ProtDiseOMIM{i,1}{1,j};
    end
end
clear count i j 
OMIM = unique(OMIM);

ProtDiseAssoMatr = zeros(length(ProtID),length(OMIM));
for i = 1 : length(ProtDiseOMIM)
    for j = 1 : length(ProtDiseOMIM{i,1})
        LineDise = find(strcmp(OMIM,ProtDiseOMIM{i,1}{1,j}));
        ProtDiseAssoMatr(i,LineDise) = 1;
    end
end
clear i j LineDise
    
save UniprotProtDiseAssoResults 

function tf = matchstart(string,pattern)
%MATCHES start of string with pattern

tf = ~isempty(regexp(string,['^',pattern],'once'));

function [data, outLine] = extractfield(embltext,ln,lineID)
%Extracts a field from the embltext cellstr

startLn = ln;
while matchstart(embltext{ln},lineID)
    ln=ln+1;
end
data = char(embltext{startLn:ln-1});
data = data(:,6:end);
outLine = ln;






