function structarray = cellofstructs2structarray(cellofstructs)

    tmp = cellofstructs;
    fields = fieldnames(tmp{1});

    clear allmasks
    structarray(length(tmp)) = struct();

    for zi = 1:length(tmp) 
        for fi = 1:length(fields)
            structarray(zi).(fields{fi}) = tmp{zi}.(fields{fi});
        end
    end
end
    