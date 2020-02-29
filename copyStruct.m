function opt2 = copyStruct(opt1,opt2,varargin)
% copy fields from structure opt1 to structure opt2
% if addnew=1, then add additional fields into opt2 from opt1

addnew = 0;

for vc = 1:2:length(varargin)
    switch (varargin{vc})
        case('addnew')
            addnew = varargin{vc+1};
    end
end

inputopt = fieldnames(opt1);
for c = 1:length(inputopt)
    s = inputopt(c); s=s{1};
    if (addnew || isfield(opt2,s))
        opt2.(s) = opt1.(s);
    end
end


end