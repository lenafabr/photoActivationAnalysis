function bound = rect2bound(rect)
% from rectangle representation to boundary coordinates
% does not include overlap of last coordinate

bound = [rect(1) rect(2); rect(1)+rect(3) rect(2); rect(1)+rect(3) rect(2)+rect(4); rect(1) rect(2)+rect(4)];

end