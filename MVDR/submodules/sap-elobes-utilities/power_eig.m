function d = power_eig(rxx)
rxx = rxx/rxx(1);
rxx1 = rxx(:,1);

while 1
	rxx = rxx*rxx;
	rxx = rxx/rxx(1);
    
	if real(rxx(:,1)'*rxx1/(rxx(:,1)'*rxx(:,1)))>0.999
		break;
    end
    rxx1 = rxx(:,1);
%     disp(rxx1)
%     dd=0;
end
d = rxx(:,1);