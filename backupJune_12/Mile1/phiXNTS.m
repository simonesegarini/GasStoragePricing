function values = phiXNTS(us, params)
%   ADD COMMENTS


alpha = params(1); sigma = params(3);
k = params(4); theta = params(5);


values = (1-alpha)/(k*alpha)*(1 - (1 - 1i.*k/(1-alpha)*(theta.*us + 1i.*us.^2.*sigma^2/2)).^alpha);

end