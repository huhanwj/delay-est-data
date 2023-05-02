function y=raisedcosine(t,alpha)
y=sinc(t).*cos(alpha*pi*t)./(1-(2*alpha*t).^2);
index=isnan(y)| isinf(y);
y(index)=pi/4.*sinc(1/(2*alpha));
end