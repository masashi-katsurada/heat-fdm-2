function a=testptrilu(n)
	a=zeros(n,n);
	for i=1:n
	 a(i,i)=4;
	end
	for i=2:n
	 a(i,i-1)=-1;
	end
	for i=1:n-1
	 a(i,i+1)=-1;
	end
	a(1,n)=-1;
	a(n,1)=-1;
