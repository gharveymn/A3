function ret=mergesort(arr)
	
	n = max(size(arr));
	if(n == 1)
		ret = arr;
		return
	else
		arr1 = mergesort(arr(1:ceil(n/2)));
		arr2 = mergesort(arr(ceil(n/2)+1:end));
	end
	
	p = size(arr1);
	q = size(arr2);
	
	ret = zeros(1,n);
	piv1 = 1;
	piv2 = 1;
	
	for i=1:n
		
		if(arr1(piv1) > arr2(piv2))
			ret(i) = arr2(piv2);
			piv2 = piv2+1;
		else
			ret(i) = arr1(piv1);
			piv1 = piv1+1;
		end
		
		if(piv1 > p)
			ret(i+1:end) = arr2(piv2:end);
			break
		elseif(piv2 > q)
			ret(i+1:end) = arr1(piv1:end);
			break
		end
	end
	
end