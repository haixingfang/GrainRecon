% sysabs checks whether a reflection is systematic absent
% 
% sysabs = sysabs(hkl,syscond)
%
% INPUT: hkl = [h k l] 
%         syscond: [1x23] with condition for systematic absences in this
%         space group, X in syscond should given as shown below
% OUTPUT: sysbs: if 1 the reflection is systematic absent 
%                if 0 its not
%
%syscond:
%class        systematic abs               sysconditions(i)
%HKL          H+K=XN                            1
%             H+L=XN                            2
%             K+L=XN                            3
%             H+K,H+L,K+L = XN                  4
%             H+K+L=XN                          5
%             -H+K+L=XN                         6 
%HHL          H=XN                              7
%             L=XN                              8
%             H+L=XN                            9
%             2H+L=XN                          10
%0KL          K=XN                             11
%             L=XN                             12
%             K+L=XN                           13
%H0L          H=XN                             14
%             L=XN                             15
%             H+L=XN                           16
%HK0          H=XN                             17
%             K=XN                             18
%             H+K=XN                           19
%HH0          H=XN                             20
%H00          H=XN                             21
%0K0          K=XN                             22
%00L          L=XN                             23

function sysabs = sysabs(HKL,syscond)

h = HKL(1);
k = HKL(2);
l = HKL(3);
sysabs=0;

%HKL class
if syscond(1) ~= 0
   x = syscond(1);
   if mod(abs(h+k),x) ~=0
     sysabs=1;
   end 
end
if syscond(2) ~= 0
   x = syscond(2);
   if mod(abs(h+l),x) ~=0
     sysabs=1;
   end 
end
if syscond(3) ~= 0
   x = syscond(3);
   if mod(abs(k+l),x) ~=0
     sysabs=1;
   end 
end
if syscond(4) ~= 0
   sysabs=1;
   x = syscond(4);
   if mod(abs(h+k),x) ==0
     if mod(abs(h+l),x) ==0
       if  mod(abs(k+l),x) ==0
          sysabs=0;
       end 
     end
   end
end
if syscond(5) ~= 0
   x = syscond(5);
   if mod(abs(h+k+l),x) ~=0
       sysabs=1;
   end
end
if syscond(6) ~= 0
   x = syscond(6);
   if mod(abs(-h+k+l),x) ~=0
       sysabs=1;
   end
end

%HHL class
if (h-k)==0
  if syscond(7) ~= 0
     x = syscond(7);
     if mod(abs(h),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(8) ~= 0
     x = syscond(8);
     if mod(abs(l),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(9) ~= 0
     x = syscond(9);
     if mod(abs(h+l),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(10) ~= 0
     x = syscond(10);
     if mod(abs(h+h+l),x) ~= 0
        sysabs = 1;
     end
  end
end

%0KL class
if h==0
  if syscond(11) ~= 0
     x = syscond(11);
     if mod(abs(k),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(12) ~= 0
     x = syscond(12);
     if mod(abs(l),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(13) ~= 0
     x = syscond(13);
     if mod(abs(k+l),x) ~= 0
        sysabs = 1;
     end
  end
end

%H0L class
if k==0
  if syscond(14) ~= 0
     x = syscond(14);
     if mod(abs(h),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(15) ~= 0
     x = syscond(15);
     if mod(abs(l),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(16) ~= 0
     x = syscond(16);
     if mod(abs(h+l),x) ~= 0
        sysabs = 1;
     end
  end
end


%HK0 class
if l==0
  if syscond(17) ~= 0
     x = syscond(17);
     if mod(abs(h),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(18) ~= 0
     x = syscond(18);
     if mod(abs(k),x) ~= 0
        sysabs = 1;
     end
  end
  if syscond(19) ~= 0
     x = syscond(19);
     if mod(abs(h+k),x) ~= 0
        sysabs = 1;
     end
  end
end

%HH0 class
if l==0
if h-k==0
  if syscond(20) ~= 0
     x = syscond(20);
     if mod(abs(h),x) ~= 0
        sysabs = 1;
     end
  end
end
end

%H00 class
if abs(k)+abs(l)==0
  if syscond(21) ~= 0
     x = syscond(21);
     if mod(abs(h),x) ~= 0
        sysabs = 1;
     end
  end
end

%0K0 class
if abs(h)+abs(l)==0
  if syscond(22) ~= 0
     x = syscond(22);
     if mod(abs(k),x) ~= 0
        sysabs = 1;
     end
  end
end

%00L class
if abs(h)+abs(k)==0
  if syscond(23) ~= 0
     x = syscond(23);
     if mod(abs(l),x) ~= 0
        sysabs = 1;
     end
  end
end

