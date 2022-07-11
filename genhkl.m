% Generate reflections up to maximum sin(theta)/lambda (sintlmax)
% The program follows the method described in: 
% Le Page and Gabe (1979) J. Appl. Cryst., 12, 464-466

function HALL = genhkl(cell,sysconditions,sintlmax)

SEG(:,:,1) = [ 0 0  0;  1 0 0; 0 1 0; 0 0  1];
SEG(:,:,2) = [-1 0  1; -1 0 0; 0 1 0; 0 0  1];
SEG(:,:,3) = [-1 1  0; -1 0 0; 0 1 0; 0 0 -1];
SEG(:,:,4) = [ 0 1 -1;  1 0 0; 0 1 0; 0 0 -1];

nref = 0;
H =[];     % Data of half sphere
HALL =[];  % Friedel mates included
STLALL = [];
sintlH = 0.0;

for i=1:4
 segn = i;
 % initialize the identifiers
 htest = 0;
 ktest = 0;
 ltest = 0;
 HLAST =SEG(1,:,segn);
 HSAVE= HLAST;
 sintlH=sintl(cell,HSAVE);



 while ltest == 0
   while ktest == 0
      while htest == 0
        nref = nref+1;
         if nref ~= 1
            if sysabs(HLAST,sysconditions) == 0
               H=[H; HLAST];
               HALL =[HALL; HLAST; -HLAST];
               STLALL = [STLALL sintlH sintlH];
            else
               nref=nref-1;
            end
         end
         HNEW = HLAST+SEG(2,:,segn);
         sintlH=sintl(cell,HNEW);
         if sintlH <= sintlmax
        	 HLAST = HNEW;
         else 
	         htest = 1;
         end 
      end
      
      HLAST(1) = HSAVE(1);
      HLAST = HLAST + SEG(3,:,segn);
      HNEW = HLAST;
      sintlH=sintl(cell,HNEW);
      if sintlH > sintlmax
         ktest = 1;
      end
      htest = 0;
   end

   HLAST(2) = HSAVE(2);
   HLAST = HLAST + SEG(4,:,segn);
   HNEW = HLAST;
   sintlH=sintl(cell,HNEW);
   if sintlH > sintlmax
      ltest = 1;
   end
   ktest = 0;
 end
end

nref = (nref -1)*2;

HALL = [HALL STLALL'];


