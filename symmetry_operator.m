function sym = symmetry_operator(symmetry_name,use_miller_bravais)
% symmetry_name: cubic, hexagonal, orthorhombic, tetragonal, triclinic
% adapt from pymicro.crystal.lattice.py
% Jan 5, 2021
% def symmetry_operators(self, use_miller_bravais=False)
% Define the equivalent crystal symmetries.
% Those come from Randle & Engler, 2000. For instance in the cubic
% crystal struture, for instance there are 24 equivalent cube orientations.
% :returns array: A numpy array of shape (n, 3, 3) where n is the \
% number of symmetries of the given crystal structure.
if nargin<2
    use_miller_bravais=0;
end
switch symmetry_name
case 'cubic'
    sym=zeros(3,3,24);
    sym(:,:,1) = [1 0 0;0 1 0;0 0 1];
    sym(:,:,2) = [0 0 -1;0 -1 0;-1 0 0];
    sym(:,:,3) = [0 0 -1;0 1 0;1 0 0];
    sym(:,:,4) = [-1 0 0;0 1 0;0 0 -1];
    sym(:,:,5) = [0 0 1;0 1 0;-1 0 0];
    sym(:,:,6) = [1 0 0;0 0 -1;0 1 0];
    sym(:,:,7) = [1 0 0;0 -1 0;0 0 -1];
    sym(:,:,8) = [1 0 0;0 0 1;0 -1 0];
    sym(:,:,9) = [0 -1 0;1 0 0;0 0 1];
    sym(:,:,10) = [-1 0 0;0 -1 0;0 0 1];
    sym(:,:,11) = [0 1 0;-1 0 0;0 0 1];
    sym(:,:,12) = [0 0 1;1 0 0;0 1 0];
    sym(:,:,13) = [0 1 0;0 0 1;1 0 0];
    sym(:,:,14) = [0 0 -1;-1 0 0;0 1 0];
    sym(:,:,15) = [0 -1 0;0 0 1;-1 0 0];
    sym(:,:,16) = [0 1 0;0 0 -1;-1 0 0];
    sym(:,:,17) = [0 0 -1;1 0 0;0 -1 0];
    sym(:,:,18) = [0 0 1;-1 0 0;0 -1 0];
    sym(:,:,19) = [0 -1 0;0 0 -1;1 0 0];
    sym(:,:,20) = [0 1 0;1 0 0;0 0 -1];
    sym(:,:,21) = [-1 0 0;0 0 1;0 1 0];
    sym(:,:,22) = [0 0 1;0 -1 0;1 0 0];
    sym(:,:,23) = [0 -1 0;-1 0 0;0 0 -1];
    sym(:,:,24) = [-1 0 0;0 0 -1;0 -1 0];
case 'hexagonal'
    if use_miller_bravais==1 % using the Miller-Bravais representation here
        sym=zeros(4,4,12);
        sym(:,:,1) = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
        sym(:,:,2) = [0 0 1 0;1 0 0 0;0 0 1 0;0 0 0 1];
        sym(:,:,3) = [0 1 0 0;0 0 1 0;1 0 0 0;0 0 0 1];
        sym(:,:,4) = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 -1];
        sym(:,:,5) = [0 0 1 0;1 0 0 0;0 1 0 0;0 0 0 -1];
        sym(:,:,6) = [0 1 0 0;0 0 1 0;1 0 0 0;0 0 0 -1];
        sym(:,:,7) = [-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1];
        sym(:,:,8) = [0 0 -1 0;-1 0 0 0;0 -1 0 0;0 0 0 1];
        sym(:,:,9) = [0 -1 0 0;0 0 -1 0;-1 0 0 0;0 0 0 1];
        sym(:,:,10) = [-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1];
        sym(:,:,11) = [0 0 -1 0;-1 0 0 0;0 -1 0 0;0 0 0 -1];
        sym(:,:,12) = [0 -1 0 0;0 0 -1 0;-1 0 0 0;0 0 0 -1];
    else
        sym=zeros(3,3,12);
        s60 = sin(60 * pi / 180);
        sym(:,:,1) = [1 0 0;0 1 0;0 0 1];
        sym(:,:,2) = [0.5 s60 0;-s60 0.5 0;0 0 1];
        sym(:,:,3) = [-0.5 s60 0;-s60 -0.5 0;0 0 1];
        sym(:,:,4) = [-1 0 0;0 -1 0;0 0 1];
        sym(:,:,5) = [-0.5 -s60 0;s60 -0.5 0;0 0 1];
        sym(:,:,6) = [0.5 -s60 0;s60 0.5 0;0 0 1];
        sym(:,:,7) = [1 0 0;0 -1 0;0 0 -1];
        sym(:,:,8) = [0.5 s60 0;s60 -0.5 0;0 0 -1];
        sym(:,:,9) = [-0.5 s60 0;s60 0.5 0;0 0 -1];
        sym(:,:,10) = [-1 0 0;0 1 0;0 0 -1];
        sym(:,:,11) = [-0.5 -s60 0;-s60 0.5 0;0 0 -1];
        sym(:,:,12) = [0.5 -s60 0;-s60 -0.5 0;0 0 -1];
    end
case 'orthorhombic'
    sym=zeros(3,3,4);
    sym(:,:,1) = [1 0 0;0 1 0;0 0 1];
    sym(:,:,2) = [1 0 0;0 -1 0;0 0 -1];
    sym(:,:,3) = [-1 0 -1;0 1 0;0 0 -1];
    sym(:,:,4) = [-1 0 0;0 -1 0;0 0 1];
case 'tetragonal'
    sym=zeros(3,3,8);
    sym(:,:,1) = [1 0 0;0 1 0;0 0 1];
    sym(:,:,2) = [0 -1 0;1 0 0;0 0 1];
    sym(:,:,3) = [-1 0 0;0 -1 0;0 0 1];
    sym(:,:,4) = [0 1 0;-1 0 0;0 0 1]; 
    sym(:,:,5) = [1 0 0;0 -1 0;0 0 -1];
    sym(:,:,6) = [-1 0 0;0 1 0;0 0 -1];
    sym(:,:,7) = [0 1 0;1 0 0;0 0 -1];
    sym(:,:,8) = [0 -1 0;-1 0 0;0 0 -1];
case 'triclinic'
    sym=zeros(3,3);
    sym = [1 0 0;0 1 0;0 0 1];
otherwise
    sym=zeros(3,3);
    sym = [1 0 0;0 1 0;0 0 1];
    warning('symmetry not supported');
end

