function [cell,sgno,atomparam,space_group_IT_number]=input_ni_fun()

cell = [3.5238 3.5238 3.5238 90 90 90];
sgno = 225;
atomparam(1).name = 'Ni';
atomparam(1).atomno = 28;
atomparam(1).pos = [0.00000 0.0000 0.0000];
atomparam(1).adp =0.01;
atomparam(1).occ = 1;
atomparam(1).symmulti = 48;
space_group_IT_number = 225;