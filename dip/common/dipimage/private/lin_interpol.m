function out = lin_interpol(in, vals)

size_in = size(in);
in = in(:);
out = zeros(size_in(1), size_in(2));

for i=1:size( in, 1)
	io = find( vals(:, 1) <= in(i));
	io2 = find( vals(:, 1) >= in(i));
   if( ~isempty(io))
	   io_low = io(end);
   else
      io_low = find( vals(:, 1) == min( vals(:, 1)));
   end
   if( ~isempty(io2))
	   io_high = io2(1);
   else
      io_high = find( vals(:, 1) == max( vals(:, 1)));
   end  

	if( io_low ~= io_high)
      val_low = vals(io_low, 2);
      val_high = vals(io_high, 2);
		dist_low = in(i)-vals(io_low, 1);
		dist_high = vals(io_high, 1)-in(i);
		out(i) = (dist_low.* val_high + dist_high .* val_low)./( dist_high + dist_low);
   else
		out(i) = vals(io_high, 2);
	end
end
