function ar=di_cast_to_mwsize(ar)
  mwsizeof=di_mwtypes_sizeof();
  switch mwsizeof(1)
    case 4
      ar=int32(ar);
    case 8
      ar=int64(ar);
    otherwise
      error('cannot cast to mwSize');
  end
end

