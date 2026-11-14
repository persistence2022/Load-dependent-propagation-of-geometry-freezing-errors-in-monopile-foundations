function h=hash_bytes(b)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(uint8(b));raw=typecast(md.digest(),'uint8');h=lower(reshape(dec2hex(raw,2).',1,[]));
end
