function uv_access_writable(path, mode = 2)
    local ret
    req = Libc.malloc(Base._sizeof_uv_fs)
	try
		if path == ""
			path = "."
		end
        ret = ccall(:uv_fs_access, Int32, (Ptr{Nothing}, Ptr{Nothing}, Cstring, Int64, Ptr{Nothing}), Base.eventloop(), req, path, mode, C_NULL)
        ccall(:uv_fs_req_cleanup, Nothing, (Ptr{Nothing},), req)
    finally
        Libc.free(req)
    end
    return ret == 0
end
