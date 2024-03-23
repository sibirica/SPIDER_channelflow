function [lib_names, krnl_info, tex_names] = defineLibrary(term_names, mode)
%{
    OBJECTIVE:
    Get useful variables for library construction & visualization from keys in "term_names"

    INPUTS:
    term_names - string vector containing keys for different library terms

    OUTPUTS:
    lib_names - string keys used for library construction & making struct containing coefficient results
    krnl_str - string input for MATLAB's native eval() function, for construting the integral kernel
    tex_names - labels for library terms in Latex format for plotting purposes
%}

fullLibrary = [];
fullKernels = {};

% Different libraries by case
if mode == "BC"
    %%%% Linear %%%%
    u = ["u", "${\bf u}$"];
    int_struct = {};
    int_struct.n_pieces = ["wy.*V"];
    int_struct.n_modes = [0, 0, 0, 0];
    int_struct.t_pieces = ["wx.*U"; "wz.*W"];
    int_struct.t_modes = [[0, 0, 0, 0]; ...
                          [0, 0, 0, 0]];
    int_struct.size_string = "mean_u";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; u]; 

    %%%% Inhomogeneous (normal) %%%%
    n = ["n", "$n$"];
    int_struct = {};
    int_struct.n_pieces = ["wy"];
    int_struct.n_modes = [0, 0, 0, 0];
    int_struct.t_pieces = ["0"];
    int_struct.t_modes = [0, 0, 0, 0];
    int_struct.t_skip = 1;
    int_struct.size_string = "1";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; n];

    %%%% Normal projection of squared velocity %%%%
    unu = ["unu", "$(u \cdot n)u$"];
    int_struct = {};
    int_struct.n_pieces = ["wy.*V.^2"];
    int_struct.n_modes = [0, 0, 0, 0];
    int_struct.t_pieces = ["wx.*U.*V"; "wz.*V.*W"];
    int_struct.t_modes = [[0, 0, 0, 0]; ...
                          [0, 0, 0, 0]];
    int_struct.size_string = "mean_u^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; unu];
    
    %%%% Normal with u %%%%
    unn = ["unn", "$(u \cdot n)n$"];
    int_struct = {};
    int_struct.n_pieces = ["wy.*V"];
    int_struct.n_modes = [0, 0, 0, 0];
    int_struct.t_pieces = [0];
    int_struct.t_skip = 1;
    int_struct.size_string = "mean_u";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; unn];

    %%%% Linear with quadratic velocity %%%%
    u3 = ["u3", "$(u^2 + v^2 + w^2){\bf u}$"];
    int_struct = {};
    int_struct.n_pieces = ["wy.*V.*(U.^2 + V.^2 + W.^2)"];
    int_struct.n_modes = [0, 0, 0, 0];
    int_struct.t_pieces = ["wx.*U.*(U.^2 + V.^2 + W.^2)"; ...
                           "wz.*W.*(U.^2 + V.^2 + W.^2)"];
    int_struct.t_modes = [[0, 0, 0, 0]; ...
                          [0, 0, 0, 0]];
    int_struct.size_string = "mean_u^3";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; u3];

    %%%% Normal derivative of velocity %%%%
    dnu = ["dnu", "$d_n u$"];
    int_struct = {};
    int_struct.n_pieces = ["wy.*dV"];
    int_struct.n_modes = [0, 0, 0, 0];
    int_struct.t_pieces = ["wx.*dU"; ...
                           "wz.*dW"];
    int_struct.t_modes = [[0, 0, 0, 0]; ...
                          [0, 0, 0, 0]];
    int_struct.size_string = "std_u/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; dnu];

    %%%% Gradient of normal velocity %%%%
    du_n = ["du_n", "$\nabla(u_n)$"];
    int_struct = {};
    int_struct.n_pieces = ["wy.*dV"];
    int_struct.n_modes = [0, 0, 0, 0];
    int_struct.n_skip = 1;
    int_struct.t_pieces = ["-wx.*V*S_x"; ...
                           "-wz.*V*S_z"];
    int_struct.t_modes = [[1, 0, 0, 0]; ...
                          [0, 0, 1, 0]];
    int_struct.size_string = "std_u/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; du_n];

    %%%% Advection %%%%
    uu_x = ["uu_x", "$({\bf u}\cdot\nabla){\bf u}$"];
    int_struct = {};
    int_struct.n_pieces = ["-U.*V.*wy*S_x"; ...
                           "-V.*W.*wy*S_z"; ...
                           "2*V.*wy.*dV"];
    int_struct.n_modes = [[1, 0, 0, 0]; ...
                          [0, 0, 1, 0]; ...
                          [0, 0, 0, 0]];
    int_struct.t_pieces = ["-U.*U.*wx*S_x"; ...
                           "-U.*W.*wx*S_z"; ...
                           "U.*dV.*wx"; ...
                           "V.*dU.*wx"; ...
                           "-U.*W.*wz*S_x"; ...
                           "-W.*W.*wz*S_z"; ...
                           "V.*dW.*wz"; ...
                           "W.*dV.*wz"];
    int_struct.t_modes = [[1, 0, 0, 0]; ...
                          [0, 0, 1, 0]; ...
                          [0, 0, 0, 0]; ...
                          [0, 0, 0, 0]; ...
                          [1, 0, 0, 0]; ...
                          [0, 0, 1, 0]; ...
                          [0, 0, 0, 0]; ...
                          [0, 0, 0, 0]];
    int_struct.size_string = "mean_u*std_u/dx";
    fullKernels(end + 1) = {int_struct};                
    fullLibrary = [fullLibrary; uu_x];
    
    %%%% Linear in n with velocity square prefactor %%%%
    u2n = ["u2n", "{\bf u}^2{\bf n}"];
    int_struct = {};
    int_struct.n_pieces = ["(U.*U+V.*V+W.*W).*wy"];
    int_struct.n_modes = [[0, 0, 0, 0]];
    int_struct.t_pieces = ["0"];
    int_struct.t_modes = [[0, 0, 0, 0]];
    int_struct.t_skip = 1;
    int_struct.size_string = "mean_u^2";
    fullKernels(end + 1) = {int_struct};                
    fullLibrary = [fullLibrary; u2n];
    
    
    %%%% Gradient of velocity squared %%%%
    du2 = ["du2", "\nabla({\bf u})"];
    int_struct = {};
    int_struct.n_pieces = ["(U.*dU+V.*dV+W.*dW).*wy"];
    int_struct.n_modes = [[0, 0, 0, 0]];
    int_struct.t_pieces = ["-(U.*U+V.*V+W.*W).*wx";
                           "-(U.*U+V.*V+W.*W).*wz"];
    int_struct.t_modes = [[1, 0, 0, 0]
                          [0, 0, 1, 0]];
    int_struct.size_string = "mean_u*std_u/dx";
    fullKernels(end + 1) = {int_struct};                
    fullLibrary = [fullLibrary; du2];
    
    %%%% Laplacian of u %%%%
    u_xx = ["u_xx", "$\nabla^2{\bf u}$"];
    int_struct = {};
    int_struct.n_pieces = ["wy.*V*S_x^2"; "wy.*dyyV"; "wy.*V*S_z^2"];
    int_struct.n_modes = [[2, 0, 0, 0]
                          [0, 0, 0, 0]
                          [0, 0, 2, 0]];
    int_struct.t_pieces = ["wx.*U*S_x^2"; "wx.*dyyU"; "wx.*U*S_z^2";
                           "wz.*W*S_x^2"; "wz.*dyyW"; "wz.*W*S_z^2"];
    int_struct.t_modes = [[2, 0, 0, 0]
                          [0, 0, 0, 0]
                          [0, 0, 2, 0]
                          [2, 0, 0, 0]
                          [0, 0, 0, 0]
                          [0, 0, 2, 0]];
    int_struct.size_string = "mean_u/dx^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; u_xx];
    
elseif mode == "NS"   
    %%%% Time derivative %%%%
    u_t = ["u_t", "$\partial_t{\bf u}$"];       % lib_name, tex_name       
    fullLibrary = [fullLibrary; u_t];
    int_struct.pieces = ["-U.*w.*S_t"];
    int_struct.modes = [[0, 0, 0, 1]];
    int_struct.size_string = "std_u/dt";
    fullKernels(end + 1) = {int_struct};

    %%%% Second time derivative %%%%
    u_tt = ["u_tt", "$\partial^2_t{\bf u}$"];       % lib_name, tex_name       
    fullLibrary = [fullLibrary; u_tt];
    int_struct.pieces = ["U.*w.*S_t^2"];
    int_struct.modes = [[0, 0, 0, 2]];
    int_struct.size_string = "std_u/dt^2";
    fullKernels(end + 1) = {int_struct};

    %%%% Advection %%%%
    uu_x = ["uu_x", "$({\bf u}\cdot\nabla){\bf u}$"];
    int_struct.pieces = ["-U.*U.*w*S_x"; ...
                         "-U.*V.*w*S_y"; ...
                         "-U.*W.*w*S_z"];
    int_struct.modes = [[1, 0, 0, 0]; ...
                        [0, 1, 0, 0]; ...
                        [0, 0, 1, 0]];
    int_struct.size_string = "mean_u*std_u/dx";
    fullKernels(end + 1) = {int_struct};                
    fullLibrary = [fullLibrary; uu_x];

    %%%% Laplacian %%%%
    u_xx = ["u_xx", "$\nabla^2{\bf u}$"];
    int_struct.pieces = ["U.*w*S_x^2"; ...
                         "U.*w*S_y^2"; ...
                         "U.*w*S_z^2"];
    int_struct.modes = [[2, 0, 0, 0]; ...
                        [0, 2, 0, 0]; ...
                        [0, 0, 2, 0]];
    int_struct.size_string = "std_u/dx^2";
    fullKernels(end + 1) = {int_struct};            
    fullLibrary = [fullLibrary; u_xx];

    %%%% Linear %%%%
    u = ["u", "${\bf u}$"];
    int_struct.pieces = ["U.*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_u";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; u];

    %%%% Linear with quadratic velocity %%%%
    u3 = ["u3", "$(u^2 + v^2 + w^2){\bf u}$"];
    int_struct.pieces = ["(U.^2 + V.^2 + W.^2).*U.*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_u^3";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; u3];
    
    %%%% Linear in u with pressure square pre-factor %%%%
    p2u = ["p2u", "$p^2{\bf u}$"];
    int_struct.pieces = ["w.*U.*p.*p"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_u*mean_p^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; p2u];

    %%%% Gradient of pressure %%%%
    dp = ["dp", "$\nabla p$"];
    int_struct.pieces = ["-p.*w.*S_x"];
    int_struct.modes = [[1, 0, 0, 0]];
    int_struct.size_string = "std_p/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; dp];

    %%%% Gradient of pressure squared %%%%
    pdp = ["pdp", "$p \nabla p$"];
    int_struct.pieces = ["-p.*p.*w.*S_x/2"];
    int_struct.modes = [[1, 0, 0, 0]];
    int_struct.size_string = "mean_p*std_p/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; pdp];
    
    %%%% Velocity squared gradient %%%%
    dvv = ["dvv", "$\nabla ({\bf u}\cdot {\bf u})$"];
    int_struct.pieces = ["-(U.*U+V.*V+W.*W).*p.*w.*S_x/2"];
    int_struct.modes = [[1, 0, 0, 0]];
    int_struct.size_string = "mean_u*std_u/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; dvv];

    %%%% p times u %%%%
    pu = ["pu", "$pu$"];
    int_struct.pieces = ["p.*U.*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_p*mean_u";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; pu];
    
    %%%% p times Laplacian of u %%%%
    pu_xx = ["pu_xx", "$p\nabla^2 u$"];
    int_struct.pieces = ["p.*dxxU.*w"; ...
                         "p.*dyyU.*w"; ...
                         "p.*dzzU.*w"];
    int_struct.modes = [[0, 0, 0, 0]; ...
                        [0, 0, 0, 0]; ...
                        [0, 0, 0, 0]];
    int_struct.size_string = "mean_p*std_u/dx^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; pu_xx];
    
    %%%% u times Laplacian of p %%%%
    up_xx = ["up_xx", "$u\nabla^2 p$$"];
    int_struct.pieces = ["U.*dxxp.*w"; ...
                         "U.*dyyp.*w"; ...
                         "U.*dzzp.*w"];
    int_struct.modes = [[0, 0, 0, 0]; ...
                        [0, 0, 0, 0]; ...
                        [0, 0, 0, 0]];
    int_struct.size_string = "mean_u*std_p/dx^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; up_xx];
    
    %%%% u times time derivative of p %%%%
    up_t = ["up_t", "$u\partial_t p$$"];
    int_struct.pieces = ["U.*dtp.*w"; ...
                         "U.*dtp.*w"; ...
                         "U.*dtp.*w"];
    int_struct.modes = [[0, 0, 0, 0]; ...
                        [0, 0, 0, 0]; ...
                        [0, 0, 0, 0]];
    int_struct.size_string = "mean_u*std_p/dt";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; up_t];
    
    %%%% p times time derivative of u %%%%
    pu_t = ["pu_t", "$p\partial_t u$"];
    int_struct.pieces = ["p.*dtU.*w"; ...
                         "p.*dtU.*w"; ...
                         "p.*dtU.*w"];
    int_struct.modes = [[0, 0, 0, 0]; ...
                        [0, 0, 0, 0]; ...
                        [0, 0, 0, 0]];
    int_struct.size_string = "mean_p*std_u/dt";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; pu_t];
elseif mode == "div"
    %%%% Divergence of u %%%%
    div_u = ["div_u", "$\nabla \cdot u$"];      
    int_struct.pieces = ["-U.*w*S_x"; ...
                         "-V.*w*S_y"; ...
                         "-W.*w*S_z"];
    int_struct.modes = [[1, 0, 0, 0]; ...
                        [0, 1, 0, 0]; ...
                        [0, 0, 1, 0]];
    int_struct.size_string = "std_u/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; div_u];

    %%%% Linear %%%%
    p = ["p", "${\bf p}$"];
    int_struct.pieces = ["p.*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_p";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; p];

    %%%% Time derivative %%%%
    p_t = ["p_t", "$\partial_tp$"];    
    int_struct.pieces = ["-p.*w.*S_t"];
    int_struct.modes = [[0, 0, 0, 1]];
    int_struct.size_string = "std_p/dt";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; p_t];
    
    %%%% Second time derivative %%%%
    p_tt = ["p_tt", "$\partial_t^2p$"];    
    int_struct.pieces = ["p.*w.*S_t^2"];
    int_struct.modes = [[0, 0, 0, 2]];
    int_struct.size_string = "std_p/dt^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; p_tt];

    %%%% Quadratic in p %%%%
    p2 = ["p2", "$p^2$"];
    int_struct.pieces = ["p.*p.*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_p^2";
    fullKernels(end + 1) = {int_struct};            
    fullLibrary = [fullLibrary; p2];

    %%%% Quadratic in u %%%%
    u2 = ["u2", "$u^2 + v^2 + w^2$"];
    int_struct.pieces = ["(U.^2 + V.^2 + W.^2).*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_u^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; u2];

    %%%% Advection-like term %%%%
    udp = ["udp", "$(u \cdot \nabla)p$"];
    int_struct.pieces = ["U.*dxp.*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_u*std_p/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; udp];
    
    %%%% Laplacian of p %%%%
    p_xx = ["p_xx", "$\nabla^2 p$"];
    int_struct.pieces = ["p.*w*S_x^2"; ...
                         "p.*w*S_y^2"; ...
                         "p.*w*S_z^2"];
    int_struct.modes = [[2, 0, 0, 0]; ...
                        [0, 2, 0, 0]; ...
                        [0, 0, 2, 0]];
    int_struct.size_string = "std_p/dx^2";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; p_xx];
    
    %%%% p times time derivative of p %%%%
    pp_t = ["pp_t", "$p \partial_t p$"];
    int_struct.pieces = ["p.*dtp.*w"];
    int_struct.size_string = "std_p*mean_p/dt";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; pp_t];
    
    %%%% p times divergence of u %%%%
    pdiv_u = ["pdiv_u", "$p(\nabla \cdot {\bf u})$"];
    int_struct.pieces = ["p.*dxU.*w"; ...
                         "p.*dyV.*w"; ...
                         "p.*dzW.*w"];
    int_struct.modes = [[0, 0, 0, 0]; ...
                        [0, 0, 0, 0]; ...
                        [0, 0, 0, 0]];
    int_struct.size_string = "std_u*mean_p/dx";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; pdiv_u];
    
    %%%% u^2 times p %%%%
    u2p = ["u2p", "$u^2p$"];
    int_struct.pieces = ["p.*(U.^2+V.^2+W.^2).*w"];
    int_struct.modes = [[0, 0, 0, 0]];
    int_struct.size_string = "mean_u^2*mean_p";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; u2p];
    
    %%%% Time derivative of velocity squared %%%%
    uu_t = ["uu_t", "$u \cdot \partial_t u$"];
    int_struct.pieces = ["(U.^2+V.^2+W.^2).*w*S_t/2"];
    int_struct.modes = [[0, 0, 0, 1]];
    int_struct.size_string = "mean_u*std_u/dt";
    fullKernels(end + 1) = {int_struct};
    fullLibrary = [fullLibrary; uu_t];
end

krnl_info = containers.Map;
tex_names = [];
lib_names = [];

while length(term_names) >= 1
    ind = fullLibrary(:,1) == term_names(1);
    if sum(ind) == 0
        error("Library key > " + term_names(1) + " < not found in Library.")
    end
    krnl_info(term_names(1)) = fullKernels{ind};
    lib_names = [lib_names, term_names(1)];
    tex_names = [tex_names, fullLibrary(ind, 2)];
    term_names = term_names(2:end);
end

end