import sys
import numpy as np

# --- Generates qdldl files ---
def gen_qdldl_solve(fname, L_p, L_i, P, path=''):
    n = len(L_p) - 1

    ident = '    '

    # --- Generates qdldl file ---
    qdldl_c = '#include \"{:}\"\n'.format(fname + '.h')

    # Generates lsolve function
    lsolve_name = 'QDLDL_Lsolve_st'
    lsolve_sig = 'static inline void ' + lsolve_name + '(const QDLDL_float Lx[{:}], QDLDL_float x[{:}]);'.format(len(L_i), n)
    lsolve = lsolve_sig[:-1] + '{\n'
    lsolve += '\n' + ident + 'QDLDL_float val;\n'

    temp = ''
    for i in range(n - 1):
        temp += '\n' + ident[:-1] + '/*\n' +\
                ident + '* i = {:}\n'.format(i) +\
                ident + '* L_p[{:}] = {:}, L_p[{:}] = {:}\n'.format(i, L_p[i], i + 1, L_p[i + 1]) +\
                ident + '*/\n'
        temp += ident + 'val = x[{:}];\n'.format(i)
        for j in range(L_p[i], L_p[i + 1]):
            temp += ident + 'x[{:}] -= Lx[{:}] * val;\n'.format(L_i[j], j)

    lsolve += temp + '}\n'

    # Generates ltsolve function
    ltsolve_name = 'QDLDL_Ltsolve_st'
    ltsolve_sig = 'static inline void ' + ltsolve_name + '(const QDLDL_float Lx[{:}], QDLDL_float x[{:}]);'.format(len(L_i), n)
    ltsolve = ltsolve_sig[:-1] + '{\n'
    ltsolve += '\n' + ident + 'QDLDL_float val;\n'

    temp = ''
    for i in reversed(range(n - 1)):
        temp += '\n' + ident[:-1] + '/*\n' +\
                ident + '* i = {:}\n'.format(i) +\
                ident + '* L_p[{:}] = {:}, L_p[{:}] = {:}\n'.format(i, L_p[i], i + 1, L_p[i + 1]) +\
                ident + '*/\n'
        temp += ident + 'val = x[{:}];\n'.format(i)
        for j in range(L_p[i], L_p[i + 1]):
            temp += ident + 'val -= Lx[{:}] * x[{:}];\n'.format(j, L_i[j])
        temp += ident + 'x[{:}] = val;\n'.format(i)

    ltsolve += temp + '}\n'

    # Generates QDLDL solve function
    solve_name = 'QDLDL_solve_st'
    solve_sig = 'void ' + solve_name + '(const QDLDL_float Lx[{:}], const QDLDL_float Dinv[{:}], QDLDL_float x[{:}]);'.format(len(L_i), n, n)
    solve = solve_sig[:-1] + '{\n'
    #solve += '\n' + ident + 'QDLDL_int i;\n'
    solve += '\n' + ident + lsolve_name + '(Lx, x);\n'
    #solve += ident + 'for (i = 0; i < {:}; i++) x[i] *= Dinv[i];'.format(n)
    solve += '\n'
    for i in range(n):
        solve += ident + 'x[{:}] *= Dinv[{:}];\n'.format(i, i)
    solve += '\n' + ident + ltsolve_name + '(Lx, x);\n'
    solve += '}\n'

    # Generates permute_function
    permute_x_name = 'permute_x_st'

    permute_x_sig = 'static inline void ' + permute_x_name + '(QDLDL_float x[{:}], const QDLDL_float b[{:}]);'.format(n, n)

    permute_x = permute_x_sig[:-1] + '{\n'
    permute_x += ident + '/* Permute x = P*b using P */\n'
    for i in range(n):
        permute_x += ident + 'x[{:}] = b[{:}];\n'.format(i, P[i])
    permute_x += '}\n'

    # Generates permutet_function
    permutet_x_name = 'permutet_x_st'

    permutet_x_sig = 'static inline void ' + permutet_x_name + '(QDLDL_float x[{:}], const QDLDL_float b[{:}]);'.format(n, n)

    permutet_x = permutet_x_sig[:-1] + '{\n'
    permutet_x += ident + '/* Permute x = P\'*b using P */\n'
    for i in range(n):
        permutet_x += ident + 'x[{:}] = b[{:}];\n'.format(P[i], i)
    permutet_x += '}\n'
    
    # Generates LDLSolve function
    ldl_solve_name = 'LDLSolve_st'
    ldl_solve_sig = 'void ' + ldl_solve_name + '(QDLDL_float *x, QDLDL_float *b, const csc *L, const c_float *Dinv, c_float *bp);'

    ldl_solve = ldl_solve_sig[:-1] + '{\n'
    ldl_solve += ident + '/* Solves P\'LDL\'P x = b for x */\n'
    ldl_solve += ident + permute_x_name + '(bp, b);\n'
    ldl_solve += ident + solve_name + '(L->x, Dinv, bp);\n'
    ldl_solve += ident + permutet_x_name + '(x, bp);\n'
    ldl_solve += '}\n'

    qdldl_c += '\n' + lsolve + '\n' + ltsolve + '\n' + solve
    qdldl_c += '\n' + permute_x + '\n' + permutet_x + '\n' + ldl_solve
    
    # Generates the header file
    qdldl_h = '#ifndef QDLDL_ST_H\n#define QDLDL_ST_H\n'
    qdldl_h += '\n#include \"qdldl_types.h\"\n'
    qdldl_h += '#include \"types.h\"\n'
    #qdldl_h += '\n' + lsolve_sig + '\n'
    #qdldl_h += ltsolve_sig + '\n'
    qdldl_h +=  '\n' + solve_sig + '\n' + ldl_solve_sig + '\n'
    qdldl_h += '\n#endif /* QDLDL_ST_H */\n'

    # Saves files
    if len(path) > 0:
        if path[-1] != '\\' or path[-1] != '/':
            path = path + '/'
    with open(path + fname + '.c', mode='w') as f:
        f.write(qdldl_c)
    with open(path + fname + '.h', mode='w') as f:
        f.write(qdldl_h)


def gen_qdldl_solve_v2(fname, L_p, L_i, P, D_inv, L_x, path=''):

    def np_array_to_c(arr, arr_name):

        if arr.ndim == 1:
            n = arr.shape[0]
            m = 1
        else:
            if (arr.shape[0] == 1) or (arr.shape[1] == 1):
                arr = arr.flatten()
                n = arr.shape[0]
                m = 1
            else:
                n, m = arr.shape

        arr_str = np.array2string(arr, separator=',')
        arr_str = arr_str.replace('[', '{')
        arr_str = arr_str.replace(']', '}')

        if m == 1:
            arr_str = '{:}[{:}] = {:};'.format(arr_name, n, arr_str)
        else:
            arr_str = '{:}[{:}][{:}] = {:};'.format(arr_name, n, m, arr_str)

        return arr_str

    #np.set_printoptions(precision=50, threshold=sys.maxsize)
    np.set_printoptions(floatmode='unique', threshold=sys.maxsize)
    
    n = len(L_p) - 1

    ident = '    '

    # --- Generates qdldl file ---
    qdldl_c = '#include \"{:}\"\n'.format(fname + '.h')

    # Generates lsolve function
    lsolve_name = 'QDLDL_Lsolve_st'
    lsolve_sig = 'static inline void ' + lsolve_name + '(const QDLDL_float Lx[{:}], QDLDL_float x[{:}]);'.format(len(L_i), n)
    lsolve = lsolve_sig[:-1] + '{\n'
    lsolve += '\n' + ident + 'QDLDL_float val;\n'

    temp = ''
    for i in range(n - 1):
        temp += '\n' + ident[:-1] + '/*\n' +\
                ident + '* i = {:}\n'.format(i) +\
                ident + '* L_p[{:}] = {:}, L_p[{:}] = {:}\n'.format(i, L_p[i], i + 1, L_p[i + 1]) +\
                ident + '*/\n'
        temp += ident + 'val = x[{:}];\n'.format(i)
        for j in range(L_p[i], L_p[i + 1]):
            temp += ident + 'x[{:}] -= Lx[{:}] * val;\n'.format(L_i[j], j)

    lsolve += temp + '}\n'

    # Generates ltsolve function
    ltsolve_name = 'QDLDL_Ltsolve_st'
    ltsolve_sig = 'static inline void ' + ltsolve_name + '(const QDLDL_float Lx[{:}], QDLDL_float x[{:}]);'.format(len(L_i), n)
    ltsolve = ltsolve_sig[:-1] + '{\n'
    ltsolve += '\n' + ident + 'QDLDL_float val;\n'

    temp = ''
    for i in reversed(range(n - 1)):
        temp += '\n' + ident[:-1] + '/*\n' +\
                ident + '* i = {:}\n'.format(i) +\
                ident + '* L_p[{:}] = {:}, L_p[{:}] = {:}\n'.format(i, L_p[i], i + 1, L_p[i + 1]) +\
                ident + '*/\n'
        temp += ident + 'val = x[{:}];\n'.format(i)
        for j in range(L_p[i], L_p[i + 1]):
            temp += ident + 'val -= Lx[{:}] * x[{:}];\n'.format(j, L_i[j])
        temp += ident + 'x[{:}] = val;\n'.format(i)

    ltsolve += temp + '}\n'

    # Generates QDLDL solve function
    solve_name = 'QDLDL_solve_st'
    solve_sig = 'void ' + solve_name + '(QDLDL_float x[{:}]);'.format(n)
    solve = solve_sig[:-1] + '{\n'
    #solve += '\n' + ident + 'QDLDL_int i;\n'
    solve += '\n' + ident + np_array_to_c(D_inv, 'static const QDLDL_float Dinv') + '\n'
    solve += '\n' + ident + np_array_to_c(L_x, 'static const QDLDL_float Lx') + '\n'
    solve += '\n' + ident + lsolve_name + '(Lx, x);\n'
    #solve += ident + 'for (i = 0; i < {:}; i++) x[i] *= Dinv[i];'.format(n)
    solve += '\n'
    for i in range(n):
        solve += ident + 'x[{:}] *= Dinv[{:}];\n'.format(i, i)
    solve += '\n' + ident + ltsolve_name + '(Lx, x);\n'
    solve += '}\n'

    # Generates permute_function
    permute_x_name = 'permute_x_st'

    permute_x_sig = 'static inline void ' + permute_x_name + '(QDLDL_float x[{:}], const QDLDL_float b[{:}]);'.format(n, n)

    permute_x = permute_x_sig[:-1] + '{\n'
    permute_x += ident + '/* Permute x = P*b using P */\n'
    for i in range(n):
        permute_x += ident + 'x[{:}] = b[{:}];\n'.format(i, P[i])
    permute_x += '}\n'

    # Generates permutet_function
    permutet_x_name = 'permutet_x_st'

    permutet_x_sig = 'static inline void ' + permutet_x_name + '(QDLDL_float x[{:}], const QDLDL_float b[{:}]);'.format(n, n)

    permutet_x = permutet_x_sig[:-1] + '{\n'
    permutet_x += ident + '/* Permute x = P\'*b using P */\n'
    for i in range(n):
        permutet_x += ident + 'x[{:}] = b[{:}];\n'.format(P[i], i)
    permutet_x += '}\n'
    
    # Generates LDLSolve function
    ldl_solve_name = 'LDLSolve_st'
    ldl_solve_sig = 'void ' + ldl_solve_name + '(QDLDL_float *x, QDLDL_float *b, c_float *bp);'

    ldl_solve = ldl_solve_sig[:-1] + '{\n'
    ldl_solve += ident + '/* Solves P\'LDL\'P x = b for x */\n'
    ldl_solve += ident + permute_x_name + '(bp, b);\n'
    ldl_solve += ident + solve_name + '(bp);\n'
    ldl_solve += ident + permutet_x_name + '(x, bp);\n'
    ldl_solve += '}\n'

    qdldl_c += '\n' + lsolve + '\n' + ltsolve + '\n' + solve
    qdldl_c += '\n' + permute_x + '\n' + permutet_x + '\n' + ldl_solve
    
    # Generates the header file
    qdldl_h = '#ifndef QDLDL_ST_H\n#define QDLDL_ST_H\n'
    qdldl_h += '\n#include \"qdldl_types.h\"\n'
    qdldl_h += '#include \"types.h\"\n'
    #qdldl_h += '\n' + lsolve_sig + '\n'
    #qdldl_h += ltsolve_sig + '\n'
    qdldl_h +=  '\n' + solve_sig + '\n' + ldl_solve_sig + '\n'
    qdldl_h += '\n#endif /* QDLDL_ST_H */\n'

    # Saves files
    if len(path) > 0:
        if path[-1] != '\\' or path[-1] != '/':
            path = path + '/'
    with open(path + fname + '.c', mode='w') as f:
        f.write(qdldl_c)
    with open(path + fname + '.h', mode='w') as f:
        f.write(qdldl_h)

    np.set_printoptions(floatmode='fixed', threshold=1000)
    #np.set_printoptions(precision=8, threshold=1000)

        
# --- Generatex auxil files ---
def gen_auxil(fname, n, m, alpha, path=''):

    ident = '    '
    
    auxil_c = '#include \"{:}\"\n'.format(fname + '.h')
    auxil_c += '#include \"stdint.h\"\n'

    # Generates update_xzy function
    update_xzy_name = 'update_xzy_st'
    x_args = '\n' + ident + 'float x[{:}], float x_prev[{:}], float delta_x[{:}],\n'.format(n, n, n)
    z_args = ident + 'float z[{:}], float z_prev[{:}],\n'.format(m, m)
    y_args = ident + 'float y[{:}], float delta_y[{:}],\n'.format(m, m)
    xz_tilde_args = ident + 'float xz_tilde[{:}],\n'.format(n + m)
    rho_args = ident + 'float rho_vec[{:}], float rho_inv_vec[{:}],\n'.format(m, m)
    bounds_args = ident + 'float l[{:}], float u[{:}]'.format(m, m)
    
    update_xzy_sig = 'void ' + update_xzy_name + '('
    update_xzy_sig += x_args + z_args + y_args + xz_tilde_args + rho_args + bounds_args
    update_xzy_sig += '\n' + ident + ');'

    update_xzy = '\n' + update_xzy_sig[:-1] + '{\n'
    
    update_xzy += '\n' + ident + 'int32_t i;\n'
    update_xzy += ident + 'float alpha = (float){:};\n'.format(alpha)

    update_x = '\n' + ident + '/* Updates x*/\n'
    update_x += ident + 'for(i = 0; i < {:}; i++){{\n'.format(n)
    update_x += 2 * ident + 'x[i] = alpha * xz_tilde[i] + ( ((float)1.0f) - alpha) * x_prev[i];\n'
    update_x += 2 * ident + 'delta_x[i] = x[i] - x_prev[i];\n'
    update_x += ident + '}\n'

    update_z = '\n' + ident + '/* Updates z*/\n'
    update_z += ident + 'for(i = 0; i < {:}; i++){{\n'.format(m)
    update_z += 2 * ident + 'z[i] = alpha * xz_tilde[i + {:}] + ( ((float)1.0f) - alpha) * z_prev[i] + rho_inv_vec[i] * y[i];\n'.format(n)
    update_z += 2 * ident + 'if( z[i] > u[i] ) z[i] = u[i];\n'
    update_z += 2 * ident + 'else if( z[i] < l[i] ) z[i] = l[i];\n'
    #update_z += 2 * ident + 'delta_x[i] = x[i] - x_prev[i];\n'
    update_z += ident + '}\n'

    update_y = '\n' + ident + '/* Updates y*/\n'
    update_y += ident + 'for(i = 0; i < {:}; i++){{\n'.format(m)
    update_y += 2 * ident + 'delta_y[i] = rho_vec[i] * (alpha * xz_tilde[i + {:}] + ( ((float)1.0f) - alpha) * z_prev[i] - z[i]);\n'.format(n)
    update_y += 2 * ident + 'y[i] += delta_y[i];\n'
    update_y += ident + '}\n'    

    #print(update_xzy_sig)
    #print(update_x)
    #print(update_z)
    #print(update_y)

    update_xzy += update_x + update_z + update_y + '}\n'

    auxil_c += update_xzy

    # Generates the header file
    auxil_h = '#ifndef AUXIL_ST_H\n#define AUXIL_ST_H\n'
    #qdldl_h += '\n' + lsolve_sig + '\n'
    #qdldl_h += ltsolve_sig + '\n'
    auxil_h +=  '\n' + update_xzy_sig + '\n'
    auxil_h += '\n#endif /* AUXIL_ST_H */\n'

    # Saves file
    if len(path) > 0:
        if path[-1] != '\\' or path[-1] != '/':
            path = path + '/'
    with open(path + fname + '.c', mode='w') as f:
        f.write(auxil_c)
    with open(path + fname + '.h', mode='w') as f:
        f.write(auxil_h)

if __name__ == '__main__':
    
    L_p = [0, 1, 2, 4, 6, 13, 19, 24, 28, 31, 33, 34, 34]
    L_i = [4, 11, 4, 11, 4, 11, 5, 6, 7, 8, 9, 10, 11, 6, 7, 8, 9, 10, 11,
          7, 8, 9, 10, 11, 8, 9, 10, 11, 9, 10, 11, 10, 11, 11]
    P = [8, 9, 11, 10, 0, 2, 3, 4, 5, 6, 7, 1]

    qdldl_fname = 'qdldl_st'
    auxil_fname = 'auxil_st'
    
    gen_qdldl_solve(qdldl_fname, L_p, L_i, P)
    gen_auxil(auxil_fname, 8, 4, 1.6)
