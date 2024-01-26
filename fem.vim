
function s:linspace(x0, x1, n)
    let xs = []
    let h = (a:x1 - a:x0) / (a:n - 1)

    for i in range(a:n)
        call add(xs, a:x0 + i * h)
    endfor

    return xs
endfunction


function s:zeros(n)
    return map(range(a:n), 0)
endfunction


function s:zeros2d(n, m)
    return map(range(a:n), {k -> s:zeros(a:m)})
endfunction


function s:solve(matrix, rhs)
    let n = len(a:matrix)
    for i in range(n - 1)
        let v = a:matrix[i][1]
        let a:matrix[i][1] = 1.0
        let a:matrix[i][2] /= v
        let a:rhs[i] /= v
        let a = a:matrix[i + 1][0]
        let a:matrix[i + 1][0] = 0.0
        let a:matrix[i + 1][1] -= a * a:matrix[i][2]
        let a:rhs[i + 1] -= a * a:rhs[i]
    endfor
    let a:rhs[n - 1] /= a:matrix[n - 1][1]
    let a:matrix[n - 1][1] = 1.0

    for i in range(n - 2, 0, -1)
        let a:rhs[i] -= a:matrix[i][2] * a:rhs[i + 1]
    endfor
endfunction


function s:basis_eval(x) dict
    let k = self.idx
    let xs = self.mesh
    if k > 0 && a:x < xs[k - 1]
        return 0.0
    elseif k < len(xs) - 1 && a:x >= xs[k + 1]
        return 0.0
    elseif a:x <= xs[k]
        return (a:x - xs[k - 1]) / (xs[k] - xs[k - 1])
    else
        return (xs[k + 1] - a:x) / (xs[k + 1] - xs[k])
    endif
endfunction


function s:basis(mesh, k)
    return {'idx': a:k, 'mesh': a:mesh, 'eval': function("s:basis_eval")}
endfunction


function s:basis_eval_der(x) dict
    let k = self.idx
    let xs = self.mesh
    if k > 0 && a:x < xs[k - 1]
        return 0.0
    elseif k < len(xs) - 1 && a:x >= xs[k + 1]
        return 0.0
    elseif a:x <= xs[k]
        return 1.0 / (xs[k] - xs[k - 1])
    else
        return -1.0 / (xs[k + 1] - xs[k])
    endif
endfunction


function s:basis_der(mesh, k)
    return {'idx': a:k, 'mesh': a:mesh, 'eval': function("s:basis_eval_der")}
endfunction


function s:find_interval(mesh, x)
    let n = len(a:mesh)
    for i in range(n)
        if a:mesh[i] > a:x
            return i
        endif
    endfor
    return n - 1
endfunction


let s:quad_data = [
\ [-0.86113631152, 0.34785484513],
\ [-0.33998104350, 0.65214515486],
\ [ 0.33998104350, 0.65214515486],
\ [ 0.86113631152, 0.34785484513]
\ ]


function s:integrate(a, b, f)
    let val = 0.0
    for [t, w] in s:quad_data
        let s = (1 + t) / 2.0
        let x = a:a + s * (a:b - a:a)
        let val += a:f(x) * w
    endfor
    return val
endfunction


function s:combination_eval(x) dict
    let idx = s:find_interval(self.mesh, a:x)
    let val = 0
    for i in range(idx - 1, idx)
        let e_i = s:basis(self.mesh, i)
        let val += self.coeffs[i] * e_i.eval(a:x)
    endfor
    return val
endfunction

function s:fun_from_coeffs(coeffs, mesh)
    return {'coeffs': a:coeffs, 'mesh': a:mesh, 'eval': function('s:combination_eval')}
endfunction


function s:fem(N)
    let pi = 3.1415926535
    let Forcing = {x -> 4 * pi * pi * sin(2*pi*x) + 6}
    let u_left = -0.7
    let u_right = 0.3

    let dim = a:N + 1
    let mesh = s:linspace(0.0, 1.0, dim)

    function B(i, j) closure
        let de_i = s:basis_der(mesh, a:i)
        let de_j = s:basis_der(mesh, a:j)
        return {x -> de_i.eval(x) * de_j.eval(x)}
    endfunction

    function L(i) closure
        let e_i = s:basis(mesh, a:i)
        return {x -> e_i.eval(x) * Forcing(x)}
    endfunction

    let matrix = s:zeros2d(dim, 3)
    let rhs = s:zeros(dim)

    let matrix[0][1] = 1.0
    let rhs[0] = u_left

    for i in range(1, a:N - 1)
        let e_i = s:basis(mesh, i)

        let matrix[i][0] = s:integrate(mesh[i - 1], mesh[i], B(i, i - 1))
        let matrix[i][1] = s:integrate(mesh[i - 1], mesh[i], B(i, i))
                        \+ s:integrate(mesh[i], mesh[i + 1], B(i, i))
        let matrix[i][2] = s:integrate(mesh[i], mesh[i + 1], B(i, i + 1))

        let rhs[i] = s:integrate(mesh[i - 1], mesh[i], L(i))
                 \ + s:integrate(mesh[i], mesh[i + 1], L(i))
    endfor

    let rhs[a:N] = u_right
    let matrix[a:N][1] = 1.0

    call s:solve(matrix, rhs)

    " solution function
    let u_h = s:fun_from_coeffs(rhs, mesh)

    " Visualization
    let s:res_x = 250
    let s:res_y = 65

    let s:max_y = 1.2
    let s:min_y = -1.2
    let s:offset_y = 5

    set virtualedit=all
    setlocal buftype=nofile

    " prepare enough empty lines
    for i in range(s:res_y + s:offset_y)
        call append(0, "")
    endfor

    function s:y_to_line(y)
        let y_norm = (s:max_y - a:y) / (s:max_y - s:min_y)
        return s:offset_y + float2nr(y_norm * s:res_y)
    endfunction

    let &textwidth = s:res_x
    for i in range(s:res_x)
        let x = (1.0 * i) / s:res_x
        let y = u_h.eval(x)
        let j = s:y_to_line(y)
        call cursor(j, i + 1)
        normal r@
    endfor

    call cursor(s:y_to_line(0.0), 1)
endfunction

call s:fem(100)
