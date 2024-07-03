

def folding(h_input, response_matrix, h_output):
    for a in range(h_output.GetNbinsX()):
        for b in range(h_output.GetNbinsY()):
            val = 0.0
            err = 0.0
            index_x_out = a + h_output.GetNbinsX()*b
            for k in range(h_input.GetNbinsX()):
                for l in range(h_input.GetNbinsY()):
                    index_x_in = k + h_input.GetNbinsX()*l
                    val += h_input.GetBinContent(k+1, l+1) * response_matrix(index_x_out, index_x_in)
                    err += h_input.GetBinError(k+1, l+1)**2 * response_matrix(index_x_out, index_x_in)**2
            h_output.SetBinContent(a+1, b+1, val)
            h_output.SetBinError(a+1, b+1, math.sqrt(err))
    return h_output
