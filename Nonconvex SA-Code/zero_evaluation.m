G = [0 ; 0];
final_loss =0;
z = final_solution;
for k =1 : data.vali
    zero_SGradient;
    G = G + GR;
    final_loss = final_loss +obj;
end
G = G/data.vali;
gradient_final = norm(G)^2;
final_loss = final_loss/data.vali;