reshaped_DCE_LD_array = zeros(320,320,24,30);
for i = 1:24
    for j = 1:30
        reshaped_DCE_LD_array(:,:,i,j) = imresize(DCE_LD(:,:,i,j),0.8);
    end
end
