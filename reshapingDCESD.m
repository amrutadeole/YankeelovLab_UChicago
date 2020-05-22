reshaped_DCE_SD_array = zeros(320,320,24,65);
for i = 1:24
    for j = 1:65
        reshaped_DCE_SD_array(:,:,i,j) = imresize(DCE_SD(:,:,i,j),0.8);
    end
end