
    a = imread('cameraman.tif');
    figure('Name','Input image');
    imshow(a);


if isgray(a)
    
    dvalue=double(a);
    
    
    if isa(a,'uint8')
        img_dct=dct2(dvalue);
        img_pow=(img_dct).^2;
        img_pow=img_pow(:);
        [B,index]=sort(img_pow);
        B=flipud(B);
        index=flipud(index);
        
        compressed_dct=zeros(size(dvalue));
        for ii=1:coeff
            compressed_dct(index(ii))=img_dct(index(ii));
        end
        im=idct2(compressed_dct);
        im=uint8(im);
    end
end

    figure('Name','Output image');
    imshow(im);
     im = imread(im);