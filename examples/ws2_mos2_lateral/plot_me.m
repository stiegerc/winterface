h3 = read_sparse_mat('H_3.bin');
h4 = read_sparse_mat('H_4.bin');
h5 = read_sparse_mat('H_5.bin');

close all
figure

subplot(1,3,1), spy(h3), title(nnz(h3)/numel(h3))
subplot(1,3,2), spy(h4), title(nnz(h4)/numel(h4))
subplot(1,3,3), spy(h5), title(nnz(h5)/numel(h5))