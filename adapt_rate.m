function l_new = adapt_rate(l_old, grad_new, grad_old, alpha, beta, l_max, l_min)

product_grad = grad_new.*grad_old;
l_max = l_max.*ones(size(l_old));
l_min = l_min.*ones(size(l_old));
l_new = zeros(size(l_old));

max_ls = min([alpha*l_old l_max], [], 2);
min_ls = max([beta*l_old l_min], [], 2);

l_new(product_grad > 0) = max_ls(product_grad > 0);
l_new(product_grad < 0) = min_ls(product_grad < 0);
l_new(product_grad == 0) = l_old(product_grad == 0);

end