
if ~exist('A_new','var')
    [A_new, b_vec_new] = remove_null_row(A,b_vec);
end

figure(10)
hold on;
for i = 1:size(A_new,1)
    if norm(A_new(i,:))<1e-3
        break;
    end
    
    if abs(A_new(i,1))<1e-3
        plot([-5,5],[(5*A_new(i,1)-b_vec_new(i))/A_new(i,2),(-5*A_new(i,1)-b_vec_new(i))/A_new(i,2)]);
        plot([0,0+A_new(i,1)],[-b_vec_new(i)/A_new(i,2),-b_vec_new(i)/A_new(i,2)+A_new(i,2)],'k->')
    else
        plot([(5*A_new(i,2)-b_vec_new(i))/A_new(i,1),(-5*A_new(i,2)-b_vec_new(i))/A_new(i,1)],[-5,5]);
        plot([-b_vec_new(i)/A_new(i,1),-b_vec_new(i)/A_new(i,1)+A_new(i,1)],[0,0+A_new(i,2)],'k->')
    end
    
end
legend;
gain = 1.2;
plot(gain*[-1,1],gain*[-1,-1],'k--');
plot(gain*[-1,1],gain*[1,1],'k--');
plot(gain*[-1,-1],gain*[-1,1],'k--');
plot(gain*[1,1],gain*[-1,1],'k--');   
axis equal
% xlim([-2,2]);
% ylim([-2,2]);


% if exist('A_new','var')
%     clear A_new;
%     clear b_vec_new;
% end

% function [A_new, b_vec_new] = remove_null_row(A,b_vec)
%     for i = 1:size(A,1)
%         if norm(A(i,:))<1e-3 && b_vec(i)<-1e-3
%             error('The %th row is problamatic',i)
%         end
%     end
%     row_A =vecnorm(A,2,2)<1e-3.*ones(size(A,1),1);
%     row_b =abs(b_vec) < 1e-3;
%     row_to_be_removed = row_A;
%    
%     A(row_to_be_removed,:) = [];
%     b_vec(row_to_be_removed,:) = [];
%     
%     A_new = A;
%     b_vec_new = b_vec;
% end
