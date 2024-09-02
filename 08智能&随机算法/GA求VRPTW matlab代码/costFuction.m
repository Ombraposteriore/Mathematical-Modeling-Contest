
%% 计算当前解的成本函数
%输入curr_vc                  每辆车所经过的顾客
%输入a,b                      顾客时间窗结束时间[a[i],b[i]]
%输入s                        对每个顾客的服务时间
%输入L                        仓库时间窗结束时间
%输入dist                     距离矩阵
%输入demands                  各个顾客需求量
%输入cap                      车辆最大载货量
%输出cost                      成本函数 f=TD+alpha*q+belta*w
function [cost]=costFuction(curr_vc,a,b,s,L,dist,demands,cap,alpha,belta)
[TD] = travel_distance(curr_vc,dist);
[q]=violateLoad(curr_vc,demands,cap);
[w]=violateTW(curr_vc,a,b,s,L,dist);
cost=TD+alpha*q+belta*w;
end

