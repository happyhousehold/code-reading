function [average_cost2 , average_cost_diff] = simulation(z1, z2, data)
s1 = z1(1);
big_s1 = z1(2);
s2 = z2(1);
big_s2 = z2(2);
invt1 = 60 ;
invt2 = 60 ;
%total_time = 300;
me = 0.1;
h = 1;
bck_cst = 5;
fix_ord_cst = 32;
var_ord_cst = 3;
total_cost1 = 0 ;
total_cost2 = 0 ;
cust_time = exprnd(me);
next_time = 0;
clock = 0;
lead_time = 1000;
order_time = 1;
next_event = '';
order1 = 0;
order2 = 0;
while clock <= data.time
    if strcmp(next_event, 'arriving_customer') == 1
        total_cost1 = total_cost1 + next_time*(h*max(invt1, 0) + bck_cst* max(-invt1, 0));
        total_cost2 = total_cost2 + next_time*(h*max(invt2, 0) + bck_cst* max(-invt2, 0));
        rnd_demand = rand;
        demand = 4;
        if rnd_demand <= 0.167
            demand = 1;
        elseif  rnd_demand <= .5
            demand = 2;
        elseif rnd_demand <= .833
            demand = 3;
        end
        invt1 = invt1 - demand;
        invt2 = invt2 - demand;
        cust_time = exprnd(me);
        %creating amount of next demand
    elseif strcmp(next_event, 'arriving_order') == 1
        total_cost1 = total_cost1 + next_time*(h*max(invt1, 0) + bck_cst* max(-invt1, 0));
        total_cost2 = total_cost2 + next_time*(h*max(invt2, 0) + bck_cst* max(-invt2, 0));
        invt1 = invt1 + order1;
        invt2 = invt2 + order2;
        lead_time = 1000;
    elseif strcmp(next_event, 'placing_order') == 1
        total_cost1 = total_cost1 + next_time*(h*max(invt1, 0) + bck_cst* max(-invt1, 0));
        total_cost2 = total_cost2 + next_time*(h*max(invt2, 0) + bck_cst* max(-invt2, 0));
        if (invt1 < s1) && (invt2 < s2)
            order1 = big_s1 - invt1;
            order2 = big_s2 - invt2;
            total_cost1 = total_cost1 +fix_ord_cst+var_ord_cst*order1;
            total_cost2 = total_cost2 +fix_ord_cst+var_ord_cst*order2;
            lead_time = unifrnd(0.5 , 1);
        elseif  invt1 < s1
            order1 = big_s1 - invt1;
            order2 = 0 ;
            total_cost1 = total_cost1 +fix_ord_cst+var_ord_cst*order1;
            lead_time = unifrnd(0.5 , 1);
        elseif invt2 < s2
            order2 = big_s2 - invt2;
            order1 = 0 ;
            total_cost2 = total_cost2 +fix_ord_cst+var_ord_cst*order2;
            lead_time = unifrnd(0.5 , 1);
        end
        order_time = 1;
    end
    next_time = min(cust_time , order_time);
    next_time = min(next_time , lead_time);
    if next_time == cust_time
        next_event = 'arriving_customer';
    elseif next_time == order_time
        next_event = 'placing_order';
    else
        next_event = 'arriving_order';
    end
    cust_time = cust_time - next_time;
    lead_time = lead_time - next_time;
    order_time = order_time - next_time;
    %plot(clock , total_cost2/clock);
    %plot(clock , invt);
    %hold all
    %plot(clock , total_cost1/clock);
    clock = clock + next_time;
end
average_cost1 = total_cost1/ data.time;
average_cost2 = total_cost2/ data.time;
average_cost_diff = (average_cost1-average_cost2)/data.mu;