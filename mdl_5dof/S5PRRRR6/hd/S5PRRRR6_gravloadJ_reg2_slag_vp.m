% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t34 = g(1) * t23 + g(2) * t22;
t21 = qJ(2) + qJ(3);
t16 = sin(t21);
t18 = cos(t21);
t7 = -g(3) * t18 + t16 * t34;
t25 = sin(qJ(2));
t52 = pkin(2) * t25;
t51 = pkin(3) * t16;
t50 = pkin(7) * t18;
t47 = g(3) * t16;
t20 = qJ(4) + qJ(5);
t15 = sin(t20);
t45 = t22 * t15;
t17 = cos(t20);
t44 = t22 * t17;
t24 = sin(qJ(4));
t43 = t22 * t24;
t26 = cos(qJ(4));
t42 = t22 * t26;
t41 = t23 * t15;
t40 = t23 * t17;
t39 = t23 * t24;
t38 = t23 * t26;
t37 = t18 * pkin(3) + t16 * pkin(7);
t14 = t26 * pkin(4) + pkin(3);
t28 = -pkin(8) - pkin(7);
t36 = t18 * t14 - t16 * t28;
t35 = -t51 - t52;
t33 = t14 * t16 + t18 * t28;
t27 = cos(qJ(2));
t30 = -g(3) * t27 + t25 * t34;
t29 = -g(1) * (-t18 * t39 + t42) - g(2) * (-t18 * t43 - t38) + t24 * t47;
t19 = t27 * pkin(2);
t11 = t23 * t50;
t10 = t22 * t50;
t8 = t18 * t34 + t47;
t6 = t7 * t26;
t5 = t7 * t24;
t4 = t7 * t17;
t3 = t7 * t15;
t2 = -g(1) * (-t18 * t40 - t45) - g(2) * (-t18 * t44 + t41) + t17 * t47;
t1 = -g(1) * (-t18 * t41 + t44) - g(2) * (-t18 * t45 - t40) + t15 * t47;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t25 + t27 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t30 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t23 * t35 + t11) - g(2) * (t22 * t35 + t10) - g(3) * (t19 + t37), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * (t19 + t36) + t34 * (t33 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-t23 * t51 + t11) - g(2) * (-t22 * t51 + t10) - g(3) * t37, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t36 + t34 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -g(1) * (-t18 * t38 - t43) - g(2) * (-t18 * t42 + t39) + t26 * t47, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
