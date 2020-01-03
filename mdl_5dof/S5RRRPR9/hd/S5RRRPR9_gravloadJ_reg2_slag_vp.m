% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t39 = sin(qJ(1));
t42 = cos(qJ(1));
t21 = g(1) * t42 + g(2) * t39;
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t53 = t42 * t40;
t41 = cos(qJ(2));
t61 = t39 * t41;
t13 = t37 * t61 + t53;
t54 = t42 * t37;
t15 = t39 * t40 - t41 * t54;
t38 = sin(qJ(2));
t65 = g(3) * t38;
t71 = -g(1) * t15 + g(2) * t13 + t37 * t65;
t11 = -g(3) * t41 + t21 * t38;
t69 = g(1) * t39;
t63 = t37 * pkin(3);
t62 = t38 * t42;
t60 = t41 * t42;
t35 = qJ(3) + pkin(9);
t26 = sin(t35);
t19 = pkin(4) * t26 + t63;
t59 = t42 * t19;
t28 = qJ(5) + t35;
t23 = sin(t28);
t58 = t42 * t23;
t24 = cos(t28);
t57 = t42 * t24;
t56 = t42 * t26;
t27 = cos(t35);
t55 = t42 * t27;
t36 = -qJ(4) - pkin(7);
t30 = t40 * pkin(3);
t20 = pkin(4) * t27 + t30;
t52 = t42 * pkin(1) + t39 * pkin(6);
t50 = t41 * pkin(2) + t38 * pkin(7);
t48 = -g(2) * t42 + t69;
t18 = pkin(2) + t20;
t34 = -pkin(8) + t36;
t47 = t41 * t18 - t38 * t34;
t25 = t30 + pkin(2);
t45 = t41 * t25 - t38 * t36;
t31 = t42 * pkin(6);
t17 = t48 * t38;
t16 = t39 * t37 + t41 * t53;
t14 = -t40 * t61 + t54;
t12 = t21 * t41 + t65;
t10 = t39 * t26 + t41 * t55;
t9 = t39 * t27 - t41 * t56;
t8 = -t27 * t61 + t56;
t7 = t26 * t61 + t55;
t6 = t39 * t23 + t41 * t57;
t5 = t39 * t24 - t41 * t58;
t4 = -t24 * t61 + t58;
t3 = t23 * t61 + t57;
t2 = g(1) * t6 - g(2) * t4 + t24 * t65;
t1 = -g(1) * t5 + g(2) * t3 + t23 * t65;
t22 = [0, 0, 0, 0, 0, 0, t48, t21, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t41, -t17, -t21, -g(1) * (-t39 * pkin(1) + t31) - g(2) * t52, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t17, -g(1) * t31 - g(2) * (t50 * t42 + t52) - (-pkin(1) - t50) * t69, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t17, -g(1) * (pkin(3) * t54 + t31) - g(2) * (t25 * t60 - t36 * t62 + t52) + (-g(1) * (-pkin(1) - t45) - g(2) * t63) * t39, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t17, -g(1) * (t31 + t59) - g(2) * (t18 * t60 - t34 * t62 + t52) + (-g(1) * (-pkin(1) - t47) - g(2) * t19) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t40, -t11 * t37, -t12, -g(3) * t50 + t21 * (pkin(2) * t38 - pkin(7) * t41), 0, 0, 0, 0, 0, 0, t11 * t27, -t11 * t26, -t12, -g(3) * t45 + t21 * (t25 * t38 + t36 * t41), 0, 0, 0, 0, 0, 0, t11 * t24, -t11 * t23, -t12, -g(3) * t47 + t21 * (t18 * t38 + t34 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, g(1) * t16 - g(2) * t14 + t40 * t65, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t26 * t65, g(1) * t10 - g(2) * t8 + t27 * t65, 0, t71 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t39 * t20 - t41 * t59) - g(2) * (-t19 * t61 - t42 * t20) + t19 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t22;
