% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t37 = cos(qJ(4));
t25 = t37 * pkin(4) + pkin(3);
t33 = qJ(2) + qJ(3);
t28 = sin(t33);
t30 = cos(t33);
t40 = -pkin(9) - pkin(8);
t75 = t30 * t25 - t28 * t40;
t74 = t30 * pkin(3) + t28 * pkin(8);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t19 = g(1) * t39 + g(2) * t36;
t54 = t39 * t37;
t34 = sin(qJ(4));
t59 = t36 * t34;
t14 = t30 * t59 + t54;
t55 = t39 * t34;
t58 = t36 * t37;
t16 = -t30 * t55 + t58;
t64 = g(3) * t28;
t73 = -g(1) * t16 + g(2) * t14 + t34 * t64;
t7 = -g(3) * t30 + t19 * t28;
t35 = sin(qJ(2));
t72 = pkin(2) * t35;
t71 = pkin(3) * t28;
t70 = pkin(8) * t30;
t38 = cos(qJ(2));
t31 = t38 * pkin(2);
t26 = t31 + pkin(1);
t20 = t39 * t26;
t66 = g(2) * t20;
t32 = qJ(4) + qJ(5);
t27 = sin(t32);
t61 = t36 * t27;
t29 = cos(t32);
t60 = t36 * t29;
t57 = t39 * t27;
t56 = t39 * t29;
t41 = -pkin(7) - pkin(6);
t51 = pkin(4) * t34 - t41;
t49 = -t71 - t72;
t47 = g(1) * t36 - g(2) * t39;
t45 = t25 * t28 + t30 * t40;
t42 = -g(3) * t38 + t19 * t35;
t22 = t39 * t70;
t21 = t36 * t70;
t17 = t30 * t54 + t59;
t15 = -t30 * t58 + t55;
t13 = t47 * t28;
t12 = t30 * t56 + t61;
t11 = -t30 * t57 + t60;
t10 = -t30 * t60 + t57;
t9 = t30 * t61 + t56;
t8 = t19 * t30 + t64;
t6 = t7 * t37;
t5 = t7 * t34;
t4 = t7 * t29;
t3 = t7 * t27;
t2 = g(1) * t12 - g(2) * t10 + t29 * t64;
t1 = -g(1) * t11 + g(2) * t9 + t27 * t64;
t18 = [0, 0, 0, 0, 0, 0, t47, t19, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t38, -t47 * t35, -t19, -g(1) * (-t36 * pkin(1) + t39 * pkin(6)) - g(2) * (t39 * pkin(1) + t36 * pkin(6)), 0, 0, 0, 0, 0, 0, t47 * t30, -t13, -t19, -g(1) * (-t36 * t26 - t39 * t41) - g(2) * (-t36 * t41 + t20), 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -t66 + (g(1) * t41 - g(2) * t74) * t39 + (-g(1) * (-t26 - t74) + g(2) * t41) * t36, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t66 + (-g(1) * t51 - g(2) * t75) * t39 + (-g(1) * (-t26 - t75) - g(2) * t51) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, g(3) * t35 + t19 * t38, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t42 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t49 * t39 + t22) - g(2) * (t49 * t36 + t21) - g(3) * (t31 + t74), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * (t31 + t75) + t19 * (t45 + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-t39 * t71 + t22) - g(2) * (-t36 * t71 + t21) - g(3) * t74, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t75 + t19 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, g(1) * t17 - g(2) * t15 + t37 * t64, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t73 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t18;
