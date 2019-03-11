% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t39 = qJ(2) + qJ(3);
t35 = qJ(4) + t39;
t32 = qJ(5) + t35;
t26 = sin(t32);
t27 = cos(t32);
t69 = t27 * pkin(5) + t26 * pkin(11);
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t22 = g(1) * t45 + g(2) * t42;
t3 = -g(3) * t27 + t22 * t26;
t46 = -pkin(8) - pkin(7);
t33 = sin(t39);
t68 = pkin(3) * t33;
t29 = sin(t35);
t67 = pkin(4) * t29;
t66 = pkin(5) * t26;
t65 = pkin(11) * t27;
t64 = g(3) * t26;
t41 = sin(qJ(2));
t62 = t41 * pkin(2);
t40 = sin(qJ(6));
t61 = t42 * t40;
t43 = cos(qJ(6));
t60 = t42 * t43;
t59 = t45 * t40;
t58 = t45 * t43;
t34 = cos(t39);
t28 = pkin(3) * t34;
t44 = cos(qJ(2));
t37 = t44 * pkin(2);
t56 = t28 + t37;
t38 = -pkin(9) + t46;
t30 = cos(t35);
t25 = pkin(4) * t30;
t55 = t25 + t69;
t54 = t25 + t56;
t17 = -t67 - t68;
t16 = t17 - t62;
t53 = t16 - t66;
t52 = t17 - t66;
t51 = t28 + t55;
t50 = -t66 - t67;
t48 = g(1) * t42 - g(2) * t45;
t5 = -g(3) * t30 + t22 * t29;
t7 = -g(3) * t34 + t22 * t33;
t47 = -g(3) * t44 + t22 * t41;
t36 = -pkin(10) + t38;
t31 = t37 + pkin(1);
t20 = t45 * t65;
t19 = t42 * t65;
t18 = pkin(1) + t56;
t15 = pkin(1) + t54;
t14 = t27 * t58 + t61;
t13 = -t27 * t59 + t60;
t12 = -t27 * t60 + t59;
t11 = t27 * t61 + t58;
t10 = t45 * t15;
t9 = t48 * t26;
t8 = g(3) * t33 + t22 * t34;
t6 = g(3) * t29 + t22 * t30;
t4 = t22 * t27 + t64;
t2 = t3 * t43;
t1 = t3 * t40;
t21 = [0, 0, 0, 0, 0, 0, t48, t22, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t44, -t48 * t41, -t22, -g(1) * (-t42 * pkin(1) + t45 * pkin(7)) - g(2) * (t45 * pkin(1) + t42 * pkin(7)) 0, 0, 0, 0, 0, 0, t48 * t34, -t48 * t33, -t22, -g(1) * (-t42 * t31 - t45 * t46) - g(2) * (t45 * t31 - t42 * t46) 0, 0, 0, 0, 0, 0, t48 * t30, -t48 * t29, -t22, -g(1) * (-t42 * t18 - t45 * t38) - g(2) * (t45 * t18 - t42 * t38) 0, 0, 0, 0, 0, 0, t48 * t27, -t9, -t22, -g(1) * (-t42 * t15 - t45 * t36) - g(2) * (-t42 * t36 + t10) 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t9, -g(2) * t10 + (g(1) * t36 - g(2) * t69) * t45 + (-g(1) * (-t15 - t69) + g(2) * t36) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, g(3) * t41 + t22 * t44, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t47 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t56 - t22 * (-t62 - t68) 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t54 - t22 * t16, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t53 * t45 + t20) - g(2) * (t53 * t42 + t19) - g(3) * (t37 + t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * (t25 + t28) - t22 * t17, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t52 * t45 + t20) - g(2) * (t52 * t42 + t19) - g(3) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t50 * t45 + t20) - g(2) * (t50 * t42 + t19) - g(3) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t45 * t66 + t20) - g(2) * (-t42 * t66 + t19) - g(3) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t11 + t40 * t64, g(1) * t14 - g(2) * t12 + t43 * t64, 0, 0;];
taug_reg  = t21;
