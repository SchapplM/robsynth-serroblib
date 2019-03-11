% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = qJ(2) + qJ(3);
t37 = qJ(4) + t41;
t30 = pkin(11) + t37;
t26 = sin(t30);
t27 = cos(t30);
t50 = t27 * pkin(5) + t26 * pkin(10);
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t23 = g(1) * t47 + g(2) * t44;
t3 = -g(3) * t27 + t23 * t26;
t48 = -pkin(8) - pkin(7);
t34 = sin(t41);
t68 = pkin(3) * t34;
t31 = sin(t37);
t67 = pkin(4) * t31;
t66 = pkin(5) * t26;
t65 = pkin(10) * t27;
t64 = g(3) * t26;
t43 = sin(qJ(2));
t62 = t43 * pkin(2);
t42 = sin(qJ(6));
t61 = t44 * t42;
t45 = cos(qJ(6));
t60 = t44 * t45;
t59 = t47 * t42;
t58 = t47 * t45;
t35 = cos(t41);
t29 = pkin(3) * t35;
t46 = cos(qJ(2));
t38 = t46 * pkin(2);
t57 = t29 + t38;
t40 = -pkin(9) + t48;
t32 = cos(t37);
t28 = pkin(4) * t32;
t56 = t28 + t50;
t55 = t28 + t57;
t17 = -t67 - t68;
t16 = t17 - t62;
t54 = t16 - t66;
t53 = t17 - t66;
t52 = t29 + t56;
t51 = -t66 - t67;
t22 = g(1) * t44 - g(2) * t47;
t5 = -g(3) * t32 + t23 * t31;
t7 = -g(3) * t35 + t23 * t34;
t49 = -g(3) * t46 + t23 * t43;
t36 = -qJ(5) + t40;
t33 = t38 + pkin(1);
t20 = t47 * t65;
t19 = t44 * t65;
t18 = pkin(1) + t57;
t15 = pkin(1) + t55;
t14 = t27 * t58 + t61;
t13 = -t27 * t59 + t60;
t12 = -t27 * t60 + t59;
t11 = t27 * t61 + t58;
t10 = t47 * t15;
t9 = t22 * t26;
t8 = g(3) * t34 + t23 * t35;
t6 = g(3) * t31 + t23 * t32;
t4 = t23 * t27 + t64;
t2 = t3 * t45;
t1 = t3 * t42;
t21 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t46, -t22 * t43, -t23, -g(1) * (-t44 * pkin(1) + t47 * pkin(7)) - g(2) * (t47 * pkin(1) + t44 * pkin(7)) 0, 0, 0, 0, 0, 0, t22 * t35, -t22 * t34, -t23, -g(1) * (-t44 * t33 - t47 * t48) - g(2) * (t47 * t33 - t44 * t48) 0, 0, 0, 0, 0, 0, t22 * t32, -t22 * t31, -t23, -g(1) * (-t44 * t18 - t47 * t40) - g(2) * (t47 * t18 - t44 * t40) 0, 0, 0, 0, 0, 0, t22 * t27, -t9, -t23, -g(1) * (-t44 * t15 - t47 * t36) - g(2) * (-t44 * t36 + t10) 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t9, -g(2) * t10 + (g(1) * t36 - g(2) * t50) * t47 + (-g(1) * (-t15 - t50) + g(2) * t36) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, g(3) * t43 + t23 * t46, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t49 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t57 - t23 * (-t62 - t68) 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t55 - t23 * t16, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t54 * t47 + t20) - g(2) * (t54 * t44 + t19) - g(3) * (t38 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * (t28 + t29) - t23 * t17, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t53 * t47 + t20) - g(2) * (t53 * t44 + t19) - g(3) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t51 * t47 + t20) - g(2) * (t51 * t44 + t19) - g(3) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t11 + t42 * t64, g(1) * t14 - g(2) * t12 + t45 * t64, 0, 0;];
taug_reg  = t21;
