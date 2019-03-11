% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t37 = qJ(2) + qJ(3);
t34 = qJ(4) + t37;
t28 = sin(t34);
t29 = cos(t34);
t42 = cos(qJ(5));
t30 = t42 * pkin(5) + pkin(4);
t38 = -qJ(6) - pkin(10);
t80 = -t28 * t38 + t29 * t30;
t79 = t29 * pkin(4) + t28 * pkin(10);
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t24 = g(1) * t44 + g(2) * t41;
t64 = t44 * t42;
t39 = sin(qJ(5));
t67 = t41 * t39;
t12 = t29 * t67 + t64;
t65 = t44 * t39;
t66 = t41 * t42;
t14 = -t29 * t65 + t66;
t70 = g(3) * t28;
t1 = -g(1) * t14 + g(2) * t12 + t39 * t70;
t7 = -g(3) * t29 + t24 * t28;
t45 = -pkin(8) - pkin(7);
t32 = sin(t37);
t78 = pkin(3) * t32;
t77 = pkin(4) * t28;
t76 = pkin(10) * t29;
t33 = cos(t37);
t27 = pkin(3) * t33;
t43 = cos(qJ(2));
t35 = t43 * pkin(2);
t62 = t27 + t35;
t20 = pkin(1) + t62;
t16 = t44 * t20;
t72 = g(2) * t16;
t60 = t27 + t79;
t36 = -pkin(9) + t45;
t59 = pkin(5) * t39 - t36;
t22 = t41 * t76;
t57 = -t41 * t77 + t22;
t23 = t44 * t76;
t56 = -t44 * t77 + t23;
t55 = t27 + t80;
t54 = -t77 - t78;
t52 = g(1) * t41 - g(2) * t44;
t50 = t28 * t30 + t29 * t38;
t49 = t50 * t41;
t48 = t50 * t44;
t9 = -g(3) * t33 + t24 * t32;
t40 = sin(qJ(2));
t46 = -g(3) * t43 + t24 * t40;
t31 = t35 + pkin(1);
t21 = -t40 * pkin(2) - t78;
t18 = t44 * t21;
t17 = t41 * t21;
t15 = t29 * t64 + t67;
t13 = -t29 * t66 + t65;
t11 = t52 * t28;
t10 = g(3) * t32 + t24 * t33;
t8 = t24 * t29 + t70;
t6 = t7 * t42;
t5 = t7 * t39;
t4 = -g(1) * t13 - g(2) * t15;
t3 = -g(1) * t12 - g(2) * t14;
t2 = g(1) * t15 - g(2) * t13 + t42 * t70;
t19 = [0, 0, 0, 0, 0, 0, t52, t24, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t43, -t52 * t40, -t24, -g(1) * (-t41 * pkin(1) + t44 * pkin(7)) - g(2) * (t44 * pkin(1) + t41 * pkin(7)) 0, 0, 0, 0, 0, 0, t52 * t33, -t52 * t32, -t24, -g(1) * (-t41 * t31 - t44 * t45) - g(2) * (t44 * t31 - t41 * t45) 0, 0, 0, 0, 0, 0, t52 * t29, -t11, -t24, -g(1) * (-t41 * t20 - t44 * t36) - g(2) * (-t41 * t36 + t16) 0, 0, 0, 0, 0, 0, t4, t3, t11, -t72 + (g(1) * t36 - g(2) * t79) * t44 + (-g(1) * (-t20 - t79) + g(2) * t36) * t41, 0, 0, 0, 0, 0, 0, t4, t3, t11, -t72 + (-g(1) * t59 - g(2) * t80) * t44 + (-g(1) * (-t20 - t80) - g(2) * t59) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(3) * t40 + t24 * t43, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t46 * pkin(2), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(3) * t62 - t24 * t21, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t18 + t56) - g(2) * (t17 + t57) - g(3) * (t35 + t60) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t18 - t48) - g(2) * (t17 - t49) - g(3) * (t35 + t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t9 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t54 * t44 + t23) - g(2) * (t54 * t41 + t22) - g(3) * t60, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t55 + t24 * (t50 + t78); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t56 - g(2) * t57 - g(3) * t79, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, g(1) * t48 + g(2) * t49 - g(3) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t19;
