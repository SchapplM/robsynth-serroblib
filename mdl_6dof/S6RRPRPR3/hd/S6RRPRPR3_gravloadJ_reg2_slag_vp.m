% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t22 = g(1) * t48 + g(2) * t45;
t40 = qJ(2) + pkin(10);
t33 = cos(t40);
t46 = cos(qJ(4));
t60 = t48 * t46;
t43 = sin(qJ(4));
t67 = t45 * t43;
t14 = t33 * t67 + t60;
t61 = t48 * t43;
t66 = t45 * t46;
t16 = -t33 * t61 + t66;
t31 = sin(t40);
t75 = g(3) * t31;
t82 = -g(1) * t16 + g(2) * t14 + t43 * t75;
t3 = -g(3) * t33 + t22 * t31;
t44 = sin(qJ(2));
t81 = pkin(2) * t44;
t47 = cos(qJ(2));
t36 = t47 * pkin(2);
t29 = t36 + pkin(1);
t23 = t48 * t29;
t77 = g(2) * t23;
t73 = t43 * pkin(4);
t39 = qJ(4) + pkin(11);
t30 = sin(t39);
t19 = pkin(5) * t30 + t73;
t72 = t19 * t33;
t34 = qJ(6) + t39;
t26 = sin(t34);
t71 = t45 * t26;
t27 = cos(t34);
t70 = t45 * t27;
t69 = t45 * t30;
t32 = cos(t39);
t68 = t45 * t32;
t65 = t48 * t26;
t64 = t48 * t27;
t63 = t48 * t30;
t62 = t48 * t32;
t41 = -qJ(5) - pkin(8);
t42 = -qJ(3) - pkin(7);
t59 = t19 - t42;
t35 = t46 * pkin(4);
t20 = pkin(5) * t32 + t35;
t57 = -t42 + t73;
t56 = pkin(3) * t33 + pkin(8) * t31;
t21 = g(1) * t45 - g(2) * t48;
t18 = pkin(3) + t20;
t38 = -pkin(9) + t41;
t55 = t18 * t33 - t31 * t38;
t28 = t35 + pkin(3);
t54 = t28 * t33 - t31 * t41;
t49 = -g(3) * t47 + t22 * t44;
t17 = t33 * t60 + t67;
t15 = -t33 * t66 + t61;
t13 = t21 * t31;
t12 = t33 * t62 + t69;
t11 = -t33 * t63 + t68;
t10 = -t33 * t68 + t63;
t9 = t33 * t69 + t62;
t8 = t33 * t64 + t71;
t7 = -t33 * t65 + t70;
t6 = -t33 * t70 + t65;
t5 = t33 * t71 + t64;
t4 = t22 * t33 + t75;
t2 = g(1) * t8 - g(2) * t6 + t27 * t75;
t1 = -g(1) * t7 + g(2) * t5 + t26 * t75;
t24 = [0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t47, -t21 * t44, -t22, -g(1) * (-pkin(1) * t45 + pkin(7) * t48) - g(2) * (pkin(1) * t48 + pkin(7) * t45) 0, 0, 0, 0, 0, 0, t21 * t33, -t13, -t22, -g(1) * (-t45 * t29 - t48 * t42) - g(2) * (-t45 * t42 + t23) 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -t77 + (g(1) * t42 - g(2) * t56) * t48 + (-g(1) * (-t29 - t56) + g(2) * t42) * t45, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t77 + (-g(1) * t57 - g(2) * t54) * t48 + (-g(1) * (-t29 - t54) - g(2) * t57) * t45, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t13, -t77 + (-g(1) * t59 - g(2) * t55) * t48 + (-g(1) * (-t29 - t55) - g(2) * t59) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, g(3) * t44 + t22 * t47, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t49 * pkin(2), 0, 0, 0, 0, 0, 0, t3 * t46, -t3 * t43, -t4, -g(3) * (t36 + t56) + t22 * (pkin(3) * t31 - pkin(8) * t33 + t81) 0, 0, 0, 0, 0, 0, t3 * t32, -t3 * t30, -t4, -g(3) * (t36 + t54) + t22 * (t28 * t31 + t33 * t41 + t81) 0, 0, 0, 0, 0, 0, t3 * t27, -t3 * t26, -t4, -g(3) * (t36 + t55) + t22 * (t18 * t31 + t33 * t38 + t81); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, g(1) * t17 - g(2) * t15 + t46 * t75, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t30 * t75, g(1) * t12 - g(2) * t10 + t32 * t75, 0, t82 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t20 * t45 - t48 * t72) - g(2) * (-t20 * t48 - t45 * t72) + t19 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t24;
