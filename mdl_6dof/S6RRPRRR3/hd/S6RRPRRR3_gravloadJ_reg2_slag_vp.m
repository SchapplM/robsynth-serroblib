% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t23 = g(1) * t47 + g(2) * t44;
t38 = qJ(2) + pkin(11);
t31 = cos(t38);
t45 = cos(qJ(4));
t62 = t47 * t45;
t42 = sin(qJ(4));
t69 = t44 * t42;
t15 = t31 * t69 + t62;
t63 = t47 * t42;
t68 = t44 * t45;
t17 = -t31 * t63 + t68;
t30 = sin(t38);
t77 = g(3) * t30;
t86 = -g(1) * t17 + g(2) * t15 + t42 * t77;
t40 = qJ(4) + qJ(5);
t33 = cos(t40);
t64 = t47 * t33;
t32 = sin(t40);
t71 = t44 * t32;
t10 = t31 * t71 + t64;
t65 = t47 * t32;
t70 = t44 * t33;
t12 = -t31 * t65 + t70;
t3 = -g(1) * t12 + g(2) * t10 + t32 * t77;
t51 = -g(3) * t31 + t23 * t30;
t48 = -pkin(9) - pkin(8);
t43 = sin(qJ(2));
t85 = pkin(2) * t43;
t46 = cos(qJ(2));
t36 = t46 * pkin(2);
t29 = t36 + pkin(1);
t24 = t47 * t29;
t79 = g(2) * t24;
t75 = t42 * pkin(4);
t20 = pkin(5) * t32 + t75;
t74 = t20 * t31;
t34 = qJ(6) + t40;
t26 = sin(t34);
t73 = t44 * t26;
t27 = cos(t34);
t72 = t44 * t27;
t67 = t47 * t26;
t66 = t47 * t27;
t41 = -qJ(3) - pkin(7);
t61 = t20 - t41;
t35 = t45 * pkin(4);
t21 = pkin(5) * t33 + t35;
t58 = -t41 + t75;
t57 = t31 * pkin(3) + t30 * pkin(8);
t22 = g(1) * t44 - g(2) * t47;
t19 = pkin(3) + t21;
t39 = -pkin(10) + t48;
t56 = t31 * t19 - t30 * t39;
t28 = t35 + pkin(3);
t55 = t31 * t28 - t30 * t48;
t49 = -g(3) * t46 + t23 * t43;
t18 = t31 * t62 + t69;
t16 = -t31 * t68 + t63;
t14 = t22 * t30;
t13 = t31 * t64 + t71;
t11 = -t31 * t70 + t65;
t9 = t31 * t66 + t73;
t8 = -t31 * t67 + t72;
t7 = -t31 * t72 + t67;
t6 = t31 * t73 + t66;
t5 = t23 * t31 + t77;
t4 = g(1) * t13 - g(2) * t11 + t33 * t77;
t2 = g(1) * t9 - g(2) * t7 + t27 * t77;
t1 = -g(1) * t8 + g(2) * t6 + t26 * t77;
t25 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t46, -t22 * t43, -t23, -g(1) * (-t44 * pkin(1) + t47 * pkin(7)) - g(2) * (t47 * pkin(1) + t44 * pkin(7)) 0, 0, 0, 0, 0, 0, t22 * t31, -t14, -t23, -g(1) * (-t44 * t29 - t47 * t41) - g(2) * (-t44 * t41 + t24) 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t14, -t79 + (g(1) * t41 - g(2) * t57) * t47 + (-g(1) * (-t29 - t57) + g(2) * t41) * t44, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, t14, -t79 + (-g(1) * t58 - g(2) * t55) * t47 + (-g(1) * (-t29 - t55) - g(2) * t58) * t44, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t14, -t79 + (-g(1) * t61 - g(2) * t56) * t47 + (-g(1) * (-t29 - t56) - g(2) * t61) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, g(3) * t43 + t23 * t46, 0, 0, 0, 0, 0, 0, 0, 0, t51, t5, 0, t49 * pkin(2), 0, 0, 0, 0, 0, 0, t51 * t45, -t51 * t42, -t5, -g(3) * (t36 + t57) + t23 * (pkin(3) * t30 - pkin(8) * t31 + t85) 0, 0, 0, 0, 0, 0, t51 * t33, -t51 * t32, -t5, -g(3) * (t36 + t55) + t23 * (t28 * t30 + t31 * t48 + t85) 0, 0, 0, 0, 0, 0, t51 * t27, -t51 * t26, -t5, -g(3) * (t36 + t56) + t23 * (t19 * t30 + t31 * t39 + t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, g(1) * t18 - g(2) * t16 + t45 * t77, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t86 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t44 * t21 - t47 * t74) - g(2) * (-t47 * t21 - t44 * t74) + t20 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t25;
