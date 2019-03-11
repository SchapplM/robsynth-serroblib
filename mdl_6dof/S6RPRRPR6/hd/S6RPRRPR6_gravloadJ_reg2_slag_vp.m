% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t23 = g(1) * t47 + g(2) * t45;
t38 = pkin(10) + qJ(3);
t32 = cos(t38);
t46 = cos(qJ(4));
t58 = t47 * t46;
t44 = sin(qJ(4));
t65 = t45 * t44;
t14 = t32 * t65 + t58;
t59 = t47 * t44;
t64 = t45 * t46;
t16 = -t32 * t59 + t64;
t30 = sin(t38);
t73 = g(3) * t30;
t79 = -g(1) * t16 + g(2) * t14 + t44 * t73;
t3 = -g(3) * t32 + t23 * t30;
t41 = cos(pkin(10));
t28 = t41 * pkin(2) + pkin(1);
t21 = t47 * t28;
t75 = g(2) * t21;
t71 = t44 * pkin(4);
t39 = qJ(4) + pkin(11);
t31 = sin(t39);
t19 = pkin(5) * t31 + t71;
t70 = t19 * t32;
t34 = qJ(6) + t39;
t26 = sin(t34);
t69 = t45 * t26;
t27 = cos(t34);
t68 = t45 * t27;
t67 = t45 * t31;
t33 = cos(t39);
t66 = t45 * t33;
t63 = t47 * t26;
t62 = t47 * t27;
t61 = t47 * t31;
t60 = t47 * t33;
t42 = -qJ(5) - pkin(8);
t43 = -pkin(7) - qJ(2);
t57 = t19 - t43;
t35 = t46 * pkin(4);
t20 = pkin(5) * t33 + t35;
t55 = -t43 + t71;
t54 = t32 * pkin(3) + t30 * pkin(8);
t22 = g(1) * t45 - g(2) * t47;
t18 = pkin(3) + t20;
t37 = -pkin(9) + t42;
t52 = t32 * t18 - t30 * t37;
t29 = t35 + pkin(3);
t50 = t32 * t29 - t30 * t42;
t17 = t32 * t58 + t65;
t15 = -t32 * t64 + t59;
t13 = t22 * t30;
t12 = t32 * t60 + t67;
t11 = -t32 * t61 + t66;
t10 = -t32 * t66 + t61;
t9 = t32 * t67 + t60;
t8 = t32 * t62 + t69;
t7 = -t32 * t63 + t68;
t6 = -t32 * t68 + t63;
t5 = t32 * t69 + t62;
t4 = t23 * t32 + t73;
t2 = g(1) * t8 - g(2) * t6 + t27 * t73;
t1 = -g(1) * t7 + g(2) * t5 + t26 * t73;
t24 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t41, -t22 * sin(pkin(10)) -t23, -g(1) * (-t45 * pkin(1) + t47 * qJ(2)) - g(2) * (t47 * pkin(1) + t45 * qJ(2)) 0, 0, 0, 0, 0, 0, t22 * t32, -t13, -t23, -g(1) * (-t45 * t28 - t47 * t43) - g(2) * (-t45 * t43 + t21) 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -t75 + (g(1) * t43 - g(2) * t54) * t47 + (-g(1) * (-t28 - t54) + g(2) * t43) * t45, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t75 + (-g(1) * t55 - g(2) * t50) * t47 + (-g(1) * (-t28 - t50) - g(2) * t55) * t45, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t13, -t75 + (-g(1) * t57 - g(2) * t52) * t47 + (-g(1) * (-t28 - t52) - g(2) * t57) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t46, -t3 * t44, -t4, -g(3) * t54 + t23 * (pkin(3) * t30 - pkin(8) * t32) 0, 0, 0, 0, 0, 0, t3 * t33, -t3 * t31, -t4, -g(3) * t50 + t23 * (t29 * t30 + t32 * t42) 0, 0, 0, 0, 0, 0, t3 * t27, -t3 * t26, -t4, -g(3) * t52 + t23 * (t18 * t30 + t32 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, g(1) * t17 - g(2) * t15 + t46 * t73, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t31 * t73, g(1) * t12 - g(2) * t10 + t33 * t73, 0, t79 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t45 * t20 - t47 * t70) - g(2) * (-t47 * t20 - t45 * t70) + t19 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t24;
