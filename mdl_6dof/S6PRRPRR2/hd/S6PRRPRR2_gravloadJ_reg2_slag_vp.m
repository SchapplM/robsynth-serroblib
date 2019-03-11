% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t62 = cos(pkin(6));
t36 = sin(pkin(6));
t40 = sin(qJ(2));
t70 = t36 * t40;
t83 = -t39 * t70 + t62 * t42;
t43 = cos(qJ(2));
t35 = sin(pkin(11));
t58 = t35 * t62;
t61 = cos(pkin(11));
t20 = -t40 * t58 + t61 * t43;
t69 = t36 * t42;
t82 = -t20 * t39 + t35 * t69;
t28 = t42 * pkin(3) + pkin(2);
t68 = t36 * t43;
t21 = t28 * t68;
t81 = g(3) * t21;
t80 = g(3) * t36;
t51 = t62 * t61;
t18 = t35 * t43 + t40 * t51;
t38 = sin(qJ(5));
t79 = t18 * t38;
t78 = t20 * t38;
t33 = qJ(3) + pkin(12);
t30 = cos(t33);
t34 = qJ(5) + qJ(6);
t31 = sin(t34);
t76 = t30 * t31;
t32 = cos(t34);
t75 = t30 * t32;
t74 = t30 * t38;
t41 = cos(qJ(5));
t73 = t30 * t41;
t72 = t30 * t43;
t71 = t35 * t36;
t37 = -qJ(4) - pkin(8);
t67 = t37 * t40;
t66 = t38 * t43;
t65 = t41 * t43;
t17 = t35 * t40 - t43 * t51;
t64 = -t17 * t28 - t18 * t37;
t19 = t61 * t40 + t43 * t58;
t63 = -t19 * t28 - t20 * t37;
t57 = t36 * t61;
t55 = t82 * pkin(3);
t29 = sin(t33);
t54 = pkin(4) * t30 + pkin(9) * t29;
t27 = t41 * pkin(5) + pkin(4);
t44 = -pkin(10) - pkin(9);
t53 = t27 * t30 - t29 * t44;
t52 = t83 * pkin(3);
t13 = -t29 * t70 + t62 * t30;
t7 = -t18 * t29 - t30 * t57;
t9 = -t20 * t29 + t30 * t71;
t50 = g(1) * t9 + g(2) * t7 + g(3) * t13;
t10 = t20 * t30 + t29 * t71;
t14 = t62 * t29 + t30 * t70;
t8 = t18 * t30 - t29 * t57;
t49 = g(1) * t10 + g(2) * t8 + g(3) * t14;
t48 = -t18 * t39 - t42 * t57;
t5 = -g(1) * t19 - g(2) * t17 + g(3) * t68;
t47 = g(1) * t20 + g(2) * t18 + g(3) * t70;
t46 = t48 * pkin(3);
t45 = -g(1) * (-t10 * t38 + t19 * t41) - g(2) * (t17 * t41 - t8 * t38) - g(3) * (-t14 * t38 - t36 * t65);
t4 = t5 * t29;
t2 = -g(1) * (-t10 * t32 - t19 * t31) - g(2) * (-t17 * t31 - t8 * t32) - g(3) * (-t14 * t32 + t31 * t68);
t1 = -g(1) * (-t10 * t31 + t19 * t32) - g(2) * (t17 * t32 - t8 * t31) - g(3) * (-t14 * t31 - t32 * t68);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t47, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t42, t5 * t39, -t47, -g(1) * (-t19 * pkin(2) + t20 * pkin(8)) - g(2) * (-t17 * pkin(2) + t18 * pkin(8)) - (pkin(2) * t43 + pkin(8) * t40) * t80, 0, 0, 0, 0, 0, 0, -t5 * t30, t4, -t47, -g(1) * t63 - g(2) * t64 - g(3) * (-t36 * t67 + t21) 0, 0, 0, 0, 0, 0, -g(1) * (-t19 * t73 + t78) - g(2) * (-t17 * t73 + t79) - (t30 * t65 + t38 * t40) * t80, -g(1) * (t19 * t74 + t20 * t41) - g(2) * (t17 * t74 + t18 * t41) - (-t30 * t66 + t40 * t41) * t80, -t4, -g(1) * (-t54 * t19 + t63) - g(2) * (-t54 * t17 + t64) - t81 - (t54 * t43 - t67) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (-t19 * t75 + t20 * t31) - g(2) * (-t17 * t75 + t18 * t31) - (t31 * t40 + t32 * t72) * t80, -g(1) * (t19 * t76 + t20 * t32) - g(2) * (t17 * t76 + t18 * t32) - (-t31 * t72 + t32 * t40) * t80, -t4, -g(1) * (pkin(5) * t78 - t53 * t19 + t63) - g(2) * (pkin(5) * t79 - t53 * t17 + t64) - t81 - (t53 * t43 + (pkin(5) * t38 - t37) * t40) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t48 - g(3) * t83, -g(1) * (-t20 * t42 - t39 * t71) - g(2) * (-t18 * t42 + t39 * t57) - g(3) * (-t62 * t39 - t40 * t69) 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49, 0, -g(1) * t55 - g(2) * t46 - g(3) * t52, 0, 0, 0, 0, 0, 0, -t50 * t41, t50 * t38, -t49, -g(1) * (t9 * pkin(4) + t10 * pkin(9) + t55) - g(2) * (t7 * pkin(4) + t8 * pkin(9) + t46) - g(3) * (t13 * pkin(4) + t14 * pkin(9) + t52) 0, 0, 0, 0, 0, 0, -t50 * t32, t50 * t31, -t49, -g(1) * (-t10 * t44 + t9 * t27 + t55) - g(2) * (t7 * t27 - t8 * t44 + t46) - g(3) * (t13 * t27 - t14 * t44 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -g(1) * (-t10 * t41 - t19 * t38) - g(2) * (-t17 * t38 - t8 * t41) - g(3) * (-t14 * t41 + t36 * t66) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t45 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
