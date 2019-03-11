% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t27 = g(1) * t52 + g(2) * t49;
t48 = sin(qJ(2));
t90 = t27 * t48;
t47 = sin(qJ(6));
t71 = pkin(10) + qJ(5);
t37 = sin(t71);
t51 = cos(qJ(2));
t67 = cos(t71);
t15 = t48 * t37 + t51 * t67;
t58 = t48 * t67;
t16 = -t51 * t37 + t58;
t5 = t16 * t49;
t78 = t51 * t52;
t7 = t37 * t78 - t52 * t58;
t55 = g(1) * t7 - g(2) * t5 + g(3) * t15;
t92 = t55 * t47;
t50 = cos(qJ(6));
t91 = t55 * t50;
t38 = t48 * qJ(3);
t74 = t51 * pkin(2) + t38;
t88 = pkin(2) * t48;
t87 = g(1) * t49;
t84 = g(3) * t16;
t83 = t51 * pkin(3);
t45 = cos(pkin(10));
t34 = t45 * pkin(4) + pkin(3);
t82 = pkin(2) + t34;
t44 = sin(pkin(10));
t81 = t48 * t44;
t79 = t51 * t44;
t72 = qJ(3) * t51;
t30 = t49 * t72;
t70 = pkin(4) * t79;
t77 = t49 * t70 + t30;
t32 = t52 * t72;
t76 = t52 * t70 + t32;
t41 = t52 * pkin(7);
t46 = -pkin(8) - qJ(4);
t75 = t52 * t46 + t41;
t73 = t52 * pkin(1) + t49 * pkin(7);
t29 = pkin(4) * t81;
t69 = t51 * t34 + t29 + t74;
t68 = pkin(2) * t78 + t52 * t38 + t73;
t6 = t15 * t49;
t66 = t5 * pkin(5) + t6 * pkin(9);
t8 = t15 * t52;
t65 = -t7 * pkin(5) + t8 * pkin(9);
t64 = -g(1) * t5 - g(2) * t7;
t63 = -t15 * pkin(5) + t16 * pkin(9);
t26 = -g(2) * t52 + t87;
t62 = t6 * t47 - t52 * t50;
t61 = t52 * t47 + t6 * t50;
t60 = -t48 * t45 + t79;
t59 = t51 * t45 + t81;
t57 = -pkin(1) - t74;
t56 = t52 * t29 + t34 * t78 + t49 * t46 + t68;
t2 = g(1) * t8 + g(2) * t6 + t84;
t54 = t82 * t90;
t53 = (-pkin(1) - t82 * t51 + (-pkin(4) * t44 - qJ(3)) * t48) * t87;
t18 = t26 * t51;
t17 = t26 * t48;
t14 = t59 * t52;
t13 = t60 * t52;
t12 = t59 * t49;
t11 = t60 * t49;
t10 = g(3) * t48 + t27 * t51;
t9 = -g(3) * t51 + t90;
t4 = -t49 * t47 + t8 * t50;
t3 = -t8 * t47 - t49 * t50;
t1 = [0, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t27, -g(1) * (-t49 * pkin(1) + t41) - g(2) * t73, 0, 0, 0, 0, 0, 0, t18, -t27, t17, -g(1) * t41 - g(2) * t68 - t57 * t87, 0, 0, 0, 0, 0, 0, g(1) * t12 - g(2) * t14, -g(1) * t11 + g(2) * t13, t27, -g(1) * (-t52 * qJ(4) + t41) - g(2) * (pkin(3) * t78 + t68) + (-g(1) * (t57 - t83) + g(2) * qJ(4)) * t49, 0, 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, -t64, t27, -g(1) * t75 - g(2) * t56 - t53, 0, 0, 0, 0, 0, 0, g(1) * t61 - g(2) * t4, -g(1) * t62 - g(2) * t3, t64, -g(1) * (-t6 * pkin(5) + t5 * pkin(9) + t75) - g(2) * (t8 * pkin(5) + t7 * pkin(9) + t56) - t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t10, -g(1) * (-t52 * t88 + t32) - g(2) * (-t49 * t88 + t30) - g(3) * t74, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t59, -g(1) * t14 - g(2) * t12 + g(3) * t60, 0, -g(1) * t32 - g(2) * t30 - g(3) * (t74 + t83) + (pkin(2) + pkin(3)) * t90, 0, 0, 0, 0, 0, 0, -t55, -t2, 0, -g(1) * t76 - g(2) * t77 - g(3) * t69 + t54, 0, 0, 0, 0, 0, 0, -t91, t92, t2, -g(1) * (-t65 + t76) - g(2) * (-t66 + t77) - g(3) * (-t63 + t69) + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t2, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t92, -t2, -g(1) * t65 - g(2) * t66 - g(3) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t62 + t47 * t84, g(1) * t4 + g(2) * t61 + t50 * t84, 0, 0;];
taug_reg  = t1;
