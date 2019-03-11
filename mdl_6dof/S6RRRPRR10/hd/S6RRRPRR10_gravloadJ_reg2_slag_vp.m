% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t51 = sin(qJ(1));
t55 = cos(qJ(1));
t28 = g(1) * t55 + g(2) * t51;
t50 = sin(qJ(2));
t120 = t28 * t50;
t54 = cos(qJ(2));
t128 = g(3) * t54 - t120;
t49 = sin(qJ(3));
t52 = cos(qJ(5));
t48 = sin(qJ(5));
t53 = cos(qJ(3));
t98 = t48 * t53;
t67 = -t49 * t52 + t98;
t63 = g(3) * t67;
t90 = t55 * t49;
t25 = -t51 * t53 + t54 * t90;
t89 = t55 * t53;
t26 = t51 * t49 + t54 * t89;
t7 = t25 * t52 - t26 * t48;
t93 = t51 * t54;
t23 = t49 * t93 + t89;
t92 = t53 * t54;
t24 = t51 * t92 - t90;
t80 = -t23 * t52 + t24 * t48;
t127 = -g(1) * t7 + g(2) * t80 + t50 * t63;
t105 = g(3) * t50;
t47 = qJ(5) + qJ(6);
t40 = sin(t47);
t41 = cos(t47);
t3 = t25 * t41 - t26 * t40;
t69 = t40 * t53 - t41 * t49;
t81 = -t23 * t41 + t24 * t40;
t1 = -g(1) * t3 + g(2) * t81 + t69 * t105;
t99 = t48 * t49;
t66 = t52 * t53 + t99;
t102 = t23 * t48;
t70 = t24 * t52 + t102;
t8 = t25 * t48 + t26 * t52;
t119 = g(1) * t8 + g(2) * t70 + t66 * t105;
t4 = t25 * t40 + t26 * t41;
t68 = t40 * t49 + t41 * t53;
t71 = t23 * t40 + t24 * t41;
t2 = g(1) * t4 + g(2) * t71 + t68 * t105;
t113 = -pkin(3) - pkin(4);
t112 = g(1) * t51;
t42 = t50 * pkin(8);
t44 = t54 * pkin(2);
t39 = t52 * pkin(5) + pkin(4);
t103 = -pkin(3) - t39;
t97 = t49 * t50;
t96 = t50 * t51;
t95 = t50 * t53;
t94 = t50 * t55;
t91 = t54 * t55;
t56 = -pkin(10) - pkin(9);
t88 = t55 * t56;
t87 = t44 + t42;
t86 = t55 * pkin(1) + t51 * pkin(7);
t85 = qJ(4) * t49;
t84 = -pkin(1) - t44;
t83 = -pkin(2) - t85;
t82 = pkin(5) * t48 + qJ(4);
t19 = t23 * pkin(3);
t79 = t24 * qJ(4) - t19;
t21 = t25 * pkin(3);
t78 = t26 * qJ(4) - t21;
t77 = pkin(3) * t92 + t54 * t85 + t87;
t76 = pkin(2) * t91 + pkin(8) * t94 + t86;
t45 = t55 * pkin(7);
t75 = -t24 * pkin(3) - t23 * qJ(4) + t45;
t74 = t26 * pkin(3) + t76;
t73 = g(1) * t23 - g(2) * t25;
t72 = -g(2) * t55 + t112;
t61 = t25 * qJ(4) + t74;
t60 = (t84 - t42) * t112;
t6 = g(1) * t25 + g(2) * t23 + g(3) * t97;
t59 = g(1) * t26 + g(2) * t24 + g(3) * t95;
t34 = pkin(8) * t91;
t31 = pkin(8) * t93;
t29 = qJ(4) * t95;
t27 = g(1) * t96 - g(2) * t94;
t15 = t28 * t54 + t105;
t11 = t128 * t53;
t10 = t128 * t49;
t9 = g(1) * t24 - g(2) * t26;
t5 = [0, 0, 0, 0, 0, 0, t72, t28, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t54, -t27, -t28, -g(1) * (-t51 * pkin(1) + t45) - g(2) * t86, 0, 0, 0, 0, 0, 0, t9, -t73, t27, -g(1) * t45 - g(2) * t76 - t60, 0, 0, 0, 0, 0, 0, t9, t27, t73, -g(1) * t75 - g(2) * t61 - t60, 0, 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t8, -g(1) * t80 - g(2) * t7, -t27, -g(1) * (-t24 * pkin(4) + t75) - g(2) * (t26 * pkin(4) - pkin(9) * t94 + t61) - ((-pkin(8) + pkin(9)) * t50 + t84) * t112, 0, 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t4, -g(1) * t81 - g(2) * t3, -t27, -g(1) * (-pkin(5) * t102 - t24 * t39 + t75) - g(2) * (t25 * t82 + t26 * t39 + t50 * t88 + t74) - ((-pkin(8) - t56) * t50 + t84) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, t15, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, -t15, -g(1) * (-pkin(2) * t94 + t34) - g(2) * (-pkin(2) * t96 + t31) - g(3) * t87, 0, 0, 0, 0, 0, 0, -t11, -t15, -t10, -g(1) * t34 - g(2) * t31 - g(3) * t77 + (pkin(3) * t53 - t83) * t120, 0, 0, 0, 0, 0, 0, -t128 * t66, -t120 * t67 + t54 * t63, t15, -g(1) * (-pkin(9) * t91 + t34) - g(2) * (-pkin(9) * t93 + t31) - g(3) * (pkin(4) * t92 + t77) + (g(3) * pkin(9) + t28 * (-t113 * t53 - t83)) * t50, 0, 0, 0, 0, 0, 0, -t128 * t68, t128 * t69, t15, -g(1) * (t54 * t88 + t34) - g(2) * (t56 * t93 + t31) - g(3) * (pkin(5) * t54 * t99 + t39 * t92 + t77) + (-g(3) * t56 + t28 * (-t103 * t53 + t49 * t82 + pkin(2))) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t59, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, -t59, -g(1) * t78 - g(2) * t79 - g(3) * (-pkin(3) * t97 + t29) 0, 0, 0, 0, 0, 0, -t127, -t119, 0, -g(1) * (-t25 * pkin(4) + t78) - g(2) * (-t23 * pkin(4) + t79) - g(3) * (t113 * t97 + t29) 0, 0, 0, 0, 0, 0, -t1, -t2, 0, -g(1) * (-t25 * t39 + t26 * t82 - t21) - g(2) * (-t23 * t39 + t24 * t82 - t19) - g(3) * t29 - (pkin(5) * t98 + t103 * t49) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t119, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t127 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t5;
