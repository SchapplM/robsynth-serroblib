% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t110 = cos(qJ(2));
t98 = cos(pkin(6));
t84 = t98 * t110;
t113 = -t59 * t58 + t61 * t84;
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t38 = -t110 * t95 - t58 * t97;
t77 = t98 * t95;
t78 = t98 * t97;
t62 = t110 * t78 - t58 * t77;
t18 = t59 * t38 + t61 * t62;
t53 = pkin(12) + qJ(6);
t51 = sin(t53);
t52 = cos(t53);
t100 = -t110 * t77 - t58 * t78;
t37 = -t110 * t97 + t58 * t95;
t17 = t61 * t100 + t59 * t37;
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t96 = sin(pkin(6));
t90 = t61 * t96;
t8 = -t17 * t60 - t57 * t90;
t112 = t18 * t52 + t8 * t51;
t111 = -t18 * t51 + t8 * t52;
t22 = t59 * t100 - t61 * t37;
t107 = t51 * t60;
t106 = t52 * t60;
t54 = sin(pkin(12));
t105 = t54 * t60;
t55 = cos(pkin(12));
t104 = t55 * t60;
t75 = t96 * t95;
t76 = t97 * t96;
t29 = -t110 * t76 + t58 * t75;
t83 = t110 * t96;
t46 = pkin(2) * t83;
t99 = -t29 * pkin(3) + t46;
t94 = pkin(5) * t54 + pkin(9);
t92 = t58 * t98;
t31 = pkin(2) * t92 + (-pkin(8) - qJ(3)) * t96;
t50 = t110 * pkin(2) + pkin(1);
t93 = -t59 * t31 + t61 * t50;
t91 = t59 * t96;
t89 = t113 * pkin(2);
t30 = t110 * t75 + t58 * t76;
t88 = t30 * pkin(9) + t99;
t87 = t22 * pkin(3) + t93;
t11 = t22 * t57 - t60 * t91;
t7 = -t17 * t57 + t60 * t90;
t86 = -g(1) * t7 + g(2) * t11;
t21 = t61 * t38 - t59 * t62;
t85 = g(1) * t18 - g(2) * t21;
t82 = pkin(4) * t60 + qJ(5) * t57;
t81 = -t61 * t31 - t59 * t50;
t49 = t55 * pkin(5) + pkin(4);
t56 = -pkin(10) - qJ(5);
t80 = t49 * t60 - t56 * t57;
t79 = t18 * pkin(3) + t89;
t73 = pkin(3) * t17 + t81;
t72 = -t21 * pkin(9) + t87;
t23 = t30 * t57 - t98 * t60;
t71 = g(1) * t11 + g(2) * t7 + g(3) * t23;
t12 = t22 * t60 + t57 * t91;
t24 = t30 * t60 + t98 * t57;
t70 = g(1) * t12 + g(2) * t8 + g(3) * t24;
t69 = -g(1) * t22 + g(2) * t17 - g(3) * t30;
t68 = g(1) * t21 + g(2) * t18 - g(3) * t29;
t67 = -t17 * pkin(9) + t79;
t66 = t18 * pkin(9) + t73;
t34 = -t61 * t58 - t59 * t84;
t65 = t34 * pkin(2);
t64 = t21 * pkin(3) + t65;
t63 = pkin(9) * t22 + t64;
t36 = -g(1) * t90 - g(2) * t91;
t35 = t61 * t110 - t59 * t92;
t33 = -t59 * t110 - t61 * t92;
t28 = -g(1) * t91 + g(2) * t90 - g(3) * t98;
t6 = t12 * t52 - t21 * t51;
t5 = -t12 * t51 - t21 * t52;
t3 = t68 * t57;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * t61, g(1) * t61 + g(2) * t59, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t33 - g(2) * t35, g(1) * t113 - g(2) * t34, t36, -g(1) * (-t59 * pkin(1) + pkin(8) * t90) - g(2) * (t61 * pkin(1) + pkin(8) * t91) 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t22, t85, t36, -g(1) * t81 - g(2) * t93, 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t12, t86, -t85, -g(1) * t66 - g(2) * t72, 0, 0, 0, 0, 0, 0, -g(1) * (t18 * t54 - t55 * t8) - g(2) * (t12 * t55 - t21 * t54) -g(1) * (t18 * t55 + t54 * t8) - g(2) * (-t12 * t54 - t21 * t55) -t86, -g(1) * (-pkin(4) * t8 - qJ(5) * t7 + t66) - g(2) * (t12 * pkin(4) + t11 * qJ(5) + t72) 0, 0, 0, 0, 0, 0, g(1) * t111 - g(2) * t6, -g(1) * t112 - g(2) * t5, -t86, -g(1) * (t94 * t18 - t49 * t8 + t56 * t7 + t73) - g(2) * (-t11 * t56 + t12 * t49 - t94 * t21 + t87); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t34 - g(2) * t113 - g(3) * t83, g(3) * t96 * t58 + g(1) * t35 - g(2) * t33, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t69, 0, -g(1) * t65 - g(2) * t89 - g(3) * t46, 0, 0, 0, 0, 0, 0, -t68 * t60, t3, t69, -g(1) * t63 - g(2) * t67 - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (t21 * t104 + t22 * t54) - g(2) * (t18 * t104 - t17 * t54) - g(3) * (-t29 * t104 + t30 * t54) -g(1) * (-t21 * t105 + t22 * t55) - g(2) * (-t18 * t105 - t17 * t55) - g(3) * (t29 * t105 + t30 * t55) -t3, -g(1) * (t82 * t21 + t63) - g(2) * (t82 * t18 + t67) - g(3) * (-t82 * t29 + t88) 0, 0, 0, 0, 0, 0, -g(1) * (t21 * t106 + t22 * t51) - g(2) * (t18 * t106 - t17 * t51) - g(3) * (-t29 * t106 + t30 * t51) -g(1) * (-t21 * t107 + t22 * t52) - g(2) * (-t18 * t107 - t17 * t52) - g(3) * (t29 * t107 + t30 * t52) -t3, -g(1) * (t80 * t21 + t22 * t94 + t64) - g(2) * (-t94 * t17 + t80 * t18 + t79) - g(3) * (-t80 * t29 + t94 * t30 + t99); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t70, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t55, -t71 * t54, -t70, -g(1) * (-t11 * pkin(4) + t12 * qJ(5)) - g(2) * (-t7 * pkin(4) + t8 * qJ(5)) - g(3) * (-t23 * pkin(4) + t24 * qJ(5)) 0, 0, 0, 0, 0, 0, t71 * t52, -t71 * t51, -t70, -g(1) * (-t11 * t49 - t12 * t56) - g(2) * (-t7 * t49 - t8 * t56) - g(3) * (-t23 * t49 - t24 * t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t112 - g(3) * (-t24 * t51 + t29 * t52) g(1) * t6 + g(2) * t111 - g(3) * (-t24 * t52 - t29 * t51) 0, 0;];
taug_reg  = t1;
