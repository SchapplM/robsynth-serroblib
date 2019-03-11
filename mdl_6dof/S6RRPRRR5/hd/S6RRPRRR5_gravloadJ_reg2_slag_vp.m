% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t103 = cos(pkin(6));
t117 = cos(qJ(2));
t87 = t103 * t117;
t124 = -t60 * t59 + t63 * t87;
t100 = sin(pkin(12));
t81 = t103 * t100;
t102 = cos(pkin(12));
t82 = t103 * t102;
t105 = -t117 * t81 - t59 * t82;
t40 = t100 * t59 - t102 * t117;
t20 = t105 * t63 + t60 * t40;
t58 = sin(qJ(4));
t62 = cos(qJ(4));
t101 = sin(pkin(6));
t94 = t63 * t101;
t11 = -t20 * t62 - t58 * t94;
t41 = -t100 * t117 - t102 * t59;
t65 = t117 * t82 - t59 * t81;
t21 = t60 * t41 + t63 * t65;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t123 = t11 * t57 + t21 * t61;
t122 = t11 * t61 - t21 * t57;
t56 = qJ(5) + qJ(6);
t54 = sin(t56);
t55 = cos(t56);
t121 = t11 * t54 + t21 * t55;
t120 = t11 * t55 - t21 * t54;
t79 = t101 * t100;
t80 = t102 * t101;
t33 = t117 * t79 + t59 * t80;
t27 = t103 * t58 + t33 * t62;
t32 = -t117 * t80 + t59 * t79;
t25 = t105 * t60 - t63 * t40;
t95 = t60 * t101;
t15 = t25 * t62 + t58 * t95;
t24 = t63 * t41 - t60 * t65;
t8 = -t15 * t57 - t24 * t61;
t119 = g(2) * t123 - g(3) * (-t27 * t57 + t32 * t61) - g(1) * t8;
t112 = t54 * t62;
t111 = t55 * t62;
t110 = t57 * t62;
t107 = t61 * t62;
t86 = t117 * t101;
t50 = pkin(2) * t86;
t104 = -t32 * pkin(3) + t50;
t99 = pkin(5) * t57 + pkin(9);
t98 = t20 * t58 - t62 * t94;
t96 = t59 * t103;
t34 = pkin(2) * t96 + (-pkin(8) - qJ(3)) * t101;
t53 = pkin(2) * t117 + pkin(1);
t97 = -t60 * t34 + t63 * t53;
t93 = t124 * pkin(2);
t92 = t33 * pkin(9) + t104;
t91 = t25 * pkin(3) + t97;
t90 = pkin(4) * t62 + pkin(10) * t58;
t14 = t25 * t58 - t62 * t95;
t89 = g(1) * t98 + g(2) * t14;
t88 = g(1) * t21 - g(2) * t24;
t85 = -t63 * t34 - t60 * t53;
t52 = t61 * pkin(5) + pkin(4);
t64 = -pkin(11) - pkin(10);
t84 = t52 * t62 - t58 * t64;
t83 = t21 * pkin(3) + t93;
t77 = pkin(3) * t20 + t85;
t76 = -t24 * pkin(9) + t91;
t26 = t103 * t62 - t33 * t58;
t75 = g(1) * t14 - g(2) * t98 - g(3) * t26;
t74 = g(1) * t15 + g(2) * t11 + g(3) * t27;
t73 = -g(1) * t25 + g(2) * t20 - g(3) * t33;
t72 = g(1) * t24 + g(2) * t21 - g(3) * t32;
t71 = -t20 * pkin(9) + t83;
t70 = t21 * pkin(9) + t77;
t37 = -t63 * t59 - t60 * t87;
t69 = t37 * pkin(2);
t68 = t24 * pkin(3) + t69;
t66 = pkin(9) * t25 + t68;
t39 = -g(1) * t94 - g(2) * t95;
t38 = t117 * t63 - t60 * t96;
t36 = -t117 * t60 - t63 * t96;
t31 = -g(1) * t95 + g(2) * t94 - g(3) * t103;
t9 = t15 * t61 - t24 * t57;
t7 = t15 * t55 - t24 * t54;
t6 = -t15 * t54 - t24 * t55;
t4 = t72 * t58;
t2 = g(1) * t7 + g(2) * t120 - g(3) * (-t27 * t55 - t32 * t54);
t1 = -g(1) * t6 + g(2) * t121 - g(3) * (-t27 * t54 + t32 * t55);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t63, g(1) * t63 + g(2) * t60, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t38, g(1) * t124 - g(2) * t37, t39, -g(1) * (-t60 * pkin(1) + pkin(8) * t94) - g(2) * (t63 * pkin(1) + pkin(8) * t95) 0, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t25, t88, t39, -g(1) * t85 - g(2) * t97, 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t15, t89, -t88, -g(1) * t70 - g(2) * t76, 0, 0, 0, 0, 0, 0, g(1) * t122 - g(2) * t9, -g(1) * t123 - g(2) * t8, -t89, -g(1) * (-pkin(4) * t11 + pkin(10) * t98 + t70) - g(2) * (t15 * pkin(4) + t14 * pkin(10) + t76) 0, 0, 0, 0, 0, 0, g(1) * t120 - g(2) * t7, -g(1) * t121 - g(2) * t6, -t89, -g(1) * (-t11 * t52 + t21 * t99 - t64 * t98 + t77) - g(2) * (-t14 * t64 + t15 * t52 - t24 * t99 + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t37 - g(2) * t124 - g(3) * t86, g(3) * t101 * t59 + g(1) * t38 - g(2) * t36, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t73, 0, -g(1) * t69 - g(2) * t93 - g(3) * t50, 0, 0, 0, 0, 0, 0, -t72 * t62, t4, t73, -g(1) * t66 - g(2) * t71 - g(3) * t92, 0, 0, 0, 0, 0, 0, -g(1) * (t107 * t24 + t25 * t57) - g(2) * (t107 * t21 - t20 * t57) - g(3) * (-t107 * t32 + t33 * t57) -g(1) * (-t110 * t24 + t25 * t61) - g(2) * (-t110 * t21 - t20 * t61) - g(3) * (t110 * t32 + t33 * t61) -t4, -g(1) * (t24 * t90 + t66) - g(2) * (t21 * t90 + t71) - g(3) * (-t32 * t90 + t92) 0, 0, 0, 0, 0, 0, -g(1) * (t111 * t24 + t25 * t54) - g(2) * (t111 * t21 - t20 * t54) - g(3) * (-t111 * t32 + t33 * t54) -g(1) * (-t112 * t24 + t25 * t55) - g(2) * (-t112 * t21 - t20 * t55) - g(3) * (t112 * t32 + t33 * t55) -t4, -g(1) * (t84 * t24 + t25 * t99 + t68) - g(2) * (-t99 * t20 + t84 * t21 + t83) - g(3) * (-t84 * t32 + t99 * t33 + t104); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t61, -t75 * t57, -t74, -g(1) * (-t14 * pkin(4) + t15 * pkin(10)) - g(2) * (pkin(4) * t98 + t11 * pkin(10)) - g(3) * (t26 * pkin(4) + t27 * pkin(10)) 0, 0, 0, 0, 0, 0, t75 * t55, -t75 * t54, -t74, -g(1) * (-t14 * t52 - t15 * t64) - g(2) * (-t11 * t64 + t52 * t98) - g(3) * (t26 * t52 - t27 * t64); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, g(1) * t9 + g(2) * t122 - g(3) * (-t27 * t61 - t32 * t57) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t119 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
