% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t123 = pkin(8) + qJ(4);
t113 = sin(t123) / 0.2e1;
t124 = pkin(8) - qJ(4);
t118 = sin(t124);
t103 = t113 + t118 / 0.2e1;
t115 = cos(t123) / 0.2e1;
t120 = cos(t124);
t104 = t120 / 0.2e1 + t115;
t78 = pkin(7) + pkin(14);
t70 = sin(t78) / 0.2e1;
t79 = pkin(7) - pkin(14);
t76 = sin(t79);
t122 = t70 + t76 / 0.2e1;
t83 = sin(pkin(6));
t117 = t83 * t122;
t125 = pkin(6) + qJ(2);
t116 = cos(t125) / 0.2e1;
t126 = pkin(6) - qJ(2);
t121 = cos(t126);
t105 = t121 / 0.2e1 + t116;
t91 = sin(qJ(2));
t92 = sin(qJ(1));
t97 = cos(qJ(1));
t55 = -t97 * t105 + t92 * t91;
t114 = sin(t125) / 0.2e1;
t119 = sin(t126);
t67 = t114 - t119 / 0.2e1;
t96 = cos(qJ(2));
t56 = t97 * t67 + t92 * t96;
t71 = cos(t79) / 0.2e1;
t77 = cos(t78);
t63 = t71 + t77 / 0.2e1;
t80 = sin(pkin(14));
t39 = t97 * t117 + t55 * t63 + t56 * t80;
t130 = t83 * t97;
t62 = t70 - t76 / 0.2e1;
t64 = t71 - t77 / 0.2e1;
t84 = cos(pkin(14));
t40 = t64 * t130 + t55 * t62 - t56 * t84;
t82 = sin(pkin(7));
t86 = cos(pkin(7));
t53 = t86 * t130 - t55 * t82;
t90 = sin(qJ(4));
t15 = -t53 * t103 - t39 * t104 + t40 * t90;
t65 = t113 - t118 / 0.2e1;
t68 = t115 - t120 / 0.2e1;
t95 = cos(qJ(4));
t16 = t39 * t65 + t40 * t95 - t53 * t68;
t81 = sin(pkin(8));
t85 = cos(pkin(8));
t31 = -t39 * t81 + t53 * t85;
t89 = sin(qJ(5));
t94 = cos(qJ(5));
t5 = t16 * t94 + t31 * t89;
t88 = sin(qJ(6));
t93 = cos(qJ(6));
t146 = -t15 * t93 + t5 * t88;
t145 = t15 * t88 + t5 * t93;
t142 = t16 * t89 - t31 * t94;
t131 = t83 * t92;
t58 = -t92 * t105 - t97 * t91;
t111 = t86 * t131 - t58 * t82;
t139 = -g(1) * t53 - g(2) * t111;
t133 = t82 * t68;
t132 = t82 * t85;
t129 = t88 * t94;
t128 = t93 * t94;
t127 = qJ(3) * t82;
t66 = t114 + t119 / 0.2e1;
t87 = cos(pkin(6));
t112 = t66 * t82 - t87 * t86;
t60 = t92 * t67 - t97 * t96;
t69 = t116 - t121 / 0.2e1;
t106 = t87 * t122 + t66 * t63 + t69 * t80;
t48 = t66 * t62 + t87 * t64 - t69 * t84;
t100 = t106 * t65 + t112 * t68 + t48 * t95;
t37 = -t106 * t81 - t112 * t85;
t102 = -t92 * t117 - t58 * t63 - t60 * t80;
t32 = t102 * t81 + t111 * t85;
t41 = t64 * t131 + t58 * t62 - t60 * t84;
t98 = -t102 * t65 - t111 * t68 + t41 * t95;
t6 = t32 * t94 - t89 * t98;
t110 = g(1) * t6 + g(2) * t142 + g(3) * (-t100 * t89 + t37 * t94);
t17 = t102 * t104 - t111 * t103 + t41 * t90;
t25 = t112 * t103 - t106 * t104 + t48 * t90;
t109 = g(1) * t17 - g(2) * t15 + g(3) * t25;
t108 = g(1) * t60 - g(2) * t56 + g(3) * t69;
t101 = t82 * t103;
t52 = t69 * t62 + t66 * t84;
t51 = t69 * t63 - t66 * t80;
t47 = t58 * t84 + t60 * t62;
t46 = -t58 * t80 + t60 * t63;
t45 = -t55 * t84 - t56 * t62;
t44 = t55 * t80 - t56 * t63;
t43 = -t69 * t132 - t51 * t81;
t34 = -t60 * t132 - t46 * t81;
t33 = t132 * t56 - t44 * t81;
t29 = t69 * t133 + t51 * t65 + t52 * t95;
t28 = t69 * t101 - t51 * t104 + t52 * t90;
t24 = t60 * t133 + t46 * t65 + t47 * t95;
t23 = t60 * t101 - t46 * t104 + t47 * t90;
t22 = -t133 * t56 + t44 * t65 + t45 * t95;
t21 = -t101 * t56 - t44 * t104 + t45 * t90;
t20 = t29 * t94 + t43 * t89;
t11 = t100 * t94 + t37 * t89;
t9 = t24 * t94 + t34 * t89;
t8 = t22 * t94 + t33 * t89;
t7 = t32 * t89 + t94 * t98;
t2 = t17 * t88 + t7 * t93;
t1 = t17 * t93 - t7 * t88;
t3 = [0, g(1) * t92 - g(2) * t97, g(1) * t97 + g(2) * t92, 0, 0, 0, 0, 0, g(1) * t56 + g(2) * t60, -g(1) * t55 - g(2) * t58, -g(1) * t40 - g(2) * t41, -g(1) * t39 + g(2) * t102, t139, -g(1) * (-t92 * pkin(1) - t56 * pkin(2) + pkin(10) * t130) - g(2) * (t97 * pkin(1) - pkin(2) * t60 + pkin(10) * t131) + t139 * qJ(3), 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t98, g(1) * t15 + g(2) * t17, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, g(1) * t142 - g(2) * t6, 0, 0, 0, 0, 0, -g(1) * t145 - g(2) * t2, g(1) * t146 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t58 + g(2) * t55 - g(3) * t66, -t108, -g(1) * t47 - g(2) * t45 - g(3) * t52, -g(1) * t46 - g(2) * t44 - g(3) * t51, t108 * t82, -g(1) * (t58 * pkin(2) - t60 * t127) - g(2) * (-t55 * pkin(2) + t127 * t56) - g(3) * (t66 * pkin(2) - t69 * t127) 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t22 - g(3) * t29, g(1) * t23 + g(2) * t21 + g(3) * t28, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t20, -g(1) * (-t24 * t89 + t34 * t94) - g(2) * (-t22 * t89 + t33 * t94) - g(3) * (-t29 * t89 + t43 * t94) 0, 0, 0, 0, 0, -g(1) * (t23 * t88 + t9 * t93) - g(2) * (t21 * t88 + t8 * t93) - g(3) * (t20 * t93 + t28 * t88) -g(1) * (t23 * t93 - t9 * t88) - g(2) * (t21 * t93 - t8 * t88) - g(3) * (-t20 * t88 + t28 * t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t111 + g(2) * t53 + g(3) * t112, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, g(1) * t98 - g(2) * t16 + g(3) * t100, 0, 0, 0, 0, 0, t109 * t94, -t109 * t89, 0, 0, 0, 0, 0, -g(1) * (-t17 * t128 + t88 * t98) - g(2) * (t128 * t15 - t16 * t88) - g(3) * (t100 * t88 - t25 * t128) -g(1) * (t17 * t129 + t93 * t98) - g(2) * (-t129 * t15 - t16 * t93) - g(3) * (t100 * t93 + t25 * t129); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, g(1) * t7 - g(2) * t5 + g(3) * t11, 0, 0, 0, 0, 0, -t110 * t93, t110 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t146 - g(3) * (-t11 * t88 + t25 * t93) g(1) * t2 - g(2) * t145 - g(3) * (-t11 * t93 - t25 * t88);];
taug_reg  = t3;
