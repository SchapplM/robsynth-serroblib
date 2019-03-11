% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t122 = sin(qJ(1));
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t123 = cos(qJ(1));
t98 = cos(pkin(6));
t93 = t98 * t123;
t53 = t122 * t77 - t80 * t93;
t78 = cos(qJ(4));
t116 = t53 * t78;
t54 = t122 * t80 + t77 * t93;
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t73 = sin(pkin(6));
t96 = t73 * t123;
t28 = t54 * t79 - t76 * t96;
t75 = sin(qJ(4));
t121 = t28 * t75;
t128 = t116 - t121;
t74 = -qJ(5) - pkin(10);
t127 = -t74 * t76 + pkin(2);
t68 = pkin(4) * t78 + pkin(3);
t126 = -t68 * t79 - t127;
t117 = t53 * t75;
t125 = t28 * t78 + t117;
t72 = qJ(4) + pkin(11);
t69 = sin(t72);
t70 = cos(t72);
t4 = t28 * t69 - t53 * t70;
t5 = t28 * t70 + t53 * t69;
t124 = g(3) * t73;
t92 = t98 * t122;
t56 = t123 * t80 - t77 * t92;
t95 = t73 * t122;
t32 = t56 * t79 + t76 * t95;
t120 = t32 * t75;
t115 = t54 * t75;
t55 = t123 * t77 + t80 * t92;
t114 = t55 * t75;
t113 = t55 * t78;
t112 = t56 * t75;
t110 = t69 * t79;
t109 = t70 * t79;
t108 = t73 * t77;
t107 = t73 * t80;
t105 = t75 * t77;
t104 = t75 * t79;
t103 = t78 * t79;
t102 = t79 * t80;
t27 = t54 * t76 + t79 * t96;
t101 = -t27 * t68 - t28 * t74;
t31 = t56 * t76 - t79 * t95;
t100 = -t31 * t68 - t32 * t74;
t51 = t108 * t76 - t79 * t98;
t52 = t108 * t79 + t76 * t98;
t99 = -t51 * t68 - t52 * t74;
t97 = t78 * t107;
t94 = -g(1) * t27 + g(2) * t31;
t91 = -pkin(5) * t70 - qJ(6) * t69;
t90 = -t52 * t75 - t97;
t21 = t107 * t70 + t52 * t69;
t8 = t32 * t69 - t55 * t70;
t89 = g(1) * t8 + g(2) * t4 + g(3) * t21;
t88 = g(1) * t31 + g(2) * t27 + g(3) * t51;
t87 = g(1) * t32 + g(2) * t28 + g(3) * t52;
t86 = -g(1) * t55 - g(2) * t53 + g(3) * t107;
t85 = pkin(1) * t123 + pkin(2) * t56 + pkin(4) * t114 + pkin(8) * t95 + pkin(9) * t55 - t31 * t74 + t32 * t68;
t84 = pkin(9) * t108 + (pkin(4) * t105 + t102 * t68) * t73 + t127 * t107;
t83 = pkin(4) * t115 + t54 * pkin(9) + t126 * t53;
t82 = pkin(4) * t112 + t56 * pkin(9) + t126 * t55;
t81 = -pkin(1) * t122 - pkin(2) * t54 - pkin(4) * t117 + pkin(8) * t96 - pkin(9) * t53 + t27 * t74 - t28 * t68;
t44 = pkin(4) * t113;
t40 = pkin(4) * t116;
t34 = (t102 * t70 + t69 * t77) * t73;
t33 = (t102 * t69 - t70 * t77) * t73;
t22 = -t107 * t69 + t52 * t70;
t16 = -t109 * t55 + t56 * t69;
t15 = -t110 * t55 - t56 * t70;
t14 = -t109 * t53 + t54 * t69;
t13 = -t110 * t53 - t54 * t70;
t12 = t86 * t76;
t11 = t32 * t78 + t114;
t10 = t113 - t120;
t9 = t32 * t70 + t55 * t69;
t1 = [0, g(1) * t122 - g(2) * t123, g(1) * t123 + g(2) * t122, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t56, -g(1) * t53 + g(2) * t55, 0, 0, 0, 0, 0, g(1) * t28 - g(2) * t32, t94, 0, 0, 0, 0, 0, g(1) * t125 - g(2) * t11, g(1) * t128 - g(2) * t10, -t94, -g(1) * t81 - g(2) * t85, g(1) * t5 - g(2) * t9, -t94, g(1) * t4 - g(2) * t8, -g(1) * (-pkin(5) * t5 - qJ(6) * t4 + t81) - g(2) * (pkin(5) * t9 + qJ(6) * t8 + t85); 0, 0, 0, 0, 0, 0, 0, 0, -t86, g(1) * t56 + g(2) * t54 + g(3) * t108, 0, 0, 0, 0, 0, -t86 * t79, t12, 0, 0, 0, 0, 0, -g(1) * (-t103 * t55 + t112) - g(2) * (-t103 * t53 + t115) - (t102 * t78 + t105) * t124, -g(1) * (t104 * t55 + t56 * t78) - g(2) * (t104 * t53 + t54 * t78) - (-t102 * t75 + t77 * t78) * t124, -t12, -g(1) * t82 - g(2) * t83 - g(3) * t84, -g(1) * t16 - g(2) * t14 - g(3) * t34, -t12, -g(1) * t15 - g(2) * t13 - g(3) * t33, -g(1) * (pkin(5) * t16 + qJ(6) * t15 + t82) - g(2) * (pkin(5) * t14 + qJ(6) * t13 + t83) - g(3) * (pkin(5) * t34 + qJ(6) * t33 + t84); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t87, 0, 0, 0, 0, 0, t88 * t78, -t88 * t75, -t87, -g(1) * t100 - g(2) * t101 - g(3) * t99, t88 * t70, -t87, t88 * t69, -g(1) * (t31 * t91 + t100) - g(2) * (t27 * t91 + t101) - g(3) * (t51 * t91 + t99); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t128 - g(3) * t90, g(1) * t11 + g(2) * t125 - g(3) * (t107 * t75 - t52 * t78) 0, -g(1) * t44 - g(2) * t40 + (g(3) * t97 + t75 * t87) * pkin(4), t89, 0, -g(1) * t9 - g(2) * t5 - g(3) * t22, -g(1) * (-pkin(4) * t120 - pkin(5) * t8 + qJ(6) * t9 + t44) - g(2) * (-pkin(4) * t121 - pkin(5) * t4 + qJ(6) * t5 + t40) - g(3) * (pkin(4) * t90 - t21 * pkin(5) + t22 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, 0, 0, 0, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89;];
taug_reg  = t1;
