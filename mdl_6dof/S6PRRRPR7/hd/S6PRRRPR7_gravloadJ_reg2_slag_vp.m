% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t113 = cos(pkin(12));
t114 = cos(pkin(7));
t110 = sin(pkin(12));
t121 = sin(qJ(2));
t123 = cos(qJ(2));
t115 = cos(pkin(6));
t93 = t115 * t113;
t69 = t110 * t121 - t123 * t93;
t111 = sin(pkin(7));
t112 = sin(pkin(6));
t89 = t112 * t111;
t126 = t113 * t89 + t69 * t114;
t91 = t115 * t110;
t70 = t113 * t121 + t123 * t91;
t88 = t112 * t110;
t125 = -t111 * t88 + t70 * t114;
t90 = t114 * t112;
t124 = t115 * t111 + t123 * t90;
t122 = cos(qJ(3));
t57 = pkin(13) + qJ(6);
t55 = sin(t57);
t63 = cos(qJ(4));
t120 = t55 * t63;
t56 = cos(t57);
t119 = t56 * t63;
t58 = sin(pkin(13));
t118 = t58 * t63;
t59 = cos(pkin(13));
t117 = t59 * t63;
t75 = t121 * t89;
t98 = t123 * t112;
t116 = pkin(2) * t98 + pkin(9) * t75;
t62 = sin(qJ(3));
t77 = t121 * t90;
t40 = t122 * t98 - t62 * t77;
t109 = t40 * pkin(3) + t116;
t108 = pkin(5) * t58 + pkin(10);
t44 = t110 * t123 + t121 * t93;
t15 = t126 * t122 + t44 * t62;
t13 = t15 * pkin(3);
t16 = t44 * t122 - t126 * t62;
t107 = t16 * pkin(10) - t13;
t45 = t113 * t123 - t121 * t91;
t17 = t125 * t122 + t45 * t62;
t14 = t17 * pkin(3);
t18 = t45 * t122 - t125 * t62;
t106 = t18 * pkin(10) - t14;
t97 = t112 * t121;
t32 = -t124 * t122 + t62 * t97;
t31 = t32 * pkin(3);
t33 = t122 * t97 + t124 * t62;
t105 = t33 * pkin(10) - t31;
t104 = t44 * t111;
t103 = t45 * t111;
t61 = sin(qJ(4));
t102 = t61 * t111;
t101 = t62 * t114;
t100 = t63 * t111;
t99 = t114 * t122;
t96 = -pkin(4) * t63 - qJ(5) * t61;
t54 = t59 * pkin(5) + pkin(4);
t60 = -pkin(11) - qJ(5);
t95 = -t54 * t63 + t60 * t61;
t39 = t122 * t77 + t62 * t98;
t94 = t39 * pkin(10) + t109;
t87 = -t69 * pkin(2) + pkin(9) * t104;
t86 = -t70 * pkin(2) + pkin(9) * t103;
t68 = t115 * t114 - t123 * t89;
t19 = t33 * t61 - t68 * t63;
t64 = t69 * t111 - t113 * t90;
t5 = t16 * t61 - t64 * t63;
t65 = t70 * t111 + t114 * t88;
t7 = t18 * t61 - t65 * t63;
t85 = g(1) * t7 + g(2) * t5 + g(3) * t19;
t20 = t33 * t63 + t68 * t61;
t6 = t16 * t63 + t64 * t61;
t8 = t18 * t63 + t65 * t61;
t84 = g(1) * t8 + g(2) * t6 + g(3) * t20;
t24 = -t44 * t101 - t69 * t122;
t83 = t24 * pkin(3) + t87;
t26 = -t45 * t101 - t70 * t122;
t82 = t26 * pkin(3) + t86;
t11 = -t45 * t100 + t26 * t61;
t27 = t40 * t61 - t63 * t75;
t9 = -t44 * t100 + t24 * t61;
t81 = g(1) * t11 + g(2) * t9 + g(3) * t27;
t80 = g(1) * t17 + g(2) * t15 + g(3) * t32;
t79 = g(1) * t18 + g(2) * t16 + g(3) * t33;
t23 = t44 * t99 - t69 * t62;
t25 = t45 * t99 - t70 * t62;
t78 = g(1) * t25 + g(2) * t23 + g(3) * t39;
t72 = t23 * pkin(10) + t83;
t71 = t25 * pkin(10) + t82;
t28 = t40 * t63 + t61 * t75;
t12 = t45 * t102 + t26 * t63;
t10 = t44 * t102 + t24 * t63;
t4 = t80 * t61;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t70 + g(2) * t69 - g(3) * t98, g(1) * t45 + g(2) * t44 + g(3) * t97, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t24 - g(3) * t40, t78, -g(1) * t103 - g(2) * t104 - g(3) * t75, -g(1) * t86 - g(2) * t87 - g(3) * t116, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10 - g(3) * t28, t81, -t78, -g(1) * t71 - g(2) * t72 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t59 + t25 * t58) - g(2) * (t10 * t59 + t23 * t58) - g(3) * (t28 * t59 + t39 * t58) -g(1) * (-t12 * t58 + t25 * t59) - g(2) * (-t10 * t58 + t23 * t59) - g(3) * (-t28 * t58 + t39 * t59) -t81, -g(1) * (t12 * pkin(4) + t11 * qJ(5) + t71) - g(2) * (t10 * pkin(4) + t9 * qJ(5) + t72) - g(3) * (t28 * pkin(4) + t27 * qJ(5) + t94) 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t56 + t25 * t55) - g(2) * (t10 * t56 + t23 * t55) - g(3) * (t28 * t56 + t39 * t55) -g(1) * (-t12 * t55 + t25 * t56) - g(2) * (-t10 * t55 + t23 * t56) - g(3) * (-t28 * t55 + t39 * t56) -t81, -g(1) * (t108 * t25 - t11 * t60 + t12 * t54 + t82) - g(2) * (t10 * t54 + t108 * t23 - t9 * t60 + t83) - g(3) * (t108 * t39 - t27 * t60 + t28 * t54 + t109); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t79, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t63, -t4, -t79, -g(1) * t106 - g(2) * t107 - g(3) * t105, 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t117 + t18 * t58) - g(2) * (-t15 * t117 + t16 * t58) - g(3) * (-t32 * t117 + t33 * t58) -g(1) * (t17 * t118 + t18 * t59) - g(2) * (t15 * t118 + t16 * t59) - g(3) * (t32 * t118 + t33 * t59) t4, -g(1) * (t96 * t17 + t106) - g(2) * (t96 * t15 + t107) - g(3) * (t96 * t32 + t105) 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t119 + t18 * t55) - g(2) * (-t15 * t119 + t16 * t55) - g(3) * (-t32 * t119 + t33 * t55) -g(1) * (t17 * t120 + t18 * t56) - g(2) * (t15 * t120 + t16 * t56) - g(3) * (t32 * t120 + t33 * t56) t4, -g(1) * (t108 * t18 + t95 * t17 - t14) - g(2) * (t108 * t16 + t95 * t15 - t13) - g(3) * (t108 * t33 + t95 * t32 - t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t84, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t59, -t85 * t58, -t84, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - g(3) * (-t19 * pkin(4) + t20 * qJ(5)) 0, 0, 0, 0, 0, 0, t85 * t56, -t85 * t55, -t84, -g(1) * (-t7 * t54 - t8 * t60) - g(2) * (-t5 * t54 - t6 * t60) - g(3) * (-t19 * t54 - t20 * t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t56 - t8 * t55) - g(2) * (t15 * t56 - t6 * t55) - g(3) * (-t20 * t55 + t32 * t56) -g(1) * (-t17 * t55 - t8 * t56) - g(2) * (-t15 * t55 - t6 * t56) - g(3) * (-t20 * t56 - t32 * t55) 0, 0;];
taug_reg  = t1;
