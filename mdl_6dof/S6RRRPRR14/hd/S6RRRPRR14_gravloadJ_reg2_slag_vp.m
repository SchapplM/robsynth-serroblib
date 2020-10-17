% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 16:43:13
% EndTime: 2019-05-07 16:43:17
% DurationCPUTime: 0.96s
% Computational Cost: add. (664->173), mult. (1570->248), div. (0->0), fcn. (1919->12), ass. (0->96)
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t63 = cos(qJ(2));
t114 = cos(qJ(1));
t93 = cos(pkin(6));
t75 = t93 * t114;
t36 = t59 * t75 + t60 * t63;
t58 = sin(qJ(3));
t62 = cos(qJ(3));
t56 = sin(pkin(6));
t89 = t56 * t114;
t18 = t36 * t58 + t62 * t89;
t35 = t60 * t59 - t63 * t75;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t127 = t18 * t57 + t35 * t61;
t126 = t18 * t61 - t35 * t57;
t55 = qJ(5) + qJ(6);
t52 = sin(t55);
t53 = cos(t55);
t125 = t18 * t52 + t35 * t53;
t124 = t18 * t53 - t35 * t52;
t109 = t35 * t62;
t94 = qJ(4) * t58;
t123 = -pkin(3) * t109 - t35 * t94;
t84 = t60 * t93;
t37 = t114 * t59 + t63 * t84;
t108 = t37 * t62;
t122 = -pkin(3) * t108 - t37 * t94;
t102 = t56 * t63;
t105 = t56 * t59;
t33 = t58 * t105 - t93 * t62;
t103 = t56 * t62;
t38 = t114 * t63 - t59 * t84;
t22 = -t60 * t103 + t38 * t58;
t8 = t22 * t61 - t37 * t57;
t121 = -g(2) * t126 - g(3) * (t57 * t102 + t33 * t61) - g(1) * t8;
t119 = pkin(4) + pkin(9);
t118 = g(3) * t56;
t117 = t35 * pkin(9);
t116 = t37 * pkin(9);
t51 = t61 * pkin(5) + pkin(4);
t115 = pkin(9) + t51;
t107 = t52 * t58;
t106 = t53 * t58;
t104 = t56 * t60;
t101 = t57 * t58;
t100 = t58 * t61;
t99 = t58 * t63;
t98 = t61 * t63;
t97 = t62 * t63;
t96 = pkin(2) * t102 + pkin(9) * t105;
t95 = t114 * pkin(1) + pkin(8) * t104;
t29 = t35 * pkin(2);
t92 = -t29 + t123;
t31 = t37 * pkin(2);
t91 = -t31 + t122;
t90 = t38 * pkin(2) + t95;
t88 = -t60 * pkin(1) + pkin(8) * t89;
t87 = t36 * pkin(9) - t29;
t86 = t38 * pkin(9) - t31;
t85 = pkin(5) * t57 + qJ(4);
t19 = t36 * t62 - t58 * t89;
t14 = t18 * pkin(3);
t83 = t19 * qJ(4) - t14;
t16 = t22 * pkin(3);
t23 = t58 * t104 + t38 * t62;
t82 = t23 * qJ(4) - t16;
t28 = t33 * pkin(3);
t34 = t59 * t103 + t93 * t58;
t81 = t34 * qJ(4) - t28;
t80 = t23 * pkin(3) + t90;
t79 = t56 * pkin(3) * t97 + t94 * t102 + t96;
t78 = -t36 * pkin(2) + t88;
t77 = -g(1) * t18 + g(2) * t22;
t76 = -g(1) * t19 + g(2) * t23;
t13 = g(1) * t35 - g(2) * t37;
t74 = g(3) * t79;
t73 = -pkin(3) * t19 + t78;
t64 = -pkin(11) - pkin(10);
t72 = pkin(5) * t101 - t62 * t64;
t71 = g(1) * t114 + g(2) * t60;
t70 = t22 * qJ(4) + t80;
t4 = g(1) * t22 + g(2) * t18 + g(3) * t33;
t69 = g(1) * t23 + g(2) * t19 + g(3) * t34;
t68 = -qJ(4) * t18 + t73;
t67 = -g(1) * t37 - g(2) * t35 + g(3) * t102;
t66 = g(1) * t38 + g(2) * t36 + g(3) * t105;
t11 = t67 * t62;
t10 = t67 * t58;
t9 = t22 * t57 + t37 * t61;
t7 = t22 * t52 + t37 * t53;
t6 = t22 * t53 - t37 * t52;
t2 = g(1) * t7 + g(2) * t125 - g(3) * (t53 * t102 - t33 * t52);
t1 = -g(1) * t6 - g(2) * t124 - g(3) * (t52 * t102 + t33 * t53);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t114, t71, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t36 - g(2) * t38, -t13, -t71 * t56, -g(1) * t88 - g(2) * t95, 0, 0, 0, 0, 0, 0, -t76, t77, t13, -g(1) * (t78 - t117) - g(2) * (t90 + t116) 0, 0, 0, 0, 0, 0, t13, t76, -t77, -g(1) * (t68 - t117) - g(2) * (t70 + t116) 0, 0, 0, 0, 0, 0, g(1) * t127 - g(2) * t9, g(1) * t126 - g(2) * t8, -t76, -g(1) * (-pkin(10) * t19 - t119 * t35 + t68) - g(2) * (t23 * pkin(10) + t119 * t37 + t70) 0, 0, 0, 0, 0, 0, g(1) * t125 - g(2) * t7, g(1) * t124 - g(2) * t6, -t76, -g(1) * (-t115 * t35 - t18 * t85 + t19 * t64 + t73) - g(2) * (t115 * t37 + t85 * t22 - t23 * t64 + t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t66, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, -t66, -g(1) * t86 - g(2) * t87 - g(3) * t96, 0, 0, 0, 0, 0, 0, -t66, t11, -t10, -g(1) * (t86 + t122) - g(2) * (t87 + t123) - t74, 0, 0, 0, 0, 0, 0, -g(1) * (-t37 * t101 + t38 * t61) - g(2) * (-t35 * t101 + t36 * t61) - (t57 * t99 + t59 * t61) * t118, -g(1) * (-t37 * t100 - t38 * t57) - g(2) * (-t35 * t100 - t36 * t57) - (-t57 * t59 + t58 * t98) * t118, -t11, -g(1) * (-pkin(10) * t108 + t119 * t38 + t91) - g(2) * (-pkin(10) * t109 + t119 * t36 + t92) - g(3) * ((pkin(4) * t59 + pkin(10) * t97) * t56 + t79) 0, 0, 0, 0, 0, 0, -g(1) * (-t37 * t107 + t38 * t53) - g(2) * (-t35 * t107 + t36 * t53) - (t52 * t99 + t53 * t59) * t118, -g(1) * (-t37 * t106 - t38 * t52) - g(2) * (-t35 * t106 - t36 * t52) - (-t52 * t59 + t53 * t99) * t118, -t11, -g(1) * (t115 * t38 - t72 * t37 + t91) - g(2) * (t115 * t36 - t72 * t35 + t92) - t74 - (t51 * t59 + t72 * t63) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t69, -g(1) * t82 - g(2) * t83 - g(3) * t81, 0, 0, 0, 0, 0, 0, -t69 * t57, -t69 * t61, t4, -g(1) * (-t22 * pkin(10) + t82) - g(2) * (-t18 * pkin(10) + t83) - g(3) * (-t33 * pkin(10) + t81) 0, 0, 0, 0, 0, 0, -t69 * t52, -t69 * t53, t4, -g(1) * (t22 * t64 + t85 * t23 - t16) - g(2) * (t18 * t64 + t85 * t19 - t14) - g(3) * (t33 * t64 + t85 * t34 - t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, g(1) * t9 + g(2) * t127 - g(3) * (-t33 * t57 + t56 * t98) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t121 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
