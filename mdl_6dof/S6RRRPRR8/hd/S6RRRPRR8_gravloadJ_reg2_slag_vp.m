% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:24:43
% EndTime: 2019-05-07 12:24:49
% DurationCPUTime: 1.10s
% Computational Cost: add. (812->173), mult. (1399->271), div. (0->0), fcn. (1686->14), ass. (0->90)
t57 = sin(pkin(6));
t61 = sin(qJ(2));
t101 = t57 * t61;
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t91 = cos(pkin(6));
t124 = -t101 * t60 + t64 * t91;
t115 = cos(qJ(1));
t65 = cos(qJ(2));
t62 = sin(qJ(1));
t83 = t62 * t91;
t34 = t115 * t65 - t61 * t83;
t99 = t57 * t64;
t19 = -t34 * t60 + t62 * t99;
t77 = t91 * t115;
t32 = t61 * t77 + t62 * t65;
t55 = qJ(3) + pkin(12);
t50 = sin(t55);
t51 = cos(t55);
t87 = t57 * t115;
t14 = t32 * t51 - t50 * t87;
t31 = t61 * t62 - t65 * t77;
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t123 = t14 * t59 - t31 * t63;
t112 = t31 * t59;
t122 = t14 * t63 + t112;
t56 = qJ(5) + qJ(6);
t52 = sin(t56);
t53 = cos(t56);
t121 = t14 * t52 - t31 * t53;
t120 = t14 * t53 + t31 * t52;
t26 = t101 * t51 + t50 * t91;
t100 = t57 * t62;
t18 = t100 * t50 + t34 * t51;
t33 = t115 * t61 + t65 * t83;
t8 = -t18 * t59 + t33 * t63;
t95 = t63 * t65;
t119 = g(2) * t123 - g(3) * (-t26 * t59 - t57 * t95) - g(1) * t8;
t49 = pkin(3) * t64 + pkin(2);
t98 = t57 * t65;
t35 = t49 * t98;
t117 = g(3) * t35;
t116 = g(3) * t57;
t110 = t32 * t59;
t109 = t33 * t59;
t108 = t34 * t59;
t106 = t51 * t52;
t105 = t51 * t53;
t104 = t51 * t59;
t103 = t51 * t63;
t102 = t51 * t65;
t58 = -qJ(4) - pkin(9);
t97 = t58 * t61;
t96 = t59 * t65;
t94 = -t31 * t49 - t32 * t58;
t93 = -t33 * t49 - t34 * t58;
t92 = pkin(1) * t115 + pkin(8) * t100;
t89 = t60 * t100;
t86 = -pkin(1) * t62 + pkin(8) * t87;
t85 = -t32 * t50 - t51 * t87;
t43 = t60 * t87;
t84 = t32 * t64 - t43;
t81 = t19 * pkin(3);
t80 = pkin(3) * t89 - t33 * t58 + t34 * t49 + t92;
t79 = pkin(4) * t51 + pkin(10) * t50;
t17 = -t100 * t51 + t34 * t50;
t78 = g(1) * t85 + g(2) * t17;
t12 = g(1) * t31 - g(2) * t33;
t48 = pkin(5) * t63 + pkin(4);
t66 = -pkin(11) - pkin(10);
t76 = t48 * t51 - t50 * t66;
t75 = t124 * pkin(3);
t74 = g(1) * t115 + g(2) * t62;
t73 = pkin(3) * t43 + t31 * t58 - t32 * t49 + t86;
t25 = -t101 * t50 + t51 * t91;
t72 = g(1) * t17 - g(2) * t85 - g(3) * t25;
t71 = g(1) * t18 + g(2) * t14 + g(3) * t26;
t70 = t32 * t60 + t64 * t87;
t10 = -g(1) * t33 - g(2) * t31 + g(3) * t98;
t69 = g(1) * t34 + g(2) * t32 + g(3) * t101;
t68 = t70 * pkin(3);
t20 = t34 * t64 + t89;
t9 = t18 * t63 + t109;
t7 = t10 * t50;
t6 = t18 * t53 + t33 * t52;
t5 = -t18 * t52 + t33 * t53;
t2 = g(1) * t6 + g(2) * t120 - g(3) * (-t26 * t53 + t52 * t98);
t1 = -g(1) * t5 + g(2) * t121 - g(3) * (-t26 * t52 - t53 * t98);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t62 - g(2) * t115, t74, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t32 - g(2) * t34, -t12, -t74 * t57, -g(1) * t86 - g(2) * t92, 0, 0, 0, 0, 0, 0, g(1) * t84 - g(2) * t20, -g(1) * t70 - g(2) * t19, t12, -g(1) * (-pkin(2) * t32 - pkin(9) * t31 + t86) - g(2) * (pkin(2) * t34 + pkin(9) * t33 + t92) 0, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t78, t12, -g(1) * t73 - g(2) * t80, 0, 0, 0, 0, 0, 0, g(1) * t122 - g(2) * t9, -g(1) * t123 - g(2) * t8, -t78, -g(1) * (-pkin(4) * t14 + pkin(10) * t85 + t73) - g(2) * (pkin(4) * t18 + pkin(10) * t17 + t80) 0, 0, 0, 0, 0, 0, g(1) * t120 - g(2) * t6, -g(1) * t121 - g(2) * t5, -t78, -g(1) * (-pkin(5) * t112 - t14 * t48 - t66 * t85 + t73) - g(2) * (pkin(5) * t109 - t17 * t66 + t18 * t48 + t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t69, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t64, t10 * t60, -t69, -g(1) * (-pkin(2) * t33 + pkin(9) * t34) - g(2) * (-pkin(2) * t31 + pkin(9) * t32) - (pkin(2) * t65 + pkin(9) * t61) * t116, 0, 0, 0, 0, 0, 0, -t10 * t51, t7, -t69, -g(1) * t93 - g(2) * t94 - g(3) * (-t57 * t97 + t35) 0, 0, 0, 0, 0, 0, -g(1) * (-t103 * t33 + t108) - g(2) * (-t103 * t31 + t110) - (t51 * t95 + t59 * t61) * t116, -g(1) * (t104 * t33 + t34 * t63) - g(2) * (t104 * t31 + t32 * t63) - (-t51 * t96 + t61 * t63) * t116, -t7, -g(1) * (-t33 * t79 + t93) - g(2) * (-t31 * t79 + t94) - t117 - (t65 * t79 - t97) * t116, 0, 0, 0, 0, 0, 0, -g(1) * (-t105 * t33 + t34 * t52) - g(2) * (-t105 * t31 + t32 * t52) - (t102 * t53 + t52 * t61) * t116, -g(1) * (t106 * t33 + t34 * t53) - g(2) * (t106 * t31 + t32 * t53) - (-t102 * t52 + t53 * t61) * t116, -t7, -g(1) * (pkin(5) * t108 - t33 * t76 + t93) - g(2) * (pkin(5) * t110 - t31 * t76 + t94) - t117 - (t76 * t65 + (pkin(5) * t59 - t58) * t61) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t19 + g(2) * t70 - g(3) * t124, g(1) * t20 + g(2) * t84 - g(3) * (-t60 * t91 - t61 * t99) 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, 0, -g(1) * t81 + g(2) * t68 - g(3) * t75, 0, 0, 0, 0, 0, 0, t72 * t63, -t72 * t59, -t71, -g(1) * (-pkin(4) * t17 + pkin(10) * t18 + t81) - g(2) * (pkin(4) * t85 + t14 * pkin(10) - t68) - g(3) * (pkin(4) * t25 + pkin(10) * t26 + t75) 0, 0, 0, 0, 0, 0, t72 * t53, -t72 * t52, -t71, -g(1) * (-t17 * t48 - t18 * t66 + t81) - g(2) * (-t14 * t66 + t48 * t85 - t68) - g(3) * (t25 * t48 - t26 * t66 + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, g(1) * t9 + g(2) * t122 - g(3) * (-t26 * t63 + t57 * t96) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t119 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
