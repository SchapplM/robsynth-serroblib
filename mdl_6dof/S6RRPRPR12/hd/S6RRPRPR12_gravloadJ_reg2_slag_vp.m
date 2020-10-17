% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:18:56
% EndTime: 2019-05-06 16:18:59
% DurationCPUTime: 0.80s
% Computational Cost: add. (527->140), mult. (1059->206), div. (0->0), fcn. (1260->12), ass. (0->79)
t57 = sin(pkin(6));
t66 = cos(qJ(1));
t101 = t57 * t66;
t61 = sin(qJ(2));
t62 = sin(qJ(1));
t65 = cos(qJ(2));
t94 = cos(pkin(6));
t87 = t66 * t94;
t35 = t61 * t62 - t65 * t87;
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t75 = t101 * t60 + t35 * t64;
t103 = t57 * t62;
t88 = t62 * t94;
t37 = t61 * t66 + t65 * t88;
t14 = -t103 * t60 + t37 * t64;
t116 = pkin(4) * t60;
t36 = t61 * t87 + t62 * t65;
t58 = -qJ(5) - pkin(9);
t121 = t116 * t36 + t35 * t58;
t38 = -t61 * t88 + t65 * t66;
t120 = t116 * t38 + t37 * t58;
t119 = -g(1) * t38 - g(2) * t36;
t59 = sin(qJ(6));
t63 = cos(qJ(6));
t56 = qJ(4) + pkin(11);
t53 = sin(t56);
t54 = cos(t56);
t77 = t101 * t54 - t35 * t53;
t118 = t36 * t63 + t59 * t77;
t117 = -t36 * t59 + t63 * t77;
t113 = g(3) * t57;
t112 = t35 * t60;
t108 = t37 * t60;
t106 = t53 * t59;
t105 = t53 * t63;
t104 = t57 * t61;
t102 = t57 * t65;
t100 = t58 * t65;
t99 = t59 * t61;
t98 = t61 * t63;
t97 = t75 * pkin(4);
t96 = pkin(2) * t102 + qJ(3) * t104;
t95 = pkin(1) * t66 + pkin(8) * t103;
t91 = t104 * t116 + t96;
t90 = -pkin(1) * t62 + pkin(8) * t101;
t89 = -t101 * t64 + t112;
t31 = t35 * pkin(2);
t86 = qJ(3) * t36 - t31;
t33 = t37 * pkin(2);
t85 = qJ(3) * t38 - t33;
t76 = t101 * t53 + t35 * t54;
t8 = t103 * t53 - t37 * t54;
t84 = g(1) * t76 + g(2) * t8;
t83 = pkin(5) * t53 - pkin(10) * t54;
t82 = g(1) * t35 - g(2) * t37;
t7 = g(1) * t36 - g(2) * t38;
t81 = g(1) * t66 + g(2) * t62;
t80 = t14 * pkin(4);
t79 = pkin(2) * t38 + qJ(3) * t37 + t95;
t19 = -t102 * t54 - t53 * t94;
t74 = g(1) * t8 - g(2) * t76 - g(3) * t19;
t20 = -t102 * t53 + t54 * t94;
t9 = t103 * t54 + t37 * t53;
t73 = g(1) * t9 - g(2) * t77 + g(3) * t20;
t72 = -pkin(2) * t36 - qJ(3) * t35 + t90;
t71 = -t102 * t64 - t60 * t94;
t4 = -g(1) * t37 - g(2) * t35 + g(3) * t102;
t70 = g(3) * t104 - t119;
t52 = pkin(4) * t64 + pkin(3);
t69 = pkin(4) * t108 + t103 * t52 - t38 * t58 + t79;
t68 = t71 * pkin(4);
t67 = -pkin(4) * t112 + t101 * t52 + t36 * t58 + t72;
t39 = t81 * t57;
t15 = t103 * t64 + t108;
t3 = t38 * t59 + t63 * t9;
t2 = t38 * t63 - t59 * t9;
t1 = t70 * t54;
t5 = [0, 0, 0, 0, 0, 0, g(1) * t62 - g(2) * t66, t81, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t82, -t39, -g(1) * t90 - g(2) * t95, 0, 0, 0, 0, 0, 0, -t39, -t7, t82, -g(1) * t72 - g(2) * t79, 0, 0, 0, 0, 0, 0, g(1) * t89 - g(2) * t15, g(1) * t75 - g(2) * t14, t7, -g(1) * (pkin(3) * t101 - pkin(9) * t36 + t72) - g(2) * (pkin(3) * t103 + pkin(9) * t38 + t79) 0, 0, 0, 0, 0, 0, -g(1) * t77 - g(2) * t9, t84, t7, -g(1) * t67 - g(2) * t69, 0, 0, 0, 0, 0, 0, -g(1) * t117 - g(2) * t3, g(1) * t118 - g(2) * t2, -t84, -g(1) * (pkin(5) * t77 + pkin(10) * t76 + t67) - g(2) * (pkin(5) * t9 + pkin(10) * t8 + t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t70, -g(1) * t85 - g(2) * t86 - g(3) * t96, 0, 0, 0, 0, 0, 0, -t70 * t60, -t70 * t64, -t4, -g(1) * (-pkin(9) * t37 + t85) - g(2) * (-pkin(9) * t35 + t86) - g(3) * (pkin(9) * t102 + t96) 0, 0, 0, 0, 0, 0, -t70 * t53, -t1, -t4, -g(1) * (t85 + t120) - g(2) * (t86 + t121) - g(3) * (-t100 * t57 + t91) 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t38 - t37 * t59) - g(2) * (t105 * t36 - t35 * t59) - (t53 * t98 + t59 * t65) * t113, -g(1) * (-t106 * t38 - t37 * t63) - g(2) * (-t106 * t36 - t35 * t63) - (-t53 * t99 + t63 * t65) * t113, t1, -g(1) * (-t33 + t120) - g(2) * (-t31 + t121) - g(3) * t91 - (t61 * t83 - t100) * t113 + t119 * (qJ(3) + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t75 - g(3) * t71, g(1) * t15 + g(2) * t89 - g(3) * (t102 * t60 - t64 * t94) 0, 0, 0, 0, 0, 0, 0, 0, t74, t73, 0, -g(1) * t80 - g(2) * t97 - g(3) * t68, 0, 0, 0, 0, 0, 0, t74 * t63, -t74 * t59, -t73, -g(1) * (-pkin(5) * t8 + pkin(10) * t9 + t80) - g(2) * (pkin(5) * t76 - pkin(10) * t77 + t97) - g(3) * (t19 * pkin(5) + t20 * pkin(10) + t68); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t118 - g(3) * (-t20 * t59 + t57 * t98) g(1) * t3 - g(2) * t117 - g(3) * (-t20 * t63 - t57 * t99) 0, 0;];
taug_reg  = t5;
