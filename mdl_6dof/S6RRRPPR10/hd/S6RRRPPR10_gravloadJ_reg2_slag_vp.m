% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:05:34
% EndTime: 2019-05-07 07:05:38
% DurationCPUTime: 1.02s
% Computational Cost: add. (606->169), mult. (1445->235), div. (0->0), fcn. (1760->12), ass. (0->88)
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t118 = -pkin(3) * t60 - qJ(4) * t57;
t106 = cos(qJ(2));
t107 = cos(qJ(1));
t59 = sin(qJ(1));
t58 = sin(qJ(2));
t92 = cos(pkin(6));
t80 = t58 * t92;
t33 = t59 * t106 + t107 * t80;
t54 = sin(pkin(6));
t88 = t54 * t107;
t15 = t33 * t57 + t60 * t88;
t70 = t92 * t106;
t32 = -t107 * t70 + t59 * t58;
t52 = pkin(11) + qJ(6);
t49 = sin(t52);
t50 = cos(t52);
t117 = t15 * t49 + t32 * t50;
t116 = t15 * t50 - t32 * t49;
t115 = t118 * t32;
t34 = t107 * t58 + t59 * t70;
t114 = t118 * t34;
t113 = pkin(4) + pkin(9);
t111 = g(3) * t54;
t110 = t32 * pkin(9);
t109 = t34 * pkin(9);
t55 = cos(pkin(11));
t48 = t55 * pkin(5) + pkin(4);
t108 = pkin(9) + t48;
t103 = t49 * t57;
t102 = t50 * t57;
t53 = sin(pkin(11));
t101 = t53 * t57;
t100 = t54 * t58;
t99 = t54 * t59;
t98 = t54 * t60;
t97 = t55 * t57;
t87 = t54 * t106;
t96 = pkin(2) * t87 + pkin(9) * t100;
t95 = t107 * pkin(1) + pkin(8) * t99;
t93 = qJ(5) * t60;
t26 = t32 * pkin(2);
t91 = -t26 + t115;
t28 = t34 * pkin(2);
t90 = -t28 + t114;
t35 = t107 * t106 - t59 * t80;
t89 = t35 * pkin(2) + t95;
t86 = t57 * t106;
t85 = t60 * t106;
t84 = -t59 * pkin(1) + pkin(8) * t88;
t83 = t33 * pkin(9) - t26;
t82 = t35 * pkin(9) - t28;
t81 = pkin(5) * t53 + qJ(4);
t16 = t33 * t60 - t57 * t88;
t11 = t15 * pkin(3);
t79 = t16 * qJ(4) - t11;
t19 = t35 * t57 - t59 * t98;
t13 = t19 * pkin(3);
t20 = t35 * t60 + t57 * t99;
t78 = t20 * qJ(4) - t13;
t30 = t57 * t100 - t92 * t60;
t25 = t30 * pkin(3);
t31 = t92 * t57 + t58 * t98;
t77 = t31 * qJ(4) - t25;
t76 = t20 * pkin(3) + t89;
t75 = t96 + (pkin(3) * t85 + qJ(4) * t86) * t54;
t74 = t53 * t86;
t73 = -t33 * pkin(2) + t84;
t72 = -g(1) * t15 + g(2) * t19;
t71 = -g(1) * t16 + g(2) * t20;
t10 = g(1) * t32 - g(2) * t34;
t69 = g(3) * t75;
t68 = -pkin(3) * t16 + t73;
t56 = -pkin(10) - qJ(5);
t67 = -pkin(5) * t101 + t56 * t60;
t66 = g(1) * t107 + g(2) * t59;
t65 = t19 * qJ(4) + t76;
t2 = g(1) * t19 + g(2) * t15 + g(3) * t30;
t64 = g(1) * t20 + g(2) * t16 + g(3) * t31;
t63 = -qJ(4) * t15 + t68;
t62 = g(1) * t35 + g(2) * t33 + g(3) * t100;
t61 = -g(1) * t34 - g(2) * t32 + g(3) * t87;
t8 = t61 * t60;
t7 = t61 * t57;
t6 = t19 * t49 + t34 * t50;
t5 = t19 * t50 - t34 * t49;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * t107, t66, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t33 - g(2) * t35, -t10, -t66 * t54, -g(1) * t84 - g(2) * t95, 0, 0, 0, 0, 0, 0, -t71, t72, t10, -g(1) * (t73 - t110) - g(2) * (t89 + t109) 0, 0, 0, 0, 0, 0, t10, t71, -t72, -g(1) * (t63 - t110) - g(2) * (t65 + t109) 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t53 - t32 * t55) - g(2) * (t19 * t53 + t34 * t55) -g(1) * (-t15 * t55 + t32 * t53) - g(2) * (t19 * t55 - t34 * t53) -t71, -g(1) * (-qJ(5) * t16 - t113 * t32 + t63) - g(2) * (t20 * qJ(5) + t113 * t34 + t65) 0, 0, 0, 0, 0, 0, g(1) * t117 - g(2) * t6, g(1) * t116 - g(2) * t5, -t71, -g(1) * (-t108 * t32 - t15 * t81 + t16 * t56 + t68) - g(2) * (t108 * t34 + t81 * t19 - t20 * t56 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t62, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, -t62, -g(1) * t82 - g(2) * t83 - g(3) * t96, 0, 0, 0, 0, 0, 0, -t62, t8, -t7, -g(1) * (t82 + t114) - g(2) * (t83 + t115) - t69, 0, 0, 0, 0, 0, 0, -g(1) * (-t34 * t101 + t35 * t55) - g(2) * (-t32 * t101 + t33 * t55) - (t55 * t58 + t74) * t111, -g(1) * (-t34 * t97 - t35 * t53) - g(2) * (-t32 * t97 - t33 * t53) - (-t53 * t58 + t55 * t86) * t111, -t8, -g(1) * (t113 * t35 - t34 * t93 + t90) - g(2) * (t113 * t33 - t32 * t93 + t91) - g(3) * ((pkin(4) * t58 + qJ(5) * t85) * t54 + t75) 0, 0, 0, 0, 0, 0, -g(1) * (-t34 * t103 + t35 * t50) - g(2) * (-t32 * t103 + t33 * t50) - (t49 * t86 + t50 * t58) * t111, -g(1) * (-t34 * t102 - t35 * t49) - g(2) * (-t32 * t102 - t33 * t49) - (-t49 * t58 + t50 * t86) * t111, -t8, -g(1) * (t108 * t35 + t67 * t34 + t90) - g(2) * (t108 * t33 + t67 * t32 + t91) - t69 - (pkin(5) * t74 + t48 * t58 - t56 * t85) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t64, -g(1) * t78 - g(2) * t79 - g(3) * t77, 0, 0, 0, 0, 0, 0, -t64 * t53, -t64 * t55, t2, -g(1) * (-t19 * qJ(5) + t78) - g(2) * (-t15 * qJ(5) + t79) - g(3) * (-t30 * qJ(5) + t77) 0, 0, 0, 0, 0, 0, -t64 * t49, -t64 * t50, t2, -g(1) * (t19 * t56 + t81 * t20 - t13) - g(2) * (t15 * t56 + t81 * t16 - t11) - g(3) * (t30 * t56 + t81 * t31 - t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t116 - g(3) * (t30 * t50 + t49 * t87) g(1) * t6 + g(2) * t117 - g(3) * (-t30 * t49 + t50 * t87) 0, 0;];
taug_reg  = t1;
