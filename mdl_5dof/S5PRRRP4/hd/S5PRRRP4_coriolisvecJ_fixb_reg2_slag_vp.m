% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:20
% EndTime: 2019-12-05 16:46:24
% DurationCPUTime: 0.84s
% Computational Cost: add. (1109->157), mult. (2335->192), div. (0->0), fcn. (1504->6), ass. (0->112)
t69 = cos(qJ(3));
t108 = qJD(3) * t69;
t101 = pkin(2) * t108;
t70 = cos(qJ(2));
t103 = t70 * qJD(1);
t67 = sin(qJ(2));
t110 = qJD(1) * t67;
t66 = sin(qJ(3));
t55 = t66 * t110;
t41 = t69 * t103 - t55;
t133 = t101 - t41;
t138 = 2 * qJD(4);
t109 = qJD(3) * t66;
t54 = qJD(2) * pkin(2) + t103;
t45 = t66 * t70 + t69 * t67;
t77 = t45 * qJD(2);
t73 = (t67 * t108 + t77) * qJD(1);
t17 = t54 * t109 + t73;
t65 = sin(qJ(4));
t63 = t65 ^ 2;
t68 = cos(qJ(4));
t64 = t68 ^ 2;
t111 = t63 + t64;
t62 = qJD(2) + qJD(3);
t113 = t62 * t55;
t91 = qJD(2) * t103;
t16 = t69 * (qJD(3) * t54 + t91) - t113;
t139 = t111 * t16;
t93 = t111 * t62;
t34 = t69 * t110 + t66 * t54;
t28 = t62 * pkin(7) + t34;
t121 = t65 * t28;
t13 = t68 * t16;
t4 = t13 + (qJD(5) - t121) * qJD(4);
t106 = qJD(4) * t68;
t11 = t65 * t16;
t7 = t28 * t106 + t11;
t137 = t4 * t68 + t7 * t65;
t44 = t66 * t67 - t69 * t70;
t136 = t44 * t62;
t134 = t111 * t28;
t112 = -t63 + t64;
t132 = t112 * t62 * t138;
t71 = qJD(4) ^ 2;
t131 = pkin(7) * t71;
t130 = t62 * pkin(3);
t129 = t69 * pkin(2);
t128 = t17 * t44;
t22 = t45 * qJD(3) + t77;
t127 = t22 * t62;
t126 = t34 * t62;
t50 = -t68 * pkin(4) - t65 * qJ(5) - pkin(3);
t125 = t50 * t62;
t56 = t66 * pkin(2) + pkin(7);
t124 = t56 * t71;
t123 = t62 * t65;
t122 = t62 * t68;
t120 = t68 * t28;
t33 = t69 * t54 - t55;
t27 = -t33 - t130;
t118 = t27 * t106 + t17 * t65;
t107 = qJD(4) * t65;
t117 = t33 * t107 + t34 * t122;
t100 = pkin(2) * t109;
t102 = qJD(4) * qJ(5);
t104 = t65 * qJD(5);
t38 = pkin(4) * t107 - t68 * t102 - t104;
t30 = t38 + t100;
t40 = t45 * qJD(1);
t116 = t30 - t40;
t115 = t41 * t107 + t40 * t122;
t114 = t34 - t38;
t20 = t102 + t120;
t105 = t20 * qJD(4);
t97 = t62 * t107;
t85 = pkin(4) * t65 - qJ(5) * t68;
t8 = (t85 * qJD(4) - t104) * t62 + t17;
t96 = -t8 - t131;
t95 = -t8 - t124;
t94 = t136 * t111;
t89 = -qJD(4) * pkin(4) + qJD(5);
t88 = t68 * t97;
t87 = t33 * t93;
t86 = -t40 + t100;
t19 = t89 + t121;
t84 = t19 * t65 + t20 * t68;
t83 = t45 * t71 + t127;
t82 = (-pkin(2) * t62 - t54) * qJD(3);
t81 = t136 * t138;
t42 = t50 - t129;
t80 = t42 * t62 - t101;
t57 = -pkin(3) - t129;
t79 = t57 * t62 - t101;
t78 = -t65 * t105 + t19 * t106 + t137;
t76 = t133 * t93;
t74 = (t19 * t68 - t20 * t65) * qJD(4) + t137;
t72 = qJD(2) ^ 2;
t61 = t62 ^ 2;
t60 = t71 * t68;
t59 = t71 * t65;
t51 = t65 * t61 * t68;
t49 = -0.2e1 * t88;
t48 = 0.2e1 * t88;
t43 = t112 * t61;
t39 = t85 * t62;
t24 = t27 * t107;
t15 = -t33 + t125;
t9 = t15 * t107;
t6 = t62 * t94;
t2 = t65 * t81 - t83 * t68;
t1 = t83 * t65 + t68 * t81;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t67, -t72 * t70, 0, 0, 0, 0, 0, 0, 0, 0, -t127, t136 * t62, 0, -t136 * t34 + t16 * t45 - t33 * t22 + t128, 0, 0, 0, 0, 0, 0, t2, t1, -t6, t139 * t45 + t27 * t22 - t28 * t94 + t128, 0, 0, 0, 0, 0, 0, t2, -t6, -t1, -t136 * t84 + t15 * t22 + t8 * t44 + t74 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t62 + t66 * t82 - t73, t41 * t62 + (t82 - t91) * t69 + t113, 0, t33 * t40 - t34 * t41 + (t16 * t66 - t17 * t69 + (-t33 * t66 + t34 * t69) * qJD(3)) * pkin(2), t48, t132, t60, t49, -t59, 0, t24 + t79 * t107 + (-t62 * t100 - t124 - t17) * t68 + t115, (t86 * t62 + t124) * t65 + (t41 + t79) * t106 + t118, t76 + t139, t133 * t134 + t139 * t56 + t17 * t57 + t86 * t27, t48, t60, -t132, 0, t59, t49, t9 + t80 * t107 + (-t30 * t62 + t95) * t68 + t115, t76 + t78, (-t116 * t62 + t95) * t65 + (-t15 - t41 - t80) * t106, t116 * t15 + t133 * t84 + t8 * t42 + t74 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126 - t17, t33 * t62 - t16, 0, 0, t48, t132, t60, t49, -t59, 0, -pkin(3) * t97 + t24 + (-t17 - t131) * t68 + t117, (-t126 + t131) * t65 + (t33 - t130) * t106 + t118, -t87 + t139, -t17 * pkin(3) + pkin(7) * t139 - t33 * t134 - t27 * t34, t48, t60, -t132, 0, t59, t49, t50 * t97 + t9 + (-t38 * t62 + t96) * t68 + t117, -t87 + t78, (-t15 - t33 - t125) * t106 + (t114 * t62 + t96) * t65, t74 * pkin(7) - t114 * t15 - t84 * t33 + t8 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t43, 0, t51, 0, 0, -t27 * t123 - t11, -t27 * t122 - t13, 0, 0, -t51, 0, t43, 0, 0, t51, -t11 + (-t15 * t65 + t39 * t68) * t62, ((t20 - t102) * t65 + (-t19 + t89) * t68) * t62, qJD(5) * t138 + t13 + (t15 * t68 + t39 * t65) * t62, -t19 * t120 - t7 * pkin(4) + t4 * qJ(5) - t15 * t39 + (qJD(5) + t121) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, 0, -t63 * t61 - t71, t15 * t123 - t105 + t7;];
tauc_reg = t3;
