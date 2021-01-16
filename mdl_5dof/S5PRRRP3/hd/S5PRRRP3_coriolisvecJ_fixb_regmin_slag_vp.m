% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP3
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
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:27
% EndTime: 2021-01-15 16:23:31
% DurationCPUTime: 0.70s
% Computational Cost: add. (1086->149), mult. (2736->200), div. (0->0), fcn. (1814->4), ass. (0->106)
t126 = cos(qJ(4));
t70 = sin(qJ(4));
t71 = sin(qJ(3));
t72 = cos(qJ(3));
t50 = t126 * t71 + t70 * t72;
t109 = qJD(2) * t50;
t110 = t109 * qJ(5);
t106 = t71 * qJD(1);
t129 = pkin(6) + pkin(7);
t55 = t129 * t72;
t39 = qJD(2) * t55 + t106;
t33 = t70 * t39;
t113 = qJD(3) * pkin(3);
t99 = qJD(2) * t129;
t38 = t72 * qJD(1) - t71 * t99;
t36 = t38 + t113;
t95 = t126 * t36 - t33;
t6 = -t110 + t95;
t67 = qJD(3) + qJD(4);
t104 = qJD(2) * qJD(3);
t131 = -0.2e1 * t104;
t130 = t109 ^ 2;
t5 = t67 * pkin(4) + t6;
t128 = t5 - t6;
t127 = pkin(3) * t67;
t100 = t126 * t72;
t121 = t70 * t71;
t88 = t67 * t121;
t93 = t126 * qJD(4);
t23 = -qJD(3) * t100 - t72 * t93 + t88;
t21 = t23 * t67;
t108 = qJD(2) * t71;
t101 = t70 * t108;
t89 = qJD(2) * t100;
t41 = -t89 + t101;
t125 = t41 * t67;
t124 = t109 * t41;
t65 = -t72 * pkin(3) - pkin(2);
t53 = t65 * qJD(2);
t122 = t53 * t109;
t74 = qJD(2) ^ 2;
t120 = t72 * t74;
t73 = qJD(3) ^ 2;
t119 = t73 * t71;
t118 = t73 * t72;
t24 = t67 * t50;
t20 = t24 * qJD(2);
t117 = -t50 * t20 + t23 * t41;
t116 = t126 * t38 - t33;
t115 = t67 * t89;
t114 = t71 ^ 2 - t72 ^ 2;
t19 = qJD(2) * t88 - t115;
t112 = t19 * qJ(5);
t111 = t41 * qJ(5);
t107 = qJD(4) * t70;
t92 = t41 * pkin(4) + qJD(5);
t25 = t53 + t92;
t105 = qJD(5) + t25;
t103 = t71 * t113;
t102 = pkin(3) * t108;
t35 = t126 * t39;
t97 = t71 * t104;
t62 = pkin(3) * t97;
t15 = t20 * pkin(4) + t62;
t98 = qJD(3) * t129;
t30 = t38 * qJD(3);
t31 = (-t72 * t99 - t106) * qJD(3);
t96 = t126 * t31 - t70 * t30;
t94 = -t70 * t38 - t35;
t91 = pkin(2) * t131;
t49 = -t100 + t121;
t87 = t109 * t24 - t49 * t19;
t86 = -t70 * t36 - t35;
t54 = t129 * t71;
t85 = -t126 * t55 + t70 * t54;
t51 = t71 * t98;
t52 = t72 * t98;
t84 = -t55 * t107 - t126 * t51 - t70 * t52 - t54 * t93;
t83 = -t67 * t101 + t115;
t82 = t86 * qJD(4) + t96;
t81 = t85 * qJD(4) - t126 * t52 + t70 * t51;
t80 = -t39 * t107 + t126 * t30 + t70 * t31 + t36 * t93;
t79 = t82 + t112;
t78 = t53 * t41 - t80;
t77 = -t20 * qJ(5) + t80;
t76 = (-t35 + (-t36 - t127) * t70) * qJD(4) + t96;
t75 = t105 * t41 - t77;
t64 = t126 * pkin(3) + pkin(4);
t56 = t93 * t127;
t40 = t41 ^ 2;
t32 = t49 * pkin(4) + t65;
t27 = pkin(4) * t109 + t102;
t22 = t24 * t67;
t18 = t24 * pkin(4) + t103;
t17 = -t49 * qJ(5) - t85;
t16 = -t50 * qJ(5) - t126 * t54 - t70 * t55;
t13 = -t40 + t130;
t10 = t83 + t125;
t9 = -t110 + t116;
t8 = t94 + t111;
t7 = -t86 - t111;
t4 = t23 * qJ(5) - t50 * qJD(5) + t81;
t3 = -t24 * qJ(5) - t49 * qJD(5) + t84;
t2 = -qJD(5) * t109 + t79;
t1 = -t41 * qJD(5) + t77;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t118, 0, 0, 0, 0, 0, -t22, t21, -t22, t21, t87 + t117, t1 * t50 - t2 * t49 - t7 * t23 - t5 * t24; 0, 0, 0, 0, 0.2e1 * t72 * t97, t114 * t131, t118, -t119, 0, -pkin(6) * t118 + t71 * t91, pkin(6) * t119 + t72 * t91, -t109 * t23 - t19 * t50, -t87 + t117, -t21, -t22, 0, t41 * t103 + t65 * t20 + t53 * t24 + t49 * t62 + t81 * t67, 0.2e1 * t109 * t103 - t65 * t19 - t53 * t23 - t84 * t67, t15 * t49 + t18 * t41 + t32 * t20 + t25 * t24 + t4 * t67, t109 * t18 + t15 * t50 - t32 * t19 - t25 * t23 - t3 * t67, -t1 * t49 - t109 * t4 + t16 * t19 - t17 * t20 - t2 * t50 + t5 * t23 - t7 * t24 - t3 * t41, t1 * t17 + t15 * t32 + t2 * t16 + t25 * t18 + t7 * t3 + t5 * t4; 0, 0, 0, 0, -t71 * t120, t114 * t74, 0, 0, 0, t74 * pkin(2) * t71, pkin(2) * t120, t124, t13, t10, 0, 0, -t41 * t102 - t94 * t67 - t122 + t76, -t102 * t109 + t116 * t67 - t56 + t78, -t105 * t109 - t27 * t41 - t8 * t67 + t112 + t76, -t109 * t27 + t9 * t67 - t56 + t75, t64 * t19 + (t7 + t8) * t109 + (-t5 + t9) * t41 + (-t20 * t70 + (t109 * t70 - t126 * t41) * qJD(4)) * pkin(3), t2 * t64 - t25 * t27 - t5 * t8 - t7 * t9 + (t1 * t70 + (t126 * t7 - t5 * t70) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t13, t10, 0, 0, -t86 * t67 - t122 + t82, t95 * t67 + t78, t7 * t67 + (-t25 - t92) * t109 + t79, -t130 * pkin(4) + t6 * t67 + t75, t19 * pkin(4) - t128 * t41, t128 * t7 + (-t109 * t25 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t109 + t20, t83 - t125, -t40 - t130, t109 * t5 + t7 * t41 + t15;];
tauc_reg = t11;
