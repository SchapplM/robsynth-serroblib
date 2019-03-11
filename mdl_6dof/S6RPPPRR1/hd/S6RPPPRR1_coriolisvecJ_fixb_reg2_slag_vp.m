% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:14
% EndTime: 2019-03-09 01:30:18
% DurationCPUTime: 1.32s
% Computational Cost: add. (1701->206), mult. (3219->291), div. (0->0), fcn. (1727->6), ass. (0->127)
t104 = qJD(3) * qJD(1);
t54 = sin(qJ(5));
t41 = sin(pkin(9)) * pkin(1) + qJ(3);
t150 = qJD(1) * t41;
t34 = qJD(4) + t150;
t26 = -qJD(1) * pkin(7) + t34;
t56 = cos(qJ(5));
t70 = t54 * qJD(2) - t56 * t26;
t11 = -t70 * qJD(5) + t54 * t104;
t82 = pkin(5) * t56 + pkin(8) * t54;
t29 = t82 * qJD(5) + qJD(4);
t23 = t29 * qJD(1);
t53 = sin(qJ(6));
t55 = cos(qJ(6));
t21 = t56 * qJD(2) + t54 * t26;
t16 = qJD(5) * pkin(8) + t21;
t39 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t25 = t54 * pkin(5) - t56 * pkin(8) + t39;
t17 = t25 * qJD(1) - qJD(3);
t75 = t53 * t16 - t55 * t17;
t1 = -t75 * qJD(6) + t55 * t11 + t53 * t23;
t108 = t54 * qJD(1);
t40 = qJD(6) + t108;
t153 = t75 * t40 + t1;
t109 = t53 * qJD(5);
t111 = qJD(6) * t56;
t152 = -t54 * t109 + t55 * t111;
t119 = qJD(1) * t56;
t92 = t55 * t119;
t32 = t92 + t109;
t135 = t32 * t40;
t115 = qJD(6) * t32;
t105 = qJD(1) * qJD(5);
t90 = t54 * t105;
t19 = -t53 * t90 + t115;
t151 = t19 - t135;
t6 = t55 * t16 + t53 * t17;
t2 = -qJD(6) * t6 - t53 * t11 + t55 * t23;
t149 = -t6 * t40 - t2;
t107 = t55 * qJD(5);
t94 = t53 * t111;
t64 = t54 * t107 + t94;
t18 = t64 * qJD(1) - qJD(6) * t107;
t30 = t53 * t119 - t107;
t68 = t30 * t40;
t148 = -t18 + t68;
t117 = qJD(5) * t56;
t123 = t32 * t117 - t18 * t54;
t50 = t56 ^ 2;
t120 = qJD(1) * t50;
t147 = -t107 * (-t40 * t54 + t120) + t40 * t94;
t42 = t56 * t104;
t12 = t21 * qJD(5) - t42;
t78 = -t53 * t75 - t55 * t6;
t146 = qJD(5) * t78 + t12;
t140 = t12 * t56;
t60 = -(-t21 * t56 - t54 * t70) * qJD(5) + t11 * t54 - t140;
t142 = t12 * t53;
t141 = t12 * t55;
t15 = -qJD(5) * pkin(5) + t70;
t139 = t15 * t53;
t138 = t15 * t55;
t137 = t30 * t56;
t136 = t32 * t30;
t134 = t32 * t56;
t133 = t40 * t53;
t132 = t40 * t55;
t49 = t54 ^ 2;
t58 = qJD(1) ^ 2;
t131 = t49 * t58;
t130 = t53 * t54;
t129 = t54 * t19;
t128 = t54 * t55;
t127 = t56 * t18;
t126 = t56 * t19;
t57 = qJD(5) ^ 2;
t125 = t57 * t54;
t124 = t57 * t56;
t122 = -t57 - t58;
t121 = qJD(1) * t39;
t118 = qJD(5) * t54;
t116 = qJD(6) * t30;
t114 = qJD(6) * t53;
t113 = qJD(6) * t54;
t112 = qJD(6) * t55;
t27 = -qJD(3) + t121;
t110 = t27 * qJD(1);
t106 = qJD(3) - t27;
t103 = t40 * t130;
t102 = t40 * t128;
t101 = t56 * t58 * t54;
t100 = t53 * t120;
t98 = t30 * t118;
t97 = t32 * t118;
t96 = t32 * t111;
t95 = t40 * t114;
t91 = 0.2e1 * qJD(4) * qJD(1);
t89 = t56 * t105;
t88 = t27 + t121;
t87 = t40 + t108;
t86 = -t18 + t116;
t84 = qJD(1) + t113;
t83 = t54 * t89;
t38 = -pkin(7) + t41;
t81 = -t38 * t113 + t29;
t80 = t40 * t112 + t53 * t89;
t79 = t53 * t6 - t55 * t75;
t74 = t21 * t54 - t56 * t70;
t72 = qJD(3) + t88;
t67 = -t38 * t57 + t91;
t66 = t152 * t40;
t65 = -pkin(8) * t117 + t15 * t54;
t63 = qJD(3) * t54 + qJD(6) * t25 + t38 * t117;
t62 = t78 * qJD(6) - t1 * t53 - t2 * t55;
t61 = -t79 * qJD(6) + t1 * t55 - t2 * t53;
t59 = qJD(5) * t15 + t61;
t47 = 0.2e1 * t104;
t46 = t50 * t58;
t33 = t82 * qJD(1);
t13 = t55 * t126;
t10 = t38 * t128 + t53 * t25;
t9 = -t38 * t130 + t55 * t25;
t8 = t53 * t33 - t55 * t70;
t7 = t55 * t33 + t53 * t70;
t4 = -t63 * t53 + t81 * t55;
t3 = t81 * t53 + t63 * t55;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0.2e1 * t150 * qJD(3), 0, 0, 0, 0, 0, 0, 0, t47, t91, t34 * qJD(3) + t27 * qJD(4) + (qJD(3) * t41 + qJD(4) * t39) * qJD(1), -0.2e1 * t83, 0.2e1 * (t49 - t50) * t105, -t125, 0.2e1 * t83, -t124, 0, t117 * t72 + t54 * t67, -t118 * t72 + t56 * t67 (-t49 - t50) * t104 - t60, qJD(3) * t74 + qJD(4) * t88 + t38 * t60, -t127 * t55 - t32 * t64, -t13 + (-t96 + t98) * t55 + (t97 + (t18 + t116) * t56) * t53, t123 - t147, t53 * t126 + t152 * t30, -t129 + (-t100 - t137) * qJD(5) - t66, t87 * t117, t4 * t40 + (t2 + (t30 * t38 - t139) * qJD(5)) * t54 + (t15 * t112 - qJD(3) * t30 + t142 - t19 * t38 + (qJD(1) * t9 - t75) * qJD(5)) * t56, -t3 * t40 + (-t1 + (t32 * t38 - t138) * qJD(5)) * t54 + (-t15 * t114 - qJD(3) * t32 + t141 + t18 * t38 + (-qJD(1) * t10 - t6) * qJD(5)) * t56, -t10 * t19 + t118 * t79 + t9 * t18 - t3 * t30 - t4 * t32 + t56 * t62, -t38 * t140 + t1 * t10 + t2 * t9 + t6 * t3 - t75 * t4 + (-qJD(3) * t56 + t118 * t38) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t125, 0, -qJD(5) * t74 + t11 * t56 + t12 * t54, 0, 0, 0, 0, 0, 0, t129 + (-t100 + t137) * qJD(5) - t66, t123 + t147, -t13 + (t96 + t98) * t55 + (t56 * t86 - t97) * t53, t146 * t54 + t59 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t150 * qJD(1), 0, 0, 0, 0, 0, 0, 0, -t58, 0 (-qJD(4) - t34) * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t89, 0.2e1 * t90, t46 + t131 (-qJD(4) - t74) * qJD(1), 0, 0, 0, 0, 0, 0, t95 + (t103 + (t30 - t107) * t56) * qJD(1) (t102 + t134) * qJD(1) + t80, t148 * t55 + t151 * t53 (t15 * t56 + t54 * t78) * qJD(1) + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t106 * qJD(1), 0, 0, 0, 0, 0, 0, t122 * t54, t122 * t56, 0, t60 - t110, 0, 0, 0, 0, 0, 0, -t126 - t84 * t132 + (-t53 * t56 * t87 + t30 * t54) * qJD(5), t127 + t84 * t133 + (-t56 * t132 + (t32 - t92) * t54) * qJD(5) (-t117 * t30 + t32 * t84 - t129) * t55 + (t30 * t84 + t123) * t53, -t79 * qJD(1) - t146 * t56 + t59 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t46 - t131, 0, -t101, 0, 0, -t110 * t56 + t42, -t106 * t108, 0, 0, t132 * t32 - t18 * t53 (-t18 - t68) * t55 + (-t19 - t135) * t53 (t102 - t134) * qJD(1) + t80, -t19 * t55 + t53 * t68, -t95 + (-t103 + (t30 + t107) * t56) * qJD(1), -t40 * t119, -pkin(5) * t19 - t141 - t21 * t30 - t7 * t40 + (-pkin(8) * t132 + t139) * qJD(6) + (t53 * t65 + t56 * t75) * qJD(1), pkin(5) * t18 + t142 - t21 * t32 + t8 * t40 + (pkin(8) * t133 + t138) * qJD(6) + (t55 * t65 + t56 * t6) * qJD(1), t8 * t30 + t7 * t32 + ((-t19 + t115) * pkin(8) + t153) * t55 + (pkin(8) * t86 + t149) * t53, -t12 * pkin(5) + pkin(8) * t61 - t15 * t21 - t6 * t8 + t7 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t30 ^ 2 + t32 ^ 2, t148, -t136, -t151, t89, -t15 * t32 - t149, t15 * t30 - t153, 0, 0;];
tauc_reg  = t5;
