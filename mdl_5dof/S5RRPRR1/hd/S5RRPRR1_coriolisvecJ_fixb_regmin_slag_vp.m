% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:18
% EndTime: 2021-01-15 21:13:25
% DurationCPUTime: 1.52s
% Computational Cost: add. (1335->198), mult. (3398->301), div. (0->0), fcn. (2303->6), ass. (0->132)
t162 = cos(qJ(4));
t119 = qJD(1) * t162;
t89 = sin(qJ(2));
t141 = qJD(1) * t89;
t88 = sin(qJ(4));
t91 = cos(qJ(2));
t51 = -t91 * t119 + t88 * t141;
t83 = qJD(2) + qJD(4);
t154 = t51 * t83;
t58 = t162 * t89 + t88 * t91;
t34 = t83 * t58;
t28 = t34 * qJD(1);
t50 = qJD(5) + t51;
t166 = qJD(5) - t50;
t92 = pkin(1) + pkin(2);
t165 = t50 ^ 2;
t140 = qJD(1) * t91;
t53 = -t119 * t89 - t140 * t88;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t106 = t90 * t53 - t87 * t83;
t13 = -qJD(5) * t106 - t154 * t87;
t144 = pkin(3) + qJ(3);
t64 = t144 * t89;
t65 = t144 * t91;
t100 = -t162 * t64 - t88 * t65;
t116 = qJD(2) * t144;
t47 = t91 * qJD(3) - t116 * t89;
t134 = t89 * qJD(3);
t48 = -t116 * t91 - t134;
t10 = t100 * qJD(4) + t162 * t47 + t88 * t48;
t59 = qJD(1) * t64;
t86 = qJD(2) * pkin(1);
t46 = qJD(2) * pkin(2) - t59 + t86;
t125 = t162 * t46;
t60 = qJD(1) * t65;
t147 = t88 * t60;
t25 = -t125 + t147;
t66 = t92 * t91;
t56 = -qJD(1) * t66 + qJD(3);
t32 = t53 * pkin(4) + t56;
t99 = t162 * t91 - t88 * t89;
t33 = t83 * t99;
t40 = t162 * t65 - t88 * t64;
t42 = -t58 * pkin(4) - t66;
t118 = t162 * qJD(4);
t137 = qJD(4) * t88;
t44 = t48 * qJD(1);
t45 = t47 * qJD(1);
t98 = -t60 * t137 + t162 * t45 + t88 * t44;
t6 = t118 * t46 + t98;
t120 = -t162 * t44 + t88 * t45;
t124 = t162 * t60;
t26 = t88 * t46 + t124;
t7 = qJD(4) * t26 + t120;
t164 = -(qJD(5) * t42 + t10) * t50 + (qJD(5) * t32 + t6) * t99 + t25 * t33 - t40 * t28 + t7 * t58;
t85 = t91 ^ 2;
t163 = 0.2e1 * t85;
t135 = qJD(5) * t90;
t136 = qJD(5) * t87;
t12 = t83 * t135 + t136 * t53 - t154 * t90;
t161 = t12 * t87;
t160 = t25 * t58;
t36 = -t87 * t53 - t90 * t83;
t159 = t36 * t50;
t158 = t106 * t50;
t157 = t106 * t53;
t156 = t42 * t28;
t155 = t50 * t53;
t153 = t53 * t36;
t152 = t53 * t51;
t151 = t53 * t83;
t150 = t56 * t51;
t148 = t87 * t28;
t24 = t90 * t28;
t94 = qJD(1) ^ 2;
t146 = t91 * t94;
t93 = qJD(2) ^ 2;
t145 = t93 * t91;
t131 = qJD(1) * qJD(2);
t121 = t89 * t131;
t76 = pkin(1) * t121;
t54 = pkin(2) * t121 + t76;
t61 = t92 * t141;
t139 = qJD(2) * t89;
t62 = t92 * t139;
t84 = t89 ^ 2;
t143 = t84 - t85;
t142 = t89 * qJ(3);
t138 = qJD(2) * t91;
t130 = pkin(1) * t140;
t70 = qJD(3) - t130;
t133 = -qJD(3) + t70;
t132 = qJ(3) * qJD(1);
t63 = -t132 * t89 + t86;
t129 = t63 * t138;
t127 = 0.2e1 * t131;
t126 = t92 * t162;
t123 = qJ(3) * t138;
t122 = t25 * (-t50 + t51);
t115 = t90 * t50;
t71 = t88 * t92 + pkin(4);
t113 = t51 * pkin(4) + qJD(5) * t71 + t61;
t112 = (-qJD(3) - t70) * qJD(1);
t111 = t91 * t127;
t20 = t83 * pkin(4) + t26;
t9 = t90 * t20 + t87 * t32;
t110 = t25 * t135 - t9 * t53 + t7 * t87;
t30 = -t88 * t59 + t124;
t109 = t137 * t92 - t30;
t108 = t87 * t20 - t90 * t32;
t107 = t25 * t51 - t71 * t28;
t105 = t24 + (-t51 * t87 - t136) * t50;
t104 = -t108 * t53 + t25 * t136 - t7 * t90;
t103 = t56 * t53 - t120;
t31 = -t162 * t59 - t147;
t102 = -t118 * t92 + t31;
t101 = -t136 * t58 + t90 * t33;
t96 = -t165 * t90 - t148;
t95 = qJ(3) ^ 2;
t55 = (-t123 - t134) * qJD(1);
t19 = -t33 * pkin(4) + t62;
t18 = -t51 ^ 2 + t53 ^ 2;
t17 = pkin(4) * t154 + t54;
t16 = t90 * t17;
t15 = -t151 - t28;
t11 = t40 * qJD(4) - t162 * t48 + t88 * t47;
t4 = t115 * t50 + t148 - t157;
t3 = t105 - t153;
t2 = -t106 * t115 + t161;
t1 = (t12 - t159) * t90 + (-t13 + t158) * t87;
t5 = [0, 0, 0, t89 * t111, -t143 * t127, t145, -t93 * t89, 0, 0, 0, -qJ(3) * t145 + (-0.3e1 * t130 + t133) * t139, t93 * t142 + (t133 * t91 + (0.2e1 * t84 - t85) * qJD(1) * pkin(1)) * qJD(2), -t129 - t55 * t89 + (-0.2e1 * t89 * t123 + (t84 + t163) * qJD(3)) * qJD(1), (qJD(1) * qJD(3) * t163 - t129) * qJ(3) + (-qJ(3) * t55 - qJD(3) * t63 + (-0.2e1 * t95 * t140 + (t70 - t130) * pkin(1)) * qJD(2)) * t89, -t154 * t58 - t53 * t33, -t154 * t99 - t58 * t28 - t33 * t51 + t53 * t34, t33 * t83, -t34 * t83, 0, -t11 * t83 - t66 * t28 + t56 * t34 + t62 * t51 - t54 * t99, -t10 * t83 + t154 * t66 + t56 * t33 - t62 * t53 + t54 * t58, t12 * t90 * t58 - t101 * t106, (t106 * t87 - t36 * t90) * t33 + (-t161 - t13 * t90 + (t106 * t90 + t36 * t87) * qJD(5)) * t58, t101 * t50 - t106 * t34 - t12 * t99 + t58 * t24, -t58 * t148 + t13 * t99 - t36 * t34 + (-t135 * t58 - t87 * t33) * t50, -t28 * t99 + t50 * t34, t11 * t36 - t100 * t13 - t16 * t99 - t108 * t34 + (t19 * t50 + t156 + (t20 * t99 - t40 * t50 + t160) * qJD(5)) * t90 + t164 * t87, -t11 * t106 - t100 * t12 - t9 * t34 + (-(-qJD(5) * t40 + t19) * t50 - t156 + (-qJD(5) * t20 + t17) * t99 - qJD(5) * t160) * t87 + t164 * t90; 0, 0, 0, -t89 * t146, t143 * t94, 0, 0, 0, 0, 0, (pkin(1) * t146 + t112) * t89, -t84 * t94 * pkin(1) + t112 * t91, (t94 * t142 + (t63 - t86) * qJD(1)) * t91, (t89 * t94 * t95 + t132 * t63) * t91 + (-t141 * t70 + t55) * pkin(1), -t152, t18, 0, t15, 0, t30 * t83 - t61 * t51 + (-t124 + (-t83 * t92 - t46) * t88) * qJD(4) + t103, t31 * t83 + t150 + t61 * t53 + (-t126 * t83 - t125) * qJD(4) - t98, t2, t1, t4, t3, t155, -t13 * t126 + t107 * t87 + t109 * t36 + (t102 * t87 - t113 * t90) * t50 + t104, -t12 * t126 + t107 * t90 - t109 * t106 + (t102 * t90 + t113 * t87) * t50 + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121, t111, (-t84 - t85) * t94, -t85 * t94 * qJ(3) + t63 * t141 + t76, 0, 0, 0, 0, 0, t28 - t151, -0.2e1 * t154, 0, 0, 0, 0, 0, t105 + t153, t96 - t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t18, 0, t15, 0, t103 + (-qJD(4) + t83) * t26, -t25 * t83 + t150 - t6, t2, t1, t4, t3, t155, pkin(4) * t96 + t122 * t87 - t26 * t36 + t104, t26 * t106 + t90 * t122 + (t165 * t87 - t24) * pkin(4) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t36, t106 ^ 2 - t36 ^ 2, t12 + t159, -t13 - t158, t28, t25 * t106 - t166 * t9 - t87 * t6 + t16, t166 * t108 - t87 * t17 + t25 * t36 - t90 * t6;];
tauc_reg = t5;
