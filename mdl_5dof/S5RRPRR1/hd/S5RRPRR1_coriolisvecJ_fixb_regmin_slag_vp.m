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
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
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
% StartTime: 2019-07-18 17:22:26
% EndTime: 2019-07-18 17:22:32
% DurationCPUTime: 1.43s
% Computational Cost: add. (1320->186), mult. (3340->282), div. (0->0), fcn. (2280->6), ass. (0->126)
t157 = cos(qJ(4));
t117 = qJD(1) * t157;
t89 = sin(qJ(2));
t137 = qJD(1) * t89;
t88 = sin(qJ(4));
t91 = cos(qJ(2));
t162 = t91 * t117 - t88 * t137;
t83 = qJD(2) + qJD(4);
t163 = t162 * t83;
t58 = t157 * t89 + t88 * t91;
t34 = t83 * t58;
t28 = t34 * qJD(1);
t50 = qJD(5) - t162;
t161 = qJD(5) - t50;
t92 = pkin(1) + pkin(2);
t160 = t50 ^ 2;
t136 = qJD(1) * t91;
t53 = -t89 * t117 - t88 * t136;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t106 = t90 * t53 - t87 * t83;
t13 = -t106 * qJD(5) + t163 * t87;
t140 = pkin(3) + qJ(3);
t64 = t140 * t89;
t65 = t140 * t91;
t100 = -t157 * t64 - t88 * t65;
t114 = qJD(2) * t140;
t47 = t91 * qJD(3) - t89 * t114;
t130 = t89 * qJD(3);
t48 = -t91 * t114 - t130;
t10 = t100 * qJD(4) + t157 * t47 + t88 * t48;
t59 = qJD(1) * t64;
t86 = qJD(2) * pkin(1);
t46 = qJD(2) * pkin(2) - t59 + t86;
t123 = t157 * t46;
t60 = qJD(1) * t65;
t142 = t88 * t60;
t25 = -t123 + t142;
t66 = t92 * t91;
t56 = -qJD(1) * t66 + qJD(3);
t32 = t53 * pkin(4) + t56;
t99 = t157 * t91 - t88 * t89;
t33 = t83 * t99;
t40 = t157 * t65 - t88 * t64;
t42 = -t58 * pkin(4) - t66;
t116 = t157 * qJD(4);
t133 = qJD(4) * t88;
t44 = t47 * qJD(1);
t45 = t48 * qJD(1);
t98 = -t60 * t133 + t157 * t44 + t88 * t45;
t6 = t46 * t116 + t98;
t118 = -t157 * t45 + t88 * t44;
t122 = t157 * t60;
t26 = t88 * t46 + t122;
t7 = t26 * qJD(4) + t118;
t159 = -(qJD(5) * t42 + t10) * t50 + (qJD(5) * t32 + t6) * t99 + t25 * t33 - t40 * t28 + t7 * t58;
t85 = t91 ^ 2;
t158 = 0.2e1 * t85;
t131 = qJD(5) * t90;
t132 = qJD(5) * t87;
t12 = t83 * t131 + t53 * t132 + t163 * t90;
t156 = t12 * t87;
t155 = t25 * t58;
t36 = -t87 * t53 - t90 * t83;
t154 = t36 * t50;
t153 = t106 * t50;
t152 = t106 * t53;
t151 = t42 * t28;
t150 = t50 * t53;
t148 = t53 * t36;
t147 = t53 * t162;
t146 = t53 * t83;
t145 = t56 * t162;
t143 = t87 * t28;
t94 = qJD(1) ^ 2;
t141 = t89 * t94;
t24 = t90 * t28;
t128 = qJD(1) * qJD(2);
t119 = t89 * t128;
t76 = pkin(1) * t119;
t54 = pkin(2) * t119 + t76;
t61 = t92 * qJD(2) * t89;
t62 = t92 * t137;
t84 = t89 ^ 2;
t139 = t84 - t85;
t138 = qJ(3) * t94;
t134 = qJD(2) * t91;
t129 = qJ(3) * qJD(1);
t127 = pkin(1) * t136;
t63 = -t89 * t129 + t86;
t126 = t63 * t134;
t124 = t92 * t157;
t121 = qJ(3) * t134;
t120 = t25 * (-t50 - t162);
t113 = t90 * t50;
t71 = t88 * t92 + pkin(4);
t111 = -pkin(4) * t162 + qJD(5) * t71 + t62;
t20 = t83 * pkin(4) + t26;
t9 = t90 * t20 + t87 * t32;
t110 = t25 * t131 - t9 * t53 + t7 * t87;
t30 = -t88 * t59 + t122;
t109 = t92 * t133 - t30;
t108 = t87 * t20 - t90 * t32;
t107 = -t162 * t25 - t71 * t28;
t105 = t24 + (t162 * t87 - t132) * t50;
t104 = -t108 * t53 + t25 * t132 - t7 * t90;
t103 = t56 * t53 - t118;
t31 = -t157 * t59 - t142;
t102 = -t92 * t116 + t31;
t101 = -t58 * t132 + t90 * t33;
t96 = -t160 * t90 - t143;
t95 = qJ(3) ^ 2;
t93 = qJD(2) ^ 2;
t70 = qJD(3) - t127;
t55 = (-t121 - t130) * qJD(1);
t19 = -t33 * pkin(4) + t61;
t18 = -t162 ^ 2 + t53 ^ 2;
t17 = -pkin(4) * t163 + t54;
t16 = t90 * t17;
t15 = -t146 - t28;
t11 = t40 * qJD(4) - t157 * t48 + t88 * t47;
t4 = t50 * t113 + t143 - t152;
t3 = t105 - t148;
t2 = -t106 * t113 + t156;
t1 = (t12 - t154) * t90 + (-t13 + t153) * t87;
t5 = [0, 0, 0, 0.2e1 * t91 * t119, -0.2e1 * t139 * t128, t93 * t91, -t93 * t89, 0, 0, 0, -t126 - t55 * t89 + (-0.2e1 * t89 * t121 + (t84 + t158) * qJD(3)) * qJD(1), (qJD(1) * qJD(3) * t158 - t126) * qJ(3) + (-qJ(3) * t55 - qJD(3) * t63 + (-0.2e1 * t95 * t136 + (t70 - t127) * pkin(1)) * qJD(2)) * t89, t163 * t58 - t53 * t33, t162 * t33 + t163 * t99 - t58 * t28 + t53 * t34, t33 * t83, -t34 * t83, 0, -t11 * t83 - t162 * t61 - t66 * t28 + t56 * t34 - t54 * t99, -t10 * t83 - t163 * t66 + t56 * t33 - t61 * t53 + t54 * t58, t12 * t90 * t58 - t101 * t106, (t106 * t87 - t36 * t90) * t33 + (-t156 - t13 * t90 + (t106 * t90 + t36 * t87) * qJD(5)) * t58, t101 * t50 - t106 * t34 - t12 * t99 + t58 * t24, -t58 * t143 + t13 * t99 - t36 * t34 + (-t58 * t131 - t87 * t33) * t50, -t28 * t99 + t50 * t34, t11 * t36 - t100 * t13 - t16 * t99 - t108 * t34 + (t19 * t50 + t151 + (t20 * t99 - t40 * t50 + t155) * qJD(5)) * t90 + t159 * t87, -t11 * t106 - t100 * t12 - t9 * t34 + (-(-qJD(5) * t40 + t19) * t50 - t151 + (-qJD(5) * t20 + t17) * t99 - qJD(5) * t155) * t87 + t159 * t90; 0, 0, 0, -t91 * t141, t139 * t94, 0, 0, 0, 0, 0, (t89 * t138 + (t63 - t86) * qJD(1)) * t91, (t63 * t129 + t95 * t141) * t91 + (-t70 * t137 + t55) * pkin(1), t147, t18, 0, t15, 0, t30 * t83 + t62 * t162 + (-t122 + (-t83 * t92 - t46) * t88) * qJD(4) + t103, t31 * t83 - t145 + t62 * t53 + (-t83 * t124 - t123) * qJD(4) - t98, t2, t1, t4, t3, t150, -t13 * t124 + t107 * t87 + t109 * t36 + (t102 * t87 - t111 * t90) * t50 + t104, -t12 * t124 + t107 * t90 - t109 * t106 + (t102 * t90 + t111 * t87) * t50 + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t84 - t85) * t94, t63 * t137 - t85 * t138 + t76, 0, 0, 0, 0, 0, t28 - t146, 0.2e1 * t163, 0, 0, 0, 0, 0, t105 + t148, t96 - t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t18, 0, t15, 0, t103 + (-qJD(4) + t83) * t26, -t25 * t83 - t145 - t6, t2, t1, t4, t3, t150, t96 * pkin(4) + t87 * t120 - t26 * t36 + t104, t26 * t106 + t90 * t120 + (t160 * t87 - t24) * pkin(4) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t36, t106 ^ 2 - t36 ^ 2, t12 + t154, -t13 - t153, t28, t25 * t106 - t161 * t9 - t87 * t6 + t16, t161 * t108 - t87 * t17 + t25 * t36 - t90 * t6;];
tauc_reg  = t5;
