% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:18
% EndTime: 2021-01-15 21:13:28
% DurationCPUTime: 2.27s
% Computational Cost: add. (1760->279), mult. (4019->384), div. (0->0), fcn. (2792->10), ass. (0->172)
t113 = qJD(2) + qJD(4);
t120 = sin(qJ(4));
t124 = cos(qJ(2));
t222 = cos(qJ(4));
t166 = qJD(1) * t222;
t121 = sin(qJ(2));
t187 = qJD(1) * t121;
t230 = -t120 * t187 + t124 * t166;
t202 = t230 * t113;
t116 = qJ(2) + qJ(4);
t109 = sin(t116);
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t148 = g(1) * t125 + g(2) * t122;
t231 = t148 * t109;
t108 = t121 * qJDD(1);
t158 = qJDD(1) * t222;
t64 = t120 * t124 + t222 * t121;
t38 = t113 * t64;
t22 = qJD(1) * t38 + t108 * t120 - t124 * t158;
t126 = pkin(1) + pkin(2);
t57 = qJD(5) - t230;
t157 = t57 ^ 2;
t229 = pkin(4) * t157;
t110 = cos(t116);
t217 = g(3) * t110;
t117 = qJDD(2) * pkin(1);
t205 = pkin(3) + qJ(3);
t159 = qJD(2) * t205;
t55 = -t121 * qJD(3) - t124 * t159;
t72 = t205 * t121;
t29 = qJDD(2) * pkin(2) + qJD(1) * t55 - qJDD(1) * t72 + t117;
t54 = t124 * qJD(3) - t121 * t159;
t73 = t205 * t124;
t33 = qJD(1) * t54 + qJDD(1) * t73;
t160 = t120 * t33 - t222 * t29;
t66 = qJD(1) * t73;
t172 = t222 * t66;
t118 = qJD(2) * pkin(1);
t65 = qJD(1) * t72;
t49 = qJD(2) * pkin(2) + t118 - t65;
t31 = t120 * t49 + t172;
t7 = qJD(4) * t31 + t160;
t228 = t7 + t217;
t67 = t126 * t187;
t87 = t120 * t126 + pkin(4);
t226 = (-pkin(4) * t230 + qJD(5) * t87 + t67) * t57;
t112 = qJDD(2) + qJDD(4);
t119 = sin(qJ(5));
t123 = cos(qJ(5));
t176 = t124 * qJDD(1);
t21 = t120 * t176 + t121 * t158 + t202;
t186 = qJD(1) * t124;
t60 = -t120 * t186 - t121 * t166;
t44 = t119 * t113 - t123 * t60;
t9 = qJD(5) * t44 - t123 * t112 + t119 * t21;
t102 = g(3) * t109;
t74 = t126 * t124;
t62 = -qJD(1) * t74 + qJD(3);
t225 = t148 * t110 - t230 * t62 + t102;
t136 = -t120 * t121 + t222 * t124;
t138 = -t120 * t73 - t222 * t72;
t16 = t138 * qJD(4) + t120 * t55 + t222 * t54;
t184 = qJD(4) * t120;
t149 = -t66 * t184 + t222 * t33;
t165 = qJD(4) * t222;
t131 = t120 * t29 + t165 * t49 + t149;
t36 = t60 * pkin(4) + t62;
t162 = t112 * pkin(4) + qJD(5) * t36 + t131;
t20 = qJDD(5) + t22;
t173 = t222 * t49;
t203 = t120 * t66;
t30 = -t173 + t203;
t37 = t113 * t136;
t45 = -t120 * t72 + t222 * t73;
t47 = -t64 * pkin(4) - t74;
t224 = -(qJD(5) * t47 + t16) * t57 + t162 * t136 - t45 * t20 + t30 * t37 + t7 * t64;
t115 = t124 ^ 2;
t223 = 0.2e1 * t115;
t221 = g(1) * t122;
t218 = g(2) * t125;
t216 = g(3) * t124;
t215 = t30 * t64;
t41 = -t123 * t113 - t119 * t60;
t214 = t41 * t57;
t213 = t44 * t57;
t212 = t44 * t60;
t211 = t47 * t20;
t210 = t57 * t60;
t209 = t60 * t41;
t208 = t60 * t230;
t182 = qJD(5) * t123;
t183 = qJD(5) * t119;
t8 = t119 * t112 + t113 * t182 + t123 * t21 + t183 * t60;
t206 = t8 * t119;
t204 = t119 * t20;
t201 = t60 * t113;
t200 = t121 * t122;
t199 = t121 * t125;
t128 = qJD(1) ^ 2;
t198 = t121 * t128;
t197 = t122 * t119;
t196 = t122 * t123;
t195 = t122 * t124;
t194 = t124 * t125;
t193 = t124 * t128;
t192 = t125 * t119;
t191 = t125 * t123;
t175 = pkin(1) * t186;
t83 = qJD(3) - t175;
t190 = -qJD(3) + t83;
t68 = t126 * qJD(2) * t121;
t114 = t121 ^ 2;
t189 = t114 - t115;
t188 = t114 + t223;
t180 = qJD(1) * qJD(2);
t164 = t121 * t180;
t181 = pkin(1) * t164 + qJDD(3);
t179 = qJD(1) * qJD(3);
t178 = qJDD(1) * t115;
t177 = qJDD(2) * t124;
t174 = g(1) * t194 + g(2) * t195 + g(3) * t121;
t171 = t64 * t183;
t170 = pkin(1) * t176;
t168 = t126 * t222;
t163 = t124 * t180;
t127 = qJD(2) ^ 2;
t61 = -t170 + t181;
t156 = -qJ(3) * t127 - t61;
t40 = pkin(2) * t164 - qJDD(1) * t74 + t181;
t13 = -t21 * pkin(4) + t40;
t25 = t113 * pkin(4) + t31;
t155 = qJD(5) * t25 - t13;
t152 = t123 * t57;
t151 = g(1) * t199 + g(2) * t200 - t216;
t150 = t121 * t163;
t147 = -t218 + t221;
t146 = t20 * t64 + t37 * t57;
t12 = t119 * t36 + t123 * t25;
t144 = t228 * t119 - t12 * t60 + t30 * t182;
t11 = -t119 * t25 + t123 * t36;
t143 = t11 * t60 + t123 * t231 + t30 * t183;
t142 = t123 * t20 + (t119 * t230 - t183) * t57;
t141 = -t162 + t102;
t140 = -pkin(4) * t20 + (-t57 - t230) * t30;
t139 = -t170 - t221;
t134 = -qJ(3) * qJDD(1) + (-qJD(3) - t83) * qJD(1);
t71 = -qJ(3) * t187 + t118;
t133 = -t71 * qJD(2) * t124 - t148;
t132 = t62 * t60 - t160 - t217 + t231;
t35 = -t222 * t65 - t203;
t130 = -t87 * t20 - t30 * t230 + (-t126 * t165 + t35) * t57;
t129 = qJ(3) ^ 2;
t98 = g(1) * t195;
t97 = g(2) * t199;
t53 = t110 * t191 + t197;
t52 = -t110 * t192 + t196;
t51 = -t110 * t196 + t192;
t50 = t110 * t197 + t191;
t46 = -t121 * t179 + t117 + (-t163 - t108) * qJ(3);
t34 = -t120 * t65 + t172;
t24 = -t37 * pkin(4) + t68;
t23 = -t230 ^ 2 + t60 ^ 2;
t17 = t45 * qJD(4) + t120 * t54 - t222 * t55;
t15 = -t201 - t22;
t14 = t21 - t202;
t10 = t123 * t13;
t4 = t152 * t57 + t204 + t212;
t3 = t142 - t209;
t2 = t152 * t44 + t206;
t1 = (t8 - t214) * t123 + (-t9 - t213) * t119;
t5 = [qJDD(1), t147, t148, t114 * qJDD(1) + 0.2e1 * t150, 0.2e1 * t121 * t176 - 0.2e1 * t180 * t189, qJDD(2) * t121 + t127 * t124, -t127 * t121 + t177, 0, -g(2) * t194 + t98, -g(1) * t200 + t97, pkin(1) * t178 + t98 + (t156 - t218) * t124 + (-qJ(3) * qJDD(2) + (-0.2e1 * t175 + t190) * qJD(2)) * t121, -qJ(3) * t177 + t97 + (t139 - t156) * t121 + (t189 * qJD(1) * pkin(1) + t190 * t124) * qJD(2), -t46 * t121 + t188 * t179 + (qJDD(1) * t188 - 0.2e1 * t150) * qJ(3) + t133, t129 * t178 + (t147 - t61) * t124 * pkin(1) + (t179 * t223 + t133) * qJ(3) + (-t46 * qJ(3) - t71 * qJD(3) + (pkin(1) * t83 - 0.2e1 * t129 * t186) * qJD(2)) * t121, t21 * t64 - t60 * t37, t136 * t21 - t64 * t22 + t230 * t37 + t60 * t38, t64 * t112 + t37 * t113, t112 * t136 - t38 * t113, 0, t110 * t147 + t112 * t138 - t17 * t113 - t136 * t40 - t74 * t22 - t230 * t68 + t62 * t38, -t109 * t147 - t45 * t112 - t16 * t113 - t74 * t21 + t62 * t37 + t40 * t64 - t68 * t60, -t44 * t171 + (t37 * t44 + t64 * t8) * t123, (-t119 * t44 - t123 * t41) * t37 + (-t206 - t123 * t9 + (t119 * t41 - t123 * t44) * qJD(5)) * t64, t123 * t146 - t136 * t8 - t171 * t57 + t44 * t38, -t182 * t57 * t64 - t119 * t146 + t136 * t9 - t41 * t38, -t136 * t20 + t57 * t38, -g(1) * t51 - g(2) * t53 - t10 * t136 + t11 * t38 + t17 * t41 - t138 * t9 + (t211 + t24 * t57 + (t136 * t25 - t45 * t57 + t215) * qJD(5)) * t123 + t224 * t119, -g(1) * t50 - g(2) * t52 - t12 * t38 + t17 * t44 - t138 * t8 + (-(-qJD(5) * t45 + t24) * t57 - t211 - t155 * t136 - qJD(5) * t215) * t119 + t224 * t123; 0, 0, 0, -t121 * t193, t189 * t128, t108, t176, qJDD(2), t151, t174, 0.2e1 * t117 + (pkin(1) * t193 + t134) * t121 + t151, -t114 * t128 * pkin(1) + t124 * t134 + t174, -pkin(1) * t108 + (qJ(3) * t198 + (t71 - t118) * qJD(1)) * t124, (qJ(3) * qJD(1) * t71 + t129 * t198) * t124 + (-t216 + t46 + (-qJD(1) * t83 + t148) * t121) * pkin(1), t208, t23, t14, t15, t112, t112 * t168 + t34 * t113 + t67 * t230 + (-t172 + (-t113 * t126 - t49) * t120) * qJD(4) + t132, t35 * t113 + t67 * t60 + (-t112 * t126 - t29) * t120 + (-t113 * t168 - t173) * qJD(4) - t149 + t225, t2, t1, t4, t3, t210, -t34 * t41 + (t41 * t184 - t222 * t9) * t126 + (-t228 - t226) * t123 + t130 * t119 + t143, -t34 * t44 + (t44 * t184 - t222 * t8) * t126 + t130 * t123 + (-t231 + t226) * t119 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t164 - t176, t108 + 0.2e1 * t163, (-t114 - t115) * t128, -t115 * t128 * qJ(3) + t187 * t71 + t139 + t181 + t218, 0, 0, 0, 0, 0, t22 - t201, t21 + t202, 0, 0, 0, 0, 0, t142 + t209, -t123 * t157 - t204 + t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t23, t14, t15, t112, t132 + (-qJD(4) + t113) * t31, -t30 * t113 - t131 + t225, t2, t1, t4, t3, t210, -t31 * t41 + t140 * t119 + (-t228 - t229) * t123 + t143, -t31 * t44 + t140 * t123 + (-t231 + t229) * t119 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t41, -t41 ^ 2 + t44 ^ 2, t8 + t214, t213 - t9, t20, -g(1) * t52 + g(2) * t50 + t119 * t141 + t12 * t57 - t182 * t25 - t30 * t44 + t10, g(1) * t53 - g(2) * t51 + t11 * t57 + t119 * t155 + t123 * t141 + t30 * t41;];
tau_reg = t5;
