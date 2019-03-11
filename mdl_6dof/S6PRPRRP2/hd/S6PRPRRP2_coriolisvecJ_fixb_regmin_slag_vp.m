% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:11
% EndTime: 2019-03-08 20:03:16
% DurationCPUTime: 2.05s
% Computational Cost: add. (2559->303), mult. (6357->411), div. (0->0), fcn. (4737->10), ass. (0->164)
t102 = sin(pkin(11));
t111 = cos(qJ(2));
t103 = sin(pkin(6));
t169 = qJD(1) * t103;
t149 = t111 * t169;
t104 = cos(pkin(11));
t108 = sin(qJ(2));
t150 = t108 * t169;
t88 = t104 * t150;
t61 = t102 * t149 + t88;
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t135 = pkin(4) * t107 - pkin(9) * t110;
t84 = t135 * qJD(4);
t212 = t61 - t84;
t167 = qJD(2) * t110;
t94 = -qJD(5) + t167;
t160 = qJD(2) * qJD(4);
t142 = t107 * t160;
t139 = pkin(5) * t142;
t106 = sin(qJ(5));
t109 = cos(qJ(5));
t86 = qJD(2) * pkin(2) + t149;
t55 = t102 * t86 + t88;
t52 = qJD(2) * pkin(8) + t55;
t105 = cos(pkin(6));
t93 = qJD(1) * t105 + qJD(3);
t208 = -t107 * t52 + t110 * t93;
t66 = (t102 * t108 - t104 * t111) * t103;
t63 = qJD(2) * t66;
t57 = qJD(1) * t63;
t13 = qJD(4) * t208 - t110 * t57;
t162 = qJD(5) * t109;
t163 = qJD(5) * t106;
t37 = t107 * t93 + t110 * t52;
t32 = qJD(4) * pkin(9) + t37;
t126 = -pkin(4) * t110 - pkin(9) * t107 - pkin(3);
t87 = t102 * t150;
t54 = t104 * t86 - t87;
t40 = qJD(2) * t126 - t54;
t67 = (t102 * t111 + t104 * t108) * t103;
t41 = (qJD(1) * t67 + t84) * qJD(2);
t141 = -t106 * t13 + t109 * t41 - t32 * t162 - t40 * t163;
t2 = -t139 - t141;
t9 = t106 * t40 + t109 * t32;
t7 = -qJ(6) * t94 + t9;
t211 = t7 * t94 + t2;
t174 = t106 * t110;
t64 = t104 * t149 - t87;
t210 = t212 * t109 - t64 * t174;
t173 = t109 * t110;
t205 = pkin(2) * t104;
t76 = t126 - t205;
t209 = -t212 * t106 + t76 * t162 - t64 * t173;
t146 = t107 * t162;
t159 = qJD(4) * qJD(5);
t164 = qJD(4) * t110;
t59 = qJD(2) * (t106 * t164 + t146) + t106 * t159;
t166 = qJD(4) * t106;
t168 = qJD(2) * t107;
t80 = t109 * t168 + t166;
t207 = t80 ^ 2;
t165 = qJD(4) * t107;
t14 = -t107 * t57 + t52 * t164 + t93 * t165;
t144 = t109 * t160;
t147 = t107 * t163;
t58 = qJD(2) * t147 - t109 * t159 - t110 * t144;
t5 = pkin(5) * t59 + qJ(6) * t58 - qJD(6) * t80 + t14;
t204 = t106 * t5;
t203 = t109 * t5;
t31 = -qJD(4) * pkin(4) - t208;
t161 = t109 * qJD(4);
t78 = t106 * t168 - t161;
t12 = pkin(5) * t78 - qJ(6) * t80 + t31;
t202 = t12 * t80;
t201 = t64 * t78;
t200 = t64 * t80;
t199 = t78 * t94;
t198 = t80 * t78;
t197 = t80 * t94;
t96 = pkin(2) * t102 + pkin(8);
t183 = t109 * t96;
t196 = (-t96 * t163 - qJD(6)) * t110 + (qJ(6) - t183) * t165 + t209;
t188 = t106 * t96;
t145 = pkin(5) + t188;
t190 = t106 * t76 + t96 * t173;
t195 = t190 * qJD(5) - t145 * t165 + t210;
t131 = pkin(5) * t106 - qJ(6) * t109;
t194 = t106 * qJD(6) + t94 * t131 + t37;
t83 = t135 * qJD(2);
t193 = t106 * t83 + t109 * t208;
t148 = t110 * t161;
t192 = -t109 * t107 * t59 - t78 * t148;
t189 = t106 * t31;
t187 = t109 * t31;
t186 = t109 * t76;
t184 = t109 * t94;
t182 = t110 * t58;
t181 = t110 * t59;
t179 = t14 * t106;
t178 = t14 * t109;
t8 = -t106 * t32 + t109 * t40;
t177 = qJD(6) - t8;
t176 = qJD(5) * t78;
t113 = qJD(2) ^ 2;
t175 = t103 * t113;
t112 = qJD(4) ^ 2;
t172 = t112 * t107;
t171 = t112 * t110;
t100 = t107 ^ 2;
t170 = -t110 ^ 2 + t100;
t158 = pkin(9) * t106 * t94;
t157 = pkin(9) * t184;
t156 = t106 * t41 + t109 * t13 + t40 * t162;
t155 = pkin(9) * t165;
t154 = pkin(9) * t161;
t153 = t80 * t164;
t152 = t94 * t163;
t151 = t94 * t162;
t140 = -t58 + t176;
t138 = t80 * t146;
t137 = t94 * t146;
t136 = qJ(6) * t142;
t6 = pkin(5) * t94 + t177;
t134 = -t106 * t7 + t109 * t6;
t133 = t106 * t6 + t109 * t7;
t132 = pkin(5) * t109 + qJ(6) * t106;
t130 = -t106 * t208 + t109 * t83;
t50 = t105 * t107 + t110 * t67;
t24 = t106 * t66 + t109 * t50;
t23 = t106 * t50 - t109 * t66;
t128 = t105 * t110 - t107 * t67;
t91 = t100 * t144;
t125 = -t94 * t147 - t91;
t124 = -t9 * t94 + t141;
t123 = t131 + t96;
t62 = qJD(2) * t67;
t56 = qJD(1) * t62;
t122 = qJD(2) * t61 - t112 * t96 - t56;
t51 = -qJD(2) * pkin(3) - t54;
t121 = qJD(4) * (qJD(2) * (-pkin(3) - t205) + t51 + t64);
t120 = t32 * t163 - t156;
t119 = (-qJD(2) * t100 + t110 * t94) * t106;
t21 = qJD(4) * t50 - t63 * t107;
t22 = qJD(4) * t128 - t63 * t110;
t3 = qJD(5) * t24 + t22 * t106 - t62 * t109;
t117 = -t128 * t59 - t142 * t23 + t21 * t78 + t3 * t94;
t1 = -qJD(6) * t94 - t120 + t136;
t116 = qJD(5) * t134 + t1 * t109 + t106 * t2;
t4 = -qJD(5) * t23 + t62 * t106 + t22 * t109;
t115 = -t128 * t58 + t142 * t24 - t21 * t80 - t4 * t94;
t114 = qJD(4) * t119 + t78 * t165 + t137 - t181;
t90 = -pkin(4) - t132;
t77 = t94 * t148;
t70 = t80 * t165;
t60 = t123 * t107;
t48 = pkin(5) * t80 + qJ(6) * t78;
t43 = t110 * t145 - t186;
t42 = -qJ(6) * t110 + t190;
t35 = -t58 - t199;
t28 = (qJD(5) * t132 - qJD(6) * t109) * t107 + t123 * t164;
t18 = -pkin(5) * t168 - t130;
t17 = qJ(6) * t168 + t193;
t10 = [0, 0, -t108 * t175, -t111 * t175, -t54 * t62 - t55 * t63 + t56 * t66 - t57 * t67, 0, 0, 0, 0, 0, -t21 * qJD(4) + (-t110 * t62 + t66 * t165) * qJD(2), -t22 * qJD(4) + (t107 * t62 + t66 * t164) * qJD(2), 0, 0, 0, 0, 0, t117, -t115, t117, -t23 * t58 - t24 * t59 + t3 * t80 - t4 * t78, t115, t1 * t24 + t12 * t21 - t128 * t5 + t2 * t23 + t3 * t6 + t4 * t7; 0, 0, 0, 0, t54 * t61 - t55 * t64 + (-t102 * t57 - t104 * t56) * pkin(2), 0.2e1 * t110 * t142, -0.2e1 * t170 * t160, t171, -t172, 0, t107 * t121 + t110 * t122, -t107 * t122 + t110 * t121, t80 * t148 + (-t109 * t58 - t80 * t163) * t107, -t138 + (-t153 + (t58 + t176) * t107) * t106 + t192, -t125 + t70 - t77 + t182, t137 + t181 + (-t107 * t78 + t119) * qJD(4) (-t94 - t167) * t165 (t76 * t163 + t210) * t94 + (t96 * t151 + (t78 * t96 + t189) * qJD(4) - t141) * t110 + (t31 * t162 + t179 + t96 * t59 - t201 + (-t94 * t188 + (-t96 * t174 + t186) * qJD(2) + t8) * qJD(4)) * t107, t209 * t94 + ((-t94 * t96 - t32) * t163 + (t80 * t96 + t187) * qJD(4) + t156) * t110 + (-t31 * t163 + t178 - t96 * t58 - t200 + (-t190 * qJD(2) - t94 * t183 - t9) * qJD(4)) * t107, t28 * t78 + t59 * t60 + t195 * t94 + (t12 * t166 + t2) * t110 + (t12 * t162 + t204 - t201 + (-qJD(2) * t43 - t6) * qJD(4)) * t107, -t42 * t59 - t43 * t58 + t195 * t80 - t196 * t78 + t134 * t164 + (-qJD(5) * t133 - t1 * t106 + t109 * t2) * t107, -t28 * t80 + t58 * t60 - t196 * t94 + (-t12 * t161 - t1) * t110 + (t12 * t163 - t203 + t200 + (qJD(2) * t42 + t7) * qJD(4)) * t107, t1 * t42 + t2 * t43 + t5 * t60 + t196 * t7 + t195 * t6 + (-t107 * t64 + t28) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t171, 0, 0, 0, 0, 0, t114, t70 + (t94 * t161 + t58) * t110 + t125, t114, t138 + (t107 * t140 + t153) * t106 + t192, -t182 - t77 + t91 + (-qJD(4) * t80 + t152) * t107 (qJD(4) * t133 - t5) * t110 + (qJD(4) * t12 + t116) * t107; 0, 0, 0, 0, 0, -t107 * t113 * t110, t170 * t113, 0, 0, 0, qJD(4) * t37 - t51 * t168 - t14 (-qJD(2) * t51 + t57) * t110, -t58 * t106 - t80 * t184 (-t58 + t199) * t109 + (-t59 + t197) * t106, -t151 + (t94 * t173 + (-t80 + t166) * t107) * qJD(2), t152 + (-t94 * t174 + (t78 + t161) * t107) * qJD(2), t94 * t168, -pkin(4) * t59 - t178 + t130 * t94 - t37 * t78 + (t157 + t189) * qJD(5) + (-t8 * t107 + (-t110 * t31 - t155) * t106) * qJD(2), pkin(4) * t58 + t179 - t193 * t94 - t37 * t80 + (-t158 + t187) * qJD(5) + (-t31 * t173 + (t9 - t154) * t107) * qJD(2), -t203 - t18 * t94 + t59 * t90 - t194 * t78 + (t106 * t12 + t157) * qJD(5) + (t107 * t6 + (-t110 * t12 - t155) * t106) * qJD(2), t17 * t78 - t18 * t80 + (t1 - t94 * t6 + (qJD(5) * t80 - t59) * pkin(9)) * t109 + (pkin(9) * t140 + t211) * t106, -t204 + t17 * t94 + t58 * t90 + t194 * t80 + (-t109 * t12 + t158) * qJD(5) + (t12 * t173 + (-t7 + t154) * t107) * qJD(2), t116 * pkin(9) - t194 * t12 - t17 * t7 - t18 * t6 + t5 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, -t78 ^ 2 + t207, t35, -t197 - t59, t142, -t31 * t80 + t124, t31 * t78 - t8 * t94 + t120, -t48 * t78 + t124 + 0.2e1 * t139 - t202, pkin(5) * t58 - qJ(6) * t59 + (t7 - t9) * t80 + (t6 - t177) * t78, 0.2e1 * t136 - t12 * t78 + t48 * t80 + (-0.2e1 * qJD(6) + t8) * t94 - t120, -pkin(5) * t2 + qJ(6) * t1 - t12 * t48 + t177 * t7 - t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142 + t198, t35, -t94 ^ 2 - t207, t202 + t211;];
tauc_reg  = t10;
