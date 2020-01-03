% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:24
% EndTime: 2019-12-31 21:09:31
% DurationCPUTime: 2.09s
% Computational Cost: add. (2045->339), mult. (4940->442), div. (0->0), fcn. (2860->4), ass. (0->166)
t195 = pkin(3) + qJ(5);
t105 = sin(qJ(2));
t159 = qJD(1) * qJD(2);
t92 = t105 * t159;
t214 = t195 * t92;
t104 = sin(qJ(3));
t106 = cos(qJ(3));
t161 = qJD(3) * t106;
t149 = t105 * t161;
t107 = cos(qJ(2));
t163 = qJD(2) * t107;
t213 = t104 * t163 + t149;
t69 = -pkin(2) * t107 - pkin(7) * t105 - pkin(1);
t51 = t69 * qJD(1);
t166 = qJD(1) * t107;
t95 = pkin(6) * t166;
t73 = qJD(2) * pkin(7) + t95;
t28 = t104 * t73 - t106 * t51;
t165 = qJD(2) * t104;
t167 = qJD(1) * t105;
t61 = t106 * t167 + t165;
t124 = t61 * pkin(4) + t28;
t170 = qJD(4) + t124;
t212 = -0.2e1 * t159;
t57 = t61 ^ 2;
t87 = -qJD(3) + t166;
t84 = t87 ^ 2;
t211 = -t57 - t84;
t76 = qJD(4) * t87;
t85 = qJ(4) * t92;
t210 = t85 - t76;
t29 = t104 * t51 + t106 * t73;
t22 = qJ(4) * t87 - t29;
t134 = pkin(3) * t92;
t135 = pkin(6) * t92;
t162 = qJD(3) * t104;
t129 = pkin(2) * t105 - pkin(7) * t107;
t67 = t129 * qJD(2);
t52 = qJD(1) * t67;
t141 = -t104 * t135 - t106 * t52 + t73 * t161 + t51 * t162;
t9 = -t134 + t141;
t209 = -t22 * t87 + t9;
t208 = t104 * qJD(4) + t95 + (t104 * t166 - t162) * pkin(3);
t158 = qJD(2) * qJD(3);
t38 = t213 * qJD(1) + t104 * t158;
t160 = t106 * qJD(2);
t59 = t104 * t167 - t160;
t207 = t59 ^ 2;
t206 = pkin(4) + pkin(7);
t205 = pkin(4) * t59;
t10 = t195 * t87 + t170;
t204 = t10 * t87;
t72 = -qJD(2) * pkin(2) + pkin(6) * t167;
t121 = -t61 * qJ(4) + t72;
t13 = t195 * t59 + t121;
t203 = t13 * t61;
t24 = pkin(3) * t59 + t121;
t201 = t24 * t61;
t200 = t59 * t87;
t147 = t107 * t159;
t151 = t105 * t162;
t37 = qJD(1) * t151 + (-t147 - t158) * t106;
t117 = pkin(6) * t147 + t37 * qJ(4) - t61 * qJD(4);
t6 = pkin(3) * t38 + t117;
t199 = t6 * t104;
t198 = t6 * t106;
t197 = t61 * t59;
t196 = t87 * t61;
t153 = -pkin(6) * t104 - pkin(3);
t175 = t106 * t107;
t112 = pkin(4) * t175 + (-qJ(5) + t153) * t105;
t64 = t129 * qJD(1);
t184 = t106 * t64;
t75 = t206 * t106;
t194 = -qJD(1) * t112 + qJD(3) * t75 + t184;
t144 = pkin(6) * t106 - qJ(4);
t177 = t104 * t107;
t113 = -pkin(4) * t177 - t144 * t105;
t48 = t104 * t64;
t193 = qJD(1) * t113 + t206 * t162 + t48;
t179 = qJ(4) * t106;
t126 = qJ(5) * t104 - t179;
t116 = t126 * t107;
t192 = qJD(1) * t116 - t126 * qJD(3) + t106 * qJD(5) + t208;
t191 = qJ(4) * t161 - t166 * t179 + t208;
t190 = t104 * t67 + t69 * t161;
t178 = t104 * t105;
t188 = pkin(3) * t178 + t105 * pkin(6);
t91 = pkin(6) * t175;
t187 = t104 * t69 + t91;
t186 = qJ(4) * t38;
t185 = t104 * t72;
t183 = t106 * t72;
t182 = t106 * t87;
t181 = t37 * t104;
t180 = t59 * qJ(4);
t176 = t105 * t106;
t109 = qJD(1) ^ 2;
t174 = t107 * t109;
t108 = qJD(2) ^ 2;
t173 = t108 * t105;
t172 = t108 * t107;
t171 = -qJD(4) - t28;
t17 = t29 - t205;
t169 = -qJD(5) - t17;
t101 = t105 ^ 2;
t168 = -t107 ^ 2 + t101;
t164 = qJD(2) * t105;
t157 = pkin(7) * t104 * t87;
t156 = pkin(7) * t182;
t90 = pkin(6) * t177;
t155 = pkin(7) * t164;
t154 = pkin(7) * t160;
t150 = t107 * t162;
t148 = -t92 + t197;
t145 = -t104 * qJ(4) - pkin(2);
t143 = t106 * t69 - t90;
t142 = -t104 * t52 + t106 * t135 - t51 * t161 + t73 * t162;
t140 = t213 * pkin(3) + pkin(6) * t163 + qJ(4) * t151;
t139 = pkin(1) * t212;
t138 = -t61 + t165;
t137 = t59 + t160;
t133 = -qJD(3) * t91 + t106 * t67 - t69 * t162;
t132 = t76 + t142;
t35 = qJ(4) * t107 - t187;
t131 = -t107 * qJD(4) + t190;
t130 = t153 * t105;
t128 = pkin(6) * (-t87 + t166);
t21 = pkin(3) * t87 - t171;
t127 = t104 * t22 + t106 * t21;
t125 = qJD(1) * t101 - t107 * t87;
t123 = -t37 * pkin(4) + t141;
t122 = -pkin(4) * t38 - t142;
t120 = -t29 * t87 - t141;
t1 = t59 * qJD(5) + t195 * t38 + t117;
t119 = t1 * t104 + t13 * t161;
t118 = -t1 * t106 + t13 * t162;
t114 = -t13 * t59 + t122;
t18 = -t37 - t200;
t111 = t123 - t214;
t110 = -t196 - t38;
t100 = t107 * pkin(3);
t83 = 0.2e1 * t85;
t74 = t206 * t104;
t68 = -pkin(3) * t106 + t145;
t50 = -t195 * t106 + t145;
t42 = -qJ(4) * t176 + t188;
t36 = t100 - t143;
t34 = t126 * t105 + t188;
t32 = pkin(3) * t61 + t180;
t31 = qJD(1) * t130 - t184;
t30 = t144 * t167 - t48;
t26 = -pkin(4) * t178 - t35;
t25 = t107 * qJ(5) + t100 + t90 + (pkin(4) * t105 - t69) * t106;
t19 = t195 * t61 + t180;
t15 = (-qJ(4) * t163 - qJD(4) * t105) * t106 + t140;
t14 = qJD(2) * t130 - t133;
t12 = qJD(5) - t22 - t205;
t11 = -qJ(4) * t164 + (t105 * t160 + t150) * pkin(6) - t131;
t8 = qJD(2) * t116 + (qJD(5) * t104 + (qJ(5) * qJD(3) - qJD(4)) * t106) * t105 + t140;
t7 = (-pkin(4) * t176 - t90) * qJD(3) + t113 * qJD(2) + t131;
t5 = -t85 + t132;
t4 = -pkin(4) * t151 + qJD(2) * t112 + t107 * qJD(5) - t133;
t3 = t122 + t210;
t2 = qJD(5) * t87 + t111;
t16 = [0, 0, 0, 0.2e1 * t107 * t92, t168 * t212, t172, -t173, 0, -pkin(6) * t172 + t105 * t139, pkin(6) * t173 + t107 * t139, t61 * t107 * t160 + (-t37 * t106 - t61 * t162) * t105, (-t104 * t61 - t106 * t59) * t163 + (t181 - t106 * t38 + (t104 * t59 - t106 * t61) * qJD(3)) * t105, t87 * t151 + t37 * t107 + (t105 * t61 + t106 * t125) * qJD(2), t87 * t149 + t38 * t107 + (-t104 * t125 - t105 * t59) * qJD(2), (-t87 - t166) * t164, -t133 * t87 + t141 * t107 + (pkin(6) * t38 + t72 * t161) * t105 + ((pkin(6) * t59 + t185) * t107 + (qJD(1) * t143 + t104 * t128 - t28) * t105) * qJD(2), (-pkin(6) * t150 + t190) * t87 - t142 * t107 + (-pkin(6) * t37 - t162 * t72) * t105 + ((pkin(6) * t61 + t183) * t107 + (-qJD(1) * t187 + t106 * t128 - t29) * t105) * qJD(2), t11 * t59 + t14 * t61 + t35 * t38 - t36 * t37 + t127 * t163 + (t104 * t5 + t106 * t9 + (-t104 * t21 + t106 * t22) * qJD(3)) * t105, -t14 * t87 - t15 * t59 - t42 * t38 + (-t165 * t24 - t9) * t107 + (-t24 * t161 - t199 + (qJD(1) * t36 + t21) * qJD(2)) * t105, t11 * t87 - t15 * t61 + t42 * t37 + (-t160 * t24 + t5) * t107 + (t24 * t162 - t198 + (-qJD(1) * t35 - t22) * qJD(2)) * t105, t11 * t22 + t14 * t21 + t15 * t24 + t35 * t5 + t36 * t9 + t42 * t6, -t25 * t37 - t26 * t38 + t4 * t61 - t59 * t7 + (t10 * t106 - t104 * t12) * t163 + (-t104 * t3 + t106 * t2 + (-t10 * t104 - t106 * t12) * qJD(3)) * t105, t34 * t37 - t8 * t61 - t7 * t87 + (-t13 * t160 - t3) * t107 + ((qJD(1) * t26 + t12) * qJD(2) + t118) * t105, t34 * t38 + t4 * t87 + t8 * t59 + (t13 * t165 + t2) * t107 + ((-qJD(1) * t25 - t10) * qJD(2) + t119) * t105, t1 * t34 + t10 * t4 + t12 * t7 + t13 * t8 + t2 * t25 + t26 * t3; 0, 0, 0, -t105 * t174, t168 * t109, 0, 0, 0, t109 * pkin(1) * t105, pkin(1) * t174, -t61 * t182 - t181, (-t37 + t200) * t106 + (-t38 + t196) * t104, -t87 * t161 + (t105 * t138 + t87 * t175) * qJD(1), t87 * t162 + (t105 * t137 - t87 * t177) * qJD(1), t87 * t167, t64 * t182 - pkin(2) * t38 + (t156 + t185) * qJD(3) + (t28 * t105 + (-t107 * t72 - t155) * t104 + (-t107 * t137 + t178 * t87) * pkin(6)) * qJD(1), pkin(2) * t37 - t48 * t87 + (-t157 + t183) * qJD(3) + (-t72 * t175 + (t29 - t154) * t105 + (t107 * t138 + t176 * t87) * pkin(6)) * qJD(1), -t30 * t59 - t31 * t61 + (-t5 - t87 * t21 + (qJD(3) * t61 - t38) * pkin(7)) * t106 + ((qJD(3) * t59 - t37) * pkin(7) + t209) * t104, t198 + t31 * t87 - t68 * t38 + t191 * t59 + (-t104 * t24 - t156) * qJD(3) + (-t105 * t21 + (t107 * t24 + t155) * t104) * qJD(1), -t199 - t30 * t87 + t68 * t37 + t191 * t61 + (-t106 * t24 + t157) * qJD(3) + (t24 * t175 + (t22 + t154) * t105) * qJD(1), -t21 * t31 - t22 * t30 + t6 * t68 - t191 * t24 + (qJD(3) * t127 + t9 * t104 - t5 * t106) * pkin(7), -t37 * t74 - t38 * t75 + t194 * t61 + t193 * t59 + (t3 - t204) * t106 + (t12 * t87 + t2) * t104, t37 * t50 + t193 * t87 + t192 * t61 + (t13 * t175 + (qJD(2) * t75 - t12) * t105) * qJD(1) - t119, t38 * t50 + t194 * t87 - t192 * t59 + (-t13 * t177 + (-qJD(2) * t74 + t10) * t105) * qJD(1) + t118, t1 * t50 + t10 * t194 - t12 * t193 - t13 * t192 + t2 * t74 + t3 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t57 - t207, t18, t110, t92, -t61 * t72 + t120, t28 * t87 + t59 * t72 + t142, pkin(3) * t37 - t186 + (-t22 - t29) * t61 + (t21 + t171) * t59, t32 * t59 - t120 - 0.2e1 * t134 + t201, t171 * t87 - t24 * t59 + t32 * t61 - t132 + t83, -t9 * pkin(3) - t5 * qJ(4) + t171 * t22 - t21 * t29 - t24 * t32, -t186 + t195 * t37 + (t12 + t169) * t61 + (t10 - t170) * t59, -t124 * t87 + t19 * t61 + t114 - 0.2e1 * t76 + t83, -t203 - t19 * t59 + (-0.2e1 * qJD(5) - t17) * t87 + 0.2e1 * t214 - t123, t3 * qJ(4) + t10 * t169 + t12 * t170 - t13 * t19 - t195 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t148, t211, t201 + t209, t18, t211, t148, t203 + (qJD(5) + t12) * t87 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t92 + t197, -t84 - t207, t114 - t204 + t210;];
tauc_reg = t16;
