% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:51
% EndTime: 2019-03-08 20:48:57
% DurationCPUTime: 2.38s
% Computational Cost: add. (2076->329), mult. (4848->484), div. (0->0), fcn. (3525->10), ass. (0->187)
t116 = cos(pkin(6));
t123 = cos(qJ(4));
t119 = sin(qJ(4));
t195 = qJD(1) * t119;
t125 = -pkin(2) - pkin(8);
t124 = cos(qJ(2));
t115 = sin(pkin(6));
t196 = qJD(1) * t115;
t169 = t124 * t196;
t146 = qJD(3) - t169;
t73 = t125 * qJD(2) + t146;
t50 = -t116 * t195 + t123 * t73;
t37 = -qJD(4) * pkin(4) - t50;
t118 = sin(qJ(5));
t122 = cos(qJ(5));
t181 = t122 * qJD(4);
t191 = qJD(2) * t123;
t83 = t118 * t191 - t181;
t24 = pkin(5) * t83 + t37;
t117 = sin(qJ(6));
t121 = cos(qJ(6));
t166 = t122 * t191;
t189 = qJD(4) * t118;
t85 = t166 + t189;
t30 = t117 * t85 + t121 * t83;
t244 = t24 * t30;
t142 = t117 * t83 - t121 * t85;
t243 = t142 * t30;
t234 = qJD(5) + qJD(6);
t86 = t117 * t118 - t121 * t122;
t242 = t234 * t86;
t87 = t117 * t122 + t118 * t121;
t128 = t234 * t87;
t241 = t142 ^ 2 - t30 ^ 2;
t192 = qJD(2) * t119;
t107 = qJD(5) + t192;
t104 = qJD(6) + t107;
t182 = qJD(6) * t121;
t183 = qJD(6) * t117;
t184 = qJD(5) * t123;
t165 = t118 * t184;
t235 = -t119 * t181 - t165;
t54 = t235 * qJD(2) + qJD(5) * t181;
t162 = t118 * t192;
t55 = -qJD(4) * t162 + qJD(5) * t85;
t6 = -t117 * t55 + t121 * t54 - t83 * t182 - t85 * t183;
t240 = t104 * t30 + t6;
t180 = qJD(2) * qJD(4);
t108 = t123 * t180;
t120 = sin(qJ(2));
t193 = qJD(2) * t115;
t168 = t120 * t193;
t187 = qJD(4) * t123;
t25 = t73 * t187 + (-qJD(4) * t116 + t168) * t195;
t148 = pkin(4) * t123 + pkin(9) * t119;
t79 = qJD(4) * t148 + qJD(3);
t52 = (t79 + t169) * qJD(2);
t157 = -t118 * t25 + t122 * t52;
t170 = t120 * t196;
t92 = pkin(4) * t119 - pkin(9) * t123 + qJ(3);
t59 = qJD(2) * t92 + t170;
t220 = t118 * t59;
t206 = t116 * t123;
t106 = qJD(1) * t206;
t51 = t119 * t73 + t106;
t38 = qJD(4) * pkin(9) + t51;
t18 = t122 * t38 + t220;
t131 = -qJD(5) * t18 + t157;
t2 = pkin(5) * t108 - pkin(10) * t54 + t131;
t185 = qJD(5) * t122;
t179 = -t118 * t52 - t122 * t25 - t59 * t185;
t186 = qJD(5) * t118;
t139 = -t38 * t186 - t179;
t5 = -pkin(10) * t55 + t139;
t173 = -t117 * t5 + t121 * t2;
t11 = -pkin(10) * t83 + t18;
t223 = t11 * t121;
t17 = -t118 * t38 + t122 * t59;
t10 = -pkin(10) * t85 + t17;
t8 = pkin(5) * t107 + t10;
t4 = t117 * t8 + t223;
t239 = -qJD(6) * t4 + t24 * t142 + t173;
t130 = qJD(6) * t142 - t117 * t54 - t121 * t55;
t238 = -t104 * t142 + t130;
t205 = t118 * t120;
t237 = t122 * t79 - (-t119 * t205 + t122 * t124) * t196;
t200 = t120 * t122;
t236 = (t118 * t124 + t119 * t200) * t196 - t123 * t125 * t181 - t118 * t79 - t92 * t185;
t160 = qJD(6) * t8 + t5;
t9 = t11 * t183;
t233 = t117 * t2 + t121 * t160 - t9;
t232 = pkin(9) + pkin(10);
t201 = t119 * t125;
t103 = t122 * t201;
t203 = t118 * t125;
t159 = pkin(5) - t203;
t202 = t119 * t122;
t177 = pkin(10) * t202;
t231 = (-t103 + (pkin(10) * t123 - t92) * t118) * qJD(5) + (t159 * t123 + t177) * qJD(4) + t237;
t188 = qJD(4) * t119;
t163 = t118 * t188;
t164 = t122 * t184;
t133 = t163 - t164;
t174 = t118 * t201;
t230 = -pkin(10) * t133 + qJD(5) * t174 + t236;
t88 = t148 * qJD(2);
t229 = t118 * t88 + t122 * t50;
t137 = t86 * t119;
t228 = -qJD(2) * t137 - t242;
t135 = t87 * qJD(2);
t227 = t119 * t135 + t128;
t226 = qJD(2) * pkin(2);
t225 = t107 * t83;
t224 = t107 * t85;
t222 = t118 * t37;
t221 = t118 * t54;
t219 = t122 * t37;
t218 = t123 * t54;
t194 = qJD(2) * qJ(3);
t91 = t170 + t194;
t216 = t124 * t91;
t152 = t123 * t168;
t213 = -qJD(4) * t106 - t73 * t188;
t26 = -qJD(1) * t152 - t213;
t215 = t26 * t118;
t214 = t26 * t122;
t212 = t118 * t92 + t103;
t211 = t107 * t118;
t210 = t107 * t119;
t209 = t107 * t122;
t208 = t115 * t124;
t127 = qJD(2) ^ 2;
t207 = t115 * t127;
t204 = t118 * t123;
t199 = t122 * t123;
t114 = t123 ^ 2;
t198 = t119 ^ 2 - t114;
t126 = qJD(4) ^ 2;
t197 = -t126 - t127;
t190 = qJD(4) * t104;
t176 = t120 * t207;
t175 = t124 * t207;
t171 = qJD(5) * t232;
t167 = t124 * t193;
t156 = -t118 * t50 + t122 * t88;
t155 = t107 * t125 + t38;
t154 = qJD(5) * t119 + qJD(2);
t153 = t119 * t168;
t151 = -t51 + (t162 + t186) * pkin(5);
t99 = t232 * t122;
t150 = qJD(6) * t99 + (pkin(5) * t123 + t177) * qJD(2) + t156 + t122 * t171;
t98 = t232 * t118;
t149 = pkin(10) * t162 + qJD(6) * t98 + t118 * t171 + t229;
t147 = -t91 + t170;
t78 = t122 * t92;
t28 = -pkin(10) * t199 + t159 * t119 + t78;
t41 = -pkin(10) * t204 + t212;
t145 = t117 * t28 + t121 * t41;
t72 = -t119 * t208 + t206;
t44 = t115 * t200 - t118 * t72;
t45 = t115 * t205 + t122 * t72;
t144 = -t117 * t45 + t121 * t44;
t143 = t117 * t44 + t121 * t45;
t141 = qJD(2) * t114 - t210;
t140 = -pkin(9) * t187 + t119 * t37;
t71 = t116 * t119 + t123 * t208;
t138 = t147 * qJD(2);
t136 = qJD(2) * t86;
t134 = t147 - t194;
t80 = (qJD(3) + t169) * qJD(2);
t129 = qJD(2) * t146 - t125 * t126 + t80;
t111 = -pkin(5) * t122 - pkin(4);
t102 = t119 * t108;
t82 = t146 - t226;
t81 = (pkin(5) * t118 - t125) * t123;
t65 = t86 * t123;
t64 = t87 * t123;
t56 = -pkin(5) * t133 + t125 * t188;
t43 = qJD(4) * t72 - t152;
t42 = -qJD(4) * t71 + t153;
t20 = -t183 * t204 + (t234 * t199 - t163) * t121 + t235 * t117;
t19 = qJD(4) * t137 - t123 * t128;
t16 = pkin(5) * t55 + t26;
t13 = qJD(5) * t44 + t118 * t167 + t122 * t42;
t12 = -qJD(5) * t45 - t118 * t42 + t122 * t167;
t3 = -t11 * t117 + t121 * t8;
t1 = [0, 0, -t176, -t175, t176, t175 (t120 * t80 + (t216 + (t82 - t169) * t120) * qJD(2)) * t115, 0, 0, 0, 0, 0, t119 * t175 + (-t43 + t152) * qJD(4), t123 * t175 + (-t42 - t153) * qJD(4), 0, 0, 0, 0, 0, t107 * t12 + t108 * t44 + t43 * t83 + t55 * t71, -t107 * t13 - t108 * t45 + t43 * t85 + t54 * t71, 0, 0, 0, 0, 0 (-qJD(6) * t143 - t117 * t13 + t12 * t121) * t104 + t144 * t108 + t43 * t30 - t71 * t130 -(qJD(6) * t144 + t117 * t12 + t121 * t13) * t104 - t143 * t108 - t43 * t142 + t71 * t6; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), qJ(3) * t80 + qJD(3) * t91 + (-t216 + (-t82 - t226) * t120) * t196, -0.2e1 * t102, 0.2e1 * t198 * t180, -t126 * t119, -t126 * t123, 0, t119 * t129 - t134 * t187, t123 * t129 + t134 * t188, -t85 * t165 + (-t188 * t85 + t218) * t122 (t118 * t85 + t122 * t83) * t188 + (-t221 - t122 * t55 + (t118 * t83 - t122 * t85) * qJD(5)) * t123, -t107 * t165 + t119 * t54 + (t122 * t141 + t123 * t85) * qJD(4), -t107 * t164 - t119 * t55 + (-t118 * t141 - t123 * t83) * qJD(4), t107 * t187 + t102 (-t186 * t92 + t237) * t107 + ((t125 * t83 - t222) * qJD(4) + (-t122 * t155 - t220) * qJD(5) + t157) * t119 + (t83 * t170 + t37 * t185 + t215 - t125 * t55 + (-t107 * t203 + (t78 - t174) * qJD(2) + t17) * qJD(4)) * t123, t236 * t107 + (t155 * t186 + (t125 * t85 - t219) * qJD(4) + t179) * t119 + (t85 * t170 - t37 * t186 + t214 - t125 * t54 + (-qJD(2) * t212 - t18) * qJD(4)) * t123, -t142 * t19 - t6 * t65, -t130 * t65 + t142 * t20 - t19 * t30 - t6 * t64, t104 * t19 + t119 * t6 + (-qJD(2) * t65 - t142) * t187, -t104 * t20 + t119 * t130 + (-qJD(2) * t64 - t30) * t187, t104 * t187 + t102, t173 * t119 + t56 * t30 - t81 * t130 + t16 * t64 + t24 * t20 + (t230 * t117 + t231 * t121) * t104 + (-t104 * t145 - t119 * t4) * qJD(6) + (t30 * t170 + ((-t117 * t41 + t121 * t28) * qJD(2) + t3) * qJD(4)) * t123, -t233 * t119 - t56 * t142 + t81 * t6 - t16 * t65 + t24 * t19 + ((-qJD(6) * t28 + t230) * t121 + (qJD(6) * t41 - t231) * t117) * t104 + (-t142 * t170 + (-qJD(2) * t145 - t4) * qJD(4)) * t123; 0, 0, 0, 0, 0, -t127, t138, 0, 0, 0, 0, 0, t197 * t119, t197 * t123, 0, 0, 0, 0, 0, -t123 * t55 - t154 * t209 + (t119 * t83 + (-t107 - t192) * t204) * qJD(4), -t218 + t154 * t211 + (-t107 * t199 + (t85 - t166) * t119) * qJD(4), 0, 0, 0, 0, 0, t104 * t136 + (-t190 * t87 + t130) * t123 + ((-t191 * t87 + t30) * qJD(4) + t104 * t242) * t119, t104 * t135 + (t190 * t86 - t6) * t123 + (t128 * t104 + (t123 * t136 - t142) * qJD(4)) * t119; 0, 0, 0, 0, 0, 0, 0, t123 * t127 * t119, -t198 * t127, 0, 0, 0, qJD(4) * t51 + t123 * t138 + t213, -t147 * t192, t209 * t85 + t221 (t54 - t225) * t122 + (-t55 - t224) * t118, t107 * t185 + (t107 * t202 + (-t85 + t189) * t123) * qJD(2), -t107 * t186 + (-t118 * t210 + (t83 + t181) * t123) * qJD(2), -t107 * t191, -pkin(4) * t55 - t214 - t156 * t107 - t51 * t83 + (-pkin(9) * t209 + t222) * qJD(5) + (t118 * t140 - t17 * t123) * qJD(2), -pkin(4) * t54 + t215 + t229 * t107 - t51 * t85 + (pkin(9) * t211 + t219) * qJD(5) + (t122 * t140 + t18 * t123) * qJD(2), -t142 * t228 + t6 * t87, t130 * t87 + t142 * t227 - t228 * t30 - t6 * t86, t228 * t104 + (qJD(4) * t87 + t142) * t191, -t227 * t104 + (-qJD(4) * t86 + t30) * t191, -t104 * t191, -t111 * t130 + t16 * t86 + t151 * t30 + t227 * t24 + (t117 * t149 - t121 * t150) * t104 + ((-t117 * t99 - t121 * t98) * qJD(4) - t3) * t191, t111 * t6 + t16 * t87 - t151 * t142 + t228 * t24 + (t117 * t150 + t121 * t149) * t104 + (-(-t117 * t98 + t121 * t99) * qJD(4) + t4) * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t83, -t83 ^ 2 + t85 ^ 2, t54 + t225, t224 - t55, t108, t107 * t18 - t37 * t85 + t131, t107 * t17 + t37 * t83 - t139, -t243, t241, t240, t238, t108 -(-t10 * t117 - t223) * t104 + (-t104 * t183 + t108 * t121 - t30 * t85) * pkin(5) + t239, t244 + t9 + (-t104 * t11 - t2) * t117 + (t10 * t104 - t160) * t121 + (-t104 * t182 - t108 * t117 + t142 * t85) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t241, t240, t238, t108, t104 * t4 + t239, t104 * t3 - t233 + t244;];
tauc_reg  = t1;
