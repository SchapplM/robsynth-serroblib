% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:22:01
% EndTime: 2019-12-31 20:22:10
% DurationCPUTime: 3.27s
% Computational Cost: add. (3660->363), mult. (8664->503), div. (0->0), fcn. (6603->14), ass. (0->197)
t160 = sin(pkin(9));
t165 = sin(qJ(2));
t218 = qJD(1) * t165;
t161 = cos(pkin(9));
t169 = cos(qJ(2));
t225 = t161 * t169;
t115 = qJD(1) * t225 - t160 * t218;
t256 = qJD(4) + qJD(5);
t266 = t115 - t256;
t106 = qJD(4) - t115;
t144 = pkin(2) * t160 + pkin(7);
t156 = qJ(2) + pkin(9);
t149 = sin(t156);
t150 = cos(t156);
t166 = sin(qJ(1));
t170 = cos(qJ(1));
t192 = g(1) * t170 + g(2) * t166;
t178 = -g(3) * t150 + t149 * t192;
t162 = -qJ(3) - pkin(6);
t200 = qJD(2) * t162;
t114 = -qJD(3) * t165 + t169 * t200;
t137 = t162 * t165;
t75 = qJDD(2) * pkin(2) + t114 * qJD(1) + qJDD(1) * t137;
t113 = qJD(3) * t169 + t165 * t200;
t138 = t162 * t169;
t83 = t113 * qJD(1) - qJDD(1) * t138;
t30 = -t160 * t83 + t161 * t75;
t28 = -qJDD(2) * pkin(3) - t30;
t265 = -qJD(4) * t144 * t106 + t178 - t28;
t164 = sin(qJ(4));
t167 = cos(qJ(5));
t163 = sin(qJ(5));
t168 = cos(qJ(4));
t224 = t163 * t168;
t130 = t164 * t167 + t224;
t247 = t266 * t130;
t127 = t160 * t169 + t161 * t165;
t117 = t127 * qJD(1);
t213 = t168 * qJD(2);
t94 = t117 * t164 - t213;
t96 = qJD(2) * t164 + t117 * t168;
t187 = t163 * t94 - t167 * t96;
t39 = t163 * t96 + t167 * t94;
t264 = t187 * t39;
t217 = qJD(4) * t164;
t235 = t115 * t164;
t263 = t217 - t235;
t262 = t187 ^ 2 - t39 ^ 2;
t104 = qJD(5) + t106;
t214 = qJD(5) * t167;
t215 = qJD(5) * t163;
t212 = qJD(1) * qJD(2);
t203 = t169 * t212;
t204 = t165 * t212;
t82 = qJDD(1) * t127 - t160 * t204 + t161 * t203;
t35 = qJD(4) * t213 + t164 * qJDD(2) - t117 * t217 + t168 * t82;
t36 = qJD(4) * t96 - t168 * qJDD(2) + t164 * t82;
t7 = -t163 * t36 + t167 * t35 - t94 * t214 - t96 * t215;
t261 = t104 * t39 + t7;
t159 = qJ(4) + qJ(5);
t155 = cos(t159);
t227 = t155 * t170;
t154 = sin(t159);
t230 = t154 * t166;
t100 = t150 * t227 + t230;
t148 = pkin(2) * t169 + pkin(1);
t135 = -t148 * qJD(1) + qJD(3);
t48 = -pkin(3) * t115 - pkin(7) * t117 + t135;
t131 = qJD(1) * t137;
t246 = qJD(2) * pkin(2);
t123 = t131 + t246;
t132 = qJD(1) * t138;
t226 = t161 * t132;
t80 = t160 * t123 - t226;
t69 = qJD(2) * pkin(7) + t80;
t23 = t164 * t48 + t168 * t69;
t17 = -pkin(8) * t94 + t23;
t14 = t17 * t215;
t254 = g(3) * t149;
t120 = t160 * t132;
t79 = t123 * t161 + t120;
t68 = -qJD(2) * pkin(3) - t79;
t37 = pkin(4) * t94 + t68;
t228 = t155 * t166;
t229 = t154 * t170;
t98 = -t150 * t228 + t229;
t260 = g(1) * t100 - g(2) * t98 + t155 * t254 + t37 * t39 + t14;
t179 = pkin(2) * t204 - t148 * qJDD(1) + qJDD(3);
t210 = t169 * qJDD(1);
t211 = t165 * qJDD(1);
t188 = -t160 * t211 + t161 * t210;
t81 = -qJD(2) * t117 + t188;
t24 = -pkin(3) * t81 - pkin(7) * t82 + t179;
t21 = t168 * t24;
t31 = t160 * t75 + t161 * t83;
t29 = qJDD(2) * pkin(7) + t31;
t116 = t127 * qJD(2);
t76 = qJD(1) * t116 + qJDD(4) - t188;
t2 = pkin(4) * t76 - pkin(8) * t35 - qJD(4) * t23 - t164 * t29 + t21;
t216 = qJD(4) * t168;
t184 = t164 * t24 + t168 * t29 + t48 * t216 - t69 * t217;
t3 = -pkin(8) * t36 + t184;
t208 = -t163 * t3 + t167 * t2;
t22 = -t164 * t69 + t168 * t48;
t16 = -pkin(8) * t96 + t22;
t11 = pkin(4) * t106 + t16;
t236 = t167 * t17;
t5 = t11 * t163 + t236;
t97 = t150 * t230 + t227;
t99 = -t150 * t229 + t228;
t259 = -g(1) * t99 + g(2) * t97 - qJD(5) * t5 + t154 * t254 + t37 * t187 + t208;
t174 = qJD(5) * t187 - t163 * t35 - t167 * t36;
t258 = -t104 * t187 + t174;
t64 = t130 * t127;
t129 = t163 * t164 - t167 * t168;
t248 = t266 * t129;
t71 = qJDD(5) + t76;
t255 = -t248 * t104 - t130 * t71;
t252 = g(3) * t169;
t251 = pkin(8) + t144;
t59 = pkin(2) * t218 + pkin(3) * t117 - pkin(7) * t115;
t85 = t131 * t161 + t120;
t250 = t164 * t59 + t168 * t85;
t126 = t160 * t165 - t225;
t78 = pkin(3) * t126 - pkin(7) * t127 - t148;
t93 = t137 * t160 - t138 * t161;
t86 = t168 * t93;
t249 = t164 * t78 + t86;
t245 = t106 * t94;
t244 = t106 * t96;
t243 = t117 * t39;
t242 = t117 * t187;
t241 = t117 * t94;
t240 = t117 * t96;
t238 = t164 * t35;
t237 = t164 * t76;
t119 = t126 * qJD(2);
t234 = t119 * t164;
t233 = t119 * t168;
t232 = t127 * t164;
t231 = t127 * t168;
t223 = t164 * t166;
t222 = t164 * t170;
t221 = t166 * t168;
t220 = t168 * t170;
t157 = t165 ^ 2;
t219 = -t169 ^ 2 + t157;
t209 = t165 * t246;
t145 = -pkin(2) * t161 - pkin(3);
t207 = t127 * t217;
t206 = t127 * t216;
t202 = qJD(5) * t11 + t3;
t199 = qJD(4) * t251;
t198 = -qJD(4) * t48 - t29;
t57 = t113 * t160 - t161 * t114;
t84 = t131 * t160 - t226;
t92 = -t161 * t137 - t138 * t160;
t196 = t106 * t168;
t195 = t247 * t104 - t129 * t71;
t194 = t263 * pkin(4) - t84;
t193 = -t69 * t216 + t21;
t191 = g(1) * t166 - g(2) * t170;
t124 = t251 * t164;
t190 = -pkin(8) * t235 + qJD(5) * t124 + t164 * t199 + t250;
t125 = t251 * t168;
t52 = t168 * t59;
t189 = pkin(4) * t117 + qJD(5) * t125 - t164 * t85 + t52 + (-pkin(8) * t115 + t199) * t168;
t186 = -t263 * t106 + t168 * t76;
t185 = -0.2e1 * pkin(1) * t212 - pkin(6) * qJDD(2);
t58 = t113 * t161 + t114 * t160;
t60 = pkin(3) * t116 + pkin(7) * t119 + t209;
t183 = t164 * t60 + t168 * t58 + t78 * t216 - t93 * t217;
t182 = t206 - t234;
t180 = t106 * t68 - t144 * t76;
t171 = qJD(2) ^ 2;
t176 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t171 + t191;
t172 = qJD(1) ^ 2;
t175 = pkin(1) * t172 - pkin(6) * qJDD(1) + t192;
t136 = -pkin(4) * t168 + t145;
t110 = t150 * t220 + t223;
t109 = -t150 * t222 + t221;
t108 = -t150 * t221 + t222;
t107 = t150 * t223 + t220;
t67 = t168 * t78;
t65 = t129 * t127;
t56 = pkin(4) * t232 + t92;
t53 = t168 * t60;
t27 = pkin(4) * t182 + t57;
t25 = -pkin(8) * t232 + t249;
t18 = pkin(4) * t126 - pkin(8) * t231 - t164 * t93 + t67;
t13 = -t119 * t224 - t163 * t207 - t215 * t232 + (t256 * t231 - t234) * t167;
t12 = t129 * t119 - t256 * t64;
t10 = pkin(4) * t36 + t28;
t9 = -pkin(8) * t182 + t183;
t6 = pkin(8) * t233 + pkin(4) * t116 - t164 * t58 + t53 + (-t86 + (pkin(8) * t127 - t78) * t164) * qJD(4);
t4 = t11 * t167 - t163 * t17;
t1 = [qJDD(1), t191, t192, qJDD(1) * t157 + 0.2e1 * t165 * t203, 0.2e1 * t165 * t210 - 0.2e1 * t219 * t212, qJDD(2) * t165 + t169 * t171, qJDD(2) * t169 - t165 * t171, 0, t165 * t185 + t169 * t176, -t165 * t176 + t169 * t185, t115 * t58 - t116 * t80 + t117 * t57 + t119 * t79 - t126 * t31 - t127 * t30 + t81 * t93 + t82 * t92 - t192, t31 * t93 + t80 * t58 - t30 * t92 - t79 * t57 - t179 * t148 + t135 * t209 - g(1) * (-t148 * t166 - t162 * t170) - g(2) * (t148 * t170 - t162 * t166), -t96 * t233 + (t168 * t35 - t96 * t217) * t127, -(-t164 * t96 - t168 * t94) * t119 + (-t238 - t168 * t36 + (t164 * t94 - t168 * t96) * qJD(4)) * t127, t76 * t231 + t116 * t96 + t126 * t35 + (-t207 - t233) * t106, -t106 * t182 - t116 * t94 - t126 * t36 - t76 * t232, t106 * t116 + t126 * t76, (-t93 * t216 + t53) * t106 + t67 * t76 + t193 * t126 + t22 * t116 + t57 * t94 + t92 * t36 + t68 * t206 - g(1) * t108 - g(2) * t110 + ((-qJD(4) * t78 - t58) * t106 - t93 * t76 + t198 * t126 + t28 * t127 - t68 * t119) * t164, -t183 * t106 - t249 * t76 - t184 * t126 - t23 * t116 + t57 * t96 + t92 * t35 - t68 * t233 - g(1) * t107 - g(2) * t109 + (t28 * t168 - t68 * t217) * t127, -t12 * t187 - t65 * t7, -t12 * t39 + t13 * t187 - t174 * t65 - t64 * t7, t104 * t12 - t116 * t187 + t126 * t7 - t65 * t71, -t104 * t13 - t116 * t39 + t126 * t174 - t64 * t71, t104 * t116 + t126 * t71, (-t163 * t9 + t167 * t6) * t104 + (-t163 * t25 + t167 * t18) * t71 + t208 * t126 + t4 * t116 + t27 * t39 - t56 * t174 + t10 * t64 + t37 * t13 - g(1) * t98 - g(2) * t100 + ((-t163 * t18 - t167 * t25) * t104 - t5 * t126) * qJD(5), -g(1) * t97 - g(2) * t99 - t10 * t65 - t5 * t116 + t37 * t12 + t14 * t126 - t27 * t187 + t56 * t7 + (-(-qJD(5) * t25 + t6) * t104 - t18 * t71 - t2 * t126) * t163 + (-(qJD(5) * t18 + t9) * t104 - t25 * t71 - t202 * t126) * t167; 0, 0, 0, -t165 * t172 * t169, t219 * t172, t211, t210, qJDD(2), t165 * t175 - t252, g(3) * t165 + t169 * t175, (t80 - t84) * t117 + (t79 - t85) * t115 + (t160 * t81 - t161 * t82) * pkin(2), t79 * t84 - t80 * t85 + (-t252 + t160 * t31 + t161 * t30 + (-qJD(1) * t135 + t192) * t165) * pkin(2), t196 * t96 + t238, (t35 - t245) * t168 + (-t36 - t244) * t164, t106 * t196 + t237 - t240, t186 + t241, -t106 * t117, -t52 * t106 - t22 * t117 + t145 * t36 - t84 * t94 + (t85 * t106 + t180) * t164 + t265 * t168, t250 * t106 + t23 * t117 + t145 * t35 - t265 * t164 + t180 * t168 - t84 * t96, t130 * t7 - t187 * t248, -t129 * t7 + t130 * t174 - t187 * t247 - t248 * t39, t242 - t255, t195 + t243, -t104 * t117, (-t124 * t167 - t125 * t163) * t71 - t136 * t174 + t10 * t129 - t4 * t117 + t194 * t39 - t247 * t37 + (t163 * t190 - t167 * t189) * t104 + t178 * t155, -(-t124 * t163 + t125 * t167) * t71 + t136 * t7 + t10 * t130 + t5 * t117 - t194 * t187 + t248 * t37 + (t163 * t189 + t167 * t190) * t104 - t178 * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 ^ 2 - t117 ^ 2, -t115 * t80 + t117 * t79 + t179 - t191, 0, 0, 0, 0, 0, t186 - t241, -t106 ^ 2 * t168 - t237 - t240, 0, 0, 0, 0, 0, t195 - t243, t242 + t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t94, -t94 ^ 2 + t96 ^ 2, t35 + t245, t244 - t36, t76, -g(1) * t109 + g(2) * t107 + t106 * t23 - t68 * t96 + (t198 + t254) * t164 + t193, g(1) * t110 - g(2) * t108 + t106 * t22 + t168 * t254 + t68 * t94 - t184, -t264, t262, t261, t258, t71, -(-t16 * t163 - t236) * t104 + (-t104 * t215 + t167 * t71 - t39 * t96) * pkin(4) + t259, (-t104 * t17 - t2) * t163 + (t104 * t16 - t202) * t167 + (-t104 * t214 - t163 * t71 + t187 * t96) * pkin(4) + t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t264, t262, t261, t258, t71, t104 * t5 + t259, t104 * t4 - t163 * t2 - t167 * t202 + t260;];
tau_reg = t1;
