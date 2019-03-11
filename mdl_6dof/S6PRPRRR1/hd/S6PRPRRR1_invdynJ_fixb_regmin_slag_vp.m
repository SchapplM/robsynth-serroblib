% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:09
% EndTime: 2019-03-08 20:25:17
% DurationCPUTime: 3.54s
% Computational Cost: add. (3353->340), mult. (7607->490), div. (0->0), fcn. (6444->16), ass. (0->205)
t152 = qJD(4) + qJD(5);
t163 = sin(qJ(5));
t167 = cos(qJ(4));
t270 = cos(qJ(5));
t220 = qJD(2) * t270;
t164 = sin(qJ(4));
t237 = qJD(2) * t164;
t283 = -t163 * t237 + t167 * t220;
t284 = t283 * t152;
t101 = qJD(6) - t283;
t282 = t101 - qJD(6);
t155 = qJ(4) + qJ(5);
t149 = sin(t155);
t150 = cos(t155);
t161 = cos(pkin(6));
t156 = sin(pkin(12));
t165 = sin(qJ(2));
t159 = cos(pkin(12));
t168 = cos(qJ(2));
t242 = t168 * t159;
t108 = t156 * t165 - t242;
t157 = sin(pkin(11));
t160 = cos(pkin(11));
t192 = t156 * t168 + t159 * t165;
t98 = t192 * t161;
t195 = -t108 * t157 + t160 * t98;
t196 = -t108 * t160 - t157 * t98;
t158 = sin(pkin(6));
t247 = t158 * t160;
t248 = t157 * t158;
t97 = t192 * t158;
t185 = -g(3) * (-t149 * t97 + t150 * t161) - g(2) * (-t149 * t195 - t150 * t247) - g(1) * (-t149 * t196 + t150 * t248);
t151 = qJDD(4) + qJDD(5);
t130 = t161 * qJDD(1) + qJDD(3);
t118 = t167 * t130;
t238 = qJD(2) * t158;
t218 = qJD(1) * t238;
t231 = qJDD(1) * t158;
t280 = t165 * t231 + t168 * t218;
t129 = t168 * t231;
t95 = qJDD(2) * pkin(2) - t165 * t218 + t129;
t61 = t156 * t95 + t280 * t159;
t54 = qJDD(2) * pkin(8) + t61;
t210 = pkin(9) * qJDD(2) + t54;
t137 = qJD(1) * t161 + qJD(3);
t239 = qJD(1) * t158;
t223 = t168 * t239;
t114 = qJD(2) * pkin(2) + t223;
t224 = t165 * t239;
t120 = t159 * t224;
t84 = t156 * t114 + t120;
t215 = t84 + (pkin(8) + pkin(9)) * qJD(2);
t58 = t137 * t164 + t215 * t167;
t15 = qJDD(4) * pkin(4) - qJD(4) * t58 - t210 * t164 + t118;
t57 = t167 * t137 - t215 * t164;
t16 = qJD(4) * t57 + t164 * t130 + t210 * t167;
t226 = t270 * t58;
t50 = qJD(4) * pkin(4) + t57;
t22 = t163 * t50 + t226;
t277 = t22 * qJD(5) - t270 * t15 + t163 * t16;
t3 = -t151 * pkin(5) + t277;
t181 = t185 - t3;
t189 = (g(1) * t157 - g(2) * t160) * t158;
t281 = -t130 + t189;
t182 = t108 * t161;
t63 = -t157 * t192 - t160 * t182;
t66 = t157 * t182 - t160 * t192;
t246 = t158 * t165;
t96 = t156 * t246 - t158 * t242;
t184 = g(1) * t66 + g(2) * t63 - g(3) * t96;
t179 = t184 * t150;
t211 = qJDD(2) * t270;
t230 = t164 * qJDD(2);
t200 = t163 * t230 - t167 * t211;
t243 = t163 * t167;
t111 = t270 * t164 + t243;
t78 = t152 * t111;
t56 = qJD(2) * t78 + t200;
t51 = qJDD(6) + t56;
t225 = -pkin(4) * t167 - pkin(3);
t269 = pkin(2) * t159;
t122 = t225 - t269;
t186 = -t163 * t164 + t270 * t167;
t69 = -pkin(5) * t186 - pkin(10) * t111 + t122;
t279 = t69 * t51 - t179;
t162 = sin(qJ(6));
t79 = t161 * t167 - t164 * t97;
t80 = t161 * t164 + t167 * t97;
t35 = t163 * t79 + t270 * t80;
t166 = cos(qJ(6));
t90 = t96 * t166;
t278 = -t162 * t35 + t90;
t236 = qJD(4) * t164;
t91 = t156 * t223 + t120;
t205 = pkin(4) * t236 - t91;
t234 = qJD(6) * t162;
t255 = t166 * t51;
t77 = t152 * t186;
t276 = -t101 * (t111 * t234 - t166 * t77) + t111 * t255;
t256 = t163 * t58;
t21 = t270 * t50 - t256;
t19 = -t152 * pkin(5) - t21;
t219 = qJD(5) * t270;
t235 = qJD(5) * t163;
t177 = t163 * t15 + t270 * t16 + t50 * t219 - t58 * t235;
t2 = t151 * pkin(10) + t177;
t143 = pkin(2) * t156 + pkin(8);
t267 = pkin(9) + t143;
t212 = qJD(4) * t267;
t100 = t167 * t212;
t106 = t267 * t164;
t107 = t267 * t167;
t187 = -t270 * t106 - t163 * t107;
t119 = t156 * t224;
t94 = t159 * t223 - t119;
t99 = t164 * t212;
t265 = -t187 * qJD(5) + t163 * t100 + t186 * t94 + t270 * t99;
t274 = -g(1) * t196 - g(2) * t195;
t104 = -qJD(2) * t243 - t164 * t220;
t83 = t114 * t159 - t119;
t73 = t225 * qJD(2) - t83;
t33 = -pkin(5) * t283 + pkin(10) * t104 + t73;
t72 = -t163 * t106 + t270 * t107;
t275 = (qJD(6) * t33 + t2) * t186 + t19 * t77 + t3 * t111 + (-qJD(6) * t69 + t265) * t101 - g(3) * t97 - t72 * t51 + t274;
t193 = t104 * t166 - t152 * t162;
t229 = t167 * qJDD(2);
t55 = t163 * t229 + t164 * t211 + t284;
t27 = -t193 * qJD(6) - t166 * t151 + t162 * t55;
t233 = qJD(6) * t166;
t26 = t104 * t234 + t162 * t151 + t152 * t233 + t166 * t55;
t266 = -t186 * t26 - t193 * t78;
t264 = t72 * qJD(5) + t270 * t100 - t111 * t94 - t163 * t99;
t263 = pkin(5) * t78 - pkin(10) * t77 + t205;
t85 = -t104 * t162 - t166 * t152;
t262 = t101 * t85;
t261 = t101 * t193;
t260 = t162 * t26;
t258 = t162 * t51;
t257 = t162 * t96;
t254 = t166 * t193;
t253 = t19 * t283;
t252 = t19 * t111;
t251 = t101 * t104;
t250 = t104 * t283;
t245 = t161 * t165;
t244 = t161 * t168;
t241 = qJDD(1) - g(3);
t153 = t164 ^ 2;
t240 = -t167 ^ 2 + t153;
t232 = qJD(2) * qJD(4);
t20 = t152 * pkin(10) + t22;
t198 = t162 * t20 - t166 * t33;
t221 = -t104 * t198 + t19 * t234;
t217 = t164 * t232;
t146 = pkin(4) * t163 + pkin(10);
t74 = -pkin(5) * t104 - pkin(10) * t283;
t208 = pkin(4) * t237 + qJD(6) * t146 + t74;
t207 = t101 * t166;
t23 = t163 * t57 + t226;
t204 = pkin(4) * t235 - t23;
t60 = -t280 * t156 + t159 * t95;
t202 = t186 * t27 - t78 * t85;
t199 = -t146 * t51 - t253;
t8 = t162 * t33 + t166 * t20;
t197 = t166 * t35 + t257;
t194 = t111 * t151 + t152 * t77;
t24 = t270 * t57 - t256;
t191 = -pkin(4) * t219 + t24;
t190 = -t163 * t80 + t270 * t79;
t183 = -t8 * t104 - t181 * t162 + t19 * t233;
t81 = -qJD(2) * pkin(3) - t83;
t180 = -qJD(2) * t81 - t274 - t54;
t144 = -pkin(3) - t269;
t178 = -qJDD(4) * t143 + (qJD(2) * t144 + t81 + t94) * qJD(4);
t36 = pkin(4) * t217 + t225 * qJDD(2) - t60;
t176 = (-t111 * t233 - t162 * t77) * t101 - t111 * t258;
t169 = qJD(4) ^ 2;
t175 = -qJD(2) * t91 + t143 * t169 + t184 - t60 + (-pkin(3) + t144) * qJDD(2);
t174 = -g(1) * (-t157 * t244 - t160 * t165) - g(2) * (-t157 * t165 + t160 * t244) - g(3) * t158 * t168;
t172 = t73 * t104 + t185 - t277;
t42 = -t149 * t247 + t150 * t195;
t44 = t149 * t248 + t150 * t196;
t76 = t149 * t161 + t150 * t97;
t171 = g(1) * t44 + g(2) * t42 + g(3) * t76 - t283 * t73 - t177;
t170 = qJD(2) ^ 2;
t147 = -t270 * pkin(4) - pkin(5);
t124 = qJDD(4) * t167 - t164 * t169;
t123 = qJDD(4) * t164 + t167 * t169;
t93 = t108 * t238;
t92 = qJD(2) * t97;
t59 = t104 ^ 2 - t283 ^ 2;
t52 = t151 * t186 - t152 * t78;
t40 = t79 * qJD(4) - t167 * t93;
t39 = -t80 * qJD(4) + t164 * t93;
t31 = -t200 + (-qJD(2) * t111 - t104) * t152;
t30 = t55 - t284;
t13 = t101 * t207 - t104 * t193 + t258;
t12 = -t101 ^ 2 * t162 - t104 * t85 + t255;
t11 = pkin(5) * t56 - pkin(10) * t55 + t36;
t10 = -t193 * t207 + t260;
t9 = t166 * t11;
t6 = t35 * qJD(5) + t163 * t40 - t270 * t39;
t5 = t190 * qJD(5) + t163 * t39 + t270 * t40;
t4 = (t26 - t262) * t166 + (-t27 + t261) * t162;
t1 = [t241, 0 (qJDD(2) * t168 - t165 * t170) * t158 (-qJDD(2) * t165 - t168 * t170) * t158, t130 * t161 - t60 * t96 + t61 * t97 - t83 * t92 - t84 * t93 - g(3), 0, 0, 0, 0, 0, -t96 * t229 + qJD(4) * t39 + qJDD(4) * t79 + (-t167 * t92 + t96 * t236) * qJD(2), t96 * t230 - qJD(4) * t40 - qJDD(4) * t80 + (qJD(4) * t167 * t96 + t164 * t92) * qJD(2), 0, 0, 0, 0, 0, t151 * t190 - t152 * t6 - t283 * t92 + t56 * t96, -t104 * t92 - t151 * t35 - t152 * t5 + t55 * t96, 0, 0, 0, 0, 0 (-qJD(6) * t197 - t162 * t5 + t166 * t92) * t101 + t278 * t51 + t6 * t85 - t190 * t27 -(qJD(6) * t278 + t162 * t92 + t166 * t5) * t101 - t197 * t51 - t6 * t193 - t190 * t26; 0, qJDD(2), t129 + t174, -g(1) * (t157 * t245 - t160 * t168) - g(2) * (-t157 * t168 - t160 * t245) - t241 * t246, t83 * t91 - t84 * t94 + (t61 * t156 + t60 * t159 + t174) * pkin(2), qJDD(2) * t153 + 0.2e1 * t167 * t217, 0.2e1 * t164 * t229 - 0.2e1 * t240 * t232, t123, t124, 0, t164 * t178 - t167 * t175, t164 * t175 + t167 * t178, -t104 * t77 + t111 * t55, t104 * t78 - t111 * t56 + t186 * t55 + t283 * t77, t194, t52, 0, t122 * t56 + t151 * t187 - t264 * t152 - t186 * t36 - t205 * t283 + t73 * t78 - t179, -t205 * t104 + t111 * t36 + t122 * t55 + t184 * t149 - t151 * t72 + t265 * t152 + t73 * t77, -t77 * t254 + (t166 * t26 + t193 * t234) * t111 (t162 * t193 - t166 * t85) * t77 + (-t260 - t166 * t27 + (t162 * t85 + t254) * qJD(6)) * t111, t266 + t276, t176 + t202, t101 * t78 - t186 * t51, -t9 * t186 - t187 * t27 - t198 * t78 + t264 * t85 + (t263 * t101 + (-t72 * t101 + t186 * t20 + t252) * qJD(6) + t279) * t166 + t275 * t162, -t187 * t26 - t8 * t78 - t264 * t193 + ((-qJD(6) * t20 + t11) * t186 - qJD(6) * t252 + (qJD(6) * t72 - t263) * t101 - t279) * t162 + t275 * t166; 0, 0, 0, 0, -g(3) * t161 - t281, 0, 0, 0, 0, 0, t124, -t123, 0, 0, 0, 0, 0, t52, -t194, 0, 0, 0, 0, 0, t176 - t202, t266 - t276; 0, 0, 0, 0, 0, -t164 * t170 * t167, t240 * t170, t230, t229, qJDD(4), -g(3) * t79 + t180 * t164 - t167 * t189 + t118, g(3) * t80 + t281 * t164 + t180 * t167, t250, t59, t30, t31, t151, t23 * t152 + (t270 * t151 - t152 * t235 + t237 * t283) * pkin(4) + t172, t24 * t152 + (t104 * t237 - t151 * t163 - t152 * t219) * pkin(4) + t171, t10, t4, t13, t12, t251, t147 * t27 + t204 * t85 + (t101 * t191 + t199) * t162 + (-t101 * t208 + t181) * t166 + t221, t147 * t26 - t204 * t193 + t199 * t166 + (t162 * t208 + t166 * t191) * t101 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, t59, t30, t31, t151, t22 * t152 + t172, t21 * t152 + t171, t10, t4, t13, t12, t251, -pkin(5) * t27 - t22 * t85 + (-pkin(10) * t51 + t21 * t101 - t253) * t162 + ((-pkin(10) * qJD(6) - t74) * t101 + t181) * t166 + t221, -pkin(5) * t26 + (t162 * t74 + t166 * t21) * t101 + t22 * t193 - t166 * t253 + (t101 * t234 - t255) * pkin(10) + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193 * t85, t193 ^ 2 - t85 ^ 2, t26 + t262, -t27 - t261, t51, -t162 * t2 + t9 + t19 * t193 - g(1) * (-t162 * t44 - t166 * t66) - g(2) * (-t162 * t42 - t166 * t63) - g(3) * (-t162 * t76 + t90) + t282 * t8, -t166 * t2 - t162 * t11 + t19 * t85 - g(1) * (t162 * t66 - t166 * t44) - g(2) * (t162 * t63 - t166 * t42) - g(3) * (-t166 * t76 - t257) - t282 * t198;];
tau_reg  = t1;
