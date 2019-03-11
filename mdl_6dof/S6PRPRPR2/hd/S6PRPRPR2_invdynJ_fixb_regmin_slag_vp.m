% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:57
% EndTime: 2019-03-08 19:33:06
% DurationCPUTime: 3.63s
% Computational Cost: add. (3193->417), mult. (7613->608), div. (0->0), fcn. (6475->16), ass. (0->212)
t179 = sin(qJ(4));
t182 = cos(qJ(4));
t216 = pkin(4) * t179 - qJ(5) * t182;
t110 = qJD(4) * t216 - t179 * qJD(5);
t175 = cos(pkin(11));
t180 = sin(qJ(2));
t173 = sin(pkin(6));
t251 = qJD(1) * t173;
t234 = t180 * t251;
t137 = t175 * t234;
t171 = sin(pkin(11));
t183 = cos(qJ(2));
t233 = t183 * t251;
t92 = t171 * t233 + t137;
t299 = t110 - t92;
t248 = qJD(2) * t182;
t154 = -qJD(6) + t248;
t170 = sin(pkin(12));
t174 = cos(pkin(12));
t242 = t174 * qJD(4);
t249 = qJD(2) * t179;
t119 = t170 * t249 - t242;
t247 = qJD(4) * t170;
t121 = t174 * t249 + t247;
t178 = sin(qJ(6));
t181 = cos(qJ(6));
t208 = t119 * t178 - t121 * t181;
t298 = t154 * t208;
t250 = qJD(2) * t173;
t229 = qJD(1) * t250;
t240 = qJDD(1) * t173;
t297 = t180 * t240 + t183 * t229;
t158 = pkin(2) * t171 + pkin(8);
t246 = qJD(4) * t179;
t231 = t158 * t246;
t268 = t170 * t182;
t136 = t171 * t234;
t95 = t175 * t233 - t136;
t281 = t170 * t231 + t299 * t174 + t268 * t95;
t260 = t174 * t182;
t296 = t299 * t170 - t260 * t95;
t177 = cos(pkin(6));
t153 = qJD(1) * t177 + qJD(3);
t133 = qJD(2) * pkin(2) + t233;
t83 = t171 * t133 + t137;
t81 = qJD(2) * pkin(8) + t83;
t295 = t153 * t182 - t179 * t81;
t166 = t182 * qJDD(2);
t241 = qJD(2) * qJD(4);
t294 = t179 * t241 - t166;
t172 = sin(pkin(10));
t176 = cos(pkin(10));
t256 = t177 * t183;
t293 = -t172 * t256 - t176 * t180;
t292 = -t154 - qJD(6);
t291 = -qJDD(4) * pkin(4) + qJDD(5);
t255 = t183 * t175;
t264 = t173 * t180;
t102 = t171 * t264 - t173 * t255;
t125 = t171 * t180 - t255;
t193 = t125 * t177;
t207 = t171 * t183 + t175 * t180;
t60 = -t172 * t207 - t176 * t193;
t63 = t172 * t193 - t176 * t207;
t290 = -g(1) * t63 - g(2) * t60 + g(3) * t102;
t257 = t177 * t180;
t253 = -t171 * t256 - t175 * t257;
t64 = -t176 * t125 + t172 * t253;
t59 = t172 * t125 + t176 * t253;
t161 = t174 * qJDD(4);
t227 = t182 * t241;
t238 = t179 * qJDD(2);
t192 = t227 + t238;
t88 = t170 * t192 - t161;
t239 = qJDD(4) * t170;
t89 = t174 * t192 + t239;
t17 = -qJD(6) * t208 + t178 * t89 + t181 * t88;
t151 = t177 * qJDD(1) + qJDD(3);
t273 = t151 * t179;
t150 = t183 * t240;
t98 = qJDD(2) * pkin(2) - t180 * t229 + t150;
t47 = t171 * t98 + t297 * t175;
t41 = qJDD(2) * pkin(8) + t47;
t10 = qJDD(4) * qJ(5) + t273 + t182 * t41 + (qJD(5) + t295) * qJD(4);
t204 = pkin(4) * t182 + qJ(5) * t179 + pkin(3);
t46 = -t297 * t171 + t175 * t98;
t22 = t110 * qJD(2) - qJDD(2) * t204 - t46;
t7 = t174 * t10 + t170 * t22;
t287 = pkin(2) * t175;
t285 = pkin(9) + qJ(5);
t205 = pkin(5) * t179 - pkin(9) * t260;
t284 = -qJD(4) * t205 - t281;
t261 = t174 * t179;
t283 = (-pkin(9) * t268 - t158 * t261) * qJD(4) + t296;
t53 = t179 * t153 + t182 * t81;
t49 = qJD(4) * qJ(5) + t53;
t82 = t133 * t175 - t136;
t58 = -qJD(2) * t204 - t82;
t15 = t170 * t58 + t174 * t49;
t126 = t170 * t181 + t174 * t178;
t104 = t126 * t179;
t123 = qJDD(6) + t294;
t195 = t126 * t182;
t243 = qJD(6) * t181;
t244 = qJD(6) * t179;
t269 = t170 * t178;
t57 = qJD(4) * t195 + t243 * t261 - t244 * t269;
t282 = -t104 * t123 + t57 * t154;
t280 = -t174 * t231 + t296;
t274 = t121 * t178;
t69 = t181 * t119 + t274;
t279 = t154 * t69;
t45 = -qJD(4) * pkin(4) + qJD(5) - t295;
t278 = t179 * t45;
t124 = -t181 * t174 + t269;
t194 = t124 * t182;
t277 = qJD(2) * t194 - t124 * qJD(6);
t276 = -qJD(2) * t195 + t126 * qJD(6);
t129 = t216 * qJD(2);
t26 = t170 * t129 + t174 * t295;
t167 = pkin(12) + qJ(6);
t164 = sin(t167);
t271 = t164 * t182;
t165 = cos(t167);
t270 = t165 * t182;
t266 = t172 * t180;
t265 = t173 * t179;
t263 = t173 * t182;
t262 = t173 * t183;
t254 = qJDD(1) - g(3);
t117 = -t204 - t287;
t77 = t170 * t117 + t158 * t260;
t168 = t179 ^ 2;
t252 = -t182 ^ 2 + t168;
t245 = qJD(4) * t182;
t236 = t176 * t256;
t6 = -t10 * t170 + t174 * t22;
t2 = t294 * pkin(5) - pkin(9) * t89 + t6;
t5 = -pkin(9) * t88 + t7;
t235 = -t178 * t5 + t181 * t2;
t232 = t170 * t248;
t230 = qJ(5) * t166;
t226 = t170 * t238;
t225 = t174 * t238;
t223 = pkin(5) * t170 + t158;
t16 = -qJD(6) * t274 - t119 * t243 - t178 * t88 + t181 * t89;
t222 = -t16 * t182 - t208 * t246;
t14 = -t170 * t49 + t174 * t58;
t25 = t174 * t129 - t170 * t295;
t220 = qJD(2) * t252;
t218 = -t170 * t6 + t174 * t7;
t217 = t178 * t2 + t181 * t5;
t12 = -pkin(5) * t248 - pkin(9) * t121 + t14;
t13 = -pkin(9) * t119 + t15;
t3 = t12 * t181 - t13 * t178;
t4 = t12 * t178 + t13 * t181;
t215 = -t14 * t170 + t15 * t174;
t103 = t207 * t173;
t79 = t103 * t182 + t177 * t179;
t30 = t102 * t174 - t170 * t79;
t31 = t102 * t170 + t174 * t79;
t214 = -t178 * t31 + t181 * t30;
t213 = t178 * t30 + t181 * t31;
t107 = t174 * t117;
t55 = -pkin(9) * t261 + t107 + (-t158 * t170 - pkin(5)) * t182;
t66 = -pkin(9) * t170 * t179 + t77;
t212 = -t178 * t66 + t181 * t55;
t211 = t178 * t55 + t181 * t66;
t105 = t124 * t179;
t56 = -qJD(4) * t194 - t126 * t244;
t210 = t105 * t123 + t154 * t56;
t78 = t103 * t179 - t177 * t182;
t206 = -t151 * t182 + t153 * t246 + t179 * t41 + t81 * t245;
t203 = t17 * t182 - t246 * t69;
t142 = t285 * t170;
t202 = pkin(9) * t232 + qJD(5) * t174 - qJD(6) * t142 - t26;
t143 = t285 * t174;
t201 = qJD(2) * t205 + qJD(5) * t170 + qJD(6) * t143 + t25;
t200 = g(1) * (-t172 * t263 + t179 * t64) + g(2) * (t176 * t263 - t179 * t59) + g(3) * t78;
t36 = -t176 * t265 - t182 * t59;
t38 = t172 * t265 + t182 * t64;
t199 = g(1) * t38 + g(2) * t36 + g(3) * t79;
t197 = -g(1) * t64 + g(2) * t59 - g(3) * t103;
t196 = -g(3) * t177 + (-g(1) * t172 + g(2) * t176) * t173;
t11 = t206 + t291;
t191 = -t11 + t200;
t190 = -qJ(5) * t246 + (qJD(5) - t45) * t182;
t159 = -pkin(3) - t287;
t80 = -qJD(2) * pkin(3) - t82;
t189 = -qJDD(4) * t158 + (qJD(2) * t159 + t80 + t95) * qJD(4);
t188 = -g(1) * t293 - g(3) * t262;
t187 = t200 - t206;
t184 = qJD(4) ^ 2;
t186 = -qJD(2) * t92 + t158 * t184 - t290 - t46 + (-pkin(3) + t159) * qJDD(2);
t185 = qJD(2) ^ 2;
t160 = -pkin(5) * t174 - pkin(4);
t144 = pkin(2) * t236;
t140 = qJDD(4) * t182 - t179 * t184;
t139 = qJDD(4) * t179 + t182 * t184;
t109 = t223 * t179;
t101 = t223 * t245;
t94 = t125 * t250;
t93 = qJD(2) * t103;
t76 = -t158 * t268 + t107;
t39 = pkin(5) * t232 + t53;
t29 = -qJD(4) * t78 - t94 * t182;
t28 = qJD(4) * t79 - t94 * t179;
t27 = pkin(5) * t119 + t45;
t19 = t170 * t93 + t174 * t29;
t18 = -t170 * t29 + t174 * t93;
t8 = pkin(5) * t88 + t11;
t1 = [t254, 0 (qJDD(2) * t183 - t180 * t185) * t173 (-qJDD(2) * t180 - t183 * t185) * t173, -t102 * t46 + t103 * t47 + t151 * t177 - t82 * t93 - t83 * t94 - g(3), 0, 0, 0, 0, 0, -t102 * t166 - qJD(4) * t28 - qJDD(4) * t78 + (t102 * t246 - t182 * t93) * qJD(2), t102 * t238 - qJD(4) * t29 - qJDD(4) * t79 + (t102 * t245 + t179 * t93) * qJD(2), -t30 * t166 + t119 * t28 + t78 * t88 + (-t18 * t182 + t246 * t30) * qJD(2), t31 * t166 + t121 * t28 + t78 * t89 + (t182 * t19 - t246 * t31) * qJD(2), -t119 * t19 - t121 * t18 - t30 * t89 - t31 * t88, t11 * t78 + t14 * t18 + t15 * t19 + t28 * t45 + t30 * t6 + t31 * t7 - g(3), 0, 0, 0, 0, 0 -(-qJD(6) * t213 - t178 * t19 + t18 * t181) * t154 + t214 * t123 + t28 * t69 + t78 * t17 (qJD(6) * t214 + t178 * t18 + t181 * t19) * t154 - t213 * t123 - t28 * t208 + t78 * t16; 0, qJDD(2), t150 - g(2) * (t236 - t266) + t188, -g(1) * (t172 * t257 - t176 * t183) - g(2) * (-t172 * t183 - t176 * t257) - t254 * t264, -g(2) * t144 + t82 * t92 - t83 * t95 + (g(2) * t266 + t47 * t171 + t46 * t175 + t188) * pkin(2), qJDD(2) * t168 + 0.2e1 * t179 * t227, -0.2e1 * qJD(4) * t220 + 0.2e1 * t166 * t179, t139, t140, 0, t179 * t189 - t182 * t186, t179 * t186 + t182 * t189, t197 * t170 + (t11 * t170 - t95 * t119 + t158 * t88 + (qJD(2) * t76 + t14) * qJD(4)) * t179 + (-t76 * qJDD(2) - t6 + (t119 * t158 + t170 * t45) * qJD(4) - t281 * qJD(2) + t290 * t174) * t182, t197 * t174 + (t11 * t174 - t95 * t121 + t158 * t89 + (-qJD(2) * t77 - t15) * qJD(4)) * t179 + (t77 * qJDD(2) + t7 + (t121 * t158 + t174 * t45) * qJD(4) + t280 * qJD(2) - t290 * t170) * t182, -t76 * t89 - t77 * t88 - t281 * t121 - t280 * t119 + (-t14 * t174 - t15 * t170) * t245 + (-t170 * t7 - t174 * t6 + t290) * t179, t7 * t77 + t6 * t76 - t95 * t278 - g(1) * (t293 * pkin(2) + pkin(8) * t64) - g(2) * (-pkin(2) * t266 - pkin(8) * t59 + t144) - g(3) * (pkin(2) * t262 + pkin(8) * t103) + (t11 * t179 + t245 * t45) * t158 + t280 * t15 + t281 * t14 + t290 * t204, -t105 * t16 - t208 * t56, -t104 * t16 + t105 * t17 + t208 * t57 - t56 * t69, -t210 + t222, t203 + t282, -t123 * t182 - t154 * t246, t212 * t123 - t235 * t182 + t101 * t69 + t109 * t17 + t8 * t104 + t27 * t57 - g(1) * (t164 * t64 + t270 * t63) - g(2) * (-t164 * t59 + t270 * t60) - g(3) * (-t102 * t270 + t103 * t164) + (qJD(4) * t3 - t69 * t95) * t179 + (t283 * t178 + t284 * t181) * t154 + (t154 * t211 + t182 * t4) * qJD(6), -t211 * t123 + t217 * t182 - t101 * t208 + t109 * t16 - t8 * t105 + t27 * t56 - g(1) * (t165 * t64 - t271 * t63) - g(2) * (-t165 * t59 - t271 * t60) - g(3) * (t102 * t271 + t103 * t165) + (-qJD(4) * t4 + t208 * t95) * t179 + (-t284 * t178 + t283 * t181) * t154 + (t154 * t212 + t182 * t3) * qJD(6); 0, 0, 0, 0, t196 + t151, 0, 0, 0, 0, 0, t140, -t139 (-t88 + t226) * t182 + (t119 * t179 - t170 * t220) * qJD(4) (-t89 + t225) * t182 + (t121 * t179 - t174 * t220) * qJD(4) (t170 * t89 - t174 * t88) * t179 + (-t119 * t174 + t121 * t170) * t245, -t11 * t182 + t218 * t179 + (t182 * t215 + t278) * qJD(4) + t196, 0, 0, 0, 0, 0, -t203 + t282, t210 + t222; 0, 0, 0, 0, 0, -t179 * t185 * t182, t252 * t185, t238, t166, qJDD(4), qJD(4) * t53 - t249 * t80 + t187, -t273 + (-qJD(2) * t80 - t41) * t182 + t199, t170 * t230 - pkin(4) * t88 - t119 * t53 + t191 * t174 + (-t14 * t179 + t170 * t190 + t182 * t25) * qJD(2), t174 * t230 - pkin(4) * t89 - t121 * t53 - t191 * t170 + (t15 * t179 + t174 * t190 - t182 * t26) * qJD(2), t119 * t26 + t121 * t25 + (-qJ(5) * t88 - qJD(5) * t119 + t14 * t248 + t7) * t174 + (qJ(5) * t89 + qJD(5) * t121 + t15 * t248 - t6) * t170 - t199, -t14 * t25 - t15 * t26 - t45 * t53 + t215 * qJD(5) + t191 * pkin(4) + (-t199 + t218) * qJ(5), t16 * t126 - t208 * t277, -t124 * t16 - t126 * t17 + t208 * t276 - t277 * t69, t123 * t126 - t277 * t154 + t208 * t249, -t123 * t124 + t276 * t154 + t69 * t249, t154 * t249 (-t142 * t181 - t143 * t178) * t123 + t160 * t17 + t8 * t124 - t3 * t249 - t39 * t69 + t276 * t27 + (t178 * t202 + t181 * t201) * t154 + t200 * t165 -(-t142 * t178 + t143 * t181) * t123 + t160 * t16 + t8 * t126 + t4 * t249 + t39 * t208 + t277 * t27 + (-t178 * t201 + t181 * t202) * t154 - t200 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226 - t161 + (-t121 + t247) * t248, t225 + t239 + (t119 + t242) * t248, -t119 ^ 2 - t121 ^ 2, t119 * t15 + t121 * t14 - t187 + t291, 0, 0, 0, 0, 0, t17 + t298, t16 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208 * t69, t208 ^ 2 - t69 ^ 2, t16 - t279, -t17 + t298, t123, t27 * t208 - g(1) * (-t164 * t38 - t165 * t63) - g(2) * (-t164 * t36 - t165 * t60) - g(3) * (t102 * t165 - t164 * t79) + t235 + t292 * t4, t27 * t69 - g(1) * (t164 * t63 - t165 * t38) - g(2) * (t164 * t60 - t165 * t36) - g(3) * (-t102 * t164 - t165 * t79) - t217 + t292 * t3;];
tau_reg  = t1;
