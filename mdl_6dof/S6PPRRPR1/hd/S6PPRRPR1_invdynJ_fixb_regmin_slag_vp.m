% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:31
% EndTime: 2019-03-08 18:47:44
% DurationCPUTime: 6.06s
% Computational Cost: add. (4051->430), mult. (10463->647), div. (0->0), fcn. (9799->18), ass. (0->209)
t163 = sin(pkin(7));
t274 = cos(pkin(6));
t230 = t163 * t274;
t284 = cos(qJ(3));
t164 = sin(pkin(6));
t162 = sin(pkin(12));
t169 = sin(qJ(3));
t166 = cos(pkin(12));
t273 = cos(pkin(7));
t227 = t166 * t273;
t189 = -t162 * t169 + t284 * t227;
t295 = t164 * t189;
t301 = t284 * t230 + t295;
t146 = qJD(1) * t274 + qJD(2);
t239 = t163 * t284;
t178 = -qJD(1) * t295 - t146 * t239;
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t217 = pkin(4) * t168 - qJ(5) * t171;
t108 = qJD(4) * t217 - t168 * qJD(5);
t190 = t284 * t162 + t169 * t227;
t184 = t190 * t164;
t262 = t163 * t169;
t69 = qJD(1) * t184 + t146 * t262;
t300 = t108 - t69;
t255 = qJD(3) * t171;
t148 = -qJD(6) + t255;
t161 = sin(pkin(13));
t165 = cos(pkin(13));
t249 = t165 * qJD(4);
t256 = qJD(3) * t168;
t123 = t161 * t256 - t249;
t254 = qJD(4) * t161;
t125 = t165 * t256 + t254;
t167 = sin(qJ(6));
t170 = cos(qJ(6));
t208 = t123 * t167 - t125 * t170;
t299 = t148 * t208;
t253 = qJD(4) * t168;
t244 = pkin(9) * t253;
t264 = t161 * t171;
t280 = t161 * t244 + t300 * t165 - t178 * t264;
t260 = t165 * t171;
t297 = t300 * t161 + t178 * t260;
t66 = qJD(3) * pkin(9) + t69;
t242 = t164 * t166 * t163;
t97 = -qJD(1) * t242 + t146 * t273;
t296 = -t168 * t66 + t171 * t97;
t144 = t274 * qJDD(1) + qJDD(2);
t236 = qJD(3) * t262;
t248 = qJD(1) * qJD(3);
t176 = -qJDD(1) * t295 - t144 * t239 + t146 * t236 + t248 * t184;
t272 = cos(pkin(11));
t219 = t274 * t272;
t271 = sin(pkin(11));
t109 = t162 * t219 + t166 * t271;
t182 = t162 * t271 - t166 * t219;
t229 = t164 * t272;
t293 = t163 * t229 + t182 * t273;
t50 = t109 * t169 + t293 * t284;
t218 = t274 * t271;
t110 = -t162 * t218 + t166 * t272;
t183 = t162 * t272 + t166 * t218;
t228 = t164 * t271;
t292 = -t163 * t228 + t183 * t273;
t52 = t110 * t169 + t292 * t284;
t199 = g(1) * t52 + g(2) * t50 - g(3) * t301;
t294 = t69 * qJD(3) - t176 + t199;
t156 = t171 * qJDD(3);
t247 = qJD(3) * qJD(4);
t232 = t168 * t247;
t291 = t232 - t156;
t290 = -t148 - qJD(6);
t289 = qJDD(1) * t190;
t288 = -qJDD(4) * pkin(4) + qJDD(5);
t151 = t165 * qJDD(4);
t245 = t168 * qJDD(3);
t193 = t171 * t247 + t245;
t95 = t161 * t193 - t151;
t246 = qJDD(4) * t161;
t96 = t165 * t193 + t246;
t35 = -qJD(6) * t208 + t167 * t96 + t170 * t95;
t94 = -qJDD(1) * t242 + t144 * t273;
t276 = t168 * t94;
t235 = qJD(3) * t284;
t39 = qJDD(3) * pkin(9) + (t144 * t169 + t146 * t235) * t163 + (t189 * t248 + t289) * t164;
t10 = qJDD(4) * qJ(5) + t276 + t171 * t39 + (qJD(5) + t296) * qJD(4);
t204 = pkin(4) * t171 + qJ(5) * t168 + pkin(3);
t20 = qJD(3) * t108 - qJDD(3) * t204 + t176;
t5 = t165 * t10 + t161 * t20;
t283 = pkin(10) + qJ(5);
t205 = pkin(5) * t168 - pkin(10) * t260;
t282 = -qJD(4) * t205 - t280;
t261 = t165 * t168;
t281 = (-pkin(9) * t261 - pkin(10) * t264) * qJD(4) + t297;
t47 = t168 * t97 + t171 * t66;
t45 = qJD(4) * qJ(5) + t47;
t56 = -qJD(3) * t204 + t178;
t17 = t161 * t56 + t165 * t45;
t279 = -t165 * t244 + t297;
t278 = qJD(3) * pkin(3);
t268 = t125 * t167;
t72 = t170 * t123 + t268;
t277 = t148 * t72;
t130 = t217 * qJD(3);
t32 = t161 * t130 + t165 * t296;
t158 = pkin(13) + qJ(6);
t154 = sin(t158);
t267 = t154 * t171;
t155 = cos(t158);
t266 = t155 * t171;
t265 = t161 * t167;
t128 = -t170 * t165 + t265;
t194 = t128 * t171;
t259 = qJD(3) * t194 - t128 * qJD(6);
t129 = t161 * t170 + t165 * t167;
t195 = t129 * t171;
t258 = -qJD(3) * t195 + t129 * qJD(6);
t99 = pkin(9) * t260 - t161 * t204;
t159 = t168 ^ 2;
t257 = -t171 ^ 2 + t159;
t252 = qJD(4) * t171;
t251 = qJD(6) * t168;
t250 = qJD(6) * t170;
t241 = pkin(5) * t161 + pkin(9);
t4 = -t10 * t161 + t165 * t20;
t2 = t291 * pkin(5) - pkin(10) * t96 + t4;
t3 = -pkin(10) * t95 + t5;
t240 = -t167 * t3 + t170 * t2;
t237 = t161 * t255;
t234 = qJ(5) * t156;
t16 = -t161 * t45 + t165 * t56;
t31 = t165 * t130 - t161 * t296;
t225 = t163 * t235;
t224 = qJD(4) * t235;
t220 = t167 * t2 + t170 * t3;
t12 = -pkin(5) * t255 - pkin(10) * t125 + t16;
t13 = -pkin(10) * t123 + t17;
t6 = t12 * t170 - t13 * t167;
t7 = t12 * t167 + t13 * t170;
t187 = t273 * t274 - t242;
t79 = t169 * t230 + t184;
t55 = t168 * t187 + t79 * t171;
t25 = -t161 * t55 - t165 * t301;
t26 = -t161 * t301 + t165 * t55;
t216 = -t167 * t26 + t170 * t25;
t215 = t167 * t25 + t170 * t26;
t122 = t165 * t204;
t77 = -pkin(10) * t261 - t122 + (-pkin(9) * t161 - pkin(5)) * t171;
t86 = -pkin(10) * t161 * t168 + t99;
t214 = -t167 * t86 + t170 * t77;
t213 = t167 * t77 + t170 * t86;
t115 = t168 * t273 + t171 * t262;
t82 = -t161 * t115 - t165 * t239;
t83 = t165 * t115 - t161 * t239;
t212 = -t167 * t83 + t170 * t82;
t211 = t167 * t82 + t170 * t83;
t209 = t168 * t39 - t171 * t94 + t66 * t252 + t97 * t253;
t140 = t283 * t161;
t203 = pkin(10) * t237 + qJD(5) * t165 - qJD(6) * t140 - t32;
t141 = t283 * t165;
t202 = qJD(3) * t205 + qJD(5) * t161 + qJD(6) * t141 + t31;
t174 = t163 * t182 - t229 * t273;
t175 = t163 * t183 + t228 * t273;
t51 = t109 * t284 - t293 * t169;
t53 = t110 * t284 - t292 * t169;
t54 = t79 * t168 - t171 * t187;
t201 = g(1) * (t168 * t53 - t171 * t175) + g(2) * (t168 * t51 - t171 * t174) + g(3) * t54;
t28 = t168 * t174 + t51 * t171;
t30 = t168 * t175 + t53 * t171;
t200 = g(1) * t30 + g(2) * t28 + g(3) * t55;
t198 = g(1) * t53 + g(2) * t51 + g(3) * t79;
t34 = -qJD(6) * t268 - t123 * t250 - t167 * t95 + t170 * t96;
t173 = qJD(3) ^ 2;
t196 = t284 * qJDD(3) - t169 * t173;
t43 = -qJD(4) * pkin(4) + qJD(5) - t296;
t11 = t209 + t288;
t192 = -t11 + t201;
t114 = t168 * t262 - t171 * t273;
t191 = -qJ(5) * t253 + (qJD(5) - t43) * t171;
t65 = t178 - t278;
t188 = -pkin(9) * qJDD(4) + (t65 - t178 - t278) * qJD(4);
t185 = -g(1) * t228 + g(2) * t229 - g(3) * t274;
t181 = t201 - t209;
t172 = qJD(4) ^ 2;
t177 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t172 + t294;
t150 = -pkin(5) * t165 - pkin(4);
t134 = t241 * t168;
t127 = qJDD(6) + t291;
t119 = t241 * t252;
t106 = t128 * t168;
t105 = t129 * t168;
t98 = -pkin(9) * t264 - t122;
t85 = qJD(4) * t115 + t168 * t225;
t84 = -qJD(4) * t114 + t171 * t225;
t71 = t79 * qJD(3);
t70 = t301 * qJD(3);
t64 = qJD(4) * t195 + t250 * t261 - t251 * t265;
t63 = -qJD(4) * t194 - t129 * t251;
t62 = t161 * t236 + t165 * t84;
t61 = -t161 * t84 + t165 * t236;
t41 = pkin(5) * t237 + t47;
t33 = pkin(5) * t123 + t43;
t24 = -qJD(4) * t54 + t70 * t171;
t23 = qJD(4) * t55 + t70 * t168;
t15 = t161 * t71 + t165 * t24;
t14 = -t161 * t24 + t165 * t71;
t8 = pkin(5) * t95 + t11;
t1 = [qJDD(1) - g(3), t144 * t274 - g(3) + (t162 ^ 2 + t166 ^ 2) * t164 ^ 2 * qJDD(1), 0, -qJD(3) * t71 + qJDD(3) * t301, -qJD(3) * t70 - qJDD(3) * t79, 0, 0, 0, 0, 0, t301 * t156 - qJD(4) * t23 - qJDD(4) * t54 + (-t171 * t71 - t253 * t301) * qJD(3), -t301 * t245 - qJD(4) * t24 - qJDD(4) * t55 + (t168 * t71 - t252 * t301) * qJD(3), -t25 * t156 + t123 * t23 + t54 * t95 + (-t14 * t171 + t25 * t253) * qJD(3), t26 * t156 + t125 * t23 + t54 * t96 + (t15 * t171 - t253 * t26) * qJD(3), -t123 * t15 - t125 * t14 - t25 * t96 - t26 * t95, t11 * t54 + t14 * t16 + t15 * t17 + t23 * t43 + t25 * t4 + t26 * t5 - g(3), 0, 0, 0, 0, 0 -(-qJD(6) * t215 + t14 * t170 - t15 * t167) * t148 + t216 * t127 + t23 * t72 + t54 * t35 (qJD(6) * t216 + t14 * t167 + t15 * t170) * t148 - t215 * t127 - t23 * t208 + t54 * t34; 0, t185 + t144, 0, t196 * t163 (-qJDD(3) * t169 - t284 * t173) * t163, 0, 0, 0, 0, 0, -t85 * qJD(4) - t114 * qJDD(4) + (-t168 * t224 + t171 * t196) * t163, -t84 * qJD(4) - t115 * qJDD(4) + (-t168 * t196 - t171 * t224) * t163, -t82 * t156 + t114 * t95 + t123 * t85 + (-t171 * t61 + t253 * t82) * qJD(3), t83 * t156 + t114 * t96 + t125 * t85 + (t171 * t62 - t253 * t83) * qJD(3), -t123 * t62 - t125 * t61 - t82 * t96 - t83 * t95, t11 * t114 + t16 * t61 + t17 * t62 + t4 * t82 + t43 * t85 + t5 * t83 + t185, 0, 0, 0, 0, 0 -(-qJD(6) * t211 - t167 * t62 + t170 * t61) * t148 + t212 * t127 + t85 * t72 + t114 * t35 (qJD(6) * t212 + t167 * t61 + t170 * t62) * t148 - t211 * t127 - t85 * t208 + t114 * t34; 0, 0, qJDD(3), t294, -t144 * t262 - t164 * t289 + t198, qJDD(3) * t159 + 0.2e1 * t171 * t232, 0.2e1 * t156 * t168 - 0.2e1 * t247 * t257, qJDD(4) * t168 + t171 * t172, qJDD(4) * t171 - t168 * t172, 0, t168 * t188 + t171 * t177, -t168 * t177 + t171 * t188, -t198 * t161 + (pkin(9) * t95 + t11 * t161 + t178 * t123 + (qJD(3) * t98 + t16) * qJD(4)) * t168 + (-t98 * qJDD(3) - t4 + (pkin(9) * t123 + t161 * t43) * qJD(4) - t280 * qJD(3) + t199 * t165) * t171, -t198 * t165 + (pkin(9) * t96 + t11 * t165 + t178 * t125 + (-qJD(3) * t99 - t17) * qJD(4)) * t168 + (t99 * qJDD(3) + t5 + (pkin(9) * t125 + t165 * t43) * qJD(4) + t279 * qJD(3) - t199 * t161) * t171, -t95 * t99 - t96 * t98 - t280 * t125 - t279 * t123 + (-t16 * t165 - t161 * t17) * t252 + (-t161 * t5 - t165 * t4 + t199) * t168, t43 * t168 * t178 + t4 * t98 + t5 * t99 + t279 * t17 + t280 * t16 + (t11 * t168 + t252 * t43 - t198) * pkin(9) + t199 * t204, -t106 * t34 - t208 * t63, -t105 * t34 + t106 * t35 + t208 * t64 - t63 * t72, -t106 * t127 - t148 * t63 - t171 * t34 - t208 * t253, -t105 * t127 + t148 * t64 + t171 * t35 - t253 * t72, -t127 * t171 - t148 * t253, t214 * t127 - t240 * t171 + t119 * t72 + t134 * t35 + t8 * t105 + t33 * t64 - g(1) * (t154 * t53 - t266 * t52) - g(2) * (t154 * t51 - t266 * t50) - g(3) * (t154 * t79 + t266 * t301) + (t6 * qJD(4) + t178 * t72) * t168 + (t281 * t167 + t282 * t170) * t148 + (t148 * t213 + t171 * t7) * qJD(6), -t213 * t127 + t220 * t171 - t119 * t208 + t134 * t34 - t8 * t106 + t33 * t63 - g(1) * (t155 * t53 + t267 * t52) - g(2) * (t155 * t51 + t267 * t50) - g(3) * (t155 * t79 - t267 * t301) + (-t7 * qJD(4) - t178 * t208) * t168 + (-t282 * t167 + t281 * t170) * t148 + (t148 * t214 + t171 * t6) * qJD(6); 0, 0, 0, 0, 0, -t171 * t173 * t168, t257 * t173, t245, t156, qJDD(4), qJD(4) * t47 - t256 * t65 + t181, -t276 + (-qJD(3) * t65 - t39) * t171 + t200, t161 * t234 - pkin(4) * t95 - t123 * t47 + t192 * t165 + (-t16 * t168 + t161 * t191 + t171 * t31) * qJD(3), t165 * t234 - pkin(4) * t96 - t125 * t47 - t192 * t161 + (t165 * t191 + t17 * t168 - t171 * t32) * qJD(3), t123 * t32 + t125 * t31 + (-qJ(5) * t95 - qJD(5) * t123 + t16 * t255 + t5) * t165 + (qJ(5) * t96 + qJD(5) * t125 + t17 * t255 - t4) * t161 - t200, -t16 * t31 - t17 * t32 - t43 * t47 + (-t16 * t161 + t165 * t17) * qJD(5) + t192 * pkin(4) + (-t161 * t4 + t165 * t5 - t200) * qJ(5), t129 * t34 - t208 * t259, -t128 * t34 - t129 * t35 + t208 * t258 - t259 * t72, t127 * t129 - t148 * t259 + t208 * t256, -t127 * t128 + t148 * t258 + t256 * t72, t148 * t256 (-t140 * t170 - t141 * t167) * t127 + t150 * t35 + t8 * t128 - t6 * t256 - t41 * t72 + t258 * t33 + (t167 * t203 + t170 * t202) * t148 + t201 * t155 -(-t140 * t167 + t141 * t170) * t127 + t150 * t34 + t8 * t129 + t7 * t256 + t41 * t208 + t259 * t33 + (-t167 * t202 + t170 * t203) * t148 - t201 * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161 * t245 - t151 + (-t125 + t254) * t255, t165 * t245 + t246 + (t123 + t249) * t255, -t123 ^ 2 - t125 ^ 2, t123 * t17 + t125 * t16 - t181 + t288, 0, 0, 0, 0, 0, t35 + t299, t34 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208 * t72, t208 ^ 2 - t72 ^ 2, t34 - t277, -t35 + t299, t127, t33 * t208 - g(1) * (-t154 * t30 + t155 * t52) - g(2) * (-t154 * t28 + t155 * t50) - g(3) * (-t154 * t55 - t155 * t301) + t240 + t290 * t7, t33 * t72 - g(1) * (-t154 * t52 - t155 * t30) - g(2) * (-t154 * t50 - t155 * t28) - g(3) * (t154 * t301 - t155 * t55) - t220 + t290 * t6;];
tau_reg  = t1;
