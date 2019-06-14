% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPPR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 02:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:46:37
% EndTime: 2019-05-05 02:46:52
% DurationCPUTime: 7.01s
% Computational Cost: add. (12365->393), mult. (27442->534), div. (0->0), fcn. (19561->12), ass. (0->221)
t207 = sin(pkin(6));
t209 = cos(pkin(6));
t213 = sin(qJ(3));
t214 = sin(qJ(2));
t216 = cos(qJ(3));
t217 = cos(qJ(2));
t249 = qJD(2) * qJD(3);
t241 = t216 * t249;
t247 = t213 * qJDD(2);
t182 = t241 + t247;
t206 = sin(pkin(11));
t208 = cos(pkin(11));
t242 = t213 * t249;
t246 = t216 * qJDD(2);
t230 = -t242 + t246;
t141 = t208 * t182 + t206 * t230;
t253 = qJD(2) * t213;
t174 = -t208 * t216 * qJD(2) + t206 * t253;
t165 = qJD(3) * t174;
t291 = t165 - t141;
t218 = qJD(3) ^ 2;
t256 = t213 * t208;
t176 = (t216 * t206 + t256) * qJD(2);
t281 = t176 ^ 2;
t157 = t281 + t218;
t264 = t176 * t174;
t303 = qJDD(3) + t264;
t322 = t206 * t303;
t97 = t208 * t157 + t322;
t320 = t208 * t303;
t99 = -t206 * t157 + t320;
t61 = t213 * t97 - t216 * t99;
t348 = t207 * (t214 * t61 + t217 * t291) - t209 * (t213 * t99 + t216 * t97);
t347 = pkin(8) * t61;
t344 = pkin(3) * t97;
t343 = qJ(4) * t97;
t342 = qJ(4) * t99;
t282 = t174 ^ 2;
t283 = -t282 - t281;
t140 = t206 * t182 - t208 * t230;
t251 = t176 * qJD(3);
t106 = t140 - t251;
t285 = t165 + t141;
t323 = -t208 * t106 + t206 * t285;
t324 = -t206 * t106 - t208 * t285;
t331 = -t213 * t324 + t216 * t323;
t341 = -pkin(2) * t283 + pkin(8) * t331;
t158 = t281 - t218;
t304 = qJDD(3) - t264;
t319 = t208 * t304;
t321 = t206 * t304;
t339 = t213 * (t206 * t158 + t319) - t216 * (t208 * t158 - t321);
t338 = t209 * (t213 * t323 + t216 * t324) + (t214 * t331 - t217 * t283) * t207;
t337 = pkin(3) * t324;
t336 = qJ(4) * t324;
t333 = -pkin(3) * t283 + qJ(4) * t323;
t305 = t140 + t251;
t121 = -t218 - t282;
t85 = t206 * t121 + t319;
t88 = -t208 * t121 + t321;
t52 = t213 * t85 + t216 * t88;
t332 = (t214 * t52 + t217 * t305) * t207 + t209 * (t213 * t88 - t216 * t85);
t153 = t282 - t218;
t330 = t213 * (-t208 * t153 + t322) - t216 * (t206 * t153 + t320);
t328 = pkin(8) * t52;
t316 = t291 * qJ(5);
t315 = t213 * (t206 * t291 - t208 * t305) + t216 * (-t206 * t305 - t208 * t291);
t314 = pkin(3) * t85;
t313 = qJ(4) * t85;
t312 = qJ(4) * t88;
t212 = sin(qJ(6));
t215 = cos(qJ(6));
t145 = t212 * qJD(3) - t215 * t174;
t147 = t215 * qJD(3) + t212 * t174;
t104 = t147 * t145;
t134 = qJDD(6) + t141;
t287 = -t104 + t134;
t295 = t212 * t287;
t294 = t215 * t287;
t267 = sin(pkin(10));
t268 = cos(pkin(10));
t229 = g(1) * t267 - g(2) * t268;
t254 = -g(3) + qJDD(1);
t290 = t207 * t254 + t209 * t229;
t289 = 2 * qJD(4);
t169 = qJD(6) + t176;
t117 = t169 * t145;
t95 = -t145 * qJD(6) + t215 * qJDD(3) + t212 * t140;
t288 = -t117 + t95;
t219 = qJD(2) ^ 2;
t194 = t213 * t219 * t216;
t188 = qJDD(3) + t194;
t186 = -g(1) * t268 - g(2) * t267;
t119 = t217 * t186 + t290 * t214;
t114 = -t219 * pkin(2) + qJDD(2) * pkin(8) + t119;
t152 = -t207 * t229 + t209 * t254;
t91 = t213 * t114 - t216 * t152;
t78 = (-t182 + t241) * qJ(4) + t188 * pkin(3) - t91;
t187 = qJD(3) * pkin(3) - qJ(4) * t253;
t203 = t216 ^ 2;
t199 = t203 * t219;
t92 = t216 * t114 + t213 * t152;
t79 = -pkin(3) * t199 + qJ(4) * t230 - qJD(3) * t187 + t92;
t44 = -0.2e1 * qJD(4) * t174 + t206 * t78 + t208 * t79;
t284 = t281 - t282;
t143 = t145 ^ 2;
t144 = t147 ^ 2;
t167 = t169 ^ 2;
t280 = 2 * qJD(5);
t279 = -pkin(4) - pkin(9);
t278 = pkin(4) * t206;
t277 = pkin(4) * t208;
t234 = t214 * t186 - t290 * t217;
t113 = -qJDD(2) * pkin(2) - t219 * pkin(8) + t234;
t84 = -t230 * pkin(3) - qJ(4) * t199 + t187 * t253 + qJDD(4) + t113;
t275 = t206 * t84;
t274 = t208 * t84;
t151 = t176 * pkin(5) - qJD(3) * pkin(9);
t120 = t174 * pkin(4) - t176 * qJ(5);
t232 = -t218 * pkin(4) - t174 * t120 + t44;
t248 = qJDD(3) * qJ(5);
t28 = t248 - t140 * pkin(5) - t282 * pkin(9) + (t280 + t151) * qJD(3) + t232;
t273 = t212 * t28;
t83 = t104 + t134;
t272 = t212 * t83;
t239 = t206 * t79 - t208 * t78;
t43 = t176 * t289 + t239;
t23 = t206 * t44 - t208 * t43;
t271 = t213 * t23;
t270 = t215 * t28;
t269 = t215 * t83;
t266 = t169 * t212;
t265 = t169 * t215;
t257 = t213 * t188;
t189 = qJDD(3) - t194;
t255 = t216 * t189;
t250 = t289 + t120;
t245 = -t144 - t167;
t244 = t206 * t104;
t243 = t208 * t104;
t240 = qJ(5) * t206 + pkin(3);
t24 = t206 * t43 + t208 * t44;
t231 = -qJDD(3) * pkin(4) - t218 * qJ(5) + qJDD(5) + t239;
t27 = -qJDD(3) * pkin(9) + t285 * pkin(5) + (pkin(9) * t174 + t250) * t176 + t231;
t221 = t140 * pkin(4) + t316 + t84;
t237 = pkin(4) * qJD(3) - (2 * qJD(5));
t35 = t221 + (-t151 + t237) * t176 + t140 * pkin(9) - t282 * pkin(5);
t16 = t212 * t35 - t215 * t27;
t55 = t213 * t91 + t216 * t92;
t236 = t212 * qJDD(3) - t215 * t140;
t17 = t212 * t27 + t215 * t35;
t7 = -t215 * t16 + t212 * t17;
t8 = t212 * t16 + t215 * t17;
t226 = qJD(3) * t280 + t232;
t33 = t226 + t248;
t34 = t176 * t250 + t231;
t19 = t206 * t33 - t208 * t34;
t20 = t206 * t34 + t208 * t33;
t9 = -t213 * t19 + t216 * t20;
t183 = -0.2e1 * t242 + t246;
t228 = (-qJD(6) + t169) * t147 - t236;
t225 = t213 * (t208 * t141 - t206 * t251) + t216 * (t206 * t141 + t208 * t251);
t224 = t213 * (t206 * t140 + t165 * t208) + t216 * (-t208 * t140 + t165 * t206);
t222 = (t213 * (-t174 * t208 + t176 * t206) + t216 * (-t174 * t206 - t176 * t208)) * qJD(3);
t220 = -t176 * t280 + t221;
t202 = t213 ^ 2;
t198 = t202 * t219;
t193 = -t199 - t218;
t192 = -t198 - t218;
t185 = t198 + t199;
t184 = (t202 + t203) * qJDD(2);
t181 = 0.2e1 * t241 + t247;
t149 = -t213 * t192 - t255;
t148 = t216 * t193 - t257;
t116 = -t144 + t167;
t115 = t143 - t167;
t101 = t144 - t143;
t94 = -t147 * qJD(6) - t236;
t93 = -t167 - t143;
t90 = -t143 - t144;
t80 = (t145 * t212 + t147 * t215) * t169;
t69 = t117 + t95;
t65 = (qJD(6) + t169) * t147 + t236;
t63 = -t147 * t265 - t212 * t95;
t62 = -t145 * t266 - t215 * t94;
t59 = -t212 * t115 - t269;
t58 = -t215 * t116 - t295;
t57 = -t212 * t245 - t269;
t56 = t215 * t245 - t272;
t54 = t215 * t93 - t295;
t53 = t212 * t93 + t294;
t50 = t176 * t237 + t221;
t47 = t212 * t69 + t215 * t228;
t46 = t212 * t65 - t215 * t288;
t45 = t212 * t228 - t215 * t69;
t41 = t220 + (t305 + t251) * pkin(4);
t40 = -pkin(4) * t251 - t220 - t316;
t39 = t206 * t56 + t208 * t288;
t38 = t206 * t288 - t208 * t56;
t37 = t206 * t53 + t208 * t65;
t36 = t206 * t65 - t208 * t53;
t32 = t206 * t45 + t208 * t90;
t31 = t206 * t90 - t208 * t45;
t30 = -qJ(5) * t283 + t34;
t29 = -pkin(4) * t283 + t33;
t25 = pkin(5) * t45 - qJ(5) * t47;
t22 = -t213 * t38 + t216 * t39;
t21 = -t213 * t36 + t216 * t37;
t18 = -t213 * t31 + t216 * t32;
t14 = pkin(5) * t288 + t279 * t57 - t273;
t13 = pkin(5) * t65 + t279 * t54 + t270;
t12 = t216 * t24 - t271;
t11 = pkin(5) * t56 - qJ(5) * t57 - t17;
t10 = pkin(5) * t53 - qJ(5) * t54 - t16;
t6 = t206 * t7 + t208 * t28;
t5 = t206 * t28 - t208 * t7;
t4 = pkin(5) * t90 + t279 * t47 - t8;
t3 = pkin(5) * t7 - qJ(5) * t8;
t2 = pkin(5) * t28 + t279 * t8;
t1 = -t213 * t5 + t216 * t6;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t254, 0, 0, 0, 0, 0, 0, (qJDD(2) * t217 - t214 * t219) * t207, (-qJDD(2) * t214 - t217 * t219) * t207, 0, t209 * t152 + (t119 * t214 - t217 * t234) * t207, 0, 0, 0, 0, 0, 0, t209 * (t216 * t188 + t213 * t193) + (t214 * t148 + t217 * t183) * t207, t209 * (-t213 * t189 + t216 * t192) + (t214 * t149 - t217 * t181) * t207, (t184 * t214 + t185 * t217) * t207, t209 * (t213 * t92 - t216 * t91) + (-t217 * t113 + t214 * t55) * t207, 0, 0, 0, 0, 0, 0, -t332, t348, t338, t209 * (t213 * t24 + t216 * t23) + (t214 * t12 - t217 * t84) * t207, 0, 0, 0, 0, 0, 0, t338, t332, -t348, t209 * (t216 * t19 + t213 * t20) + (t214 * t9 - t217 * t50) * t207, 0, 0, 0, 0, 0, 0, t209 * (t213 * t37 + t216 * t36) + (t214 * t21 - t217 * t54) * t207, t209 * (t213 * t39 + t216 * t38) + (t214 * t22 - t217 * t57) * t207, t209 * (t213 * t32 + t216 * t31) + (t214 * t18 - t217 * t47) * t207, t209 * (t213 * t6 + t216 * t5) + (t214 * t1 - t217 * t8) * t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t234, -t119, 0, 0, (t182 + t241) * t213, t216 * t181 + t213 * t183, t257 + t216 * (-t198 + t218), t183 * t216, t213 * (t199 - t218) + t255, 0, pkin(2) * t183 + pkin(8) * t148 - t216 * t113, -pkin(2) * t181 + pkin(8) * t149 + t213 * t113, pkin(2) * t185 + pkin(8) * t184 + t55, -pkin(2) * t113 + pkin(8) * t55, t225, t315, t339, t224, -t330, t222, t213 * (t275 - t313) + t216 * (-pkin(3) * t305 - t274 - t312) - pkin(2) * t305 - t328, t213 * (t274 + t343) + t216 * (pkin(3) * t291 + t275 - t342) + pkin(2) * t291 + t347, t213 * (-t23 - t336) + t216 * (t24 + t333) + t341, -qJ(4) * t271 + t216 * (-pkin(3) * t84 + qJ(4) * t24) - pkin(2) * t84 + pkin(8) * t12, t222, -t339, t330, t225, t315, t224, t213 * (-t206 * t29 + t208 * t30 - t336) + t216 * (t206 * t30 + t208 * t29 + t333) + t341, t213 * (-t206 * t41 + t313) + t216 * (t208 * t41 + t312) + t328 + (qJ(5) * t256 + t216 * t240 + pkin(2)) * t305, t213 * (t208 * t40 - t343) + t216 * (t206 * t40 + t342) - t347 - (-t213 * t278 + t216 * (pkin(3) + t277) + pkin(2)) * t291, (t213 * (-qJ(5) * t208 + t278) + t216 * (-t240 - t277) - pkin(2)) * t50 + (pkin(8) + qJ(4)) * t9, t213 * (-t206 * t63 + t243) + t216 * (t208 * t63 + t244), t213 * (t208 * t101 - t206 * t46) + t216 * (t206 * t101 + t208 * t46), t213 * (-t206 * t58 + t208 * t69) + t216 * (t206 * t69 + t208 * t58), t213 * (-t206 * t62 - t243) + t216 * (t208 * t62 - t244), t213 * (-t206 * t59 + t208 * t228) + t216 * (t206 * t228 + t208 * t59), t213 * (t208 * t134 - t206 * t80) + t216 * (t206 * t134 + t208 * t80), t213 * (-qJ(4) * t36 + t208 * t10 - t206 * t13) + t216 * (-pkin(3) * t54 + qJ(4) * t37 + t206 * t10 + t208 * t13) - pkin(2) * t54 + pkin(8) * t21, t213 * (-qJ(4) * t38 + t208 * t11 - t206 * t14) + t216 * (-pkin(3) * t57 + qJ(4) * t39 + t206 * t11 + t208 * t14) - pkin(2) * t57 + pkin(8) * t22, t213 * (-qJ(4) * t31 - t206 * t4 + t208 * t25) + t216 * (-pkin(3) * t47 + qJ(4) * t32 + t206 * t25 + t208 * t4) - pkin(2) * t47 + pkin(8) * t18, t213 * (-qJ(4) * t5 - t206 * t2 + t208 * t3) + t216 * (-pkin(3) * t8 + qJ(4) * t6 + t208 * t2 + t206 * t3) - pkin(2) * t8 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, t198 - t199, t247, t194, t246, qJDD(3), -t91, -t92, 0, 0, t264, t284, t285, -t264, -t106, qJDD(3), -t43 + t314, -t44 - t344, t337, pkin(3) * t23, qJDD(3), -t285, t106, t264, t284, -t264, -pkin(4) * t285 - qJ(5) * t106 + t337, -pkin(4) * t304 - qJ(5) * t121 - t314 + t34, t344 + pkin(4) * t157 + (qJDD(3) + t303) * qJ(5) + t226, pkin(3) * t19 - pkin(4) * t34 + qJ(5) * t33, -t147 * t266 + t215 * t95, -t212 * t288 - t215 * t65, -t212 * t116 + t294, t145 * t265 - t212 * t94, t215 * t115 - t272, (-t145 * t215 + t147 * t212) * t169, pkin(3) * t36 + qJ(5) * t65 + t279 * t53 + t273, pkin(3) * t38 + qJ(5) * t288 + t279 * t56 + t270, pkin(3) * t31 + qJ(5) * t90 + t279 * t45 - t7, pkin(3) * t5 + qJ(5) * t28 + t279 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, -t291, t283, t84, 0, 0, 0, 0, 0, 0, t283, -t305, t291, t50, 0, 0, 0, 0, 0, 0, t54, t57, t47, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, t304, -t157, t34, 0, 0, 0, 0, 0, 0, t53, t56, t45, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t101, t69, -t104, t228, t134, -t16, -t17, 0, 0;];
tauJ_reg  = t15;
