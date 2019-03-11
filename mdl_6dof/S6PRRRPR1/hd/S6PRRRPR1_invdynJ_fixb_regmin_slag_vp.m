% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:26
% EndTime: 2019-03-08 23:03:36
% DurationCPUTime: 4.34s
% Computational Cost: add. (6425->435), mult. (15446->611), div. (0->0), fcn. (12484->18), ass. (0->248)
t211 = sin(qJ(4));
t212 = sin(qJ(3));
t215 = cos(qJ(4));
t216 = cos(qJ(3));
t151 = t211 * t216 + t212 * t215;
t346 = pkin(8) + pkin(9);
t282 = qJD(3) * t346;
t156 = t212 * t282;
t157 = t216 * t282;
t206 = sin(pkin(6));
t217 = cos(qJ(2));
t309 = t206 * t217;
t281 = qJD(1) * t309;
t165 = t346 * t212;
t166 = t346 * t216;
t302 = -t211 * t165 + t215 * t166;
t362 = -qJD(4) * t302 + t151 * t281 + t156 * t211 - t215 * t157;
t150 = t211 * t212 - t215 * t216;
t294 = qJD(4) * t215;
t295 = qJD(4) * t211;
t361 = -t150 * t281 + t215 * t156 + t211 * t157 + t165 * t294 + t166 * t295;
t296 = qJD(2) * t216;
t275 = t215 * t296;
t298 = qJD(2) * t212;
t277 = t211 * t298;
t143 = -t275 + t277;
t145 = -t211 * t296 - t215 * t298;
t204 = sin(pkin(12));
t207 = cos(pkin(12));
t257 = -t207 * t143 + t145 * t204;
t353 = qJD(6) - t257;
t354 = t353 - qJD(6);
t210 = sin(qJ(6));
t200 = qJD(3) + qJD(4);
t290 = qJD(2) * qJD(3);
t272 = t216 * t290;
t287 = t216 * qJDD(2);
t288 = t212 * qJDD(2);
t73 = qJD(4) * t275 - t200 * t277 + t211 * t287 + (t272 + t288) * t215;
t112 = t200 * t151;
t249 = t211 * t288 - t215 * t287;
t74 = qJD(2) * t112 + t249;
t40 = -t204 * t73 - t207 * t74;
t38 = qJDD(6) - t40;
t324 = t210 * t38;
t214 = cos(qJ(6));
t357 = t214 * t353;
t360 = -t353 * t357 - t324;
t245 = -t143 * t204 - t207 * t145;
t75 = -t214 * t200 + t210 * t245;
t359 = t353 * t75;
t358 = -qJ(5) * t112 - qJD(5) * t150 - t361;
t111 = t200 * t150;
t356 = qJ(5) * t111 - qJD(5) * t151 + t362;
t205 = sin(pkin(11));
t208 = cos(pkin(11));
t209 = cos(pkin(6));
t213 = sin(qJ(2));
t307 = t209 * t213;
t137 = t205 * t217 + t208 * t307;
t139 = -t205 * t307 + t208 * t217;
t203 = qJ(3) + qJ(4);
t194 = pkin(12) + t203;
t184 = sin(t194);
t185 = cos(t194);
t198 = qJDD(3) + qJDD(4);
t299 = qJD(1) * t213;
t276 = t206 * t299;
t267 = qJD(2) * t346 + t276;
t300 = qJD(1) * t209;
t116 = t212 * t300 + t216 * t267;
t108 = t215 * t116;
t115 = -t267 * t212 + t216 * t300;
t327 = qJD(3) * pkin(3);
t109 = t115 + t327;
t246 = -t211 * t109 - t108;
t289 = t209 * qJDD(1);
t171 = t216 * t289;
t291 = qJD(1) * qJD(2);
t123 = qJDD(2) * pkin(8) + (qJDD(1) * t213 + t217 * t291) * t206;
t261 = pkin(9) * qJDD(2) + t123;
t64 = qJDD(3) * pkin(3) - qJD(3) * t116 - t261 * t212 + t171;
t65 = qJD(3) * t115 + t212 * t289 + t261 * t216;
t226 = qJD(4) * t246 - t211 * t65 + t215 * t64;
t15 = pkin(4) * t198 - qJ(5) * t73 + qJD(5) * t145 + t226;
t347 = -(qJD(4) * t109 + t65) * t215 + t116 * t295 - t211 * t64;
t17 = -qJ(5) * t74 - qJD(5) * t143 - t347;
t5 = t15 * t207 - t17 * t204;
t3 = -pkin(5) * t198 - t5;
t311 = t206 * t213;
t312 = t206 * t208;
t313 = t205 * t206;
t355 = t3 + g(3) * (-t184 * t311 + t185 * t209) + g(2) * (-t137 * t184 - t185 * t312) + g(1) * (-t139 * t184 + t185 * t313);
t197 = t216 * pkin(3);
t332 = pkin(2) + t197;
t293 = qJD(6) * t210;
t242 = -t214 * t38 + t293 * t353;
t134 = t145 * qJ(5);
t106 = t211 * t116;
t265 = t215 * t109 - t106;
t59 = t134 + t265;
t352 = pkin(4) * t74 + qJDD(5);
t195 = sin(t203);
t196 = cos(t203);
t351 = -g(3) * (-t195 * t311 + t196 * t209) - g(1) * (-t139 * t195 + t196 * t313) - g(2) * (-t137 * t195 - t196 * t312);
t104 = t207 * t150 + t151 * t204;
t105 = -t150 * t204 + t151 * t207;
t253 = g(1) * t139 + g(2) * t137;
t230 = -g(3) * t311 - t253;
t317 = qJ(5) * t143;
t60 = -t246 - t317;
t326 = t204 * t60;
t49 = pkin(4) * t200 + t59;
t26 = t207 * t49 - t326;
t24 = -pkin(5) * t200 - t26;
t330 = t204 * t356 + t207 * t358;
t135 = -qJD(2) * t332 - t281;
t95 = pkin(4) * t143 + qJD(5) + t135;
t35 = -pkin(5) * t257 - pkin(10) * t245 + t95;
t6 = t204 * t15 + t207 * t17;
t4 = pkin(10) * t198 + t6;
t255 = pkin(4) * t150 - t332;
t50 = pkin(5) * t104 - pkin(10) * t105 + t255;
t256 = -t215 * t165 - t166 * t211;
t237 = -qJ(5) * t151 + t256;
t84 = -qJ(5) * t150 + t302;
t54 = t204 * t237 + t207 * t84;
t69 = -t111 * t207 - t112 * t204;
t350 = -(qJD(6) * t35 + t4) * t104 + t24 * t69 + t3 * t105 + (-qJD(6) * t50 - t330) * t353 - t54 * t38 + t230;
t306 = t209 * t217;
t136 = t205 * t213 - t208 * t306;
t138 = t205 * t306 + t208 * t213;
t254 = g(1) * t138 + g(2) * t136;
t231 = g(3) * t309 - t254;
t349 = -t231 * t185 + t50 * t38;
t273 = t213 * t291;
t167 = t206 * t273;
t218 = qJD(3) ^ 2;
t270 = qJDD(1) * t309;
t348 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t218 + t206 * (-g(3) * t217 + t273) - t167 + t254 + t270;
t338 = t145 * pkin(4);
t337 = t24 * t257;
t335 = t75 * t245;
t77 = t200 * t210 + t214 * t245;
t334 = t77 * t245;
t333 = t353 * t245;
t331 = t204 * t358 - t207 * t356;
t51 = t207 * t60;
t27 = t204 * t49 + t51;
t329 = pkin(3) * qJD(4);
t328 = qJD(2) * pkin(2);
t292 = qJD(6) * t214;
t41 = -t204 * t74 + t207 * t73;
t31 = t210 * t198 + t200 * t292 + t214 * t41 - t245 * t293;
t325 = t210 * t31;
t323 = t210 * t353;
t322 = t214 * t77;
t320 = t24 * t105;
t258 = -t115 * t211 - t108;
t238 = t258 + t317;
t308 = t207 * t211;
t303 = t215 * t115 - t106;
t62 = t134 + t303;
t319 = -t204 * t62 + t207 * t238 + (t204 * t215 + t308) * t329;
t314 = t204 * t211;
t318 = -t204 * t238 - t207 * t62 + (t207 * t215 - t314) * t329;
t315 = t145 * t143;
t310 = t206 * t216;
t304 = qJDD(1) - g(3);
t190 = pkin(3) * t215 + pkin(4);
t133 = pkin(3) * t308 + t204 * t190;
t161 = pkin(4) * t196 + t197;
t201 = t212 ^ 2;
t301 = -t216 ^ 2 + t201;
t297 = qJD(2) * t213;
t193 = t212 * t327;
t192 = pkin(3) * t298;
t285 = t210 * t309;
t284 = t214 * t309;
t279 = t206 * t297;
t278 = qJD(2) * t309;
t25 = pkin(10) * t200 + t27;
t248 = t210 * t25 - t214 * t35;
t274 = t24 * t293 + t245 * t248;
t271 = t217 * t290;
t269 = pkin(4) * t112 + t193;
t264 = -t214 * t198 + t210 * t41;
t263 = t206 * t304;
t129 = pkin(10) + t133;
t47 = pkin(5) * t245 - pkin(10) * t257 - t338;
t259 = qJD(6) * t129 + t192 + t47;
t252 = g(1) * t205 - g(2) * t208;
t68 = -t111 * t204 + t207 * t112;
t251 = pkin(5) * t68 - pkin(10) * t69 + t269 - t276;
t250 = t245 * t27 + t257 * t26;
t12 = t210 * t35 + t214 * t25;
t141 = t209 * t216 - t212 * t311;
t142 = t209 * t212 + t213 * t310;
t85 = t141 * t215 - t142 * t211;
t86 = t141 * t211 + t142 * t215;
t244 = t257 * t323 - t242;
t219 = qJD(2) ^ 2;
t243 = qJDD(2) * t217 - t213 * t219;
t132 = -pkin(3) * t314 + t190 * t207;
t57 = t204 * t85 + t207 * t86;
t241 = -t210 * t57 - t284;
t240 = -t214 * t57 + t285;
t236 = t12 * t245 + t210 * t355 + t24 * t292;
t234 = -t276 + t193;
t232 = qJD(3) * t192 - qJDD(2) * t332 + t167;
t159 = -t281 - t328;
t229 = -qJD(2) * t159 - t123 + t253;
t228 = -t129 * t38 - t318 * t353 - t337;
t98 = t232 - t270;
t224 = -pkin(8) * qJDD(3) + (t159 + t281 - t328) * qJD(3);
t55 = t98 + t352;
t222 = -g(1) * (-t139 * t196 - t195 * t313) - g(2) * (-t137 * t196 + t195 * t312) - g(3) * (-t195 * t209 - t196 * t311) + t135 * t143 + t347;
t221 = t135 * t145 + t226 + t351;
t199 = -qJ(5) - t346;
t187 = -pkin(4) * t207 - pkin(5);
t186 = pkin(4) * t204 + pkin(10);
t160 = -pkin(3) * t212 - pkin(4) * t195;
t155 = pkin(2) + t161;
t128 = -pkin(5) - t132;
t121 = t184 * t209 + t185 * t311;
t114 = -qJD(3) * t142 - t212 * t278;
t113 = qJD(3) * t141 + t216 * t278;
t94 = t139 * t185 + t184 * t313;
t92 = t137 * t185 - t184 * t312;
t78 = -t143 ^ 2 + t145 ^ 2;
t67 = -t249 + (-qJD(2) * t151 - t145) * t200;
t66 = t143 * t200 + t73;
t56 = t204 * t86 - t207 * t85;
t53 = t204 * t84 - t207 * t237;
t44 = -qJD(4) * t86 - t113 * t211 + t114 * t215;
t43 = qJD(4) * t85 + t113 * t215 + t114 * t211;
t32 = qJD(6) * t77 + t264;
t29 = t207 * t59 - t326;
t28 = t204 * t59 + t51;
t19 = t204 * t44 + t207 * t43;
t18 = t204 * t43 - t207 * t44;
t13 = -pkin(5) * t40 - pkin(10) * t41 + t55;
t10 = t214 * t13;
t9 = t357 * t77 + t325;
t8 = -t334 - t360;
t7 = t244 + t335;
t1 = (t31 - t359) * t214 + (-t353 * t77 - t32) * t210;
t2 = [t304, 0, t243 * t206 (-qJDD(2) * t213 - t217 * t219) * t206, 0, 0, 0, 0, 0, qJD(3) * t114 + qJDD(3) * t141 + (-t212 * t271 + t216 * t243) * t206, -qJD(3) * t113 - qJDD(3) * t142 + (-t212 * t243 - t216 * t271) * t206, 0, 0, 0, 0, 0, t198 * t85 + t200 * t44 + (t143 * t297 - t217 * t74) * t206, -t198 * t86 - t200 * t43 + (-t145 * t297 - t217 * t73) * t206, t18 * t245 + t19 * t257 + t40 * t57 + t41 * t56, -t18 * t26 + t19 * t27 - t5 * t56 + t57 * t6 - g(3) + (-t217 * t55 + t297 * t95) * t206, 0, 0, 0, 0, 0 (qJD(6) * t240 - t19 * t210 + t214 * t279) * t353 + t241 * t38 + t18 * t75 + t56 * t32 -(qJD(6) * t241 + t19 * t214 + t210 * t279) * t353 + t240 * t38 + t18 * t77 + t56 * t31; 0, qJDD(2), t304 * t309 + t254, -t213 * t263 + t253, qJDD(2) * t201 + 0.2e1 * t212 * t272, 0.2e1 * t212 * t287 - 0.2e1 * t301 * t290, qJDD(3) * t212 + t216 * t218, qJDD(3) * t216 - t212 * t218, 0, t224 * t212 + t216 * t348, -t212 * t348 + t224 * t216, t111 * t145 + t151 * t73, t111 * t143 + t112 * t145 - t150 * t73 - t151 * t74, -t111 * t200 + t151 * t198, -t112 * t200 - t150 * t198, 0, t135 * t112 + t234 * t143 + t98 * t150 - t231 * t196 + t256 * t198 + t200 * t362 - t332 * t74, -t135 * t111 - t234 * t145 + t98 * t151 + t231 * t195 - t302 * t198 + t200 * t361 - t332 * t73, -t104 * t6 - t105 * t5 + t245 * t331 + t257 * t330 - t26 * t69 - t27 * t68 + t40 * t54 + t41 * t53 + t230, t6 * t54 - t5 * t53 + t55 * t255 + t95 * t269 - g(1) * (-t138 * t155 - t139 * t199) - g(2) * (-t136 * t155 - t137 * t199) + t330 * t27 - t331 * t26 + (-t95 * t299 - g(3) * (t155 * t217 - t199 * t213)) * t206, t69 * t322 + (t214 * t31 - t293 * t77) * t105 (-t210 * t77 - t214 * t75) * t69 + (-t325 - t214 * t32 + (t210 * t75 - t322) * qJD(6)) * t105, t104 * t31 - t105 * t242 + t357 * t69 + t68 * t77, -t69 * t323 - t104 * t32 - t68 * t75 + (-t292 * t353 - t324) * t105, t104 * t38 + t353 * t68, t10 * t104 - t248 * t68 + t53 * t32 + t331 * t75 + (t251 * t353 + (-t25 * t104 - t353 * t54 + t320) * qJD(6) + t349) * t214 + t350 * t210, -t12 * t68 + t53 * t31 + t331 * t77 + (-(-qJD(6) * t25 + t13) * t104 - qJD(6) * t320 + (qJD(6) * t54 - t251) * t353 - t349) * t210 + t350 * t214; 0, 0, 0, 0, -t212 * t219 * t216, t301 * t219, t288, t287, qJDD(3), -g(3) * t141 + t212 * t229 - t252 * t310 + t171, g(3) * t142 + (t206 * t252 - t289) * t212 + t229 * t216, -t315, t78, t66, t67, t198, -t258 * t200 + (-t143 * t298 + t215 * t198 - t200 * t295) * pkin(3) + t221, t303 * t200 + (t145 * t298 - t211 * t198 - t200 * t294) * pkin(3) + t222, -t132 * t41 + t133 * t40 + t245 * t319 + t257 * t318 + t250, t6 * t133 + t5 * t132 - t95 * (t192 - t338) - g(1) * (t139 * t160 + t161 * t313) - g(2) * (t137 * t160 - t161 * t312) - g(3) * (t160 * t311 + t161 * t209) + t318 * t27 - t319 * t26, t9, t1, t8, t7, -t333, t128 * t32 + t319 * t75 + t228 * t210 + (-t259 * t353 - t355) * t214 + t274, t128 * t31 + t214 * t228 + t259 * t323 + t319 * t77 + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t315, t78, t66, t67, t198, -t200 * t246 + t221, t200 * t265 + t222, -t28 * t245 - t29 * t257 + (t204 * t40 - t207 * t41) * pkin(4) + t250, t26 * t28 - t27 * t29 + (t95 * t145 + t6 * t204 + t5 * t207 + t351) * pkin(4), t9, t1, t8, t7, -t333, t187 * t32 - t28 * t75 + (-t186 * t38 + t29 * t353 - t337) * t210 + ((-qJD(6) * t186 - t47) * t353 - t355) * t214 + t274, t187 * t31 + (t210 * t47 + t214 * t29) * t353 - t28 * t77 - t214 * t337 + t242 * t186 + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245 ^ 2 - t257 ^ 2, -t217 * t263 + t245 * t26 - t257 * t27 + t232 - t254 + t352, 0, 0, 0, 0, 0, t244 - t335, -t334 + t360; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75, -t75 ^ 2 + t77 ^ 2, t31 + t359, t354 * t77 - t264, t38, -t210 * t4 + t10 - t24 * t77 - g(1) * (t138 * t214 - t210 * t94) - g(2) * (t136 * t214 - t210 * t92) - g(3) * (-t121 * t210 - t284) + t354 * t12, -t214 * t4 - t210 * t13 + t24 * t75 - g(1) * (-t138 * t210 - t214 * t94) - g(2) * (-t136 * t210 - t214 * t92) - g(3) * (-t121 * t214 + t285) - t354 * t248;];
tau_reg  = t2;
