% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:05
% EndTime: 2019-12-31 22:06:20
% DurationCPUTime: 5.89s
% Computational Cost: add. (14825->376), mult. (29801->467), div. (0->0), fcn. (20690->8), ass. (0->249)
t229 = sin(qJ(2));
t233 = cos(qJ(2));
t228 = sin(qJ(3));
t232 = cos(qJ(3));
t278 = qJD(1) * t229;
t198 = -t232 * qJD(2) + t228 * t278;
t218 = t229 * qJDD(1);
t274 = qJD(1) * qJD(2);
t269 = t233 * t274;
t203 = t218 + t269;
t248 = -t228 * qJDD(2) - t232 * t203;
t170 = -t198 * qJD(3) - t248;
t200 = t228 * qJD(2) + t232 * t278;
t227 = sin(qJ(4));
t231 = cos(qJ(4));
t175 = t231 * t198 + t227 * t200;
t249 = t232 * qJDD(2) - t228 * t203;
t240 = -t200 * qJD(3) + t249;
t106 = -t175 * qJD(4) + t231 * t170 + t227 * t240;
t215 = t233 * qJD(1) - qJD(3);
t211 = -qJD(4) + t215;
t295 = t175 * t211;
t330 = t106 + t295;
t177 = -t227 * t198 + t231 * t200;
t136 = t177 * t175;
t217 = t229 * t274;
t273 = t233 * qJDD(1);
t204 = -t217 + t273;
t197 = -qJDD(3) + t204;
t194 = -qJDD(4) + t197;
t329 = -t136 + t194;
t288 = t227 * t329;
t174 = t177 ^ 2;
t315 = t211 ^ 2;
t327 = -t174 - t315;
t65 = -t231 * t327 - t288;
t282 = t231 * t329;
t89 = -t227 * t327 + t282;
t53 = t228 * t89 - t232 * t65;
t55 = t228 * t65 + t232 * t89;
t377 = -pkin(6) * (t229 * t330 + t233 * t55) + pkin(1) * t53;
t375 = pkin(2) * t53;
t374 = pkin(7) * t53;
t373 = pkin(2) * t330 - pkin(7) * t55;
t159 = t177 * t211;
t266 = -t227 * t170 + t231 * t240;
t244 = t177 * qJD(4) - t266;
t79 = t159 + t244;
t317 = t175 ^ 2;
t155 = t317 - t315;
t95 = t227 * t155 - t282;
t99 = t231 * t155 + t288;
t371 = t229 * (t228 * t95 - t232 * t99) - t233 * t79;
t326 = t174 - t317;
t328 = -t159 + t244;
t46 = -t227 * t328 + t231 * t330;
t304 = t227 * t330;
t48 = t231 * t328 + t304;
t370 = t229 * (t228 * t46 + t232 * t48) + t233 * t326;
t369 = pkin(3) * t65;
t368 = pkin(8) * t65;
t367 = pkin(8) * t89;
t366 = t228 * t48 - t232 * t46;
t364 = t228 * t99 + t232 * t95;
t156 = -t174 + t315;
t121 = t194 + t136;
t287 = t227 * t121;
t355 = t231 * t156 - t287;
t117 = t231 * t121;
t356 = -t227 * t156 - t117;
t362 = t228 * t356 + t232 * t355;
t325 = -t295 + t106;
t361 = t229 * (-t228 * t355 + t232 * t356) - t233 * t325;
t324 = -t315 - t317;
t333 = t231 * t324 + t287;
t334 = t227 * t324 - t117;
t343 = t228 * t333 + t232 * t334;
t360 = pkin(2) * t343;
t359 = pkin(7) * t343;
t344 = -t228 * t334 + t232 * t333;
t354 = -pkin(2) * t328 + pkin(7) * t344;
t353 = pkin(6) * (t229 * t328 + t233 * t344) - pkin(1) * t343;
t111 = -t317 - t174;
t352 = pkin(2) * t111;
t351 = pkin(3) * t111;
t350 = pkin(3) * t334;
t349 = pkin(8) * t333;
t348 = pkin(8) * t334;
t346 = qJ(5) * t330;
t235 = qJD(1) ^ 2;
t230 = sin(qJ(1));
t234 = cos(qJ(1));
t256 = t234 * g(1) + t230 * g(2);
t296 = qJDD(1) * pkin(6);
t193 = -t235 * pkin(1) - t256 + t296;
t257 = -t233 * pkin(2) - t229 * pkin(7);
t263 = t235 * t257 + t193;
t309 = t233 * g(3);
t314 = qJD(2) ^ 2;
t151 = -qJDD(2) * pkin(2) - t314 * pkin(7) + t229 * t263 + t309;
t183 = -t215 * pkin(3) - t200 * pkin(8);
t316 = t198 ^ 2;
t91 = -t240 * pkin(3) - t316 * pkin(8) + t200 * t183 + t151;
t347 = pkin(4) * t244 - t346 + t91;
t345 = t229 * t111;
t293 = t200 * t198;
t242 = -t197 - t293;
t332 = t228 * t242;
t331 = t232 * t242;
t186 = t198 * t215;
t144 = t170 - t186;
t140 = (qJD(3) + t215) * t200 - t249;
t243 = (t175 * t227 + t177 * t231) * t211;
t292 = t211 * t227;
t153 = t177 * t292;
t291 = t211 * t231;
t272 = t175 * t291;
t254 = -t153 + t272;
t321 = t228 * t254 + t232 * t243;
t245 = t227 * t244 - t272;
t255 = -t175 * t292 - t231 * t244;
t320 = t228 * t245 + t232 * t255;
t319 = t229 * (-t228 * t243 + t232 * t254) + t233 * t194;
t271 = t233 * t136;
t318 = t229 * (-t228 * t255 + t232 * t245) + t271;
t196 = t200 ^ 2;
t213 = t215 ^ 2;
t267 = t230 * g(1) - t234 * g(2);
t192 = qJDD(1) * pkin(1) + t235 * pkin(6) + t267;
t252 = -t204 + t217;
t253 = t203 + t269;
t139 = pkin(2) * t252 - pkin(7) * t253 - t192;
t310 = t229 * g(3);
t152 = -t314 * pkin(2) + qJDD(2) * pkin(7) + t233 * t263 - t310;
t108 = -t232 * t139 + t228 * t152;
t62 = t242 * pkin(3) - t144 * pkin(8) - t108;
t109 = t228 * t139 + t232 * t152;
t64 = -t316 * pkin(3) + pkin(8) * t240 + t215 * t183 + t109;
t36 = t227 * t64 - t231 * t62;
t37 = t227 * t62 + t231 * t64;
t18 = t227 * t37 - t231 * t36;
t313 = pkin(3) * t18;
t75 = t231 * t325;
t81 = (-qJD(4) - t211) * t177 + t266;
t47 = t227 * t81 - t75;
t312 = pkin(3) * t47;
t311 = pkin(4) * t231;
t277 = qJD(5) * t211;
t206 = -0.2e1 * t277;
t132 = t175 * pkin(4) - t177 * qJ(5);
t264 = -t194 * qJ(5) - t175 * t132 + t37;
t246 = -pkin(4) * t315 + t264;
t29 = t206 + t246;
t31 = t194 * pkin(4) - qJ(5) * t315 + t177 * t132 + qJDD(5) + t36;
t308 = -pkin(4) * t31 + qJ(5) * t29;
t307 = -pkin(4) * t325 - qJ(5) * t79;
t305 = t227 * t325;
t303 = t227 * t91;
t302 = t228 * t18;
t300 = t231 * t91;
t299 = t232 * t18;
t297 = qJ(5) * t231;
t290 = t215 * t228;
t289 = t215 * t232;
t285 = t228 * t151;
t164 = t197 - t293;
t284 = t228 * t164;
t214 = t233 * t235 * t229;
t283 = t229 * (qJDD(2) + t214);
t281 = t232 * t151;
t280 = t232 * t164;
t279 = t233 * (-t214 + qJDD(2));
t276 = qJD(3) - t215;
t270 = t233 * t293;
t268 = -qJ(5) * t227 - pkin(3);
t19 = t227 * t36 + t231 * t37;
t59 = t228 * t108 + t232 * t109;
t181 = t229 * t193 + t309;
t182 = t233 * t193 - t310;
t265 = t229 * t181 + t233 * t182;
t13 = t227 * t29 - t231 * t31;
t262 = pkin(3) * t13 + t308;
t45 = -t227 * t79 - t75;
t261 = pkin(3) * t45 + t307;
t260 = -t37 - t369;
t73 = t227 * t106 - t177 * t291;
t74 = t231 * t106 + t153;
t258 = t229 * (-t228 * t73 + t232 * t74) - t271;
t251 = t232 * t108 - t228 * t109;
t250 = -pkin(1) + t257;
t247 = -t36 + t350;
t241 = -pkin(4) * t327 - qJ(5) * t329 + t246;
t239 = t241 + t369;
t238 = -pkin(4) * t121 + qJ(5) * t324 - t31;
t237 = t238 + t350;
t236 = 0.2e1 * qJD(5) * t177 - t347;
t224 = t233 ^ 2;
t223 = t229 ^ 2;
t221 = t224 * t235;
t219 = t223 * t235;
t205 = -0.2e1 * t217 + t273;
t202 = t218 + 0.2e1 * t269;
t185 = -t196 + t213;
t184 = -t213 + t316;
t179 = t196 - t316;
t178 = -t196 - t213;
t171 = -t213 - t316;
t163 = t196 + t316;
t145 = t276 * t198 + t248;
t143 = t170 + t186;
t141 = -t276 * t200 + t249;
t131 = -t228 * t178 + t280;
t130 = t232 * t178 + t284;
t126 = t232 * t171 - t332;
t125 = t228 * t171 + t331;
t101 = -t140 * t232 + t228 * t144;
t57 = t300 + t368;
t52 = t303 - t348;
t51 = t231 * t81 + t305;
t49 = -t231 * t79 + t305;
t43 = t228 * t74 + t232 * t73;
t34 = -pkin(3) * t330 + t303 + t367;
t33 = -pkin(3) * t328 - t300 + t349;
t32 = (-pkin(4) * t211 - 0.2e1 * qJD(5)) * t177 + t347;
t27 = (-t328 + t159) * pkin(4) + t236;
t26 = pkin(4) * t159 + t236 + t346;
t25 = -t228 * t47 + t232 * t51;
t24 = -t228 * t45 + t232 * t49;
t23 = t228 * t51 + t232 * t47;
t22 = t228 * t49 + t232 * t45;
t21 = -qJ(5) * t111 + t31;
t20 = t206 + (-t111 - t315) * pkin(4) + t264;
t17 = -t227 * t27 - t297 * t328 - t348;
t16 = -pkin(3) * t91 + pkin(8) * t19;
t15 = -pkin(4) * t304 + t231 * t26 - t368;
t14 = t227 * t31 + t231 * t29;
t12 = t231 * t27 + t268 * t328 + t349;
t11 = -t367 + t227 * t26 + (pkin(3) + t311) * t330;
t10 = -pkin(8) * t47 - t18;
t9 = pkin(8) * t51 + t19 - t351;
t8 = -pkin(8) * t45 - t227 * t20 + t231 * t21;
t7 = t232 * t19 - t302;
t6 = t228 * t19 + t299;
t5 = pkin(8) * t49 + t231 * t20 + t227 * t21 - t351;
t4 = -pkin(8) * t13 + (pkin(4) * t227 - t297) * t32;
t3 = -t228 * t13 + t232 * t14;
t2 = t232 * t13 + t228 * t14;
t1 = pkin(8) * t14 + (t268 - t311) * t32;
t28 = [0, 0, 0, 0, 0, qJDD(1), t267, t256, 0, 0, t253 * t229, t233 * t202 + t229 * t205, t283 + t233 * (-t219 + t314), -t252 * t233, t229 * (t221 - t314) + t279, 0, t233 * t192 + pkin(1) * t205 + pkin(6) * (t233 * (-t221 - t314) - t283), -t229 * t192 - pkin(1) * t202 + pkin(6) * (-t279 - t229 * (-t219 - t314)), pkin(1) * (t219 + t221) + (t223 + t224) * t296 + t265, pkin(1) * t192 + pkin(6) * t265, t229 * (t232 * t170 + t200 * t290) - t270, t229 * (t232 * t141 - t228 * t143) - t233 * t179, t229 * (-t228 * t185 + t331) - t233 * t144, t229 * (-t198 * t289 - t228 * t240) + t270, t229 * (t232 * t184 + t284) + t233 * t140, t233 * t197 + t229 * (t198 * t232 - t200 * t228) * t215, t229 * (-pkin(7) * t125 + t285) + t233 * (-pkin(2) * t125 + t108) - pkin(1) * t125 + pkin(6) * (t233 * t126 - t229 * t141), t229 * (-pkin(7) * t130 + t281) + t233 * (-pkin(2) * t130 + t109) - pkin(1) * t130 + pkin(6) * (t233 * t131 - t229 * t145), t229 * t251 + pkin(6) * (t233 * t101 - t229 * t163) + t250 * (-t140 * t228 - t232 * t144), pkin(6) * (t229 * t151 + t233 * t59) - t250 * t251, t258, -t370, t361, t318, -t371, t319, t229 * (-t228 * t33 + t232 * t52 - t359) + t233 * (-t247 - t360) + t353, t229 * (-t228 * t34 + t232 * t57 - t374) + t233 * (-t260 - t375) - t377, t229 * (-pkin(7) * t23 + t232 * t10 - t228 * t9) + t233 * (-pkin(2) * t23 - t312) - pkin(1) * t23 + pkin(6) * (t233 * t25 + t345), t229 * (-pkin(7) * t6 - pkin(8) * t299 - t228 * t16) + t233 * (-pkin(2) * t6 - t313) - pkin(1) * t6 + pkin(6) * (t229 * t91 + t233 * t7), t258, t361, t370, t319, t371, t318, t229 * (-t228 * t12 + t232 * t17 - t359) + t233 * (-t237 - t360) + t353, t229 * (-pkin(7) * t22 - t228 * t5 + t232 * t8) + t233 * (-pkin(2) * t22 - t261) - pkin(1) * t22 + pkin(6) * (t233 * t24 + t345), t229 * (-t228 * t11 + t232 * t15 + t374) + t233 * (-t239 + 0.2e1 * t277 + t375) + t377, t229 * (-pkin(7) * t2 - t228 * t1 + t232 * t4) + t233 * (-pkin(2) * t2 - t262) - pkin(1) * t2 + pkin(6) * (t229 * t32 + t233 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, t219 - t221, t218, t214, t273, qJDD(2), -t181, -t182, 0, 0, t228 * t170 - t200 * t289, t228 * t141 + t232 * t143, t232 * t185 + t332, -t198 * t290 + t232 * t240, t228 * t184 - t280, (t198 * t228 + t200 * t232) * t215, pkin(2) * t141 + pkin(7) * t126 - t281, pkin(2) * t145 + pkin(7) * t131 + t285, pkin(2) * t163 + pkin(7) * t101 + t59, -pkin(2) * t151 + pkin(7) * t59, t43, -t366, t362, t320, t364, t321, t228 * t52 + t232 * t33 + t354, t228 * t57 + t232 * t34 - t373, pkin(7) * t25 + t228 * t10 + t232 * t9 - t352, -pkin(2) * t91 + pkin(7) * t7 - pkin(8) * t302 + t232 * t16, t43, t362, t366, t321, -t364, t320, t232 * t12 + t228 * t17 + t354, pkin(7) * t24 + t228 * t8 + t232 * t5 - t352, t232 * t11 + t228 * t15 + t373, -pkin(2) * t32 + pkin(7) * t3 + t232 * t1 + t228 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, t179, t144, -t293, -t140, -t197, -t108, -t109, 0, 0, t136, t326, t325, -t136, -t79, -t194, t247, t260, t312, t313, t136, t325, -t326, -t194, t79, -t136, t237, t261, t206 + t239, t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t326, t325, -t136, -t79, -t194, -t36, -t37, 0, 0, t136, t325, -t326, -t194, t79, -t136, t238, t307, t206 + t241, t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t325, t327, t31;];
tauJ_reg = t28;
