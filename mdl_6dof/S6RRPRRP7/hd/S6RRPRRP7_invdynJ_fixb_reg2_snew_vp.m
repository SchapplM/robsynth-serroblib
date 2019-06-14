% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 18:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:16:19
% EndTime: 2019-05-06 18:16:40
% DurationCPUTime: 7.54s
% Computational Cost: add. (19437->441), mult. (40016->517), div. (0->0), fcn. (25937->8), ass. (0->278)
t232 = sin(qJ(4));
t233 = sin(qJ(2));
t236 = cos(qJ(4));
t237 = cos(qJ(2));
t189 = (t233 * t232 + t237 * t236) * qJD(1);
t186 = qJD(5) + t189;
t358 = t186 ^ 2;
t306 = t233 * qJD(1);
t191 = -t232 * t237 * qJD(1) + t236 * t306;
t231 = sin(qJ(5));
t235 = cos(qJ(5));
t303 = qJD(2) - qJD(4);
t167 = t191 * t231 + t235 * t303;
t359 = t167 ^ 2;
t147 = t359 - t358;
t305 = qJD(1) * qJD(2);
t295 = t237 * t305;
t304 = t233 * qJDD(1);
t200 = t295 + t304;
t223 = t237 * qJDD(1);
t296 = t233 * t305;
t201 = t223 - t296;
t289 = t232 * t200 + t236 * t201;
t142 = -t191 * qJD(4) - t289;
t140 = qJDD(5) - t142;
t169 = t235 * t191 - t231 * t303;
t323 = t169 * t167;
t100 = -t323 - t140;
t344 = t100 * t231;
t76 = -t147 * t235 - t344;
t151 = t186 * t169;
t143 = -t189 * qJD(4) + t200 * t236 - t201 * t232;
t302 = qJDD(2) - qJDD(4);
t291 = -t143 * t231 - t235 * t302;
t260 = qJD(5) * t169 - t291;
t88 = -t151 + t260;
t424 = t233 * (t232 * t88 + t236 * t76) - t237 * (t232 * t76 - t236 * t88);
t165 = t169 ^ 2;
t374 = -t165 - t358;
t64 = t235 * t374 + t344;
t423 = pkin(4) * t64;
t422 = pkin(9) * t64;
t338 = t100 * t235;
t66 = -t231 * t374 + t338;
t421 = pkin(9) * t66;
t420 = t232 * t66;
t419 = t236 * t66;
t335 = qJ(3) * t233;
t355 = pkin(2) + pkin(3);
t255 = t237 * t355 + pkin(1) + t335;
t418 = t255 * t64;
t373 = t165 - t359;
t254 = -t235 * t143 + t231 * t302;
t251 = -t167 * qJD(5) - t254;
t324 = t167 * t186;
t380 = t324 - t251;
t345 = t231 * t380;
t375 = t151 + t260;
t50 = t235 * t375 - t345;
t415 = t233 * (-t232 * t373 + t236 * t50) - t237 * (t232 * t50 + t236 * t373);
t369 = -t323 + t140;
t343 = t231 * t369;
t366 = -t358 - t359;
t378 = t235 * t366 - t343;
t394 = t232 * t375 + t236 * t378;
t412 = pkin(8) * t394;
t395 = t232 * t378 - t236 * t375;
t411 = pkin(8) * t395;
t320 = t191 * t189;
t398 = -t302 - t320;
t410 = t232 * t398;
t409 = t236 * t398;
t408 = t380 * qJ(6);
t72 = -t147 * t231 + t338;
t407 = qJ(3) * t394 - t355 * t395;
t337 = t235 * t369;
t379 = t231 * t366 + t337;
t406 = t255 * t379 + pkin(7) * (t233 * t395 + t237 * t394);
t367 = t324 + t251;
t148 = -t165 + t358;
t396 = -t148 * t231 + t337;
t405 = t233 * (t232 * t367 + t236 * t396) + t237 * (-t232 * t396 + t236 * t367);
t404 = pkin(4) * t379;
t403 = pkin(9) * t378;
t402 = pkin(9) * t379;
t397 = t148 * t235 + t343;
t372 = t165 + t359;
t393 = pkin(4) * t372;
t238 = qJD(2) ^ 2;
t228 = t233 ^ 2;
t239 = qJD(1) ^ 2;
t317 = t228 * t239;
t209 = t238 + t317;
t213 = t237 * t239 * t233;
t207 = qJDD(2) - t213;
t313 = t237 * t207;
t392 = pkin(7) * (-t209 * t233 + t313);
t389 = t232 * t372;
t384 = t236 * t372;
t286 = t189 * t303;
t124 = t143 + t286;
t381 = t201 - t296;
t48 = -t231 * t375 - t235 * t380;
t301 = t303 ^ 2;
t377 = pkin(7) - pkin(8);
t376 = (t167 * t231 + t169 * t235) * t186;
t229 = t237 ^ 2;
t316 = t229 * t239;
t371 = t313 + (-t238 + t316) * t233;
t131 = pkin(5) * t167 - qJ(6) * t169;
t208 = -qJD(2) * pkin(3) - pkin(8) * t306;
t234 = sin(qJ(1));
t354 = cos(qJ(1));
t292 = t234 * g(1) - t354 * g(2);
t192 = qJDD(1) * pkin(1) + t239 * pkin(7) + t292;
t248 = t381 * pkin(2) + t192;
t308 = qJD(2) * t237;
t357 = 2 * qJD(3);
t101 = t200 * qJ(3) + t201 * pkin(3) - pkin(8) * t316 + (qJ(3) * t308 + (t357 + t208) * t233) * qJD(1) + t248;
t285 = t303 * t191;
t56 = -t124 * pkin(9) + (-t142 - t285) * pkin(4) + t101;
t158 = pkin(4) * t189 - pkin(9) * t191;
t263 = g(1) * t354 + t234 * g(2);
t333 = qJDD(1) * pkin(7);
t250 = -t239 * pkin(1) - t263 + t333;
t183 = t237 * t250;
t351 = t237 * pkin(2);
t271 = -t335 - t351;
t258 = t239 * t271;
t284 = qJDD(2) * qJ(3) + qJD(2) * t357 + t237 * t258 + t183;
t352 = t233 * g(3);
t259 = t284 - t352;
t349 = t238 * pkin(2);
t139 = t259 - t349;
t115 = -pkin(3) * t316 - pkin(8) * t201 + qJD(2) * t208 + t139;
t206 = qJDD(2) + t213;
t247 = t233 * t250;
t242 = -qJDD(2) * pkin(2) - t238 * qJ(3) + qJDD(3) + t247;
t350 = t237 * g(3);
t240 = t350 - t200 * pkin(8) - t206 * pkin(3) + (pkin(8) * t308 + t271 * t306) * qJD(1) + t242;
t69 = t236 * t115 + t232 * t240;
t59 = -pkin(4) * t301 - pkin(9) * t302 - t189 * t158 + t69;
t29 = t231 * t56 + t235 * t59;
t288 = t140 * qJ(6) - t167 * t131 + t29;
t365 = -pkin(5) * (t374 + t358) - qJ(6) * t100 + t288;
t68 = t115 * t232 - t236 * t240;
t58 = t302 * pkin(4) - t301 * pkin(9) + t158 * t191 + t68;
t249 = t260 * pkin(5) + t408 + t58;
t26 = (pkin(5) * t186 - (2 * qJD(6))) * t169 + t249;
t293 = qJ(6) * t231 + pkin(4);
t353 = pkin(5) * t235;
t307 = qJD(6) * t186;
t179 = 0.2e1 * t307;
t277 = t179 + t288;
t23 = -pkin(5) * t358 + t277;
t28 = t231 * t59 - t235 * t56;
t24 = -t140 * pkin(5) - qJ(6) * t358 + t131 * t169 + qJDD(6) + t28;
t7 = t23 * t235 + t231 * t24;
t364 = -t26 * (t293 + t353) + pkin(9) * t7;
t246 = 0.2e1 * qJD(6) * t169 - t249;
t20 = (-t375 - t151) * pkin(5) + t246;
t363 = t20 * t235 - t293 * t375 + t403;
t19 = -pkin(5) * t151 + t246 - t408;
t362 = -t421 + t19 * t231 - (pkin(4) + t353) * t380;
t321 = t186 * t235;
t298 = t167 * t321;
t262 = t231 * t260 + t298;
t299 = t236 * t323;
t300 = t232 * t323;
t361 = t233 * (t236 * t262 - t300) + t237 * (-t232 * t262 - t299);
t322 = t186 * t231;
t145 = t169 * t322;
t278 = t145 - t298;
t360 = t233 * (t140 * t232 + t236 * t278) + t237 * (t140 * t236 - t232 * t278);
t187 = t189 ^ 2;
t188 = t191 ^ 2;
t348 = t231 * t58;
t346 = t231 * t367;
t341 = t235 * t58;
t339 = t235 * t367;
t334 = qJ(6) * t235;
t332 = t101 * t232;
t331 = t101 * t236;
t155 = -t320 + t302;
t326 = t155 * t232;
t325 = t155 * t236;
t202 = t223 - 0.2e1 * t296;
t315 = t233 * t202;
t314 = t233 * t206;
t211 = -t238 - t316;
t311 = pkin(7) * (t211 * t237 - t314) + pkin(1) * t202;
t309 = t228 + t229;
t204 = t309 * t239;
t310 = pkin(1) * t204 + t309 * t333;
t294 = pkin(4) * t232 + qJ(3);
t16 = t231 * t28 + t235 * t29;
t175 = t247 + t350;
t176 = t183 - t352;
t290 = t175 * t233 + t237 * t176;
t287 = pkin(4) * t236 + t355;
t283 = -pkin(4) * t58 + pkin(9) * t16;
t282 = -pkin(5) * t24 + qJ(6) * t23;
t281 = -pkin(5) * t367 - qJ(6) * t88;
t275 = t232 * t286;
t274 = t232 * t285;
t273 = t236 * t286;
t272 = t236 * t285;
t270 = t200 + t295;
t269 = t231 * t29 - t235 * t28;
t36 = t232 * t69 - t236 * t68;
t37 = t232 * t68 + t236 * t69;
t199 = 0.2e1 * t295 + t304;
t266 = t237 * t199 + t315;
t82 = t167 * t322 - t235 * t260;
t121 = t191 * qJD(2) + t289;
t95 = (qJD(5) + t186) * t167 + t254;
t257 = pkin(4) * t95 + t348 + t421;
t256 = pkin(4) * t375 + t341 - t403;
t90 = (-qJD(5) + t186) * t169 + t291;
t53 = t235 * t90 + t346;
t253 = pkin(9) * t53 + t16 + t393;
t17 = (t372 - t358) * pkin(5) + t277;
t18 = qJ(6) * t372 + t24;
t51 = -t235 * t88 + t346;
t252 = pkin(9) * t51 + t17 * t235 + t18 * t231 + t393;
t245 = pkin(5) * t369 + qJ(6) * t366 - t24;
t86 = t235 * t251 - t145;
t244 = t233 * (t236 * t86 + t300) + t237 * (-t232 * t86 + t299);
t243 = t306 * t357 + t248;
t241 = t233 * t258 + t242;
t141 = t241 + t350;
t205 = (t228 - t229) * t239;
t178 = -t188 + t301;
t177 = t187 - t301;
t174 = -t188 - t301;
t173 = t314 + t237 * (t238 - t317);
t172 = t270 * t233;
t171 = t381 * t237;
t160 = t188 - t187;
t154 = -t301 - t187;
t127 = -t174 * t232 + t325;
t126 = t174 * t236 + t326;
t125 = t143 - t286;
t120 = (0.2e1 * qJD(4) - qJD(2)) * t191 + t289;
t117 = t154 * t236 - t410;
t116 = t154 * t232 + t409;
t85 = t169 * t321 + t231 * t251;
t79 = -t121 * t236 + t125 * t232;
t78 = -t121 * t232 - t125 * t236;
t49 = t231 * t90 - t339;
t47 = -t231 * t88 - t339;
t44 = -t232 * t95 + t419;
t42 = t236 * t95 + t420;
t40 = t232 * t380 - t419;
t38 = -t236 * t380 - t420;
t35 = t236 * t53 - t389;
t34 = t236 * t51 - t389;
t33 = t232 * t53 + t384;
t32 = t232 * t51 + t384;
t31 = t341 - t422;
t30 = t348 - t402;
t25 = -pkin(4) * t47 - t281;
t22 = t29 - t423;
t21 = t28 - t404;
t14 = -t245 - t404;
t13 = -t20 * t231 - t334 * t375 - t402;
t12 = pkin(5) * t345 + t19 * t235 + t422;
t11 = -0.2e1 * t307 - t365 + t423;
t10 = t16 * t236 + t232 * t58;
t9 = t16 * t232 - t236 * t58;
t8 = -pkin(9) * t49 - t269;
t6 = t23 * t231 - t235 * t24;
t5 = -pkin(9) * t47 - t17 * t231 + t18 * t235;
t4 = t232 * t26 + t236 * t7;
t3 = t232 * t7 - t236 * t26;
t2 = -pkin(9) * t6 + (pkin(5) * t231 - t334) * t26;
t1 = -pkin(4) * t6 - t282;
t15 = [0, 0, 0, 0, 0, qJDD(1), t292, t263, 0, 0, t172, t266, t173, t171, t371, 0, t192 * t237 + t311, -pkin(1) * t199 - t233 * t192 - t392, t290 + t310, pkin(1) * t192 + pkin(7) * t290, t172, t173, -t266, 0, -t371, t171, t237 * (pkin(2) * t202 + t243) + (t237 * t270 + t315) * qJ(3) + t311, t233 * (qJ(3) * t204 + t241) + t237 * (pkin(2) * t204 + t284 - t349) + t310, t233 * t243 + t392 + (pkin(1) + t351) * t199 + (t199 + t270) * t335, pkin(7) * (t139 * t237 + t141 * t233) + (pkin(1) - t271) * (qJ(3) * t270 + t243), t233 * (t236 * t143 + t274) + t237 * (-t232 * t143 + t272), t233 * (-t120 * t236 - t124 * t232) + t237 * (t120 * t232 - t124 * t236), t233 * (-t178 * t232 + t409) + t237 * (-t178 * t236 - t410), t233 * (-t232 * t142 - t273) + t237 * (-t236 * t142 + t275), t233 * (t177 * t236 + t326) + t237 * (-t177 * t232 + t325), t233 * (t273 - t274) + t237 * (-t275 - t272), t233 * (-pkin(8) * t116 + t332) + t237 * (-pkin(8) * t117 + t331) + pkin(7) * (t116 * t233 + t117 * t237) + t255 * t120, t233 * (-pkin(8) * t126 + t331) + t237 * (-pkin(8) * t127 - t332) + pkin(7) * (t126 * t233 + t127 * t237) + t255 * t124, t233 * (-pkin(8) * t78 - t36) + t237 * (-pkin(8) * t79 - t37) + pkin(7) * (t233 * t78 + t237 * t79) + t255 * (-t187 - t188), t101 * t255 + t377 * (t233 * t36 + t237 * t37), t244, -t415, t405, t361, -t424, t360, t233 * (-t21 * t232 + t236 * t30 - t411) + t237 * (-t21 * t236 - t232 * t30 - t412) + t406, t233 * (-pkin(8) * t42 - t22 * t232 + t236 * t31) + t237 * (-pkin(8) * t44 - t22 * t236 - t232 * t31) + pkin(7) * (t233 * t42 + t237 * t44) + t418, t233 * (-pkin(8) * t33 + t236 * t8) + t237 * (-pkin(8) * t35 - t232 * t8) + pkin(7) * (t233 * t33 + t237 * t35) + (t233 * t294 + t237 * t287 + pkin(1)) * t49, (t233 * (-pkin(9) * t236 + t294) + t237 * (pkin(9) * t232 + t287) + pkin(1)) * t269 + t377 * (t10 * t237 + t233 * t9), t244, t405, t415, t360, t424, t361, t233 * (t13 * t236 - t14 * t232 - t411) + t237 * (-t13 * t232 - t14 * t236 - t412) + t406, t233 * (-pkin(8) * t32 - t232 * t25 + t236 * t5) + t237 * (-pkin(8) * t34 - t232 * t5 - t236 * t25) + pkin(7) * (t233 * t32 + t237 * t34) + t255 * t47, t233 * (-pkin(8) * t38 - t11 * t232 + t12 * t236) + t237 * (-pkin(8) * t40 - t11 * t236 - t12 * t232) + pkin(7) * (t233 * t38 + t237 * t40) - t418, t233 * (-pkin(8) * t3 - t1 * t232 + t2 * t236) + t237 * (-pkin(8) * t4 - t1 * t236 - t2 * t232) + pkin(7) * (t233 * t3 + t237 * t4) + t255 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, t205, t304, t213, t223, qJDD(2), -t175, -t176, 0, 0, -t213, t304, -t205, qJDD(2), -t223, t213, pkin(2) * t206 + qJ(3) * t211 - t141, (-t233 * pkin(2) + qJ(3) * t237) * qJDD(1), qJ(3) * t207 + (t209 - t238) * pkin(2) + t259, -pkin(2) * t141 + qJ(3) * t139, -t320, -t160, -t125, t320, t121, t302, qJ(3) * t117 - t116 * t355 + t68, qJ(3) * t127 - t126 * t355 + t69, qJ(3) * t79 - t355 * t78, qJ(3) * t37 - t355 * t36, -t85, -t48, -t397, -t82, t72, t376, t256 + t407, qJ(3) * t44 - t355 * t42 - t257, qJ(3) * t35 - t33 * t355 - t253, qJ(3) * t10 - t355 * t9 - t283, -t85, -t397, t48, t376, -t72, -t82, -t363 + t407, qJ(3) * t34 - t32 * t355 - t252, qJ(3) * t40 - t355 * t38 - t362, qJ(3) * t4 - t3 * t355 - t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, t304, -t209, t141, 0, 0, 0, 0, 0, 0, t116, t126, t78, t36, 0, 0, 0, 0, 0, 0, t395, t42, t33, t9, 0, 0, 0, 0, 0, 0, t395, t32, t38, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, t160, t125, -t320, -t121, -t302, -t68, -t69, 0, 0, t85, t48, t397, t82, -t72, -t376, -t256, t257, t253, t283, t85, t397, -t48, -t376, t72, t82, t363, t252, t362, t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323, t373, t367, -t323, -t88, t140, -t28, -t29, 0, 0, t323, t367, -t373, t140, t88, -t323, t245, t281, t179 + t365, t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t369, t367, t374, t24;];
tauJ_reg  = t15;
