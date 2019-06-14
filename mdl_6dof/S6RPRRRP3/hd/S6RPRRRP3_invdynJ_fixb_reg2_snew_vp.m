% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:20:08
% EndTime: 2019-05-06 01:20:29
% DurationCPUTime: 8.66s
% Computational Cost: add. (21012->426), mult. (40939->546), div. (0->0), fcn. (27900->10), ass. (0->275)
t250 = sin(pkin(10));
t251 = cos(pkin(10));
t256 = sin(qJ(3));
t260 = cos(qJ(3));
t255 = sin(qJ(4));
t259 = cos(qJ(4));
t309 = qJD(1) * t256;
t215 = -t259 * qJD(3) + t255 * t309;
t305 = qJD(1) * qJD(3);
t297 = t260 * t305;
t304 = t256 * qJDD(1);
t221 = t297 + t304;
t277 = -qJDD(3) * t255 - t221 * t259;
t187 = -qJD(4) * t215 - t277;
t216 = qJD(3) * t255 + t259 * t309;
t254 = sin(qJ(5));
t258 = cos(qJ(5));
t192 = t258 * t215 + t216 * t254;
t278 = t259 * qJDD(3) - t255 * t221;
t268 = -qJD(4) * t216 + t278;
t264 = -t192 * qJD(5) + t258 * t187 + t254 * t268;
t237 = qJD(1) * t260 - qJD(4);
t231 = -qJD(5) + t237;
t320 = t192 * t231;
t356 = t264 + t320;
t194 = -t215 * t254 + t216 * t258;
t154 = t194 * t192;
t241 = t256 * t305;
t303 = t260 * qJDD(1);
t222 = -t241 + t303;
t214 = -qJDD(4) + t222;
t211 = -qJDD(5) + t214;
t360 = -t154 + t211;
t328 = t360 * t254;
t191 = t194 ^ 2;
t345 = t231 ^ 2;
t358 = -t191 - t345;
t102 = t258 * t358 + t328;
t327 = t360 * t258;
t104 = -t254 * t358 + t327;
t49 = t102 * t255 - t104 * t259;
t40 = -t256 * t356 + t260 * t49;
t65 = t102 * t259 + t104 * t255;
t417 = pkin(1) * (t250 * t40 + t251 * t65) + pkin(2) * t65 + pkin(7) * t40;
t413 = pkin(3) * t65;
t412 = pkin(8) * t65;
t410 = pkin(3) * t356 + pkin(8) * t49;
t409 = t256 * t49 + t260 * t356;
t347 = t192 ^ 2;
t170 = t347 - t345;
t110 = t170 * t254 - t327;
t114 = t170 * t258 + t328;
t174 = t194 * t231;
t293 = -t187 * t254 + t258 * t268;
t273 = qJD(5) * t194 - t293;
t93 = t174 + t273;
t408 = t256 * (t110 * t255 - t114 * t259) - t260 * t93;
t357 = t191 - t347;
t359 = -t174 + t273;
t58 = -t254 * t359 + t258 * t356;
t333 = t254 * t356;
t60 = t258 * t359 + t333;
t407 = t256 * (t255 * t58 + t259 * t60) + t260 * t357;
t406 = pkin(4) * t102;
t405 = pkin(9) * t102;
t404 = pkin(9) * t104;
t403 = t255 * t60 - t259 * t58;
t401 = t110 * t259 + t114 * t255;
t171 = -t191 + t345;
t135 = t211 + t154;
t326 = t135 * t254;
t389 = t258 * t171 - t326;
t129 = t258 * t135;
t390 = -t171 * t254 - t129;
t399 = t255 * t390 + t259 * t389;
t355 = -t320 + t264;
t398 = t256 * (-t255 * t389 + t259 * t390) - t260 * t355;
t354 = -t345 - t347;
t363 = t258 * t354 + t326;
t364 = t254 * t354 - t129;
t375 = t255 * t363 + t259 * t364;
t376 = -t255 * t364 + t259 * t363;
t386 = t256 * t359 + t260 * t376;
t397 = pkin(1) * (t250 * t386 - t251 * t375) + pkin(7) * t386 - pkin(2) * t375;
t395 = pkin(3) * t375;
t394 = pkin(8) * t375;
t388 = -pkin(3) * t359 + pkin(8) * t376;
t387 = t256 * t376 - t260 * t359;
t123 = -t347 - t191;
t385 = pkin(3) * t123;
t384 = pkin(4) * t123;
t383 = pkin(4) * t364;
t382 = pkin(9) * t363;
t381 = pkin(9) * t364;
t379 = qJ(6) * t356;
t311 = -g(3) + qJDD(2);
t240 = t260 * t311;
t262 = qJD(1) ^ 2;
t257 = sin(qJ(1));
t261 = cos(qJ(1));
t295 = t257 * g(1) - g(2) * t261;
t217 = qJDD(1) * pkin(1) + t295;
t284 = g(1) * t261 + g(2) * t257;
t218 = -pkin(1) * t262 - t284;
t310 = t250 * t217 + t251 * t218;
t184 = -pkin(2) * t262 + qJDD(1) * pkin(7) + t310;
t285 = -t260 * pkin(3) - t256 * pkin(8);
t290 = t262 * t285 + t184;
t344 = qJD(3) ^ 2;
t148 = -qJDD(3) * pkin(3) - t344 * pkin(8) + t290 * t256 - t240;
t199 = -pkin(4) * t237 - pkin(9) * t216;
t346 = t215 ^ 2;
t86 = -t268 * pkin(4) - t346 * pkin(9) + t199 * t216 + t148;
t380 = pkin(5) * t273 - t379 + t86;
t378 = t123 * t256;
t377 = t123 * t260;
t318 = t216 * t215;
t270 = -t214 - t318;
t362 = t255 * t270;
t361 = t259 * t270;
t202 = t215 * t237;
t159 = t187 - t202;
t155 = (qJD(4) + t237) * t216 - t278;
t272 = (t192 * t254 + t194 * t258) * t231;
t317 = t231 * t254;
t168 = t194 * t317;
t316 = t231 * t258;
t302 = t192 * t316;
t282 = -t168 + t302;
t351 = t255 * t282 + t259 * t272;
t274 = t254 * t273 - t302;
t283 = -t192 * t317 - t258 * t273;
t350 = t255 * t274 + t259 * t283;
t349 = t256 * (-t255 * t272 + t259 * t282) + t260 * t211;
t301 = t260 * t154;
t348 = t256 * (-t255 * t283 + t259 * t274) + t301;
t213 = t216 ^ 2;
t235 = t237 ^ 2;
t292 = t251 * t217 - t250 * t218;
t183 = -qJDD(1) * pkin(2) - t262 * pkin(7) - t292;
t280 = -t222 + t241;
t281 = t221 + t297;
t141 = t280 * pkin(3) - t281 * pkin(8) + t183;
t294 = t256 * t311;
t149 = -t344 * pkin(3) + qJDD(3) * pkin(8) + t290 * t260 + t294;
t88 = -t259 * t141 + t149 * t255;
t71 = t270 * pkin(4) - t159 * pkin(9) - t88;
t89 = t255 * t141 + t259 * t149;
t73 = -t346 * pkin(4) + pkin(9) * t268 + t237 * t199 + t89;
t38 = t254 * t73 - t258 * t71;
t39 = t254 * t71 + t258 * t73;
t18 = t254 * t39 - t258 * t38;
t343 = pkin(4) * t18;
t87 = t258 * t355;
t95 = (-qJD(5) - t231) * t194 + t293;
t59 = t254 * t95 - t87;
t342 = pkin(4) * t59;
t341 = pkin(5) * t258;
t308 = qJD(6) * t231;
t224 = -0.2e1 * t308;
t150 = pkin(5) * t192 - qJ(6) * t194;
t291 = -t211 * qJ(6) - t192 * t150 + t39;
t275 = -pkin(5) * t345 + t291;
t29 = t224 + t275;
t35 = t211 * pkin(5) - qJ(6) * t345 + t150 * t194 + qJDD(6) + t38;
t340 = -pkin(5) * t35 + qJ(6) * t29;
t339 = -pkin(5) * t355 - qJ(6) * t93;
t338 = t18 * t255;
t337 = t18 * t259;
t336 = t254 * t86;
t334 = t254 * t355;
t332 = t258 * t86;
t329 = qJ(6) * t258;
t325 = t148 * t255;
t324 = t148 * t259;
t179 = t214 - t318;
t322 = t179 * t255;
t321 = t179 * t259;
t315 = t237 * t255;
t314 = t237 * t259;
t236 = t260 * t262 * t256;
t229 = qJDD(3) + t236;
t313 = t256 * t229;
t228 = -t236 + qJDD(3);
t312 = t260 * t228;
t307 = qJD(4) - t237;
t300 = t260 * t318;
t298 = pkin(1) * t250 + pkin(7);
t296 = -qJ(6) * t254 - pkin(4);
t19 = t254 * t38 + t258 * t39;
t55 = t255 * t88 + t259 * t89;
t165 = t184 * t256 - t240;
t166 = t260 * t184 + t294;
t121 = t256 * t165 + t260 * t166;
t13 = t254 * t29 - t258 * t35;
t289 = pkin(4) * t13 + t340;
t57 = -t254 * t93 - t87;
t288 = pkin(4) * t57 + t339;
t287 = -t39 + t406;
t83 = -t194 * t316 + t254 * t264;
t84 = t258 * t264 + t168;
t286 = t256 * (-t255 * t83 + t259 * t84) - t301;
t279 = t255 * t89 - t259 * t88;
t276 = -t38 + t383;
t271 = -pkin(1) * t251 - pkin(2) + t285;
t269 = -pkin(5) * t358 - qJ(6) * t360 + t275;
t267 = t269 - t406;
t266 = -pkin(5) * t135 + qJ(6) * t354 - t35;
t265 = t266 + t383;
t263 = 0.2e1 * qJD(6) * t194 - t380;
t247 = t260 ^ 2;
t246 = t256 ^ 2;
t244 = t247 * t262;
t242 = t246 * t262;
t234 = -t244 - t344;
t233 = -t242 - t344;
t226 = t242 + t244;
t225 = (t246 + t247) * qJDD(1);
t223 = -0.2e1 * t241 + t303;
t220 = 0.2e1 * t297 + t304;
t201 = -t213 + t235;
t200 = -t235 + t346;
t198 = -t233 * t256 - t312;
t197 = t234 * t260 - t313;
t196 = t213 - t346;
t195 = -t213 - t235;
t188 = -t235 - t346;
t178 = t213 + t346;
t160 = t307 * t215 + t277;
t158 = t187 + t202;
t156 = -t307 * t216 + t278;
t146 = -t195 * t255 + t321;
t145 = t195 * t259 + t322;
t140 = t188 * t259 - t362;
t139 = t188 * t255 + t361;
t116 = -t155 * t259 + t159 * t255;
t106 = t146 * t260 - t160 * t256;
t101 = t140 * t260 - t156 * t256;
t64 = t332 - t405;
t63 = t258 * t95 + t334;
t61 = -t258 * t93 + t334;
t53 = t336 - t381;
t52 = t255 * t84 + t259 * t83;
t43 = -pkin(4) * t356 + t336 + t404;
t42 = -pkin(4) * t359 - t332 + t382;
t36 = (-pkin(5) * t231 - 0.2e1 * qJD(6)) * t194 + t380;
t33 = -t255 * t59 + t259 * t63;
t32 = -t255 * t57 + t259 * t61;
t31 = t255 * t63 + t259 * t59;
t30 = t255 * t61 + t259 * t57;
t27 = (-t359 + t174) * pkin(5) + t263;
t26 = pkin(5) * t174 + t263 + t379;
t25 = -qJ(6) * t123 + t35;
t24 = t224 + (-t123 - t345) * pkin(5) + t291;
t23 = t260 * t33 + t378;
t22 = t260 * t32 + t378;
t21 = -t254 * t27 - t329 * t359 - t381;
t20 = -pkin(5) * t333 + t258 * t26 + t405;
t17 = t258 * t27 + t296 * t359 + t382;
t16 = -t404 + t254 * t26 + (pkin(4) + t341) * t356;
t15 = -pkin(4) * t86 + pkin(9) * t19;
t14 = t254 * t35 + t258 * t29;
t12 = -pkin(9) * t59 - t18;
t11 = pkin(9) * t63 + t19 - t384;
t10 = -pkin(9) * t57 - t24 * t254 + t25 * t258;
t9 = pkin(9) * t61 + t24 * t258 + t25 * t254 - t384;
t8 = t19 * t259 - t338;
t7 = t19 * t255 + t337;
t6 = t256 * t86 + t260 * t8;
t5 = -pkin(9) * t13 + (pkin(5) * t254 - t329) * t36;
t4 = -t13 * t255 + t14 * t259;
t3 = t13 * t259 + t14 * t255;
t2 = pkin(9) * t14 + (t296 - t341) * t36;
t1 = t256 * t36 + t260 * t4;
t28 = [0, 0, 0, 0, 0, qJDD(1), t295, t284, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t251 - t250 * t262) + t292, pkin(1) * (-qJDD(1) * t250 - t251 * t262) - t310, 0, pkin(1) * (t250 * t310 + t251 * t292), t281 * t256, t220 * t260 + t223 * t256, t313 + t260 * (-t242 + t344), -t280 * t260, t256 * (t244 - t344) + t312, 0, -t260 * t183 + pkin(2) * t223 + pkin(7) * t197 + pkin(1) * (t197 * t250 + t223 * t251), t256 * t183 - pkin(2) * t220 + pkin(7) * t198 + pkin(1) * (t198 * t250 - t220 * t251), pkin(2) * t226 + pkin(7) * t225 + pkin(1) * (t225 * t250 + t226 * t251) + t121, -pkin(2) * t183 + pkin(7) * t121 + pkin(1) * (t121 * t250 - t183 * t251), t256 * (t187 * t259 + t216 * t315) - t300, t256 * (t156 * t259 - t158 * t255) - t260 * t196, t256 * (-t201 * t255 + t361) - t260 * t159, t256 * (-t215 * t314 - t255 * t268) + t300, t256 * (t200 * t259 + t322) + t260 * t155, t260 * t214 + t256 * (t215 * t259 - t216 * t255) * t237, t256 * (-pkin(8) * t139 + t325) + t260 * (-pkin(3) * t139 + t88) - pkin(2) * t139 + pkin(7) * t101 + pkin(1) * (t101 * t250 - t251 * t139), t256 * (-pkin(8) * t145 + t324) + t260 * (-pkin(3) * t145 + t89) - pkin(2) * t145 + pkin(7) * t106 + pkin(1) * (t250 * t106 - t251 * t145), -t256 * t279 + t298 * (t116 * t260 - t178 * t256) + t271 * (-t155 * t255 - t159 * t259), t298 * (t148 * t256 + t260 * t55) + t271 * t279, t286, -t407, t398, t348, -t408, t349, t256 * (-t255 * t42 + t259 * t53 - t394) + t260 * (-t276 - t395) + t397, t256 * (-t255 * t43 + t259 * t64 - t412) + t260 * (-t287 - t413) - t417, t256 * (-pkin(8) * t31 - t11 * t255 + t12 * t259) + t260 * (-pkin(3) * t31 - t342) - pkin(2) * t31 + pkin(7) * t23 + pkin(1) * (t23 * t250 - t251 * t31), t256 * (-pkin(8) * t7 - pkin(9) * t337 - t15 * t255) + t260 * (-pkin(3) * t7 - t343) - pkin(2) * t7 + pkin(7) * t6 + pkin(1) * (t250 * t6 - t251 * t7), t286, t398, t407, t349, t408, t348, t256 * (-t17 * t255 + t21 * t259 - t394) + t260 * (-t265 - t395) + t397, t256 * (-pkin(8) * t30 + t10 * t259 - t255 * t9) + t260 * (-pkin(3) * t30 - t288) - pkin(2) * t30 + pkin(7) * t22 + pkin(1) * (t22 * t250 - t251 * t30), t256 * (-t16 * t255 + t20 * t259 + t412) + t260 * (-t267 + 0.2e1 * t308 + t413) + t417, t256 * (-pkin(8) * t3 - t2 * t255 + t259 * t5) + t260 * (-pkin(3) * t3 - t289) - pkin(2) * t3 + pkin(7) * t1 + pkin(1) * (t1 * t250 - t251 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, 0, 0, 0, 0, 0, 0, t229 * t260 + t234 * t256, -t228 * t256 + t233 * t260, 0, -t165 * t260 + t166 * t256, 0, 0, 0, 0, 0, 0, t140 * t256 + t156 * t260, t146 * t256 + t160 * t260, t116 * t256 + t178 * t260, -t148 * t260 + t256 * t55, 0, 0, 0, 0, 0, 0, t387, -t409, t256 * t33 - t377, t256 * t8 - t260 * t86, 0, 0, 0, 0, 0, 0, t387, t256 * t32 - t377, t409, t256 * t4 - t260 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, t242 - t244, t304, t236, t303, qJDD(3), -t165, -t166, 0, 0, t187 * t255 - t216 * t314, t156 * t255 + t158 * t259, t201 * t259 + t362, -t215 * t315 + t259 * t268, t200 * t255 - t321, (t215 * t255 + t216 * t259) * t237, pkin(3) * t156 + pkin(8) * t140 - t324, pkin(3) * t160 + pkin(8) * t146 + t325, pkin(3) * t178 + pkin(8) * t116 + t55, -pkin(3) * t148 + pkin(8) * t55, t52, -t403, t399, t350, t401, t351, t255 * t53 + t259 * t42 + t388, t255 * t64 + t259 * t43 - t410, pkin(8) * t33 + t11 * t259 + t12 * t255 - t385, -pkin(3) * t86 + pkin(8) * t8 - pkin(9) * t338 + t15 * t259, t52, t399, t403, t351, -t401, t350, t17 * t259 + t21 * t255 + t388, pkin(8) * t32 + t10 * t255 + t259 * t9 - t385, t16 * t259 + t20 * t255 + t410, -pkin(3) * t36 + pkin(8) * t4 + t2 * t259 + t255 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t196, t159, -t318, -t155, -t214, -t88, -t89, 0, 0, t154, t357, t355, -t154, -t93, -t211, t276, t287, t342, t343, t154, t355, -t357, -t211, t93, -t154, t265, t288, t224 + t267, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t357, t355, -t154, -t93, -t211, -t38, -t39, 0, 0, t154, t355, -t357, -t211, t93, -t154, t266, t339, t224 + t269, t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t355, t358, t35;];
tauJ_reg  = t28;
