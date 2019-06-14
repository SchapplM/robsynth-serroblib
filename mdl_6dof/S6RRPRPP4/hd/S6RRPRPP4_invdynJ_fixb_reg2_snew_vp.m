% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRPP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:38:49
% EndTime: 2019-05-06 12:39:17
% DurationCPUTime: 8.23s
% Computational Cost: add. (21226->416), mult. (45128->502), div. (0->0), fcn. (28347->8), ass. (0->270)
t259 = sin(qJ(4));
t260 = sin(qJ(2));
t262 = cos(qJ(4));
t263 = cos(qJ(2));
t327 = t260 * qJ(3);
t362 = -pkin(2) - pkin(8);
t276 = t263 * t362 - pkin(1) - t327;
t317 = qJD(1) * t263;
t220 = qJD(2) * t262 - t259 * t317;
t313 = qJD(1) * qJD(2);
t303 = t260 * t313;
t310 = t263 * qJDD(1);
t225 = -t303 + t310;
t294 = t259 * qJDD(2) + t262 * t225;
t180 = -qJD(4) * t220 - t294;
t218 = qJD(2) * t259 + t262 * t317;
t282 = t262 * qJDD(2) - t259 * t225;
t181 = -qJD(4) * t218 + t282;
t256 = sin(pkin(9));
t257 = cos(pkin(9));
t137 = t180 * t256 + t181 * t257;
t186 = t257 * t218 + t220 * t256;
t318 = qJD(1) * t260;
t242 = qJD(4) + t318;
t339 = t186 * t242;
t382 = t137 - t339;
t245 = t263 * t313;
t247 = t260 * qJDD(1);
t224 = t247 + t245;
t214 = qJDD(4) + t224;
t188 = -t218 * t256 + t220 * t257;
t338 = t188 * t186;
t131 = -t338 - t214;
t345 = t131 * t256;
t185 = t188 ^ 2;
t364 = t242 ^ 2;
t379 = -t185 - t364;
t80 = t257 * t379 + t345;
t344 = t131 * t257;
t82 = -t256 * t379 + t344;
t50 = t259 * t82 + t262 * t80;
t444 = -pkin(7) * (t260 * t50 + t263 * t382) + t276 * (t259 * t80 - t262 * t82);
t443 = pkin(3) * t50;
t441 = -qJ(3) * t382 - t362 * t50;
t169 = t242 * t188;
t297 = -t257 * t180 + t181 * t256;
t104 = -t297 + t169;
t365 = t186 ^ 2;
t163 = t365 - t364;
t87 = -t163 * t256 + t344;
t91 = -t163 * t257 - t345;
t439 = t260 * t104 + t263 * (t259 * t91 + t262 * t87);
t438 = pkin(4) * t80;
t437 = qJ(5) * t80;
t436 = qJ(5) * t82;
t144 = t185 - t365;
t380 = t169 + t297;
t63 = -t380 * t256 + t257 * t382;
t348 = t382 * t256;
t65 = t380 * t257 + t348;
t435 = t260 * t144 + t263 * (t259 * t65 - t262 * t63);
t432 = t259 * t87 - t262 * t91;
t381 = t137 + t339;
t401 = t104 * t257 + t381 * t256;
t402 = t104 * t256 - t257 * t381;
t420 = t259 * t401 + t262 * t402;
t431 = pkin(3) * t420;
t429 = t259 * t63 + t262 * t65;
t113 = -t365 - t185;
t428 = qJ(3) * t113 + t362 * t420;
t427 = t276 * (-t259 * t402 + t262 * t401) + pkin(7) * (t113 * t263 + t260 * t420);
t375 = -t338 + t214;
t343 = t375 * t256;
t374 = -t364 - t365;
t384 = t257 * t374 - t343;
t120 = t257 * t375;
t385 = t256 * t374 + t120;
t400 = t259 * t384 + t262 * t385;
t426 = pkin(3) * t400;
t60 = pkin(4) * t402;
t423 = qJ(5) * t402;
t422 = -pkin(4) * t113 + qJ(5) * t401;
t421 = qJ(3) * t380 + t362 * t400;
t165 = -t185 + t364;
t404 = t257 * t165 + t343;
t405 = -t165 * t256 + t120;
t419 = -t259 * t404 + t262 * t405;
t418 = t276 * (-t259 * t385 + t262 * t384) + pkin(7) * (t260 * t400 + t263 * t380);
t416 = t260 * t381 + t263 * (-t259 * t405 - t262 * t404);
t415 = pkin(3) * t113;
t412 = qJ(5) * t384;
t411 = qJ(5) * t385;
t410 = qJ(6) * t382;
t315 = qJD(5) * t188;
t403 = pkin(4) * t385 - 0.2e1 * t315;
t398 = 2 * qJD(3);
t397 = pkin(3) * t380;
t396 = pkin(3) * t382;
t253 = t260 ^ 2;
t266 = qJD(1) ^ 2;
t248 = t253 * t266;
t265 = qJD(2) ^ 2;
t236 = -t248 - t265;
t323 = t260 * t266;
t304 = t263 * t323;
t231 = -qJDD(2) + t304;
t321 = t263 * t231;
t395 = pkin(7) * (-t236 * t260 + t321);
t193 = t220 * t218;
t378 = -t193 + t214;
t390 = t259 * t378;
t388 = t262 * t378;
t288 = t224 + t245;
t383 = t288 * qJ(3);
t201 = t242 * t218;
t155 = t181 + t201;
t261 = sin(qJ(1));
t264 = cos(qJ(1));
t292 = g(1) * t264 + g(2) * t261;
t351 = qJDD(1) * pkin(7);
t206 = -pkin(1) * t266 - t292 + t351;
t360 = pkin(2) * t263;
t289 = -t327 - t360;
t221 = t289 * qJD(1);
t250 = t260 * g(3);
t274 = (qJD(1) * t221 + t206) * t263 - t265 * pkin(2) - t250;
t272 = qJD(2) * t398 + t274;
t233 = pkin(3) * t318 - qJD(2) * pkin(8);
t241 = pkin(2) * t303;
t302 = qJD(3) * t318;
t244 = -0.2e1 * t302;
t254 = t263 ^ 2;
t300 = t261 * g(1) - t264 * g(2);
t280 = -qJDD(1) * pkin(1) - t300;
t115 = -t233 * t318 + t241 + t244 + (-pkin(3) * t254 - pkin(7)) * t266 + t362 * t225 - t383 + t280;
t325 = t260 * t206;
t275 = -qJDD(2) * pkin(2) - t265 * qJ(3) + t221 * t318 + qJDD(3) + t325;
t135 = t224 * pkin(3) - qJDD(2) * pkin(8) + (-pkin(3) * t313 - pkin(8) * t323 + g(3)) * t263 + t275;
t69 = t259 * t115 - t262 * t135;
t54 = pkin(4) * t378 - qJ(5) * t155 - t69;
t196 = pkin(4) * t242 - qJ(5) * t220;
t212 = t218 ^ 2;
t70 = t115 * t262 + t135 * t259;
t56 = -pkin(4) * t212 + qJ(5) * t180 - t196 * t242 + t70;
t30 = -0.2e1 * qJD(5) * t186 + t256 * t54 + t257 * t56;
t329 = t254 * t266;
t377 = t321 - (-t265 + t329) * t260;
t373 = pkin(5) * t297 - t410;
t311 = qJDD(2) * qJ(3);
t372 = t225 * pkin(3) - pkin(8) * t329 + t311;
t226 = -0.2e1 * t303 + t310;
t237 = t265 + t329;
t230 = qJDD(2) + t304;
t336 = t230 * t260;
t371 = pkin(7) * (t237 * t263 + t336) - pkin(1) * t226;
t370 = -t180 * pkin(4) - t212 * qJ(5) + t220 * t196 + qJDD(5);
t277 = (-t186 * t256 - t188 * t257) * t242;
t333 = t242 * t256;
t162 = t188 * t333;
t332 = t242 * t257;
t307 = t186 * t332;
t290 = t162 - t307;
t369 = -t259 * t277 + t262 * t290;
t279 = t256 * t297 + t307;
t291 = t186 * t333 - t257 * t297;
t368 = -t259 * t291 + t262 * t279;
t203 = t260 * t214;
t367 = t203 + t263 * (-t259 * t290 - t262 * t277);
t306 = t260 * t338;
t366 = -t306 + t263 * (-t259 * t279 - t262 * t291);
t213 = t220 ^ 2;
t363 = 2 * qJD(6);
t359 = pkin(5) * t257;
t358 = t263 * g(3);
t299 = t256 * t56 - t257 * t54;
t29 = t299 + 0.2e1 * t315;
t16 = t256 * t30 - t257 * t29;
t356 = t16 * t259;
t355 = t16 * t262;
t133 = (t398 + t233) * qJD(2) + t274 + t372;
t71 = t133 + t370;
t354 = t256 * t71;
t353 = t257 * t71;
t352 = qJ(6) * t257;
t273 = (-qJD(4) + t242) * t220 - t294;
t110 = -t155 * t262 + t259 * t273;
t346 = t110 * t260;
t172 = t193 + t214;
t341 = t172 * t259;
t340 = t172 * t262;
t331 = t242 * t259;
t330 = t242 * t262;
t328 = t259 * t133;
t324 = t260 * t226;
t322 = t262 * t133;
t228 = t248 + t329;
t320 = (t253 + t254) * t351 + pkin(1) * t228;
t314 = qJD(4) + t242;
t143 = pkin(5) * t186 - qJ(6) * t188;
t281 = t214 * qJ(6) - t186 * t143 + t242 * t363 + t30;
t25 = -pkin(5) * t364 + t281;
t278 = -t214 * pkin(5) - qJ(6) * t364 + qJDD(6) + t299;
t27 = (0.2e1 * qJD(5) + t143) * t188 + t278;
t11 = t25 * t256 - t257 * t27;
t309 = pkin(4) * t11 - pkin(5) * t27 + qJ(6) * t25;
t308 = -pkin(5) * t381 + qJ(6) * t104 + t60;
t305 = t260 * t193;
t301 = -qJ(6) * t256 - pkin(4);
t17 = t256 * t29 + t257 * t30;
t298 = -t30 + t438;
t197 = t325 + t358;
t198 = t263 * t206 - t250;
t296 = t197 * t260 + t263 * t198;
t96 = t137 * t256 + t188 * t332;
t97 = t137 * t257 - t162;
t293 = t263 * (-t259 * t97 - t262 * t96) + t306;
t40 = t259 * t70 - t262 * t69;
t287 = t259 * t69 + t262 * t70;
t285 = (-t248 + t265) * t263 + t336;
t283 = -t299 + t403;
t205 = t266 * pkin(7) - t280;
t160 = t275 + t358;
t271 = -pkin(5) * t379 - qJ(6) * t131 + t25 - t438;
t270 = t225 * pkin(2) + t205 - t241;
t269 = pkin(5) * t375 + qJ(6) * t374 - t143 * t188 - t278 + t403;
t157 = t272 + t311;
t268 = t270 + 0.2e1 * t302;
t267 = -qJD(2) * t233 + t188 * t363 - t272 - t370 - t372 - t373;
t229 = t248 - t329;
t223 = t247 + 0.2e1 * t245;
t200 = -t213 + t364;
t199 = t212 - t364;
t195 = t288 * t260;
t194 = (t225 - t303) * t263;
t191 = t213 - t212;
t190 = -t213 - t364;
t189 = t223 * t263 + t324;
t182 = -t364 - t212;
t170 = -t212 - t213;
t154 = t181 - t201;
t153 = -t314 * t218 + t282;
t150 = t314 * t220 + t294;
t141 = t190 * t262 - t341;
t125 = t182 * t259 + t388;
t58 = -t259 * t96 + t262 * t97;
t47 = t353 - t437;
t42 = t354 - t411;
t39 = -pkin(4) * t382 + t354 + t436;
t38 = (pkin(5) * t242 - (2 * qJD(6))) * t188 + t71 + t373;
t37 = -pkin(4) * t380 - t353 + t412;
t32 = (-t380 - t169) * pkin(5) + t267;
t31 = -pkin(5) * t169 + t267 + t410;
t23 = -qJ(6) * t113 + t27;
t22 = (-t113 - t364) * pkin(5) + t281;
t21 = -t256 * t32 - t352 * t380 - t411;
t20 = -pkin(5) * t348 + t257 * t31 + t437;
t19 = t257 * t32 + t301 * t380 + t412;
t18 = -t436 + t256 * t31 + (pkin(4) + t359) * t382;
t15 = pkin(4) * t16;
t14 = -pkin(4) * t71 + qJ(5) * t17;
t13 = -t16 - t423;
t12 = t25 * t257 + t256 * t27;
t9 = t17 + t422;
t8 = -t22 * t256 + t23 * t257 - t423;
t7 = t22 * t257 + t23 * t256 + t422;
t6 = -qJ(5) * t11 + (pkin(5) * t256 - t352) * t38;
t4 = t17 * t259 + t355;
t3 = qJ(5) * t12 + (t301 - t359) * t38;
t1 = t11 * t262 + t12 * t259;
t2 = [0, 0, 0, 0, 0, qJDD(1), t300, t292, 0, 0, t195, t189, t285, t194, -t377, 0, t263 * t205 - t371, -pkin(1) * t223 - t260 * t205 + t395, t296 + t320, pkin(1) * t205 + pkin(7) * t296, 0, -t285, t377, t195, t189, t194, t260 * (qJ(3) * t228 + t275) + (pkin(2) * t228 + t157 + t250) * t263 + t320, t263 * (-pkin(2) * t226 + t244 - t270) + (-t263 * t288 - t324) * qJ(3) + t371, t260 * t268 - t395 + (pkin(1) + t360) * t223 + (t223 + t288) * t327, pkin(7) * (t157 * t263 + t160 * t260) + (pkin(1) - t289) * (t268 + t383), t305 + t263 * (-t181 * t259 - t220 * t330), t260 * t191 + t263 * (t150 * t259 - t154 * t262), t260 * t155 + t263 * (-t200 * t262 - t390), -t305 + t263 * (-t180 * t262 - t218 * t331), t260 * t273 + t263 * (-t199 * t259 - t340), t203 + t263 * (t218 * t259 + t220 * t262) * t242, t260 * (pkin(3) * t125 - t69) + t263 * (pkin(3) * t150 + t322) + pkin(7) * (t125 * t260 + t150 * t263) + t276 * (t182 * t262 - t390), t260 * (pkin(3) * t141 - t70) + t263 * (pkin(3) * t153 - t328) + pkin(7) * (t141 * t260 + t153 * t263) + t276 * (-t190 * t259 - t340), pkin(3) * t346 + t263 * (pkin(3) * t170 - t287) + pkin(7) * (t170 * t263 + t346) + t276 * (t155 * t259 + t262 * t273), t276 * t287 + (pkin(3) + pkin(7)) * (t133 * t263 + t260 * t40), t293, t435, t416, t366, t439, t367, t260 * (t283 + t426) + t263 * (-t259 * t42 - t262 * t37 + t397) + t418, t260 * (t298 + t443) + t263 * (-t259 * t47 - t262 * t39 + t396) - t444, t260 * (t60 + t431) + t263 * (-t259 * t13 - t262 * t9 + t415) + t427, t260 * (pkin(3) * t4 + t15) + t263 * (pkin(3) * t71 + qJ(5) * t356 - t262 * t14) + pkin(7) * (t260 * t4 + t263 * t71) + t276 * (t17 * t262 - t356), t293, t416, -t435, t367, -t439, t366, t260 * (t269 + t426) + t263 * (-t262 * t19 - t259 * t21 + t397) + t418, t260 * (t308 + t431) + t263 * (-t259 * t8 - t262 * t7 + t415) + t427, t260 * (t271 - t443) + t263 * (-t262 * t18 - t259 * t20 - t396) + t444, t260 * (pkin(3) * t1 + t309) + t263 * (pkin(3) * t38 - t259 * t6 - t262 * t3) + pkin(7) * (t1 * t260 + t263 * t38) + t276 * (-t11 * t259 + t12 * t262); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t304, t229, t247, t304, t310, qJDD(2), -t197, -t198, 0, 0, qJDD(2), -t247, -t310, -t304, t229, t304, (-pkin(2) * t260 + qJ(3) * t263) * qJDD(1), -pkin(2) * t230 + qJ(3) * t237 + t160, -pkin(2) * t236 + (qJDD(2) - t231) * qJ(3) + t272, -pkin(2) * t160 + qJ(3) * t157, t181 * t262 - t220 * t331, -t150 * t262 - t154 * t259, -t200 * t259 + t388, -t180 * t259 + t218 * t330, t199 * t262 - t341, (-t218 * t262 + t220 * t259) * t242, qJ(3) * t150 + t125 * t362 + t328, qJ(3) * t153 + t141 * t362 + t322, qJ(3) * t170 + t110 * t362 - t40, qJ(3) * t133 + t362 * t40, t58, -t429, t419, t368, t432, t369, -t259 * t37 + t262 * t42 + t421, -t259 * t39 + t262 * t47 - t441, t262 * t13 - t259 * t9 + t428, qJ(3) * t71 - qJ(5) * t355 - t259 * t14 + t362 * t4, t58, t419, t429, t369, -t432, t368, -t259 * t19 + t262 * t21 + t421, -t259 * t7 + t262 * t8 + t428, -t259 * t18 + t262 * t20 + t441, qJ(3) * t38 + t1 * t362 - t259 * t3 + t262 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, t230, t236, t160, 0, 0, 0, 0, 0, 0, t125, t141, t110, t40, 0, 0, 0, 0, 0, 0, t400, t50, t420, t4, 0, 0, 0, 0, 0, 0, t400, t420, -t50, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t191, t155, -t193, t273, t214, -t69, -t70, 0, 0, t338, t144, t381, -t338, t104, t214, t283, t298, t60, t15, t338, t381, -t144, t214, -t104, -t338, t269, t308, t271, t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t380, t382, t113, t71, 0, 0, 0, 0, 0, 0, t380, t113, -t382, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t375, t381, t379, t27;];
tauJ_reg  = t2;
