% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:40
% EndTime: 2019-12-05 18:14:51
% DurationCPUTime: 8.89s
% Computational Cost: add. (33385->383), mult. (48058->556), div. (0->0), fcn. (28562->10), ass. (0->257)
t410 = qJD(1) + qJD(3);
t407 = qJD(4) + t410;
t404 = t407 ^ 2;
t421 = cos(qJ(4));
t409 = qJDD(1) + qJDD(3);
t405 = qJDD(4) + t409;
t417 = sin(qJ(4));
t453 = t417 * t405;
t367 = t421 * t404 + t453;
t446 = t421 * t405;
t370 = t417 * t404 - t446;
t418 = sin(qJ(3));
t422 = cos(qJ(3));
t311 = t422 * t367 - t418 * t370;
t413 = g(1) - qJDD(2);
t348 = pkin(7) * t367 - t421 * t413;
t478 = pkin(7) * t370 - t417 * t413;
t260 = pkin(6) * t311 + t422 * t348 - t418 * t478;
t316 = t418 * t367 + t422 * t370;
t414 = sin(pkin(9));
t415 = cos(pkin(9));
t269 = t415 * t311 - t414 * t316;
t489 = pkin(6) * t316 + t418 * t348 + t422 * t478;
t204 = qJ(2) * t269 + t415 * t260 - t414 * t489;
t419 = sin(qJ(1));
t423 = cos(qJ(1));
t273 = t414 * t311 + t415 * t316;
t500 = qJ(2) * t273 + t414 * t260 + t415 * t489;
t502 = t419 * t269 + t423 * t273;
t509 = pkin(5) * t502 + t419 * t204 + t423 * t500;
t226 = t423 * t269 - t419 * t273;
t508 = pkin(5) * t226 + t423 * t204 - t419 * t500;
t397 = t423 * g(2) + t419 * g(3);
t384 = qJDD(1) * pkin(1) + t397;
t396 = t419 * g(2) - t423 * g(3);
t425 = qJD(1) ^ 2;
t385 = -t425 * pkin(1) + t396;
t331 = -t415 * t384 + t414 * t385;
t327 = qJDD(1) * pkin(2) - t331;
t332 = t414 * t384 + t415 * t385;
t328 = -t425 * pkin(2) + t332;
t288 = t418 * t327 + t422 * t328;
t408 = t410 ^ 2;
t280 = -t408 * pkin(3) + t288;
t428 = t422 * t327 - t418 * t328;
t426 = t409 * pkin(3) + t428;
t236 = t417 * t280 - t421 * t426;
t237 = t421 * t280 + t417 * t426;
t437 = t417 * t236 + t421 * t237;
t195 = t421 * t236 - t417 * t237;
t444 = t422 * t195;
t169 = -t418 * t437 + t444;
t451 = t418 * t195;
t481 = t422 * t437 + t451;
t156 = t415 * t169 - t414 * t481;
t496 = t414 * t169 + t415 * t481;
t142 = t419 * t156 + t423 * t496;
t141 = t423 * t156 - t419 * t496;
t378 = t422 * t408 + t418 * t409;
t381 = t418 * t408 - t422 * t409;
t320 = t415 * t378 - t414 * t381;
t354 = pkin(6) * t378 - t422 * t413;
t479 = pkin(6) * t381 - t418 * t413;
t268 = qJ(2) * t320 + t415 * t354 - t414 * t479;
t324 = t414 * t378 + t415 * t381;
t488 = qJ(2) * t324 + t414 * t354 + t415 * t479;
t490 = t419 * t320 + t423 * t324;
t503 = pkin(5) * t490 + t419 * t268 + t423 * t488;
t281 = t423 * t320 - t419 * t324;
t501 = pkin(5) * t281 + t423 * t268 - t419 * t488;
t436 = t422 * t288 - t418 * t428;
t241 = -t418 * t288 - t422 * t428;
t458 = t415 * t241;
t199 = -t414 * t436 + t458;
t459 = t414 * t241;
t480 = t415 * t436 + t459;
t174 = t419 * t199 + t423 * t480;
t173 = t423 * t199 - t419 * t480;
t435 = t414 * t331 + t415 * t332;
t291 = t415 * t331 - t414 * t332;
t442 = t423 * t291;
t243 = -t419 * t435 + t442;
t449 = t419 * t291;
t244 = t423 * t435 + t449;
t386 = t414 * qJDD(1) + t415 * t425;
t359 = qJ(2) * t386 - t415 * t413;
t387 = t415 * qJDD(1) - t414 * t425;
t427 = t419 * t386 - t423 * t387;
t429 = -qJ(2) * t387 - t414 * t413;
t476 = pkin(5) * t427 + t419 * t359 + t423 * t429;
t340 = t423 * t386 + t419 * t387;
t475 = pkin(5) * t340 + t423 * t359 - t419 * t429;
t462 = pkin(1) * t413;
t461 = pkin(2) * t413;
t416 = sin(qJ(5));
t411 = t416 ^ 2;
t460 = t411 * t404;
t233 = -t405 * pkin(4) - t404 * pkin(8) + t236;
t457 = t416 * t233;
t420 = cos(qJ(5));
t392 = t420 * t404 * t416;
t382 = qJDD(5) + t392;
t456 = t416 * t382;
t383 = qJDD(5) - t392;
t455 = t416 * t383;
t454 = t416 * t405;
t448 = t420 * t233;
t447 = t420 * t383;
t398 = t420 * t405;
t234 = -t404 * pkin(4) + t405 * pkin(8) + t237;
t225 = t420 * t234 - t416 * t413;
t412 = t420 ^ 2;
t441 = t411 + t412;
t440 = qJD(5) * t407;
t439 = t416 * t440;
t438 = t420 * t440;
t224 = t416 * t234 + t420 * t413;
t189 = t416 * t224 + t420 * t225;
t432 = t417 * t392;
t431 = t421 * t392;
t393 = t419 * qJDD(1) + t423 * t425;
t430 = pkin(5) * t393 - t423 * g(1);
t188 = t420 * t224 - t416 * t225;
t350 = t423 * t396 - t419 * t397;
t349 = -t419 * t396 - t423 * t397;
t424 = qJD(5) ^ 2;
t399 = t412 * t404;
t394 = -t423 * qJDD(1) + t419 * t425;
t391 = -t399 - t424;
t390 = t399 - t424;
t389 = -t424 - t460;
t388 = t424 - t460;
t374 = t420 * t382;
t373 = -pkin(5) * t394 + t419 * g(1);
t372 = t399 - t460;
t371 = t399 + t460;
t366 = t441 * t405;
t364 = t398 - 0.2e1 * t439;
t363 = t398 - t439;
t362 = t438 + t454;
t361 = 0.2e1 * t438 + t454;
t360 = t441 * t440;
t344 = t417 * qJDD(5) + t421 * t360;
t343 = -t421 * qJDD(5) + t417 * t360;
t338 = -t416 * t389 - t447;
t337 = -t416 * t388 + t374;
t336 = t420 * t391 - t456;
t335 = t420 * t390 - t455;
t334 = t420 * t389 - t455;
t333 = t416 * t391 + t374;
t330 = t420 * t362 - t411 * t440;
t329 = -t416 * t363 - t412 * t440;
t314 = t421 * t366 - t417 * t371;
t310 = t417 * t366 + t421 * t371;
t309 = -t416 * t361 + t420 * t364;
t308 = t421 * t337 + t416 * t453;
t307 = t421 * t335 + t417 * t398;
t306 = t417 * t337 - t416 * t446;
t305 = t417 * t335 - t420 * t446;
t304 = t421 * t330 - t432;
t303 = t421 * t329 + t432;
t302 = t417 * t330 + t431;
t301 = t417 * t329 - t431;
t300 = t421 * t338 + t417 * t361;
t299 = t421 * t336 - t417 * t364;
t298 = t417 * t338 - t421 * t361;
t297 = t417 * t336 + t421 * t364;
t296 = -t418 * t343 + t422 * t344;
t295 = t422 * t343 + t418 * t344;
t294 = t421 * t309 - t417 * t372;
t293 = t417 * t309 + t421 * t372;
t286 = qJ(2) * t435 + t462;
t276 = -t418 * t310 + t422 * t314;
t275 = t422 * t310 + t418 * t314;
t264 = -t418 * t306 + t422 * t308;
t263 = -t418 * t305 + t422 * t307;
t262 = t422 * t306 + t418 * t308;
t261 = t422 * t305 + t418 * t307;
t256 = -t418 * t302 + t422 * t304;
t255 = -t418 * t301 + t422 * t303;
t254 = t422 * t302 + t418 * t304;
t253 = t422 * t301 + t418 * t303;
t252 = -t418 * t298 + t422 * t300;
t251 = -t418 * t297 + t422 * t299;
t250 = t422 * t298 + t418 * t300;
t249 = t422 * t297 + t418 * t299;
t248 = -t414 * t295 + t415 * t296;
t247 = t415 * t295 + t414 * t296;
t246 = -t418 * t293 + t422 * t294;
t245 = t422 * t293 + t418 * t294;
t238 = pkin(6) * t436 + t461;
t231 = -t414 * t275 + t415 * t276;
t230 = t415 * t275 + t414 * t276;
t222 = -t414 * t262 + t415 * t264;
t221 = -t414 * t261 + t415 * t263;
t220 = t415 * t262 + t414 * t264;
t219 = t415 * t261 + t414 * t263;
t218 = -t414 * t254 + t415 * t256;
t217 = -t414 * t253 + t415 * t255;
t216 = t415 * t254 + t414 * t256;
t215 = t415 * t253 + t414 * t255;
t214 = -pkin(8) * t334 + t448;
t213 = -pkin(8) * t333 + t457;
t212 = -pkin(4) * t334 + t225;
t211 = -pkin(4) * t333 + t224;
t210 = -t414 * t250 + t415 * t252;
t209 = -t414 * t249 + t415 * t251;
t208 = t415 * t250 + t414 * t252;
t207 = t415 * t249 + t414 * t251;
t206 = -t414 * t245 + t415 * t246;
t205 = t415 * t245 + t414 * t246;
t192 = pkin(3) * t413 + pkin(7) * t437;
t191 = -t419 * t230 + t423 * t231;
t190 = -t423 * t230 - t419 * t231;
t186 = -pkin(7) * t310 + t421 * t188;
t185 = pkin(7) * t314 + t417 * t188;
t184 = -t419 * t208 + t423 * t210;
t183 = -t419 * t207 + t423 * t209;
t182 = -t423 * t208 - t419 * t210;
t181 = -t423 * t207 - t419 * t209;
t180 = -pkin(7) * t298 - t417 * t212 + t421 * t214;
t179 = -pkin(7) * t297 - t417 * t211 + t421 * t213;
t178 = -pkin(3) * t334 + pkin(7) * t300 + t421 * t212 + t417 * t214;
t177 = -pkin(3) * t333 + pkin(7) * t299 + t421 * t211 + t417 * t213;
t176 = t421 * t189 + t417 * t233;
t175 = t417 * t189 - t421 * t233;
t172 = pkin(6) * t458 + qJ(2) * t199 - t414 * t238;
t171 = pkin(6) * t459 + qJ(2) * t480 + t415 * t238 + t462;
t166 = -pkin(6) * t275 - t418 * t185 + t422 * t186;
t165 = pkin(6) * t276 + t422 * t185 + t418 * t186;
t164 = -t418 * t175 + t422 * t176;
t163 = t422 * t175 + t418 * t176;
t162 = -pkin(6) * t250 - t418 * t178 + t422 * t180;
t161 = -pkin(6) * t249 - t418 * t177 + t422 * t179;
t160 = -pkin(2) * t334 + pkin(6) * t252 + t422 * t178 + t418 * t180;
t159 = -pkin(2) * t333 + pkin(6) * t251 + t422 * t177 + t418 * t179;
t158 = -pkin(7) * t175 - (pkin(4) * t417 - pkin(8) * t421) * t188;
t153 = pkin(6) * t169 + pkin(7) * t444 - t418 * t192;
t152 = pkin(6) * t481 + pkin(7) * t451 + t422 * t192 + t461;
t151 = pkin(7) * t176 - (-pkin(4) * t421 - pkin(8) * t417 - pkin(3)) * t188;
t150 = -qJ(2) * t230 - t414 * t165 + t415 * t166;
t149 = qJ(2) * t231 + t415 * t165 + t414 * t166;
t148 = -t414 * t163 + t415 * t164;
t147 = t415 * t163 + t414 * t164;
t146 = -qJ(2) * t208 - t414 * t160 + t415 * t162;
t145 = -qJ(2) * t207 - t414 * t159 + t415 * t161;
t144 = -pkin(1) * t334 + qJ(2) * t210 + t415 * t160 + t414 * t162;
t143 = -pkin(1) * t333 + qJ(2) * t209 + t415 * t159 + t414 * t161;
t140 = qJ(2) * t156 - t414 * t152 + t415 * t153;
t139 = qJ(2) * t496 + t415 * t152 + t414 * t153 + t462;
t138 = -pkin(6) * t163 - t418 * t151 + t422 * t158;
t137 = -t419 * t147 + t423 * t148;
t136 = -t423 * t147 - t419 * t148;
t135 = pkin(2) * t188 + pkin(6) * t164 + t422 * t151 + t418 * t158;
t134 = -qJ(2) * t147 - t414 * t135 + t415 * t138;
t133 = pkin(1) * t188 + qJ(2) * t148 + t415 * t135 + t414 * t138;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t413, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t413, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t413, 0, 0, 0, 0, 0, 0, t333, t334, 0, -t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t394, t393, 0, t349, 0, 0, 0, 0, 0, 0, t427, t340, 0, t243, 0, 0, 0, 0, 0, 0, t490, t281, 0, t173, 0, 0, 0, 0, 0, 0, t502, t226, 0, t141, 0, 0, 0, 0, 0, 0, t181, t182, t190, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t393, t394, 0, t350, 0, 0, 0, 0, 0, 0, -t340, t427, 0, t244, 0, 0, 0, 0, 0, 0, -t281, t490, 0, t174, 0, 0, 0, 0, 0, 0, -t226, t502, 0, t142, 0, 0, 0, 0, 0, 0, t183, t184, t191, t137; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), t397, -t396, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t387 - t331, -pkin(1) * t386 - t332, 0, -pkin(1) * t291, 0, 0, 0, 0, 0, t409, -pkin(1) * t324 - pkin(2) * t381 + t428, -pkin(1) * t320 - pkin(2) * t378 - t288, 0, -pkin(1) * t199 - pkin(2) * t241, 0, 0, 0, 0, 0, t405, -pkin(1) * t273 - pkin(2) * t316 - pkin(3) * t370 - t236, -pkin(1) * t269 - pkin(2) * t311 - pkin(3) * t367 - t237, 0, -pkin(1) * t156 - pkin(2) * t169 - pkin(3) * t195, (t362 + t438) * t416, t420 * t361 + t416 * t364, t420 * t388 + t456, (t363 - t439) * t420, t416 * t390 + t447, 0, pkin(1) * t207 + pkin(2) * t249 + pkin(3) * t297 + pkin(4) * t364 + pkin(8) * t336 - t448, pkin(1) * t208 + pkin(2) * t250 + pkin(3) * t298 - pkin(4) * t361 + pkin(8) * t338 + t457, pkin(1) * t230 + pkin(2) * t275 + pkin(3) * t310 + pkin(4) * t371 + pkin(8) * t366 + t189, pkin(1) * t147 + pkin(2) * t163 + pkin(3) * t175 - pkin(4) * t233 + pkin(8) * t189; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t393, 0, t394, 0, t430, t373, -t350, -pkin(5) * t350, 0, 0, -t340, 0, t427, 0, t475, -t476, -t244, -pkin(5) * t244 - qJ(2) * t449 - t423 * t286, 0, 0, -t281, 0, t490, 0, t501, -t503, -t174, -pkin(5) * t174 - t423 * t171 - t419 * t172, 0, 0, -t226, 0, t502, 0, t508, -t509, -t142, -pkin(5) * t142 - t423 * t139 - t419 * t140, -t423 * t216 - t419 * t218, -t423 * t205 - t419 * t206, -t423 * t220 - t419 * t222, -t423 * t215 - t419 * t217, -t423 * t219 - t419 * t221, -t423 * t247 - t419 * t248, -pkin(5) * t183 - t423 * t143 - t419 * t145, -pkin(5) * t184 - t423 * t144 - t419 * t146, -pkin(5) * t191 - t423 * t149 - t419 * t150, -pkin(5) * t137 - t423 * t133 - t419 * t134; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t394, 0, -t393, 0, -t373, t430, t349, pkin(5) * t349, 0, 0, -t427, 0, -t340, 0, t476, t475, t243, pkin(5) * t243 + qJ(2) * t442 - t419 * t286, 0, 0, -t490, 0, -t281, 0, t503, t501, t173, pkin(5) * t173 - t419 * t171 + t423 * t172, 0, 0, -t502, 0, -t226, 0, t509, t508, t141, pkin(5) * t141 - t419 * t139 + t423 * t140, -t419 * t216 + t423 * t218, -t419 * t205 + t423 * t206, -t419 * t220 + t423 * t222, -t419 * t215 + t423 * t217, -t419 * t219 + t423 * t221, -t419 * t247 + t423 * t248, pkin(5) * t181 - t419 * t143 + t423 * t145, pkin(5) * t182 - t419 * t144 + t423 * t146, pkin(5) * t190 - t419 * t149 + t423 * t150, pkin(5) * t136 - t419 * t133 + t423 * t134;];
tauB_reg = t1;
