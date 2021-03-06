% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:04
% EndTime: 2019-12-05 15:17:15
% DurationCPUTime: 6.73s
% Computational Cost: add. (20638->486), mult. (38902->764), div. (0->0), fcn. (28142->10), ass. (0->336)
t441 = sin(qJ(5));
t444 = cos(qJ(5));
t445 = cos(qJ(4));
t442 = sin(qJ(4));
t470 = qJD(3) * t442;
t391 = -t444 * t445 * qJD(3) + t441 * t470;
t393 = (t441 * t445 + t442 * t444) * qJD(3);
t356 = t393 * t391;
t462 = qJDD(4) + qJDD(5);
t502 = -t356 + t462;
t504 = t441 * t502;
t503 = t444 * t502;
t437 = sin(pkin(8));
t439 = cos(pkin(8));
t412 = t439 * g(1) + t437 * g(2);
t434 = g(3) - qJDD(1);
t436 = sin(pkin(9));
t438 = cos(pkin(9));
t383 = -t436 * t412 + t438 * t434;
t363 = t438 * t383;
t384 = -t438 * t412 - t436 * t434;
t319 = -t436 * t384 + t363;
t431 = qJD(4) + qJD(5);
t388 = t431 * t391;
t467 = qJD(3) * qJD(4);
t457 = t445 * t467;
t465 = t442 * qJDD(3);
t402 = t457 + t465;
t425 = t445 * qJDD(3);
t458 = t442 * t467;
t403 = t425 - t458;
t448 = t391 * qJD(5) - t444 * t402 - t441 * t403;
t501 = -t388 - t448;
t453 = t441 * t402 - t444 * t403;
t301 = (qJD(5) - t431) * t393 + t453;
t389 = t391 ^ 2;
t390 = t393 ^ 2;
t429 = t431 ^ 2;
t500 = qJD(3) ^ 2;
t499 = t431 * t441;
t498 = t431 * t444;
t497 = t436 * t383;
t443 = sin(qJ(3));
t446 = cos(qJ(3));
t463 = t446 * qJDD(3);
t408 = -t443 * t500 + t463;
t495 = t436 * t408;
t411 = t437 * g(1) - t439 * g(2);
t405 = -qJDD(2) + t411;
t494 = t437 * t405;
t493 = t437 * t434;
t492 = t437 * t438;
t491 = t438 * t405;
t490 = t438 * t439;
t489 = t439 * t405;
t464 = t443 * qJDD(3);
t407 = t446 * t500 + t464;
t488 = t439 * t407;
t487 = t439 * t408;
t486 = t439 * t434;
t335 = t443 * t384 + t446 * t405;
t327 = -qJDD(3) * pkin(3) - t500 * pkin(6) + t335;
t415 = qJD(4) * pkin(4) - pkin(7) * t470;
t433 = t445 ^ 2;
t427 = t433 * t500;
t281 = -t403 * pkin(4) - pkin(7) * t427 + t415 * t470 + t327;
t485 = t441 * t281;
t349 = t356 + t462;
t484 = t441 * t349;
t336 = t446 * t384 - t443 * t405;
t328 = -t500 * pkin(3) + qJDD(3) * pkin(6) + t336;
t292 = t442 * t328 - t445 * t383;
t420 = t442 * t500 * t445;
t413 = qJDD(4) + t420;
t268 = (-t402 + t457) * pkin(7) + t413 * pkin(4) - t292;
t293 = t445 * t328 + t442 * t383;
t269 = -pkin(4) * t427 + t403 * pkin(7) - qJD(4) * t415 + t293;
t213 = -t444 * t268 + t441 * t269;
t214 = t441 * t268 + t444 * t269;
t177 = -t444 * t213 + t441 * t214;
t483 = t442 * t177;
t482 = t442 * t327;
t481 = t442 * t413;
t414 = qJDD(4) - t420;
t480 = t442 * t414;
t479 = t443 * t383;
t478 = t444 * t281;
t477 = t444 * t349;
t476 = t445 * t177;
t475 = t445 * t327;
t474 = t445 * t413;
t473 = t445 * t414;
t472 = t446 * t383;
t432 = t442 ^ 2;
t471 = t432 + t433;
t469 = t432 * t500;
t466 = qJDD(3) * t436;
t461 = t443 * t356;
t460 = t446 * t356;
t459 = pkin(1) * t436 + pkin(5);
t241 = t442 * t292 + t445 * t293;
t406 = t471 * qJDD(3);
t409 = t427 + t469;
t359 = t443 * t406 + t446 * t409;
t360 = t446 * t406 - t443 * t409;
t456 = -pkin(2) * t359 - pkin(3) * t409 - pkin(6) * t406 + qJ(2) * t360 - t241;
t455 = pkin(2) * t408 + qJ(2) * t407 - t335;
t454 = -pkin(2) * t407 + qJ(2) * t408 - t336;
t178 = t441 * t213 + t444 * t214;
t320 = t438 * t384 + t497;
t362 = -t437 * t411 - t439 * t412;
t452 = t443 * t420;
t451 = t446 * t420;
t240 = t445 * t292 - t442 * t293;
t277 = t446 * t335 - t443 * t336;
t278 = t443 * t335 + t446 * t336;
t361 = t439 * t411 - t437 * t412;
t450 = t437 * t407 + t438 * t487;
t449 = t408 * t492 - t488;
t447 = qJD(4) ^ 2;
t419 = -t427 - t447;
t418 = t427 - t447;
t417 = -t447 - t469;
t416 = t447 - t469;
t410 = t427 - t469;
t404 = t425 - 0.2e1 * t458;
t401 = 0.2e1 * t457 + t465;
t400 = t471 * t467;
t397 = t436 * t407;
t386 = -t390 + t429;
t385 = t389 - t429;
t382 = t443 * qJDD(4) + t446 * t400;
t381 = t445 * t402 - t432 * t467;
t380 = t446 * qJDD(4) - t443 * t400;
t379 = -t442 * t403 - t433 * t467;
t375 = -t390 - t429;
t374 = -t442 * t417 - t473;
t373 = -t442 * t416 + t474;
t372 = t445 * t419 - t481;
t371 = t445 * t418 - t480;
t370 = t445 * t417 - t480;
t369 = -t445 * t416 - t481;
t368 = t442 * t419 + t474;
t367 = -t442 * t418 - t473;
t366 = (-t402 - t457) * t442;
t365 = (-t403 + t458) * t445;
t358 = -t442 * t401 + t445 * t404;
t357 = -t445 * t401 - t442 * t404;
t355 = t437 * t408 - t438 * t488;
t353 = -t407 * t492 - t487;
t351 = -t390 + t389;
t347 = t446 * t381 - t452;
t346 = t446 * t379 + t452;
t345 = -t443 * t381 - t451;
t344 = -t443 * t379 + t451;
t343 = t446 * t373 + t442 * t464;
t342 = t446 * t371 + t443 * t425;
t341 = -t443 * t373 + t442 * t463;
t340 = -t443 * t371 + t445 * t463;
t339 = -t429 - t389;
t338 = pkin(5) * t407 + t472;
t337 = -pkin(5) * t408 + t479;
t334 = t446 * t374 + t443 * t401;
t333 = t446 * t372 - t443 * t404;
t332 = t443 * t374 - t446 * t401;
t331 = t443 * t372 + t446 * t404;
t330 = (-t391 * t444 + t393 * t441) * t431;
t329 = (-t391 * t441 - t393 * t444) * t431;
t326 = -t389 - t390;
t325 = t446 * t358 - t443 * t410;
t324 = -t443 * t358 - t446 * t410;
t322 = -t393 * qJD(5) - t453;
t315 = t459 * t407 + t472;
t314 = t459 * t408 - t479;
t313 = t444 * t385 - t484;
t312 = -t441 * t386 + t503;
t311 = t441 * t385 + t477;
t310 = t444 * t386 + t504;
t309 = t437 * t359 + t360 * t490;
t308 = -t439 * t359 + t360 * t492;
t307 = -t441 * t375 - t477;
t306 = t444 * t375 - t484;
t305 = -t388 + t448;
t300 = (qJD(5) + t431) * t393 + t453;
t299 = t438 * t347 - t436 * t366;
t298 = t438 * t346 - t436 * t365;
t297 = t438 * t343 - t436 * t369;
t296 = t438 * t342 - t436 * t367;
t295 = t439 * t320 - t494;
t294 = t437 * t320 + t489;
t291 = -t393 * t499 - t444 * t448;
t290 = t393 * t498 - t441 * t448;
t289 = -t441 * t322 + t391 * t498;
t288 = t444 * t322 + t391 * t499;
t287 = -pkin(6) * t370 + t475;
t286 = -pkin(6) * t368 + t482;
t285 = t438 * t334 + t436 * t370;
t284 = t438 * t333 + t436 * t368;
t283 = t436 * t334 - t438 * t370;
t282 = t436 * t333 - t438 * t368;
t280 = t444 * t339 - t504;
t279 = t441 * t339 + t503;
t274 = t438 * t325 - t436 * t357;
t273 = -t442 * t329 + t445 * t330;
t272 = -t445 * t329 - t442 * t330;
t271 = -pkin(3) * t370 + t293;
t270 = -pkin(3) * t368 + t292;
t266 = t446 * t273 + t443 * t462;
t265 = -t443 * t273 + t446 * t462;
t262 = t438 * t338 + t454 * t436;
t261 = t438 * t337 + t455 * t436;
t260 = t438 * t278 + t497;
t259 = t436 * t278 - t363;
t258 = t439 * t285 + t437 * t332;
t257 = t439 * t284 + t437 * t331;
t256 = t437 * t285 - t439 * t332;
t255 = t437 * t284 - t439 * t331;
t254 = -pkin(2) * t331 - pkin(3) * t404 - pkin(6) * t372 + t475;
t253 = -pkin(2) * t332 + pkin(3) * t401 - pkin(6) * t374 - t482;
t252 = -t442 * t311 + t445 * t313;
t251 = -t442 * t310 + t445 * t312;
t250 = -t445 * t311 - t442 * t313;
t249 = -t445 * t310 - t442 * t312;
t248 = -t442 * t306 + t445 * t307;
t247 = t445 * t306 + t442 * t307;
t246 = -t301 * t444 - t441 * t305;
t245 = -t444 * t300 - t441 * t501;
t244 = -t301 * t441 + t444 * t305;
t243 = -t441 * t300 + t444 * t501;
t242 = -pkin(7) * t306 + t478;
t238 = -t442 * t290 + t445 * t291;
t237 = -t442 * t288 + t445 * t289;
t236 = -t445 * t290 - t442 * t291;
t235 = -t445 * t288 - t442 * t289;
t234 = -pkin(7) * t279 + t485;
t233 = -t442 * t279 + t445 * t280;
t232 = t445 * t279 + t442 * t280;
t231 = -pkin(5) * t359 + t446 * t240;
t230 = t446 * t238 + t461;
t229 = t446 * t237 - t461;
t228 = -t443 * t238 + t460;
t227 = -t443 * t237 - t460;
t226 = t446 * t241 + t443 * t327;
t225 = t443 * t241 - t446 * t327;
t223 = t446 * t252 - t443 * t301;
t222 = t446 * t251 - t443 * t305;
t221 = -t443 * t252 - t446 * t301;
t220 = -t443 * t251 - t446 * t305;
t219 = t439 * t260 - t277 * t437;
t218 = t437 * t260 + t277 * t439;
t217 = t438 * t266 - t436 * t272;
t216 = t446 * t248 + t443 * t501;
t215 = t443 * t248 - t446 * t501;
t211 = -pkin(4) * t501 + pkin(7) * t307 + t485;
t210 = -t443 * t240 - t459 * t360;
t209 = -pkin(5) * t332 - t443 * t271 + t446 * t287;
t208 = -pkin(5) * t331 - t443 * t270 + t446 * t286;
t207 = -pkin(1) * t259 + pkin(2) * t383 - pkin(5) * t278;
t206 = -pkin(4) * t300 + pkin(7) * t280 - t478;
t205 = t446 * t233 + t443 * t300;
t204 = t443 * t233 - t446 * t300;
t203 = -t442 * t244 + t445 * t246;
t202 = -t442 * t243 + t445 * t245;
t201 = t445 * t244 + t442 * t246;
t200 = -t445 * t243 - t442 * t245;
t199 = -qJ(2) * t259 - (pkin(2) * t436 - pkin(5) * t438) * t277;
t198 = -pkin(1) * t283 + pkin(2) * t370 - pkin(5) * t334 - t446 * t271 - t443 * t287;
t197 = -pkin(1) * t282 + pkin(2) * t368 - pkin(5) * t333 - t446 * t270 - t443 * t286;
t196 = t446 * t202 - t443 * t351;
t195 = -t443 * t202 - t446 * t351;
t194 = t446 * t203 + t443 * t326;
t193 = t443 * t203 - t446 * t326;
t192 = t438 * t223 - t436 * t250;
t191 = t438 * t222 - t436 * t249;
t190 = t438 * t230 - t436 * t236;
t189 = t438 * t229 - t436 * t235;
t188 = t438 * t216 + t436 * t247;
t187 = t436 * t216 - t438 * t247;
t186 = t438 * t226 - t240 * t436;
t185 = t436 * t226 + t240 * t438;
t184 = t438 * t205 + t436 * t232;
t183 = t436 * t205 - t438 * t232;
t182 = -pkin(2) * t225 + pkin(3) * t327 - pkin(6) * t241;
t181 = -pkin(3) * t247 - pkin(4) * t306 + t214;
t180 = t438 * t231 - t456 * t436;
t179 = -pkin(3) * t201 - pkin(4) * t244;
t176 = -qJ(2) * t283 + t438 * t209 - t436 * t253;
t175 = -qJ(2) * t282 + t438 * t208 - t436 * t254;
t174 = -pkin(3) * t232 - pkin(4) * t279 + t213;
t173 = -pkin(6) * t247 - t442 * t211 + t445 * t242;
t172 = -pkin(4) * t281 + pkin(7) * t178;
t171 = -pkin(5) * t225 - (pkin(3) * t443 - pkin(6) * t446) * t240;
t170 = -pkin(6) * t232 - t442 * t206 + t445 * t234;
t169 = -pkin(7) * t244 - t177;
t168 = t439 * t186 + t437 * t225;
t167 = t437 * t186 - t439 * t225;
t166 = t439 * t188 + t437 * t215;
t165 = t437 * t188 - t439 * t215;
t164 = t438 * t196 - t436 * t200;
t163 = -pkin(4) * t326 + pkin(7) * t246 + t178;
t162 = t438 * t194 + t436 * t201;
t161 = t436 * t194 - t438 * t201;
t160 = t439 * t184 + t437 * t204;
t159 = t437 * t184 - t439 * t204;
t158 = -pkin(2) * t215 + pkin(3) * t501 - pkin(6) * t248 - t445 * t211 - t442 * t242;
t157 = -pkin(2) * t204 + pkin(3) * t300 - pkin(6) * t233 - t445 * t206 - t442 * t234;
t156 = t445 * t178 - t483;
t155 = t442 * t178 + t476;
t154 = t446 * t156 + t443 * t281;
t153 = t443 * t156 - t446 * t281;
t152 = t439 * t162 + t437 * t193;
t151 = t437 * t162 - t439 * t193;
t150 = -pkin(1) * t185 - pkin(5) * t226 - (pkin(3) * t446 + pkin(6) * t443 + pkin(2)) * t240;
t149 = -pkin(5) * t215 + t446 * t173 - t443 * t181;
t148 = -pkin(5) * t204 + t446 * t170 - t443 * t174;
t147 = -qJ(2) * t185 + t438 * t171 - t436 * t182;
t146 = -pkin(3) * t155 - pkin(4) * t177;
t145 = -pkin(6) * t201 - t442 * t163 + t445 * t169;
t144 = -pkin(1) * t187 + pkin(2) * t247 - pkin(5) * t216 - t443 * t173 - t446 * t181;
t143 = -pkin(1) * t183 + pkin(2) * t232 - pkin(5) * t205 - t443 * t170 - t446 * t174;
t142 = -pkin(6) * t155 - pkin(7) * t476 - t442 * t172;
t141 = t438 * t154 + t436 * t155;
t140 = t436 * t154 - t438 * t155;
t139 = -pkin(2) * t193 + pkin(3) * t326 - pkin(6) * t203 - t445 * t163 - t442 * t169;
t138 = -pkin(5) * t193 + t446 * t145 - t443 * t179;
t137 = -qJ(2) * t187 + t438 * t149 - t436 * t158;
t136 = -qJ(2) * t183 + t438 * t148 - t436 * t157;
t135 = t439 * t141 + t437 * t153;
t134 = t437 * t141 - t439 * t153;
t133 = -pkin(2) * t153 + pkin(3) * t281 - pkin(6) * t156 + pkin(7) * t483 - t445 * t172;
t132 = -pkin(1) * t161 + pkin(2) * t201 - pkin(5) * t194 - t443 * t145 - t446 * t179;
t131 = -pkin(5) * t153 + t446 * t142 - t443 * t146;
t130 = -qJ(2) * t161 + t438 * t138 - t436 * t139;
t129 = -pkin(1) * t140 + pkin(2) * t155 - pkin(5) * t154 - t443 * t142 - t446 * t146;
t128 = -qJ(2) * t140 + t438 * t131 - t436 * t133;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t362, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, 0, 0, 0, 0, 0, 0, t355, -t450, 0, t219, 0, 0, 0, 0, 0, 0, t257, t258, t309, t168, 0, 0, 0, 0, 0, 0, t160, t166, t152, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t361, 0, 0, 0, 0, 0, 0, 0, 0, 0, t294, 0, 0, 0, 0, 0, 0, t353, -t449, 0, t218, 0, 0, 0, 0, 0, 0, t255, t256, t308, t167, 0, 0, 0, 0, 0, 0, t159, t165, t151, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t434, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t319, 0, 0, 0, 0, 0, 0, -t397, -t495, 0, t259, 0, 0, 0, 0, 0, 0, t282, t283, t436 * t360, t185, 0, 0, 0, 0, 0, 0, t183, t187, t161, t140; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t493, -t486, -t361, -qJ(1) * t361, 0, 0, 0, 0, 0, 0, -t437 * t383 - t436 * t489, -t437 * t384 - t438 * t489, t439 * t319, -qJ(1) * t294 - (pkin(1) * t437 - qJ(2) * t439) * t319, 0, 0, t450, 0, t355, t439 * t466, -qJ(1) * t353 + t439 * t261 - t437 * t315, qJ(1) * t449 + t439 * t262 - t437 * t314, t277 * t490 + t278 * t437, -qJ(1) * t218 + t439 * t199 - t437 * t207, t439 * t299 - t437 * t345, t439 * t274 - t437 * t324, t439 * t297 - t437 * t341, t439 * t298 - t437 * t344, t439 * t296 - t437 * t340, -t437 * t380 + t382 * t490, -qJ(1) * t255 + t439 * t175 - t437 * t197, -qJ(1) * t256 + t439 * t176 - t437 * t198, -qJ(1) * t308 + t439 * t180 - t437 * t210, -qJ(1) * t167 + t439 * t147 - t437 * t150, t439 * t190 - t437 * t228, t439 * t164 - t437 * t195, t439 * t191 - t437 * t220, t439 * t189 - t437 * t227, t439 * t192 - t437 * t221, t439 * t217 - t437 * t265, -qJ(1) * t159 + t439 * t136 - t437 * t143, -qJ(1) * t165 + t439 * t137 - t437 * t144, -qJ(1) * t151 + t439 * t130 - t437 * t132, -qJ(1) * t134 + t439 * t128 - t437 * t129; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t486, -t493, t362, qJ(1) * t362, 0, 0, 0, 0, 0, 0, t439 * t383 - t436 * t494, t439 * t384 - t437 * t491, t437 * t319, qJ(1) * t295 - (-pkin(1) * t439 - qJ(2) * t437) * t319, 0, 0, t449, 0, t353, t437 * t466, qJ(1) * t355 + t437 * t261 + t439 * t315, -qJ(1) * t450 + t437 * t262 + t439 * t314, t277 * t492 - t278 * t439, qJ(1) * t219 + t437 * t199 + t439 * t207, t437 * t299 + t439 * t345, t437 * t274 + t439 * t324, t437 * t297 + t439 * t341, t437 * t298 + t439 * t344, t437 * t296 + t439 * t340, t439 * t380 + t382 * t492, qJ(1) * t257 + t437 * t175 + t439 * t197, qJ(1) * t258 + t437 * t176 + t439 * t198, qJ(1) * t309 + t437 * t180 + t439 * t210, qJ(1) * t168 + t437 * t147 + t439 * t150, t437 * t190 + t439 * t228, t437 * t164 + t439 * t195, t437 * t191 + t439 * t220, t437 * t189 + t439 * t227, t437 * t192 + t439 * t221, t437 * t217 + t439 * t265, qJ(1) * t160 + t437 * t136 + t439 * t143, qJ(1) * t166 + t437 * t137 + t439 * t144, qJ(1) * t152 + t437 * t130 + t439 * t132, qJ(1) * t135 + t437 * t128 + t439 * t129; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t411, t412, 0, 0, 0, 0, 0, 0, 0, 0, t491, -t436 * t405, t320, pkin(1) * t405 + qJ(2) * t320, 0, 0, t495, 0, -t397, -t438 * qJDD(3), -pkin(1) * t408 + t436 * t337 - t438 * t455, pkin(1) * t407 + t436 * t338 - t438 * t454, t436 * t277, qJ(2) * t260 - (-pkin(2) * t438 - pkin(5) * t436 - pkin(1)) * t277, t436 * t347 + t438 * t366, t436 * t325 + t438 * t357, t436 * t343 + t438 * t369, t436 * t346 + t438 * t365, t436 * t342 + t438 * t367, t436 * t382, -pkin(1) * t331 + qJ(2) * t284 + t436 * t208 + t438 * t254, -pkin(1) * t332 + qJ(2) * t285 + t436 * t209 + t438 * t253, -pkin(1) * t359 + t436 * t231 + t438 * t456, -pkin(1) * t225 + qJ(2) * t186 + t436 * t171 + t438 * t182, t436 * t230 + t438 * t236, t436 * t196 + t438 * t200, t436 * t222 + t438 * t249, t436 * t229 + t438 * t235, t436 * t223 + t438 * t250, t436 * t266 + t438 * t272, -pkin(1) * t204 + qJ(2) * t184 + t436 * t148 + t438 * t157, -pkin(1) * t215 + qJ(2) * t188 + t436 * t149 + t438 * t158, -pkin(1) * t193 + qJ(2) * t162 + t436 * t138 + t438 * t139, -pkin(1) * t153 + qJ(2) * t141 + t436 * t131 + t438 * t133;];
tauB_reg = t1;
