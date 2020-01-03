% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:10
% EndTime: 2019-12-31 19:06:18
% DurationCPUTime: 5.41s
% Computational Cost: add. (21162->417), mult. (28284->596), div. (0->0), fcn. (15183->8), ass. (0->269)
t438 = sin(qJ(3));
t432 = -qJD(1) + qJD(3);
t428 = t432 ^ 2;
t430 = qJDD(1) - qJDD(3);
t442 = cos(qJ(3));
t464 = t442 * t430;
t451 = t428 * t438 + t464;
t372 = pkin(6) * t451 + g(3) * t438;
t439 = sin(qJ(1));
t443 = cos(qJ(1));
t471 = t438 * t430;
t454 = -t428 * t442 + t471;
t487 = t439 * t454 + t443 * t451;
t488 = -pkin(6) * t454 + g(3) * t442;
t497 = -pkin(5) * t487 + t372 * t443 - t439 * t488;
t345 = t439 * t451 - t443 * t454;
t496 = -pkin(5) * t345 + t372 * t439 + t443 * t488;
t436 = sin(qJ(5));
t440 = cos(qJ(5));
t441 = cos(qJ(4));
t437 = sin(qJ(4));
t480 = t432 * t437;
t378 = -t432 * t440 * t441 + t436 * t480;
t380 = (t436 * t441 + t437 * t440) * t432;
t342 = t380 * t378;
t460 = qJDD(4) + qJDD(5);
t490 = -t342 + t460;
t495 = t436 * t490;
t494 = t440 * t490;
t433 = qJDD(1) * qJ(2);
t413 = t443 * g(1) + t439 * g(2);
t450 = (2 * qJD(2) * qJD(1)) - t413;
t448 = t433 + t450;
t485 = (pkin(1) + pkin(2));
t486 = qJD(1) ^ 2;
t371 = -(t485 * t486) + t448;
t412 = t439 * g(1) - t443 * g(2);
t449 = -qJDD(2) + t412;
t446 = -(qJ(2) * t486) - t449;
t445 = -(qJDD(1) * t485) + t446;
t325 = t438 * t371 - t442 * t445;
t326 = t371 * t442 + t438 * t445;
t282 = t325 * t442 - t326 * t438;
t283 = t325 * t438 + t326 * t442;
t243 = t282 * t443 + t283 * t439;
t491 = t282 * t439 - t283 * t443;
t431 = qJD(4) + qJD(5);
t367 = t431 * t378;
t462 = qJD(4) * t432;
t456 = t441 * t462;
t472 = t437 * t430;
t390 = t456 - t472;
t416 = t441 * t430;
t457 = t437 * t462;
t391 = -t416 - t457;
t447 = qJD(5) * t378 - t390 * t440 - t391 * t436;
t489 = -t367 - t447;
t455 = t436 * t390 - t391 * t440;
t296 = (qJD(5) - t431) * t380 + t455;
t376 = t378 ^ 2;
t377 = t380 ^ 2;
t427 = t431 ^ 2;
t483 = qJDD(1) * pkin(1);
t482 = t431 * t436;
t481 = t431 * t440;
t434 = t437 ^ 2;
t479 = t434 * t428;
t435 = t441 ^ 2;
t418 = t435 * t428;
t317 = pkin(3) * t430 - pkin(7) * t428 + t325;
t403 = qJD(4) * pkin(4) - pkin(8) * t480;
t284 = -pkin(4) * t391 - pkin(8) * t418 + t403 * t480 + t317;
t478 = t436 * t284;
t333 = t342 + t460;
t477 = t436 * t333;
t318 = -pkin(3) * t428 - pkin(7) * t430 + t326;
t309 = -g(3) * t441 + t437 * t318;
t411 = t441 * t428 * t437;
t401 = qJDD(4) + t411;
t278 = (-t390 + t456) * pkin(8) + t401 * pkin(4) - t309;
t310 = g(3) * t437 + t318 * t441;
t279 = -pkin(4) * t418 + pkin(8) * t391 - qJD(4) * t403 + t310;
t241 = -t278 * t440 + t279 * t436;
t242 = t278 * t436 + t279 * t440;
t211 = -t241 * t440 + t242 * t436;
t476 = t437 * t211;
t475 = t437 * t317;
t474 = t437 * t401;
t402 = qJDD(4) - t411;
t473 = t437 * t402;
t470 = t440 * t284;
t469 = t440 * t333;
t468 = t441 * t211;
t467 = t441 * t317;
t466 = t441 * t401;
t465 = t441 * t402;
t463 = t434 + t435;
t459 = t438 * t342;
t458 = t442 * t342;
t212 = t241 * t436 + t242 * t440;
t382 = -(pkin(1) * t486) + t448;
t385 = -t446 + t483;
t336 = t382 * t443 - t385 * t439;
t362 = -t412 * t439 - t413 * t443;
t453 = t438 * t411;
t452 = t442 * t411;
t404 = qJDD(1) * t439 + t443 * t486;
t388 = -pkin(5) * t404 + g(3) * t443;
t405 = qJDD(1) * t443 - t439 * t486;
t387 = pkin(5) * t405 + g(3) * t439;
t264 = t309 * t441 - t310 * t437;
t265 = t437 * t309 + t441 * t310;
t335 = t382 * t439 + t385 * t443;
t361 = t412 * t443 - t413 * t439;
t444 = qJD(4) ^ 2;
t409 = -t418 - t444;
t408 = t418 - t444;
t407 = -t444 - t479;
t406 = t444 - t479;
t399 = t418 - t479;
t398 = t418 + t479;
t393 = t463 * t430;
t392 = -t416 - 0.2e1 * t457;
t389 = 0.2e1 * t456 - t472;
t386 = t463 * t462;
t365 = -t377 + t427;
t364 = t376 - t427;
t360 = qJDD(4) * t438 + t386 * t442;
t359 = qJDD(4) * t442 - t386 * t438;
t358 = t390 * t441 - t434 * t462;
t357 = -t391 * t437 - t435 * t462;
t356 = -t377 - t427;
t355 = -t407 * t437 - t465;
t354 = -t406 * t437 + t466;
t353 = t409 * t441 - t474;
t352 = t408 * t441 - t473;
t351 = t407 * t441 - t473;
t350 = t409 * t437 + t466;
t347 = -t393 * t442 - t398 * t438;
t344 = -t393 * t438 + t398 * t442;
t343 = -t389 * t437 + t392 * t441;
t341 = -t377 + t376;
t340 = t354 * t442 - t437 * t471;
t339 = t352 * t442 - t416 * t438;
t338 = -t354 * t438 - t437 * t464;
t337 = -t352 * t438 - t441 * t464;
t331 = t358 * t442 - t453;
t330 = t357 * t442 + t453;
t329 = -t358 * t438 - t452;
t328 = -t357 * t438 + t452;
t327 = -t427 - t376;
t324 = t355 * t442 + t389 * t438;
t323 = t353 * t442 - t392 * t438;
t322 = t355 * t438 - t389 * t442;
t321 = t353 * t438 + t392 * t442;
t320 = (-t378 * t440 + t380 * t436) * t431;
t319 = (-t378 * t436 - t380 * t440) * t431;
t316 = -t376 - t377;
t315 = t343 * t442 - t399 * t438;
t314 = -t343 * t438 - t399 * t442;
t312 = -qJD(5) * t380 - t455;
t308 = t364 * t440 - t477;
t307 = -t365 * t436 + t494;
t306 = t364 * t436 + t469;
t305 = t365 * t440 + t495;
t304 = -t356 * t436 - t469;
t303 = t356 * t440 - t477;
t302 = t344 * t439 + t347 * t443;
t301 = -t344 * t443 + t347 * t439;
t300 = -t367 + t447;
t295 = (qJD(5) + t431) * t380 + t455;
t294 = -t380 * t482 - t440 * t447;
t293 = t380 * t481 - t436 * t447;
t292 = -t312 * t436 + t378 * t481;
t291 = t312 * t440 + t378 * t482;
t290 = -pkin(7) * t351 + t467;
t289 = -pkin(7) * t350 + t475;
t288 = t327 * t440 - t495;
t287 = t327 * t436 + t494;
t286 = -pkin(3) * t351 + t310;
t285 = -pkin(3) * t350 + t309;
t276 = t322 * t439 + t324 * t443;
t275 = t321 * t439 + t323 * t443;
t274 = -t322 * t443 + t324 * t439;
t273 = -t321 * t443 + t323 * t439;
t272 = pkin(6) * t282 + qJ(2) * g(3);
t269 = -pkin(6) * t283 + g(3) * t485;
t268 = -t319 * t437 + t320 * t441;
t267 = t268 * t442 + t438 * t460;
t266 = -t268 * t438 + t442 * t460;
t262 = -t306 * t437 + t308 * t441;
t261 = -t305 * t437 + t307 * t441;
t260 = -t303 * t437 + t304 * t441;
t259 = t303 * t441 + t304 * t437;
t258 = -t296 * t440 - t300 * t436;
t257 = -t295 * t440 - t436 * t489;
t256 = -t296 * t436 + t300 * t440;
t255 = -t295 * t436 + t440 * t489;
t254 = -pkin(6) * t344 + t264 * t442;
t253 = -pkin(6) * t347 - t264 * t438;
t252 = -pkin(8) * t303 + t470;
t251 = -t293 * t437 + t294 * t441;
t250 = -t291 * t437 + t292 * t441;
t249 = -t287 * t437 + t288 * t441;
t248 = t287 * t441 + t288 * t437;
t247 = -pkin(8) * t287 + t478;
t246 = t265 * t442 + t317 * t438;
t245 = t265 * t438 - t317 * t442;
t239 = t251 * t442 + t459;
t238 = t250 * t442 - t459;
t237 = -t251 * t438 + t458;
t236 = -t250 * t438 - t458;
t235 = t262 * t442 - t296 * t438;
t234 = t261 * t442 - t300 * t438;
t233 = -t262 * t438 - t296 * t442;
t232 = -t261 * t438 - t300 * t442;
t231 = t260 * t442 + t438 * t489;
t230 = t260 * t438 - t442 * t489;
t229 = -pkin(4) * t489 + pkin(8) * t304 + t478;
t228 = -pkin(6) * t322 + qJ(2) * t351 - t286 * t438 + t290 * t442;
t227 = -pkin(6) * t321 + qJ(2) * t350 - t285 * t438 + t289 * t442;
t226 = t249 * t442 + t295 * t438;
t225 = t249 * t438 - t295 * t442;
t224 = -pkin(4) * t295 + pkin(8) * t288 - t470;
t223 = -pkin(6) * t324 - t442 * t286 - t438 * t290 + t351 * t485;
t222 = -pkin(6) * t323 - t442 * t285 - t438 * t289 + t350 * t485;
t221 = -t256 * t437 + t258 * t441;
t220 = -t255 * t437 + t257 * t441;
t219 = t256 * t441 + t258 * t437;
t218 = t220 * t442 - t341 * t438;
t217 = -t220 * t438 - t341 * t442;
t216 = t221 * t442 + t316 * t438;
t215 = t221 * t438 - t316 * t442;
t214 = t245 * t439 + t246 * t443;
t213 = -t245 * t443 + t246 * t439;
t210 = -pkin(3) * t259 - pkin(4) * t303 + t242;
t209 = -pkin(3) * t248 - pkin(4) * t287 + t241;
t208 = t230 * t439 + t231 * t443;
t207 = -t230 * t443 + t231 * t439;
t206 = -pkin(3) * t219 - pkin(4) * t256;
t205 = -pkin(4) * t284 + pkin(8) * t212;
t204 = t225 * t439 + t226 * t443;
t203 = -t225 * t443 + t226 * t439;
t202 = -pkin(8) * t256 - t211;
t201 = -pkin(7) * t259 - t229 * t437 + t252 * t441;
t200 = -pkin(6) * t245 - (pkin(3) * t438 - pkin(7) * t442 + qJ(2)) * t264;
t199 = -pkin(4) * t316 + pkin(8) * t258 + t212;
t198 = -pkin(7) * t248 - t224 * t437 + t247 * t441;
t197 = -pkin(6) * t246 - (pkin(3) * t442 + pkin(7) * t438 + t485) * t264;
t196 = t215 * t439 + t216 * t443;
t195 = -t215 * t443 + t216 * t439;
t194 = t212 * t441 - t476;
t193 = t212 * t437 + t468;
t192 = t194 * t442 + t284 * t438;
t191 = t194 * t438 - t284 * t442;
t190 = -pkin(3) * t193 - pkin(4) * t211;
t189 = -pkin(6) * t230 + qJ(2) * t259 + t201 * t442 - t210 * t438;
t188 = -pkin(7) * t219 - t199 * t437 + t202 * t441;
t187 = -pkin(6) * t231 - t438 * t201 - t442 * t210 + t259 * t485;
t186 = -pkin(6) * t225 + qJ(2) * t248 + t198 * t442 - t209 * t438;
t185 = -pkin(6) * t226 - t438 * t198 - t442 * t209 + t248 * t485;
t184 = -pkin(7) * t193 - pkin(8) * t468 - t205 * t437;
t183 = t191 * t439 + t192 * t443;
t182 = -t191 * t443 + t192 * t439;
t181 = -pkin(6) * t215 + qJ(2) * t219 + t188 * t442 - t206 * t438;
t180 = -pkin(6) * t216 - t438 * t188 - t442 * t206 + t219 * t485;
t179 = -pkin(6) * t191 + qJ(2) * t193 + t184 * t442 - t190 * t438;
t178 = -pkin(6) * t192 - t438 * t184 - t442 * t190 + t193 * t485;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t404, -t405, 0, t362, 0, 0, 0, 0, 0, 0, -t404, 0, t405, t336, 0, 0, 0, 0, 0, 0, -t345, t487, 0, -t491, 0, 0, 0, 0, 0, 0, t275, t276, t302, t214, 0, 0, 0, 0, 0, 0, t204, t208, t196, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t405, -t404, 0, t361, 0, 0, 0, 0, 0, 0, t405, 0, t404, t335, 0, 0, 0, 0, 0, 0, t487, t345, 0, t243, 0, 0, 0, 0, 0, 0, t273, t274, t301, t213, 0, 0, 0, 0, 0, 0, t203, t207, t195, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t350, -t351, 0, t264, 0, 0, 0, 0, 0, 0, -t248, -t259, -t219, -t193; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t405, 0, -t404, 0, -t387, -t388, -t361, -pkin(5) * t361, 0, t405, 0, 0, t404, 0, -t387, -t335, t388, -pkin(5) * t335 + (-pkin(1) * t439 + qJ(2) * t443) * g(3), 0, 0, -t487, 0, -t345, 0, t497, t496, t243, -pkin(5) * t243 - t269 * t439 + t272 * t443, -t329 * t439 + t331 * t443, -t314 * t439 + t315 * t443, -t338 * t439 + t340 * t443, -t328 * t439 + t330 * t443, -t337 * t439 + t339 * t443, -t359 * t439 + t360 * t443, -pkin(5) * t273 - t222 * t439 + t227 * t443, -pkin(5) * t274 - t223 * t439 + t228 * t443, -pkin(5) * t301 - t253 * t439 + t254 * t443, -pkin(5) * t213 - t197 * t439 + t200 * t443, -t237 * t439 + t239 * t443, -t217 * t439 + t218 * t443, -t232 * t439 + t234 * t443, -t236 * t439 + t238 * t443, -t233 * t439 + t235 * t443, -t266 * t439 + t267 * t443, -pkin(5) * t203 - t185 * t439 + t186 * t443, -pkin(5) * t207 - t187 * t439 + t189 * t443, -pkin(5) * t195 - t180 * t439 + t181 * t443, -pkin(5) * t182 - t178 * t439 + t179 * t443; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t404, 0, t405, 0, t388, -t387, t362, pkin(5) * t362, 0, t404, 0, 0, -t405, 0, t388, t336, t387, pkin(5) * t336 + (pkin(1) * t443 + qJ(2) * t439) * g(3), 0, 0, -t345, 0, t487, 0, t496, -t497, t491, -pkin(5) * t491 + t269 * t443 + t272 * t439, t329 * t443 + t331 * t439, t314 * t443 + t315 * t439, t338 * t443 + t340 * t439, t328 * t443 + t330 * t439, t337 * t443 + t339 * t439, t359 * t443 + t360 * t439, pkin(5) * t275 + t222 * t443 + t227 * t439, pkin(5) * t276 + t223 * t443 + t228 * t439, pkin(5) * t302 + t253 * t443 + t254 * t439, pkin(5) * t214 + t197 * t443 + t200 * t439, t237 * t443 + t239 * t439, t217 * t443 + t218 * t439, t232 * t443 + t234 * t439, t236 * t443 + t238 * t439, t233 * t443 + t235 * t439, t266 * t443 + t267 * t439, pkin(5) * t204 + t185 * t443 + t186 * t439, pkin(5) * t208 + t187 * t443 + t189 * t439, pkin(5) * t196 + t180 * t443 + t181 * t439, pkin(5) * t183 + t178 * t443 + t179 * t439; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t412, t413, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t449 + (2 * t483), 0, 0.2e1 * t433 + t450, pkin(1) * t385 + qJ(2) * t382, 0, 0, 0, 0, 0, t430, qJ(2) * t454 + t451 * t485 + t325, qJ(2) * t451 - t454 * t485 + t326, 0, qJ(2) * t283 + t282 * t485, (-t390 - t456) * t437, -t389 * t441 - t392 * t437, -t406 * t441 - t474, (-t391 + t457) * t441, -t408 * t437 - t465, 0, -pkin(3) * t392 - pkin(7) * t353 + qJ(2) * t323 - t321 * t485 + t467, pkin(3) * t389 - pkin(7) * t355 + qJ(2) * t324 - t322 * t485 - t475, -pkin(3) * t398 + pkin(7) * t393 + qJ(2) * t347 - t344 * t485 - t265, pkin(3) * t317 - pkin(7) * t265 + qJ(2) * t246 - t245 * t485, -t293 * t441 - t294 * t437, -t255 * t441 - t257 * t437, -t305 * t441 - t307 * t437, -t291 * t441 - t292 * t437, -t306 * t441 - t308 * t437, -t319 * t441 - t320 * t437, pkin(3) * t295 - pkin(7) * t249 + qJ(2) * t226 - t441 * t224 - t225 * t485 - t437 * t247, pkin(3) * t489 - pkin(7) * t260 + qJ(2) * t231 - t441 * t229 - t230 * t485 - t437 * t252, pkin(3) * t316 - pkin(7) * t221 + qJ(2) * t216 - t441 * t199 - t437 * t202 - t215 * t485, pkin(3) * t284 - pkin(7) * t194 + pkin(8) * t476 + qJ(2) * t192 - t191 * t485 - t441 * t205;];
tauB_reg = t1;