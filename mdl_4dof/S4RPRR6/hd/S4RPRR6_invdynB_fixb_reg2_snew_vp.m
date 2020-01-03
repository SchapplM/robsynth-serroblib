% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRR6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:46
% EndTime: 2019-12-31 16:52:51
% DurationCPUTime: 4.99s
% Computational Cost: add. (18614->400), mult. (45878->613), div. (0->0), fcn. (33364->8), ass. (0->282)
t447 = sin(qJ(3));
t445 = cos(pkin(7));
t450 = cos(qJ(3));
t444 = sin(pkin(7));
t491 = t444 * t447;
t416 = (t445 * t450 - t491) * qJD(1);
t459 = t444 * t450 + t445 * t447;
t418 = t459 * qJD(1);
t497 = t418 * t416;
t506 = qJDD(3) + t497;
t511 = t447 * t506;
t510 = t450 * t506;
t448 = sin(qJ(1));
t451 = cos(qJ(1));
t430 = g(1) * t451 + g(2) * t448;
t453 = qJD(1) ^ 2;
t507 = -pkin(1) * t453 + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t430;
t446 = sin(qJ(4));
t449 = cos(qJ(4));
t374 = -t416 * t449 + t418 * t446;
t376 = t416 * t446 + t418 * t449;
t338 = t376 * t374;
t470 = qJDD(3) + qJDD(4);
t505 = -t338 + t470;
t509 = t446 * t505;
t508 = t449 * t505;
t443 = qJD(3) + qJD(4);
t368 = t443 * t374;
t473 = qJDD(1) * t445;
t414 = qJDD(1) * t491 - t450 * t473;
t476 = t418 * qJD(3);
t385 = -t414 - t476;
t477 = t416 * qJD(3);
t502 = t459 * qJDD(1);
t387 = t502 + t477;
t456 = qJD(4) * t374 - t385 * t446 - t387 * t449;
t504 = -t368 - t456;
t441 = t444 ^ 2;
t442 = t445 ^ 2;
t503 = t441 + t442;
t501 = t453 * t503;
t463 = -t385 * t449 + t446 * t387;
t297 = (qJD(4) - t443) * t376 + t463;
t372 = t374 ^ 2;
t373 = t376 ^ 2;
t412 = t416 ^ 2;
t413 = t418 ^ 2;
t439 = t443 ^ 2;
t500 = pkin(2) * t445;
t499 = t445 * g(3);
t498 = qJDD(1) * pkin(1);
t496 = t441 * t453;
t436 = t442 * t453;
t495 = t443 * t446;
t494 = t443 * t449;
t364 = -t499 + (-pkin(5) * qJDD(1) + t453 * t500 - t507) * t444;
t392 = -t444 * g(3) + t445 * t507;
t369 = -pkin(2) * t436 + pkin(5) * t473 + t392;
t326 = -t364 * t450 + t447 * t369;
t327 = t364 * t447 + t369 * t450;
t277 = -t326 * t450 + t327 * t447;
t493 = t444 * t277;
t492 = t444 * t445;
t490 = t445 * t277;
t429 = t448 * g(1) - g(2) * t451;
t461 = -qJDD(2) + t429;
t380 = (pkin(1) + t500) * qJDD(1) + (pkin(5) * t503 + qJ(2)) * t453 + t461;
t460 = qJD(3) * pkin(3) - pkin(6) * t418;
t317 = t385 * pkin(3) + t412 * pkin(6) - t418 * t460 + t380;
t489 = t446 * t317;
t335 = t338 + t470;
t488 = t446 * t335;
t283 = (-t387 + t477) * pkin(6) + t506 * pkin(3) - t326;
t289 = -t412 * pkin(3) + t385 * pkin(6) - qJD(3) * t460 + t327;
t247 = -t283 * t449 + t289 * t446;
t248 = t283 * t446 + t289 * t449;
t222 = -t247 * t449 + t248 * t446;
t487 = t447 * t222;
t486 = t447 * t380;
t382 = qJDD(3) - t497;
t485 = t447 * t382;
t410 = t453 * qJ(2) + t461 + t498;
t484 = t448 * t410;
t483 = t449 * t317;
t482 = t449 * t335;
t481 = t450 * t222;
t480 = t450 * t380;
t479 = t450 * t382;
t478 = t451 * t410;
t472 = t448 * qJDD(1);
t471 = t451 * qJDD(1);
t468 = t448 * t338;
t467 = t448 * t497;
t466 = t451 * t338;
t465 = t451 * t497;
t464 = t410 + t498;
t223 = t247 * t446 + t248 * t449;
t278 = t326 * t447 + t327 * t450;
t391 = t444 * t507 + t499;
t347 = t391 * t444 + t392 * t445;
t400 = -t429 * t448 - t430 * t451;
t428 = -t448 * t453 + t471;
t462 = -pkin(4) * t428 - g(3) * t448;
t346 = t391 * t445 - t392 * t444;
t399 = t429 * t451 - t430 * t448;
t427 = t451 * t453 + t472;
t421 = t445 * t501;
t396 = -t421 * t448 + t445 * t471;
t458 = t421 * t451 + t445 * t472;
t452 = qJD(3) ^ 2;
t435 = t442 * qJDD(1);
t434 = t441 * qJDD(1);
t426 = t436 - t496;
t425 = t436 + t496;
t424 = t435 - t434;
t423 = t435 + t434;
t420 = t444 * t501;
t411 = -pkin(4) * t427 + g(3) * t451;
t405 = -t413 - t452;
t404 = -t413 + t452;
t403 = t412 - t452;
t402 = t428 * t492;
t401 = t427 * t492;
t397 = t420 * t451 + t444 * t472;
t395 = t420 * t448 - t444 * t471;
t394 = t423 * t451 - t425 * t448;
t393 = t423 * t448 + t425 * t451;
t389 = -t413 + t412;
t386 = t502 + 0.2e1 * t477;
t384 = t414 + 0.2e1 * t476;
t379 = -t452 - t412;
t371 = (t416 * t450 + t418 * t447) * qJD(3);
t370 = (t416 * t447 - t418 * t450) * qJD(3);
t366 = -t373 + t439;
t365 = t372 - t439;
t362 = -t373 - t439;
t359 = -t412 - t413;
t357 = t387 * t450 - t447 * t476;
t356 = t387 * t447 + t450 * t476;
t355 = -t385 * t447 - t450 * t477;
t354 = t385 * t450 - t447 * t477;
t353 = -t405 * t447 - t479;
t352 = -t404 * t447 + t510;
t351 = t403 * t450 - t485;
t350 = t405 * t450 - t485;
t349 = t404 * t450 + t511;
t348 = t403 * t447 + t479;
t344 = -t384 * t450 - t386 * t447;
t343 = -t414 * t450 + t447 * t502;
t342 = -t384 * t447 + t386 * t450;
t341 = -t414 * t447 - t450 * t502;
t340 = t379 * t450 - t511;
t339 = t379 * t447 + t510;
t337 = -t373 + t372;
t333 = -t439 - t372;
t332 = -t370 * t444 + t371 * t445;
t331 = (-t374 * t449 + t376 * t446) * t443;
t330 = (-t374 * t446 - t376 * t449) * t443;
t329 = t347 * t451 - t484;
t328 = t347 * t448 + t478;
t324 = -pkin(5) * t350 - t480;
t322 = -qJD(4) * t376 - t463;
t321 = -pkin(5) * t339 - t486;
t320 = -t372 - t373;
t319 = -t356 * t444 + t357 * t445;
t318 = -t354 * t444 + t355 * t445;
t316 = -t350 * t444 + t353 * t445;
t315 = -t349 * t444 + t352 * t445;
t314 = -t348 * t444 + t351 * t445;
t313 = t350 * t445 + t353 * t444;
t312 = t365 * t449 - t488;
t311 = -t366 * t446 + t508;
t310 = t365 * t446 + t482;
t309 = t366 * t449 + t509;
t308 = -pkin(2) * t386 + pkin(5) * t353 - t486;
t307 = -t362 * t446 - t482;
t306 = t362 * t449 - t488;
t305 = -pkin(2) * t384 + pkin(5) * t340 + t480;
t304 = -t342 * t444 + t344 * t445;
t303 = -t341 * t444 + t343 * t445;
t302 = t341 * t445 + t343 * t444;
t301 = -t368 + t456;
t296 = (qJD(4) + t443) * t376 + t463;
t295 = -t339 * t444 + t340 * t445;
t294 = t339 * t445 + t340 * t444;
t293 = -t376 * t495 - t449 * t456;
t292 = t376 * t494 - t446 * t456;
t291 = -t322 * t446 + t374 * t494;
t290 = t322 * t449 + t374 * t495;
t288 = t316 * t451 + t386 * t448;
t287 = t316 * t448 - t386 * t451;
t286 = t333 * t449 - t509;
t285 = t333 * t446 + t508;
t280 = -t330 * t447 + t331 * t450;
t279 = t330 * t450 + t331 * t447;
t276 = t295 * t451 + t384 * t448;
t275 = t295 * t448 - t384 * t451;
t274 = t303 * t451 + t359 * t448;
t273 = t303 * t448 - t359 * t451;
t272 = pkin(2) * t380 + pkin(5) * t278;
t271 = -pkin(1) * t302 - pkin(2) * t341;
t270 = -pkin(6) * t306 - t483;
t269 = -pkin(5) * t341 - t277;
t268 = -t310 * t447 + t312 * t450;
t267 = -t309 * t447 + t311 * t450;
t266 = t310 * t450 + t312 * t447;
t265 = t309 * t450 + t311 * t447;
t264 = -t306 * t447 + t307 * t450;
t263 = t306 * t450 + t307 * t447;
t262 = -pkin(1) * t313 - pkin(2) * t350 + t327;
t261 = -pkin(6) * t285 - t489;
t260 = -pkin(2) * t359 + pkin(5) * t343 + t278;
t259 = -t297 * t449 - t301 * t446;
t258 = -t296 * t449 - t446 * t504;
t257 = -t297 * t446 + t301 * t449;
t256 = -t296 * t446 + t449 * t504;
t255 = -pkin(1) * t294 - pkin(2) * t339 + t326;
t254 = -t292 * t447 + t293 * t450;
t253 = -t290 * t447 + t291 * t450;
t252 = t292 * t450 + t293 * t447;
t251 = t290 * t450 + t291 * t447;
t250 = -t285 * t447 + t286 * t450;
t249 = t285 * t450 + t286 * t447;
t245 = -qJ(2) * t313 - t308 * t444 + t324 * t445;
t244 = -t279 * t444 + t280 * t445;
t243 = t278 * t445 - t493;
t242 = t278 * t444 + t490;
t241 = -pkin(3) * t504 + pkin(6) * t307 - t489;
t240 = -qJ(2) * t294 - t305 * t444 + t321 * t445;
t239 = t243 * t451 - t380 * t448;
t238 = t243 * t448 + t380 * t451;
t237 = -pkin(3) * t296 + pkin(6) * t286 + t483;
t236 = -t266 * t444 + t268 * t445;
t235 = -t265 * t444 + t267 * t445;
t234 = -t263 * t444 + t264 * t445;
t233 = t263 * t445 + t264 * t444;
t232 = -pkin(1) * t242 - pkin(2) * t277;
t231 = -t257 * t447 + t259 * t450;
t230 = -t256 * t447 + t258 * t450;
t229 = t257 * t450 + t259 * t447;
t228 = t256 * t450 + t258 * t447;
t227 = -t252 * t444 + t254 * t445;
t226 = -t251 * t444 + t253 * t445;
t225 = -t249 * t444 + t250 * t445;
t224 = t249 * t445 + t250 * t444;
t221 = t234 * t451 + t448 * t504;
t220 = t234 * t448 - t451 * t504;
t219 = -qJ(2) * t302 - t260 * t444 + t269 * t445;
t218 = -pkin(5) * t490 - qJ(2) * t242 - t272 * t444;
t217 = pkin(3) * t317 + pkin(6) * t223;
t216 = t225 * t451 + t296 * t448;
t215 = t225 * t448 - t296 * t451;
t214 = -pkin(5) * t263 - t241 * t447 + t270 * t450;
t213 = -pkin(6) * t257 - t222;
t212 = -pkin(5) * t249 - t237 * t447 + t261 * t450;
t211 = -pkin(2) * t504 + pkin(5) * t264 + t241 * t450 + t270 * t447;
t210 = -pkin(3) * t320 + pkin(6) * t259 + t223;
t209 = -pkin(2) * t296 + pkin(5) * t250 + t237 * t450 + t261 * t447;
t208 = -t229 * t444 + t231 * t445;
t207 = -t228 * t444 + t230 * t445;
t206 = t229 * t445 + t231 * t444;
t205 = -pkin(1) * t233 - pkin(2) * t263 - pkin(3) * t306 + t248;
t204 = t208 * t451 + t320 * t448;
t203 = t208 * t448 - t320 * t451;
t202 = t223 * t450 - t487;
t201 = t223 * t447 + t481;
t200 = -pkin(1) * t224 - pkin(2) * t249 - pkin(3) * t285 + t247;
t199 = -pkin(1) * t206 - pkin(2) * t229 - pkin(3) * t257;
t198 = -qJ(2) * t233 - t211 * t444 + t214 * t445;
t197 = -pkin(5) * t229 - t210 * t447 + t213 * t450;
t196 = -t201 * t444 + t202 * t445;
t195 = t201 * t445 + t202 * t444;
t194 = -pkin(2) * t320 + pkin(5) * t231 + t210 * t450 + t213 * t447;
t193 = -pkin(5) * t201 - pkin(6) * t481 - t217 * t447;
t192 = t196 * t451 - t317 * t448;
t191 = t196 * t448 + t317 * t451;
t190 = -qJ(2) * t224 - t209 * t444 + t212 * t445;
t189 = pkin(2) * t317 + pkin(5) * t202 - pkin(6) * t487 + t217 * t450;
t188 = -pkin(1) * t195 - pkin(2) * t201 - pkin(3) * t222;
t187 = -qJ(2) * t206 - t194 * t444 + t197 * t445;
t186 = -qJ(2) * t195 - t189 * t444 + t193 * t445;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t427, -t428, 0, t400, 0, 0, 0, 0, 0, 0, -t458, t397, t394, t329, 0, 0, 0, 0, 0, 0, t276, t288, t274, t239, 0, 0, 0, 0, 0, 0, t216, t221, t204, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t428, -t427, 0, t399, 0, 0, 0, 0, 0, 0, t396, t395, t393, t328, 0, 0, 0, 0, 0, 0, t275, t287, t273, t238, 0, 0, 0, 0, 0, 0, t215, t220, t203, t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t346, 0, 0, 0, 0, 0, 0, t294, t313, t302, t242, 0, 0, 0, 0, 0, 0, t224, t233, t206, t195; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t428, 0, -t427, 0, t462, -t411, -t399, -pkin(4) * t399, t402, t424 * t451 - t426 * t448, t397, -t402, t458, 0, -pkin(4) * t396 - t391 * t448 - t444 * t478, -pkin(4) * t395 - t392 * t448 - t445 * t478, -pkin(4) * t393 + t346 * t451, -pkin(4) * t328 - (pkin(1) * t448 - qJ(2) * t451) * t346, t319 * t451 - t467, t304 * t451 - t389 * t448, t315 * t451 + t448 * t502, t318 * t451 + t467, t314 * t451 - t414 * t448, qJDD(3) * t448 + t332 * t451, -pkin(4) * t275 + t240 * t451 - t255 * t448, -pkin(4) * t287 + t245 * t451 - t262 * t448, -pkin(4) * t273 + t219 * t451 - t271 * t448, -pkin(4) * t238 + t218 * t451 - t232 * t448, t227 * t451 + t468, t207 * t451 - t337 * t448, t235 * t451 - t301 * t448, t226 * t451 - t468, t236 * t451 - t297 * t448, t244 * t451 + t448 * t470, -pkin(4) * t215 + t190 * t451 - t200 * t448, -pkin(4) * t220 + t198 * t451 - t205 * t448, -pkin(4) * t203 + t187 * t451 - t199 * t448, -pkin(4) * t191 + t186 * t451 - t188 * t448; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t427, 0, t428, 0, t411, t462, t400, pkin(4) * t400, t401, t424 * t448 + t426 * t451, t395, -t401, -t396, 0, -pkin(4) * t458 + t391 * t451 - t444 * t484, pkin(4) * t397 + t392 * t451 - t445 * t484, pkin(4) * t394 + t346 * t448, pkin(4) * t329 - (-pkin(1) * t451 - qJ(2) * t448) * t346, t319 * t448 + t465, t304 * t448 + t389 * t451, t315 * t448 - t451 * t502, t318 * t448 - t465, t314 * t448 + t414 * t451, -qJDD(3) * t451 + t332 * t448, pkin(4) * t276 + t240 * t448 + t255 * t451, pkin(4) * t288 + t245 * t448 + t262 * t451, pkin(4) * t274 + t219 * t448 + t271 * t451, pkin(4) * t239 + t218 * t448 + t232 * t451, t227 * t448 - t466, t207 * t448 + t337 * t451, t235 * t448 + t301 * t451, t226 * t448 + t466, t236 * t448 + t297 * t451, t244 * t448 - t451 * t470, pkin(4) * t216 + t190 * t448 + t200 * t451, pkin(4) * t221 + t198 * t448 + t205 * t451, pkin(4) * t204 + t187 * t448 + t199 * t451, pkin(4) * t192 + t186 * t448 + t188 * t451; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t429, t430, 0, 0, t434, 0.2e1 * t444 * t473, 0, t435, 0, 0, -qJ(2) * t421 + t445 * t464, qJ(2) * t420 - t444 * t464, pkin(1) * t425 + qJ(2) * t423 + t347, pkin(1) * t410 + qJ(2) * t347, t356 * t445 + t357 * t444, t342 * t445 + t344 * t444, t349 * t445 + t352 * t444, t354 * t445 + t355 * t444, t348 * t445 + t351 * t444, t370 * t445 + t371 * t444, -pkin(1) * t384 + qJ(2) * t295 + t305 * t445 + t321 * t444, -pkin(1) * t386 + qJ(2) * t316 + t308 * t445 + t324 * t444, -pkin(1) * t359 + qJ(2) * t303 + t260 * t445 + t269 * t444, pkin(1) * t380 - pkin(5) * t493 + qJ(2) * t243 + t272 * t445, t252 * t445 + t254 * t444, t228 * t445 + t230 * t444, t265 * t445 + t267 * t444, t251 * t445 + t253 * t444, t266 * t445 + t268 * t444, t279 * t445 + t280 * t444, -pkin(1) * t296 + qJ(2) * t225 + t209 * t445 + t212 * t444, -pkin(1) * t504 + qJ(2) * t234 + t211 * t445 + t214 * t444, -pkin(1) * t320 + qJ(2) * t208 + t194 * t445 + t197 * t444, pkin(1) * t317 + qJ(2) * t196 + t189 * t445 + t193 * t444;];
tauB_reg = t1;