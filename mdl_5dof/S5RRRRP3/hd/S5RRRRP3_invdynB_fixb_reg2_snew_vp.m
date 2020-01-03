% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRRP3
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
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRRP3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:29
% EndTime: 2019-12-31 21:49:39
% DurationCPUTime: 7.03s
% Computational Cost: add. (21071->393), mult. (25393->520), div. (0->0), fcn. (14266->8), ass. (0->283)
t468 = qJD(4) ^ 2;
t455 = qJD(1) + qJD(2);
t449 = qJD(3) + t455;
t447 = t449 ^ 2;
t459 = sin(qJ(4));
t457 = t459 ^ 2;
t527 = t457 * t447;
t432 = t468 + t527;
t463 = cos(qJ(4));
t437 = t463 * t447 * t459;
t430 = qJDD(4) - t437;
t511 = t463 * t430;
t384 = -t459 * t432 + t511;
t501 = qJD(4) * t449;
t495 = t463 * t501;
t499 = qJDD(1) + qJDD(2);
t448 = qJDD(3) + t499;
t519 = t459 * t448;
t408 = 0.2e1 * t495 + t519;
t460 = sin(qJ(3));
t464 = cos(qJ(3));
t333 = t460 * t384 + t464 * t408;
t336 = t464 * t384 - t460 * t408;
t461 = sin(qJ(2));
t465 = cos(qJ(2));
t282 = t465 * t333 + t461 * t336;
t286 = t461 * t333 - t465 * t336;
t462 = sin(qJ(1));
t466 = cos(qJ(1));
t234 = t466 * t282 - t462 * t286;
t586 = pkin(5) * t234;
t238 = t462 * t282 + t466 * t286;
t585 = pkin(5) * t238;
t584 = pkin(6) * t282;
t583 = pkin(1) * t282 + pkin(2) * t333 + pkin(8) * t384;
t522 = t459 * t430;
t377 = t463 * t432 + t522;
t582 = -pkin(1) * t377 - pkin(6) * t286;
t368 = t459 * t408;
t496 = t459 * t501;
t510 = t463 * t448;
t409 = -0.2e1 * t496 + t510;
t369 = t463 * t409;
t349 = -t369 + t368;
t458 = t463 ^ 2;
t418 = (t457 - t458) * t447;
t325 = t460 * t349 + t464 * t418;
t327 = t464 * t349 - t460 * t418;
t268 = t465 * t325 + t461 * t327;
t269 = t461 * t325 - t465 * t327;
t581 = t466 * t268 - t462 * t269;
t580 = t462 * t268 + t466 * t269;
t526 = t458 * t447;
t434 = -t468 + t526;
t382 = -t463 * t434 + t522;
t507 = t464 * t448;
t343 = t460 * t382 + t463 * t507;
t346 = t464 * t382 - t460 * t510;
t295 = t465 * t343 + t461 * t346;
t297 = t461 * t343 - t465 * t346;
t579 = t466 * t295 - t462 * t297;
t578 = t462 * t295 + t466 * t297;
t516 = t460 * t448;
t413 = t464 * t447 + t516;
t416 = t460 * t447 - t507;
t353 = t465 * t413 - t461 * t416;
t395 = pkin(7) * t413 - t464 * g(3);
t551 = pkin(7) * t416 - t460 * g(3);
t303 = pkin(6) * t353 + t465 * t395 - t461 * t551;
t358 = t461 * t413 + t465 * t416;
t314 = t462 * t353 + t466 * t358;
t563 = pkin(6) * t358 + t461 * t395 + t465 * t551;
t577 = pkin(5) * t314 + t462 * t303 + t466 * t563;
t550 = t466 * t353 - t462 * t358;
t576 = pkin(5) * t550 + t466 * t303 - t462 * t563;
t575 = pkin(7) * t333;
t571 = -pkin(2) * t377 + pkin(7) * t336;
t444 = t466 * g(1) + t462 * g(2);
t541 = qJD(1) ^ 2;
t473 = -t541 * pkin(1) - t444;
t443 = t462 * g(1) - t466 * g(2);
t478 = qJDD(1) * pkin(1) + t443;
t375 = t461 * t478 + t465 * t473;
t454 = t455 ^ 2;
t363 = -t454 * pkin(2) + t375;
t470 = -t461 * t473 + t465 * t478;
t469 = t499 * pkin(2) + t470;
t317 = t460 * t363 - t464 * t469;
t318 = t464 * t363 + t460 * t469;
t491 = t460 * t317 + t464 * t318;
t260 = t464 * t317 - t460 * t318;
t506 = t465 * t260;
t221 = -t461 * t491 + t506;
t515 = t461 * t260;
t558 = t465 * t491 + t515;
t201 = t462 * t221 + t466 * t558;
t570 = t466 * t221 - t462 * t558;
t425 = t465 * t454 + t461 * t499;
t427 = -t461 * t454 + t465 * t499;
t366 = t462 * t425 - t466 * t427;
t399 = pkin(6) * t425 - t465 * g(3);
t557 = -pkin(6) * t427 - t461 * g(3);
t569 = pkin(5) * t366 + t462 * t399 + t466 * t557;
t480 = t466 * t425 + t462 * t427;
t568 = pkin(5) * t480 + t466 * t399 - t462 * t557;
t490 = t465 * t375 - t461 * t470;
t322 = -t461 * t375 - t465 * t470;
t504 = t466 * t322;
t559 = -t462 * t490 + t504;
t513 = t462 * t322;
t266 = t466 * t490 + t513;
t556 = 2 * qJD(5);
t553 = pkin(3) * t377;
t552 = pkin(8) * t377;
t540 = pkin(4) * t463;
t483 = -qJ(5) * t459 - t540;
t407 = t483 * t449;
t307 = -t447 * pkin(3) + t448 * pkin(8) + t318;
t492 = t459 * g(3) - t463 * t307;
t475 = t463 * t449 * t407 + qJDD(4) * qJ(5) + (qJD(4) * t556) - t492;
t543 = t459 * t434 + t511;
t528 = t449 * t459;
t542 = t407 * t528 + qJDD(5);
t435 = -t468 - t526;
t429 = qJDD(4) + t437;
t523 = t459 * t429;
t381 = t463 * t435 - t523;
t332 = t460 * t381 + t464 * t409;
t335 = t464 * t381 - t460 * t409;
t281 = t465 * t332 + t461 * t335;
t284 = -t461 * t332 + t465 * t335;
t233 = t466 * t281 + t462 * t284;
t539 = pkin(5) * t233;
t502 = t457 + t458;
t411 = t502 * t448;
t417 = t502 * t447;
t352 = t460 * t411 + t464 * t417;
t356 = t464 * t411 - t460 * t417;
t310 = t465 * t352 + t461 * t356;
t313 = -t461 * t352 + t465 * t356;
t253 = t466 * t310 + t462 * t313;
t538 = pkin(5) * t253;
t537 = pkin(6) * t281;
t536 = pkin(6) * t310;
t535 = pkin(7) * t332;
t534 = pkin(7) * t352;
t420 = t463 * t429;
t376 = t459 * t435 + t420;
t533 = pkin(8) * t376;
t306 = -t448 * pkin(3) - t447 * pkin(8) + t317;
t525 = t459 * t306;
t524 = t459 * t409;
t512 = t463 * t306;
t292 = t463 * g(3) + t459 * t307;
t503 = t417 - t468;
t494 = -pkin(1) * t376 + pkin(6) * t284;
t493 = -pkin(2) * t376 + pkin(7) * t335;
t247 = t459 * t292 - t463 * t492;
t391 = -t462 * t443 - t466 * t444;
t488 = t460 * t437;
t487 = t464 * t437;
t486 = pkin(1) * t281 + pkin(2) * t332 + pkin(3) * t409 + pkin(8) * t381;
t485 = pkin(1) * t310 + pkin(2) * t352 + pkin(3) * t417 + pkin(8) * t411;
t275 = -pkin(3) * t376 + t292;
t439 = t466 * qJDD(1) - t462 * t541;
t484 = -pkin(5) * t439 - t462 * g(3);
t482 = pkin(4) * t459 - qJ(5) * t463;
t246 = t463 * t292 + t459 * t492;
t481 = t463 * t408 + t524;
t390 = t466 * t443 - t462 * t444;
t477 = t495 + t519;
t476 = -t496 + t510;
t474 = -qJDD(4) * pkin(4) + t292 + t542;
t472 = -t476 * pkin(4) + t306 + (-t477 - t495) * qJ(5);
t471 = t528 * t556 - t472;
t467 = pkin(1) * g(3);
t438 = t462 * qJDD(1) + t466 * t541;
t433 = t468 - t527;
t419 = -pkin(5) * t438 + t466 * g(3);
t405 = t502 * t501;
t404 = t482 * t448;
t387 = t460 * qJDD(4) + t464 * t405;
t386 = -t464 * qJDD(4) + t460 * t405;
t383 = -t459 * t433 + t420;
t379 = t463 * t433 + t523;
t373 = -t457 * t501 + t463 * t477;
t372 = -t458 * t501 - t459 * t476;
t350 = pkin(7) * t356;
t347 = t464 * t383 + t459 * t516;
t344 = t460 * t383 - t459 * t507;
t341 = t464 * t373 - t488;
t340 = t464 * t372 + t488;
t339 = t460 * t373 + t487;
t338 = t460 * t372 - t487;
t329 = -t461 * t386 + t465 * t387;
t328 = t465 * t386 + t461 * t387;
t319 = pkin(6) * t490 + t467;
t308 = pkin(6) * t313;
t299 = -t461 * t344 + t465 * t347;
t296 = t465 * t344 + t461 * t347;
t290 = -t461 * t339 + t465 * t341;
t289 = -t461 * t338 + t465 * t340;
t288 = t465 * t339 + t461 * t341;
t287 = t465 * t338 + t461 * t340;
t278 = t512 + t552;
t277 = t525 - t533;
t276 = -t492 + t553;
t274 = -t462 * t328 + t466 * t329;
t273 = t466 * t328 + t462 * t329;
t272 = t468 * qJ(5) - t474;
t271 = -t468 * pkin(4) + t475;
t264 = t503 * qJ(5) + t474;
t263 = t503 * pkin(4) + t475;
t262 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t528 + t472;
t257 = pkin(2) * g(3) + pkin(7) * t491;
t256 = (t409 - t496) * pkin(4) + t471;
t255 = -pkin(4) * t496 + qJ(5) * t408 + t471;
t254 = -t462 * t310 + t466 * t313;
t252 = pkin(5) * t254;
t251 = (-t435 - t468) * qJ(5) + (-qJDD(4) - t429) * pkin(4) + t275 + t542;
t250 = -t553 - qJ(5) * t430 + (-t432 + t468) * pkin(4) - t475;
t249 = -t462 * t296 + t466 * t299;
t248 = t466 * t296 + t462 * t299;
t244 = -t462 * t288 + t466 * t290;
t243 = -t462 * t287 + t466 * t289;
t242 = t466 * t288 + t462 * t290;
t241 = t466 * t287 + t462 * t289;
t240 = -pkin(4) * t368 + t463 * t255 - t552;
t239 = qJ(5) * t369 - t459 * t256 - t533;
t236 = -t462 * t281 + t466 * t284;
t232 = pkin(5) * t236;
t231 = t464 * t246 - t534;
t230 = t460 * t246 + t350;
t229 = t463 * t271 - t459 * t272;
t228 = t459 * t271 + t463 * t272;
t227 = t464 * t247 + t460 * t306;
t226 = t460 * t247 - t464 * t306;
t225 = -t459 * t263 + t463 * t264;
t224 = -t460 * t276 + t464 * t278 + t575;
t223 = -t460 * t275 + t464 * t277 - t535;
t218 = t464 * t276 + t460 * t278 - t571;
t217 = t464 * t275 + t460 * t277 + t493;
t216 = t464 * t225 - t460 * t404 - t534;
t215 = t460 * t225 + t464 * t404 + t350;
t214 = t464 * t229 + t460 * t262;
t213 = t460 * t229 - t464 * t262;
t212 = t464 * t239 - t460 * t251 - t535;
t211 = t464 * t240 - t460 * t250 - t575;
t210 = t460 * t239 + t464 * t251 + t493;
t209 = t460 * t240 + t464 * t250 + t571;
t208 = -pkin(3) * t228 - pkin(4) * t272 - qJ(5) * t271;
t207 = -t461 * t230 + t465 * t231 - t536;
t206 = t465 * t230 + t461 * t231 + t308;
t205 = -pkin(8) * t228 + t482 * t262;
t204 = -t461 * t226 + t465 * t227;
t203 = t465 * t226 + t461 * t227;
t202 = pkin(6) * t221 + pkin(7) * t506 - t461 * t257;
t199 = pkin(6) * t558 + pkin(7) * t515 + t465 * t257 + t467;
t198 = -pkin(7) * t226 - (pkin(3) * t460 - pkin(8) * t464) * t246;
t197 = -t461 * t218 + t465 * t224 + t584;
t196 = -t461 * t217 + t465 * t223 - t537;
t195 = t465 * t218 + t461 * t224 - t582;
t194 = t465 * t217 + t461 * t223 + t494;
t193 = -t461 * t215 + t465 * t216 - t536;
t192 = t465 * t215 + t461 * t216 + t308;
t191 = -t461 * t213 + t465 * t214;
t190 = t465 * t213 + t461 * t214;
t189 = pkin(7) * t227 - (-pkin(3) * t464 - pkin(8) * t460 - pkin(2)) * t246;
t188 = -t461 * t210 + t465 * t212 - t537;
t187 = -t461 * t209 + t465 * t211 - t584;
t186 = t465 * t210 + t461 * t212 + t494;
t185 = t465 * t209 + t461 * t211 + t582;
t184 = -t462 * t203 + t466 * t204;
t183 = t466 * t203 + t462 * t204;
t182 = -pkin(7) * t213 + t464 * t205 - t460 * t208;
t181 = -t462 * t190 + t466 * t191;
t180 = t466 * t190 + t462 * t191;
t179 = -pkin(2) * t228 + pkin(7) * t214 + t460 * t205 + t464 * t208;
t178 = -pkin(6) * t203 - t461 * t189 + t465 * t198;
t177 = pkin(1) * t246 + pkin(6) * t204 + t465 * t189 + t461 * t198;
t176 = -pkin(6) * t190 - t461 * t179 + t465 * t182;
t175 = -pkin(1) * t228 + pkin(6) * t191 + t465 * t179 + t461 * t182;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t438, -t439, 0, t391, 0, 0, 0, 0, 0, 0, -t480, t366, 0, t266, 0, 0, 0, 0, 0, 0, -t550, t314, 0, t201, 0, 0, 0, 0, 0, 0, t236, t238, t254, t184, 0, 0, 0, 0, 0, 0, t236, t254, -t238, t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t439, -t438, 0, t390, 0, 0, 0, 0, 0, 0, -t366, -t480, 0, -t559, 0, 0, 0, 0, 0, 0, -t314, -t550, 0, -t570, 0, 0, 0, 0, 0, 0, t233, -t234, t253, t183, 0, 0, 0, 0, 0, 0, t233, t253, t234, t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t376, -t377, 0, -t246, 0, 0, 0, 0, 0, 0, t376, 0, t377, t228; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t439, 0, -t438, 0, t484, -t419, -t390, -pkin(5) * t390, 0, 0, -t366, 0, -t480, 0, t569, t568, t559, pkin(5) * t559 + pkin(6) * t504 - t462 * t319, 0, 0, -t314, 0, -t550, 0, t577, t576, t570, pkin(5) * t570 - t462 * t199 + t466 * t202, t244, t580, t249, t243, t578, t274, -t462 * t194 + t466 * t196 - t539, -t462 * t195 + t466 * t197 + t586, -t462 * t206 + t466 * t207 - t538, -pkin(5) * t183 - t462 * t177 + t466 * t178, t244, t249, -t580, t274, -t578, t243, -t462 * t186 + t466 * t188 - t539, -t462 * t192 + t466 * t193 - t538, -t462 * t185 + t466 * t187 - t586, -pkin(5) * t180 - t462 * t175 + t466 * t176; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t438, 0, t439, 0, t419, t484, t391, pkin(5) * t391, 0, 0, t480, 0, -t366, 0, -t568, t569, t266, pkin(5) * t266 + pkin(6) * t513 + t466 * t319, 0, 0, t550, 0, -t314, 0, -t576, t577, t201, pkin(5) * t201 + t466 * t199 + t462 * t202, t242, -t581, t248, t241, -t579, t273, t466 * t194 + t462 * t196 + t232, t466 * t195 + t462 * t197 + t585, t466 * t206 + t462 * t207 + t252, pkin(5) * t184 + t466 * t177 + t462 * t178, t242, t248, t581, t273, t579, t241, t466 * t186 + t462 * t188 + t232, t466 * t192 + t462 * t193 + t252, t466 * t185 + t462 * t187 - t585, pkin(5) * t181 + t466 * t175 + t462 * t176; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t443, t444, 0, 0, 0, 0, 0, 0, 0, t499, pkin(1) * t427 + t470, -pkin(1) * t425 - t375, 0, -pkin(1) * t322, 0, 0, 0, 0, 0, t448, -pkin(1) * t358 - pkin(2) * t416 - t317, -pkin(1) * t353 - pkin(2) * t413 - t318, 0, -pkin(1) * t221 - pkin(2) * t260, t368, t481, t379, t369, t543, 0, t486 - t512, -pkin(3) * t408 + t525 - t583, t247 + t485, pkin(1) * t203 + pkin(2) * t226 - pkin(3) * t306 + pkin(8) * t247, t368, t379, -t481, 0, -t543, t369, qJ(5) * t524 + t463 * t256 + t486, t463 * t263 + t459 * t264 + t485, t459 * t255 + (pkin(3) + t540) * t408 + t583, pkin(1) * t190 + pkin(2) * t213 + pkin(8) * t229 + (-pkin(3) + t483) * t262;];
tauB_reg = t1;
