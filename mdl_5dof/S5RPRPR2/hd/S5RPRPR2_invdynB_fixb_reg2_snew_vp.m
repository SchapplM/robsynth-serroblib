% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:10
% EndTime: 2020-01-03 11:34:25
% DurationCPUTime: 9.53s
% Computational Cost: add. (36153->449), mult. (54314->677), div. (0->0), fcn. (34125->10), ass. (0->303)
t498 = cos(qJ(3));
t486 = qJDD(1) + qJDD(3);
t495 = sin(qJ(3));
t535 = t495 * t486;
t488 = (qJD(1) + qJD(3));
t550 = t488 ^ 2;
t454 = t498 * t550 + t535;
t528 = t498 * t486;
t534 = t495 * t550;
t457 = -t528 + t534;
t491 = sin(pkin(8));
t493 = cos(pkin(8));
t400 = t493 * t454 - t491 * t457;
t489 = g(1) - qJDD(2);
t430 = pkin(6) * t454 - t498 * t489;
t568 = pkin(6) * t457 - t495 * t489;
t335 = qJ(2) * t400 + t493 * t430 - t491 * t568;
t404 = t491 * t454 + t493 * t457;
t496 = sin(qJ(1));
t499 = cos(qJ(1));
t350 = t499 * t400 - t496 * t404;
t578 = qJ(2) * t404 + t491 * t430 + t493 * t568;
t585 = -pkin(5) * t350 - t499 * t335 + t496 * t578;
t567 = t496 * t400 + t499 * t404;
t584 = pkin(5) * t567 + t496 * t335 + t499 * t578;
t471 = t499 * g(2) + t496 * g(3);
t506 = qJDD(1) * pkin(1) - t471;
t470 = t496 * g(2) - t499 * g(3);
t549 = qJD(1) ^ 2;
t507 = -t549 * pkin(1) - t470;
t409 = t491 * t506 + t493 * t507;
t407 = -t549 * pkin(2) + t409;
t505 = -t491 * t507 + t493 * t506;
t504 = qJDD(1) * pkin(2) + t505;
t353 = t495 * t407 - t498 * t504;
t354 = t498 * t407 + t495 * t504;
t514 = t495 * t353 + t498 * t354;
t304 = t498 * t353 - t495 * t354;
t539 = t493 * t304;
t255 = -t491 * t514 + t539;
t542 = t491 * t304;
t574 = t493 * t514 + t542;
t581 = t496 * t255 + t499 * t574;
t220 = -t499 * t255 + t496 * t574;
t492 = cos(pkin(9));
t490 = sin(pkin(9));
t487 = t490 ^ 2;
t502 = t492 ^ 2;
t552 = t550 * (t487 + t502);
t446 = t492 * t552;
t517 = t492 * t528;
t412 = -t495 * t446 + t517;
t414 = t498 * t446 + t492 * t535;
t369 = t493 * t412 - t491 * t414;
t372 = t491 * t412 + t493 * t414;
t315 = t499 * t369 - t496 * t372;
t317 = t496 * t369 + t499 * t372;
t570 = -t550 * pkin(3) + t486 * qJ(4) + (2 * qJD(4) * t488) + t354;
t513 = t493 * t409 - t491 * t505;
t363 = -t491 * t409 - t493 * t505;
t526 = t499 * t363;
t310 = t496 * t513 - t526;
t532 = t496 * t363;
t575 = t499 * t513 + t532;
t462 = t491 * qJDD(1) + t493 * t549;
t463 = t493 * qJDD(1) - t491 * t549;
t417 = -t499 * t462 - t496 * t463;
t440 = qJ(2) * t462 - t493 * t489;
t510 = -qJ(2) * t463 - t491 * t489;
t573 = pkin(5) * t417 - t499 * t440 + t496 * t510;
t494 = sin(qJ(5));
t497 = cos(qJ(5));
t569 = t490 * t494 - t492 * t497;
t441 = t569 * t488;
t509 = t490 * t497 + t492 * t494;
t443 = t509 * t488;
t392 = t443 * t441;
t553 = qJDD(5) - t392;
t572 = t494 * t553;
t571 = t497 * t553;
t554 = t496 * t462 - t499 * t463;
t566 = pkin(5) * t554 + t496 * t440 + t499 * t510;
t544 = t490 * t492;
t420 = t454 * t544;
t421 = t490 * t517 - t534 * t544;
t377 = t493 * t420 + t491 * t421;
t380 = t491 * t420 - t493 * t421;
t564 = t499 * t377 - t496 * t380;
t563 = t496 * t377 + t499 * t380;
t555 = t509 * t486;
t477 = t502 * t550;
t546 = t487 * t550;
t452 = t477 + t546;
t435 = t441 ^ 2;
t436 = t443 ^ 2;
t548 = pkin(1) * t489;
t480 = t486 * pkin(3);
t547 = t486 * t492;
t479 = t492 * t489;
t321 = -t479 + (pkin(4) * t550 * t492 - pkin(7) * t486 - t570) * t490;
t329 = -t490 * t489 + t570 * t492;
t322 = -pkin(4) * t477 + pkin(7) * t547 + t329;
t265 = -t497 * t321 + t494 * t322;
t266 = t494 * t321 + t497 * t322;
t236 = -t497 * t265 + t494 * t266;
t545 = t490 * t236;
t541 = t492 * t236;
t341 = -t550 * qJ(4) + qJDD(4) + t353 - t480;
t325 = -pkin(4) * t547 - t452 * pkin(7) + t341;
t538 = t494 * t325;
t385 = qJDD(5) + t392;
t537 = t494 * t385;
t536 = t495 * t341;
t531 = t497 * t325;
t530 = t497 * t385;
t529 = t498 * t341;
t522 = t441 * qJD(5);
t521 = t443 * qJD(5);
t519 = t495 * t392;
t518 = t498 * t392;
t516 = -t341 + t480;
t464 = -t496 * qJDD(1) - t499 * t549;
t515 = pkin(5) * t464 + t499 * g(1);
t237 = t494 * t265 + t497 * t266;
t328 = t570 * t490 + t479;
t281 = t490 * t328 + t492 * t329;
t422 = -t496 * t470 - t499 * t471;
t433 = t569 * t486;
t280 = t492 * t328 - t490 * t329;
t423 = t499 * t470 - t496 * t471;
t500 = qJD(5) ^ 2;
t476 = t502 * t486;
t475 = t487 * t486;
t465 = t499 * qJDD(1) - t496 * t549;
t453 = t477 - t546;
t450 = t476 - t475;
t449 = t476 + t475;
t447 = pkin(5) * t465 + t496 * g(1);
t445 = t490 * t552;
t426 = -t436 - t500;
t425 = -t436 + t500;
t424 = t435 - t500;
t413 = t498 * t445 + t490 * t535;
t410 = t495 * t445 - t490 * t528;
t399 = t498 * t450 - t495 * t453;
t398 = t498 * t449 - t495 * t452;
t397 = t495 * t450 + t498 * t453;
t396 = t495 * t449 + t498 * t452;
t391 = -t436 + t435;
t390 = t555 - t522;
t389 = t555 - 0.2e1 * t522;
t388 = -t433 - t521;
t387 = t433 + 0.2e1 * t521;
t383 = -t500 - t435;
t382 = (-t441 * t497 + t443 * t494) * qJD(5);
t381 = (-t441 * t494 - t443 * t497) * qJD(5);
t376 = -t435 - t436;
t374 = t497 * t390 - t494 * t521;
t371 = -t491 * t410 + t493 * t413;
t370 = t494 * t390 + t497 * t521;
t367 = t493 * t410 + t491 * t413;
t366 = -t494 * t388 + t497 * t522;
t365 = t497 * t388 + t494 * t522;
t360 = -t494 * t426 - t530;
t359 = -t494 * t425 + t571;
t358 = t497 * t424 - t537;
t357 = t497 * t426 - t537;
t356 = t497 * t425 + t572;
t355 = t494 * t424 + t530;
t352 = qJ(2) * t513 + t548;
t346 = -t491 * t397 + t493 * t399;
t345 = -t491 * t396 + t493 * t398;
t344 = t493 * t397 + t491 * t399;
t343 = t493 * t396 + t491 * t398;
t339 = -t497 * t387 - t494 * t389;
t338 = -t433 * t497 + t494 * t555;
t337 = -t494 * t387 + t497 * t389;
t336 = -t433 * t494 - t497 * t555;
t331 = t497 * t383 - t572;
t330 = t494 * t383 + t571;
t327 = -t490 * t381 + t492 * t382;
t324 = t495 * qJDD(5) + t498 * t327;
t323 = -t498 * qJDD(5) + t495 * t327;
t316 = t496 * t367 - t499 * t371;
t314 = t499 * t367 + t496 * t371;
t312 = -t490 * t370 + t492 * t374;
t311 = -t490 * t365 + t492 * t366;
t309 = -t490 * t357 + t492 * t360;
t308 = -t490 * t356 + t492 * t359;
t307 = -t490 * t355 + t492 * t358;
t306 = t492 * t357 + t490 * t360;
t301 = pkin(2) * t489 + pkin(6) * t514;
t300 = t496 * t343 - t499 * t345;
t299 = t499 * t343 + t496 * t345;
t298 = t498 * t308 + t495 * t555;
t297 = t498 * t307 - t495 * t433;
t296 = t495 * t308 - t498 * t555;
t295 = t495 * t307 + t498 * t433;
t294 = -pkin(7) * t357 + t531;
t293 = -t490 * t337 + t492 * t339;
t292 = -t490 * t336 + t492 * t338;
t291 = t492 * t336 + t490 * t338;
t290 = -t490 * t330 + t492 * t331;
t289 = t492 * t330 + t490 * t331;
t288 = t498 * t312 + t519;
t287 = t498 * t311 - t519;
t286 = t495 * t312 - t518;
t285 = t495 * t311 + t518;
t284 = t498 * t309 + t495 * t389;
t283 = t495 * t309 - t498 * t389;
t282 = -pkin(7) * t330 + t538;
t278 = -pkin(4) * t389 + pkin(7) * t360 + t538;
t277 = -pkin(6) * t410 - t495 * t329 + t492 * t529;
t276 = -pkin(6) * t412 - t495 * t328 + t490 * t529;
t275 = pkin(6) * t413 + t498 * t329 + t492 * t536;
t274 = -pkin(6) * t414 + t498 * t328 + t490 * t536;
t273 = -t491 * t323 + t493 * t324;
t272 = t493 * t323 + t491 * t324;
t271 = t498 * t293 - t495 * t391;
t270 = t495 * t293 + t498 * t391;
t269 = t498 * t290 + t495 * t387;
t268 = t495 * t290 - t498 * t387;
t267 = -pkin(4) * t387 + pkin(7) * t331 - t531;
t263 = t498 * t292 + t495 * t376;
t262 = t495 * t292 - t498 * t376;
t261 = -pkin(6) * t396 + t498 * t280;
t260 = pkin(6) * t398 + t495 * t280;
t259 = -pkin(3) * t291 - pkin(4) * t336;
t258 = t498 * t281 + t536;
t257 = t495 * t281 - t529;
t252 = -t491 * t296 + t493 * t298;
t251 = -t491 * t295 + t493 * t297;
t250 = t493 * t296 + t491 * t298;
t249 = t493 * t295 + t491 * t297;
t248 = -t491 * t286 + t493 * t288;
t247 = -t491 * t285 + t493 * t287;
t246 = t493 * t286 + t491 * t288;
t245 = t493 * t285 + t491 * t287;
t244 = -t491 * t283 + t493 * t284;
t243 = t493 * t283 + t491 * t284;
t242 = -pkin(3) * t306 - pkin(4) * t357 + t266;
t241 = -t491 * t270 + t493 * t271;
t240 = t493 * t270 + t491 * t271;
t239 = -t491 * t268 + t493 * t269;
t238 = t493 * t268 + t491 * t269;
t235 = -t491 * t262 + t493 * t263;
t234 = t493 * t262 + t491 * t263;
t233 = -pkin(3) * t289 - pkin(4) * t330 + t265;
t232 = -qJ(2) * t367 - t491 * t275 + t493 * t277;
t231 = -qJ(2) * t369 - t491 * t274 + t493 * t276;
t230 = qJ(2) * t371 + t493 * t275 + t491 * t277;
t229 = -qJ(2) * t372 + t493 * t274 + t491 * t276;
t228 = -qJ(4) * t306 - t490 * t278 + t492 * t294;
t227 = -pkin(7) * t336 - t236;
t226 = -qJ(2) * t343 - t491 * t260 + t493 * t261;
t225 = qJ(2) * t345 + t493 * t260 + t491 * t261;
t224 = -pkin(4) * t325 + pkin(7) * t237;
t223 = -t491 * t257 + t493 * t258;
t222 = t493 * t257 + t491 * t258;
t219 = pkin(6) * t539 + qJ(2) * t255 - t491 * t301;
t218 = -qJ(4) * t289 - t490 * t267 + t492 * t282;
t217 = pkin(6) * t542 + qJ(2) * t574 + t493 * t301 + t548;
t216 = -pkin(4) * t376 + pkin(7) * t338 + t237;
t215 = -pkin(6) * t257 - (pkin(3) * t495 - qJ(4) * t498) * t280;
t214 = t496 * t243 - t499 * t244;
t213 = t499 * t243 + t496 * t244;
t212 = pkin(6) * t258 - (-pkin(3) * t498 - qJ(4) * t495 - pkin(2)) * t280;
t211 = t496 * t238 - t499 * t239;
t210 = t499 * t238 + t496 * t239;
t209 = t492 * t237 - t545;
t208 = t490 * t237 + t541;
t207 = t496 * t234 - t499 * t235;
t206 = t499 * t234 + t496 * t235;
t205 = t498 * t209 + t495 * t325;
t204 = t495 * t209 - t498 * t325;
t203 = -pkin(6) * t283 + t498 * t228 - t495 * t242;
t202 = t496 * t222 - t499 * t223;
t201 = t499 * t222 + t496 * t223;
t200 = -pkin(2) * t306 + pkin(6) * t284 + t495 * t228 + t498 * t242;
t199 = -pkin(6) * t268 + t498 * t218 - t495 * t233;
t198 = -qJ(4) * t291 - t490 * t216 + t492 * t227;
t197 = -pkin(3) * t208 - pkin(4) * t236;
t196 = -pkin(2) * t289 + pkin(6) * t269 + t495 * t218 + t498 * t233;
t195 = -pkin(6) * t262 + t498 * t198 - t495 * t259;
t194 = -pkin(7) * t541 - qJ(4) * t208 - t490 * t224;
t193 = -pkin(2) * t291 + pkin(6) * t263 + t495 * t198 + t498 * t259;
t192 = -t491 * t204 + t493 * t205;
t191 = t493 * t204 + t491 * t205;
t190 = -qJ(2) * t222 - t491 * t212 + t493 * t215;
t189 = pkin(1) * t280 + qJ(2) * t223 + t493 * t212 + t491 * t215;
t188 = -qJ(2) * t243 - t491 * t200 + t493 * t203;
t187 = -pkin(1) * t306 + qJ(2) * t244 + t493 * t200 + t491 * t203;
t186 = -qJ(2) * t238 - t491 * t196 + t493 * t199;
t185 = -pkin(1) * t289 + qJ(2) * t239 + t493 * t196 + t491 * t199;
t184 = t496 * t191 - t499 * t192;
t183 = t499 * t191 + t496 * t192;
t182 = -qJ(2) * t234 - t491 * t193 + t493 * t195;
t181 = -pkin(6) * t204 + t498 * t194 - t495 * t197;
t180 = -pkin(1) * t291 + qJ(2) * t235 + t493 * t193 + t491 * t195;
t179 = -pkin(2) * t208 + pkin(6) * t205 + t495 * t194 + t498 * t197;
t178 = -qJ(2) * t191 - t491 * t179 + t493 * t181;
t177 = -pkin(1) * t208 + qJ(2) * t192 + t493 * t179 + t491 * t181;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t489, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t489, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, 0, 0, 0, 0, 0, 0, t289, t306, t291, t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t465, t464, 0, t422, 0, 0, 0, 0, 0, 0, -t554, t417, 0, t310, 0, 0, 0, 0, 0, 0, -t567, -t350, 0, t220, 0, 0, 0, 0, 0, 0, t315, t314, t299, t201, 0, 0, 0, 0, 0, 0, t210, t213, t206, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t464, t465, 0, t423, 0, 0, 0, 0, 0, 0, -t417, -t554, 0, -t575, 0, 0, 0, 0, 0, 0, t350, -t567, 0, -t581, 0, 0, 0, 0, 0, 0, t317, t316, t300, t202, 0, 0, 0, 0, 0, 0, t211, t214, t207, t184; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t471, t470, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t463 + t505, -pkin(1) * t462 - t409, 0, -pkin(1) * t363, 0, 0, 0, 0, 0, t486, -pkin(1) * t404 - pkin(2) * t457 - t353, -pkin(1) * t400 - pkin(2) * t454 - t354, 0, -pkin(1) * t255 - pkin(2) * t304, t475, 0.2e1 * t486 * t544, 0, t476, 0, 0, pkin(1) * t369 + pkin(2) * t412 - qJ(4) * t446 + t492 * t516, pkin(1) * t367 + pkin(2) * t410 + qJ(4) * t445 - t490 * t516, pkin(1) * t343 + pkin(2) * t396 + pkin(3) * t452 + qJ(4) * t449 + t281, pkin(1) * t222 + pkin(2) * t257 - pkin(3) * t341 + qJ(4) * t281, t492 * t370 + t490 * t374, t492 * t337 + t490 * t339, t492 * t356 + t490 * t359, t492 * t365 + t490 * t366, t492 * t355 + t490 * t358, t492 * t381 + t490 * t382, pkin(1) * t238 + pkin(2) * t268 - pkin(3) * t387 + qJ(4) * t290 + t492 * t267 + t490 * t282, pkin(1) * t243 + pkin(2) * t283 - pkin(3) * t389 + qJ(4) * t309 + t492 * t278 + t490 * t294, pkin(1) * t234 + pkin(2) * t262 - pkin(3) * t376 + qJ(4) * t292 + t492 * t216 + t490 * t227, pkin(1) * t191 + pkin(2) * t204 - pkin(3) * t325 - pkin(7) * t545 + qJ(4) * t209 + t492 * t224; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t464, 0, t465, 0, t515, -t447, -t423, -pkin(5) * t423, 0, 0, -t417, 0, -t554, 0, t573, t566, t575, pkin(5) * t575 + qJ(2) * t532 + t499 * t352, 0, 0, t350, 0, -t567, 0, t585, t584, t581, pkin(5) * t581 + t499 * t217 + t496 * t219, t564, t499 * t344 + t496 * t346, t314, -t564, -t315, 0, -pkin(5) * t317 + t499 * t229 + t496 * t231, -pkin(5) * t316 + t499 * t230 + t496 * t232, -pkin(5) * t300 + t499 * t225 + t496 * t226, -pkin(5) * t202 + t499 * t189 + t496 * t190, t499 * t246 + t496 * t248, t499 * t240 + t496 * t241, t499 * t250 + t496 * t252, t499 * t245 + t496 * t247, t499 * t249 + t496 * t251, t499 * t272 + t496 * t273, -pkin(5) * t211 + t499 * t185 + t496 * t186, -pkin(5) * t214 + t499 * t187 + t496 * t188, -pkin(5) * t207 + t499 * t180 + t496 * t182, -pkin(5) * t184 + t499 * t177 + t496 * t178; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t465, 0, -t464, 0, t447, t515, t422, pkin(5) * t422, 0, 0, t554, 0, -t417, 0, -t566, t573, t310, pkin(5) * t310 - qJ(2) * t526 + t496 * t352, 0, 0, t567, 0, t350, 0, -t584, t585, t220, pkin(5) * t220 + t496 * t217 - t499 * t219, t563, t496 * t344 - t499 * t346, t316, -t563, -t317, 0, pkin(5) * t315 + t496 * t229 - t499 * t231, pkin(5) * t314 + t496 * t230 - t499 * t232, pkin(5) * t299 + t496 * t225 - t499 * t226, pkin(5) * t201 + t496 * t189 - t499 * t190, t496 * t246 - t499 * t248, t496 * t240 - t499 * t241, t496 * t250 - t499 * t252, t496 * t245 - t499 * t247, t496 * t249 - t499 * t251, t496 * t272 - t499 * t273, pkin(5) * t210 + t496 * t185 - t499 * t186, pkin(5) * t213 + t496 * t187 - t499 * t188, pkin(5) * t206 + t496 * t180 - t499 * t182, pkin(5) * t183 + t496 * t177 - t499 * t178;];
tauB_reg = t1;
