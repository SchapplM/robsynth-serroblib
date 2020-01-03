% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRP6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:22
% EndTime: 2019-12-31 17:55:29
% DurationCPUTime: 7.20s
% Computational Cost: add. (8728->356), mult. (19880->473), div. (0->0), fcn. (12652->6), ass. (0->254)
t446 = sin(pkin(7));
t447 = cos(pkin(7));
t448 = sin(qJ(4));
t450 = cos(qJ(4));
t418 = (-t446 * t448 + t447 * t450) * qJD(1);
t413 = t418 ^ 2;
t452 = qJD(4) ^ 2;
t363 = t452 + t413;
t465 = t446 * t450 + t447 * t448;
t416 = t465 * qJD(1);
t495 = t418 * t416;
t366 = qJDD(4) + t495;
t573 = t366 * t448;
t307 = t363 * t450 + t573;
t572 = t366 * t450;
t309 = -t363 * t448 + t572;
t266 = t307 * t447 + t309 * t446;
t449 = sin(qJ(1));
t592 = t266 * t449;
t451 = cos(qJ(1));
t591 = t266 * t451;
t511 = -qJ(3) - pkin(1);
t590 = t511 * t266;
t284 = t307 * t446 - t309 * t447;
t589 = t511 * t284;
t518 = t416 ^ 2;
t396 = t518 - t452;
t314 = t396 * t448 + t572;
t321 = -t396 * t450 + t573;
t278 = t314 * t447 - t321 * t446;
t524 = t465 * qJDD(1);
t588 = t278 * t449 - t524 * t451;
t587 = t278 * t451 + t524 * t449;
t586 = pkin(2) * t266 + pkin(3) * t307 + qJ(2) * t284;
t584 = pkin(6) * t307;
t583 = pkin(6) * t309;
t582 = t314 * t446 + t321 * t447;
t397 = -t413 + t452;
t543 = qJDD(4) - t495;
t502 = t543 * t448;
t570 = t450 * t397 + t502;
t354 = t450 * t543;
t571 = -t397 * t448 + t354;
t579 = -t446 * t571 - t447 * t570;
t581 = t449 * t579;
t580 = t451 * t579;
t578 = -t446 * t570 + t447 * t571;
t342 = -t518 - t413;
t484 = qJDD(1) * t447;
t485 = qJDD(1) * t446;
t415 = -t448 * t485 + t450 * t484;
t529 = t415 * t448 - t450 * t524;
t533 = -t450 * t415 - t448 * t524;
t538 = t446 * t529 + t447 * t533;
t564 = t342 * t449 - t451 * t538;
t577 = pkin(5) * t564;
t565 = t342 * t451 + t449 * t538;
t576 = pkin(5) * t565;
t488 = t418 * qJD(4);
t369 = t524 + 0.2e1 * t488;
t525 = -t518 - t452;
t528 = t450 * t525 - t502;
t534 = t448 * t525 + t354;
t540 = t446 * t528 + t447 * t534;
t567 = t451 * t369 + t449 * t540;
t575 = pkin(5) * t567;
t568 = t449 * t369 - t451 * t540;
t574 = pkin(5) * t568;
t489 = qJD(4) * t416;
t372 = t415 - t489;
t544 = t372 - t489;
t569 = t544 * qJ(5);
t566 = qJ(2) * t369 + t511 * t540;
t541 = -t446 * t534 + t447 * t528;
t563 = pkin(2) * t369 + t511 * t541;
t562 = qJ(2) * t342 + t511 * t538;
t539 = -t446 * t533 + t447 * t529;
t561 = pkin(2) * t342 + t511 * t539;
t234 = pkin(2) * t538 + pkin(3) * t533 - qJ(2) * t539;
t560 = pkin(2) * t540 + pkin(3) * t534 - qJ(2) * t541;
t556 = pkin(6) * t528;
t555 = pkin(6) * t533;
t554 = pkin(6) * t534;
t375 = t413 - t518;
t548 = t375 * t449;
t547 = t375 * t451;
t453 = qJD(1) ^ 2;
t428 = g(1) * t449 - t451 * g(2);
t469 = qJDD(2) - t428;
t460 = -qJ(2) * t453 + t469;
t473 = -0.2e1 * qJD(1) * qJD(3) + t511 * qJDD(1) + t460;
t542 = -pkin(3) * t342 + pkin(6) * t529;
t461 = (-t416 * t448 - t418 * t450) * qJD(4);
t393 = t448 * t488;
t479 = t450 * t489;
t467 = t393 - t479;
t522 = -t446 * t467 - t447 * t461;
t532 = t449 * qJDD(4) + t451 * t522;
t531 = t451 * qJDD(4) - t449 * t522;
t480 = t451 * t495;
t370 = -t524 - t488;
t462 = -t370 * t448 + t479;
t468 = t450 * t370 + t448 * t489;
t519 = -t446 * t462 - t447 * t468;
t530 = -t449 * t519 - t480;
t481 = t449 * t495;
t527 = t451 * t519 - t481;
t440 = t446 ^ 2;
t441 = t447 ^ 2;
t490 = t440 + t441;
t526 = t490 * t453;
t361 = -g(3) * t447 + t473 * t446;
t521 = -t446 * t461 + t447 * t467;
t520 = -t446 * t468 + t447 * t462;
t517 = 2 * qJD(5);
t515 = pkin(3) * t453;
t514 = pkin(4) * t370;
t513 = pkin(4) * t450;
t512 = g(3) * t446;
t510 = qJDD(1) * pkin(1);
t341 = t512 + (-pkin(6) * qJDD(1) - t446 * t515 + t473) * t447;
t344 = -pkin(6) * t485 - t440 * t515 + t361;
t295 = -t450 * t341 + t344 * t448;
t296 = t448 * t341 + t450 * t344;
t248 = -t295 * t450 + t296 * t448;
t509 = t248 * t446;
t508 = t248 * t447;
t507 = t544 * t448;
t429 = g(1) * t451 + g(2) * t449;
t443 = qJDD(1) * qJ(2);
t464 = t429 - t443;
t486 = qJD(2) * qJD(1);
t458 = -qJDD(3) + t464 - 0.2e1 * t486;
t351 = -pkin(3) * t485 + (t490 * pkin(6) - t511) * t453 + t458;
t506 = t351 * t448;
t505 = t351 * t450;
t501 = t369 * t448;
t500 = t369 * t450;
t422 = t490 * qJDD(1);
t494 = t422 * t449;
t493 = t422 * t451;
t492 = t446 * t447;
t491 = t440 - t441;
t483 = qJDD(1) * t449;
t482 = qJDD(1) * t451;
t477 = -qJ(5) * t448 - pkin(3);
t249 = t295 * t448 + t450 * t296;
t438 = 0.2e1 * t486;
t401 = -pkin(1) * t453 + t438 - t464;
t406 = -t460 + t510;
t346 = t451 * t401 - t406 * t449;
t386 = -t428 * t449 - t451 * t429;
t388 = -t511 * t453 + t458;
t474 = -t388 + t443;
t426 = -t449 * t453 + t482;
t472 = pkin(5) * t426 + g(3) * t449;
t427 = t451 * t453 + t483;
t471 = -pkin(5) * t427 + g(3) * t451;
t330 = t372 * t448 + t450 * t488;
t331 = t372 * t450 - t393;
t288 = -t330 * t447 - t331 * t446;
t470 = t451 * t288 + t481;
t359 = pkin(4) * t416 - qJ(5) * t418;
t466 = qJDD(4) * qJ(5) + qJD(4) * t517 - t416 * t359 + t296;
t360 = t447 * t473 + t512;
t304 = t360 * t447 + t361 * t446;
t305 = -t360 * t446 + t361 * t447;
t345 = t401 * t449 + t406 * t451;
t385 = t428 * t451 - t429 * t449;
t463 = -t288 * t449 + t480;
t420 = t446 * t526;
t383 = -t420 * t449 + t446 * t482;
t381 = t420 * t451 + t446 * t483;
t459 = -qJDD(4) * pkin(4) - t452 * qJ(5) + t359 * t418 + qJDD(5) + t295;
t457 = -pkin(4) * t488 + t418 * t517 + t351;
t456 = t457 + t569;
t425 = t491 * t453;
t423 = t491 * qJDD(1);
t419 = t447 * t526;
t392 = t427 * t492;
t391 = t426 * t492;
t384 = -t419 * t449 + t447 * t482;
t382 = t419 * t451 + t447 * t483;
t380 = -t451 * t526 - t494;
t379 = -t449 * t526 + t493;
t374 = pkin(2) * t485 - t388 * t447;
t373 = pkin(2) * t484 + t388 * t446;
t371 = t415 - 0.2e1 * t489;
t339 = -pkin(2) * t420 + t360;
t338 = -pkin(2) * t419 - t361;
t336 = t372 + t489;
t313 = -t371 * t448 - t500;
t311 = t371 * t450 - t501;
t301 = -pkin(2) * t526 - t305;
t298 = t500 + t507;
t297 = -t450 * t544 + t501;
t294 = -t505 + t584;
t292 = t304 * t449 - t388 * t451;
t291 = -t304 * t451 - t388 * t449;
t290 = -t506 - t554;
t289 = -t330 * t446 + t331 * t447;
t276 = -pkin(3) * t371 - t506 - t583;
t272 = -t311 * t447 - t313 * t446;
t270 = -pkin(3) * t369 + t505 + t556;
t269 = pkin(2) * t304 - qJ(2) * t305;
t263 = -pkin(4) * t452 + t466;
t262 = t456 + t514;
t261 = t371 * t451 - t592;
t259 = t371 * t449 + t591;
t257 = -pkin(2) * t388 + t511 * t305;
t256 = -t297 * t447 - t298 * t446;
t253 = (-t369 + t370) * pkin(4) + t456;
t250 = -qJ(5) * t342 + t459;
t247 = (-t342 - t452) * pkin(4) + t466;
t246 = -t451 * t544 + t592;
t245 = -t449 * t544 - t591;
t244 = t457 + t514 + 0.2e1 * t569;
t243 = -pkin(4) * t415 - qJ(5) * t524 + t234;
t242 = pkin(3) * t351 + pkin(6) * t249;
t241 = -t248 - t555;
t240 = -qJ(5) * t500 - t253 * t448 - t554;
t239 = t263 * t450 + t448 * t459;
t238 = t263 * t448 - t450 * t459;
t237 = t249 + t542;
t236 = t253 * t450 + t477 * t369 + t556;
t235 = -pkin(4) * t507 + t244 * t450 - t584;
t233 = t583 + t244 * t448 + (pkin(3) + t513) * t544;
t232 = -t296 - t586;
t231 = t249 * t447 - t509;
t230 = t249 * t446 + t508;
t229 = -t247 * t448 + t250 * t450 - t555;
t228 = t230 * t449 - t351 * t451;
t227 = -t230 * t451 - t351 * t449;
t226 = -t295 + t560;
t225 = t247 * t450 + t250 * t448 + t542;
t224 = pkin(4) * t543 + qJ(5) * t525 - t459 + t560;
t223 = pkin(2) * t371 - t276 * t447 - t294 * t446 + t589;
t222 = qJ(5) * t366 + (t363 - t452) * pkin(4) + t466 + t586;
t221 = -t270 * t447 - t290 * t446 + t563;
t220 = -t238 * t446 + t239 * t447;
t219 = t238 * t447 + t239 * t446;
t218 = -pkin(6) * t238 + (-pkin(4) * t448 + qJ(5) * t450) * t262;
t217 = pkin(6) * t239 + (-t477 + t513) * t262;
t216 = t219 * t449 - t262 * t451;
t215 = -t219 * t451 - t262 * t449;
t214 = -t236 * t447 - t240 * t446 + t563;
t213 = -t237 * t447 - t241 * t446 + t561;
t212 = pkin(2) * t230 + pkin(3) * t248 - qJ(2) * t231;
t211 = -pkin(2) * t544 - t233 * t447 - t235 * t446 - t589;
t210 = -t225 * t447 - t229 * t446 + t561;
t209 = -pkin(2) * t351 + pkin(6) * t509 + t511 * t231 - t242 * t447;
t208 = pkin(2) * t219 + pkin(3) * t238 - pkin(4) * t459 - qJ(2) * t220 + qJ(5) * t263;
t207 = -pkin(2) * t262 - t217 * t447 - t218 * t446 + t511 * t220;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t427, -t426, 0, t386, 0, 0, 0, 0, 0, 0, 0, t427, t426, t346, 0, 0, 0, 0, 0, 0, t383, t384, t380, t292, 0, 0, 0, 0, 0, 0, t567, t261, t565, t228, 0, 0, 0, 0, 0, 0, t567, t565, t246, t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t426, -t427, 0, t385, 0, 0, 0, 0, 0, 0, 0, -t426, t427, t345, 0, 0, 0, 0, 0, 0, t381, t382, t379, t291, 0, 0, 0, 0, 0, 0, t568, t259, t564, t227, 0, 0, 0, 0, 0, 0, t568, t564, t245, t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, 0, 0, 0, 0, 0, 0, t541, t284, t539, t231, 0, 0, 0, 0, 0, 0, t541, t539, -t284, t220; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t426, 0, -t427, 0, -t472, -t471, -t385, -pkin(5) * t385, 0, -t426, t427, 0, 0, 0, -t345, t472, t471, -pkin(5) * t345 + (-pkin(1) * t449 + qJ(2) * t451) * g(3), t392, -t423 * t449 - t425 * t451, t384, -t392, -t383, 0, -pkin(5) * t381 + t339 * t451 - t374 * t449, -pkin(5) * t382 + t338 * t451 - t373 * t449, -pkin(2) * t493 - pkin(5) * t379 - t301 * t449, -pkin(5) * t291 - t257 * t449 + t269 * t451, t463, -t272 * t449 + t547, t415 * t451 - t581, t530, t588, t531, -t221 * t449 + t226 * t451 - t574, -pkin(5) * t259 - t223 * t449 + t232 * t451, -t213 * t449 + t234 * t451 - t577, -pkin(5) * t227 - t209 * t449 + t212 * t451, t463, t336 * t451 - t581, -t256 * t449 - t547, t531, -t588, t530, -t214 * t449 + t224 * t451 - t574, -t210 * t449 + t243 * t451 - t577, -pkin(5) * t245 - t211 * t449 + t222 * t451, -pkin(5) * t215 - t207 * t449 + t208 * t451; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t427, 0, t426, 0, t471, -t472, t386, pkin(5) * t386, 0, -t427, -t426, 0, 0, 0, t346, -t471, t472, pkin(5) * t346 + (pkin(1) * t451 + qJ(2) * t449) * g(3), -t391, t423 * t451 - t425 * t449, t382, t391, -t381, 0, pkin(5) * t383 + t339 * t449 + t374 * t451, pkin(5) * t384 + t338 * t449 + t373 * t451, -pkin(2) * t494 + pkin(5) * t380 + t301 * t451, pkin(5) * t292 + t257 * t451 + t269 * t449, t470, t272 * t451 + t548, t415 * t449 + t580, t527, -t587, t532, t221 * t451 + t226 * t449 + t575, pkin(5) * t261 + t223 * t451 + t232 * t449, t213 * t451 + t234 * t449 + t576, pkin(5) * t228 + t209 * t451 + t212 * t449, t470, t336 * t449 + t580, t256 * t451 - t548, t532, t587, t527, t214 * t451 + t224 * t449 + t575, t210 * t451 + t243 * t449 + t576, pkin(5) * t246 + t211 * t451 + t222 * t449, pkin(5) * t216 + t207 * t451 + t208 * t449; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t428, t429, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t469 - 0.2e1 * t510, -t429 + t438 + 0.2e1 * t443, pkin(1) * t406 + qJ(2) * t401, t441 * qJDD(1), -0.2e1 * t446 * t484, 0, t440 * qJDD(1), 0, 0, -t511 * t420 + t474 * t446, -t511 * t419 + t474 * t447, -qJ(2) * t526 - t511 * t422 - t304, -qJ(2) * t388 + t511 * t304, t289, -t311 * t446 + t313 * t447, t578, t520, -t582, t521, -t270 * t446 + t290 * t447 + t566, qJ(2) * t371 - t276 * t446 + t294 * t447 - t590, -t237 * t446 + t241 * t447 + t562, -pkin(6) * t508 - qJ(2) * t351 + t511 * t230 - t242 * t446, t289, t578, -t297 * t446 + t298 * t447, t521, t582, t520, -t236 * t446 + t240 * t447 + t566, -t225 * t446 + t229 * t447 + t562, -qJ(2) * t544 - t233 * t446 + t235 * t447 + t590, -qJ(2) * t262 - t217 * t446 + t218 * t447 + t511 * t219;];
tauB_reg = t1;
