% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPP3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:55
% EndTime: 2019-12-31 16:58:01
% DurationCPUTime: 5.66s
% Computational Cost: add. (7776->354), mult. (18481->468), div. (0->0), fcn. (11802->6), ass. (0->271)
t450 = sin(pkin(6));
t451 = cos(pkin(6));
t452 = sin(qJ(2));
t454 = cos(qJ(2));
t413 = (t450 * t454 + t451 * t452) * qJD(1);
t410 = t413 ^ 2;
t456 = qJD(2) ^ 2;
t364 = t456 + t410;
t493 = qJD(1) * t452;
t411 = -t451 * t454 * qJD(1) + t450 * t493;
t374 = t413 * t411;
t528 = qJDD(2) + t374;
t515 = t450 * t528;
t302 = t451 * t364 + t515;
t510 = t451 * t528;
t318 = t450 * t364 - t510;
t283 = t452 * t302 + t454 * t318;
t453 = sin(qJ(1));
t455 = cos(qJ(1));
t487 = qJD(1) * qJD(2);
t477 = t454 * t487;
t485 = t452 * qJDD(1);
t422 = t477 + t485;
t442 = t454 * qJDD(1);
t478 = t452 * t487;
t423 = t442 - t478;
t376 = t451 * t422 + t450 * t423;
t406 = qJD(2) * t411;
t532 = t376 - t406;
t262 = t453 * t283 - t455 * t532;
t598 = pkin(4) * t262;
t264 = t455 * t283 + t453 * t532;
t597 = pkin(4) * t264;
t269 = t454 * t302 - t452 * t318;
t596 = pkin(5) * t269;
t595 = -pkin(1) * t269 - pkin(2) * t302;
t594 = pkin(1) * t532 - pkin(5) * t283;
t475 = t450 * t422 - t451 * t423;
t492 = qJD(2) * t413;
t342 = t475 + t492;
t293 = -t450 * t342 + t451 * t532;
t511 = t451 * t342;
t517 = t450 * t532;
t295 = t511 + t517;
t253 = t452 * t293 + t454 * t295;
t523 = t411 ^ 2;
t371 = t410 - t523;
t593 = t453 * t253 + t455 * t371;
t592 = t455 * t253 - t453 * t371;
t395 = t523 - t456;
t307 = t450 * t395 + t510;
t313 = t451 * t395 - t515;
t279 = t452 * t307 - t454 * t313;
t344 = -t475 + t492;
t591 = t453 * t279 + t455 * t344;
t590 = t455 * t279 - t453 * t344;
t340 = -t523 - t410;
t533 = -t376 - t406;
t540 = t451 * t344 - t450 * t533;
t541 = t450 * t344 + t451 * t533;
t550 = -t452 * t541 + t454 * t540;
t568 = t453 * t340 + t455 * t550;
t588 = pkin(4) * t568;
t529 = qJDD(2) - t374;
t514 = t450 * t529;
t530 = -t523 - t456;
t537 = t451 * t530 - t514;
t361 = t451 * t529;
t539 = t450 * t530 + t361;
t552 = -t452 * t539 + t454 * t537;
t569 = t453 * t342 + t455 * t552;
t587 = pkin(4) * t569;
t570 = -t455 * t340 + t453 * t550;
t586 = pkin(4) * t570;
t571 = -t455 * t342 + t453 * t552;
t585 = pkin(4) * t571;
t584 = qJ(3) * t302;
t583 = qJ(3) * t318;
t582 = -t454 * t293 + t452 * t295;
t396 = -t410 + t456;
t554 = t451 * t396 + t514;
t555 = -t450 * t396 + t361;
t566 = -t452 * t554 + t454 * t555;
t581 = t453 * t566 + t455 * t533;
t580 = -t453 * t533 + t455 * t566;
t579 = t454 * t307 + t452 * t313;
t549 = t452 * t540 + t454 * t541;
t578 = pkin(5) * t549;
t551 = t452 * t537 + t454 * t539;
t577 = pkin(5) * t551;
t238 = -pkin(1) * t549 - pkin(2) * t541;
t574 = -pkin(1) * t551 - pkin(2) * t539;
t573 = -pkin(1) * t342 + pkin(5) * t552;
t572 = -pkin(1) * t340 + pkin(5) * t550;
t567 = t452 * t555 + t454 * t554;
t562 = qJ(3) * t537;
t561 = qJ(3) * t539;
t560 = qJ(3) * t541;
t553 = -pkin(2) * t340 + qJ(3) * t540;
t548 = 2 * qJD(4);
t542 = t532 * qJ(4);
t441 = t453 * qJDD(2);
t461 = (-t411 * t450 - t413 * t451) * qJD(2);
t491 = qJD(2) * t450;
t394 = t413 * t491;
t490 = qJD(2) * t451;
t479 = t411 * t490;
t466 = t394 - t479;
t526 = -t452 * t461 + t454 * t466;
t538 = t455 * t526 + t441;
t480 = t455 * t374;
t462 = t450 * t475 + t479;
t467 = t411 * t491 - t451 * t475;
t525 = -t452 * t467 + t454 * t462;
t536 = t453 * t525 + t480;
t482 = t455 * qJDD(2);
t535 = t453 * t526 - t482;
t481 = t453 * t374;
t534 = t455 * t525 - t481;
t362 = t411 * pkin(3) - t413 * qJ(4);
t457 = qJD(1) ^ 2;
t504 = t452 * t457;
t433 = t455 * g(1) + t453 * g(2);
t416 = -t457 * pkin(1) + qJDD(1) * pkin(5) - t433;
t507 = t452 * t416;
t329 = qJDD(2) * pkin(2) - t422 * qJ(3) - t507 + (pkin(2) * t504 + qJ(3) * t487 - g(3)) * t454;
t390 = -t452 * g(3) + t454 * t416;
t449 = t454 ^ 2;
t444 = t449 * t457;
t463 = qJD(2) * pkin(2) - qJ(3) * t493;
t330 = -pkin(2) * t444 + t423 * qJ(3) - qJD(2) * t463 + t390;
t497 = t450 * t329 + t451 * t330;
t531 = qJDD(2) * qJ(4) + qJD(2) * t548 - t411 * t362 + t497;
t527 = t452 * t466 + t454 * t461;
t524 = t452 * t462 + t454 * t467;
t522 = pkin(3) * t451;
t521 = t475 * pkin(3);
t448 = t452 ^ 2;
t520 = t448 * t457;
t432 = t453 * g(1) - t455 * g(2);
t465 = qJDD(1) * pkin(1) + t432;
t339 = t423 * pkin(2) - qJDD(3) - t463 * t493 + (qJ(3) * t449 + pkin(5)) * t457 + t465;
t519 = t450 * t339;
t512 = t451 * t339;
t488 = qJD(3) * t413;
t404 = 0.2e1 * t488;
t496 = -t451 * t329 + t450 * t330;
t273 = t404 + t496;
t489 = qJD(3) * t411;
t402 = -0.2e1 * t489;
t274 = t402 + t497;
t239 = -t451 * t273 + t450 * t274;
t509 = t452 * t239;
t415 = t457 * pkin(5) + t465;
t508 = t452 * t415;
t439 = t454 * t504;
t430 = qJDD(2) + t439;
t506 = t452 * t430;
t431 = qJDD(2) - t439;
t505 = t452 * t431;
t501 = t454 * t239;
t500 = t454 * t415;
t499 = t454 * t431;
t495 = -t340 - t456;
t494 = t448 + t449;
t484 = t453 * qJDD(1);
t483 = t455 * qJDD(1);
t476 = -qJ(4) * t450 - pkin(2);
t240 = t450 * t273 + t451 * t274;
t389 = t454 * g(3) + t507;
t333 = t452 * t389 + t454 * t390;
t382 = -t453 * t432 - t455 * t433;
t474 = t453 * t439;
t473 = t455 * t439;
t471 = t413 * t362 + qJDD(4) + t496;
t427 = -t453 * t457 + t483;
t470 = -pkin(4) * t427 - t453 * g(3);
t327 = t450 * t376 + t413 * t490;
t328 = t451 * t376 - t394;
t289 = -t452 * t327 + t454 * t328;
t469 = t453 * t289 - t480;
t468 = t455 * t289 + t481;
t332 = t454 * t389 - t452 * t390;
t381 = t455 * t432 - t453 * t433;
t464 = t402 + t531;
t460 = -qJDD(2) * pkin(3) + t471;
t459 = -pkin(3) * t492 + t413 * t548 + t339;
t458 = t459 + t542;
t438 = -t444 - t456;
t437 = t444 - t456;
t436 = -t456 - t520;
t435 = t456 - t520;
t429 = t444 - t520;
t428 = t444 + t520;
t426 = t455 * t457 + t484;
t425 = t494 * qJDD(1);
t424 = t442 - 0.2e1 * t478;
t421 = 0.2e1 * t477 + t485;
t419 = t454 * t430;
t418 = t494 * t487;
t407 = -pkin(4) * t426 + t455 * g(3);
t393 = t454 * t422 - t448 * t487;
t392 = -t452 * t423 - t449 * t487;
t388 = -t452 * t436 - t499;
t387 = -t452 * t435 + t419;
t386 = t454 * t438 - t506;
t385 = t454 * t437 - t505;
t384 = t454 * t436 - t505;
t383 = t452 * t438 + t419;
t379 = t455 * t425 - t453 * t428;
t378 = t453 * t425 + t455 * t428;
t377 = -t452 * t421 + t454 * t424;
t356 = t455 * t388 + t453 * t421;
t355 = t455 * t386 - t453 * t424;
t354 = t453 * t388 - t455 * t421;
t353 = t453 * t386 + t455 * t424;
t352 = -pkin(5) * t384 - t500;
t351 = -pkin(5) * t383 - t508;
t350 = -pkin(1) * t384 + t390;
t349 = -pkin(1) * t383 + t389;
t306 = t455 * t333 - t453 * t415;
t305 = t453 * t333 + t455 * t415;
t290 = -t512 + t584;
t286 = t454 * t327 + t452 * t328;
t275 = -t519 - t561;
t267 = -pkin(2) * t532 - t519 + t583;
t266 = t458 - t521;
t265 = -pkin(2) * t342 + t512 + t562;
t260 = t456 * qJ(4) - t460 - 0.2e1 * t488;
t259 = -t456 * pkin(3) + t464;
t248 = (-t342 - t475) * pkin(3) + t458;
t247 = t459 - t521 + 0.2e1 * t542;
t246 = t495 * qJ(4) + t404 + t460;
t245 = t495 * pkin(3) + t464;
t237 = t274 - t595;
t236 = pkin(2) * t339 + qJ(3) * t240;
t235 = -qJ(4) * t511 - t450 * t248 - t561;
t234 = -pkin(3) * t517 + t451 * t247 - t584;
t233 = t273 + t574;
t232 = -t239 - t560;
t231 = t451 * t248 + t342 * t476 + t562;
t230 = t451 * t259 - t450 * t260;
t229 = t450 * t259 + t451 * t260;
t228 = -t583 + t450 * t247 + (pkin(2) + t522) * t532;
t227 = -t452 * t267 + t454 * t290 + t596;
t226 = -pkin(3) * t533 - qJ(4) * t344 + t238;
t225 = t240 + t553;
t224 = t404 + (-t530 - t456) * qJ(4) + (-qJDD(2) - t529) * pkin(3) + t471 + t574;
t223 = -t452 * t265 + t454 * t275 - t577;
t222 = -qJ(4) * t528 + 0.2e1 * t489 + (-t364 + t456) * pkin(3) - t531 + t595;
t221 = t454 * t240 - t509;
t220 = t452 * t240 + t501;
t219 = -t450 * t245 + t451 * t246 - t560;
t218 = t455 * t221 - t453 * t339;
t217 = t453 * t221 + t455 * t339;
t216 = t451 * t245 + t450 * t246 + t553;
t215 = -qJ(3) * t229 + (-pkin(3) * t450 + qJ(4) * t451) * t266;
t214 = -t452 * t229 + t454 * t230;
t213 = t454 * t229 + t452 * t230;
t212 = -pkin(1) * t220 - pkin(2) * t239;
t211 = -t452 * t231 + t454 * t235 - t577;
t210 = qJ(3) * t230 + (-t476 + t522) * t266;
t209 = -t452 * t228 + t454 * t234 - t596;
t208 = t455 * t214 - t453 * t266;
t207 = t453 * t214 + t455 * t266;
t206 = -t452 * t225 + t454 * t232 - t578;
t205 = -pkin(5) * t220 - qJ(3) * t501 - t452 * t236;
t204 = -t452 * t216 + t454 * t219 - t578;
t203 = -pkin(1) * t213 - pkin(2) * t229 - pkin(3) * t260 - qJ(4) * t259;
t202 = -pkin(5) * t213 - t452 * t210 + t454 * t215;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t426, -t427, 0, t382, 0, 0, 0, 0, 0, 0, t355, t356, t379, t306, 0, 0, 0, 0, 0, 0, t569, t264, t568, t218, 0, 0, 0, 0, 0, 0, t569, t568, -t264, t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t427, -t426, 0, t381, 0, 0, 0, 0, 0, 0, t353, t354, t378, t305, 0, 0, 0, 0, 0, 0, t571, t262, t570, t217, 0, 0, 0, 0, 0, 0, t571, t570, -t262, t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t383, t384, 0, -t332, 0, 0, 0, 0, 0, 0, t551, -t269, t549, t220, 0, 0, 0, 0, 0, 0, t551, t549, t269, t213; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t427, 0, -t426, 0, t470, -t407, -t381, -pkin(4) * t381, t455 * t393 - t474, t455 * t377 - t453 * t429, t455 * t387 + t452 * t484, t455 * t392 + t474, t455 * t385 + t442 * t453, t455 * t418 + t441, -pkin(4) * t353 - t453 * t349 + t455 * t351, -pkin(4) * t354 - t453 * t350 + t455 * t352, -pkin(4) * t378 + t455 * t332, -pkin(4) * t305 - (pkin(1) * t453 - pkin(5) * t455) * t332, t468, -t592, t580, t534, -t590, t538, t455 * t223 - t453 * t233 - t585, t455 * t227 - t453 * t237 - t598, t455 * t206 - t453 * t238 - t586, -pkin(4) * t217 + t455 * t205 - t453 * t212, t468, t580, t592, t538, t590, t534, t455 * t211 - t453 * t224 - t585, t455 * t204 - t453 * t226 - t586, t455 * t209 - t453 * t222 + t598, -pkin(4) * t207 + t455 * t202 - t453 * t203; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t426, 0, t427, 0, t407, t470, t382, pkin(4) * t382, t453 * t393 + t473, t453 * t377 + t455 * t429, t453 * t387 - t452 * t483, t453 * t392 - t473, t453 * t385 - t454 * t483, t453 * t418 - t482, pkin(4) * t355 + t455 * t349 + t453 * t351, pkin(4) * t356 + t455 * t350 + t453 * t352, pkin(4) * t379 + t453 * t332, pkin(4) * t306 - (-pkin(1) * t455 - pkin(5) * t453) * t332, t469, -t593, t581, t536, -t591, t535, t453 * t223 + t455 * t233 + t587, t453 * t227 + t455 * t237 + t597, t453 * t206 + t455 * t238 + t588, pkin(4) * t218 + t453 * t205 + t455 * t212, t469, t581, t593, t535, t591, t536, t453 * t211 + t455 * t224 + t587, t453 * t204 + t455 * t226 + t588, t453 * t209 + t455 * t222 - t597, pkin(4) * t208 + t453 * t202 + t455 * t203; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t432, t433, 0, 0, (t422 + t477) * t452, t454 * t421 + t452 * t424, t454 * t435 + t506, (t423 - t478) * t454, t452 * t437 + t499, 0, pkin(1) * t424 + pkin(5) * t386 + t500, -pkin(1) * t421 + pkin(5) * t388 - t508, pkin(1) * t428 + pkin(5) * t425 + t333, pkin(1) * t415 + pkin(5) * t333, t286, -t582, t567, t524, t579, t527, t454 * t265 + t452 * t275 + t573, t454 * t267 + t452 * t290 - t594, t454 * t225 + t452 * t232 + t572, pkin(1) * t339 + pkin(5) * t221 - qJ(3) * t509 + t454 * t236, t286, t567, t582, t527, -t579, t524, t454 * t231 + t452 * t235 + t573, t454 * t216 + t452 * t219 + t572, t454 * t228 + t452 * t234 + t594, pkin(1) * t266 + pkin(5) * t214 + t454 * t210 + t452 * t215;];
tauB_reg = t1;
