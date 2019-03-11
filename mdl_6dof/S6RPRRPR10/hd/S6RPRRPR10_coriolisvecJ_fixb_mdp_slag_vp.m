% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:20
% EndTime: 2019-03-09 05:37:32
% DurationCPUTime: 6.25s
% Computational Cost: add. (2969->459), mult. (6207->624), div. (0->0), fcn. (3769->6), ass. (0->201)
t409 = sin(qJ(3));
t501 = qJD(1) * t409;
t398 = qJD(4) + t501;
t414 = -pkin(1) - pkin(7);
t408 = sin(qJ(4));
t411 = cos(qJ(4));
t533 = qJ(5) * t411;
t538 = pkin(4) + pkin(5);
t432 = t408 * t538 - t533;
t557 = t414 - t432;
t484 = t411 * qJD(3);
t412 = cos(qJ(3));
t500 = qJD(1) * t412;
t369 = t408 * t500 - t484;
t471 = t411 * t500;
t498 = qJD(3) * t408;
t371 = t471 + t498;
t407 = sin(qJ(6));
t410 = cos(qJ(6));
t322 = t369 * t407 + t371 * t410;
t439 = -t410 * t369 + t371 * t407;
t480 = qJD(1) * qJD(3);
t464 = t412 * t480;
t556 = t439 * t322 * MDP(25) + (t322 ^ 2 - t439 ^ 2) * MDP(26) - MDP(29) * t464;
t481 = qJD(6) - t398;
t555 = t322 * t481;
t451 = pkin(3) * t412 + pkin(8) * t409;
t368 = qJD(3) * t451 + qJD(2);
t493 = qJD(4) * t409;
t554 = -t414 * t493 + t368;
t545 = qJD(4) - qJD(6);
t394 = qJD(1) * t414 + qJD(2);
t521 = t394 * t412;
t361 = -qJD(3) * pkin(3) - t521;
t423 = qJ(5) * t371 - t361;
t299 = -t369 * t538 + t423;
t386 = t398 * qJD(5);
t395 = qJ(5) * t464;
t348 = t368 * qJD(1);
t379 = pkin(3) * t409 - pkin(8) * t412 + qJ(2);
t358 = t379 * qJD(1);
t378 = t409 * t394;
t360 = qJD(3) * pkin(8) + t378;
t513 = t411 * t412;
t366 = t394 * t513;
t492 = qJD(4) * t411;
t494 = qJD(4) * t408;
t428 = -qJD(3) * t366 - t408 * t348 - t358 * t492 + t360 * t494;
t285 = t386 + t395 - t428;
t497 = qJD(3) * t409;
t465 = t408 * t497;
t387 = qJD(1) * t465;
t495 = qJD(4) * t371;
t334 = -t387 + t495;
t282 = pkin(9) * t334 + t285;
t491 = qJD(4) * t412;
t468 = t408 * t491;
t424 = t484 * t409 + t468;
t333 = qJD(1) * t424 - qJD(4) * t484;
t496 = qJD(3) * t412;
t470 = t408 * t496;
t453 = -t411 * t348 + t358 * t494 + t360 * t492 + t394 * t470;
t283 = pkin(9) * t333 - t464 * t538 + t453;
t462 = t407 * t282 - t410 * t283;
t552 = t299 * t322 + t462;
t479 = 0.2e1 * qJD(1);
t406 = t412 ^ 2;
t549 = MDP(8) * (t409 ^ 2 - t406);
t548 = t481 * t439;
t454 = pkin(4) * t464;
t288 = t453 - t454;
t313 = t408 * t358 + t411 * t360;
t389 = t398 * qJ(5);
t306 = t389 + t313;
t547 = -t306 * t398 + t288;
t546 = qJD(5) * t408 + t378;
t312 = t411 * t358 - t408 * t360;
t482 = qJD(5) - t312;
t483 = -pkin(9) * t371 + t482;
t293 = -t398 * t538 + t483;
t487 = qJD(6) * t410;
t473 = t410 * t282 + t407 * t283 + t293 * t487;
t544 = -t299 * t439 + t473;
t534 = qJ(5) * t408;
t542 = -t411 * t538 - t534;
t539 = t371 ^ 2;
t537 = pkin(8) - pkin(9);
t416 = qJD(1) ^ 2;
t536 = qJ(2) * t416;
t535 = qJ(5) * t369;
t421 = -qJ(5) * t333 + qJD(5) * t371 - t394 * t497;
t289 = pkin(4) * t334 - t421;
t532 = t289 * t408;
t531 = t289 * t411;
t309 = pkin(4) * t369 - t423;
t529 = t309 * t371;
t528 = t333 * t408;
t527 = t333 * t409;
t525 = t334 * t409;
t524 = t369 * t371;
t523 = t369 * t398;
t522 = t371 * t398;
t520 = t398 * t408;
t519 = t398 * t409;
t518 = t398 * t411;
t302 = pkin(9) * t369 + t313;
t296 = t302 + t389;
t517 = t407 * t296;
t516 = t408 * t412;
t515 = t409 * t411;
t514 = t409 * t414;
t512 = t412 * t416;
t415 = qJD(3) ^ 2;
t511 = t414 * t415;
t373 = t407 * t408 + t410 * t411;
t323 = t545 * t373;
t426 = t373 * qJD(1);
t510 = t409 * t426 + t323;
t438 = t407 * t411 - t408 * t410;
t429 = t438 * t409;
t488 = qJD(6) * t407;
t509 = qJD(1) * t429 + t407 * t492 + t408 * t487 - t410 * t494 - t411 * t488;
t446 = pkin(4) * t408 - t533;
t508 = -t398 * t446 + t546;
t507 = -t398 * t432 + t546;
t375 = t451 * qJD(1);
t506 = t408 * t375 + t366;
t505 = t408 * t379 + t411 * t514;
t350 = t373 * t412;
t502 = qJD(1) * t350;
t499 = qJD(3) * t481;
t489 = qJD(5) * t411;
t485 = t361 * qJD(4);
t478 = pkin(8) * t520;
t477 = pkin(8) * t518;
t476 = pkin(8) * t496;
t385 = t537 * t411;
t472 = -t410 * t333 + t407 * t334 + t369 * t487;
t315 = qJ(5) * t500 + t506;
t331 = t409 * qJ(5) + t505;
t469 = t414 * t496;
t466 = t411 * t491;
t463 = MDP(18) * t500;
t461 = -t333 * t407 - t410 * t334;
t460 = -t361 + t521;
t459 = t375 * t411 - t394 * t516;
t391 = t408 * t514;
t458 = t379 * t411 - t391;
t457 = t369 + t484;
t456 = -t371 + t498;
t455 = qJD(1) + t493;
t388 = t409 * t464;
t425 = pkin(9) * t515 - t412 * t538;
t449 = qJD(1) * t425 - t385 * t545 - t459;
t384 = t537 * t408;
t448 = pkin(9) * t408 * t501 + qJD(6) * t384 - t537 * t494 - t315;
t447 = pkin(4) * t411 + t534;
t445 = qJ(5) * t410 - t407 * t538;
t444 = -qJ(5) * t407 - t410 * t538;
t278 = t407 * t293 + t410 * t296;
t305 = -pkin(4) * t398 + t482;
t443 = t305 * t411 - t306 * t408;
t442 = t305 * t408 + t306 * t411;
t311 = t391 + (-pkin(9) * t412 - t379) * t411 - t538 * t409;
t317 = pkin(9) * t516 + t331;
t441 = t311 * t410 - t317 * t407;
t440 = t311 * t407 + t317 * t410;
t437 = qJD(1) * t406 - t519;
t436 = -t379 * t494 - t408 * t469 + t411 * t554;
t435 = -t414 + t446;
t434 = -t309 * t409 + t476;
t433 = t361 * t409 - t476;
t430 = t313 * t398 - t453;
t286 = -t371 * t488 + t472;
t427 = t438 * qJD(1);
t422 = t379 * t492 + t408 * t554 + t411 * t469;
t420 = t312 * t398 + t428;
t295 = qJ(5) * t496 + t409 * qJD(5) + t422;
t287 = qJD(6) * t322 + t461;
t419 = qJD(4) * t443 + t285 * t411 + t288 * t408;
t418 = t545 * t438;
t380 = -pkin(3) - t447;
t364 = pkin(3) - t542;
t349 = t407 * t513 - t410 * t516;
t340 = t435 * t412;
t332 = -pkin(4) * t409 - t458;
t330 = t557 * t412;
t325 = pkin(4) * t371 + t535;
t316 = -pkin(4) * t500 - t459;
t310 = -t371 * t538 - t535;
t307 = -t333 + t523;
t304 = (qJD(4) * t447 - t489) * t412 - t435 * t497;
t300 = -pkin(4) * t496 - t436;
t298 = -t373 * t497 + t412 * t418;
t297 = qJD(3) * t429 + t323 * t412;
t294 = (t542 * qJD(4) + t489) * t412 - t557 * t497;
t291 = (-t465 + t466) * pkin(9) + t295;
t290 = pkin(9) * t468 + qJD(3) * t425 - t436;
t284 = -t334 * t538 + t421;
t277 = t293 * t410 - t517;
t1 = [(-t286 * t409 + t298 * t481 + (-t322 - t502) * t496) * MDP(27) + (-t422 * t398 + t428 * t409 + (t414 * t333 - t408 * t485) * t412 + ((-qJD(1) * t505 - t313) * t412 + (t414 * t371 + t460 * t411) * t409) * qJD(3)) * MDP(20) + (-t295 * t369 + t300 * t371 - t331 * t334 - t332 * t333 - t443 * t497 + (-qJD(4) * t442 - t285 * t408 + t288 * t411) * t412) * MDP(22) + (-t286 * t349 - t287 * t350 + t297 * t322 - t298 * t439) * MDP(26) + (t286 * t350 + t298 * t322) * MDP(25) + (t285 * t331 + t288 * t332 + t289 * t340 + t295 * t306 + t300 * t305 + t304 * t309) * MDP(24) + (t287 * t409 + t297 * t481 + (qJD(1) * t349 + t439) * t496) * MDP(28) + (-(qJD(6) * t441 + t290 * t407 + t291 * t410) * t481 + (-t296 * t488 + t473) * t409 + t294 * t322 + t330 * t286 + t284 * t350 + t299 * t298 + (qJD(1) * t440 + t278) * t496) * MDP(31) + ((t290 * t410 - t291 * t407) * t481 + t462 * t409 + t294 * t439 + t330 * t287 + t284 * t349 - t299 * t297 + (t278 * t409 - t440 * t481) * qJD(6) + (-t441 * qJD(1) - t277) * t496) * MDP(30) + (t398 * t496 + t388) * MDP(18) + (-t481 * t496 + t388) * MDP(29) + (t436 * t398 - t453 * t409 + (-t414 * t334 + t411 * t485) * t412 + ((qJD(1) * t458 + t312) * t412 + (t414 * t369 + t408 * t460) * t409) * qJD(3)) * MDP(19) - 0.2e1 * MDP(7) * t388 + (-t333 * t513 - t371 * t424) * MDP(14) + (-t412 * t511 + (-qJ(2) * t497 + qJD(2) * t412) * t479) * MDP(13) + (-t409 * t511 + (qJ(2) * t496 + qJD(2) * t409) * t479) * MDP(12) + 0.2e1 * t480 * t549 + (-t398 * t466 - t525 + (-t369 * t412 - t408 * t437) * qJD(3)) * MDP(17) + (-t398 * t468 - t527 + (t371 * t412 + t411 * t437) * qJD(3)) * MDP(16) + ((t369 * t411 + t371 * t408) * t497 + (t528 - t334 * t411 + (t369 * t408 - t371 * t411) * qJD(4)) * t412) * MDP(15) + (-t300 * t398 + t304 * t369 + t334 * t340 + (-t309 * t498 - t288) * t409 + (t309 * t492 + t532 + (-qJD(1) * t332 - t305) * qJD(3)) * t412) * MDP(21) + (t295 * t398 - t304 * t371 + t333 * t340 + (t309 * t484 + t285) * t409 + (t309 * t494 - t531 + (qJD(1) * t331 + t306) * qJD(3)) * t412) * MDP(23) + (MDP(6) * qJ(2) + MDP(5)) * qJD(2) * t479 + (-MDP(10) * t412 - MDP(9) * t409) * t415; -t416 * MDP(5) - MDP(6) * t536 + ((-t369 * t496 + t371 * t455 - t525) * t411 + (t369 * t455 + t371 * t496 - t527) * t408) * MDP(22) + (t443 * qJD(1) + (qJD(3) * t442 - t289) * t412 + (qJD(3) * t309 + t419) * t409) * MDP(24) + (t481 * t426 + (-t438 * t499 + t287) * t412 + (t323 * t481 + (t412 * t427 - t439) * qJD(3)) * t409) * MDP(30) + (-t481 * t427 + (-t373 * t499 + t286) * t412 + (-t418 * t481 + (-t322 + t502) * qJD(3)) * t409) * MDP(31) + (t409 * MDP(12) + t412 * MDP(13)) * (-t415 - t416) + (MDP(19) + MDP(21)) * (t369 * t497 + (-t334 - t387) * t412 + (-t411 * t455 - t470) * t398) + (-MDP(20) + MDP(23)) * ((t409 * (-t371 + t471) + t398 * t513) * qJD(3) - t333 * t412 - t455 * t520); -t416 * t549 - qJ(2) * MDP(12) * t512 + (t371 * t518 - t528) * MDP(14) + ((-t333 - t523) * t411 + (-t334 - t522) * t408) * MDP(15) + (t398 * t492 + (t398 * t515 + t412 * t456) * qJD(1)) * MDP(16) + (-t398 * t494 + (-t408 * t519 + t412 * t457) * qJD(1)) * MDP(17) - t398 * t463 + (-pkin(3) * t334 - t459 * t398 - t457 * t378 + (t361 * t408 - t477) * qJD(4) + (-t312 * t412 + t408 * t433) * qJD(1)) * MDP(19) + (pkin(3) * t333 + t506 * t398 + t456 * t378 + (t361 * t411 + t478) * qJD(4) + (t313 * t412 + t411 * t433) * qJD(1)) * MDP(20) + (-t531 + t316 * t398 + t334 * t380 - t508 * t369 + (t309 * t408 - t477) * qJD(4) + (t305 * t412 - t408 * t434) * qJD(1)) * MDP(21) + (t315 * t369 - t316 * t371 + (t285 + t398 * t305 + (-t334 + t495) * pkin(8)) * t411 + ((qJD(4) * t369 - t333) * pkin(8) + t547) * t408) * MDP(22) + (-t532 - t315 * t398 + t333 * t380 + t508 * t371 + (-t309 * t411 - t478) * qJD(4) + (-t306 * t412 + t411 * t434) * qJD(1)) * MDP(23) + (pkin(8) * t419 + t289 * t380 - t305 * t316 - t306 * t315 - t309 * t508) * MDP(24) + (-t286 * t438 + t322 * t510) * MDP(25) + (-t286 * t373 + t287 * t438 - t322 * t509 - t439 * t510) * MDP(26) + (t510 * t481 + (qJD(3) * t438 + t322) * t500) * MDP(27) + (-t509 * t481 + (qJD(3) * t373 - t439) * t500) * MDP(28) + t481 * MDP(29) * t500 + (t284 * t373 + t364 * t287 - (t407 * t448 + t410 * t449) * t481 + t507 * t439 + t509 * t299 + (-(t384 * t410 - t385 * t407) * qJD(3) + t277) * t500) * MDP(30) + (-t284 * t438 + t364 * t286 - (-t407 * t449 + t410 * t448) * t481 + t507 * t322 + t510 * t299 + ((t384 * t407 + t385 * t410) * qJD(3) - t278) * t500) * MDP(31) + (MDP(13) * t536 + MDP(7) * t512) * t409; MDP(14) * t524 + (-t369 ^ 2 + t539) * MDP(15) + t307 * MDP(16) + (t522 - t334) * MDP(17) + qJD(3) * t463 + (-t361 * t371 + t430) * MDP(19) + (t361 * t369 + t420) * MDP(20) + (-t325 * t369 + t430 + 0.2e1 * t454 - t529) * MDP(21) + (pkin(4) * t333 - qJ(5) * t334 + (t306 - t313) * t371 + (t305 - t482) * t369) * MDP(22) + (-t309 * t369 + t325 * t371 + 0.2e1 * t386 + 0.2e1 * t395 - t420) * MDP(23) + (-pkin(4) * t288 + qJ(5) * t285 - t305 * t313 + t306 * t482 - t309 * t325) * MDP(24) + (-t286 - t548) * MDP(27) + (t287 - t555) * MDP(28) + (-t444 * t464 - t310 * t439 - (t410 * t302 + t407 * t483) * t481 + (-t445 * t481 + t278) * qJD(6) + t552) * MDP(30) + (t445 * t464 - t310 * t322 - (-t407 * t302 + t410 * t483) * t481 + (-t444 * t481 - t517) * qJD(6) + t544) * MDP(31) - t556; (-t464 + t524) * MDP(21) + t307 * MDP(22) + (-t398 ^ 2 - t539) * MDP(23) + (t529 + t547) * MDP(24) + (-t371 * t439 - t410 * t464) * MDP(30) + (-t322 * t371 + t407 * t464) * MDP(31) - (MDP(30) * t407 + MDP(31) * t410) * t481 ^ 2; (t472 + t548) * MDP(27) + (-t461 + t555) * MDP(28) + (t278 * t481 - t552) * MDP(30) + (t277 * t481 - t544) * MDP(31) + ((-MDP(28) * t371 - MDP(30) * t296) * t410 + (-MDP(27) * t371 - MDP(28) * t369 - MDP(30) * t293 + MDP(31) * t296) * t407) * qJD(6) + t556;];
tauc  = t1;
