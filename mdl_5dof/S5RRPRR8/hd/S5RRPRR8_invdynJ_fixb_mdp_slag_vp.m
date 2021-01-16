% Calculate vector of inverse dynamics joint torques for
% S5RRPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:22
% EndTime: 2021-01-15 21:35:40
% DurationCPUTime: 6.25s
% Computational Cost: add. (4108->410), mult. (9836->533), div. (0->0), fcn. (7588->14), ass. (0->184)
t446 = sin(pkin(9));
t447 = cos(pkin(9));
t455 = cos(qJ(2));
t504 = qJD(1) * t455;
t490 = t447 * t504;
t451 = sin(qJ(2));
t505 = qJD(1) * t451;
t392 = -t446 * t505 + t490;
t454 = cos(qJ(4));
t377 = t454 * t392;
t403 = t446 * t455 + t447 * t451;
t394 = t403 * qJD(1);
t450 = sin(qJ(4));
t349 = -t394 * t450 + t377;
t453 = cos(qJ(5));
t501 = qJD(5) * t453;
t558 = -t349 * t453 + t501;
t443 = qJ(2) + pkin(9);
t440 = qJ(4) + t443;
t431 = sin(t440);
t452 = sin(qJ(1));
t456 = cos(qJ(1));
t475 = g(1) * t456 + g(2) * t452;
t557 = t475 * t431;
t471 = t392 * t450 + t454 * t394;
t393 = t403 * qJD(2);
t496 = qJDD(1) * t455;
t422 = t447 * t496;
t497 = qJDD(1) * t451;
t358 = qJD(1) * t393 + t446 * t497 - t422;
t498 = qJD(1) * qJD(2);
t489 = t451 * t498;
t420 = t446 * t489;
t488 = t455 * t498;
t359 = t403 * qJDD(1) + t447 * t488 - t420;
t503 = qJD(4) * t450;
t306 = qJD(4) * t377 - t450 * t358 + t454 * t359 - t394 * t503;
t441 = qJDD(2) + qJDD(4);
t442 = qJD(2) + qJD(4);
t449 = sin(qJ(5));
t491 = t453 * t306 + t449 * t441 + t442 * t501;
t502 = qJD(5) * t449;
t294 = -t471 * t502 + t491;
t293 = t294 * t453;
t335 = t442 * t449 + t453 * t471;
t427 = t453 * t441;
t295 = t335 * qJD(5) + t306 * t449 - t427;
t333 = -t453 * t442 + t449 * t471;
t556 = -t449 * t295 - t558 * t333 + t293;
t292 = t294 * t449;
t307 = t471 * qJD(4) + t454 * t358 + t359 * t450;
t305 = qJDD(5) + t307;
t299 = t449 * t305;
t520 = t349 * t442;
t522 = t471 * t442;
t524 = t335 * t471;
t549 = qJD(5) - t349;
t555 = t441 * MDP(19) + (-t307 + t522) * MDP(18) - t349 ^ 2 * MDP(16) + (-MDP(15) * t349 + t471 * MDP(16) - MDP(26) * t549) * t471 + (t306 - t520) * MDP(17) + (t558 * t335 + t292) * MDP(22) + (t558 * t549 + t299 - t524) * MDP(24);
t448 = -qJ(3) - pkin(6);
t418 = t448 * t455;
t410 = qJD(1) * t418;
t397 = t446 * t410;
t417 = t448 * t451;
t409 = qJD(1) * t417;
t528 = qJD(2) * pkin(2);
t401 = t409 + t528;
t356 = t447 * t401 + t397;
t533 = pkin(7) * t394;
t328 = qJD(2) * pkin(3) + t356 - t533;
t517 = t447 * t410;
t357 = t446 * t401 - t517;
t534 = pkin(7) * t392;
t332 = t357 + t534;
t309 = t328 * t454 - t332 * t450;
t301 = -pkin(4) * t442 - t309;
t554 = t301 * t349;
t486 = qJD(2) * t448;
t390 = -qJD(3) * t451 + t455 * t486;
t355 = qJDD(2) * pkin(2) + t390 * qJD(1) + qJDD(1) * t417;
t389 = qJD(3) * t455 + t451 * t486;
t362 = t389 * qJD(1) - qJDD(1) * t418;
t321 = t447 * t355 - t362 * t446;
t308 = qJDD(2) * pkin(3) - pkin(7) * t359 + t321;
t310 = t328 * t450 + t332 * t454;
t322 = t446 * t355 + t447 * t362;
t311 = -pkin(7) * t358 + t322;
t538 = t310 * qJD(4) - t454 * t308 + t450 * t311;
t284 = -pkin(4) * t441 + t538;
t432 = cos(t440);
t530 = g(3) * t432;
t552 = t284 + t530;
t435 = pkin(2) * t455 + pkin(1);
t411 = -t435 * qJD(1) + qJD(3);
t365 = -pkin(3) * t392 + t411;
t424 = g(3) * t431;
t539 = (qJD(4) * t328 + t311) * t454 + t450 * t308 - t332 * t503;
t548 = -t365 * t349 + t432 * t475 + t424 - t539;
t545 = pkin(4) * t471;
t525 = t333 * t471;
t544 = (pkin(8) * t549 + t545) * t549;
t433 = pkin(2) * t447 + pkin(3);
t536 = pkin(2) * t446;
t508 = t450 * t433 + t454 * t536;
t302 = pkin(8) * t442 + t310;
t312 = -pkin(4) * t349 - pkin(8) * t471 + t365;
t287 = -t302 * t449 + t312 * t453;
t543 = -t287 * t471 + t301 * t502 + t453 * t557;
t288 = t302 * t453 + t312 * t449;
t542 = t288 * t471 + t301 * t501 + t449 * t552;
t541 = -t365 * t471 - t530 - t538 + t557;
t339 = -t389 * t446 + t447 * t390;
t518 = t446 * t451;
t402 = -t447 * t455 + t518;
t396 = t402 * qJD(2);
t325 = pkin(7) * t396 + t339;
t340 = t447 * t389 + t446 * t390;
t326 = -pkin(7) * t393 + t340;
t366 = t447 * t417 + t418 * t446;
t341 = -pkin(7) * t403 + t366;
t367 = t446 * t417 - t447 * t418;
t342 = -pkin(7) * t402 + t367;
t472 = t341 * t454 - t342 * t450;
t289 = t472 * qJD(4) + t325 * t450 + t326 * t454;
t360 = t454 * t402 + t403 * t450;
t361 = -t402 * t450 + t403 * t454;
t371 = pkin(3) * t402 - t435;
t316 = pkin(4) * t360 - pkin(8) * t361 + t371;
t318 = t341 * t450 + t342 * t454;
t323 = -t360 * qJD(4) - t393 * t450 - t396 * t454;
t283 = pkin(8) * t441 + t539;
t479 = qJD(5) * t312 + t283;
t537 = t284 * t361 + t301 * t323 - t318 * t305 - (qJD(5) * t316 + t289) * t549 - t479 * t360;
t535 = pkin(2) * t451;
t529 = g(3) * t455;
t527 = t301 * t361;
t526 = t316 * t305;
t523 = t335 * t449;
t516 = t449 * t452;
t515 = t449 * t456;
t514 = t452 * t453;
t300 = t453 * t305;
t513 = t453 * t456;
t363 = -t409 * t446 + t517;
t336 = t363 - t534;
t364 = t447 * t409 + t397;
t337 = t364 - t533;
t469 = t433 * t454 - t450 * t536;
t510 = -t469 * qJD(4) + t336 * t450 + t337 * t454;
t509 = qJD(4) * t508 + t336 * t454 - t337 * t450;
t444 = t451 ^ 2;
t507 = -t455 ^ 2 + t444;
t495 = pkin(2) * t489 + qJDD(3);
t437 = t451 * t528;
t369 = pkin(3) * t393 + t437;
t368 = pkin(2) * t505 + pkin(3) * t394;
t481 = t549 * t449;
t386 = -t435 * qJDD(1) + t495;
t331 = pkin(3) * t358 + t386;
t286 = pkin(4) * t307 - pkin(8) * t306 + t331;
t478 = qJD(5) * t302 - t286;
t388 = pkin(8) + t508;
t476 = -pkin(8) * t349 + qJD(5) * t388 + t368 + t545;
t474 = g(1) * t452 - g(2) * t456;
t473 = -t305 * t388 - t554;
t470 = t300 - (-t349 * t449 + t502) * t549;
t467 = -0.2e1 * pkin(1) * t498 - pkin(6) * qJDD(2);
t466 = t323 * t453 - t361 * t502;
t465 = -pkin(8) * t305 + t309 * t549 - t554;
t457 = qJD(2) ^ 2;
t463 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t457 + t474;
t458 = qJD(1) ^ 2;
t462 = pkin(1) * t458 - pkin(6) * qJDD(1) + t475;
t439 = cos(t443);
t438 = sin(t443);
t387 = -pkin(4) - t469;
t385 = t432 * t513 + t516;
t384 = -t432 * t515 + t514;
t383 = -t432 * t514 + t515;
t382 = t432 * t516 + t513;
t324 = t361 * qJD(4) + t454 * t393 - t396 * t450;
t296 = pkin(4) * t324 - pkin(8) * t323 + t369;
t290 = t318 * qJD(4) - t325 * t454 + t326 * t450;
t285 = t453 * t286;
t1 = [(-g(1) * t383 - g(2) * t385 + t285 * t360 + t287 * t324 + t290 * t333 - t472 * t295 + (t296 * t549 + t526 + (-t302 * t360 - t318 * t549 + t527) * qJD(5)) * t453 + t537 * t449) * MDP(27) + (t294 * t360 + t361 * t300 + t324 * t335 + t466 * t549) * MDP(24) + (-t361 * t299 - t295 * t360 - t324 * t333 + (-t323 * t449 - t361 * t501) * t549) * MDP(25) + (t305 * t360 + t324 * t549) * MDP(26) + (-g(1) * t382 - g(2) * t384 - t288 * t324 + t290 * t335 - t472 * t294 + (-(-qJD(5) * t318 + t296) * t549 - t526 + t478 * t360 - qJD(5) * t527) * t449 + t537 * t453) * MDP(28) + (-t324 * t442 - t360 * t441) * MDP(18) + (t323 * t442 + t361 * t441) * MDP(17) + (qJDD(2) * t451 + t455 * t457) * MDP(6) + (qJDD(2) * t455 - t451 * t457) * MDP(7) + (-t290 * t442 + t307 * t371 + t324 * t365 + t331 * t360 - t349 * t369 + t474 * t432 + t441 * t472) * MDP(20) + (-t306 * t360 - t307 * t361 + t323 * t349 - t324 * t471) * MDP(16) + (t306 * t361 + t323 * t471) * MDP(15) + (-t289 * t442 + t306 * t371 - t318 * t441 + t323 * t365 + t331 * t361 + t369 * t471 - t474 * t431) * MDP(21) + (-qJDD(2) * t367 - t359 * t435 + t386 * t403 - t396 * t411 - t474 * t438 + (t394 * t535 - t340) * qJD(2)) * MDP(12) + (-t321 * t403 - t322 * t402 - t339 * t394 + t340 * t392 + t356 * t396 - t357 * t393 - t358 * t367 - t359 * t366 - t475) * MDP(13) + (qJDD(2) * t366 - t358 * t435 + t386 * t402 + t393 * t411 + t474 * t439 + (-t392 * t535 + t339) * qJD(2)) * MDP(11) + qJDD(1) * MDP(1) + (t361 * t293 + t466 * t335) * MDP(22) + ((-t333 * t453 - t523) * t323 + (-t292 - t295 * t453 + (t333 * t449 - t335 * t453) * qJD(5)) * t361) * MDP(23) + 0.2e1 * (t451 * t496 - t507 * t498) * MDP(5) + (t322 * t367 + t357 * t340 + t321 * t366 + t356 * t339 - t386 * t435 + t411 * t437 - g(1) * (-t435 * t452 - t448 * t456) - g(2) * (t435 * t456 - t448 * t452)) * MDP(14) + (qJDD(1) * t444 + 0.2e1 * t451 * t488) * MDP(4) + (t467 * t451 + t463 * t455) * MDP(9) + (-t463 * t451 + t467 * t455) * MDP(10) + t474 * MDP(2) + t475 * MDP(3); (-t368 * t471 - t508 * t441 + t510 * t442 + t548) * MDP(21) + t555 + ((t357 + t363) * t394 + (t356 - t364) * t392 + (-t358 * t446 - t359 * t447) * pkin(2)) * MDP(13) + (t387 * t295 - t552 * t453 + t473 * t449 + t509 * t333 + (t510 * t449 - t476 * t453) * t549 + t543) * MDP(27) + (t387 * t294 + t473 * t453 - t449 * t557 + t509 * t335 + (t476 * t449 + t510 * t453) * t549 + t542) * MDP(28) + (-t523 * t549 + t556) * MDP(23) + (t349 * t368 + t469 * t441 - t509 * t442 + t541) * MDP(20) + (g(3) * t451 + t462 * t455) * MDP(10) + MDP(7) * t496 + MDP(6) * t497 + (-t451 * t455 * MDP(4) + t507 * MDP(5)) * t458 + (g(3) * t438 + qJD(2) * t364 - t392 * t411 + t475 * t439 + (-qJDD(2) * t446 - t394 * t505) * pkin(2) - t322) * MDP(12) + (-g(3) * t439 - qJD(2) * t363 - t394 * t411 + t475 * t438 + (qJDD(2) * t447 + t392 * t505) * pkin(2) + t321) * MDP(11) + (t462 * t451 - t529) * MDP(9) + (-t356 * t363 - t357 * t364 + (-t529 + t321 * t447 + t322 * t446 + (-qJD(1) * t411 + t475) * t451) * pkin(2)) * MDP(14) + qJDD(2) * MDP(8) + (t470 + t525) * MDP(25); -t422 * MDP(11) - t420 * MDP(12) + (-t392 ^ 2 - t394 ^ 2) * MDP(13) + (t356 * t394 - t357 * t392 - t474 + t495) * MDP(14) + (t307 + t522) * MDP(20) + (t306 + t520) * MDP(21) + (t470 - t525) * MDP(27) + (-t453 * t549 ^ 2 - t299 - t524) * MDP(28) + (MDP(11) * t518 + t403 * MDP(12) - t435 * MDP(14)) * qJDD(1) + ((t446 * t504 + t447 * t505 + t394) * MDP(11) + (t392 + t490) * MDP(12)) * qJD(2); (t310 * t442 + t541) * MDP(20) + (t309 * t442 + t548) * MDP(21) + (-t335 * t481 + t556) * MDP(23) + (-t481 * t549 + t300 + t525) * MDP(25) + (-pkin(4) * t295 - t310 * t333 + t465 * t449 + (-t552 - t544) * t453 + t543) * MDP(27) + (-pkin(4) * t294 - t310 * t335 + t465 * t453 + (-t557 + t544) * t449 + t542) * MDP(28) + t555; t335 * t333 * MDP(22) + (-t333 ^ 2 + t335 ^ 2) * MDP(23) + (t333 * t549 + t491) * MDP(24) + (t335 * t549 + t427) * MDP(25) + t305 * MDP(26) + (-g(1) * t384 + g(2) * t382 + t288 * t549 - t301 * t335 + t285) * MDP(27) + (g(1) * t385 - g(2) * t383 + t287 * t549 + t301 * t333) * MDP(28) + ((-t283 + t424) * MDP(28) + (-MDP(25) * t471 - MDP(27) * t302 - MDP(28) * t312) * qJD(5)) * t453 + (-qJD(5) * t471 * MDP(24) + (-qJD(5) * t442 - t306) * MDP(25) + (-t479 + t424) * MDP(27) + t478 * MDP(28)) * t449;];
tau = t1;
