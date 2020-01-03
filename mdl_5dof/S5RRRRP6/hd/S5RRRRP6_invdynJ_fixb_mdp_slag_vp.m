% Calculate vector of inverse dynamics joint torques for
% S5RRRRP6
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:44
% EndTime: 2019-12-31 21:54:52
% DurationCPUTime: 4.71s
% Computational Cost: add. (3873->406), mult. (8623->519), div. (0->0), fcn. (5990->10), ass. (0->192)
t430 = cos(qJ(2));
t432 = -pkin(7) - pkin(6);
t395 = t432 * t430;
t386 = qJD(1) * t395;
t426 = sin(qJ(3));
t377 = t426 * t386;
t427 = sin(qJ(2));
t394 = t432 * t427;
t384 = qJD(1) * t394;
t540 = cos(qJ(3));
t344 = t384 * t540 + t377;
t485 = qJD(3) * t540;
t547 = -pkin(2) * t485 + t344;
t423 = qJ(2) + qJ(3);
t418 = sin(t423);
t431 = cos(qJ(1));
t519 = t418 * t431;
t428 = sin(qJ(1));
t520 = t418 * t428;
t548 = g(1) * t519 + g(2) * t520;
t470 = g(1) * t431 + g(2) * t428;
t429 = cos(qJ(4));
t532 = t429 * pkin(4);
t414 = pkin(3) + t532;
t419 = cos(t423);
t424 = -qJ(5) - pkin(8);
t479 = t419 * t414 - t418 * t424;
t516 = t426 * t430;
t383 = t427 * t540 + t516;
t531 = t430 * pkin(2);
t416 = pkin(1) + t531;
t546 = -pkin(8) * t383 - t416;
t486 = qJD(1) * t540;
t503 = qJD(1) * t427;
t374 = t426 * t503 - t430 * t486;
t425 = sin(qJ(4));
t524 = t374 * t425;
t545 = -qJ(5) * t524 + t429 * qJD(5);
t488 = t540 * t430;
t453 = -t426 * t427 + t488;
t496 = qJD(2) + qJD(3);
t347 = t496 * t453;
t437 = t347 * qJD(1);
t436 = t383 * qJDD(1) + t437;
t376 = -qJD(1) * t516 - t427 * t486;
t337 = -pkin(3) * t376 + pkin(8) * t374;
t328 = pkin(2) * t503 + t337;
t544 = t425 * t328 + t547 * t429;
t514 = t429 * t431;
t518 = t425 * t428;
t363 = t419 * t518 + t514;
t515 = t428 * t429;
t517 = t425 * t431;
t365 = -t419 * t517 + t515;
t543 = -g(1) * t365 + g(2) * t363;
t534 = g(3) * t419;
t542 = -t534 + t548;
t447 = t429 * t376 - t425 * t496;
t495 = qJDD(2) + qJDD(3);
t304 = -qJD(4) * t447 + t425 * t436 - t429 * t495;
t541 = t447 ^ 2;
t410 = g(3) * t418;
t533 = g(3) * t425;
t530 = qJD(2) * pkin(2);
t393 = t416 * qJD(1);
t327 = pkin(3) * t374 + pkin(8) * t376 - t393;
t378 = t540 * t386;
t379 = t384 + t530;
t341 = t426 * t379 - t378;
t330 = pkin(8) * t496 + t341;
t305 = t429 * t327 - t330 * t425;
t292 = qJ(5) * t447 + t305;
t370 = qJD(4) + t374;
t288 = pkin(4) * t370 + t292;
t529 = t288 * t429;
t475 = t429 * t496;
t501 = qJD(4) * t425;
t303 = -qJD(4) * t475 - t376 * t501 - t425 * t495 - t429 * t436;
t528 = t303 * t425;
t340 = t379 * t540 + t377;
t329 = -pkin(3) * t496 - t340;
t527 = t329 * t374;
t354 = -t376 * t425 - t475;
t526 = t354 * t370;
t525 = t447 * t370;
t523 = t383 * t425;
t522 = t383 * t429;
t420 = t429 * qJ(5);
t359 = t426 * t394 - t395 * t540;
t351 = t429 * t359;
t413 = pkin(2) * t426 + pkin(8);
t513 = -qJ(5) - t413;
t512 = -t292 + t288;
t511 = t425 * t337 + t429 * t340;
t478 = qJD(4) * t513;
t509 = t425 * t478 - t544 + t545;
t326 = t429 * t328;
t468 = -t376 * pkin(4) + t374 * t420;
t508 = t429 * t478 - t326 - t468 + (-qJD(5) + t547) * t425;
t339 = -pkin(3) * t453 + t546;
t507 = t425 * t339 + t351;
t480 = qJD(4) * t424;
t506 = t425 * t480 - t511 + t545;
t333 = t429 * t337;
t505 = t429 * t480 - t333 - t468 + (-qJD(5) + t340) * t425;
t421 = t427 ^ 2;
t504 = -t430 ^ 2 + t421;
t502 = qJD(3) * t426;
t500 = qJD(4) * t429;
t499 = qJD(1) * qJD(2);
t498 = qJDD(1) * t427;
t497 = qJDD(1) * t430;
t494 = t427 * t530;
t493 = qJD(4) * pkin(8) * t370;
t348 = t496 * t383;
t312 = pkin(3) * t348 - pkin(8) * t347 + t494;
t489 = qJD(2) * t432;
t385 = t427 * t489;
t387 = t430 * t489;
t454 = t394 * t540 + t426 * t395;
t315 = qJD(3) * t454 + t385 * t540 + t426 * t387;
t491 = t425 * t312 + t429 * t315 + t339 * t500;
t490 = t470 * t419 + t410;
t487 = t383 * t500;
t324 = t329 * t500;
t484 = t427 * t499;
t483 = t430 * t499;
t482 = pkin(4) * t425 - t432;
t350 = qJDD(2) * pkin(2) - t432 * (-t483 - t498);
t353 = t432 * (-t484 + t497);
t473 = -t540 * t350 - t426 * t353 + t379 * t502 - t386 * t485;
t299 = -pkin(3) * t495 + t473;
t481 = -t299 - t534;
t477 = t370 * t429;
t441 = t426 * t350 - t353 * t540 + t379 * t485 + t386 * t502;
t298 = pkin(8) * t495 + t441;
t476 = -qJD(4) * t327 - t298;
t415 = -t540 * pkin(2) - pkin(3);
t343 = t426 * t384 - t378;
t472 = t502 * pkin(2) - t343;
t471 = (t501 + t524) * pkin(4);
t469 = g(1) * t428 - g(2) * t431;
t465 = -qJDD(1) * t488 + t426 * t498;
t319 = qJD(1) * t348 + t465;
t408 = pkin(2) * t484;
t291 = t319 * pkin(3) - pkin(8) * t437 + qJDD(1) * t546 + t408;
t290 = t429 * t291;
t467 = -t330 * t500 + t290;
t318 = qJDD(4) + t319;
t466 = -pkin(8) * t318 + t527;
t306 = t327 * t425 + t330 * t429;
t293 = -qJ(5) * t354 + t306;
t463 = -t293 * t425 - t529;
t462 = -t318 * t413 + t527;
t461 = t414 * t418 + t419 * t424;
t460 = -qJ(5) * t347 - qJD(5) * t383;
t459 = t299 * t425 - t306 * t376 + t419 * t533 + t324;
t458 = t305 * t376 + t329 * t501 + t548 * t429;
t457 = t470 * t418;
t456 = -0.2e1 * pkin(1) * t499 - pkin(6) * qJDD(2);
t455 = t416 + t479;
t452 = t347 * t425 + t487;
t451 = t347 * t429 - t383 * t501;
t448 = t425 * t291 + t429 * t298 + t327 * t500 - t330 * t501;
t433 = qJD(2) ^ 2;
t444 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t433 + t469;
t434 = qJD(1) ^ 2;
t443 = pkin(1) * t434 - pkin(6) * qJDD(1) + t470;
t442 = -t393 * t376 - t473 + t542;
t282 = t304 * pkin(4) + qJDD(5) + t299;
t316 = qJD(3) * t359 + t426 * t385 - t387 * t540;
t440 = ((-t303 - t526) * t429 + (-t304 + t525) * t425) * MDP(19) + (-t447 * t477 - t528) * MDP(18) + (-t370 ^ 2 * t425 + t318 * t429 - t354 * t376) * MDP(21) + (t318 * t425 + t370 * t477 - t376 * t447) * MDP(20) + (t374 * t496 + t436) * MDP(13) + (-t465 + (-qJD(1) * t383 - t376) * t496) * MDP(14) + (-t374 ^ 2 + t376 ^ 2) * MDP(12) + t495 * MDP(15) + (-MDP(11) * t374 + t370 * MDP(22)) * t376;
t276 = pkin(4) * t318 + qJ(5) * t303 - qJD(4) * t306 + qJD(5) * t447 - t298 * t425 + t290;
t278 = -qJ(5) * t304 - qJD(5) * t354 + t448;
t439 = qJD(4) * t463 - t276 * t425 + t278 * t429 - t293 * t524 - t374 * t529 - t490;
t438 = -t393 * t374 - t441 + t490;
t392 = pkin(8) * t429 + t420;
t391 = t424 * t425;
t381 = t413 * t429 + t420;
t380 = t513 * t425;
t371 = -qJDD(1) * t416 + t408;
t366 = t419 * t514 + t518;
t364 = -t419 * t515 + t517;
t352 = t354 ^ 2;
t335 = t429 * t339;
t313 = t354 * pkin(4) + qJD(5) + t329;
t309 = t429 * t312;
t307 = -qJ(5) * t523 + t507;
t301 = -pkin(4) * t453 - t359 * t425 - t383 * t420 + t335;
t281 = -qJ(5) * t487 + (-qJD(4) * t359 + t460) * t425 + t491;
t280 = pkin(4) * t348 - t315 * t425 + t309 + t460 * t429 + (-t351 + (qJ(5) * t383 - t339) * t425) * qJD(4);
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t421 + 0.2e1 * t427 * t483) * MDP(4) + 0.2e1 * (t427 * t497 - t499 * t504) * MDP(5) + (qJDD(2) * t427 + t430 * t433) * MDP(6) + (qJDD(2) * t430 - t427 * t433) * MDP(7) + (t427 * t456 + t430 * t444) * MDP(9) + (-t427 * t444 + t430 * t456) * MDP(10) + (-t376 * t347 + t383 * t436) * MDP(11) + (-t383 * t319 - t347 * t374 + t376 * t348 + t436 * t453) * MDP(12) + (t347 * t496 + t383 * t495) * MDP(13) + (-t348 * t496 + t453 * t495) * MDP(14) + (-t316 * t496 - t416 * t319 - t393 * t348 - t371 * t453 + t374 * t494 + t419 * t469 + t454 * t495) * MDP(16) + (-g(1) * t520 + g(2) * t519 - t315 * t496 - t393 * t347 - t359 * t495 + t371 * t383 - t376 * t494 - t416 * t436) * MDP(17) + (-t303 * t522 - t447 * t451) * MDP(18) + ((-t354 * t429 + t425 * t447) * t347 + (t528 - t304 * t429 + (t354 * t425 + t429 * t447) * qJD(4)) * t383) * MDP(19) + (t303 * t453 + t318 * t522 - t348 * t447 + t370 * t451) * MDP(20) + (t304 * t453 - t318 * t523 - t348 * t354 - t370 * t452) * MDP(21) + (-t318 * t453 + t348 * t370) * MDP(22) + ((-t359 * t500 + t309) * t370 + t335 * t318 - t467 * t453 + t305 * t348 + t316 * t354 - t454 * t304 + t383 * t324 - g(1) * t364 - g(2) * t366 + ((-qJD(4) * t339 - t315) * t370 - t359 * t318 - t476 * t453 + t299 * t383 + t329 * t347) * t425) * MDP(23) + (-(-t359 * t501 + t491) * t370 - t507 * t318 + t448 * t453 - t306 * t348 - t316 * t447 + t454 * t303 + t299 * t522 - g(1) * t363 - g(2) * t365 + t451 * t329) * MDP(24) + (t280 * t447 - t281 * t354 + t301 * t303 - t304 * t307 + t469 * t418 + t463 * t347 + (-t276 * t429 - t278 * t425 + (t288 * t425 - t293 * t429) * qJD(4)) * t383) * MDP(25) + (t278 * t307 + t293 * t281 + t276 * t301 + t288 * t280 + t282 * (pkin(4) * t523 - t454) + t313 * (pkin(4) * t452 + t316) - g(1) * (-t428 * t455 + t431 * t482) - g(2) * (t428 * t482 + t431 * t455)) * MDP(26) + t469 * MDP(2) + t470 * MDP(3); (t278 * t381 + t276 * t380 + t282 * (t415 - t532) - g(3) * (t479 + t531) + (t378 + (pkin(2) * qJD(3) - t384) * t426 + t471) * t313 + t509 * t293 + t508 * t288 + t470 * (pkin(2) * t427 + t461)) * MDP(26) + t440 + qJDD(2) * MDP(8) + (-g(3) * t430 + t443 * t427) * MDP(9) + (-t415 * t303 + t462 * t429 - t425 * t457 - t472 * t447 + (t413 * t501 + t544) * t370 + t459) * MDP(24) + (t415 * t304 + t481 * t429 + t462 * t425 + t472 * t354 + (-t413 * t500 + t547 * t425 - t326) * t370 + t458) * MDP(23) + (t344 * t496 + (t376 * t503 - t426 * t495 - t485 * t496) * pkin(2) + t438) * MDP(17) + (t303 * t380 - t304 * t381 - t354 * t509 + t447 * t508 + t439) * MDP(25) + (g(3) * t427 + t430 * t443) * MDP(10) + (t343 * t496 + (-t374 * t503 + t495 * t540 - t496 * t502) * pkin(2) + t442) * MDP(16) + MDP(6) * t498 + MDP(7) * t497 + (-MDP(4) * t427 * t430 + MDP(5) * t504) * t434; (t278 * t392 + t276 * t391 - t282 * t414 - g(3) * t479 + (t471 - t341) * t313 + t506 * t293 + t505 * t288 + t470 * t461) * MDP(26) + t440 + (-pkin(3) * t304 - t333 * t370 - t341 * t354 + (t340 * t370 + t466) * t425 + (t481 - t493) * t429 + t458) * MDP(23) + (t340 * t496 + t438) * MDP(17) + (t303 * t391 - t304 * t392 - t354 * t506 + t447 * t505 + t439) * MDP(25) + (pkin(3) * t303 + t511 * t370 + t341 * t447 + t466 * t429 + (-t457 + t493) * t425 + t459) * MDP(24) + (t341 * t496 + t442) * MDP(16); -t447 * t354 * MDP(18) + (-t352 + t541) * MDP(19) + (-t303 + t526) * MDP(20) + (-t304 - t525) * MDP(21) + t318 * MDP(22) + (t306 * t370 + t329 * t447 + (t476 + t410) * t425 + t467 + t543) * MDP(23) + (g(1) * t366 - g(2) * t364 + t305 * t370 + t329 * t354 + t410 * t429 - t448) * MDP(24) + (pkin(4) * t303 - t354 * t512) * MDP(25) + (t512 * t293 + (t313 * t447 + t418 * t533 + t276 + t543) * pkin(4)) * MDP(26); (-t352 - t541) * MDP(25) + (-t288 * t447 + t293 * t354 + t282 - t542) * MDP(26);];
tau = t1;
