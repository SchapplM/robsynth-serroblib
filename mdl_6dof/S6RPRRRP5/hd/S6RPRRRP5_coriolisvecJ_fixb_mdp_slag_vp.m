% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:34
% EndTime: 2019-03-09 06:12:46
% DurationCPUTime: 6.79s
% Computational Cost: add. (9590->460), mult. (25132->569), div. (0->0), fcn. (19771->8), ass. (0->193)
t440 = cos(qJ(5));
t504 = qJD(5) * t440;
t436 = cos(pkin(10));
t441 = cos(qJ(3));
t521 = t436 * t441;
t435 = sin(pkin(10));
t439 = sin(qJ(3));
t522 = t435 * t439;
t465 = -t521 + t522;
t409 = t465 * qJD(1);
t416 = t435 * t441 + t436 * t439;
t410 = t416 * qJD(1);
t438 = sin(qJ(4));
t546 = cos(qJ(4));
t460 = -t546 * t409 - t438 * t410;
t570 = t460 * t440;
t582 = t504 - t570;
t495 = qJD(1) * qJD(3);
t483 = t441 * t495;
t425 = t436 * t483;
t484 = t439 * t495;
t405 = -t435 * t484 + t425;
t510 = t435 * t483 + t436 * t484;
t479 = t438 * t405 + t546 * t510;
t459 = -t438 * t409 + t410 * t546;
t554 = t459 * qJD(4);
t359 = t479 + t554;
t358 = t440 * t359;
t498 = qJD(5) - t460;
t437 = sin(qJ(5));
t576 = t498 * t437;
t581 = -t498 * t576 + t358;
t450 = t405 * t546 - t438 * t510;
t572 = t460 * qJD(4);
t444 = t450 + t572;
t494 = qJD(3) + qJD(4);
t476 = qJD(5) * t494;
t505 = qJD(5) * t437;
t331 = t459 * t505 + (-t444 - t476) * t440;
t329 = t331 * t437;
t330 = t331 * t440;
t384 = t437 * t494 + t440 * t459;
t356 = t437 * t359;
t512 = t498 * t504 + t356;
t462 = -t498 * t570 + t512;
t443 = t437 * t444;
t506 = qJD(5) * t384;
t332 = t443 + t506;
t382 = t437 * t459 - t440 * t494;
t463 = -t437 * t332 - t382 * t582;
t501 = t460 * qJD(3);
t502 = t459 * qJD(3);
t531 = t384 * t459;
t580 = (t462 - t531) * MDP(24) + (t502 - t479) * MDP(18) - t460 ^ 2 * MDP(16) + (-MDP(15) * t460 + MDP(16) * t459 - MDP(26) * t498) * t459 + (-t501 + t450) * MDP(17) + (t384 * t582 - t329) * MDP(22) + (-t384 * t576 - t330 + t463) * MDP(23);
t579 = t498 ^ 2;
t543 = pkin(7) + qJ(2);
t421 = t543 * t435;
t417 = qJD(1) * t421;
t422 = t543 * t436;
t418 = qJD(1) * t422;
t467 = t417 * t439 - t418 * t441;
t381 = -pkin(8) * t409 - t467;
t378 = t546 * t381;
t559 = -t441 * t417 - t418 * t439;
t380 = -pkin(8) * t410 + t559;
t379 = qJD(3) * pkin(3) + t380;
t343 = t438 * t379 + t378;
t337 = pkin(9) * t494 + t343;
t429 = -pkin(2) * t436 - pkin(1);
t419 = qJD(1) * t429 + qJD(2);
t398 = pkin(3) * t409 + t419;
t344 = -pkin(4) * t460 - pkin(9) * t459 + t398;
t311 = t337 * t440 + t344 * t437;
t300 = qJ(6) * t498 + t311;
t578 = t300 * t498;
t377 = t438 * t381;
t342 = t379 * t546 - t377;
t336 = -pkin(4) * t494 - t342;
t319 = t382 * pkin(5) - t384 * qJ(6) + t336;
t577 = t319 * t498;
t362 = pkin(4) * t459 - pkin(9) * t460;
t367 = -t510 * pkin(8) - qJD(2) * t409 + qJD(3) * t559;
t452 = t416 * qJD(2);
t449 = qJD(1) * t452;
t368 = -pkin(8) * t405 + qJD(3) * t467 - t449;
t486 = qJD(4) * t546;
t507 = qJD(4) * t438;
t447 = -t367 * t546 - t438 * t368 - t379 * t486 + t381 * t507;
t569 = -t398 * t460 + t447;
t471 = pkin(5) * t437 - qJ(6) * t440;
t568 = -pkin(5) * t505 + qJ(6) * t504 + qJD(6) * t437 + t460 * t471;
t475 = pkin(3) * t486;
t563 = pkin(5) * t459;
t562 = qJ(6) * t459;
t534 = t382 * t459;
t485 = t510 * pkin(3);
t323 = t359 * pkin(4) - pkin(9) * t444 + t485;
t454 = t437 * t323 - t337 * t505 + t344 * t504 - t440 * t447;
t541 = qJ(6) * t359;
t291 = qJD(6) * t498 + t454 + t541;
t474 = -t440 * t323 + t337 * t504 + t344 * t505 - t437 * t447;
t545 = pkin(5) * t359;
t293 = t474 - t545;
t560 = t291 * t440 + t293 * t437;
t345 = t438 * t380 + t378;
t473 = pkin(3) * t507 - t345;
t523 = t421 * t441;
t388 = -pkin(8) * t416 - t422 * t439 - t523;
t466 = t421 * t439 - t422 * t441;
t389 = -pkin(8) * t465 - t466;
t354 = t438 * t388 + t389 * t546;
t397 = t416 * t546 - t438 * t465;
t401 = pkin(3) * t465 + t429;
t458 = -t438 * t416 - t465 * t546;
t355 = -pkin(4) * t458 - pkin(9) * t397 + t401;
t513 = t440 * t354 + t437 * t355;
t558 = t498 * t505 - t358;
t557 = t546 * t388 - t438 * t389;
t310 = -t337 * t437 + t344 * t440;
t497 = qJD(6) - t310;
t299 = -pkin(5) * t498 + t497;
t555 = t299 * t437 + t300 * t440;
t553 = (t435 ^ 2 + t436 ^ 2) * (qJ(2) * MDP(7) + MDP(6));
t305 = t438 * t367 - t546 * t368 + t379 * t507 + t381 * t486;
t296 = pkin(5) * t332 + qJ(6) * t331 - qJD(6) * t384 + t305;
t552 = -t296 * t440 + t299 * t459 + t319 * t505;
t540 = t296 * t437;
t551 = -t300 * t459 - t540;
t550 = -t398 * t459 - t305;
t549 = -t305 * t440 - t310 * t459 + t336 * t505;
t548 = t305 * t437 + t311 * t459 + t336 * t504;
t547 = t384 ^ 2;
t412 = t416 * qJD(3);
t544 = t412 * pkin(3);
t537 = t311 * t498;
t536 = t332 * t440;
t430 = pkin(3) * t438 + pkin(9);
t535 = t359 * t430;
t533 = t382 * t437;
t532 = t384 * t382;
t530 = t384 * t437;
t529 = t384 * t440;
t528 = t460 * t437;
t525 = t397 * t440;
t518 = t343 + t568;
t517 = -t473 + t568;
t515 = t440 * t342 + t437 * t362;
t346 = t380 * t546 - t377;
t349 = pkin(3) * t410 + t362;
t514 = t440 * t346 + t437 * t349;
t508 = qJD(3) * t410;
t496 = qJD(1) * qJD(2);
t493 = t546 * pkin(3);
t491 = t299 * t504 + t560;
t489 = qJD(1) * t522;
t488 = t430 * t505;
t472 = t440 * pkin(5) + t437 * qJ(6);
t470 = t299 * t440 - t300 * t437;
t469 = -t336 * t460 - t535;
t468 = -t342 * t437 + t362 * t440;
t420 = -pkin(4) - t472;
t461 = t460 * t576 - t558;
t411 = t465 * qJD(3);
t365 = qJD(4) * t458 - t411 * t546 - t438 * t412;
t457 = t365 * t437 + t397 * t504;
t456 = -t365 * t440 + t397 * t505;
t455 = t319 * t384 + t474;
t448 = -qJD(3) * t523 + qJD(2) * t521 + (-qJD(2) * t435 - qJD(3) * t422) * t439;
t369 = -pkin(8) * t412 + t448;
t445 = qJD(3) * t466 - t452;
t370 = pkin(8) * t411 + t445;
t314 = qJD(4) * t557 + t546 * t369 + t438 * t370;
t366 = qJD(4) * t397 - t438 * t411 + t412 * t546;
t326 = pkin(4) * t366 - pkin(9) * t365 + t544;
t453 = t440 * t314 + t437 * t326 - t354 * t505 + t355 * t504;
t451 = t512 * pkin(9);
t446 = qJD(5) * t470 + t560;
t315 = qJD(4) * t354 + t438 * t369 - t370 * t546;
t431 = -t493 - pkin(4);
t414 = -t493 + t420;
t347 = pkin(5) * t384 + qJ(6) * t382;
t327 = t397 * t471 - t557;
t317 = pkin(5) * t458 + t354 * t437 - t355 * t440;
t316 = -qJ(6) * t458 + t513;
t312 = t382 * t498 - t331;
t309 = -t468 - t563;
t308 = t515 + t562;
t307 = t346 * t437 - t349 * t440 - t563;
t306 = t514 + t562;
t297 = t471 * t365 + (qJD(5) * t472 - qJD(6) * t440) * t397 + t315;
t295 = -pkin(5) * t366 + qJD(5) * t513 + t314 * t437 - t326 * t440;
t294 = qJ(6) * t366 - qJD(6) * t458 + t453;
t1 = [(t429 * t405 - t419 * t411) * MDP(14) + (-t397 * t359 + t365 * t460 - t366 * t459 + t444 * t458) * MDP(16) + (t365 * t459 + t397 * t444) * MDP(15) + (-t331 * t525 - t384 * t456) * MDP(22) + (-t291 * t458 + t294 * t498 - t296 * t525 - t297 * t384 + t300 * t366 + t316 * t359 + t319 * t456 + t327 * t331) * MDP(31) + (t305 * t525 - t311 * t366 + t315 * t384 + t331 * t557 - t336 * t456 - t359 * t513 - t453 * t498 + t454 * t458) * MDP(28) + (t419 * t412 + t429 * t510) * MDP(13) + (t401 * t359 + t398 * t366 + (-t412 * t460 - t458 * t510) * pkin(3)) * MDP(20) + (-t405 * t465 + t411 * t409 - t410 * t412 - t416 * t510) * MDP(9) + (t332 * t458 - t356 * t397 - t366 * t382 - t457 * t498) * MDP(25) + (-t294 * t382 + t295 * t384 - t316 * t332 - t317 * t331 + t470 * t365 + (-qJD(5) * t555 - t291 * t437 + t293 * t440) * t397) * MDP(30) + (t293 * t458 - t295 * t498 + t297 * t382 - t299 * t366 - t317 * t359 + t319 * t457 + t327 * t332 + t397 * t540) * MDP(29) + (t398 * t365 + t397 * t485 + t401 * t444 + t459 * t544) * MDP(21) + (t474 * t458 + t310 * t366 + t315 * t382 - t557 * t332 + ((-qJD(5) * t354 + t326) * t498 + t355 * t359 + t336 * qJD(5) * t397) * t440 + ((-qJD(5) * t355 - t314) * t498 - t354 * t359 + t305 * t397 + t336 * t365) * t437) * MDP(27) + ((-t382 * t440 - t530) * t365 + (t329 - t536 + (-t529 + t533) * qJD(5)) * t397) * MDP(23) + (t331 * t458 + t358 * t397 + t366 * t384 - t456 * t498) * MDP(24) + (t405 * t416 - t410 * t411) * MDP(8) + (-t359 * t458 + t366 * t498) * MDP(26) + (t291 * t316 + t293 * t317 + t294 * t300 + t295 * t299 + t296 * t327 + t297 * t319) * MDP(32) + (t365 * MDP(17) - t366 * MDP(18) - t315 * MDP(20) - t314 * MDP(21)) * t494 + 0.2e1 * t496 * t553 + (-t411 * MDP(10) - t412 * MDP(11) + MDP(13) * t445 - MDP(14) * t448) * qJD(3); (t508 + t510) * MDP(13) + (t425 + (-t409 - t489) * qJD(3)) * MDP(14) + (t502 + t479 + 0.2e1 * t554) * MDP(20) + (t501 + t450 + 0.2e1 * t572) * MDP(21) + (t461 - t534) * MDP(27) + (-t440 * t579 - t356 - t531) * MDP(28) + (-t534 + t581) * MDP(29) + (t498 * t530 + t330 + t463) * MDP(30) + (t462 + t531) * MDP(31) + (-t319 * t459 + (-t293 + t578) * t440 + (t299 * t498 + t291) * t437) * MDP(32) - qJD(1) ^ 2 * t553; (t345 * t494 + (t410 * t460 - t494 * t507) * pkin(3) + t550) * MDP(20) + (t346 * t494 + (-t410 * t459 - t486 * t494) * pkin(3) + t569) * MDP(21) + (t414 * t332 + (-t319 * t460 - t535) * t437 - t517 * t382 + (-t430 * t504 - t437 * t475 + t307) * t498 + t552) * MDP(29) + (t306 * t382 - t307 * t384 + (-t382 * t475 - t299 * t460 + (-t332 + t506) * t430) * t440 + (t384 * t475 + t300 * t460 - t331 * t430 + (t382 * t430 - t300) * qJD(5)) * t437 + t491) * MDP(30) + (-t419 * t410 - t449) * MDP(13) + (t431 * t332 + t469 * t437 + t473 * t382 + ((-qJD(5) * t430 - t349) * t440 + (-t475 + t346) * t437) * t498 + t549) * MDP(27) + (-t431 * t331 + t469 * t440 + t473 * t384 + (-t440 * t475 + t488 + t514) * t498 + t548) * MDP(28) + (t414 * t331 + (-t306 - t488) * t498 + t517 * t384 + (t475 * t498 + t535 - t577) * t440 + t551) * MDP(31) + (t425 + (t409 - t489) * qJD(3)) * MDP(10) + (t419 * t409 + t465 * t496) * MDP(14) + (t296 * t414 - t299 * t307 - t300 * t306 - t517 * t319 + t446 * t430 + t475 * t555) * MDP(32) + (t508 - t510) * MDP(11) + (t461 + t534) * MDP(25) + (-t409 ^ 2 + t410 ^ 2) * MDP(9) + t410 * t409 * MDP(8) + t580; (t343 * t494 + t550) * MDP(20) + (t342 * t494 + t569) * MDP(21) + (t534 + t581) * MDP(25) + (-pkin(4) * t332 - t336 * t528 - t343 * t382 - t468 * t498 - t451 + t549) * MDP(27) + (pkin(4) * t331 + pkin(9) * t558 - t336 * t570 - t343 * t384 + t498 * t515 + t548) * MDP(28) + (t309 * t498 - t319 * t528 + t332 * t420 - t382 * t518 - t451 + t552) * MDP(29) + (-t299 * t570 + t308 * t382 - t309 * t384 - t300 * t576 + (-t329 - t536 + (t529 + t533) * qJD(5)) * pkin(9) + t491) * MDP(30) + (t331 * t420 + (-pkin(9) * t505 - t308) * t498 + t518 * t384 + (pkin(9) * t359 - t577) * t440 + t551) * MDP(31) + (pkin(9) * t446 + t296 * t420 - t299 * t309 - t300 * t308 - t319 * t518) * MDP(32) + t580; MDP(22) * t532 + (-t382 ^ 2 + t547) * MDP(23) + t312 * MDP(24) + (t384 * t498 - t437 * t476 - t459 * t504 - t443) * MDP(25) + t359 * MDP(26) + (-t336 * t384 - t474 + t537) * MDP(27) + (t310 * t498 + t336 * t382 - t454) * MDP(28) + (-t347 * t382 - t455 + t537 + 0.2e1 * t545) * MDP(29) + (pkin(5) * t331 - qJ(6) * t332 + (t300 - t311) * t384 + (t299 - t497) * t382) * MDP(30) + (0.2e1 * t541 - t319 * t382 + t347 * t384 + (0.2e1 * qJD(6) - t310) * t498 + t454) * MDP(31) + (-pkin(5) * t293 + qJ(6) * t291 - t299 * t311 + t300 * t497 - t319 * t347) * MDP(32); (-t359 + t532) * MDP(29) + t312 * MDP(30) + (-t547 - t579) * MDP(31) + (t455 - t545 - t578) * MDP(32);];
tauc  = t1;
