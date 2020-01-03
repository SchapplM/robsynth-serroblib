% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:58
% EndTime: 2019-12-31 22:26:07
% DurationCPUTime: 5.07s
% Computational Cost: add. (4426->402), mult. (10681->548), div. (0->0), fcn. (7896->8), ass. (0->185)
t552 = qJD(4) + qJD(5);
t459 = sin(qJ(3));
t463 = cos(qJ(2));
t549 = cos(qJ(3));
t501 = qJD(1) * t549;
t460 = sin(qJ(2));
t519 = qJD(1) * t460;
t567 = -t459 * t519 + t463 * t501;
t568 = -t567 + t552;
t551 = pkin(6) + pkin(7);
t439 = t551 * t463;
t429 = qJD(1) * t439;
t415 = t459 * t429;
t437 = t551 * t460;
t427 = qJD(1) * t437;
t384 = -t427 * t549 - t415;
t500 = qJD(3) * t549;
t562 = -pkin(2) * t500 + t384;
t531 = t459 * t463;
t414 = -qJD(1) * t531 - t460 * t501;
t378 = -pkin(3) * t414 - pkin(8) * t567;
t363 = pkin(2) * t519 + t378;
t458 = sin(qJ(4));
t462 = cos(qJ(4));
t566 = -t462 * t363 + t458 * t562;
t454 = qJD(2) + qJD(3);
t393 = -t414 * t458 - t462 * t454;
t461 = cos(qJ(5));
t457 = sin(qJ(5));
t478 = t414 * t462 - t454 * t458;
t538 = t478 * t457;
t341 = t461 * t393 - t538;
t409 = qJD(4) - t567;
t404 = qJD(5) + t409;
t565 = t341 * t404;
t479 = t393 * t457 + t461 * t478;
t564 = t404 * t479;
t424 = t457 * t462 + t458 * t461;
t563 = t568 * t424;
t422 = t457 * t458 - t461 * t462;
t525 = t568 * t422;
t452 = -pkin(2) * t463 - pkin(1);
t435 = t452 * qJD(1);
t361 = -pkin(3) * t567 + pkin(8) * t414 + t435;
t416 = t549 * t429;
t546 = qJD(2) * pkin(2);
t417 = -t427 + t546;
t382 = t459 * t417 + t416;
t365 = pkin(8) * t454 + t382;
t321 = t361 * t458 + t365 * t462;
t311 = -pkin(9) * t393 + t321;
t515 = qJD(5) * t457;
t309 = t311 * t515;
t381 = t417 * t549 - t415;
t364 = -t454 * pkin(3) - t381;
t338 = t393 * pkin(4) + t364;
t561 = t338 * t341 + t309;
t376 = t567 * t454;
t516 = qJD(4) * t462;
t517 = qJD(4) * t458;
t336 = t462 * t376 + t414 * t517 + t454 * t516;
t425 = t460 * t549 + t531;
t389 = t454 * t425;
t377 = t389 * qJD(1);
t512 = qJD(1) * qJD(2);
t499 = t460 * t512;
t325 = pkin(2) * t499 + pkin(3) * t377 - pkin(8) * t376;
t323 = t462 * t325;
t506 = qJD(2) * t551;
t490 = qJD(1) * t506;
t418 = t460 * t490;
t419 = t463 * t490;
t518 = qJD(3) * t459;
t328 = t417 * t500 - t418 * t549 - t459 * t419 - t429 * t518;
t468 = -qJD(4) * t321 - t328 * t458 + t323;
t289 = pkin(4) * t377 - pkin(9) * t336 + t468;
t337 = -qJD(4) * t478 + t376 * t458;
t472 = t458 * t325 + t462 * t328 + t361 * t516 - t365 * t517;
t290 = -pkin(9) * t337 + t472;
t496 = t461 * t289 - t457 * t290;
t560 = t338 * t479 + t496;
t559 = t377 * MDP(29) + (-t341 ^ 2 + t479 ^ 2) * MDP(26) - t341 * t479 * MDP(25);
t558 = -0.2e1 * t512;
t557 = MDP(5) * (t460 ^ 2 - t463 ^ 2);
t370 = t424 * t425;
t383 = -t459 * t427 + t416;
t488 = pkin(2) * t518 - t383;
t537 = t567 * t458;
t556 = (t517 - t537) * pkin(4);
t555 = t458 * t363 + t562 * t462;
t476 = -t459 * t460 + t463 * t549;
t388 = t454 * t476;
t529 = t462 * t388;
t474 = -t425 * t517 + t529;
t554 = -t549 * t437 - t459 * t439;
t553 = qJD(1) * t425;
t495 = t336 * t457 + t461 * t337;
t299 = -qJD(5) * t479 + t495;
t550 = -pkin(8) - pkin(9);
t548 = t462 * pkin(4);
t449 = pkin(2) * t459 + pkin(8);
t547 = -pkin(9) - t449;
t320 = t462 * t361 - t365 * t458;
t310 = pkin(9) * t478 + t320;
t308 = pkin(4) * t409 + t310;
t545 = t308 * t461;
t544 = t311 * t461;
t543 = t336 * t458;
t541 = t377 * t462;
t540 = t393 * t409;
t539 = t478 * t409;
t536 = t567 * t462;
t535 = t425 * t458;
t534 = t425 * t462;
t533 = t458 * t377;
t532 = t458 * t388;
t464 = qJD(2) ^ 2;
t530 = t460 * t464;
t397 = -t459 * t437 + t439 * t549;
t392 = t462 * t397;
t528 = t463 * t464;
t465 = qJD(1) ^ 2;
t527 = t463 * t465;
t524 = t458 * t378 + t462 * t381;
t380 = -pkin(3) * t476 - pkin(8) * t425 + t452;
t522 = t458 * t380 + t392;
t521 = t556 + t488;
t514 = qJD(5) * t461;
t511 = pkin(9) * t537;
t510 = t460 * t546;
t507 = t461 * t336 - t457 * t337 - t393 * t514;
t505 = qJD(4) * t550;
t356 = t364 * t516;
t498 = qJD(4) * t547;
t497 = pkin(1) * t558;
t494 = t462 * t378 - t381 * t458;
t493 = t409 * t462;
t492 = qJD(5) * t308 + t290;
t329 = t417 * t518 - t459 * t418 + t549 * t419 + t429 * t500;
t450 = -pkin(2) * t549 - pkin(3);
t489 = -t414 * pkin(4) - pkin(9) * t536;
t487 = -t382 + t556;
t486 = -t321 * t414 + t329 * t458 + t356;
t453 = t462 * pkin(9);
t421 = t449 * t462 + t453;
t485 = qJD(5) * t421 - t462 * t498 + t489 - t566;
t438 = pkin(8) * t462 + t453;
t484 = qJD(5) * t438 - t462 * t505 + t489 + t494;
t420 = t547 * t458;
t483 = -qJD(5) * t420 - t458 * t498 - t511 + t555;
t436 = t550 * t458;
t482 = -qJD(5) * t436 - t458 * t505 - t511 + t524;
t294 = t308 * t457 + t544;
t480 = -t364 * t567 - t377 * t449;
t477 = t320 * t414 - t329 * t462 + t364 * t517;
t475 = t425 * t516 + t532;
t473 = t414 * t435 - t329;
t335 = pkin(3) * t389 - pkin(8) * t388 + t510;
t428 = t460 * t506;
t430 = t463 * t506;
t345 = qJD(3) * t554 - t549 * t428 - t459 * t430;
t471 = t458 * t335 + t462 * t345 + t380 * t516 - t397 * t517;
t298 = t478 * t515 + t507;
t293 = -t311 * t457 + t545;
t305 = pkin(4) * t337 + t329;
t470 = t293 * t414 + t305 * t422 + t338 * t563;
t469 = -t294 * t414 + t305 * t424 - t338 * t525;
t467 = -t435 * t567 - t328;
t346 = qJD(3) * t397 - t459 * t428 + t430 * t549;
t466 = (-t298 * t422 - t299 * t424 + t341 * t525 + t479 * t563) * MDP(26) + (t298 * t424 + t479 * t525) * MDP(25) + ((t336 - t540) * t462 + (-t337 + t539) * t458) * MDP(19) + (t377 * t424 - t404 * t525 - t414 * t479) * MDP(27) + (-t341 * t414 - t377 * t422 - t404 * t563) * MDP(28) + (-t478 * t493 + t543) * MDP(18) + (-t409 ^ 2 * t458 - t393 * t414 + t541) * MDP(21) + (t409 * t493 - t414 * t478 + t533) * MDP(20) + t376 * MDP(13) + (t414 ^ 2 - t567 ^ 2) * MDP(12) + (MDP(11) * t567 + t409 * MDP(22) + t404 * MDP(29)) * t414 + (-t567 * MDP(13) + (-t414 - t553) * MDP(14)) * t454;
t451 = -pkin(3) - t548;
t434 = t450 - t548;
t374 = t462 * t380;
t371 = t422 * t425;
t362 = pkin(4) * t535 - t554;
t352 = t377 * t476;
t331 = t462 * t335;
t324 = -pkin(9) * t535 + t522;
t316 = -pkin(4) * t476 - pkin(9) * t534 - t397 * t458 + t374;
t314 = pkin(4) * t475 + t346;
t303 = -t515 * t535 + (t534 * t552 + t532) * t461 + t474 * t457;
t302 = -t370 * t552 - t422 * t388;
t295 = -pkin(9) * t475 + t471;
t292 = -pkin(9) * t529 + pkin(4) * t389 - t345 * t458 + t331 + (-t392 + (pkin(9) * t425 - t380) * t458) * qJD(4);
t1 = [(t376 * t476 - t377 * t425 + t388 * t567 + t389 * t414) * MDP(12) + (t377 * t452 + t389 * t435 + (-qJD(1) * t476 - t567) * t510) * MDP(16) + (t376 * t452 + t388 * t435 + (-t414 + t553) * t510) * MDP(17) + (-pkin(6) * t528 + t460 * t497) * MDP(9) - MDP(7) * t530 + (pkin(6) * t530 + t463 * t497) * MDP(10) + t557 * t558 + ((t292 * t461 - t295 * t457) * t404 + (t316 * t461 - t324 * t457) * t377 - t496 * t476 + t293 * t389 + t314 * t341 + t362 * t299 + t305 * t370 + t338 * t303 + ((-t316 * t457 - t324 * t461) * t404 + t294 * t476) * qJD(5)) * MDP(30) + (t299 * t476 - t303 * t404 - t341 * t389 - t370 * t377) * MDP(28) + (-t294 * t389 + t362 * t298 + t338 * t302 - t305 * t371 - t309 * t476 - t314 * t479 + (-(-qJD(5) * t324 + t292) * t404 - t316 * t377 + t289 * t476) * t457 + (-(qJD(5) * t316 + t295) * t404 - t324 * t377 + t492 * t476) * t461) * MDP(31) + (-t298 * t476 + t302 * t404 - t371 * t377 - t389 * t479) * MDP(27) + (-t336 * t476 + t377 * t534 - t389 * t478 + t409 * t474) * MDP(20) + (t389 * t404 - t352) * MDP(29) + (t389 * t409 - t352) * MDP(22) + (-t321 * t389 + t329 * t534 - t336 * t554 - t346 * t478 + t364 * t474 - t377 * t522 - t409 * t471 + t472 * t476) * MDP(24) + (t346 * t393 - t554 * t337 + t425 * t356 + (-t397 * t516 + t331) * t409 + t374 * t377 - (-t365 * t516 + t323) * t476 + t320 * t389 + (t329 * t425 + t364 * t388 + (-qJD(4) * t380 - t345) * t409 - t397 * t377 - (-qJD(4) * t361 - t328) * t476) * t458) * MDP(23) + (MDP(13) * t388 - MDP(14) * t389 - MDP(16) * t346 - MDP(17) * t345) * t454 + 0.2e1 * t463 * MDP(4) * t499 + MDP(6) * t528 + (-t298 * t371 - t302 * t479) * MDP(25) + (-t298 * t370 + t299 * t371 - t302 * t341 + t303 * t479) * MDP(26) + ((-t393 * t462 + t458 * t478) * t388 + (-t543 - t337 * t462 + (t393 * t458 + t462 * t478) * qJD(4)) * t425) * MDP(19) + (t336 * t534 - t474 * t478) * MDP(18) + (t337 * t476 - t389 * t393 - t409 * t475 - t425 * t533) * MDP(21) + (t376 * t425 - t388 * t414) * MDP(11); t466 + (t384 * t454 + (t414 * t519 - t454 * t500) * pkin(2) + t467) * MDP(17) + ((t420 * t461 - t421 * t457) * t377 + t434 * t299 + (t457 * t483 - t461 * t485) * t404 + t521 * t341 + t470) * MDP(30) + (t450 * t336 + t480 * t462 - t488 * t478 + (t449 * t517 + t555) * t409 + t486) * MDP(24) + (t383 * t454 + (-t454 * t518 + t519 * t567) * pkin(2) + t473) * MDP(16) + (-(t420 * t457 + t421 * t461) * t377 + t434 * t298 + (t457 * t485 + t461 * t483) * t404 - t521 * t479 + t469) * MDP(31) + (t450 * t337 + t480 * t458 + t488 * t393 + (-t449 * t516 + t566) * t409 + t477) * MDP(23) + t465 * t557 - t460 * MDP(4) * t527 + (MDP(9) * t460 * t465 + MDP(10) * t527) * pkin(1); t466 + (-pkin(3) * t337 - t382 * t393 - t364 * t537 - t494 * t409 + (-t516 * t409 - t533) * pkin(8) + t477) * MDP(23) + (-pkin(3) * t336 + t382 * t478 - t364 * t536 + t524 * t409 + (t409 * t517 - t541) * pkin(8) + t486) * MDP(24) + (-(t436 * t457 + t438 * t461) * t377 + t451 * t298 + (t457 * t484 + t461 * t482) * t404 - t487 * t479 + t469) * MDP(31) + ((t436 * t461 - t438 * t457) * t377 + t451 * t299 + (t457 * t482 - t461 * t484) * t404 + t487 * t341 + t470) * MDP(30) + (t382 * t454 + t473) * MDP(16) + (t381 * t454 + t467) * MDP(17); -t478 * t393 * MDP(18) + (-t393 ^ 2 + t478 ^ 2) * MDP(19) + (t336 + t540) * MDP(20) + (-t337 - t539) * MDP(21) + t377 * MDP(22) + (t321 * t409 + t364 * t478 + t468) * MDP(23) + (t320 * t409 + t364 * t393 - t472) * MDP(24) + (t298 + t565) * MDP(27) + (-t299 - t564) * MDP(28) + (-(-t310 * t457 - t544) * t404 - t294 * qJD(5) + (t341 * t478 + t461 * t377 - t404 * t515) * pkin(4) + t560) * MDP(30) + ((-t311 * t404 - t289) * t457 + (t310 * t404 - t492) * t461 + (-t457 * t377 - t404 * t514 - t478 * t479) * pkin(4) + t561) * MDP(31) + t559; (t507 + t565) * MDP(27) + (-t495 - t564) * MDP(28) + (t294 * t404 + t560) * MDP(30) + (-t457 * t289 - t461 * t290 + t293 * t404 + t561) * MDP(31) + (MDP(27) * t538 + MDP(28) * t479 - MDP(30) * t294 - MDP(31) * t545) * qJD(5) + t559;];
tauc = t1;
