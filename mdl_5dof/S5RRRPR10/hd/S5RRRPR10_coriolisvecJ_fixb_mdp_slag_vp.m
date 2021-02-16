% Calculate Coriolis joint torque vector for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:34
% EndTime: 2021-01-15 23:41:56
% DurationCPUTime: 8.43s
% Computational Cost: add. (5943->463), mult. (16181->649), div. (0->0), fcn. (12543->10), ass. (0->193)
t463 = sin(pkin(5));
t471 = cos(qJ(2));
t531 = qJD(1) * t471;
t515 = t463 * t531;
t571 = qJD(3) - t515;
t469 = cos(qJ(5));
t466 = sin(qJ(5));
t465 = cos(pkin(5));
t532 = qJD(1) * t465;
t452 = qJD(2) + t532;
t470 = cos(qJ(3));
t467 = sin(qJ(3));
t468 = sin(qJ(2));
t533 = qJD(1) * t463;
t516 = t468 * t533;
t496 = t467 * t516;
t401 = -t470 * t452 + t496;
t403 = t452 * t467 + t470 * t516;
t462 = sin(pkin(10));
t464 = cos(pkin(10));
t484 = -t401 * t462 + t464 * t403;
t552 = t484 * t466;
t344 = -t469 * t571 + t552;
t363 = t464 * t401 + t403 * t462;
t564 = qJD(5) + t363;
t570 = t344 * t564;
t346 = t466 * t571 + t469 * t484;
t569 = t346 * t564;
t520 = pkin(1) * t532;
t418 = -pkin(7) * t516 + t471 * t520;
t482 = (pkin(2) * t468 - pkin(8) * t471) * t463;
t419 = qJD(1) * t482;
t501 = -t418 * t467 + t470 * t419;
t557 = qJ(4) + pkin(8);
t507 = qJD(3) * t557;
t542 = t470 * t471;
t568 = (pkin(3) * t468 - qJ(4) * t542) * t533 + t501 + qJD(4) * t467 + t470 * t507;
t495 = t467 * t515;
t536 = t470 * t418 + t467 * t419;
t567 = -qJ(4) * t495 - qJD(4) * t470 + t467 * t507 + t536;
t430 = t462 * t467 - t464 * t470;
t566 = t571 * t430;
t431 = t462 * t470 + t464 * t467;
t535 = t571 * t431;
t500 = t564 * t469;
t521 = qJD(1) * qJD(2);
t508 = t463 * t521;
t493 = t471 * t508;
t526 = qJD(3) * t470;
t376 = -qJD(3) * t496 + t452 * t526 + t470 * t493;
t528 = qJD(2) * t471;
t512 = t467 * t528;
t527 = qJD(3) * t467;
t377 = (t468 * t526 + t512) * t533 + t452 * t527;
t340 = t376 * t462 + t464 * t377;
t544 = t466 * t340;
t565 = -t500 * t564 - t544;
t459 = t463 ^ 2;
t563 = -0.2e1 * t459 * t521;
t561 = (t468 ^ 2 - t471 ^ 2) * MDP(5);
t541 = -t462 * t567 + t464 * t568;
t540 = t462 * t568 + t464 * t567;
t449 = t468 * t520;
t421 = pkin(7) * t515 + t449;
t492 = -t421 + (-t495 + t527) * pkin(3);
t545 = t463 * t471;
t559 = pkin(1) * t468;
t415 = pkin(7) * t545 + (pkin(8) + t559) * t465;
t416 = (-pkin(2) * t471 - pkin(8) * t468 - pkin(1)) * t463;
t537 = t470 * t415 + t467 * t416;
t390 = pkin(8) * t452 + t421;
t397 = qJD(1) * t416;
t356 = t390 * t470 + t397 * t467;
t420 = qJD(2) * t482;
t410 = qJD(1) * t420;
t546 = t463 * t468;
t453 = pkin(7) * t546;
t558 = pkin(1) * t471;
t422 = (t465 * t558 - t453) * qJD(2);
t411 = qJD(1) * t422;
t473 = -qJD(3) * t356 + t470 * t410 - t467 * t411;
t494 = t468 * t508;
t306 = pkin(3) * t494 - qJ(4) * t376 - qJD(4) * t403 + t473;
t481 = -t390 * t527 + t397 * t526 + t467 * t410 + t470 * t411;
t310 = -qJ(4) * t377 - qJD(4) * t401 + t481;
t289 = t464 * t306 - t310 * t462;
t287 = -pkin(4) * t494 - t289;
t456 = pkin(3) * t462 + pkin(9);
t560 = t564 * (pkin(3) * t403 + pkin(4) * t484 + pkin(9) * t363 + qJD(5) * t456) + t287;
t341 = t376 * t464 - t377 * t462;
t524 = qJD(5) * t469;
t517 = t469 * t341 + t466 * t494 + t524 * t571;
t525 = qJD(5) * t466;
t311 = -t484 * t525 + t517;
t556 = t311 * t466;
t339 = -qJ(4) * t401 + t356;
t555 = t339 * t462;
t554 = t344 * t484;
t553 = t346 * t484;
t551 = t401 * t571;
t550 = t403 * t571;
t549 = t431 * t469;
t548 = t571 * t470;
t472 = qJD(1) ^ 2;
t547 = t459 * t472;
t334 = t464 * t339;
t543 = t467 * t571;
t336 = t469 * t340;
t290 = t462 * t306 + t464 * t310;
t427 = -t465 * t470 + t467 * t546;
t513 = t463 * t528;
t381 = -qJD(3) * t427 + t470 * t513;
t428 = t465 * t467 + t470 * t546;
t474 = -qJD(3) * t537 + t470 * t420 - t422 * t467;
t530 = qJD(2) * t468;
t514 = t463 * t530;
t317 = pkin(3) * t514 - qJ(4) * t381 - qJD(4) * t428 + t474;
t380 = qJD(3) * t428 + t463 * t512;
t480 = -t415 * t527 + t416 * t526 + t467 * t420 + t470 * t422;
t321 = -qJ(4) * t380 - qJD(4) * t427 + t480;
t296 = t462 * t317 + t464 * t321;
t355 = -t390 * t467 + t470 * t397;
t338 = -qJ(4) * t403 + t355;
t332 = pkin(3) * t571 + t338;
t305 = t462 * t332 + t334;
t502 = -t415 * t467 + t470 * t416;
t343 = -pkin(3) * t545 - qJ(4) * t428 + t502;
t351 = -qJ(4) * t427 + t537;
t320 = t462 * t343 + t464 * t351;
t539 = pkin(4) * t516 + t541;
t412 = pkin(7) * t493 + qJD(2) * t449;
t423 = t465 * pkin(1) * t530 + pkin(7) * t513;
t529 = qJD(2) * t470;
t523 = qJD(2) - t452;
t518 = t466 * t545;
t458 = -pkin(3) * t470 - pkin(2);
t510 = t557 * t467;
t288 = pkin(9) * t494 + t290;
t353 = pkin(3) * t377 + t412;
t298 = pkin(4) * t340 - pkin(9) * t341 + t353;
t506 = -t288 * t466 + t469 * t298;
t505 = t341 * t466 - t469 * t494;
t504 = t466 * t566 - t469 * t516;
t503 = t466 * t516 + t469 * t566;
t497 = t459 * t468 * t471 * MDP(4);
t366 = pkin(3) * t380 + t423;
t491 = pkin(1) * t563;
t378 = pkin(4) * t430 - pkin(9) * t431 + t458;
t490 = pkin(9) * t516 - qJD(5) * t378 + t540;
t446 = t557 * t470;
t388 = t464 * t446 - t462 * t510;
t489 = -t535 * pkin(4) - pkin(9) * t566 + qJD(5) * t388 - t492;
t488 = t288 * t469 + t298 * t466;
t300 = pkin(9) * t571 + t305;
t389 = -pkin(2) * t452 - t418;
t361 = pkin(3) * t401 + qJD(4) + t389;
t313 = pkin(4) * t363 - pkin(9) * t484 + t361;
t292 = t300 * t469 + t313 * t466;
t487 = t300 * t466 - t313 * t469;
t316 = -pkin(9) * t545 + t320;
t372 = t464 * t427 + t428 * t462;
t373 = -t427 * t462 + t428 * t464;
t414 = t453 + (-pkin(2) - t558) * t465;
t375 = pkin(3) * t427 + t414;
t329 = pkin(4) * t372 - pkin(9) * t373 + t375;
t486 = t316 * t469 + t329 * t466;
t485 = -t316 * t466 + t329 * t469;
t295 = t317 * t464 - t321 * t462;
t304 = t332 * t464 - t555;
t319 = t343 * t464 - t351 * t462;
t483 = t336 + (-t363 * t466 - t525) * t564;
t358 = t373 * t466 + t469 * t545;
t479 = t431 * t524 - t504;
t478 = -t431 * t525 - t503;
t299 = -pkin(4) * t571 - t304;
t309 = t338 * t464 - t555;
t475 = -t456 * t340 + (t299 + t309) * t564;
t457 = -pkin(3) * t464 - pkin(4);
t387 = t446 * t462 + t464 * t510;
t359 = t373 * t469 - t518;
t350 = -t380 * t462 + t381 * t464;
t349 = t464 * t380 + t381 * t462;
t323 = -qJD(5) * t518 + t350 * t466 + t373 * t524 - t469 * t514;
t322 = -qJD(5) * t358 + t350 * t469 + t466 * t514;
t315 = pkin(4) * t545 - t319;
t312 = qJD(5) * t346 + t505;
t308 = t338 * t462 + t334;
t301 = pkin(4) * t349 - pkin(9) * t350 + t366;
t294 = pkin(9) * t514 + t296;
t293 = -pkin(4) * t514 - t295;
t286 = -qJD(5) * t292 + t506;
t285 = -qJD(5) * t487 + t488;
t1 = [t561 * t563 + (-t412 * t465 - t423 * t452 + t468 * t491) * MDP(9) + (-t411 * t465 - t422 * t452 + t471 * t491) * MDP(10) + (t376 * t428 + t381 * t403) * MDP(11) + (-t376 * t427 - t377 * t428 - t380 * t403 - t381 * t401) * MDP(12) + (t381 * t571 + (-t376 * t471 + (qJD(1) * t428 + t403) * t530) * t463) * MDP(13) + (-t380 * t571 + (t377 * t471 + (-qJD(1) * t427 - t401) * t530) * t463) * MDP(14) + (-t459 * t531 + t463 * t571) * MDP(15) * t530 + (t474 * t571 + t423 * t401 + t414 * t377 + t412 * t427 + t389 * t380 + (-t473 * t471 + (qJD(1) * t502 + t355) * t530) * t463) * MDP(16) + (-t480 * t571 + t423 * t403 + t414 * t376 + t412 * t428 + t389 * t381 + (t481 * t471 + (-qJD(1) * t537 - t356) * t530) * t463) * MDP(17) + (t295 * t571 + t340 * t375 + t349 * t361 + t353 * t372 + t363 * t366 + (-t289 * t471 + (qJD(1) * t319 + t304) * t530) * t463) * MDP(18) + (-t296 * t571 + t341 * t375 + t350 * t361 + t353 * t373 + t484 * t366 + (t290 * t471 + (-qJD(1) * t320 - t305) * t530) * t463) * MDP(19) + (-t289 * t373 - t290 * t372 - t295 * t484 - t296 * t363 - t304 * t350 - t305 * t349 - t319 * t341 - t320 * t340) * MDP(20) + (t289 * t319 + t290 * t320 + t295 * t304 + t296 * t305 + t353 * t375 + t361 * t366) * MDP(21) + (t311 * t359 + t322 * t346) * MDP(22) + (-t311 * t358 - t312 * t359 - t322 * t344 - t323 * t346) * MDP(23) + (t311 * t372 + t322 * t564 + t340 * t359 + t346 * t349) * MDP(24) + (-t312 * t372 - t323 * t564 - t340 * t358 - t344 * t349) * MDP(25) + (t340 * t372 + t349 * t564) * MDP(26) + ((-qJD(5) * t486 - t294 * t466 + t301 * t469) * t564 + t485 * t340 + t286 * t372 - t487 * t349 + t293 * t344 + t315 * t312 + t287 * t358 + t299 * t323) * MDP(27) + (-(qJD(5) * t485 + t294 * t469 + t301 * t466) * t564 - t486 * t340 - t285 * t372 - t292 * t349 + t293 * t346 + t315 * t311 + t287 * t359 + t299 * t322) * MDP(28) + 0.2e1 * t497 * t521 + (MDP(6) * t513 - MDP(7) * t514) * (t452 + t532); t547 * t561 + t523 * MDP(6) * t515 + (t421 * t452 + t547 * t559 - t412) * MDP(9) + (pkin(7) * t494 + t418 * t452 + (-t465 * t521 + t547) * t558) * MDP(10) + (t376 * t467 + t403 * t548) * MDP(11) + ((t376 - t551) * t470 + (-t377 - t550) * t467) * MDP(12) + (t571 * t526 + (-t571 * t542 + (qJD(2) * t467 - t403) * t468) * t533) * MDP(13) + (-t571 * t527 + (t471 * t543 + (t401 + t529) * t468) * t533) * MDP(14) + (-pkin(2) * t377 - t412 * t470 - t501 * t571 - t421 * t401 + (-pkin(8) * t548 + t389 * t467) * qJD(3) + (-t355 * t468 + (-pkin(8) * t530 - t389 * t471) * t467) * t533) * MDP(16) + (-pkin(2) * t376 + t412 * t467 + t536 * t571 - t421 * t403 + (pkin(8) * t543 + t389 * t470) * qJD(3) + (-t389 * t542 + (-pkin(8) * t529 + t356) * t468) * t533) * MDP(17) + (t340 * t458 + t353 * t430 + t535 * t361 + t492 * t363 - t541 * t571) * MDP(18) + (t341 * t458 + t353 * t431 - t361 * t566 + t484 * t492 + t540 * t571) * MDP(19) + (-t289 * t431 - t290 * t430 + t304 * t566 - t305 * t535 - t340 * t388 + t341 * t387 + t363 * t540 + t484 * t541) * MDP(20) + (-t289 * t387 + t290 * t388 - t304 * t541 - t305 * t540 + t353 * t458 + t361 * t492) * MDP(21) + (t311 * t549 + t346 * t478) * MDP(22) + (t504 * t346 + t503 * t344 + (-t556 - t312 * t469 + (t344 * t466 - t346 * t469) * qJD(5)) * t431) * MDP(23) + (t311 * t430 + t336 * t431 + t346 * t535 + t478 * t564) * MDP(24) + (-t312 * t430 - t344 * t535 - t431 * t544 - t479 * t564) * MDP(25) + (t340 * t430 + t535 * t564) * MDP(26) + ((t378 * t469 - t388 * t466) * t340 + t286 * t430 + t387 * t312 + t287 * t466 * t431 + (t466 * t490 - t469 * t489) * t564 + t539 * t344 - t535 * t487 + t479 * t299) * MDP(27) + (-(t378 * t466 + t388 * t469) * t340 - t285 * t430 + t387 * t311 + t287 * t549 + (t466 * t489 + t469 * t490) * t564 + t539 * t346 - t535 * t292 + t478 * t299) * MDP(28) - t472 * t497 + (-t523 * MDP(7) + (-qJD(2) * t387 - t304) * MDP(18) + (-qJD(2) * t388 + t305) * MDP(19) - t571 * MDP(15)) * t516; t403 * t401 * MDP(11) + (-t401 ^ 2 + t403 ^ 2) * MDP(12) + (t376 + t551) * MDP(13) + (-t377 + t550) * MDP(14) + MDP(15) * t494 + (t356 * t571 - t389 * t403 + t473) * MDP(16) + (t355 * t571 + t389 * t401 - t481) * MDP(17) + (t308 * t571 - t361 * t484 + (-t363 * t403 + t464 * t494) * pkin(3) + t289) * MDP(18) + (t309 * t571 + t361 * t363 + (-t403 * t484 - t462 * t494) * pkin(3) - t290) * MDP(19) + ((-t340 * t462 - t341 * t464) * pkin(3) + (t305 - t308) * t484 + (-t304 + t309) * t363) * MDP(20) + (t304 * t308 - t305 * t309 + (t289 * t464 + t290 * t462 - t361 * t403) * pkin(3)) * MDP(21) + (t346 * t500 + t556) * MDP(22) + ((t311 - t570) * t469 + (-t312 - t569) * t466) * MDP(23) + (-t553 - t565) * MDP(24) + (t483 + t554) * MDP(25) - t564 * t484 * MDP(26) + (-t308 * t344 + t457 * t312 + t475 * t466 - t469 * t560 + t484 * t487) * MDP(27) + (t292 * t484 - t308 * t346 + t457 * t311 + t466 * t560 + t475 * t469) * MDP(28); (t484 * t571 + t340) * MDP(18) + (-t363 * t571 + t341) * MDP(19) + (-t363 ^ 2 - t484 ^ 2) * MDP(20) + (t304 * t484 + t305 * t363 + t353) * MDP(21) + (t483 - t554) * MDP(27) + (-t553 + t565) * MDP(28); t346 * t344 * MDP(22) + (-t344 ^ 2 + t346 ^ 2) * MDP(23) + (t517 + t570) * MDP(24) + (-t505 + t569) * MDP(25) + t340 * MDP(26) + (t292 * t564 - t299 * t346 + t506) * MDP(27) + (t299 * t344 - t487 * t564 - t488) * MDP(28) + (-MDP(24) * t552 - MDP(25) * t346 - MDP(27) * t292 + MDP(28) * t487) * qJD(5);];
tauc = t1;
