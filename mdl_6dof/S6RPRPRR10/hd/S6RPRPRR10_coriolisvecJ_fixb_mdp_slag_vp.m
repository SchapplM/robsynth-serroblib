% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:30
% EndTime: 2019-03-09 04:09:43
% DurationCPUTime: 7.43s
% Computational Cost: add. (4121->452), mult. (9348->631), div. (0->0), fcn. (6799->8), ass. (0->189)
t467 = sin(qJ(3));
t533 = qJD(1) * t467;
t457 = qJD(5) + t533;
t451 = qJD(6) + t457;
t463 = sin(pkin(10));
t470 = cos(qJ(3));
t532 = qJD(1) * t470;
t514 = t463 * t532;
t464 = cos(pkin(10));
t521 = t464 * qJD(3);
t427 = t514 - t521;
t513 = t464 * t532;
t531 = qJD(3) * t463;
t429 = t513 + t531;
t466 = sin(qJ(5));
t469 = cos(qJ(5));
t370 = t469 * t427 + t429 * t466;
t468 = cos(qJ(6));
t369 = t427 * t466 - t429 * t469;
t465 = sin(qJ(6));
t556 = t369 * t465;
t577 = -t468 * t370 + t556;
t575 = t451 * t577;
t490 = pkin(3) * t470 + qJ(4) * t467;
t437 = t490 * qJD(1);
t471 = -pkin(1) - pkin(7);
t564 = qJD(1) * t471;
t452 = qJD(2) + t564;
t552 = t464 * t470;
t383 = t463 * t437 + t452 * t552;
t515 = t463 * t533;
t360 = pkin(8) * t515 + t383;
t578 = -qJD(4) * t464 + t360;
t565 = -t463 * t466 + t469 * t464;
t576 = qJD(5) * t565;
t574 = t369 * t457;
t573 = t370 * t457;
t484 = t369 * t468 + t370 * t465;
t572 = t451 * t484;
t566 = t565 * t467;
t541 = qJD(1) * t566 + t576;
t435 = t463 * t469 + t464 * t466;
t415 = t435 * qJD(1);
t563 = t435 * qJD(5);
t540 = t467 * t415 + t563;
t442 = pkin(3) * t467 - qJ(4) * t470 + qJ(2);
t420 = t442 * qJD(1);
t441 = t467 * t452;
t421 = qJD(3) * qJ(4) + t441;
t361 = t464 * t420 - t421 * t463;
t333 = pkin(4) * t533 - pkin(8) * t429 + t361;
t362 = t463 * t420 + t464 * t421;
t335 = -pkin(8) * t427 + t362;
t305 = t333 * t466 + t335 * t469;
t301 = -pkin(9) * t370 + t305;
t523 = qJD(6) * t465;
t299 = t301 * t523;
t560 = qJD(3) * pkin(3);
t506 = -qJD(4) + t560;
t555 = t452 * t470;
t412 = -t506 - t555;
t380 = pkin(4) * t427 + t412;
t320 = pkin(5) * t370 + t380;
t570 = -t320 * t577 + t299;
t520 = qJD(1) * qJD(3);
t510 = t467 * t520;
t495 = t469 * t510;
t496 = t466 * t510;
t525 = qJD(5) * t469;
t526 = qJD(5) * t466;
t330 = -t427 * t525 - t429 * t526 + t463 * t496 - t464 * t495;
t407 = qJD(3) * t490 - qJD(4) * t470 + qJD(2);
t386 = t407 * qJD(1);
t409 = (qJD(4) + t555) * qJD(3);
t340 = t464 * t386 - t409 * t463;
t518 = pkin(8) * t464 * t467;
t476 = (pkin(4) * t470 + t518) * qJD(1);
t321 = qJD(3) * t476 + t340;
t341 = t463 * t386 + t464 * t409;
t497 = t463 * t510;
t334 = pkin(8) * t497 + t341;
t504 = t469 * t321 - t334 * t466;
t474 = -t305 * qJD(5) + t504;
t509 = t470 * t520;
t292 = pkin(5) * t509 - pkin(9) * t330 + t474;
t331 = -qJD(5) * t369 - t463 * t495 - t464 * t496;
t480 = t466 * t321 + t333 * t525 + t469 * t334 - t335 * t526;
t293 = -pkin(9) * t331 + t480;
t505 = t468 * t292 - t465 * t293;
t569 = t320 * t484 + t505;
t568 = MDP(29) * t509 + (t484 ^ 2 - t577 ^ 2) * MDP(26) + t577 * MDP(25) * t484;
t519 = 0.2e1 * qJD(1);
t567 = MDP(8) * (t467 ^ 2 - t470 ^ 2);
t425 = t464 * t442;
t508 = -t463 * t471 + pkin(4);
t368 = -pkin(8) * t552 + t467 * t508 + t425;
t551 = t467 * t471;
t389 = t463 * t442 + t464 * t551;
t553 = t463 * t470;
t381 = -pkin(8) * t553 + t389;
t542 = t466 * t368 + t469 * t381;
t561 = pkin(8) + qJ(4);
t446 = t561 * t463;
t447 = t561 * t464;
t538 = -t466 * t446 + t469 * t447;
t503 = t465 * t330 + t468 * t331;
t297 = -qJD(6) * t484 + t503;
t382 = t464 * t437 - t452 * t553;
t346 = t382 + t476;
t482 = qJD(4) * t463 + qJD(5) * t447;
t562 = t446 * t525 + t578 * t469 + (t346 + t482) * t466;
t473 = qJD(1) ^ 2;
t559 = qJ(2) * t473;
t304 = t469 * t333 - t335 * t466;
t300 = pkin(9) * t369 + t304;
t298 = pkin(5) * t457 + t300;
t558 = t298 * t468;
t557 = t301 * t468;
t472 = qJD(3) ^ 2;
t548 = t471 * t472;
t375 = t435 * t465 - t468 * t565;
t547 = -qJD(6) * t375 - t465 * t540 + t468 * t541;
t376 = t435 * t468 + t465 * t565;
t546 = qJD(6) * t376 + t465 * t541 + t468 * t540;
t406 = t565 * t470;
t544 = -qJD(3) * t406 + t467 * t563 + t415;
t529 = qJD(3) * t470;
t543 = -qJD(1) * t565 - qJD(5) * t566 - t435 * t529;
t528 = qJD(3) * t471;
t511 = t470 * t528;
t379 = t463 * t407 + t464 * t511;
t530 = qJD(3) * t467;
t522 = qJD(6) * t468;
t516 = t468 * t330 - t465 * t331 - t370 * t522;
t459 = -pkin(4) * t464 - pkin(3);
t512 = t463 * t530;
t456 = t467 * t528;
t397 = -pkin(4) * t515 + t441;
t507 = pkin(5) * t540 - t397;
t392 = t464 * t407;
t344 = t392 + (t470 * t508 + t518) * qJD(3);
t358 = pkin(8) * t512 + t379;
t502 = t469 * t344 - t358 * t466;
t501 = t469 * t368 - t381 * t466;
t500 = -t469 * t446 - t447 * t466;
t426 = pkin(4) * t553 - t470 * t471;
t499 = -t429 + t531;
t498 = qJD(6) * t298 + t293;
t449 = t467 * t509;
t343 = t469 * t346;
t353 = pkin(9) * t565 + t538;
t494 = pkin(5) * t532 + pkin(9) * t541 + t435 * qJD(4) + t538 * qJD(5) + qJD(6) * t353 - t360 * t466 + t343;
t352 = -pkin(9) * t435 + t500;
t493 = pkin(9) * t540 - qJD(6) * t352 + t562;
t403 = t435 * t467;
t492 = qJD(6) * t403 + t544;
t491 = qJD(6) * t566 - t543;
t290 = t298 * t465 + t557;
t308 = pkin(5) * t467 - pkin(9) * t406 + t501;
t404 = t435 * t470;
t309 = -pkin(9) * t404 + t542;
t489 = t308 * t465 + t309 * t468;
t488 = -t340 * t463 + t341 * t464;
t487 = -t361 * t464 - t362 * t463;
t486 = -t361 * t463 + t362 * t464;
t348 = t468 * t404 + t406 * t465;
t349 = -t404 * t465 + t406 * t468;
t483 = (-t427 - t521) * t467;
t413 = -pkin(4) * t512 + t456;
t436 = t452 * t530;
t390 = -pkin(4) * t497 + t436;
t479 = t466 * t344 + t469 * t358 + t368 * t525 - t381 * t526;
t296 = t369 * t523 + t516;
t477 = -t412 + (t452 + t564) * t470;
t475 = -qJ(4) * t529 + (t412 + t506) * t467;
t396 = -pkin(5) * t565 + t459;
t388 = -t463 * t551 + t425;
t378 = -t463 * t511 + t392;
t377 = pkin(5) * t404 + t426;
t357 = -t466 * t467 * t521 - t469 * t512 + t470 * t576;
t355 = -qJD(3) * t566 - t470 * t563;
t332 = pkin(5) * t357 + t413;
t310 = pkin(5) * t331 + t390;
t303 = qJD(6) * t349 + t355 * t465 + t468 * t357;
t302 = -qJD(6) * t348 + t355 * t468 - t357 * t465;
t295 = -pkin(9) * t357 + t479;
t294 = pkin(5) * t529 - pkin(9) * t355 - qJD(5) * t542 + t502;
t289 = -t301 * t465 + t558;
t1 = [(t340 * t388 + t341 * t389 + t361 * t378 + t362 * t379 + (t412 - t555) * t456) * MDP(17) + (-t470 * t548 + (-qJ(2) * t530 + qJD(2) * t470) * t519) * MDP(13) + (-t467 * t548 + (qJ(2) * t529 + qJD(2) * t467) * t519) * MDP(12) + (-t479 * t457 - t480 * t467 - t413 * t369 + t426 * t330 + t390 * t406 + t380 * t355 + (-qJD(1) * t542 - t305) * t529) * MDP(24) + (-t378 * t429 - t379 * t427 + (-t340 * t464 - t341 * t463) * t470 + ((t388 * t464 + t389 * t463) * qJD(1) - t487) * t530) * MDP(16) + ((t294 * t468 - t295 * t465) * t451 + t505 * t467 - t332 * t577 + t377 * t297 + t310 * t348 + t320 * t303 + (-t290 * t467 - t451 * t489) * qJD(6) + ((t308 * t468 - t309 * t465) * qJD(1) + t289) * t529) * MDP(30) + (t330 * t467 + t355 * t457 + (qJD(1) * t406 - t369) * t529) * MDP(20) + (-t331 * t467 - t357 * t457 + (-qJD(1) * t404 - t370) * t529) * MDP(21) + (t296 * t467 + t302 * t451 + (qJD(1) * t349 - t484) * t529) * MDP(27) + (-t297 * t467 - t303 * t451 + (-qJD(1) * t348 + t577) * t529) * MDP(28) + (t502 * t457 + t504 * t467 + t413 * t370 + t426 * t331 + t390 * t404 + t380 * t357 + (-t305 * t467 - t457 * t542) * qJD(5) + (qJD(1) * t501 + t304) * t529) * MDP(23) + (t377 * t296 + t299 * t467 + t320 * t302 + t310 * t349 - t332 * t484 + (-(-qJD(6) * t309 + t294) * t451 - t292 * t467) * t465 + (-(qJD(6) * t308 + t295) * t451 - t498 * t467) * t468 + (-qJD(1) * t489 - t290) * t529) * MDP(31) + (t457 * t529 + t449) * MDP(22) + (t451 * t529 + t449) * MDP(29) - 0.2e1 * MDP(7) * t449 + (t296 * t349 - t302 * t484) * MDP(25) + (-t296 * t348 - t297 * t349 + t302 * t577 + t303 * t484) * MDP(26) + (t330 * t406 - t355 * t369) * MDP(18) + (-t330 * t404 - t331 * t406 - t355 * t370 + t357 * t369) * MDP(19) + ((-qJD(1) * t379 - t341) * t467 + ((-qJD(1) * t389 - t362) * t470 + (t429 * t471 + t464 * t477) * t467) * qJD(3)) * MDP(15) + ((qJD(1) * t378 + t340) * t467 + ((qJD(1) * t388 + t361) * t470 + (t427 * t471 + t463 * t477) * t467) * qJD(3)) * MDP(14) + 0.2e1 * t520 * t567 + (MDP(6) * qJ(2) + MDP(5)) * qJD(2) * t519 + (-MDP(10) * t470 - MDP(9) * t467) * t472; -t473 * MDP(5) - MDP(6) * t559 + (-t464 * t473 + (t427 - t514) * qJD(3)) * t467 * MDP(14) + (t463 * t473 + (t429 - t513) * qJD(3)) * t467 * MDP(15) + ((-t427 * t464 + t429 * t463) * t529 + (t427 * t463 + t429 * t464) * qJD(1)) * MDP(16) + (t488 * t467 + t487 * qJD(1) + (t412 * t467 + (t486 - t441) * t470) * qJD(3)) * MDP(17) + (-t331 * t470 + t543 * t457 + (t370 * t467 - t403 * t532) * qJD(3)) * MDP(23) + (-t330 * t470 + t544 * t457 + (-t369 * t467 - t532 * t566) * qJD(3)) * MDP(24) + (-t470 * t297 + (t465 * t492 - t468 * t491) * t451 + ((-t403 * t468 - t465 * t566) * t532 - t467 * t577) * qJD(3)) * MDP(30) + (-t470 * t296 + (t465 * t491 + t468 * t492) * t451 + (-(-t403 * t465 + t468 * t566) * t532 - t467 * t484) * qJD(3)) * MDP(31) + (t467 * MDP(12) + t470 * MDP(13)) * (-t472 - t473); t467 * MDP(13) * t559 + (t452 * t483 + (-t361 * t470 - t382 * t467 + t463 * t475) * qJD(1)) * MDP(14) + (t499 * t441 + (t362 * t470 + t383 * t467 + t464 * t475) * qJD(1)) * MDP(15) + (t382 * t429 + t383 * t427 + (-qJD(4) * t427 - t361 * t533 + t341) * t464 + (qJD(4) * t429 - t362 * t533 - t340) * t463) * MDP(16) + (-t361 * t382 - t362 * t383 + (-t412 - t560) * t441 + t486 * qJD(4) + t488 * qJ(4)) * MDP(17) + (t330 * t435 - t369 * t541) * MDP(18) + (t330 * t565 - t331 * t435 + t369 * t540 - t370 * t541) * MDP(19) + (t459 * t331 - t397 * t370 + t540 * t380 - t390 * t565) * MDP(23) + (t459 * t330 + t369 * t397 + t541 * t380 + t390 * t435) * MDP(24) + (t296 * t376 - t484 * t547) * MDP(25) + (-t296 * t375 - t297 * t376 + t484 * t546 + t547 * t577) * MDP(26) + (t396 * t297 + t310 * t375 + t546 * t320 - t507 * t577) * MDP(30) + (t396 * t296 + t310 * t376 + t547 * t320 - t484 * t507) * MDP(31) + (t541 * MDP(20) - t540 * MDP(21) + (-t343 - t482 * t469 + (qJD(5) * t446 + t578) * t466) * MDP(23) + t562 * MDP(24)) * t457 + (t547 * MDP(27) - t546 * MDP(28) + (t465 * t493 - t468 * t494) * MDP(30) + (t465 * t494 + t468 * t493) * MDP(31)) * t451 + ((qJD(3) * t435 + t369) * MDP(20) + (qJD(3) * t565 + t370) * MDP(21) - t457 * MDP(22) + (qJD(3) * t500 - t304) * MDP(23) + (-qJD(3) * t538 + t305) * MDP(24) + (qJD(3) * t376 + t484) * MDP(27) + (-qJD(3) * t375 - t577) * MDP(28) - t451 * MDP(29) + ((t352 * t468 - t353 * t465) * qJD(3) - t289) * MDP(30) + (-(t352 * t465 + t353 * t468) * qJD(3) + t290) * MDP(31)) * t532 + (-t567 + (-qJ(2) * MDP(12) + t467 * MDP(7)) * t470) * t473; -t499 * MDP(14) * t533 + qJD(1) * MDP(15) * t483 + (-t427 ^ 2 - t429 ^ 2) * MDP(16) + (t361 * t429 + t362 * t427 + t436) * MDP(17) + (t331 - t574) * MDP(23) + (t330 - t573) * MDP(24) + (t297 - t572) * MDP(30) + (t296 + t575) * MDP(31); -t369 * t370 * MDP(18) + (t369 ^ 2 - t370 ^ 2) * MDP(19) + (t330 + t573) * MDP(20) + (-t331 - t574) * MDP(21) + MDP(22) * t509 + (t305 * t457 + t369 * t380 + t474) * MDP(23) + (t304 * t457 + t370 * t380 - t480) * MDP(24) + (t296 - t575) * MDP(27) + (-t297 - t572) * MDP(28) + (-(-t300 * t465 - t557) * t451 - t290 * qJD(6) + (-t369 * t577 - t451 * t523 + t468 * t509) * pkin(5) + t569) * MDP(30) + ((-t301 * t451 - t292) * t465 + (t300 * t451 - t498) * t468 + (-t369 * t484 - t451 * t522 - t465 * t509) * pkin(5) + t570) * MDP(31) + t568; (t516 - t575) * MDP(27) + (-t503 - t572) * MDP(28) + (t290 * t451 + t569) * MDP(30) + (t289 * t451 - t465 * t292 - t468 * t293 + t570) * MDP(31) + (MDP(27) * t556 + MDP(28) * t484 - MDP(30) * t290 - MDP(31) * t558) * qJD(6) + t568;];
tauc  = t1;
