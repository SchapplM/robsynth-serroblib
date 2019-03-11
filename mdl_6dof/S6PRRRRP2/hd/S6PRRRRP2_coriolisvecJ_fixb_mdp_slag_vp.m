% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:29
% EndTime: 2019-03-09 00:04:39
% DurationCPUTime: 5.23s
% Computational Cost: add. (5941->434), mult. (14478->576), div. (0->0), fcn. (10826->10), ass. (0->199)
t446 = sin(qJ(3));
t564 = -pkin(9) - pkin(8);
t510 = qJD(3) * t564;
t420 = t446 * t510;
t449 = cos(qJ(3));
t421 = t449 * t510;
t445 = sin(qJ(4));
t563 = cos(qJ(4));
t425 = t564 * t446;
t426 = t564 * t449;
t567 = t563 * t425 + t445 * t426;
t354 = qJD(4) * t567 + t563 * t420 + t445 * t421;
t475 = -t445 * t446 + t449 * t563;
t450 = cos(qJ(2));
t442 = sin(pkin(6));
t525 = qJD(1) * t442;
t505 = t450 * t525;
t397 = t475 * t505;
t578 = -t354 + t397;
t448 = cos(qJ(5));
t518 = qJD(5) * t448;
t501 = qJD(2) * t563;
t522 = qJD(2) * t446;
t413 = t445 * t522 - t449 * t501;
t544 = t413 * t448;
t577 = t518 + t544;
t444 = sin(qJ(5));
t519 = qJD(5) * t444;
t545 = t413 * t444;
t576 = -t519 - t545;
t538 = t445 * t449;
t418 = t446 * t563 + t538;
t515 = qJD(3) + qJD(4);
t390 = t515 * t418;
t575 = t390 * qJD(2);
t409 = qJD(5) + t413;
t551 = t575 * t448;
t574 = pkin(10) * (t409 * t519 - t551);
t573 = MDP(5) * t449;
t572 = MDP(6) * (t446 ^ 2 - t449 ^ 2);
t447 = sin(qJ(2));
t506 = t447 * t525;
t497 = -t564 * qJD(2) + t506;
t443 = cos(pkin(6));
t524 = qJD(1) * t443;
t394 = -t497 * t446 + t449 * t524;
t523 = qJD(2) * t442;
t499 = qJD(1) * t523;
t490 = t450 * t499;
t367 = qJD(3) * t394 + t449 * t490;
t395 = t446 * t524 + t449 * t497;
t368 = -qJD(3) * t395 - t446 * t490;
t558 = qJD(3) * pkin(3);
t388 = t394 + t558;
t500 = t563 * qJD(4);
t520 = qJD(4) * t445;
t308 = t367 * t563 + t445 * t368 + t388 * t500 - t395 * t520;
t516 = qJD(2) * qJD(3);
t498 = t446 * t516;
t410 = pkin(3) * t498 + t447 * t499;
t389 = t515 * t475;
t454 = t389 * qJD(2);
t334 = pkin(4) * t575 - pkin(10) * t454 + t410;
t387 = t563 * t395;
t349 = t445 * t388 + t387;
t342 = pkin(10) * t515 + t349;
t438 = -pkin(3) * t449 - pkin(2);
t407 = qJD(2) * t438 - t505;
t415 = -qJD(2) * t538 - t446 * t501;
t363 = pkin(4) * t413 + pkin(10) * t415 + t407;
t468 = t448 * t308 + t444 * t334 - t342 * t519 + t363 * t518;
t557 = qJ(6) * t575;
t294 = qJD(6) * t409 + t468 + t557;
t493 = t444 * t308 - t448 * t334 + t342 * t518 + t363 * t519;
t562 = pkin(5) * t575;
t296 = t493 - t562;
t571 = t294 * t448 + t296 * t444;
t350 = t445 * t394 + t387;
t489 = pkin(3) * t520 - t350;
t403 = t445 * t425 - t426 * t563;
t528 = qJD(4) * t403 - t418 * t505 + t445 * t420 - t421 * t563;
t513 = t446 * t558;
t344 = pkin(4) * t390 - pkin(10) * t389 + t513;
t383 = -pkin(4) * t475 - pkin(10) * t418 + t438;
t570 = -t383 * t518 + t403 * t519 + t578 * t448 + (-t344 + t506) * t444;
t527 = t444 * t383 + t448 * t403;
t569 = pkin(5) * t576 + qJ(6) * t577 + qJD(6) * t444;
t568 = t506 - t513;
t466 = t448 * t415 - t444 * t515;
t346 = -qJD(5) * t466 + t444 * t454;
t566 = t466 ^ 2;
t565 = t409 ^ 2;
t561 = pkin(5) * t415;
t560 = pkin(3) * qJD(4);
t559 = qJD(2) * pkin(2);
t320 = t342 * t448 + t363 * t444;
t556 = t320 * t409;
t495 = t448 * t515;
t345 = -qJD(5) * t495 - t415 * t519 - t448 * t454;
t555 = t345 * t444;
t554 = t346 * t448;
t436 = pkin(3) * t445 + pkin(10);
t553 = t575 * t436;
t552 = t575 * t444;
t398 = -t415 * t444 - t495;
t550 = t398 * t409;
t549 = t398 * t444;
t548 = t466 * t398;
t547 = t466 * t409;
t546 = t466 * t448;
t543 = t418 * t444;
t542 = t418 * t448;
t541 = t442 * t447;
t540 = t442 * t450;
t452 = qJD(2) ^ 2;
t539 = t442 * t452;
t386 = t445 * t395;
t451 = qJD(3) ^ 2;
t537 = t446 * t451;
t536 = t449 * t451;
t535 = qJ(6) * t390 - qJD(6) * t475 - t570;
t370 = t397 * t444 - t448 * t506;
t534 = -pkin(5) * t390 + qJD(5) * t527 - t344 * t448 + t354 * t444 - t370;
t484 = pkin(5) * t444 - qJ(6) * t448;
t485 = t448 * pkin(5) + t444 * qJ(6);
t533 = t484 * t389 + (qJD(5) * t485 - qJD(6) * t448) * t418 + t528;
t532 = t349 + t569;
t531 = t569 - t489;
t348 = t388 * t563 - t386;
t382 = -pkin(4) * t415 + pkin(10) * t413;
t530 = t448 * t348 + t444 * t382;
t351 = t394 * t563 - t386;
t372 = pkin(3) * t522 + t382;
t529 = t448 * t351 + t444 * t372;
t521 = qJD(2) * t447;
t319 = -t342 * t444 + t363 * t448;
t517 = qJD(6) - t319;
t514 = t563 * pkin(3);
t511 = t447 * t539;
t509 = t444 * t563;
t508 = t448 * t563;
t503 = t442 * t521;
t502 = t450 * t523;
t496 = t409 * t448;
t494 = pkin(3) * t500;
t309 = t445 * t367 - t563 * t368 + t388 * t520 + t395 * t500;
t492 = t446 * t502;
t491 = t449 * t502;
t341 = -pkin(4) * t515 - t348;
t487 = t309 * t444 - t320 * t415 + t341 * t518;
t311 = -pkin(5) * t409 + t517;
t313 = qJ(6) * t409 + t320;
t482 = t311 * t448 - t313 * t444;
t481 = t341 * t413 - t553;
t480 = -t348 * t444 + t382 * t448;
t424 = -pkin(4) - t485;
t411 = t443 * t449 - t446 * t541;
t412 = t443 * t446 + t449 * t541;
t374 = t445 * t411 + t412 * t563;
t361 = t374 * t444 + t448 * t540;
t362 = t374 * t448 - t444 * t540;
t297 = pkin(5) * t346 + qJ(6) * t345 + qJD(6) * t466 + t309;
t323 = t398 * pkin(5) + qJ(6) * t466 + t341;
t479 = -t297 * t448 - t311 * t415 + t323 * t519;
t478 = -t297 * t444 + t313 * t415 - t323 * t544;
t477 = -t309 * t448 + t319 * t415 + t341 * t519;
t476 = t411 * t563 - t445 * t412;
t473 = t389 * t444 + t418 * t518;
t472 = -t389 * t448 + t418 * t519;
t471 = t559 * qJD(2);
t470 = -t323 * t466 + t493;
t469 = t407 * t415 - t309;
t464 = (-t409 * t518 - t552) * pkin(10);
t463 = -0.2e1 * qJD(3) * t559;
t462 = -t436 * t519 + t448 * t494;
t460 = t311 * t577 + t313 * t576 + t571;
t459 = qJD(5) * t482 + t571;
t458 = -t555 - t554 + (-t546 + t549) * qJD(5);
t457 = ((-t345 - t550) * t448 + (-t346 + t547) * t444) * MDP(20) + (-t466 * t496 - t555) * MDP(19) + (-t398 * t415 - t444 * t565 + t551) * MDP(22) + (t409 * t496 - t415 * t466 + t552) * MDP(21) + (t413 * t515 + t454) * MDP(14) + (-t415 * t515 - t575) * MDP(15) + (-t413 ^ 2 + t415 ^ 2) * MDP(13) + (-MDP(12) * t413 + MDP(23) * t409) * t415;
t456 = t407 * t413 - t308;
t437 = -t514 - pkin(4);
t416 = -t514 + t424;
t406 = t415 * qJ(6);
t392 = -qJD(3) * t412 - t492;
t391 = qJD(3) * t411 + t491;
t358 = -pkin(5) * t466 + qJ(6) * t398;
t353 = t418 * t484 - t567;
t336 = pkin(5) * t475 - t383 * t448 + t403 * t444;
t335 = -qJ(6) * t475 + t527;
t330 = qJD(4) * t374 + t445 * t391 - t392 * t563;
t329 = qJD(4) * t476 + t391 * t563 + t445 * t392;
t326 = -t345 + t550;
t325 = -t480 + t561;
t324 = -t406 + t530;
t322 = t351 * t444 - t372 * t448 + t561;
t321 = -t406 + t529;
t303 = qJD(5) * t362 + t329 * t444 - t448 * t503;
t302 = -qJD(5) * t361 + t329 * t448 + t444 * t503;
t1 = [-MDP(3) * t511 - t450 * MDP(4) * t539 + (-t449 * t511 + (t392 - t492) * qJD(3)) * MDP(10) + (t446 * t511 + (-t391 - t491) * qJD(3)) * MDP(11) + (-t330 * t515 + (t413 * t521 - t450 * t575) * t442) * MDP(17) + (-t329 * t515 + (-t389 * t450 - t447 * t415) * t523) * MDP(18) + (-t302 * t398 - t303 * t466 - t345 * t361 - t346 * t362) * MDP(27) + (t294 * t362 + t296 * t361 - t297 * t476 + t302 * t313 + t303 * t311 + t323 * t330) * MDP(29) + (MDP(24) + MDP(26)) * (-t303 * t409 + t330 * t398 - t346 * t476 - t361 * t575) + (-MDP(25) + MDP(28)) * (t302 * t409 + t330 * t466 - t345 * t476 + t362 * t575); 0.2e1 * t498 * t573 - 0.2e1 * t516 * t572 + MDP(7) * t536 - MDP(8) * t537 + (-pkin(8) * t536 + t446 * t463) * MDP(10) + (pkin(8) * t537 + t449 * t463) * MDP(11) + (-t415 * t389 + t418 * t454) * MDP(12) + (-t389 * t413 + t415 * t390 - t418 * t575 + t454 * t475) * MDP(13) + (t407 * t390 - t410 * t475 - t413 * t568 + t438 * t575) * MDP(17) + (t407 * t389 + t410 * t418 + t568 * t415 + t438 * t454) * MDP(18) + (-t345 * t542 + t466 * t472) * MDP(19) + ((-t398 * t448 + t444 * t466) * t389 + (t555 - t554 + (t546 + t549) * qJD(5)) * t418) * MDP(20) + (t345 * t475 - t390 * t466 - t409 * t472 + t542 * t575) * MDP(21) + (t346 * t475 - t390 * t398 - t409 * t473 - t543 * t575) * MDP(22) + (t390 * t409 - t475 * t575) * MDP(23) + (t493 * t475 + t319 * t390 - t567 * t346 + t370 * t409 + t528 * t398 + ((-qJD(5) * t403 + t344) * t409 + t383 * t575 + t341 * qJD(5) * t418) * t448 + ((-qJD(5) * t383 - t354) * t409 - t403 * t575 + t309 * t418 + t341 * t389) * t444) * MDP(24) + (t309 * t542 - t320 * t390 - t472 * t341 + t345 * t567 + t409 * t570 - t466 * t528 + t468 * t475 - t527 * t575) * MDP(25) + (t296 * t475 + t297 * t543 - t311 * t390 + t323 * t473 - t336 * t575 + t346 * t353 + t398 * t533 - t409 * t534) * MDP(26) + (-t335 * t346 - t336 * t345 - t534 * t466 - t535 * t398 + t482 * t389 + (-t294 * t444 + t296 * t448 + (-t311 * t444 - t313 * t448) * qJD(5)) * t418) * MDP(27) + (-t294 * t475 - t297 * t542 + t313 * t390 + t323 * t472 + t335 * t575 + t345 * t353 + t409 * t535 + t466 * t533) * MDP(28) + (t294 * t335 + t296 * t336 + t297 * t353 + t311 * t534 + t313 * t535 + t323 * t533) * MDP(29) + (t389 * MDP(14) - t390 * MDP(15) - t528 * MDP(17) + MDP(18) * t578) * t515; (t350 * t515 + (-t413 * t522 - t515 * t520) * pkin(3) + t469) * MDP(17) + (t297 * t416 - t311 * t322 - t313 * t321 - t531 * t323 + (t311 * t509 + t313 * t508) * t560 + t459 * t436) * MDP(29) + (t351 * t515 + (t415 * t522 - t500 * t515) * pkin(3) + t456) * MDP(18) + (t416 * t346 + (t323 * t413 - t553) * t444 - t531 * t398 + (-t436 * t518 - t444 * t494 + t322) * t409 + t479) * MDP(26) + t449 * t471 * MDP(11) + (t321 * t398 + t322 * t466 + (-t398 * t508 - t466 * t509) * t560 + t458 * t436 + t460) * MDP(27) + (t416 * t345 + (-qJD(5) * t323 + t553) * t448 - t531 * t466 + (-t321 + t462) * t409 + t478) * MDP(28) + t452 * t572 + (t437 * t346 + t481 * t444 + t489 * t398 + ((-qJD(5) * t436 - t372) * t448 + (-t494 + t351) * t444) * t409 + t477) * MDP(24) + t457 + (-t437 * t345 + t481 * t448 - t489 * t466 + (-t462 + t529) * t409 + t487) * MDP(25) + (MDP(10) * t471 - t452 * t573) * t446; (-pkin(4) * t346 + t341 * t545 - t349 * t398 - t409 * t480 + t464 + t477) * MDP(24) + (pkin(4) * t345 + t341 * t544 + t349 * t466 + t409 * t530 + t487 + t574) * MDP(25) + (t348 * t515 + t456) * MDP(18) + (pkin(10) * t459 + t297 * t424 - t311 * t325 - t313 * t324 - t323 * t532) * MDP(29) + (pkin(10) * t458 + t324 * t398 + t325 * t466 + t460) * MDP(27) + (-t323 * t518 - t324 * t409 + t345 * t424 - t466 * t532 + t478 - t574) * MDP(28) + (t323 * t545 + t325 * t409 + t346 * t424 - t398 * t532 + t464 + t479) * MDP(26) + (t349 * t515 + t469) * MDP(17) + t457; -MDP(19) * t548 + (-t398 ^ 2 + t566) * MDP(20) + t326 * MDP(21) + (-t346 - t547) * MDP(22) + t575 * MDP(23) + (t341 * t466 - t493 + t556) * MDP(24) + (t319 * t409 + t341 * t398 - t468) * MDP(25) + (-t358 * t398 - t470 + t556 + 0.2e1 * t562) * MDP(26) + (pkin(5) * t345 - qJ(6) * t346 - (t313 - t320) * t466 + (t311 - t517) * t398) * MDP(27) + (0.2e1 * t557 - t323 * t398 - t358 * t466 + (0.2e1 * qJD(6) - t319) * t409 + t468) * MDP(28) + (-pkin(5) * t296 + qJ(6) * t294 - t311 * t320 + t313 * t517 - t323 * t358) * MDP(29); t326 * MDP(27) + (-t565 - t566) * MDP(28) + (-t313 * t409 + t470 - t562) * MDP(29) + (-t548 - t575) * MDP(26);];
tauc  = t1;
