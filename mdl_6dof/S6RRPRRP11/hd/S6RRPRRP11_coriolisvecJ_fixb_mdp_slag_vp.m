% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:45
% EndTime: 2019-03-09 12:48:55
% DurationCPUTime: 5.51s
% Computational Cost: add. (4590->462), mult. (10206->609), div. (0->0), fcn. (6425->6), ass. (0->213)
t474 = cos(qJ(2));
t470 = sin(qJ(4));
t471 = sin(qJ(2));
t562 = t470 * t471;
t493 = pkin(4) * t474 - pkin(9) * t562;
t542 = qJD(1) * t471;
t460 = pkin(2) * t542;
t501 = pkin(8) * t471 - qJ(3) * t474;
t391 = qJD(1) * t501 + t460;
t541 = qJD(1) * t474;
t456 = pkin(7) * t541;
t427 = pkin(3) * t541 + t456;
t473 = cos(qJ(4));
t505 = -t391 * t470 + t473 * t427;
t536 = qJD(4) * t470;
t475 = -pkin(2) - pkin(8);
t578 = pkin(9) - t475;
t600 = qJD(1) * t493 - t578 * t536 + t505;
t432 = t578 * t473;
t523 = t473 * t542;
t547 = t473 * t391 + t470 * t427;
t597 = pkin(9) * t523 + qJD(4) * t432 + t547;
t530 = t470 * qJD(2);
t418 = -t473 * t541 - t530;
t472 = cos(qJ(5));
t520 = t470 * t541;
t539 = qJD(2) * t473;
t419 = -t520 + t539;
t469 = sin(qJ(5));
t569 = t419 * t469;
t360 = -t472 * t418 + t569;
t358 = t360 ^ 2;
t528 = qJD(1) * qJD(2);
t518 = t474 * t528;
t498 = t418 * t469 + t472 * t419;
t581 = t498 ^ 2;
t599 = MDP(26) * t518 + (-t358 + t581) * MDP(23);
t580 = pkin(3) + pkin(7);
t598 = qJ(6) * t360;
t420 = t469 * t473 + t470 * t472;
t583 = qJD(4) + qJD(5);
t366 = t583 * t420;
t489 = t420 * t471;
t383 = qJD(1) * t489;
t549 = t366 + t383;
t533 = qJD(5) * t469;
t560 = t472 * t473;
t564 = t469 * t470;
t548 = -t469 * t536 - t470 * t533 + t472 * t523 - t542 * t564 + t560 * t583;
t455 = pkin(7) * t542;
t596 = qJD(3) + t455;
t519 = t471 * t528;
t375 = t418 * qJD(4) + t470 * t519;
t535 = qJD(4) * t473;
t453 = qJD(2) * t535;
t534 = qJD(4) * t474;
t522 = t470 * t534;
t484 = t539 * t471 + t522;
t482 = qJD(1) * t484 - t453;
t532 = qJD(5) * t472;
t525 = t472 * t375 + t418 * t532 + t469 * t482;
t327 = t419 * t533 - t525;
t576 = qJ(3) * t471;
t517 = -pkin(1) - t576;
t586 = t474 * t475;
t411 = t517 + t586;
t381 = t411 * qJD(1);
t529 = pkin(3) * t542 + t596;
t386 = t475 * qJD(2) + t529;
t345 = -t381 * t470 + t473 * t386;
t337 = -pkin(9) * t419 + t345;
t450 = qJD(4) + t542;
t331 = pkin(4) * t450 + t337;
t346 = t381 * t473 + t386 * t470;
t338 = pkin(9) * t418 + t346;
t336 = t472 * t338;
t309 = t331 * t469 + t336;
t449 = pkin(2) * t519;
t537 = qJD(3) * t471;
t483 = qJD(2) * t501 - t537;
t370 = qJD(1) * t483 + t449;
t448 = pkin(7) * t518;
t410 = pkin(3) * t518 + t448;
t508 = -t470 * t370 + t473 * t410;
t480 = -t346 * qJD(4) + t508;
t313 = pkin(4) * t518 - pkin(9) * t375 + t480;
t526 = -t473 * t370 - t386 * t535 - t470 * t410;
t488 = -t381 * t536 - t526;
t317 = pkin(9) * t482 + t488;
t512 = t472 * t313 - t469 * t317;
t481 = -t309 * qJD(5) + t512;
t299 = pkin(5) * t518 + qJ(6) * t327 - qJD(6) * t498 + t481;
t507 = t375 * t469 - t472 * t482;
t328 = qJD(5) * t498 + t507;
t502 = t469 * t313 + t472 * t317 + t331 * t532 - t338 * t533;
t300 = -qJ(6) * t328 - qJD(6) * t360 + t502;
t334 = t469 * t338;
t308 = t472 * t331 - t334;
t592 = qJ(6) * t498;
t304 = t308 - t592;
t441 = qJD(5) + t450;
t303 = pkin(5) * t441 + t304;
t305 = t309 - t598;
t496 = -t560 + t564;
t595 = -t299 * t496 + t300 * t420 - t303 * t549 + t305 * t548;
t594 = -0.2e1 * t528;
t467 = t471 ^ 2;
t468 = t474 ^ 2;
t543 = t467 - t468;
t593 = MDP(5) * t543;
t466 = qJD(2) * qJ(3);
t404 = t466 + t427;
t368 = -pkin(4) * t418 + t404;
t333 = pkin(5) * t360 + qJD(6) + t368;
t589 = t333 * t498;
t588 = t368 * t498;
t587 = t441 * t498;
t392 = t496 * t474;
t524 = pkin(4) * t473 + pkin(3);
t540 = qJD(2) * t471;
t369 = (-pkin(7) - t524) * t540 - pkin(4) * t522;
t465 = qJD(2) * qJD(3);
t351 = t453 * pkin(4) + t369 * qJD(1) + t465;
t310 = t328 * pkin(5) + t351;
t436 = t580 * t471;
t423 = t473 * t436;
t515 = pkin(9) * t474 - t411;
t354 = pkin(4) * t471 + t470 * t515 + t423;
t422 = t470 * t436;
t546 = t473 * t411 + t422;
t559 = t473 * t474;
t357 = -pkin(9) * t559 + t546;
t551 = t469 * t354 + t472 * t357;
t585 = t600 * t472;
t431 = t578 * t470;
t545 = -t472 * t431 - t469 * t432;
t584 = -t431 * t533 + t432 * t532 + t600 * t469 + t597 * t472;
t582 = t327 * t496 - t498 * t549;
t544 = pkin(4) * t535 + t524 * t542 + t596;
t577 = qJD(2) * pkin(2);
t574 = t375 * t470;
t573 = t375 * t473;
t426 = t580 * t540;
t388 = -qJD(1) * t426 + t465;
t572 = t388 * t470;
t571 = t388 * t473;
t570 = t404 * t471;
t568 = t419 * t474;
t567 = t450 * t471;
t566 = t450 * t473;
t565 = t450 * t475;
t563 = t470 * t418;
t476 = qJD(2) ^ 2;
t561 = t471 * t476;
t558 = t474 * t476;
t477 = qJD(1) ^ 2;
t557 = t474 * t477;
t452 = t470 * pkin(4) + qJ(3);
t556 = t303 - t304;
t555 = -qJ(6) * t548 - qJD(6) * t420 - t584;
t554 = -pkin(5) * t541 + qJ(6) * t549 - t545 * qJD(5) + qJD(6) * t496 + t469 * t597 - t585;
t553 = t472 * t337 - t334;
t550 = -t441 * t366 - t496 * t518;
t437 = t580 * t474;
t433 = -pkin(2) * t474 + t517;
t405 = qJD(1) * t433;
t538 = qJD(2) * t474;
t531 = t498 * MDP(22);
t527 = t471 * t557;
t402 = pkin(4) * t559 + t437;
t521 = t473 * t534;
t514 = pkin(1) * t594;
t513 = qJD(3) - t577;
t459 = pkin(2) * t540;
t377 = t459 + t483;
t428 = t580 * t538;
t506 = -t377 * t470 + t473 * t428;
t324 = t493 * qJD(2) + (t473 * t515 - t422) * qJD(4) + t506;
t487 = t473 * t377 - t411 * t536 + t470 * t428 + t436 * t535;
t326 = pkin(9) * t484 + t487;
t511 = t472 * t324 - t326 * t469;
t510 = -t337 * t469 - t336;
t509 = t472 * t354 - t357 * t469;
t504 = t431 * t469 - t472 * t432;
t439 = t471 * t518;
t497 = t419 * t473 + t563;
t495 = -0.2e1 * qJD(2) * t405;
t494 = t450 * t470;
t490 = -qJ(3) * t538 - t537;
t379 = qJD(1) * t490 + t449;
t396 = t459 + t490;
t492 = pkin(7) * t476 + qJD(1) * t396 + t379;
t486 = t469 * t324 + t472 * t326 + t354 * t532 - t357 * t533;
t429 = pkin(7) * t519 - t465;
t430 = t455 + t513;
t435 = -t456 - t466;
t478 = -t429 * t474 + (t430 * t474 + (t435 + t456) * t471) * qJD(2);
t454 = pkin(4) * t472 + pkin(5);
t440 = t473 * t518;
t424 = -qJ(3) * t541 + t460;
t393 = t420 * t474;
t385 = t405 * t542;
t353 = -qJ(6) * t420 + t545;
t352 = qJ(6) * t496 + t504;
t340 = t366 * t474 - t496 * t540;
t339 = qJD(2) * t489 + t392 * t583;
t321 = qJ(6) * t392 + t551;
t320 = pkin(5) * t471 + qJ(6) * t393 + t509;
t307 = t553 - t592;
t306 = t510 + t598;
t302 = qJ(6) * t340 + qJD(6) * t392 + t486;
t301 = pkin(5) * t538 - qJ(6) * t339 - qJD(5) * t551 + qJD(6) * t393 + t511;
t1 = [(-t486 * t441 - t502 * t471 + t369 * t498 - t402 * t327 - t351 * t393 + t368 * t339 + (-qJD(1) * t551 - t309) * t538) * MDP(28) + (-t327 * t471 + t339 * t441 + (-qJD(1) * t393 + t498) * t538) * MDP(24) + (t299 * t393 + t300 * t392 - t301 * t498 - t302 * t360 - t303 * t339 + t305 * t340 + t320 * t327 - t321 * t328) * MDP(29) + (t327 * t393 + t339 * t498) * MDP(22) + (-t327 * t392 + t328 * t393 - t339 * t360 + t340 * t498) * MDP(23) + (-t487 * t450 - t426 * t419 + t437 * t375 + ((qJD(2) * t404 + qJD(4) * t381) * t470 + t526) * t471 + (-t404 * t535 - t572 + (-qJD(1) * t546 - t346) * qJD(2)) * t474) * MDP(21) + MDP(6) * t558 + (t497 * t540 + (-t470 * (t473 * t519 - t453) - t573 + (-t473 * t418 + (t419 - t520) * t470) * qJD(4)) * t474) * MDP(16) + (-t474 * t574 + (t471 * t530 - t521) * t419) * MDP(15) + (-t450 * t521 + t375 * t471 + (t568 + (-qJD(1) * t468 + t567) * t470) * qJD(2)) * MDP(17) - MDP(7) * t561 + (pkin(7) * t561 + t474 * t514) * MDP(10) + (-pkin(7) * t558 + t471 * t514) * MDP(9) + (-t328 * t471 + t340 * t441 + (qJD(1) * t392 - t360) * t538) * MDP(25) + (t450 * t538 + t439) * MDP(19) + (t441 * t538 + t439) * MDP(26) + (pkin(7) * t478 + t379 * t433 + t396 * t405) * MDP(14) + t478 * MDP(11) + (t506 * t450 + t426 * t418 + t437 * t453 + ((-qJD(1) * t437 - t404) * t539 + t508) * t471 + (-t346 * t471 - t450 * t546) * qJD(4) + (-t404 * t536 + t345 * qJD(2) + t571 + ((-t411 * t470 + t423) * qJD(2) - t437 * t536) * qJD(1)) * t474) * MDP(20) + (t300 * t321 + t305 * t302 + t299 * t320 + t303 * t301 + t310 * (-pkin(5) * t392 + t402) + (-pkin(5) * t340 + t369) * t333) * MDP(30) + (t450 * t522 + (qJD(4) * t520 - t453) * t471 + (t418 * t474 + (qJD(1) * t543 + t567) * t473) * qJD(2)) * MDP(18) + (t511 * t441 + t512 * t471 + t369 * t360 + t402 * t328 - t351 * t392 - t368 * t340 + (-t309 * t471 - t441 * t551) * qJD(5) + (qJD(1) * t509 + t308) * t538) * MDP(27) + t593 * t594 + (-t471 * t492 + t474 * t495) * MDP(13) + (t471 * t495 + t474 * t492) * MDP(12) + 0.2e1 * MDP(4) * t439; -MDP(4) * t527 + t477 * t593 + ((-t435 - t466) * t471 + (-t430 + t513) * t474) * qJD(1) * MDP(11) + t385 * MDP(12) + (0.2e1 * t465 + (t405 * t474 + t424 * t471) * qJD(1)) * MDP(13) + (-qJ(3) * t429 - qJD(3) * t435 - t405 * t424 + (-t435 * t471 + (-t430 - t577) * t474) * qJD(1) * pkin(7)) * MDP(14) + (-t419 * t494 + t573) * MDP(15) + (-t574 - t473 * t453 - t497 * qJD(4) + (t470 * t521 + (-t563 + (-t419 + t539) * t473) * t471) * qJD(1)) * MDP(16) + (-t450 * t536 + t440 + (-t450 * t562 - t568) * qJD(1)) * MDP(17) + (-t450 * t535 + (-t471 * t566 + (-t418 - t530) * t474) * qJD(1)) * MDP(18) + (qJ(3) * t453 + t572 - t505 * t450 - t529 * t418 + (t404 * t473 - t470 * t565) * qJD(4) + ((-qJ(3) * t536 - t345) * t474 + (t570 + (-t576 + t586) * qJD(2)) * t473) * qJD(1)) * MDP(20) + (qJ(3) * t375 + t571 + t547 * t450 + t529 * t419 + (-t404 * t470 - t473 * t565) * qJD(4) + (t346 * t474 + (-t475 * t538 - t570) * t470) * qJD(1)) * MDP(21) + t582 * MDP(22) + (t327 * t420 + t328 * t496 + t360 * t549 - t498 * t548) * MDP(23) + (-t383 * t441 + t550) * MDP(24) - t548 * t441 * MDP(25) + (t452 * t328 + t351 * t420 + (t431 * t532 + (qJD(5) * t432 + t597) * t469 - t585) * t441 + t548 * t368 + t544 * t360) * MDP(27) + (-t452 * t327 - t351 * t496 - t549 * t368 + t584 * t441 + t544 * t498) * MDP(28) + (t327 * t352 - t328 * t353 - t360 * t555 - t498 * t554 - t595) * MDP(29) + (t300 * t353 + t299 * t352 + t310 * (pkin(5) * t420 + t452) + (pkin(5) * t548 + t544) * t333 + t555 * t305 + t554 * t303) * MDP(30) + (-t450 * MDP(19) - t441 * MDP(26) - t424 * MDP(12) - t498 * MDP(24) + (-qJD(2) * t420 + t360) * MDP(25) + (qJD(2) * t504 - t308) * MDP(27) + (-qJD(2) * t545 + t309) * MDP(28)) * t541 + (MDP(9) * t471 * t477 + MDP(10) * t557) * pkin(1); MDP(12) * t527 + (-t467 * t477 - t476) * MDP(13) + (t385 + t448) * MDP(14) + t440 * MDP(20) + t550 * MDP(27) + (-t328 * t420 - t360 * t548 - t582) * MDP(29) + t595 * MDP(30) + (-t383 * MDP(27) - MDP(28) * t548) * t441 + (-MDP(20) * t494 - MDP(21) * t566) * t450 + (t435 * MDP(14) + t418 * MDP(20) + (-t419 - t520) * MDP(21) - t360 * MDP(27) + (-t420 * t541 - t498) * MDP(28) - t333 * MDP(30)) * qJD(2); -t419 * t418 * MDP(15) + (-t418 ^ 2 + t419 ^ 2) * MDP(16) + (-t418 * t450 + t375) * MDP(17) + (t419 * t450 + t482) * MDP(18) + MDP(19) * t518 + (t346 * t450 - t404 * t419 + t480) * MDP(20) + (t345 * t450 - t404 * t418 - t488) * MDP(21) + t360 * t531 + (t360 * t441 - t327) * MDP(24) + (-t328 + t587) * MDP(25) + (-t510 * t441 - t588 + (-t360 * t419 - t441 * t533 + t472 * t518) * pkin(4) + t481) * MDP(27) + (t553 * t441 + t368 * t360 + (-t419 * t498 - t441 * t532 - t469 * t518) * pkin(4) - t502) * MDP(28) + (-t303 * t360 + t305 * t498 + t306 * t498 + t307 * t360 + t327 * t454 + (-t328 * t469 + (-t360 * t472 + t469 * t498) * qJD(5)) * pkin(4)) * MDP(29) + (-pkin(5) * t589 + t299 * t454 - t303 * t306 - t305 * t307 + (t300 * t469 - t333 * t419 + (-t303 * t469 + t305 * t472) * qJD(5)) * pkin(4)) * MDP(30) + t599; t525 * MDP(24) + (-t507 + t587) * MDP(25) + (t309 * t441 + t512 - t588) * MDP(27) + (t308 * t441 - t502) * MDP(28) + t556 * MDP(30) * t305 + (t327 * MDP(29) + (t299 - t589) * MDP(30)) * pkin(5) + (t441 * MDP(24) + t368 * MDP(28) - MDP(29) * t556 + t531) * t360 + (-MDP(24) * t569 - t498 * MDP(25) - MDP(27) * t309) * qJD(5) + t599; (-t358 - t581) * MDP(29) + (t303 * t498 + t305 * t360 + t310) * MDP(30);];
tauc  = t1;
