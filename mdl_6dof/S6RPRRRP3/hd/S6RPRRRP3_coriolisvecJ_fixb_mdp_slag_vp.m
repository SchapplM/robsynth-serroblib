% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP3
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
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:25
% EndTime: 2019-03-09 06:05:36
% DurationCPUTime: 5.37s
% Computational Cost: add. (5438->466), mult. (12397->585), div. (0->0), fcn. (8132->8), ass. (0->197)
t467 = sin(qJ(4));
t468 = sin(qJ(3));
t544 = qJD(1) * t468;
t518 = t467 * t544;
t469 = cos(qJ(4));
t537 = qJD(3) * t469;
t419 = -t518 + t537;
t539 = qJD(3) * t467;
t420 = t469 * t544 + t539;
t466 = sin(qJ(5));
t582 = cos(qJ(5));
t368 = -t582 * t419 + t420 * t466;
t486 = t466 * t419 + t420 * t582;
t610 = t368 * t486;
t470 = cos(qJ(3));
t563 = t469 * t470;
t491 = pkin(4) * t468 - pkin(9) * t563;
t496 = pkin(3) * t468 - pkin(8) * t470;
t424 = t496 * qJD(1);
t453 = sin(pkin(10)) * pkin(1) + pkin(7);
t433 = t453 * qJD(1);
t590 = qJD(2) * t470 - t468 * t433;
t507 = t469 * t424 - t467 * t590;
t583 = -pkin(9) - pkin(8);
t525 = qJD(4) * t583;
t609 = qJD(1) * t491 - t469 * t525 + t507;
t543 = qJD(1) * t470;
t523 = t467 * t543;
t550 = t467 * t424 + t469 * t590;
t608 = pkin(9) * t523 + t467 * t525 - t550;
t534 = qJD(4) * t468;
t522 = t467 * t534;
t536 = qJD(3) * t470;
t476 = t469 * t536 - t522;
t528 = qJD(3) * qJD(4);
t475 = qJD(1) * t476 + t469 * t528;
t516 = t582 * qJD(5);
t513 = qJD(1) * t534;
t529 = qJD(1) * qJD(3);
t514 = t470 * t529;
t526 = t469 * t513 + (t514 + t528) * t467;
t532 = qJD(5) * t466;
t325 = -t419 * t516 + t420 * t532 + t466 * t526 - t582 * t475;
t451 = -qJD(4) + t543;
t440 = -qJD(5) + t451;
t313 = -t368 * t440 - t325;
t326 = qJD(5) * t486 + t466 * t475 + t582 * t526;
t512 = MDP(23) * t544;
t584 = t486 ^ 2;
t607 = qJD(3) * t512 + (-t440 * t486 - t326) * MDP(22) + MDP(19) * t610 + (-t368 ^ 2 + t584) * MDP(20) + t313 * MDP(21);
t579 = qJD(3) * pkin(3);
t388 = -t590 - t579;
t365 = -pkin(4) * t419 + t388;
t323 = pkin(5) * t368 - qJ(6) * t486 + t365;
t606 = t323 * t368;
t605 = t365 * t368;
t422 = t466 * t469 + t467 * t582;
t588 = qJD(4) + qJD(5);
t374 = t588 * t422;
t499 = t582 * t536;
t519 = t467 * t536;
t346 = t374 * t468 + t466 * t519 - t469 * t499;
t568 = t466 * t467;
t484 = t469 * t582 - t568;
t401 = t484 * t468;
t515 = t468 * t529;
t604 = t346 * t440 + t401 * t515;
t538 = qJD(3) * t468;
t603 = t325 * t470 + t486 * t538;
t589 = t582 * qJD(4) + t516;
t552 = -t469 * t589 + t484 * t543 + t568 * t588;
t551 = -t422 * t543 + t374;
t533 = qJD(4) * t469;
t520 = t468 * t533;
t601 = t519 + t520;
t339 = pkin(5) * t486 + qJ(6) * t368;
t597 = MDP(5) * t468;
t462 = t468 ^ 2;
t596 = MDP(6) * (-t470 ^ 2 + t462);
t578 = t323 * t486;
t595 = t365 * t486;
t594 = t470 * t526;
t566 = t467 * t468;
t347 = t467 * t499 - t466 * t522 - t532 * t566 + (t466 * t536 + t468 * t589) * t469;
t400 = t422 * t468;
t555 = t347 * t440 - t400 * t515;
t436 = t583 * t467;
t437 = t583 * t469;
t485 = t436 * t582 + t466 * t437;
t593 = qJD(5) * t485 - t466 * t609 + t582 * t608;
t381 = t466 * t436 - t437 * t582;
t592 = qJD(5) * t381 + t466 * t608 + t582 * t609;
t461 = t468 * qJD(2);
t397 = t470 * t433 + t461;
t389 = qJD(3) * pkin(8) + t397;
t455 = -cos(pkin(10)) * pkin(1) - pkin(2);
t414 = -pkin(3) * t470 - pkin(8) * t468 + t455;
t392 = t414 * qJD(1);
t567 = t467 * t392;
t354 = t469 * t389 + t567;
t345 = pkin(9) * t419 + t354;
t403 = t469 * t414;
t564 = t469 * t468;
t570 = t453 * t467;
t360 = -pkin(9) * t564 + t403 + (-pkin(4) - t570) * t470;
t423 = t453 * t563;
t547 = t467 * t414 + t423;
t364 = -pkin(9) * t566 + t547;
t591 = t466 * t360 + t582 * t364;
t509 = -t326 * t470 + t368 * t538;
t535 = qJD(4) * t467;
t498 = -t397 + (-t523 + t535) * pkin(4);
t477 = t491 * qJD(3);
t427 = t496 * qJD(3);
t548 = t469 * t427 + t538 * t570;
t332 = t477 + (-t423 + (pkin(9) * t468 - t414) * t467) * qJD(4) + t548;
t521 = t470 * t535;
t549 = t414 * t533 + t467 * t427;
t334 = (-t468 * t537 - t521) * t453 - t601 * pkin(9) + t549;
t586 = -qJD(5) * t591 + t332 * t582 - t466 * t334;
t581 = pkin(8) * t451;
t575 = t388 * t467;
t391 = qJD(3) * t461 + t433 * t536;
t574 = t391 * t467;
t573 = t391 * t469;
t572 = t419 * t468;
t571 = t420 * t451;
t569 = t466 * t345;
t471 = qJD(3) ^ 2;
t565 = t468 * t471;
t562 = t470 * t451;
t561 = t470 * t471;
t353 = -t389 * t467 + t469 * t392;
t344 = -pkin(9) * t420 + t353;
t312 = t344 * t582 - t569;
t560 = -pkin(4) * t516 - qJD(6) + t312;
t559 = qJ(6) * t544 - t593;
t558 = pkin(5) * t544 + t592;
t557 = -t401 * t326 + t346 * t368;
t556 = pkin(5) * t551 + qJ(6) * t552 - qJD(6) * t422 + t498;
t406 = pkin(4) * t566 + t468 * t453;
t434 = qJD(1) * t455;
t541 = qJD(3) * t485;
t540 = qJD(3) * t381;
t531 = t462 * qJD(1);
t337 = -pkin(4) * t451 + t344;
t309 = t337 * t582 - t569;
t530 = qJD(6) - t309;
t375 = pkin(4) * t601 + t453 * t536;
t460 = -pkin(4) * t469 - pkin(3);
t524 = t582 * t345;
t511 = MDP(16) * t538;
t390 = t590 * qJD(3);
t413 = qJD(1) * t427;
t508 = t467 * t390 - t469 * t413;
t506 = t451 + t543;
t505 = -t419 + t537;
t504 = qJD(4) + t543;
t503 = pkin(5) * t515;
t316 = qJD(1) * t477 - qJD(4) * t345 - t508;
t480 = -t389 * t535 + t469 * t390 + t392 * t533 + t467 * t413;
t319 = -pkin(9) * t526 + t480;
t502 = t466 * t316 + t582 * t319 + t337 * t516 - t345 * t532;
t501 = -t582 * t316 + t466 * t319 + t337 * t532 + t345 * t516;
t500 = t451 * t520;
t311 = t466 * t344 + t524;
t497 = pkin(4) * t532 - t311;
t495 = -t451 + t504;
t494 = -t325 * t400 + t347 * t486;
t492 = 0.2e1 * qJD(3) * t434;
t429 = t440 * qJD(6);
t447 = qJ(6) * t515;
t301 = t447 - t429 + t502;
t489 = t360 * t582 - t466 * t364;
t310 = t466 * t337 + t524;
t483 = t504 * t539;
t482 = -t309 * t440 - t502;
t481 = -t310 * t440 - t501;
t479 = t466 * t332 + t582 * t334 + t360 * t516 - t364 * t532;
t478 = (-t531 + t562) * t467;
t356 = pkin(4) * t526 + t391;
t302 = t501 - t503;
t459 = -pkin(4) * t582 - pkin(5);
t454 = pkin(4) * t466 + qJ(6);
t404 = t420 * t538;
t366 = -pkin(5) * t484 - qJ(6) * t422 + t460;
t350 = pkin(5) * t400 - qJ(6) * t401 + t406;
t331 = pkin(4) * t420 + t339;
t328 = t470 * pkin(5) - t489;
t327 = -qJ(6) * t470 + t591;
t308 = -t440 * qJ(6) + t310;
t307 = t440 * pkin(5) + t530;
t306 = pkin(5) * t347 + qJ(6) * t346 - qJD(6) * t401 + t375;
t305 = -pkin(5) * t538 - t586;
t304 = t326 * pkin(5) + t325 * qJ(6) - qJD(6) * t486 + t356;
t303 = qJ(6) * t538 - qJD(6) * t470 + t479;
t1 = [MDP(7) * t561 - MDP(8) * t565 + (t555 - t509) * MDP(22) + (-t440 - t543) * MDP(23) * t538 + (t309 * t538 + t406 * t326 + t365 * t347 + t356 * t400 + t375 * t368 - t440 * t586 + t501 * t470 + t489 * t515) * MDP(24) + (t479 * t440 + t502 * t470 + t375 * t486 - t406 * t325 + t356 * t401 - t365 * t346 + (-qJD(1) * t591 - t310) * t538) * MDP(25) + (t302 * t470 + t304 * t400 + t305 * t440 + t306 * t368 + t323 * t347 + t326 * t350 + (-qJD(1) * t328 - t307) * t538) * MDP(26) + (-t301 * t400 + t302 * t401 - t303 * t368 + t305 * t486 - t307 * t346 - t308 * t347 - t325 * t328 - t326 * t327) * MDP(27) + (-t301 * t470 - t303 * t440 - t304 * t401 - t306 * t486 + t323 * t346 + t325 * t350 + (qJD(1) * t327 + t308) * t538) * MDP(28) + (t301 * t327 + t302 * t328 + t303 * t308 + t304 * t350 + t305 * t307 + t306 * t323) * MDP(29) - 0.2e1 * t529 * t596 + (-t453 * t561 + t468 * t492) * MDP(10) + (t453 * t565 + t470 * t492) * MDP(11) + (t420 * t476 + t475 * t564) * MDP(12) + ((t469 * t419 - t420 * t467) * t536 + ((-t419 + t518) * t535 + (-t420 * qJD(4) - t483 - t526) * t469) * t468) * MDP(13) + (t404 + t506 * t522 + (t531 + (-t451 - t504) * t470) * t537) * MDP(14) + (t500 + t594 + (t478 + t572) * qJD(3)) * MDP(15) - t506 * t511 + (-(-t414 * t535 + t548) * t451 + ((-t453 * t419 + t575) * qJD(3) + (t567 + (t451 * t453 + t389) * t469) * qJD(4) + t508) * t470 + (t453 * t526 + t574 + t388 * t533 + ((-t470 * t570 + t403) * qJD(1) + t353) * qJD(3)) * t468) * MDP(17) + ((-t453 * t521 + t549) * t451 + t480 * t470 + (t573 + (-t453 * t544 - t388) * t535) * t468 + ((t388 * t469 + t453 * t420) * t470 + (t453 * t469 * t495 - qJD(1) * t547 - t354) * t468) * qJD(3)) * MDP(18) + (-t325 * t401 - t346 * t486) * MDP(19) + (-t494 + t557) * MDP(20) + (t603 + t604) * MDP(21) + 0.2e1 * t514 * t597; -MDP(10) * t565 - MDP(11) * t561 + (t500 - t594 + (t478 - t572) * qJD(3)) * MDP(17) + (t404 + (-t451 + t543) * t522 + (-t470 * t495 - t531) * t537) * MDP(18) + (t494 + t557) * MDP(27) + (t301 * t401 + t302 * t400 - t304 * t470 + t307 * t347 - t308 * t346 + t323 * t538) * MDP(29) + (MDP(25) - MDP(28)) * (t603 - t604) + (MDP(24) + MDP(26)) * (t509 + t555); (qJD(3) * t397 - t434 * t544 - t391) * MDP(10) - t434 * t543 * MDP(11) + (-t467 ^ 2 * t513 + (t483 - t571) * t469) * MDP(12) + ((-t526 + t571) * t467 + ((t419 + t537) * qJD(4) + (t470 * t505 - t522) * qJD(1)) * t469) * MDP(13) + (-t451 * t533 + (t469 * t562 + (-t420 + t539) * t468) * qJD(1)) * MDP(14) + (t451 * t535 + (-t467 * t562 + t468 * t505) * qJD(1)) * MDP(15) + t451 * MDP(16) * t544 + (-pkin(3) * t526 - t573 + t507 * t451 + t397 * t419 + (t469 * t581 + t575) * qJD(4) + (-t353 * t468 + (-pkin(8) * t538 - t388 * t470) * t467) * qJD(1)) * MDP(17) + (t574 - t550 * t451 - t397 * t420 + (-t467 * t581 + (t388 - t579) * t469) * qJD(4) + ((-t388 - t579) * t563 + (pkin(3) * t535 - pkin(8) * t537 + t354) * t468) * qJD(1)) * MDP(18) + (-t325 * t422 - t486 * t552) * MDP(19) + (-t325 * t484 - t326 * t422 + t368 * t552 - t486 * t551) * MDP(20) + (t552 * t440 + (qJD(3) * t422 - t486) * t544) * MDP(21) + (t551 * t440 + (qJD(3) * t484 + t368) * t544) * MDP(22) + t440 * t512 + (t460 * t326 - t356 * t484 + t592 * t440 + t498 * t368 + t551 * t365 + (-t309 + t541) * t544) * MDP(24) + (-t460 * t325 + t356 * t422 + t593 * t440 + t498 * t486 - t552 * t365 + (t310 - t540) * t544) * MDP(25) + (-t304 * t484 + t326 * t366 + t558 * t440 + t556 * t368 + t551 * t323 + (t307 + t541) * t544) * MDP(26) + (t301 * t484 + t302 * t422 - t307 * t552 - t308 * t551 + t325 * t485 - t326 * t381 + t368 * t559 + t486 * t558) * MDP(27) + (-t304 * t422 + t325 * t366 + t559 * t440 - t556 * t486 + t552 * t323 + (-t308 + t540) * t544) * MDP(28) + (t301 * t381 - t302 * t485 + t304 * t366 + t307 * t558 - t308 * t559 + t323 * t556) * MDP(29) + (-t470 * t597 + t596) * qJD(1) ^ 2; -t420 * t419 * MDP(12) + (-t419 ^ 2 + t420 ^ 2) * MDP(13) + (t419 * t451 + t475) * MDP(14) + (-t526 - t571) * MDP(15) + qJD(1) * t511 + (-t388 * t420 - t508 + (-qJD(4) - t451) * t354) * MDP(17) + (-t353 * t451 - t388 * t419 - t480) * MDP(18) + (-t311 * t440 - t595 + (-t368 * t420 + t440 * t532 + t515 * t582) * pkin(4) - t501) * MDP(24) + (-t312 * t440 + t605 + (-t420 * t486 + t440 * t516 - t466 * t515) * pkin(4) - t502) * MDP(25) + (-t578 - t331 * t368 + t497 * t440 + (pkin(5) - t459) * t515 - t501) * MDP(26) + (-t325 * t459 - t326 * t454 + (t308 + t497) * t486 + (t307 + t560) * t368) * MDP(27) + (t331 * t486 + t440 * t560 + t454 * t515 + t301 - t606) * MDP(28) + (t301 * t454 + t302 * t459 + t307 * t497 - t308 * t560 - t323 * t331) * MDP(29) + t607; (t481 - t595) * MDP(24) + (t482 + t605) * MDP(25) + (-t339 * t368 + t481 + 0.2e1 * t503 - t578) * MDP(26) + (pkin(5) * t325 - qJ(6) * t326 + (t308 - t310) * t486 + (t307 - t530) * t368) * MDP(27) + (t339 * t486 - 0.2e1 * t429 + 0.2e1 * t447 - t482 - t606) * MDP(28) + (-pkin(5) * t302 + qJ(6) * t301 - t307 * t310 + t308 * t530 - t323 * t339) * MDP(29) + t607; (-t515 + t610) * MDP(26) + t313 * MDP(27) + (-t440 ^ 2 - t584) * MDP(28) + (t308 * t440 + t302 + t578) * MDP(29);];
tauc  = t1;
