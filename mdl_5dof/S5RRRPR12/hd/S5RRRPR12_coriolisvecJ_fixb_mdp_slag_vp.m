% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR12
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
%   see S5RRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:52
% EndTime: 2019-12-31 21:41:08
% DurationCPUTime: 9.41s
% Computational Cost: add. (5663->501), mult. (15387->711), div. (0->0), fcn. (11868->10), ass. (0->195)
t502 = sin(qJ(3));
t505 = cos(qJ(3));
t500 = cos(pkin(5));
t558 = qJD(1) * t500;
t532 = qJD(2) + t558;
t503 = sin(qJ(2));
t498 = sin(pkin(5));
t559 = qJD(1) * t498;
t543 = t503 * t559;
t593 = -t502 * t543 + t505 * t532;
t430 = qJD(5) - t593;
t438 = t502 * t532 + t505 * t543;
t506 = cos(qJ(2));
t557 = qJD(1) * t506;
t542 = t498 * t557;
t483 = -qJD(3) + t542;
t497 = sin(pkin(10));
t499 = cos(pkin(10));
t395 = t438 * t497 + t499 * t483;
t504 = cos(qJ(5));
t397 = t438 * t499 - t483 * t497;
t501 = sin(qJ(5));
t582 = t397 * t501;
t596 = -t504 * t395 - t582;
t598 = t430 * t596;
t547 = pkin(1) * t558;
t458 = pkin(7) * t542 + t503 * t547;
t597 = qJD(4) * t502 + t458 + t483 * (pkin(3) * t502 - qJ(4) * t505);
t518 = t395 * t501 - t397 * t504;
t595 = t430 * t518;
t494 = t498 ^ 2;
t549 = qJD(1) * qJD(2);
t592 = -0.2e1 * t494 * t549;
t591 = MDP(5) * (t503 ^ 2 - t506 ^ 2);
t455 = -pkin(7) * t543 + t506 * t547;
t515 = (pkin(2) * t503 - pkin(8) * t506) * t498;
t456 = qJD(1) * t515;
t564 = t505 * t455 + t502 * t456;
t379 = qJ(4) * t543 + t564;
t554 = qJD(3) * t502;
t546 = pkin(8) * t554;
t569 = t597 * t499 + (-t379 - t546) * t497;
t590 = t499 * t379 + t497 * t597;
t507 = qJD(1) ^ 2;
t589 = pkin(1) * t503;
t588 = pkin(1) * t506;
t587 = pkin(4) * t502;
t586 = pkin(9) + qJ(4);
t537 = t498 * t549;
t526 = t506 * t537;
t517 = t502 * t526;
t401 = qJD(3) * t438 + t517;
t585 = qJ(4) * t401;
t527 = t503 * t537;
t419 = pkin(8) * t532 + t458;
t452 = (-pkin(2) * t506 - pkin(8) * t503 - pkin(1)) * t498;
t429 = qJD(1) * t452;
t457 = qJD(2) * t515;
t445 = qJD(1) * t457;
t574 = t498 * t503;
t488 = pkin(7) * t574;
t459 = (t500 * t588 - t488) * qJD(2);
t446 = qJD(1) * t459;
t553 = qJD(3) * t505;
t529 = t419 * t553 + t429 * t554 - t505 * t445 + t502 * t446;
t332 = -pkin(3) * t527 + t529;
t584 = t332 * t497;
t583 = t332 * t499;
t581 = t593 * t483;
t580 = t593 * t497;
t579 = t438 * t483;
t513 = t483 * t502;
t578 = t483 * t505;
t577 = t494 * t507;
t576 = t497 * t501;
t575 = t497 * t505;
t573 = t498 * t506;
t572 = t499 * t502;
t571 = t499 * t505;
t570 = t505 * t506;
t510 = t419 * t554 - t429 * t553 - t502 * t445 - t505 * t446;
t327 = qJ(4) * t527 - qJD(4) * t483 - t510;
t400 = qJD(3) * t593 + t505 * t526;
t548 = pkin(7) * t573;
t460 = (t500 * t589 + t548) * qJD(2);
t447 = qJD(1) * t460;
t331 = pkin(3) * t401 - qJ(4) * t400 - qJD(4) * t438 + t447;
t308 = t499 * t327 + t497 * t331;
t451 = t548 + (pkin(8) + t589) * t500;
t509 = -t451 * t554 + t452 * t553 + t502 * t457 + t505 * t459;
t556 = qJD(2) * t503;
t339 = (qJ(4) * t556 - qJD(4) * t506) * t498 + t509;
t465 = t500 * t502 + t505 * t574;
t540 = qJD(2) * t573;
t411 = qJD(3) * t465 + t502 * t540;
t464 = -t500 * t505 + t502 * t574;
t412 = -qJD(3) * t464 + t505 * t540;
t345 = pkin(3) * t411 - qJ(4) * t412 - qJD(4) * t465 + t460;
t311 = t499 * t339 + t497 * t345;
t568 = t499 * t546 + t590;
t418 = -pkin(2) * t532 - t455;
t354 = -pkin(3) * t593 - t438 * qJ(4) + t418;
t364 = t505 * t419 + t502 * t429;
t357 = -qJ(4) * t483 + t364;
t321 = t497 * t354 + t499 * t357;
t363 = -t502 * t419 + t429 * t505;
t384 = pkin(3) * t438 - qJ(4) * t593;
t336 = t499 * t363 + t497 * t384;
t450 = t488 + (-pkin(2) - t588) * t500;
t374 = pkin(3) * t464 - qJ(4) * t465 + t450;
t565 = t505 * t451 + t502 * t452;
t375 = -qJ(4) * t573 + t565;
t334 = t497 * t374 + t499 * t375;
t420 = -t499 * t543 + t542 * t575;
t421 = (t497 * t503 + t499 * t570) * t559;
t470 = -t504 * t499 + t576;
t471 = t497 * t504 + t499 * t501;
t552 = qJD(5) * t502;
t567 = t420 * t501 - t421 * t504 - t470 * t553 - t471 * t552;
t551 = qJD(5) * t504;
t566 = -t504 * t420 - t421 * t501 + t471 * t553 + t551 * t572 - t552 * t576;
t563 = t430 * t470;
t562 = t430 * t471;
t533 = -t502 * t455 + t456 * t505;
t380 = -pkin(3) * t543 - t533;
t544 = pkin(4) * t497 + pkin(8);
t561 = -pkin(4) * t420 + t544 * t553 - t380;
t479 = -pkin(3) * t505 - qJ(4) * t502 - pkin(2);
t435 = pkin(8) * t571 + t497 * t479;
t555 = qJD(2) * t505;
t356 = pkin(3) * t483 + qJD(4) - t363;
t550 = -qJD(4) + t356;
t377 = t400 * t497 - t499 * t527;
t378 = t400 * t499 + t497 * t527;
t545 = -t501 * t377 + t504 * t378 - t395 * t551;
t541 = t498 * t556;
t307 = -t327 * t497 + t499 * t331;
t303 = pkin(4) * t401 - pkin(9) * t378 + t307;
t304 = -pkin(9) * t377 + t308;
t536 = t504 * t303 - t304 * t501;
t310 = -t339 * t497 + t499 * t345;
t320 = t499 * t354 - t357 * t497;
t333 = t499 * t374 - t375 * t497;
t335 = -t363 * t497 + t499 * t384;
t535 = t504 * t377 + t501 * t378;
t534 = -t502 * t451 + t452 * t505;
t530 = t494 * t503 * t506 * MDP(4);
t525 = pkin(1) * t592;
t376 = pkin(3) * t573 - t534;
t413 = -pkin(9) * t497 * t502 + t435;
t524 = -pkin(9) * t421 + qJD(5) * t413 + t542 * t587 - (-pkin(9) * t571 + t587) * qJD(3) + t569;
t468 = t499 * t479;
t402 = -pkin(9) * t572 + t468 + (-pkin(8) * t497 - pkin(4)) * t505;
t523 = -pkin(9) * t420 - qJD(5) * t402 - (-pkin(8) * t572 - pkin(9) * t575) * qJD(3) + t590;
t521 = t303 * t501 + t304 * t504;
t309 = -pkin(4) * t593 - pkin(9) * t397 + t320;
t314 = -pkin(9) * t395 + t321;
t300 = t309 * t504 - t314 * t501;
t301 = t309 * t501 + t314 * t504;
t410 = t465 * t499 - t497 * t573;
t316 = pkin(4) * t464 - pkin(9) * t410 + t333;
t409 = t465 * t497 + t499 * t573;
t322 = -pkin(9) * t409 + t334;
t520 = t316 * t504 - t322 * t501;
t519 = t316 * t501 + t322 * t504;
t358 = t504 * t409 + t410 * t501;
t359 = -t409 * t501 + t410 * t504;
t516 = -t451 * t553 - t452 * t554 + t457 * t505 - t502 * t459;
t482 = t586 * t499;
t512 = -pkin(9) * t499 * t593 + pkin(4) * t438 + qJD(4) * t497 + qJD(5) * t482 + t335;
t481 = t586 * t497;
t511 = -pkin(9) * t580 - qJD(4) * t499 + qJD(5) * t481 + t336;
t312 = -qJD(5) * t582 + t545;
t508 = pkin(1) * (-t500 * t549 + t577);
t344 = -pkin(3) * t541 - t516;
t313 = -qJD(5) * t518 + t535;
t493 = -pkin(4) * t499 - pkin(3);
t472 = t544 * t502;
t454 = t470 * t502;
t453 = t471 * t502;
t434 = -pkin(8) * t575 + t468;
t389 = t412 * t499 + t497 * t541;
t388 = t412 * t497 - t499 * t541;
t350 = pkin(4) * t580 + t364;
t349 = pkin(4) * t409 + t376;
t338 = pkin(4) * t395 + t356;
t323 = pkin(4) * t388 + t344;
t318 = qJD(5) * t359 + t504 * t388 + t389 * t501;
t317 = -qJD(5) * t358 - t388 * t501 + t389 * t504;
t315 = pkin(4) * t377 + t332;
t306 = -pkin(9) * t388 + t311;
t305 = pkin(4) * t411 - pkin(9) * t389 + t310;
t299 = -qJD(5) * t301 + t536;
t298 = qJD(5) * t300 + t521;
t1 = [t591 * t592 + (-t447 * t500 - t460 * t532 + t503 * t525) * MDP(9) + (-t446 * t500 - t459 * t532 + t506 * t525) * MDP(10) + (t400 * t465 + t412 * t438) * MDP(11) + (-t400 * t464 - t401 * t465 - t411 * t438 + t412 * t593) * MDP(12) + (-t412 * t483 + (-t400 * t506 + (qJD(1) * t465 + t438) * t556) * t498) * MDP(13) + (t411 * t483 + (t401 * t506 + (-qJD(1) * t464 + t593) * t556) * t498) * MDP(14) + (-t483 * t498 - t494 * t557) * MDP(15) * t556 + (-t516 * t483 - t460 * t593 + t450 * t401 + t447 * t464 + t418 * t411 + (t529 * t506 + (qJD(1) * t534 + t363) * t556) * t498) * MDP(16) + (t509 * t483 + t460 * t438 + t450 * t400 + t447 * t465 + t418 * t412 + (-t510 * t506 + (-qJD(1) * t565 - t364) * t556) * t498) * MDP(17) + (t307 * t464 - t310 * t593 + t320 * t411 + t332 * t409 + t333 * t401 + t344 * t395 + t356 * t388 + t376 * t377) * MDP(18) + (-t308 * t464 + t311 * t593 - t321 * t411 + t332 * t410 - t334 * t401 + t344 * t397 + t356 * t389 + t376 * t378) * MDP(19) + (-t307 * t410 - t308 * t409 - t310 * t397 - t311 * t395 - t320 * t389 - t321 * t388 - t333 * t378 - t334 * t377) * MDP(20) + (t307 * t333 + t308 * t334 + t310 * t320 + t311 * t321 + t332 * t376 + t344 * t356) * MDP(21) + (t312 * t359 - t317 * t518) * MDP(22) + (-t312 * t358 - t313 * t359 + t317 * t596 + t318 * t518) * MDP(23) + (t312 * t464 + t317 * t430 + t359 * t401 - t411 * t518) * MDP(24) + (-t313 * t464 - t318 * t430 - t358 * t401 + t411 * t596) * MDP(25) + (t401 * t464 + t411 * t430) * MDP(26) + ((-qJD(5) * t519 + t305 * t504 - t306 * t501) * t430 + t520 * t401 + t299 * t464 + t300 * t411 - t323 * t596 + t349 * t313 + t315 * t358 + t338 * t318) * MDP(27) + (-t323 * t518 + t349 * t312 + t315 * t359 + t338 * t317 - (qJD(5) * t520 + t305 * t501 + t306 * t504) * t430 - t519 * t401 - t298 * t464 - t301 * t411) * MDP(28) + 0.2e1 * t530 * t549 + (MDP(6) * t540 - MDP(7) * t541) * (qJD(2) + 0.2e1 * t558); (t483 * t554 + (-t506 * t513 + (-t593 + t555) * t503) * t559) * MDP(14) + (-pkin(2) * t401 - t447 * t505 + t533 * t483 + t458 * t593 + (pkin(8) * t578 + t418 * t502) * qJD(3) + (-t363 * t503 + (-pkin(8) * t556 - t418 * t506) * t502) * t559) * MDP(16) + (-pkin(2) * t400 + t447 * t502 - t564 * t483 - t458 * t438 + (-pkin(8) * t513 + t418 * t505) * qJD(3) + (-t418 * t570 + (-pkin(8) * t555 + t364) * t503) * t559) * MDP(17) + (-t356 * t420 - t380 * t395 + t401 * t434 + t569 * t593 + (-t307 + (pkin(8) * t395 + t356 * t497) * qJD(3)) * t505 + (pkin(8) * t377 - t320 * t483 + t584) * t502) * MDP(18) + (-t356 * t421 - t380 * t397 - t401 * t435 - t568 * t593 + (t308 + (pkin(8) * t397 + t356 * t499) * qJD(3)) * t505 + (pkin(8) * t378 + t321 * t483 + t583) * t502) * MDP(19) + (t320 * t421 + t321 * t420 - t377 * t435 - t378 * t434 + (-t307 * t499 - t308 * t497) * t502 + t569 * t397 + t568 * t395 + (-t320 * t499 - t321 * t497) * t553) * MDP(20) + (t307 * t434 + t308 * t435 - t356 * t380 - t568 * t321 - t569 * t320 + (t332 * t502 + t356 * t553) * pkin(8)) * MDP(21) + (-t312 * t454 - t518 * t567) * MDP(22) + (-t312 * t453 + t313 * t454 + t518 * t566 + t567 * t596) * MDP(23) + (-t312 * t505 - t401 * t454 + t430 * t567 + t513 * t518) * MDP(24) + (t313 * t505 - t401 * t453 - t430 * t566 - t513 * t596) * MDP(25) + (-t401 * t505 - t430 * t513) * MDP(26) + ((t402 * t504 - t413 * t501) * t401 - t299 * t505 + t472 * t313 + t315 * t453 + (t501 * t523 - t504 * t524) * t430 - t561 * t596 + t566 * t338 - t300 * t513) * MDP(27) + (t472 * t312 - t315 * t454 - (t402 * t501 + t413 * t504) * t401 + t298 * t505 + (t501 * t524 + t504 * t523) * t430 - t561 * t518 + t567 * t338 + t301 * t513) * MDP(28) + t577 * t591 + (-pkin(7) * t526 + t458 * t532 + t503 * t508) * MDP(9) + (pkin(7) * t527 + t455 * t532 + t506 * t508) * MDP(10) + (t400 * t502 - t438 * t578) * MDP(11) + ((t400 - t581) * t505 + (-t401 + t579) * t502) * MDP(12) + (-t483 * t553 + (t483 * t570 + (qJD(2) * t502 - t438) * t503) * t559) * MDP(13) + t483 * MDP(15) * t543 + (-t530 + (-MDP(6) * t506 + MDP(7) * t503) * t498 * t500) * t507; -t593 ^ 2 * MDP(12) + (t400 + t581) * MDP(13) + (-t517 - t579) * MDP(14) + MDP(15) * t527 + (-t364 * t483 - t529) * MDP(16) + (-t363 * t483 - t418 * t593 + t510) * MDP(17) + (-t497 * t585 - pkin(3) * t377 - t583 - t364 * t395 - (t497 * t550 - t335) * t593) * MDP(18) + (-t499 * t585 - pkin(3) * t378 + t584 - t364 * t397 - (t499 * t550 + t336) * t593) * MDP(19) + (t335 * t397 + t336 * t395 + (-qJ(4) * t377 - qJD(4) * t395 + t320 * t593 + t308) * t499 + (qJ(4) * t378 + qJD(4) * t397 + t321 * t593 - t307) * t497) * MDP(20) + (-pkin(3) * t332 - t320 * t335 - t321 * t336 - t356 * t364 + (-t320 * t497 + t321 * t499) * qJD(4) + (-t307 * t497 + t308 * t499) * qJ(4)) * MDP(21) + (t312 * t471 + t518 * t563) * MDP(22) + (-t312 * t470 - t313 * t471 + t518 * t562 - t563 * t596) * MDP(23) + (t401 * t471 - t430 * t563) * MDP(24) + (-t401 * t470 - t430 * t562) * MDP(25) + ((-t481 * t504 - t482 * t501) * t401 + t493 * t313 + t315 * t470 + t350 * t596 + (t501 * t511 - t504 * t512) * t430 + t562 * t338) * MDP(27) + (t493 * t312 + t315 * t471 - (-t481 * t501 + t482 * t504) * t401 + t350 * t518 + (t501 * t512 + t504 * t511) * t430 - t563 * t338) * MDP(28) + (-MDP(11) * t593 + t438 * MDP(12) - MDP(14) * qJD(3) - t418 * MDP(16) - t320 * MDP(18) + t321 * MDP(19) + MDP(24) * t518 - MDP(25) * t596 - t430 * MDP(26) - t300 * MDP(27) + t301 * MDP(28)) * t438; (-t397 * t593 + t377) * MDP(18) + (t395 * t593 + t378) * MDP(19) + (-t395 ^ 2 - t397 ^ 2) * MDP(20) + (t320 * t397 + t321 * t395 + t332) * MDP(21) + (t313 - t595) * MDP(27) + (t312 + t598) * MDP(28); t518 * t596 * MDP(22) + (t518 ^ 2 - t596 ^ 2) * MDP(23) + (t545 - t598) * MDP(24) + (-t535 - t595) * MDP(25) + t401 * MDP(26) + (t301 * t430 + t338 * t518 + t536) * MDP(27) + (t300 * t430 - t338 * t596 - t521) * MDP(28) + (-MDP(24) * t582 + MDP(25) * t518 - MDP(27) * t301 - t300 * MDP(28)) * qJD(5);];
tauc = t1;
