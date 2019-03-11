% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:19
% EndTime: 2019-03-09 09:57:28
% DurationCPUTime: 5.57s
% Computational Cost: add. (5430->497), mult. (13350->622), div. (0->0), fcn. (9229->6), ass. (0->199)
t480 = sin(qJ(4));
t477 = sin(pkin(9));
t483 = cos(qJ(2));
t550 = qJD(1) * t483;
t531 = t477 * t550;
t478 = cos(pkin(9));
t482 = cos(qJ(4));
t567 = t482 * t478;
t532 = t483 * t567;
t542 = qJD(4) * t482;
t543 = qJD(4) * t480;
t582 = -t477 * t543 + t478 * t542;
t555 = -qJD(1) * t532 + t480 * t531 + t582;
t481 = sin(qJ(2));
t551 = qJD(1) * t481;
t526 = t477 * t551;
t546 = qJD(2) * t478;
t431 = -t526 + t546;
t530 = t478 * t551;
t547 = qJD(2) * t477;
t432 = t530 + t547;
t383 = -t482 * t431 + t432 * t480;
t592 = t383 ^ 2;
t465 = -qJD(4) + t550;
t591 = t383 * t465;
t508 = pkin(2) * t481 - qJ(3) * t483;
t437 = t508 * qJD(1);
t401 = pkin(7) * t526 + t478 * t437;
t569 = t478 * t483;
t503 = pkin(3) * t481 - pkin(8) * t569;
t372 = qJD(1) * t503 + t401;
t422 = t477 * t437;
t570 = t478 * t481;
t571 = t477 * t483;
t497 = -pkin(7) * t570 - pkin(8) * t571;
t388 = qJD(1) * t497 + t422;
t577 = pkin(8) + qJ(3);
t446 = t577 * t477;
t447 = t577 * t478;
t590 = qJD(3) * t567 - t482 * t388 - t446 * t542 + (-qJD(3) * t477 - qJD(4) * t447 - t372) * t480;
t435 = t477 * t482 + t478 * t480;
t421 = t435 * qJD(4);
t586 = t435 * t483;
t554 = -qJD(1) * t586 + t421;
t442 = -pkin(2) * t483 - qJ(3) * t481 - pkin(1);
t424 = t442 * qJD(1);
t471 = pkin(7) * t550;
t449 = qJD(2) * qJ(3) + t471;
t391 = t478 * t424 - t449 * t477;
t350 = -pkin(3) * t550 - pkin(8) * t432 + t391;
t392 = t477 * t424 + t478 * t449;
t358 = pkin(8) * t431 + t392;
t321 = -t482 * t350 + t358 * t480;
t386 = t431 * t480 + t432 * t482;
t500 = pkin(5) * t386 + t321;
t539 = qJD(5) + t500;
t589 = MDP(24) + MDP(27);
t537 = qJD(1) * qJD(2);
t588 = -0.2e1 * t537;
t587 = MDP(5) * (t481 ^ 2 - t483 ^ 2);
t562 = qJ(5) * t551 - t590;
t399 = -t446 * t480 + t447 * t482;
t585 = qJD(3) * t435 + qJD(4) * t399 + t372 * t482 - t480 * t388;
t425 = pkin(3) * t531 + t471;
t583 = qJ(5) * t555 + qJD(5) * t435 + t425;
t493 = t386 * qJD(4);
t581 = MDP(20) - MDP(23);
t535 = MDP(22) + MDP(26);
t580 = -MDP(21) + t589;
t524 = t483 * t537;
t514 = t477 * t524;
t515 = qJD(2) * t532;
t556 = qJD(1) * t515 + t431 * t542;
t341 = t480 * (qJD(4) * t432 + t514) - t556;
t460 = t465 ^ 2;
t579 = pkin(5) * t383;
t578 = pkin(4) + qJ(6);
t576 = qJD(2) * pkin(2);
t489 = qJD(2) * t586;
t342 = qJD(1) * t489 + t493;
t575 = qJ(5) * t342;
t574 = qJ(5) * t383;
t573 = t383 * t386;
t572 = t477 * t481;
t484 = qJD(2) ^ 2;
t568 = t481 * t484;
t566 = t483 * t484;
t485 = qJD(1) ^ 2;
t565 = t483 * t485;
t525 = t578 * t481;
t564 = pkin(5) * t555 + qJD(1) * t525 + t585;
t563 = pkin(5) * t554 + t562;
t561 = -pkin(4) * t551 - t585;
t434 = t477 * t480 - t567;
t560 = qJD(6) * t434 + t554 * t578 - t583;
t559 = -pkin(4) * t554 + t583;
t322 = t480 * t350 + t482 * t358;
t430 = t478 * t442;
t390 = -pkin(8) * t570 + t430 + (-pkin(7) * t477 - pkin(3)) * t483;
t463 = pkin(7) * t569;
t406 = t477 * t442 + t463;
t397 = -pkin(8) * t572 + t406;
t557 = t480 * t390 + t482 * t397;
t418 = qJD(2) * t508 - qJD(3) * t481;
t409 = t418 * qJD(1);
t470 = pkin(7) * t551;
t440 = (qJD(3) - t470) * qJD(2);
t371 = t477 * t409 + t478 * t440;
t545 = qJD(2) * t481;
t533 = pkin(7) * t545;
t395 = t478 * t418 + t477 * t533;
t464 = pkin(7) * t524;
t417 = pkin(3) * t514 + t464;
t544 = qJD(2) * t483;
t472 = pkin(7) * t544;
t529 = t477 * t544;
t426 = pkin(3) * t529 + t472;
t438 = pkin(3) * t572 + t481 * pkin(7);
t398 = t446 * t482 + t447 * t480;
t549 = qJD(2) * t398;
t548 = qJD(2) * t399;
t541 = t481 * MDP(19);
t540 = -qJD(5) - t321;
t313 = t322 - t579;
t538 = -qJD(6) - t313;
t534 = pkin(7) * t571;
t468 = -pkin(3) * t478 - pkin(2);
t469 = t481 * t537;
t523 = MDP(28) + t581;
t522 = pkin(1) * t588;
t521 = qJD(3) - t576;
t518 = t390 * t482 - t480 * t397;
t370 = t478 * t409 - t440 * t477;
t492 = t503 * qJD(2);
t345 = qJD(1) * t492 + t370;
t357 = -pkin(8) * t514 + t371;
t517 = t480 * t345 + t350 * t542 + t482 * t357 - t358 * t543;
t516 = t482 * t345 - t350 * t543 - t480 * t357 - t358 * t542;
t513 = qJD(2) * t525;
t441 = t470 + t521;
t512 = -t441 + t521;
t317 = qJ(5) * t465 - t322;
t331 = qJ(5) * t483 - t557;
t415 = -t480 * t572 + t481 * t567;
t510 = -qJ(5) * t415 + t438;
t332 = t483 * pkin(4) - t518;
t507 = t523 * t480;
t451 = qJD(5) * t465;
t506 = -t451 + t517;
t502 = -qJ(5) * t435 + t468;
t362 = t492 + t395;
t410 = t477 * t418;
t373 = qJD(2) * t497 + t410;
t501 = t362 * t482 - t480 * t373 - t390 * t543 - t397 * t542;
t499 = -pkin(5) * t342 + t517;
t498 = pkin(5) * t341 + t516;
t494 = t480 * t362 + t482 * t373 + t390 * t542 - t397 * t543;
t400 = -pkin(3) * t431 + t441;
t366 = t421 * t481 + t480 * t529 - t515;
t491 = qJ(5) * t366 - qJD(5) * t415 + t426;
t461 = qJ(5) * t469;
t300 = -t451 + t461 + t499;
t303 = -pkin(4) * t469 - t516;
t488 = -qJ(5) * t386 + t400;
t307 = pkin(4) * t342 + qJ(5) * t341 - qJD(5) * t386 + t417;
t487 = t480 * t580 + t523 * t482;
t301 = t342 * qJ(6) + t383 * qJD(6) + t307;
t486 = -qJD(1) * t513 - t498;
t308 = -qJ(5) * t545 + qJD(5) * t483 - t494;
t323 = -t341 - t591;
t456 = 0.2e1 * t461;
t414 = t435 * t481;
t405 = t430 - t534;
t402 = -pkin(7) * t530 + t422;
t396 = -t478 * t533 + t410;
t380 = pkin(4) * t434 + t502;
t367 = t481 * t582 + t489;
t365 = -pkin(5) * t434 + t399;
t364 = pkin(5) * t435 + t398;
t355 = t434 * t578 + t502;
t351 = pkin(4) * t414 + t510;
t334 = pkin(4) * t386 + t574;
t333 = t414 * t578 + t510;
t326 = -pkin(5) * t414 - t331;
t325 = pkin(4) * t383 + t488;
t324 = pkin(5) * t415 + qJ(6) * t483 + t332;
t320 = t386 * t578 + t574;
t316 = pkin(4) * t465 - t540;
t315 = pkin(4) * t367 + t491;
t314 = t383 * t578 + t488;
t311 = qJD(6) - t317 - t579;
t310 = t465 * t578 + t539;
t309 = -pkin(4) * t545 - t501;
t306 = qJD(6) * t414 + t367 * t578 + t491;
t305 = -pkin(5) * t367 - t308;
t304 = -pkin(5) * t366 + qJD(6) * t483 - t501 - t513;
t302 = -t461 - t506;
t299 = qJD(6) * t465 + t486;
t1 = [(t370 * t405 + t371 * t406 + t391 * t395 + t392 * t396 + (t441 + t470) * t472) * MDP(14) + ((qJD(1) * t396 + t371) * t483 + ((pkin(7) * t432 + t441 * t478) * t483 + (-t392 + (-t406 + 0.2e1 * t463) * qJD(1)) * t481) * qJD(2)) * MDP(12) + t587 * t588 + (-t395 * t432 + t396 * t431 + (-t370 * t478 - t371 * t477) * t481 + (-t391 * t478 - t392 * t477 + (-t405 * t478 - t406 * t477) * qJD(1)) * t544) * MDP(13) + (t341 * t483 + t366 * t465 + (qJD(1) * t415 + t386) * t545) * MDP(17) + (t342 * t483 + t367 * t465 + (-qJD(1) * t414 - t383) * t545) * MDP(18) + (-t501 * t465 - t516 * t483 + t426 * t383 + t438 * t342 + t417 * t414 + t400 * t367 + (qJD(1) * t518 - t321) * t545) * MDP(20) + (t302 * t483 - t307 * t415 + t308 * t465 - t315 * t386 + t325 * t366 + t341 * t351 + (-qJD(1) * t331 - t317) * t545) * MDP(24) + (-t303 * t483 - t307 * t414 - t309 * t465 - t315 * t383 - t325 * t367 - t342 * t351 + (qJD(1) * t332 + t316) * t545) * MDP(23) + (-t300 * t483 - t301 * t415 - t305 * t465 - t306 * t386 + t314 * t366 + t333 * t341 + (qJD(1) * t326 + t311) * t545) * MDP(27) + (t299 * t483 + t301 * t414 + t304 * t465 + t306 * t383 + t314 * t367 + t333 * t342 + (-qJD(1) * t324 - t310) * t545) * MDP(28) - MDP(7) * t568 + (pkin(7) * t568 + t483 * t522) * MDP(10) + (t494 * t465 + t517 * t483 + t426 * t386 - t438 * t341 + t417 * t415 - t400 * t366 + (-qJD(1) * t557 - t322) * t545) * MDP(21) + (-pkin(7) * t566 + t481 * t522) * MDP(9) + ((-qJD(1) * t395 - t370) * t483 + ((-pkin(7) * t431 + t441 * t477) * t483 + (t391 + (t405 + 0.2e1 * t534) * qJD(1)) * t481) * qJD(2)) * MDP(11) + (-t465 - t550) * qJD(2) * t541 + (t302 * t414 + t303 * t415 + t308 * t383 + t309 * t386 - t316 * t366 + t317 * t367 + t331 * t342 - t332 * t341) * MDP(22) + (t341 * t414 - t342 * t415 + t366 * t383 - t367 * t386) * MDP(16) + (t299 * t415 - t300 * t414 + t304 * t386 - t305 * t383 - t310 * t366 - t311 * t367 - t324 * t341 - t326 * t342) * MDP(26) + (-t341 * t415 - t366 * t386) * MDP(15) + MDP(6) * t566 + 0.2e1 * t483 * MDP(4) * t469 + (t299 * t324 + t300 * t326 + t301 * t333 + t304 * t310 + t305 * t311 + t306 * t314) * MDP(29) + (t302 * t331 + t303 * t332 + t307 * t351 + t308 * t317 + t309 * t316 + t315 * t325) * MDP(25); -t481 * MDP(4) * t565 + t485 * t587 + (t401 * t432 - t402 * t431 + (qJD(3) * t431 + t391 * t550 + t371) * t478 + (qJD(3) * t432 + t392 * t550 - t370) * t477) * MDP(13) + (-t391 * t401 - t392 * t402 + (-t391 * t477 + t392 * t478) * qJD(3) + (-t370 * t477 + t371 * t478) * qJ(3) + (-t441 - t576) * t471) * MDP(14) + (-t341 * t435 + t386 * t555) * MDP(15) + (t341 * t434 - t342 * t435 - t383 * t555 - t386 * t554) * MDP(16) + (-t555 * t465 + (qJD(2) * t435 - t386) * t551) * MDP(17) + (t554 * t465 + (-qJD(2) * t434 + t383) * t551) * MDP(18) + (t468 * t342 - t425 * t383 + t417 * t434 + t585 * t465 + t554 * t400 + (t321 - t549) * t551) * MDP(20) + (-t468 * t341 - t425 * t386 + t417 * t435 + t590 * t465 + t555 * t400 + (t322 - t548) * t551) * MDP(21) + (t302 * t434 + t303 * t435 + t316 * t555 + t317 * t554 - t341 * t398 - t342 * t399 + t383 * t562 - t561 * t386) * MDP(22) + (-t307 * t434 - t342 * t380 + t561 * t465 + t559 * t383 - t554 * t325 + (-t316 + t549) * t551) * MDP(23) + (-t307 * t435 + t341 * t380 + t562 * t465 + t559 * t386 - t555 * t325 + (t317 + t548) * t551) * MDP(24) + (-t302 * t399 + t303 * t398 + t307 * t380 - t316 * t561 + t317 * t562 - t325 * t559) * MDP(25) + (t299 * t435 - t300 * t434 + t310 * t555 - t311 * t554 - t341 * t364 - t342 * t365 + t383 * t563 + t386 * t564) * MDP(26) + (-t301 * t435 + t341 * t355 + t563 * t465 - t560 * t386 - t555 * t314 + (qJD(2) * t365 - t311) * t551) * MDP(27) + (t301 * t434 + t342 * t355 + t564 * t465 + t560 * t383 + t554 * t314 + (-qJD(2) * t364 + t310) * t551) * MDP(28) + (t299 * t364 + t300 * t365 + t301 * t355 + t310 * t564 - t311 * t563 + t314 * t560) * MDP(29) + (MDP(9) * t481 * t485 + MDP(10) * t565) * pkin(1) + (((-qJ(3) * t547 - t391) * t481 + (t401 + (t431 - t546) * pkin(7) + t512 * t477) * t483) * MDP(11) + ((-qJ(3) * t546 + t392) * t481 + (-t402 + (-t432 + t547) * pkin(7) + t512 * t478) * t483) * MDP(12) + t465 * t541) * qJD(1); (-t431 ^ 2 - t432 ^ 2) * MDP(13) + (t391 * t432 - t392 * t431 + t464) * MDP(14) + (-t317 * t383 + t307) * MDP(25) + (t311 * t383 + t301) * MDP(29) - t535 * t592 + (-t316 * MDP(25) - t310 * MDP(29) - t386 * t535 - t465 * t523) * t386 + (t431 * t507 + t432 * t487) * qJD(4) + (-t432 * MDP(11) - t431 * MDP(12) + ((MDP(12) + t507) * t478 + (MDP(11) + t487) * t477) * qJD(2)) * t550 + t580 * (-t556 - t591); t323 * MDP(17) - t517 * MDP(21) + (pkin(4) * t341 - t575) * MDP(22) + (t456 + t506) * MDP(24) + (-pkin(4) * t303 - qJ(5) * t302 - t316 * t322 + t317 * t540 - t325 * t334) * MDP(25) + (t341 * t578 - t575) * MDP(26) + (-0.2e1 * t451 + t456 + t499) * MDP(27) + t498 * MDP(28) + (qJ(5) * t300 - t299 * t578 + t310 * t538 + t311 * t539 - t314 * t320) * MDP(29) + (t321 * MDP(21) + t540 * MDP(24) - t500 * MDP(27) + (-0.2e1 * qJD(6) - t313) * MDP(28) - t581 * t322) * t465 - MDP(18) * t493 + (-t465 * MDP(18) - t400 * MDP(20) + (-t317 - t322) * MDP(22) + t325 * MDP(23) + t334 * MDP(24) + (t311 + t538) * MDP(26) + t320 * MDP(27) - t314 * MDP(28) + MDP(16) * t386) * t386 + (-MDP(18) * t586 + (-0.2e1 * pkin(4) * MDP(23) + 0.2e1 * t578 * MDP(28) + MDP(19)) * t481) * t537 + (t386 * MDP(15) + t400 * MDP(21) + (t316 + t540) * MDP(22) + t334 * MDP(23) - t325 * MDP(24) + (t310 - t539) * MDP(26) - t314 * MDP(27) - t320 * MDP(28) - MDP(16) * t383) * t383 + t581 * t516; (-t317 * t465 + t325 * t386 + t303) * MDP(25) + (t314 * t386 + (qJD(6) + t311) * t465 + t486) * MDP(29) + (-MDP(23) + MDP(28)) * (-t469 + t573) + t535 * t323 + t589 * (-t386 ^ 2 - t460); (-t386 * t465 - t537 * t586 - t493) * MDP(26) + (t469 + t573) * MDP(27) + (-t460 - t592) * MDP(28) + (-t310 * t465 - t314 * t383 + t300) * MDP(29);];
tauc  = t1;
