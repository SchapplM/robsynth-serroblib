% Calculate Coriolis joint torque vector for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR14_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:40
% EndTime: 2024-09-27 18:43:44
% DurationCPUTime: 1.91s
% Computational Cost: add. (4421->308), mult. (9292->429), div. (0->0), fcn. (6791->10), ass. (0->192)
t455 = sin(qJ(4));
t456 = sin(qJ(3));
t459 = cos(qJ(4));
t460 = cos(qJ(3));
t592 = -t455 * t456 + t459 * t460;
t452 = sin(pkin(5));
t591 = t592 * t452;
t461 = cos(qJ(2));
t566 = pkin(1) * qJD(1);
t449 = qJD(1) + qJD(2);
t567 = t449 * pkin(2);
t431 = t461 * t566 + t567;
t453 = cos(pkin(5));
t554 = t453 * t460;
t417 = t431 * t554;
t457 = sin(qJ(2));
t534 = t457 * t566;
t559 = t449 * t452;
t420 = pkin(8) * t559 + t534;
t498 = pkin(9) * t559 + t420;
t374 = -t498 * t456 + t417;
t558 = t449 * t453;
t436 = qJD(3) + t558;
t361 = pkin(3) * t436 + t374;
t555 = t453 * t456;
t530 = t431 * t555;
t375 = t498 * t460 + t530;
t371 = t459 * t375;
t489 = -t361 * t455 - t371;
t585 = t591 * t449;
t568 = pkin(10) * t585;
t321 = -t489 + t568;
t454 = sin(qJ(5));
t537 = qJD(5) * t454;
t319 = t321 * t537;
t458 = cos(qJ(5));
t397 = t458 * t585;
t556 = t452 * t460;
t527 = t449 * t556;
t528 = t456 * t559;
t407 = -t455 * t527 - t459 * t528;
t359 = t407 * t454 + t397;
t557 = t449 * t460;
t404 = (-pkin(3) * t557 - t431) * t452;
t366 = -pkin(4) * t585 + t404;
t590 = -t366 * t359 + t319;
t540 = qJD(3) * t460;
t521 = t453 * t540;
t548 = t460 * t461;
t551 = t456 * t457;
t589 = (-t453 * t551 + t548) * t566 - pkin(2) * t521;
t532 = t452 * (-pkin(8) - pkin(9));
t505 = t456 * t532;
t588 = -qJD(3) * t505 + t589;
t471 = pkin(1) * (-t456 * t461 - t457 * t554);
t410 = qJD(1) * t471;
t535 = pkin(2) * t555;
t587 = (t460 * t532 - t535) * qJD(3) - t410;
t576 = qJD(3) + qJD(4);
t364 = t585 * t576;
t565 = pkin(1) * qJD(2);
t531 = qJD(1) * t565;
t504 = t457 * t531;
t491 = t453 * t504;
t503 = t461 * t531;
t546 = t431 * t521 + t460 * t503;
t341 = (-t498 * qJD(3) - t491) * t456 + t546;
t470 = qJD(2) * t471;
t469 = qJD(1) * t470;
t342 = -qJD(3) * t375 + t469;
t515 = -t455 * t341 + t459 * t342;
t466 = t489 * qJD(4) + t515;
t302 = -pkin(10) * t364 + t466;
t432 = qJD(4) + t436;
t429 = qJD(5) + t432;
t586 = (-t321 * t429 - t302) * t454 + t590;
t486 = t458 * t407 - t454 * t585;
t584 = t366 * t486;
t448 = t452 ^ 2;
t560 = t448 * t460;
t583 = t456 * MDP(7) * t560 - (t456 ^ 2 - t460 ^ 2) * MDP(8) * t448;
t485 = t455 * t460 + t456 * t459;
t574 = t452 * t576;
t381 = t485 * t574;
t365 = t449 * t381;
t313 = qJD(5) * t397 + t458 * t364 - t454 * t365 + t407 * t537;
t464 = t486 * qJD(5) - t364 * t454 - t458 * t365;
t582 = t359 * MDP(21) * t486 + (-t359 ^ 2 + t486 ^ 2) * MDP(22) + (-t359 * t429 + t313) * MDP(23) + (-t429 * t486 + t464) * MDP(24);
t581 = MDP(6) * t461;
t541 = qJD(3) * t456;
t580 = t452 * (-pkin(3) * t541 + t534);
t446 = t453 * pkin(3);
t401 = pkin(2) * t554 + t446 + t505;
t442 = pkin(9) * t556;
t477 = -pkin(8) * t556 - t535;
t409 = t442 - t477;
t487 = -t401 * t455 - t409 * t459;
t578 = t487 * qJD(4) + t588 * t455 + t587 * t459;
t538 = qJD(4) * t459;
t539 = qJD(4) * t455;
t577 = -t401 * t538 + t409 * t539 - t587 * t455 + t588 * t459;
t399 = t407 * pkin(10);
t369 = t455 * t375;
t513 = t459 * t361 - t369;
t320 = t399 + t513;
t536 = qJD(3) - t436;
t575 = t536 * t420 + t491;
t573 = pkin(3) * t429;
t572 = pkin(3) * t460;
t571 = pkin(4) * t407;
t570 = pkin(4) * t591;
t380 = t592 * t574;
t569 = pkin(10) * t380;
t562 = t404 * t407;
t561 = t448 * t456;
t552 = t455 * t458;
t550 = t458 * t321;
t547 = t459 * t374 - t369;
t444 = pkin(1) * t461 + pkin(2);
t545 = t444 * t521 + t548 * t565;
t522 = t452 * t541;
t438 = pkin(3) * t522;
t408 = t449 * t438 + t452 * t504;
t533 = t457 * t565;
t529 = t444 * t555;
t415 = t485 * t452;
t378 = t415 * t458 + t454 * t591;
t323 = t378 * qJD(5) + t380 * t454 + t458 * t381;
t336 = pkin(4) * t365 + t408;
t377 = t415 * t454 - t458 * t591;
t317 = pkin(4) * t432 + t320;
t490 = -t454 * t317 - t550;
t499 = -t459 * t341 - t455 * t342 - t361 * t538 + t375 * t539;
t301 = -pkin(10) * t365 - t499;
t516 = -t454 * t301 + t458 * t302;
t467 = t490 * qJD(5) + t516;
t524 = t366 * t323 + t336 * t377 + t467 * t453;
t523 = t404 * t381 - t408 * t591 + t466 * t453;
t520 = -t431 - t567;
t519 = -pkin(4) * t429 - t317;
t367 = pkin(4) * t381 + t438;
t435 = pkin(1) * t457 + pkin(8) * t452;
t518 = -pkin(9) * t452 - t435;
t517 = t453 * pkin(4) - pkin(10) * t415;
t514 = -((-qJD(3) * t420 - t491) * t456 + t546) * t453 + t504 * t561;
t511 = -t374 * t455 - t371;
t509 = pkin(1) * t449 * t551;
t508 = pkin(3) * t528;
t507 = qJD(5) * t317 + t301;
t506 = t453 * t533;
t425 = (-pkin(2) - t572) * t452;
t497 = (-qJD(2) + t449) * t566;
t496 = (-qJD(1) - t449) * t565;
t322 = -t377 * qJD(5) + t380 * t458 - t381 * t454;
t495 = -(t454 * t302 + t507 * t458 - t319) * t453 + t366 * t322 + t336 * t378;
t494 = t404 * t380 + t408 * t415 + t453 * t499;
t379 = t381 * pkin(10);
t493 = -qJD(5) * (t401 * t459 - t409 * t455 + t517) + t379 + t577;
t412 = t591 * pkin(10);
t492 = qJD(5) * (t412 - t487) + t569 - t578;
t418 = (-t444 - t572) * t452;
t387 = t444 * t554 + t518 * t456 + t446;
t479 = -t435 * t460 - t529;
t390 = t442 - t479;
t488 = -t387 * t455 - t390 * t459;
t483 = qJD(3) * (-t444 * t449 - t431);
t482 = -t452 * t534 + t367;
t481 = t457 * t496;
t480 = t457 * t497;
t478 = t516 + t584;
t476 = -t404 * t585 + t499;
t362 = (t518 * qJD(3) - t506) * t456 + t545;
t363 = t470 + (t518 * t460 - t529) * qJD(3);
t474 = t459 * t362 + t455 * t363 + t387 * t538 - t390 * t539;
t468 = t488 * qJD(4) - t362 * t455 + t459 * t363;
t463 = (-t313 * t377 + t322 * t359 + t323 * t486 + t378 * t464) * MDP(22) + (t313 * t378 - t322 * t486) * MDP(21) + (t313 * t453 + t322 * t429) * MDP(23) + (-t323 * t429 + t453 * t464) * MDP(24) + (t364 * t591 - t365 * t415 + t380 * t585 + t381 * t407) * MDP(15) + (t364 * t415 - t380 * t407) * MDP(14) + (t364 * t453 + t380 * t432) * MDP(16) + (-t365 * t453 - t381 * t432) * MDP(17) + 0.2e1 * t583 * qJD(3) * t449 + (t452 * MDP(9) * t540 - MDP(10) * t522) * (t436 + t558);
t462 = t407 * t585 * MDP(14) + (-t432 * t585 + t364) * MDP(16) + (-t485 * t559 * t576 - t407 * t432) * MDP(17) + (t407 ^ 2 - t585 ^ 2) * MDP(15) + t582;
t443 = pkin(3) * t459 + pkin(4);
t439 = t452 * t533;
t416 = t439 + t438;
t391 = t425 - t570;
t386 = t418 - t570;
t385 = t508 - t571;
t347 = t367 + t439;
t345 = ((-t420 * t460 - t530) * qJD(3) + t469) * t453;
t331 = t412 - t488;
t330 = t387 * t459 - t390 * t455 + t517;
t325 = t399 + t547;
t324 = t511 - t568;
t309 = t468 - t569;
t308 = -t379 + t474;
t1 = [t463 + (t418 * t364 - t416 * t407 - t474 * t432 + t494) * MDP(20) + ((t479 * qJD(3) + t470) * t436 + t345 + (t456 * t483 + t460 * t481) * t448) * MDP(12) + (-(t458 * t308 + t454 * t309 + (t330 * t458 - t331 * t454) * qJD(5)) * t429 - t347 * t486 + t386 * t313 + t495) * MDP(27) + t496 * t581 + ((-t308 * t454 + t309 * t458 + (-t330 * t454 - t331 * t458) * qJD(5)) * t429 - t347 * t359 - t386 * t464 + t524) * MDP(26) + (-((-qJD(3) * t435 - t506) * t456 + t545) * t436 + (qJD(2) * t509 + t460 * t483) * t448 + t514) * MDP(13) + (t418 * t365 - t416 * t585 + t468 * t432 + t523) * MDP(19) + MDP(5) * t481; t463 + (-t391 * t464 + (t493 * t454 - t492 * t458) * t429 - t482 * t359 + t524) * MDP(26) + (t391 * t313 + (t492 * t454 + t493 * t458) * t429 - t482 * t486 + t495) * MDP(27) + (t425 * t365 + t432 * t578 + t580 * t585 + t523) * MDP(19) + (t425 * t364 + t407 * t580 + t432 * t577 + t494) * MDP(20) + ((pkin(8) * t522 + t589) * t436 + (-qJD(1) * t509 + t520 * t540) * t448 + t514) * MDP(13) + t497 * t581 + (-t410 * t436 + t345 + t480 * t560 + (t477 * t436 + t520 * t561) * qJD(3)) * MDP(12) + MDP(5) * t480; t462 + (t448 * t431 * t557 + t417 * t436 + t456 * t575 - t546) * MDP(13) + (t385 * t486 + (t324 * t429 - t302 - (-qJD(4) - qJD(5)) * t455 * t573) * t454 + ((-pkin(3) * t538 - qJD(5) * t443 + t325) * t429 - t507) * t458 + t590) * MDP(27) + (-(t324 * t458 - t325 * t454) * t429 + t385 * t359 + (-t454 * t459 - t552) * qJD(4) * t573 + ((-pkin(3) * t552 - t443 * t454) * t429 + t490) * qJD(5) + t478) * MDP(26) + t536 * MDP(9) * t527 - t536 * MDP(10) * t528 + (t547 * t432 + (t407 * t528 - t432 * t538) * pkin(3) + t476) * MDP(20) + (-t511 * t432 + t585 * t508 + t562 + (-t371 + (-pkin(3) * t432 - t361) * t455) * qJD(4) + t515) * MDP(19) + (-t575 * t460 + (-t503 + (t448 * t449 - t536 * t453) * t431) * t456) * MDP(12) - t583 * t449 ^ 2; (-t489 * t432 + t466 + t562) * MDP(19) + (t513 * t432 + t476) * MDP(20) + (-(-t320 * t454 - t550) * t429 - t359 * t571 + (t519 * t454 - t550) * qJD(5) + t478) * MDP(26) + (-t486 * t571 + (t519 * qJD(5) + t320 * t429 - t301) * t458 + t586) * MDP(27) + t462; (-t490 * t429 + t467 + t584) * MDP(26) + ((-t301 + (-qJD(5) + t429) * t317) * t458 + t586) * MDP(27) + t582;];
tauc = t1;
