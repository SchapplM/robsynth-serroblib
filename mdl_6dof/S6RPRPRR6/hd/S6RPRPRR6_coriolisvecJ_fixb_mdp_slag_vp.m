% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:53:05
% EndTime: 2019-03-09 03:53:19
% DurationCPUTime: 8.09s
% Computational Cost: add. (6862->465), mult. (18041->620), div. (0->0), fcn. (14732->10), ass. (0->186)
t513 = cos(pkin(10));
t594 = cos(qJ(3));
t549 = t594 * t513;
t502 = qJD(1) * t549;
t511 = sin(pkin(10));
t516 = sin(qJ(3));
t574 = t516 * t511;
t548 = qJD(1) * t574;
t471 = -t502 + t548;
t465 = qJD(5) + t471;
t488 = t594 * t511 + t516 * t513;
t473 = t488 * qJD(1);
t510 = sin(pkin(11));
t512 = cos(pkin(11));
t447 = -t512 * qJD(3) + t473 * t510;
t449 = qJD(3) * t510 + t473 * t512;
t515 = sin(qJ(5));
t518 = cos(qJ(5));
t398 = t518 * t447 + t449 * t515;
t517 = cos(qJ(6));
t514 = sin(qJ(6));
t601 = t447 * t515 - t449 * t518;
t583 = t601 * t514;
t345 = t517 * t398 - t583;
t459 = qJD(6) + t465;
t606 = t345 * t459;
t610 = t398 * t465;
t531 = t398 * t514 + t517 * t601;
t605 = t459 * t531;
t572 = t518 * t512;
t485 = t510 * t515 - t572;
t567 = t465 * t485;
t487 = t510 * t518 + t512 * t515;
t475 = t487 * qJD(5);
t566 = t487 * t471 + t475;
t609 = t465 * t601;
t433 = pkin(3) * t473 + qJ(4) * t471;
t592 = pkin(7) + qJ(2);
t493 = t592 * t511;
t489 = qJD(1) * t493;
t495 = t592 * t513;
t490 = qJD(1) * t495;
t527 = t594 * t489 + t516 * t490;
t382 = t510 * t433 - t512 * t527;
t578 = t471 * t510;
t366 = pkin(8) * t578 + t382;
t607 = -qJD(4) * t512 + t366;
t437 = -t485 * t514 + t487 * t517;
t568 = qJD(6) * t437 - t567 * t514 + t566 * t517;
t604 = t488 * qJD(2);
t506 = -pkin(2) * t513 - pkin(1);
t491 = qJD(1) * t506 + qJD(2);
t414 = pkin(3) * t471 - qJ(4) * t473 + t491;
t439 = -t516 * t489 + t594 * t490;
t432 = qJD(3) * qJ(4) + t439;
t371 = t512 * t414 - t432 * t510;
t338 = pkin(4) * t471 - pkin(8) * t449 + t371;
t372 = t510 * t414 + t512 * t432;
t351 = -pkin(8) * t447 + t372;
t323 = t338 * t515 + t351 * t518;
t317 = -pkin(9) * t398 + t323;
t557 = qJD(6) * t514;
t315 = t317 * t557;
t430 = -qJD(3) * pkin(3) + qJD(4) + t527;
t388 = t447 * pkin(4) + t430;
t343 = t398 * pkin(5) + t388;
t603 = t343 * t345 + t315;
t498 = qJD(3) * t502;
t460 = -qJD(3) * t548 + t498;
t558 = qJD(5) * t518;
t581 = t460 * t510;
t356 = -t447 * t558 + t460 * t572 + (-qJD(5) * t449 - t581) * t515;
t477 = t488 * qJD(3);
t461 = qJD(1) * t477;
t397 = pkin(3) * t461 - qJ(4) * t460 - qJD(4) * t473;
t526 = t549 - t574;
t522 = t526 * qJD(2);
t403 = qJD(1) * t522 + (qJD(4) - t527) * qJD(3);
t341 = t512 * t397 - t403 * t510;
t580 = t460 * t512;
t329 = pkin(4) * t461 - pkin(8) * t580 + t341;
t342 = t510 * t397 + t512 * t403;
t332 = -pkin(8) * t581 + t342;
t543 = t518 * t329 - t332 * t515;
t520 = -t323 * qJD(5) + t543;
t306 = pkin(5) * t461 - pkin(9) * t356 + t520;
t357 = -qJD(5) * t601 + t460 * t487;
t559 = qJD(5) * t515;
t525 = t515 * t329 + t518 * t332 + t338 * t558 - t351 * t559;
t307 = -pkin(9) * t357 + t525;
t544 = t517 * t306 - t514 * t307;
t602 = t343 * t531 + t544;
t600 = t461 * MDP(30) + (-t345 ^ 2 + t531 ^ 2) * MDP(27) - t345 * MDP(26) * t531;
t434 = -pkin(3) * t526 - qJ(4) * t488 + t506;
t446 = -t516 * t493 + t495 * t594;
t385 = t512 * t434 - t446 * t510;
t575 = t488 * t512;
t363 = -pkin(4) * t526 - pkin(8) * t575 + t385;
t386 = t510 * t434 + t512 * t446;
t576 = t488 * t510;
t374 = -pkin(8) * t576 + t386;
t570 = t515 * t363 + t518 * t374;
t591 = pkin(8) + qJ(4);
t492 = t591 * t510;
t494 = t591 * t512;
t564 = -t515 * t492 + t518 * t494;
t599 = -t594 * t493 - t516 * t495;
t598 = (t511 ^ 2 + t513 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t436 = t517 * t485 + t487 * t514;
t569 = -qJD(6) * t436 - t514 * t566 - t517 * t567;
t597 = -t437 * t461 - t459 * t569;
t596 = -t461 * t487 + t465 * t567;
t541 = t514 * t356 + t517 * t357;
t313 = -qJD(6) * t531 + t541;
t381 = t512 * t433 + t510 * t527;
t593 = pkin(8) * t512;
t352 = pkin(4) * t473 + t471 * t593 + t381;
t528 = qJD(4) * t510 + qJD(5) * t494;
t595 = t492 * t558 + t607 * t518 + (t352 + t528) * t515;
t466 = t471 ^ 2;
t322 = t518 * t338 - t351 * t515;
t316 = pkin(9) * t601 + t322;
t314 = pkin(5) * t465 + t316;
t589 = t314 * t517;
t588 = t317 * t517;
t587 = t345 * t473;
t586 = t531 * t473;
t585 = t398 * t473;
t584 = t601 * t473;
t476 = t526 * qJD(3);
t577 = t476 * t510;
t409 = pkin(3) * t477 - qJ(4) * t476 - qJD(4) * t488;
t415 = qJD(3) * t599 + t522;
t362 = t510 * t409 + t512 * t415;
t561 = qJD(3) * t516;
t556 = qJD(6) * t517;
t554 = qJD(1) * qJD(2);
t553 = t517 * t356 - t514 * t357 - t398 * t556;
t407 = -pkin(4) * t578 + t439;
t505 = -pkin(4) * t512 - pkin(3);
t547 = qJD(3) * t594;
t546 = pkin(5) * t566 - t407;
t361 = t512 * t409 - t415 * t510;
t335 = pkin(4) * t477 - t476 * t593 + t361;
t340 = -pkin(8) * t577 + t362;
t542 = t518 * t335 - t340 * t515;
t540 = t518 * t363 - t374 * t515;
t539 = -t518 * t492 - t494 * t515;
t538 = qJD(6) * t314 + t307;
t406 = qJD(1) * t604 - t489 * t561 + t490 * t547;
t416 = -t493 * t561 + t495 * t547 + t604;
t537 = -t436 * t461 - t459 * t568;
t536 = -t485 * t461 - t465 * t566;
t350 = t518 * t352;
t422 = -pkin(9) * t485 + t564;
t535 = pkin(5) * t473 - pkin(9) * t567 + t487 * qJD(4) + t564 * qJD(5) + qJD(6) * t422 - t366 * t515 + t350;
t421 = -pkin(9) * t487 + t539;
t534 = pkin(9) * t566 - qJD(6) * t421 + t595;
t380 = pkin(4) * t581 + t406;
t387 = pkin(4) * t577 + t416;
t309 = t314 * t514 + t588;
t532 = -t371 * t510 + t372 * t512;
t426 = t487 * t488;
t427 = t485 * t488;
t377 = t517 * t426 - t427 * t514;
t378 = -t426 * t514 - t427 * t517;
t417 = pkin(4) * t576 - t599;
t524 = t515 * t335 + t518 * t340 + t363 * t558 - t374 * t559;
t312 = t557 * t601 + t553;
t523 = t406 * t488 + t430 * t476 - t460 * t599;
t521 = -pkin(3) * t460 - qJ(4) * t461 + (-qJD(4) + t430) * t471;
t453 = pkin(5) * t485 + t505;
t440 = t461 * t526;
t379 = t426 * pkin(5) + t417;
t376 = t476 * t487 + t558 * t575 - t559 * t576;
t375 = -t475 * t488 - t476 * t485;
t330 = pkin(5) * t376 + t387;
t326 = pkin(5) * t357 + t380;
t325 = -pkin(9) * t426 + t570;
t324 = -pkin(5) * t526 + pkin(9) * t427 + t540;
t320 = qJD(6) * t378 + t375 * t514 + t517 * t376;
t319 = -qJD(6) * t377 + t375 * t517 - t376 * t514;
t311 = -pkin(9) * t376 + t524;
t310 = pkin(5) * t477 - pkin(9) * t375 - qJD(5) * t570 + t542;
t308 = -t317 * t514 + t589;
t1 = [(-t361 * t449 - t362 * t447 + (-t341 * t488 - t371 * t476 - t385 * t460) * t512 + (-t342 * t488 - t372 * t476 - t386 * t460) * t510) * MDP(17) + 0.2e1 * t554 * t598 + (-t312 * t377 - t313 * t378 - t319 * t345 + t320 * t531) * MDP(27) + (t312 * t378 - t319 * t531) * MDP(26) + (-t356 * t426 + t357 * t427 - t375 * t398 + t376 * t601) * MDP(20) + (-t356 * t427 - t375 * t601) * MDP(19) + (t460 * t526 - t461 * t488 - t471 * t476 - t473 * t477) * MDP(9) + (t313 * t526 - t320 * t459 - t345 * t477 - t377 * t461) * MDP(29) + (t357 * t526 - t376 * t465 - t398 * t477 - t426 * t461) * MDP(22) + ((t310 * t517 - t311 * t514) * t459 + (t324 * t517 - t325 * t514) * t461 - t544 * t526 + t308 * t477 + t330 * t345 + t379 * t313 + t326 * t377 + t343 * t320 + ((-t324 * t514 - t325 * t517) * t459 + t309 * t526) * qJD(6)) * MDP(31) + (t342 * t526 - t362 * t471 - t372 * t477 - t386 * t461 + t416 * t449 + t512 * t523) * MDP(16) + (-t341 * t526 + t361 * t471 + t371 * t477 + t385 * t461 + t416 * t447 + t510 * t523) * MDP(15) + (-t309 * t477 + t379 * t312 - t315 * t526 + t343 * t319 + t326 * t378 - t330 * t531 + (-(-qJD(6) * t325 + t310) * t459 - t324 * t461 + t306 * t526) * t514 + (-(qJD(6) * t324 + t311) * t459 - t325 * t461 + t538 * t526) * t517) * MDP(32) + (-t312 * t526 + t319 * t459 + t378 * t461 - t477 * t531) * MDP(28) + (-t323 * t477 + t417 * t356 + t388 * t375 - t380 * t427 - t387 * t601 - t461 * t570 - t465 * t524 + t525 * t526) * MDP(25) + (-t356 * t526 + t375 * t465 - t427 * t461 - t477 * t601) * MDP(21) + (t465 * t477 - t440) * MDP(23) + (t459 * t477 - t440) * MDP(30) + (t341 * t385 + t342 * t386 + t361 * t371 + t362 * t372 - t406 * t599 + t416 * t430) * MDP(18) + (t542 * t465 + t540 * t461 - t543 * t526 + t322 * t477 + t387 * t398 + t417 * t357 + t380 * t426 + t388 * t376 + (t323 * t526 - t465 * t570) * qJD(5)) * MDP(24) + (t461 * t506 + t477 * t491) * MDP(13) + (t460 * t506 + t476 * t491) * MDP(14) + (t460 * t488 + t473 * t476) * MDP(8) + (MDP(10) * t476 - MDP(11) * t477 - MDP(13) * t416 - MDP(14) * t415) * qJD(3); 0.2e1 * t473 * qJD(3) * MDP(13) + (t498 + (-t471 - t548) * qJD(3)) * MDP(14) + (-t447 * t473 + t461 * t512 - t466 * t510) * MDP(15) + (-t449 * t473 - t461 * t510 - t466 * t512) * MDP(16) + ((-t447 * t512 + t449 * t510) * t471 + (-t510 ^ 2 - t512 ^ 2) * t460) * MDP(17) + (t341 * t512 + t342 * t510 - t430 * t473 + t471 * t532) * MDP(18) + (t536 - t585) * MDP(24) + (t584 + t596) * MDP(25) + (t537 - t587) * MDP(31) + (t586 + t597) * MDP(32) - qJD(1) ^ 2 * t598; -t466 * MDP(9) + (t498 + (t471 - t548) * qJD(3)) * MDP(10) + (qJD(3) * t439 - t406) * MDP(13) + (t491 * t471 - t526 * t554) * MDP(14) + (-t381 * t471 - t406 * t512 - t439 * t447 + t510 * t521) * MDP(15) + (t382 * t471 + t406 * t510 - t439 * t449 + t512 * t521) * MDP(16) + (t381 * t449 + t382 * t447 + (-qJD(4) * t447 - t371 * t471 + t342) * t512 + (qJD(4) * t449 - t372 * t471 - t341) * t510) * MDP(17) + (-pkin(3) * t406 - t371 * t381 - t372 * t382 - t430 * t439 + t532 * qJD(4) + (-t341 * t510 + t342 * t512) * qJ(4)) * MDP(18) + (t356 * t487 + t567 * t601) * MDP(19) + (-t356 * t485 - t357 * t487 + t398 * t567 + t566 * t601) * MDP(20) + (t584 - t596) * MDP(21) + (t536 + t585) * MDP(22) + (t539 * t461 + t505 * t357 + t380 * t485 - t407 * t398 + (-t350 - t528 * t518 + (qJD(5) * t492 + t607) * t515) * t465 + t566 * t388) * MDP(24) + (t505 * t356 + t380 * t487 - t567 * t388 + t407 * t601 - t564 * t461 + t465 * t595) * MDP(25) + (t312 * t437 - t531 * t569) * MDP(26) + (-t312 * t436 - t313 * t437 - t345 * t569 + t531 * t568) * MDP(27) + (t586 - t597) * MDP(28) + (t537 + t587) * MDP(29) + ((t421 * t517 - t422 * t514) * t461 + t453 * t313 + t326 * t436 + (t514 * t534 - t517 * t535) * t459 + t546 * t345 + t568 * t343) * MDP(31) + (-(t421 * t514 + t422 * t517) * t461 + t453 * t312 + t326 * t437 + (t514 * t535 + t517 * t534) * t459 - t546 * t531 + t569 * t343) * MDP(32) + (-t491 * MDP(13) - t371 * MDP(15) + t372 * MDP(16) - t465 * MDP(23) - t322 * MDP(24) + t323 * MDP(25) - t459 * MDP(30) - t308 * MDP(31) + t309 * MDP(32) + MDP(8) * t471 + t473 * MDP(9)) * t473; (t449 * t471 + t581) * MDP(15) + (-t447 * t471 + t580) * MDP(16) + (-t447 ^ 2 - t449 ^ 2) * MDP(17) + (t371 * t449 + t372 * t447 + t406) * MDP(18) + (t357 - t609) * MDP(24) + (t356 - t610) * MDP(25) + (t313 - t605) * MDP(31) + (t312 - t606) * MDP(32); -t601 * t398 * MDP(19) + (-t398 ^ 2 + t601 ^ 2) * MDP(20) + (t356 + t610) * MDP(21) + (-t357 - t609) * MDP(22) + t461 * MDP(23) + (t323 * t465 + t388 * t601 + t520) * MDP(24) + (t322 * t465 + t388 * t398 - t525) * MDP(25) + (t312 + t606) * MDP(28) + (-t313 - t605) * MDP(29) + (-(-t316 * t514 - t588) * t459 - t309 * qJD(6) + (t345 * t601 - t459 * t557 + t517 * t461) * pkin(5) + t602) * MDP(31) + ((-t317 * t459 - t306) * t514 + (t316 * t459 - t538) * t517 + (-t459 * t556 - t514 * t461 - t531 * t601) * pkin(5) + t603) * MDP(32) + t600; (t553 + t606) * MDP(28) + (-t541 - t605) * MDP(29) + (t309 * t459 + t602) * MDP(31) + (-t514 * t306 - t517 * t307 + t308 * t459 + t603) * MDP(32) + (MDP(28) * t583 + MDP(29) * t531 - MDP(31) * t309 - MDP(32) * t589) * qJD(6) + t600;];
tauc  = t1;
