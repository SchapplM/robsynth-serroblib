% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:44
% EndTime: 2019-03-09 11:50:54
% DurationCPUTime: 5.66s
% Computational Cost: add. (7293->438), mult. (18267->578), div. (0->0), fcn. (13651->8), ass. (0->200)
t500 = sin(qJ(4));
t554 = qJD(4) * t500;
t497 = sin(pkin(10));
t501 = sin(qJ(2));
t551 = t501 * qJD(1);
t498 = cos(pkin(10));
t503 = cos(qJ(2));
t573 = t498 * t503;
t461 = qJD(1) * t573 - t497 * t551;
t579 = t461 * t500;
t614 = t554 - t579;
t473 = t497 * t503 + t498 * t501;
t463 = t473 * qJD(1);
t401 = pkin(2) * t551 + pkin(3) * t463 - pkin(8) * t461;
t502 = cos(qJ(4));
t394 = t502 * t401;
t590 = -qJ(3) - pkin(7);
t482 = t590 * t503;
t478 = qJD(1) * t482;
t466 = t497 * t478;
t481 = t590 * t501;
t477 = qJD(1) * t481;
t417 = t477 * t498 + t466;
t489 = pkin(2) * t497 + pkin(8);
t591 = pkin(9) + t489;
t535 = qJD(4) * t591;
t613 = pkin(4) * t463 - t417 * t500 + t394 + (-pkin(9) * t461 + t535) * t502;
t559 = t500 * t401 + t502 * t417;
t612 = -pkin(9) * t579 + t500 * t535 + t559;
t550 = qJD(1) * qJD(2);
t537 = t503 * t550;
t538 = t501 * t550;
t520 = -t497 * t538 + t498 * t537;
t609 = qJD(2) * qJD(4) + t520;
t386 = -t463 * t554 + t502 * t609;
t427 = qJD(2) * t502 - t500 * t463;
t428 = qJD(2) * t500 + t463 * t502;
t499 = sin(qJ(5));
t592 = cos(qJ(5));
t539 = t592 * qJD(5);
t553 = qJD(4) * t502;
t546 = t463 * t553 + t500 * t609;
t552 = qJD(5) * t499;
t334 = -t592 * t386 - t427 * t539 + t428 * t552 + t499 * t546;
t519 = t499 * t427 + t428 * t592;
t335 = qJD(5) * t519 + t499 * t386 + t592 * t546;
t369 = -t592 * t427 + t428 * t499;
t367 = t369 ^ 2;
t455 = qJD(4) - t461;
t448 = qJD(5) + t455;
t462 = t473 * qJD(2);
t449 = qJD(1) * t462;
t593 = t519 ^ 2;
t611 = t449 * MDP(24) + (t448 * t519 - t335) * MDP(23) + t369 * MDP(20) * t519 + (t369 * t448 - t334) * MDP(22) + (-t367 + t593) * MDP(21);
t610 = t369 * qJ(6);
t574 = t498 * t478;
t416 = t497 * t477 - t574;
t526 = t614 * pkin(4) - t416;
t544 = t592 * t500;
t476 = t499 * t502 + t544;
t597 = qJD(4) + qJD(5);
t421 = t597 * t476;
t561 = t476 * t461 - t421;
t572 = t499 * t500;
t518 = t592 * t502 - t572;
t598 = t592 * qJD(4) + t539;
t560 = t518 * t461 - t502 * t598 + t572 * t597;
t472 = t497 * t501 - t573;
t465 = t472 * qJD(2);
t541 = t473 * t553;
t608 = -t465 * t500 + t541;
t589 = qJD(2) * pkin(2);
t469 = t477 + t589;
t413 = t469 * t498 + t466;
t408 = -qJD(2) * pkin(3) - t413;
t366 = -pkin(4) * t427 + t408;
t545 = -pkin(2) * t503 - pkin(1);
t523 = t545 * qJD(1);
t479 = qJD(3) + t523;
t390 = -pkin(3) * t461 - pkin(8) * t463 + t479;
t414 = t497 * t469 - t574;
t409 = qJD(2) * pkin(8) + t414;
t361 = t390 * t500 + t409 * t502;
t487 = pkin(2) * t538;
t389 = t449 * pkin(3) - pkin(8) * t520 + t487;
t380 = t502 * t389;
t536 = qJD(2) * t590;
t459 = qJD(3) * t503 + t501 * t536;
t434 = t459 * qJD(1);
t460 = -qJD(3) * t501 + t503 * t536;
t435 = t460 * qJD(1);
t385 = t434 * t498 + t435 * t497;
t510 = -qJD(4) * t361 - t385 * t500 + t380;
t320 = pkin(4) * t449 - pkin(9) * t386 + t510;
t515 = t502 * t385 + t500 * t389 + t390 * t553 - t409 * t554;
t326 = -pkin(9) * t546 + t515;
t360 = t502 * t390 - t409 * t500;
t351 = -pkin(9) * t428 + t360;
t338 = pkin(4) * t455 + t351;
t352 = pkin(9) * t427 + t361;
t527 = t499 * t320 + t592 * t326 + t338 * t539 - t352 * t552;
t607 = t366 * t369 - t527;
t605 = -0.2e1 * t550;
t604 = MDP(4) * t501;
t603 = MDP(5) * (t501 ^ 2 - t503 ^ 2);
t602 = qJ(6) * t519;
t343 = pkin(5) * t369 + qJD(6) + t366;
t601 = t343 * t519;
t412 = pkin(3) * t472 - pkin(8) * t473 + t545;
t407 = t502 * t412;
t426 = t481 * t497 - t482 * t498;
t575 = t473 * t502;
t357 = pkin(4) * t472 - pkin(9) * t575 - t426 * t500 + t407;
t419 = t502 * t426;
t558 = t500 * t412 + t419;
t576 = t473 * t500;
t364 = -pkin(9) * t576 + t558;
t562 = t499 * t357 + t592 * t364;
t470 = t591 * t500;
t471 = t591 * t502;
t557 = -t499 * t470 + t592 * t471;
t600 = -qJD(5) * t557 + t612 * t499 - t613 * t592;
t599 = t470 * t539 + t471 * t552 + t613 * t499 + t612 * t592;
t349 = t592 * t352;
t323 = t499 * t338 + t349;
t508 = -qJD(5) * t323 + t592 * t320 - t499 * t326;
t596 = -t366 * t519 + t508;
t595 = t334 * t518 - t519 * t561;
t587 = t369 * t463;
t586 = t519 * t463;
t585 = t386 * t500;
t584 = t427 * t455;
t583 = t427 * t463;
t582 = t428 * t455;
t581 = t428 * t463;
t580 = t449 * t476;
t577 = t465 * t502;
t347 = t499 * t352;
t571 = t500 * t449;
t504 = qJD(2) ^ 2;
t570 = t501 * t504;
t431 = t502 * t449;
t569 = t503 * t504;
t505 = qJD(1) ^ 2;
t568 = t503 * t505;
t322 = t592 * t338 - t347;
t312 = t322 - t602;
t311 = pkin(5) * t448 + t312;
t567 = t311 - t312;
t566 = qJ(6) * t561 + qJD(6) * t518 - t599;
t565 = -pkin(5) * t463 + qJ(6) * t560 - t476 * qJD(6) + t600;
t564 = t592 * t351 - t347;
t548 = t501 * t589;
t490 = -pkin(2) * t498 - pkin(3);
t542 = t473 * t554;
t534 = pkin(1) * t605;
t532 = -t351 * t499 - t349;
t531 = t592 * t357 - t364 * t499;
t384 = t497 * t434 - t498 * t435;
t399 = t459 * t497 - t498 * t460;
t529 = -t592 * t470 - t471 * t499;
t425 = -t498 * t481 - t482 * t497;
t528 = t455 * t502;
t525 = -t476 * t335 + t369 * t560;
t524 = t448 * t561 + t518 * t449;
t398 = pkin(4) * t576 + t425;
t522 = t384 * t473 - t426 * t449;
t480 = -pkin(4) * t502 + t490;
t365 = pkin(4) * t608 + t399;
t521 = -t614 * t455 + t431;
t516 = -t542 - t577;
t400 = t459 * t498 + t460 * t497;
t402 = pkin(3) * t462 + pkin(8) * t465 + t548;
t514 = t502 * t400 + t500 * t402 + t412 * t553 - t426 * t554;
t395 = t502 * t402;
t330 = pkin(9) * t577 + pkin(4) * t462 - t400 * t500 + t395 + (-t419 + (pkin(9) * t473 - t412) * t500) * qJD(4);
t332 = -pkin(9) * t608 + t514;
t513 = t499 * t330 + t592 * t332 + t357 * t539 - t364 * t552;
t511 = t408 * t455 - t489 * t449;
t358 = pkin(4) * t546 + t384;
t321 = t335 * pkin(5) + t358;
t507 = -qJD(5) * t562 + t592 * t330 - t499 * t332;
t494 = pkin(4) * t592 + pkin(5);
t415 = t449 * t472;
t405 = t518 * t473;
t404 = t476 * t473;
t383 = qJ(6) * t518 + t557;
t382 = -qJ(6) * t476 + t529;
t340 = -t465 * t544 - t499 * t542 - t552 * t576 + (-t465 * t499 + t473 * t598) * t502;
t339 = t421 * t473 + t465 * t518;
t327 = -qJ(6) * t404 + t562;
t325 = pkin(5) * t472 - qJ(6) * t405 + t531;
t315 = t564 - t602;
t314 = t532 + t610;
t313 = t323 - t610;
t310 = -qJ(6) * t340 - qJD(6) * t404 + t513;
t309 = t462 * pkin(5) + t339 * qJ(6) - t405 * qJD(6) + t507;
t308 = -qJ(6) * t335 - qJD(6) * t369 + t527;
t307 = t449 * pkin(5) + t334 * qJ(6) - qJD(6) * t519 + t508;
t1 = [t603 * t605 + (-pkin(7) * t569 + t501 * t534) * MDP(9) + (pkin(7) * t570 + t503 * t534) * MDP(10) + (-t385 * t472 + t399 * t463 + t400 * t461 + t413 * t465 - t414 * t462 + t425 * t520 + t522) * MDP(11) + (t384 * t425 + t385 * t426 - t413 * t399 + t414 * t400 + (t479 + t523) * t548) * MDP(12) + (t386 * t575 + t428 * t516) * MDP(13) + (-(t427 * t502 - t428 * t500) * t465 + (-t502 * t546 - t585 + (-t500 * t427 - t428 * t502) * qJD(4)) * t473) * MDP(14) + (t386 * t472 + t428 * t462 + t431 * t473 + t455 * t516) * MDP(15) + (t427 * t462 - t455 * t608 - t472 * t546 - t473 * t571) * MDP(16) + (t455 * t462 + t415) * MDP(17) + ((-t426 * t553 + t395) * t455 + t407 * t449 + (-t409 * t553 + t380) * t472 + t360 * t462 - t399 * t427 + t425 * t546 + t408 * t541 + ((-qJD(4) * t412 - t400) * t455 + (-qJD(4) * t390 - t385) * t472 - t408 * t465 + t522) * t500) * MDP(18) + (-t361 * t462 + t384 * t575 + t425 * t386 + t399 * t428 + t408 * t516 - t449 * t558 - t455 * t514 - t472 * t515) * MDP(19) + (-t334 * t405 - t339 * t519) * MDP(20) + (t334 * t404 - t335 * t405 + t339 * t369 - t340 * t519) * MDP(21) + (-t334 * t472 - t339 * t448 + t405 * t449 + t462 * t519) * MDP(22) + (-t335 * t472 - t340 * t448 - t369 * t462 - t404 * t449) * MDP(23) + (t448 * t462 + t415) * MDP(24) + (t322 * t462 + t398 * t335 + t366 * t340 + t358 * t404 + t365 * t369 + t448 * t507 + t449 * t531 + t472 * t508) * MDP(25) + (-t323 * t462 - t398 * t334 - t366 * t339 + t358 * t405 + t365 * t519 - t448 * t513 - t449 * t562 - t472 * t527) * MDP(26) + (-t307 * t405 - t308 * t404 - t309 * t519 - t310 * t369 + t311 * t339 - t313 * t340 + t325 * t334 - t327 * t335) * MDP(27) + (t308 * t327 + t313 * t310 + t307 * t325 + t311 * t309 + t321 * (pkin(5) * t404 + t398) + t343 * (pkin(5) * t340 + t365)) * MDP(28) + MDP(6) * t569 - MDP(7) * t570 + 0.2e1 * t537 * t604; -t568 * t604 + t505 * t603 + ((t414 - t416) * t463 + (-t417 + t413) * t461 + (-t497 * t449 - t498 * t520) * pkin(2)) * MDP(11) + (t413 * t416 - t414 * t417 + (-t384 * t498 + t385 * t497 - t479 * t551) * pkin(2)) * MDP(12) + (t428 * t528 + t585) * MDP(13) + ((t386 + t584) * t502 + (-t546 - t582) * t500) * MDP(14) + (t455 * t528 + t571 - t581) * MDP(15) + (t521 - t583) * MDP(16) - t455 * t463 * MDP(17) + (t490 * t546 - t384 * t502 - t360 * t463 + t416 * t427 + (-t489 * t553 - t394) * t455 + (t417 * t455 + t511) * t500) * MDP(18) + (t361 * t463 + t384 * t500 + t490 * t386 - t416 * t428 + (t489 * t554 + t559) * t455 + t511 * t502) * MDP(19) + (-t334 * t476 - t519 * t560) * MDP(20) + (t525 - t595) * MDP(21) + (t580 - t586) * MDP(22) + (t524 + t587) * MDP(23) + (-t322 * t463 + t480 * t335 - t358 * t518 - t366 * t561 + t369 * t526 + t529 * t449) * MDP(25) + (t323 * t463 - t480 * t334 + t358 * t476 - t560 * t366 - t557 * t449 + t519 * t526) * MDP(26) + (-t307 * t476 + t308 * t518 + t311 * t560 + t313 * t561 + t334 * t382 - t335 * t383 - t369 * t566 - t519 * t565) * MDP(27) + (t308 * t383 + t307 * t382 + t321 * (-pkin(5) * t518 + t480) + (-pkin(5) * t561 + t526) * t343 + t566 * t313 + t565 * t311) * MDP(28) + (-t560 * MDP(22) - t463 * MDP(24) + MDP(25) * t600 + MDP(26) * t599) * t448 + (MDP(9) * t501 * t505 + MDP(10) * t568) * pkin(1); (-t461 ^ 2 - t463 ^ 2) * MDP(11) + (t413 * t463 - t414 * t461 + t487) * MDP(12) + (t521 + t583) * MDP(18) + (-t455 ^ 2 * t502 - t571 - t581) * MDP(19) + (t524 - t587) * MDP(25) + (t448 * t560 - t580 - t586) * MDP(26) + (t525 + t595) * MDP(27) + (t307 * t518 + t308 * t476 + t311 * t561 - t313 * t560 - t343 * t463) * MDP(28); -t428 * t427 * MDP(13) + (-t427 ^ 2 + t428 ^ 2) * MDP(14) + (t386 - t584) * MDP(15) + (-t546 + t582) * MDP(16) + t449 * MDP(17) + (t361 * t455 - t408 * t428 + t510) * MDP(18) + (t360 * t455 - t408 * t427 - t515) * MDP(19) + (-t532 * t448 + (-t428 * t369 - t448 * t552 + t449 * t592) * pkin(4) + t596) * MDP(25) + (t564 * t448 + (-t428 * t519 - t448 * t539 - t499 * t449) * pkin(4) + t607) * MDP(26) + (-t311 * t369 + t313 * t519 + t314 * t519 + t315 * t369 + t494 * t334 + (-t335 * t499 + (-t369 * t592 + t499 * t519) * qJD(5)) * pkin(4)) * MDP(27) + (-pkin(5) * t601 + t307 * t494 - t311 * t314 - t313 * t315 + (t308 * t499 - t343 * t428 + (-t311 * t499 + t313 * t592) * qJD(5)) * pkin(4)) * MDP(28) + t611; (t323 * t448 + t596) * MDP(25) + (t322 * t448 + t607) * MDP(26) + (pkin(5) * t334 - t369 * t567) * MDP(27) + (t567 * t313 + (t307 - t601) * pkin(5)) * MDP(28) + t611; (-t367 - t593) * MDP(27) + (t311 * t519 + t313 * t369 + t321) * MDP(28);];
tauc  = t1;
