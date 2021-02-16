% Calculate Coriolis joint torque vector for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:53
% EndTime: 2021-01-16 03:47:11
% DurationCPUTime: 6.73s
% Computational Cost: add. (4134->454), mult. (10718->632), div. (0->0), fcn. (8480->12), ass. (0->212)
t454 = cos(qJ(3));
t558 = cos(pkin(12));
t495 = t558 * t454;
t433 = qJD(2) * t495;
t445 = sin(pkin(12));
t450 = sin(qJ(3));
t521 = qJD(2) * t450;
t412 = t445 * t521 - t433;
t564 = qJD(5) + qJD(6);
t577 = t412 + t564;
t451 = sin(qJ(2));
t446 = sin(pkin(6));
t523 = qJD(1) * t446;
t508 = t451 * t523;
t519 = qJD(3) * t450;
t576 = pkin(3) * t519 - t508;
t560 = qJ(4) + pkin(8);
t499 = qJD(3) * t560;
t410 = qJD(4) * t454 - t450 * t499;
t464 = -qJD(4) * t450 - t454 * t499;
t358 = t558 * t410 + t445 * t464;
t468 = -t445 * t450 + t495;
t455 = cos(qJ(2));
t507 = t455 * t523;
t390 = t468 * t507;
t529 = t358 - t390;
t422 = t445 * t454 + t558 * t450;
t414 = t422 * qJD(3);
t417 = t468 * qJD(3);
t575 = pkin(4) * t414 - pkin(9) * t417 + t576;
t426 = qJD(2) * pkin(8) + t508;
t490 = qJ(4) * qJD(2) + t426;
t447 = cos(pkin(6));
t522 = qJD(1) * t447;
t506 = t450 * t522;
t386 = t490 * t454 + t506;
t374 = t445 * t386;
t434 = t454 * t522;
t385 = -t490 * t450 + t434;
t379 = qJD(3) * pkin(3) + t385;
t322 = t558 * t379 - t374;
t319 = -qJD(3) * pkin(4) - t322;
t415 = t422 * qJD(2);
t449 = sin(qJ(5));
t453 = cos(qJ(5));
t514 = t453 * qJD(3);
t393 = t415 * t449 - t514;
t315 = t393 * pkin(5) + t319;
t452 = cos(qJ(6));
t395 = qJD(3) * t449 + t415 * t453;
t448 = sin(qJ(6));
t548 = t395 * t448;
t333 = t452 * t393 + t548;
t574 = t315 * t333;
t407 = qJD(5) + t412;
t403 = qJD(6) + t407;
t573 = t333 * t403;
t476 = t393 * t448 - t452 * t395;
t572 = t403 * t476;
t425 = t448 * t453 + t449 * t452;
t526 = t577 * t425;
t518 = qJD(5) * t449;
t545 = t412 * t449;
t571 = t518 + t545;
t513 = qJD(2) * qJD(3);
t500 = t450 * t513;
t431 = t445 * t500;
t405 = qJD(3) * t433 - t431;
t346 = qJD(5) * t514 + t453 * t405 - t415 * t518;
t404 = qJD(2) * t414;
t496 = t558 * t386;
t323 = t445 * t379 + t496;
t320 = qJD(3) * pkin(9) + t323;
t441 = -pkin(3) * t454 - pkin(2);
t406 = t441 * qJD(2) + qJD(4) - t507;
t339 = pkin(4) * t412 - pkin(9) * t415 + t406;
t308 = t320 * t453 + t339 * t449;
t481 = qJD(4) + t507;
t345 = (-t426 * t450 + t434) * qJD(3) + (-qJ(4) * t519 + t481 * t454) * qJD(2);
t497 = t558 * t345;
t472 = -t426 * t454 - t506;
t503 = qJD(3) * t454 * qJ(4);
t565 = t472 * qJD(3) + (-t481 * t450 - t503) * qJD(2);
t313 = t565 * t445 + t497;
t520 = qJD(2) * t451;
t505 = t446 * t520;
t411 = pkin(3) * t500 + qJD(1) * t505;
t331 = pkin(4) * t404 - pkin(9) * t405 + t411;
t329 = t453 * t331;
t460 = -t308 * qJD(5) - t313 * t449 + t329;
t290 = pkin(5) * t404 - pkin(10) * t346 + t460;
t347 = qJD(5) * t395 + t405 * t449;
t517 = qJD(5) * t453;
t467 = t453 * t313 - t320 * t518 + t449 * t331 + t339 * t517;
t291 = -pkin(10) * t347 + t467;
t494 = t452 * t290 - t448 * t291;
t570 = t315 * t476 + t494;
t569 = t404 * MDP(27) + (-t333 ^ 2 + t476 ^ 2) * MDP(24) - t333 * MDP(23) * t476;
t366 = t425 * t422;
t568 = (t450 ^ 2 - t454 ^ 2) * MDP(6);
t567 = t390 * t449 + t575 * t453;
t530 = t410 * t445 - t422 * t507 - t558 * t464;
t372 = -pkin(4) * t468 - pkin(9) * t422 + t441;
t429 = t560 * t450;
t430 = t560 * t454;
t392 = -t445 * t429 + t558 * t430;
t566 = -t372 * t517 + t392 * t518 - t575 * t449 - t529 * t453;
t536 = t453 * t417;
t470 = -t422 * t518 + t536;
t424 = t448 * t449 - t452 * t453;
t527 = t577 * t424;
t563 = t527 * t403 - t404 * t425;
t492 = t346 * t448 + t452 * t347;
t298 = -t476 * qJD(6) + t492;
t562 = pkin(3) * t450;
t437 = pkin(3) * t445 + pkin(9);
t561 = pkin(10) + t437;
t559 = qJD(2) * pkin(2);
t307 = -t320 * t449 + t453 * t339;
t301 = -pkin(10) * t395 + t307;
t296 = pkin(5) * t407 + t301;
t557 = t296 * t452;
t302 = -pkin(10) * t393 + t308;
t556 = t302 * t452;
t555 = t333 * t415;
t554 = t476 * t415;
t553 = t346 * t449;
t552 = t393 * t407;
t551 = t393 * t415;
t550 = t395 * t407;
t549 = t395 * t415;
t544 = t422 * t449;
t543 = t422 * t453;
t542 = t446 * t451;
t541 = t446 * t455;
t457 = qJD(2) ^ 2;
t540 = t446 * t457;
t539 = t449 * t404;
t538 = t449 * t417;
t456 = qJD(3) ^ 2;
t537 = t450 * t456;
t380 = t453 * t392;
t396 = t453 * t404;
t535 = t454 * t456;
t534 = -pkin(10) * t536 + pkin(5) * t414 - t358 * t449 + (-t380 + (pkin(10) * t422 - t372) * t449) * qJD(5) + t567;
t501 = t422 * t517;
t471 = t501 + t538;
t533 = pkin(10) * t471 + t566;
t532 = pkin(5) * t471 + t530;
t327 = t558 * t385 - t374;
t512 = pkin(3) * t521;
t359 = pkin(4) * t415 + pkin(9) * t412 + t512;
t531 = t453 * t327 + t449 * t359;
t528 = t449 * t372 + t380;
t516 = qJD(6) * t448;
t515 = qJD(6) * t452;
t510 = t451 * t540;
t509 = t452 * t346 - t448 * t347 - t393 * t515;
t504 = qJD(2) * t541;
t498 = qJD(5) * t561;
t299 = t302 * t516;
t493 = t448 * t290 - t299;
t312 = t345 * t445 - t558 * t565;
t325 = t385 * t445 + t496;
t391 = t558 * t429 + t430 * t445;
t491 = t407 * t453;
t489 = qJD(6) * t296 + t291;
t488 = t454 * t504;
t487 = t450 * t504;
t438 = -t558 * pkin(3) - pkin(4);
t486 = t571 * pkin(5) - t325;
t485 = -t526 * t403 - t424 * t404;
t352 = t453 * t359;
t420 = t561 * t453;
t483 = pkin(5) * t415 + qJD(6) * t420 - t327 * t449 + t352 + (pkin(10) * t412 + t498) * t453;
t419 = t561 * t449;
t482 = pkin(10) * t545 + qJD(6) * t419 + t449 * t498 + t531;
t293 = t296 * t448 + t556;
t480 = t312 * t422 - t392 * t404;
t369 = t453 * t372;
t314 = -pkin(5) * t468 - pkin(10) * t543 - t392 * t449 + t369;
t316 = -pkin(10) * t544 + t528;
t479 = t314 * t448 + t316 * t452;
t418 = t447 * t450 + t454 * t542;
t473 = t447 * t454 - t450 * t542;
t365 = t558 * t418 + t445 * t473;
t343 = -t365 * t449 - t453 * t541;
t474 = -t365 * t453 + t449 * t541;
t478 = t343 * t452 + t448 * t474;
t477 = t343 * t448 - t452 * t474;
t475 = -t571 * t407 + t396;
t469 = qJD(2) * t559;
t297 = -t395 * t516 + t509;
t465 = t407 * t319 - t437 * t404;
t463 = t418 * qJD(3);
t461 = -0.2e1 * qJD(3) * t559;
t458 = -t463 - t487;
t428 = -t453 * pkin(5) + t438;
t384 = qJD(3) * t473 + t488;
t373 = t404 * t468;
t367 = t424 * t422;
t364 = t418 * t445 - t558 * t473;
t356 = pkin(5) * t544 + t391;
t326 = t558 * t384 + t445 * t458;
t324 = t384 * t445 - t558 * t458;
t310 = -t516 * t544 + (t564 * t543 + t538) * t452 + t470 * t448;
t309 = -t564 * t366 - t424 * t417;
t305 = qJD(5) * t474 - t326 * t449 + t453 * t505;
t304 = qJD(5) * t343 + t326 * t453 + t449 * t505;
t300 = pkin(5) * t347 + t312;
t292 = -t302 * t448 + t557;
t1 = [-MDP(3) * t510 - t455 * MDP(4) * t540 + (-t454 * t510 + (-t463 - 0.2e1 * t487) * qJD(3)) * MDP(10) + (t450 * t510 + (-t384 - t488) * qJD(3)) * MDP(11) + (-qJD(3) * t324 + (-t404 * t455 + t412 * t520) * t446) * MDP(12) + (-qJD(3) * t326 + (-t405 * t455 + t415 * t520) * t446) * MDP(13) + (t324 * t415 - t326 * t412 + t364 * t405 - t365 * t404) * MDP(14) + (t312 * t364 + t313 * t365 - t322 * t324 + t323 * t326 + (t406 * t520 - t411 * t455) * t446) * MDP(15) + (t305 * t407 + t324 * t393 + t343 * t404 + t347 * t364) * MDP(21) + (-t304 * t407 + t324 * t395 + t346 * t364 + t404 * t474) * MDP(22) + ((-t477 * qJD(6) - t304 * t448 + t305 * t452) * t403 + t478 * t404 + t324 * t333 + t364 * t298) * MDP(28) + (-(t478 * qJD(6) + t304 * t452 + t305 * t448) * t403 - t477 * t404 - t324 * t476 + t364 * t297) * MDP(29); -0.2e1 * t513 * t568 + (-pkin(8) * t535 + t450 * t461) * MDP(10) + (pkin(8) * t537 + t454 * t461) * MDP(11) + (-t412 * t508 + t404 * t441 + t406 * t414 - t411 * t468 + (t412 * t562 - t530) * qJD(3)) * MDP(12) + (-t415 * t508 + t405 * t441 + t406 * t417 + t411 * t422 + (t415 * t562 - t529) * qJD(3)) * MDP(13) + (t313 * t468 - t322 * t417 - t323 * t414 + t391 * t405 - t529 * t412 + t530 * t415 + t480) * MDP(14) + (t312 * t391 + t313 * t392 - t530 * t322 + t529 * t323 + t576 * t406 + t411 * t441) * MDP(15) + (t346 * t543 + t395 * t470) * MDP(16) + ((-t393 * t453 - t395 * t449) * t417 + (-t553 - t347 * t453 + (t393 * t449 - t395 * t453) * qJD(5)) * t422) * MDP(17) + (-t346 * t468 + t395 * t414 + t422 * t396 + t407 * t470) * MDP(18) + (t347 * t468 - t393 * t414 - t407 * t471 - t422 * t539) * MDP(19) + (t407 * t414 - t373) * MDP(20) + (t369 * t404 - (-t320 * t517 + t329) * t468 + t307 * t414 + t391 * t347 + t319 * t501 + (-t392 * t517 + t567) * t407 + t530 * t393 + ((-qJD(5) * t372 - t358) * t407 - (-qJD(5) * t339 - t313) * t468 + t319 * t417 + t480) * t449) * MDP(21) + (-t308 * t414 + t312 * t543 + t470 * t319 + t391 * t346 + t530 * t395 - t528 * t404 + t566 * t407 + t467 * t468) * MDP(22) + (-t297 * t367 - t309 * t476) * MDP(23) + (-t297 * t366 + t298 * t367 - t309 * t333 + t310 * t476) * MDP(24) + (-t297 * t468 + t309 * t403 - t367 * t404 - t414 * t476) * MDP(25) + (t298 * t468 - t310 * t403 - t333 * t414 - t366 * t404) * MDP(26) + (t403 * t414 - t373) * MDP(27) + ((t314 * t452 - t316 * t448) * t404 - t494 * t468 + t292 * t414 + t356 * t298 + t300 * t366 + t315 * t310 + (t533 * t448 + t534 * t452) * t403 + t532 * t333 + (t293 * t468 - t479 * t403) * qJD(6)) * MDP(28) + (-t479 * t404 + (t489 * t452 + t493) * t468 - t293 * t414 + t356 * t297 - t300 * t367 + t315 * t309 + ((-qJD(6) * t314 + t533) * t452 + (qJD(6) * t316 - t534) * t448) * t403 - t532 * t476) * MDP(29) + 0.2e1 * t454 * MDP(5) * t500 + MDP(7) * t535 - MDP(8) * t537; t457 * t568 + t450 * MDP(10) * t469 + (qJD(3) * t325 - t406 * t415 - t412 * t512 - t312) * MDP(12) + (-t497 + t406 * t412 + (-t445 * t472 + t327) * qJD(3) + (t445 * t503 + (-pkin(3) * t415 + t445 * t481) * t450) * qJD(2)) * MDP(13) + ((t323 - t325) * t415 + (-t322 + t327) * t412 + (-t404 * t445 - t558 * t405) * pkin(3)) * MDP(14) + (t322 * t325 - t323 * t327 + (-t558 * t312 + t313 * t445 - t406 * t521) * pkin(3)) * MDP(15) + (t395 * t491 + t553) * MDP(16) + ((t346 - t552) * t453 + (-t347 - t550) * t449) * MDP(17) + (t407 * t491 + t539 - t549) * MDP(18) + (t475 + t551) * MDP(19) - t407 * t415 * MDP(20) + (-t307 * t415 - t312 * t453 - t325 * t393 + t438 * t347 + (-t437 * t517 - t352) * t407 + (t327 * t407 + t465) * t449) * MDP(21) + (t308 * t415 + t312 * t449 - t325 * t395 + t438 * t346 + (t437 * t518 + t531) * t407 + t465 * t453) * MDP(22) + (t297 * t425 + t476 * t527) * MDP(23) + (-t297 * t424 - t298 * t425 + t527 * t333 + t476 * t526) * MDP(24) + (t554 - t563) * MDP(25) + (t485 + t555) * MDP(26) - t403 * t415 * MDP(27) + ((-t419 * t452 - t420 * t448) * t404 + t428 * t298 + t300 * t424 - t292 * t415 + (t482 * t448 - t483 * t452) * t403 + t486 * t333 + t526 * t315) * MDP(28) + (-(-t419 * t448 + t420 * t452) * t404 + t428 * t297 + t300 * t425 + t293 * t415 + (t483 * t448 + t482 * t452) * t403 - t486 * t476 - t527 * t315) * MDP(29) + (-t450 * t457 * MDP(5) + t469 * MDP(11)) * t454; 0.2e1 * t415 * qJD(3) * MDP(12) + (-t431 + (t433 - t412) * qJD(3)) * MDP(13) + (-t412 ^ 2 - t415 ^ 2) * MDP(14) + (t322 * t415 + t323 * t412 + t411) * MDP(15) + (t475 - t551) * MDP(21) + (-t407 ^ 2 * t453 - t539 - t549) * MDP(22) + (t485 - t555) * MDP(28) + (t554 + t563) * MDP(29); t395 * t393 * MDP(16) + (-t393 ^ 2 + t395 ^ 2) * MDP(17) + (t346 + t552) * MDP(18) + (-t347 + t550) * MDP(19) + t404 * MDP(20) + (t308 * t407 - t319 * t395 + t460) * MDP(21) + (t307 * t407 + t319 * t393 - t467) * MDP(22) + (t297 + t573) * MDP(25) + (-t298 - t572) * MDP(26) + (-(-t301 * t448 - t556) * t403 - t293 * qJD(6) + (-t333 * t395 - t403 * t516 + t452 * t404) * pkin(5) + t570) * MDP(28) + (t574 + t299 + (-t302 * t403 - t290) * t448 + (t301 * t403 - t489) * t452 + (t395 * t476 - t403 * t515 - t448 * t404) * pkin(5)) * MDP(29) + t569; (t509 + t573) * MDP(25) + (-t492 - t572) * MDP(26) + (t293 * t403 + t570) * MDP(28) + (-t452 * t291 + t292 * t403 - t493 + t574) * MDP(29) + (-MDP(25) * t548 + t476 * MDP(26) - t293 * MDP(28) - MDP(29) * t557) * qJD(6) + t569;];
tauc = t1;
