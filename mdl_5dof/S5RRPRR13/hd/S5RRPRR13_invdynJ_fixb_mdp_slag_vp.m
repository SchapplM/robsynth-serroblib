% Calculate vector of inverse dynamics joint torques for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR13_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR13_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:34:07
% EndTime: 2019-12-31 20:34:19
% DurationCPUTime: 8.48s
% Computational Cost: add. (4264->496), mult. (9779->658), div. (0->0), fcn. (7360->14), ass. (0->204)
t505 = cos(qJ(2));
t568 = qJD(1) * t505;
t478 = -qJD(4) + t568;
t473 = -qJD(5) + t478;
t497 = sin(pkin(9));
t501 = sin(qJ(2));
t569 = qJD(1) * t501;
t550 = t497 * t569;
t498 = cos(pkin(9));
t558 = t498 * qJD(2);
t448 = t550 - t558;
t549 = t498 * t569;
t567 = qJD(2) * t497;
t450 = t549 + t567;
t500 = sin(qJ(4));
t504 = cos(qJ(4));
t384 = t504 * t448 + t450 * t500;
t503 = cos(qJ(5));
t383 = t448 * t500 - t450 * t504;
t499 = sin(qJ(5));
t586 = t383 * t499;
t606 = -t503 * t384 + t586;
t605 = t473 * t606;
t532 = pkin(2) * t501 - qJ(3) * t505;
t459 = t532 * qJD(1);
t439 = t497 * t459;
t582 = t498 * t501;
t583 = t497 * t505;
t525 = -pkin(6) * t582 - pkin(7) * t583;
t390 = qJD(1) * t525 + t439;
t607 = qJD(3) * t498 - t390;
t555 = qJDD(1) * t501;
t483 = pkin(6) * t555;
t556 = qJD(1) * qJD(2);
t546 = t505 * t556;
t431 = -qJDD(2) * pkin(2) + pkin(6) * t546 + qJDD(3) + t483;
t502 = sin(qJ(1));
t506 = cos(qJ(1));
t536 = g(1) * t506 + g(2) * t502;
t590 = g(3) * t505;
t516 = t501 * t536 - t590;
t597 = t431 - t516;
t604 = t383 * t478;
t603 = t384 * t478;
t531 = t383 * t503 + t384 * t499;
t602 = t473 * t531;
t585 = t497 * t500;
t456 = -t504 * t498 + t585;
t522 = t456 * t505;
t574 = qJD(1) * t522 - t456 * qJD(4);
t457 = t497 * t504 + t498 * t500;
t523 = t457 * t505;
t573 = -qJD(1) * t523 + t457 * qJD(4);
t527 = pkin(2) * t505 + qJ(3) * t501 + pkin(1);
t441 = t527 * qJD(1);
t485 = pkin(6) * t568;
t468 = qJD(2) * qJ(3) + t485;
t392 = -t498 * t441 - t468 * t497;
t553 = pkin(3) * t568;
t352 = -pkin(7) * t450 + t392 - t553;
t393 = -t497 * t441 + t498 * t468;
t355 = -pkin(7) * t448 + t393;
t323 = t352 * t500 + t355 * t504;
t317 = -pkin(8) * t384 + t323;
t560 = qJD(5) * t499;
t315 = t317 * t560;
t461 = -qJD(2) * pkin(2) + pkin(6) * t569 + qJD(3);
t402 = pkin(3) * t448 + t461;
t348 = pkin(4) * t384 + t402;
t494 = pkin(9) + qJ(4);
t492 = qJ(5) + t494;
t480 = sin(t492);
t481 = cos(t492);
t580 = t502 * t505;
t410 = t480 * t506 - t481 * t580;
t579 = t505 * t506;
t412 = t480 * t502 + t481 * t579;
t591 = g(3) * t501;
t600 = g(1) * t412 - g(2) * t410 - t348 * t606 + t481 * t591 + t315;
t409 = t480 * t580 + t481 * t506;
t411 = -t480 * t579 + t481 * t502;
t487 = t498 * qJDD(2);
t521 = t546 + t555;
t407 = t497 * t521 - t487;
t554 = qJDD(2) * t497;
t408 = t498 * t521 + t554;
t561 = qJD(4) * t504;
t563 = qJD(4) * t500;
t334 = -t500 * t407 + t504 * t408 - t448 * t561 - t450 * t563;
t491 = t505 * qJDD(1);
t520 = t501 * t556 - t491;
t455 = qJDD(4) + t520;
t433 = qJD(2) * t532 - qJD(3) * t501;
t382 = qJD(1) * t433 - qJDD(1) * t527;
t421 = -pkin(6) * t520 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t344 = t498 * t382 - t421 * t497;
t327 = pkin(3) * t520 - pkin(7) * t408 + t344;
t345 = t497 * t382 + t498 * t421;
t336 = -pkin(7) * t407 + t345;
t543 = t504 * t327 - t336 * t500;
t509 = -t323 * qJD(4) + t543;
t306 = pkin(4) * t455 - pkin(8) * t334 + t509;
t335 = -qJD(4) * t383 + t504 * t407 + t408 * t500;
t519 = t500 * t327 + t504 * t336 + t352 * t561 - t355 * t563;
t307 = -pkin(8) * t335 + t519;
t544 = t503 * t306 - t499 * t307;
t599 = -g(1) * t411 + g(2) * t409 + t348 * t531 + t480 * t591 + t544;
t442 = qJDD(5) + t455;
t598 = t442 * MDP(26) + (t531 ^ 2 - t606 ^ 2) * MDP(23) + t606 * MDP(22) * t531;
t447 = t498 * t527;
t391 = -pkin(7) * t582 - t447 + (-pkin(6) * t497 - pkin(3)) * t505;
t581 = t498 * t505;
t414 = pkin(6) * t581 - t497 * t527;
t584 = t497 * t501;
t398 = -pkin(7) * t584 + t414;
t575 = t500 * t391 + t504 * t398;
t589 = pkin(7) + qJ(3);
t466 = t589 * t497;
t467 = t589 * t498;
t572 = -t500 * t466 + t504 * t467;
t403 = pkin(6) * t550 + t498 * t459;
t528 = pkin(3) * t501 - pkin(7) * t581;
t372 = qJD(1) * t528 + t403;
t529 = qJD(3) * t497 + qJD(4) * t467;
t596 = -t466 * t561 + t607 * t504 + (-t372 - t529) * t500;
t535 = g(1) * t502 - g(2) * t506;
t542 = t499 * t334 + t503 * t335;
t311 = -qJD(5) * t531 + t542;
t595 = pkin(6) * t448;
t594 = pkin(6) * t450;
t322 = t504 * t352 - t355 * t500;
t316 = pkin(8) * t383 + t322;
t314 = -pkin(4) * t478 + t316;
t588 = t314 * t503;
t587 = t317 * t503;
t388 = t503 * t456 + t457 * t499;
t578 = -qJD(5) * t388 - t499 * t573 + t503 * t574;
t389 = -t456 * t499 + t457 * t503;
t577 = qJD(5) * t389 + t499 * t574 + t503 * t573;
t566 = qJD(2) * t501;
t552 = pkin(6) * t566;
t396 = t498 * t433 + t497 * t552;
t443 = t497 * t553 + t485;
t565 = qJD(2) * t505;
t444 = (pkin(3) * t497 + pkin(6)) * t565;
t460 = pkin(3) * t584 + t501 * pkin(6);
t495 = t501 ^ 2;
t571 = -t505 ^ 2 + t495;
t562 = qJD(4) * t501;
t559 = qJD(5) * t503;
t557 = qJD(3) - t461;
t551 = t503 * t334 - t499 * t335 - t384 * t559;
t482 = -pkin(3) * t498 - pkin(2);
t548 = qJ(3) * t491;
t545 = pkin(4) * t573 - t443;
t361 = qJD(2) * t528 + t396;
t419 = t497 * t433;
t374 = qJD(2) * t525 + t419;
t541 = t504 * t361 - t374 * t500;
t540 = t504 * t391 - t398 * t500;
t538 = -t504 * t466 - t467 * t500;
t537 = qJD(5) * t314 + t307;
t364 = t504 * t372;
t366 = -pkin(8) * t456 + t572;
t534 = pkin(4) * t569 + pkin(8) * t574 + t457 * qJD(3) + t572 * qJD(4) + qJD(5) * t366 - t390 * t500 + t364;
t365 = -pkin(8) * t457 + t538;
t533 = -pkin(8) * t573 + qJD(5) * t365 + t596;
t309 = t314 * t499 + t587;
t429 = t457 * t501;
t430 = t456 * t501;
t359 = t503 * t429 - t430 * t499;
t360 = -t429 * t499 - t430 * t503;
t526 = -0.2e1 * pkin(1) * t556 - pkin(6) * qJDD(2);
t518 = t500 * t361 + t504 * t374 + t391 * t561 - t398 * t563;
t310 = t383 * t560 + t551;
t508 = qJD(1) ^ 2;
t517 = pkin(1) * t508 + t536;
t362 = pkin(3) * t407 + t431;
t514 = -t505 * t536 - t591;
t507 = qJD(2) ^ 2;
t512 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t507 + t535;
t490 = cos(t494);
t489 = sin(t494);
t425 = t489 * t502 + t490 * t579;
t424 = -t489 * t579 + t490 * t502;
t423 = t489 * t506 - t490 * t580;
t422 = t489 * t580 + t490 * t506;
t418 = pkin(4) * t456 + t482;
t413 = -pkin(6) * t583 - t447;
t404 = -pkin(6) * t549 + t439;
t397 = -t498 * t552 + t419;
t395 = pkin(4) * t429 + t460;
t368 = qJD(2) * t523 + t561 * t582 - t562 * t585;
t367 = -qJD(2) * t522 - t457 * t562;
t349 = pkin(4) * t368 + t444;
t333 = -pkin(8) * t429 + t575;
t328 = -pkin(4) * t505 + pkin(8) * t430 + t540;
t320 = pkin(4) * t335 + t362;
t319 = qJD(5) * t360 + t367 * t499 + t503 * t368;
t318 = -qJD(5) * t359 + t367 * t503 - t368 * t499;
t313 = -pkin(8) * t368 + t518;
t312 = pkin(4) * t566 - pkin(8) * t367 - qJD(4) * t575 + t541;
t308 = -t317 * t499 + t588;
t1 = [0.2e1 * (t491 * t501 - t556 * t571) * MDP(5) + (-t396 * t450 - t397 * t448 - t407 * t414 - t408 * t413 + (-t392 * t498 - t393 * t497) * t565 + (-t344 * t498 - t345 * t497 + t535) * t501) * MDP(13) + (-t455 * t505 - t478 * t566) * MDP(19) + (-t442 * t505 - t473 * t566) * MDP(26) + (t335 * t505 + t368 * t478 - t384 * t566 - t429 * t455) * MDP(18) + (qJDD(1) * t495 + 0.2e1 * t501 * t546) * MDP(4) + (-t309 * t566 - g(1) * t409 - g(2) * t411 + t395 * t310 - t315 * t505 + t348 * t318 + t320 * t360 - t349 * t531 + ((-qJD(5) * t333 + t312) * t473 - t328 * t442 + t306 * t505) * t499 + ((qJD(5) * t328 + t313) * t473 - t333 * t442 + t537 * t505) * t503) * MDP(28) + (-t310 * t505 - t318 * t473 + t360 * t442 - t531 * t566) * MDP(24) + (t310 * t360 - t318 * t531) * MDP(22) + (-g(1) * t422 - g(2) * t424 - t323 * t566 + t460 * t334 - t362 * t430 + t402 * t367 - t383 * t444 - t455 * t575 + t478 * t518 + t505 * t519) * MDP(21) + (-t334 * t505 - t367 * t478 - t383 * t566 - t430 * t455) * MDP(17) + (-t334 * t430 - t367 * t383) * MDP(15) + (-t334 * t429 + t335 * t430 - t367 * t384 + t368 * t383) * MDP(16) + (-t536 * t498 + (pkin(6) * t408 + t431 * t498 + (-qJD(1) * t414 - t393) * qJD(2)) * t501 + (t397 * qJD(1) + t414 * qJDD(1) + t345 - t535 * t497 + (t461 * t498 + t594) * qJD(2)) * t505) * MDP(12) + (-t536 * t497 + (pkin(6) * t407 + t431 * t497 + (qJD(1) * t413 + t392) * qJD(2)) * t501 + (-t396 * qJD(1) - t413 * qJDD(1) - t344 + t535 * t498 + (t461 * t497 + t595) * qJD(2)) * t505) * MDP(11) + t536 * MDP(3) + t535 * MDP(2) + (t344 * t413 + t345 * t414 + t392 * t396 + t393 * t397 + (t431 * t501 + t461 * t565 - t536) * pkin(6) + t535 * t527) * MDP(14) + (t501 * t526 + t505 * t512) * MDP(9) + (-t501 * t512 + t505 * t526) * MDP(10) + (t311 * t505 + t319 * t473 - t359 * t442 + t566 * t606) * MDP(25) + (-(t312 * t503 - t313 * t499) * t473 + (t328 * t503 - t333 * t499) * t442 - t544 * t505 + t308 * t566 - t349 * t606 + t395 * t311 + t320 * t359 + t348 * t319 - g(1) * t410 - g(2) * t412 + (-(-t328 * t499 - t333 * t503) * t473 + t309 * t505) * qJD(5)) * MDP(27) + (-t310 * t359 - t311 * t360 + t318 * t606 + t319 * t531) * MDP(23) + (-t541 * t478 + t540 * t455 - t543 * t505 + t322 * t566 + t444 * t384 + t460 * t335 + t362 * t429 + t402 * t368 - g(1) * t423 - g(2) * t425 + (t323 * t505 + t478 * t575) * qJD(4)) * MDP(20) + qJDD(1) * MDP(1) + (qJDD(2) * t501 + t505 * t507) * MDP(6) + (qJDD(2) * t505 - t501 * t507) * MDP(7); MDP(7) * t491 + MDP(6) * t555 + ((t365 * t503 - t366 * t499) * t442 + t418 * t311 + t320 * t388 + (t499 * t533 + t503 * t534) * t473 + t577 * t348 - t545 * t606 + t516 * t481) * MDP(27) + (-t388 * t442 + t473 * t577) * MDP(25) + (-(t365 * t499 + t366 * t503) * t442 + t418 * t310 + t320 * t389 + (-t499 * t534 + t503 * t533) * t473 + t578 * t348 - t545 * t531 - t516 * t480) * MDP(28) + (t310 * t389 - t531 * t578) * MDP(22) + (-t310 * t388 - t311 * t389 + t531 * t577 + t578 * t606) * MDP(23) + (t389 * t442 - t473 * t578) * MDP(24) + (t538 * t455 + t482 * t335 + t362 * t456 - t443 * t384 + (t364 + t529 * t504 + (-qJD(4) * t466 + t607) * t500) * t478 + t573 * t402 + t516 * t490) * MDP(20) + (-t455 * t456 + t478 * t573) * MDP(18) + (t334 * t457 - t383 * t574) * MDP(15) + (-t334 * t456 - t335 * t457 + t383 * t573 - t384 * t574) * MDP(16) + (t455 * t457 - t478 * t574) * MDP(17) + (t482 * t334 + t362 * t457 + t383 * t443 + t574 * t402 - t572 * t455 + t596 * t478 - t489 * t516) * MDP(21) + (t403 * t450 + t404 * t448 + (-qJ(3) * t407 - qJD(3) * t448 + t392 * t568 + t345) * t498 + (qJ(3) * t408 + qJD(3) * t450 + t393 * t568 - t344) * t497 + t514) * MDP(13) + (-t461 * t485 - t392 * t403 - t393 * t404 + (-t392 * t497 + t393 * t498) * qJD(3) - t597 * pkin(2) + (-t344 * t497 + t345 * t498 + t514) * qJ(3)) * MDP(14) + (t591 + (-pkin(6) * qJDD(1) + t517) * t505) * MDP(10) + (t498 * t548 - pkin(2) * t408 + t597 * t497 + ((-qJ(3) * t558 + t393) * t501 + (t498 * t557 - t404 - t594) * t505) * qJD(1)) * MDP(12) + (t497 * t548 - pkin(2) * t407 - t597 * t498 + ((-qJ(3) * t567 - t392) * t501 + (t497 * t557 + t403 - t595) * t505) * qJD(1)) * MDP(11) + (t501 * t517 - t483 - t590) * MDP(9) + qJDD(2) * MDP(8) + (MDP(17) * t383 + t384 * MDP(18) + t478 * MDP(19) - t322 * MDP(20) + t323 * MDP(21) + MDP(24) * t531 - MDP(25) * t606 + t473 * MDP(26) - t308 * MDP(27) + t309 * MDP(28)) * t569 + (-MDP(4) * t501 * t505 + MDP(5) * t571) * t508; (t497 * t555 - t487 + (-t450 + t567) * t568) * MDP(11) + (t498 * t555 + t554 + (t448 + t558) * t568) * MDP(12) + (-t448 ^ 2 - t450 ^ 2) * MDP(13) + (t392 * t450 + t393 * t448 + t597) * MDP(14) + (t335 + t604) * MDP(20) + (t334 + t603) * MDP(21) + (t311 + t602) * MDP(27) + (t310 - t605) * MDP(28); -t383 * t384 * MDP(15) + (t383 ^ 2 - t384 ^ 2) * MDP(16) + (t334 - t603) * MDP(17) + (-t335 + t604) * MDP(18) + t455 * MDP(19) + (-g(1) * t424 + g(2) * t422 - t323 * t478 + t383 * t402 + t489 * t591 + t509) * MDP(20) + (g(1) * t425 - g(2) * t423 - t322 * t478 + t384 * t402 + t490 * t591 - t519) * MDP(21) + (t310 + t605) * MDP(24) + (-t311 + t602) * MDP(25) + ((-t316 * t499 - t587) * t473 - t309 * qJD(5) + (-t383 * t606 + t442 * t503 + t473 * t560) * pkin(4) + t599) * MDP(27) + ((t317 * t473 - t306) * t499 + (-t316 * t473 - t537) * t503 + (-t383 * t531 - t442 * t499 + t473 * t559) * pkin(4) + t600) * MDP(28) + t598; (t551 + t605) * MDP(24) + (-t542 + t602) * MDP(25) + (-t309 * t473 + t599) * MDP(27) + (-t499 * t306 - t503 * t307 - t308 * t473 + t600) * MDP(28) + (MDP(24) * t586 + MDP(25) * t531 - MDP(27) * t309 - MDP(28) * t588) * qJD(5) + t598;];
tau = t1;
