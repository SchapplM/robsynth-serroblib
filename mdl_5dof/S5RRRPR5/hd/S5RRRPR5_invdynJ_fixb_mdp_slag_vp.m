% Calculate vector of inverse dynamics joint torques for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:58
% EndTime: 2021-01-15 23:11:14
% DurationCPUTime: 5.75s
% Computational Cost: add. (5019->427), mult. (12044->555), div. (0->0), fcn. (8877->14), ass. (0->196)
t463 = qJ(2) + qJ(3);
t454 = pkin(9) + t463;
t442 = sin(t454);
t468 = sin(qJ(1));
t472 = cos(qJ(1));
t506 = g(1) * t472 + g(2) * t468;
t575 = t506 * t442;
t455 = sin(t463);
t456 = cos(t463);
t574 = -g(3) * t456 + t506 * t455;
t460 = qJD(2) + qJD(3);
t466 = sin(qJ(3));
t470 = cos(qJ(3));
t471 = cos(qJ(2));
t532 = qJD(1) * t471;
t520 = t470 * t532;
t467 = sin(qJ(2));
t533 = qJD(1) * t467;
t521 = t466 * t533;
t525 = qJDD(1) * t471;
t527 = qJD(1) * qJD(2);
t518 = t471 * t527;
t526 = qJDD(1) * t467;
t570 = t518 + t526;
t350 = qJD(3) * t520 - t460 * t521 + t466 * t525 + t570 * t470;
t395 = -t466 * t532 - t470 * t533;
t458 = qJDD(2) + qJDD(3);
t566 = pkin(6) + pkin(7);
t371 = qJDD(2) * pkin(2) - t566 * t570;
t519 = t467 * t527;
t373 = t566 * (-t519 + t525);
t421 = t566 * t467;
t409 = qJD(1) * t421;
t558 = qJD(2) * pkin(2);
t400 = -t409 + t558;
t422 = t566 * t471;
t411 = qJD(1) * t422;
t547 = t411 * t470;
t501 = -t400 * t466 - t547;
t480 = t501 * qJD(3) + t470 * t371 - t466 * t373;
t308 = pkin(3) * t458 - qJ(4) * t350 + qJD(4) * t395 + t480;
t405 = t466 * t471 + t467 * t470;
t370 = t460 * t405;
t504 = t466 * t526 - t470 * t525;
t351 = t370 * qJD(1) + t504;
t393 = -t520 + t521;
t531 = qJD(3) * t466;
t568 = (qJD(3) * t400 + t373) * t470 + t466 * t371 - t411 * t531;
t312 = -qJ(4) * t351 - qJD(4) * t393 + t568;
t464 = sin(pkin(9));
t557 = cos(pkin(9));
t294 = t557 * t308 - t464 * t312;
t292 = -pkin(4) * t458 - t294;
t443 = cos(t454);
t562 = g(3) * t443;
t572 = t292 + t562;
t457 = t471 * pkin(2);
t560 = pkin(1) + t457;
t362 = -t557 * t393 + t395 * t464;
t492 = -t464 * t393 - t557 * t395;
t565 = pkin(3) * t395;
t332 = pkin(4) * t492 - pkin(8) * t362 - t565;
t358 = qJD(5) - t362;
t444 = pkin(3) * t464 + pkin(8);
t569 = (qJD(5) * t444 + t332) * t358;
t389 = t395 * qJ(4);
t396 = t466 * t411;
t514 = t470 * t400 - t396;
t348 = t389 + t514;
t536 = -t466 * t421 + t470 * t422;
t404 = t466 * t467 - t470 * t471;
t522 = qJD(2) * t566;
t410 = t467 * t522;
t412 = t471 * t522;
t530 = qJD(3) * t470;
t489 = -t470 * t410 - t466 * t412 - t421 * t530 - t422 * t531;
t328 = -qJ(4) * t370 - qJD(4) * t404 + t489;
t369 = t460 * t404;
t479 = -t536 * qJD(3) + t410 * t466 - t412 * t470;
t476 = qJ(4) * t369 - qJD(4) * t405 + t479;
t304 = t557 * t328 + t464 * t476;
t343 = pkin(3) * t460 + t348;
t556 = qJ(4) * t393;
t349 = -t501 - t556;
t546 = t464 * t349;
t317 = t557 * t343 - t546;
t315 = -t460 * pkin(4) - t317;
t324 = t350 * t464 + t557 * t351;
t321 = qJDD(5) + t324;
t366 = t557 * t404 + t405 * t464;
t367 = -t464 * t404 + t557 * t405;
t375 = pkin(3) * t404 - t560;
t333 = pkin(4) * t366 - pkin(8) * t367 + t375;
t359 = -qJ(4) * t404 + t536;
t499 = -t421 * t470 - t422 * t466;
t486 = -qJ(4) * t405 + t499;
t335 = t557 * t359 + t464 * t486;
t341 = -t557 * t369 - t464 * t370;
t295 = t464 * t308 + t557 * t312;
t293 = pkin(8) * t458 + t295;
t420 = t560 * qJD(1);
t372 = pkin(3) * t393 + qJD(4) - t420;
t327 = -pkin(4) * t362 - pkin(8) * t492 + t372;
t511 = qJD(5) * t327 + t293;
t567 = t292 * t367 + t315 * t341 - t335 * t321 - (qJD(5) * t333 + t304) * t358 - t511 * t366;
t430 = g(3) * t442;
t559 = pkin(2) * qJD(3);
t325 = t557 * t350 - t464 * t351;
t465 = sin(qJ(5));
t469 = cos(qJ(5));
t528 = qJD(5) * t469;
t524 = t469 * t325 + t465 * t458 + t460 * t528;
t529 = qJD(5) * t465;
t310 = -t492 * t529 + t524;
t555 = t310 * t465;
t554 = t315 * t362;
t553 = t315 * t367;
t552 = t333 * t321;
t353 = -t469 * t460 + t465 * t492;
t551 = t353 * t358;
t550 = t353 * t492;
t355 = t460 * t465 + t469 * t492;
t549 = t355 * t358;
t548 = t355 * t492;
t545 = t464 * t466;
t544 = t465 * t321;
t543 = t465 * t468;
t542 = t465 * t472;
t541 = t468 * t469;
t319 = t469 * t321;
t540 = t469 * t472;
t344 = t557 * t349;
t318 = t464 * t343 + t344;
t537 = -t470 * t409 - t396;
t356 = t389 + t537;
t500 = t409 * t466 - t547;
t487 = t500 + t556;
t516 = t557 * t466;
t539 = -t356 * t464 + t557 * t487 + (t464 * t470 + t516) * t559;
t538 = -t557 * t356 - t464 * t487 + (t557 * t470 - t545) * t559;
t450 = pkin(2) * t470 + pkin(3);
t388 = pkin(2) * t516 + t464 * t450;
t535 = pkin(3) * t456 + t457;
t461 = t467 ^ 2;
t534 = -t471 ^ 2 + t461;
t453 = t467 * t558;
t364 = pkin(3) * t370 + t453;
t513 = t358 * t469;
t390 = pkin(2) * t519 - qJDD(1) * t560;
t339 = pkin(3) * t351 + qJDD(4) + t390;
t300 = pkin(4) * t324 - pkin(8) * t325 + t339;
t316 = pkin(8) * t460 + t318;
t510 = qJD(5) * t316 - t300;
t383 = pkin(8) + t388;
t452 = pkin(2) * t533;
t508 = qJD(5) * t383 + t332 + t452;
t505 = g(1) * t468 - g(2) * t472;
t503 = -t321 * t383 - t554;
t502 = t317 * t362 + t318 * t492;
t302 = t316 * t469 + t327 * t465;
t498 = t302 * t492 + t315 * t528 + t572 * t465;
t301 = -t316 * t465 + t327 * t469;
t497 = -t301 * t492 + t315 * t529 + t469 * t575;
t496 = t319 + (t362 * t465 - t529) * t358;
t494 = -0.2e1 * pkin(1) * t527 - pkin(6) * qJDD(2);
t493 = t341 * t469 - t367 * t529;
t387 = -pkin(2) * t545 + t557 * t450;
t488 = -t362 * t372 + t506 * t443 - t295 + t430;
t323 = t557 * t348 - t546;
t484 = -t321 * t444 + t323 * t358 - t554;
t473 = qJD(2) ^ 2;
t483 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t473 + t505;
t474 = qJD(1) ^ 2;
t482 = pkin(1) * t474 - pkin(6) * qJDD(1) + t506;
t481 = -t372 * t492 + t294 - t562 + t575;
t478 = g(3) * t455 - t420 * t393 + t506 * t456 - t568;
t439 = t469 * t458;
t311 = t355 * qJD(5) + t325 * t465 - t439;
t477 = -t395 * t393 * MDP(11) - t358 * t492 * MDP(26) + ((t310 - t551) * t469 + (-t311 - t549) * t465) * MDP(23) + (t496 + t550) * MDP(25) + (t358 * t513 + t544 - t548) * MDP(24) + (t355 * t513 + t555) * MDP(22) + (t393 * t460 + t350) * MDP(13) + (-t504 + (-qJD(1) * t405 - t395) * t460) * MDP(14) + (-t393 ^ 2 + t395 ^ 2) * MDP(12) + t458 * MDP(15);
t475 = -t420 * t395 + t480 + t574;
t459 = -qJ(4) - t566;
t445 = -t557 * pkin(3) - pkin(4);
t408 = pkin(1) + t535;
t382 = -pkin(4) - t387;
t381 = t443 * t540 + t543;
t380 = -t443 * t542 + t541;
t379 = -t443 * t541 + t542;
t378 = t443 * t543 + t540;
t374 = t452 - t565;
t340 = -t369 * t464 + t557 * t370;
t334 = t359 * t464 - t557 * t486;
t322 = t348 * t464 + t344;
t307 = pkin(4) * t340 - pkin(8) * t341 + t364;
t303 = t328 * t464 - t557 * t476;
t299 = t469 * t300;
t1 = [(-t304 * t460 + t325 * t375 - t335 * t458 + t339 * t367 + t341 * t372 + t364 * t492 - t505 * t442) * MDP(19) + (-t370 * t460 - t404 * t458) * MDP(14) + (-t351 * t560 - t420 * t370 + t390 * t404 + t393 * t453 + t505 * t456 + t499 * t458 + t479 * t460) * MDP(16) + (-t294 * t367 - t295 * t366 + t303 * t492 + t304 * t362 - t317 * t341 - t318 * t340 - t324 * t335 + t325 * t334 - t506) * MDP(20) + (-t303 * t460 + t324 * t375 - t334 * t458 + t339 * t366 + t340 * t372 - t362 * t364 + t505 * t443) * MDP(18) + (-t350 * t560 + t420 * t369 + t390 * t405 - t395 * t453 - t505 * t455 - t536 * t458 - t489 * t460) * MDP(17) + (-t369 * t460 + t405 * t458) * MDP(13) + (t350 * t405 + t369 * t395) * MDP(11) + (-t350 * t404 - t351 * t405 + t369 * t393 + t370 * t395) * MDP(12) + 0.2e1 * (t467 * t525 - t534 * t527) * MDP(5) + (qJDD(2) * t467 + t471 * t473) * MDP(6) + (qJDD(2) * t471 - t467 * t473) * MDP(7) + (t295 * t335 + t318 * t304 - t294 * t334 - t317 * t303 + t339 * t375 + t372 * t364 - g(1) * (-t408 * t468 - t459 * t472) - g(2) * (t408 * t472 - t459 * t468)) * MDP(21) + (-g(1) * t378 - g(2) * t380 - t302 * t340 + t303 * t355 + t334 * t310 + (-(-qJD(5) * t335 + t307) * t358 - t552 + t510 * t366 - qJD(5) * t553) * t465 + t567 * t469) * MDP(28) + (-g(1) * t379 - g(2) * t381 + t299 * t366 + t301 * t340 + t303 * t353 + t334 * t311 + (t307 * t358 + t552 + (-t316 * t366 - t335 * t358 + t553) * qJD(5)) * t469 + t567 * t465) * MDP(27) + (t310 * t367 * t469 + t493 * t355) * MDP(22) + (t494 * t467 + t483 * t471) * MDP(9) + (-t483 * t467 + t494 * t471) * MDP(10) + ((-t353 * t469 - t355 * t465) * t341 + (-t555 - t311 * t469 + (t353 * t465 - t355 * t469) * qJD(5)) * t367) * MDP(23) + t505 * MDP(2) + t506 * MDP(3) + qJDD(1) * MDP(1) + (t321 * t366 + t340 * t358) * MDP(26) + (qJDD(1) * t461 + 0.2e1 * t467 * t518) * MDP(4) + (t310 * t366 + t367 * t319 + t340 * t355 + t493 * t358) * MDP(24) + (-t367 * t544 - t311 * t366 - t340 * t353 + (-t341 * t465 - t367 * t528) * t358) * MDP(25); (t295 * t388 + t294 * t387 - t372 * t374 - g(3) * t535 - t506 * (-pkin(2) * t467 - pkin(3) * t455) + t538 * t318 - t539 * t317) * MDP(21) + (t537 * t460 + (t395 * t533 - t466 * t458 - t460 * t530) * pkin(2) + t478) * MDP(17) + (t362 * t374 + t387 * t458 - t539 * t460 + t481) * MDP(18) + (-t374 * t492 - t388 * t458 - t538 * t460 + t488) * MDP(19) + (-t324 * t388 - t325 * t387 + t362 * t538 + t492 * t539 + t502) * MDP(20) + (-g(3) * t471 + t482 * t467) * MDP(9) + (g(3) * t467 + t482 * t471) * MDP(10) + (t382 * t310 + t503 * t469 - t465 * t575 + t539 * t355 + (t508 * t465 - t538 * t469) * t358 + t498) * MDP(28) + t477 + (-t500 * t460 + (-t393 * t533 + t470 * t458 - t460 * t531) * pkin(2) + t475) * MDP(16) + qJDD(2) * MDP(8) + MDP(6) * t526 + MDP(7) * t525 + (t382 * t311 - t572 * t469 + t503 * t465 + t539 * t353 + (-t538 * t465 - t508 * t469) * t358 + t497) * MDP(27) + (-t467 * t471 * MDP(4) + t534 * MDP(5)) * t474; (t323 * t460 + (t395 * t492 - t458 * t464) * pkin(3) + t488) * MDP(19) + (-t322 * t492 - t323 * t362 + (-t324 * t464 - t557 * t325) * pkin(3) + t502) * MDP(20) + (t445 * t311 - t322 * t353 + t484 * t465 + (-t572 - t569) * t469 + t497) * MDP(27) + (t445 * t310 - t322 * t355 + t484 * t469 + (-t575 + t569) * t465 + t498) * MDP(28) + (t317 * t322 - t318 * t323 + (t557 * t294 + t295 * t464 + t372 * t395 + t574) * pkin(3)) * MDP(21) + t477 + (-t501 * t460 + t475) * MDP(16) + (t322 * t460 + (-t362 * t395 + t557 * t458) * pkin(3) + t481) * MDP(18) + (t514 * t460 + t478) * MDP(17); (t460 * t492 + t324) * MDP(18) + (t362 * t460 + t325) * MDP(19) + (-t362 ^ 2 - t492 ^ 2) * MDP(20) + (t317 * t492 - t318 * t362 + t339 - t505) * MDP(21) + (t496 - t550) * MDP(27) + (-t358 ^ 2 * t469 - t544 - t548) * MDP(28); t355 * t353 * MDP(22) + (-t353 ^ 2 + t355 ^ 2) * MDP(23) + (t524 + t551) * MDP(24) + (t439 + t549) * MDP(25) + t321 * MDP(26) + (-g(1) * t380 + g(2) * t378 + t302 * t358 - t315 * t355 + t299) * MDP(27) + (g(1) * t381 - g(2) * t379 + t301 * t358 + t315 * t353) * MDP(28) + ((-t293 + t430) * MDP(28) + (-MDP(25) * t492 - MDP(27) * t316 - MDP(28) * t327) * qJD(5)) * t469 + (-qJD(5) * t492 * MDP(24) + (-qJD(5) * t460 - t325) * MDP(25) + (-t511 + t430) * MDP(27) + t510 * MDP(28)) * t465;];
tau = t1;
