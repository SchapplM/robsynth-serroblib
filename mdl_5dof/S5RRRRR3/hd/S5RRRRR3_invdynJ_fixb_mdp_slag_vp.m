% Calculate vector of inverse dynamics joint torques for
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:02
% EndTime: 2019-07-18 17:19:14
% DurationCPUTime: 7.24s
% Computational Cost: add. (3966->423), mult. (8791->567), div. (0->0), fcn. (7088->14), ass. (0->189)
t438 = qJ(2) + qJ(3);
t430 = sin(t438);
t443 = sin(qJ(1));
t447 = cos(qJ(1));
t476 = g(1) * t447 + g(2) * t443;
t466 = t476 * t430;
t432 = cos(t438);
t552 = g(3) * t432;
t573 = -t552 + t466;
t441 = sin(qJ(3));
t442 = sin(qJ(2));
t556 = cos(qJ(3));
t491 = qJD(1) * t556;
t446 = cos(qJ(2));
t519 = qJD(1) * t446;
t386 = -t441 * t519 - t442 * t491;
t503 = pkin(1) * t519;
t520 = qJD(1) * t442;
t576 = -t441 * t520 + t446 * t491;
t346 = -pkin(2) * t576 + pkin(5) * t386 - t503;
t555 = pkin(1) * t441;
t428 = qJD(2) * t555;
t434 = qJD(2) + qJD(3);
t394 = pkin(5) * t434 + t428;
t440 = sin(qJ(4));
t445 = cos(qJ(4));
t331 = t445 * t346 - t394 * t440;
t381 = qJD(4) - t576;
t317 = pkin(3) * t381 + t331;
t439 = sin(qJ(5));
t332 = t346 * t440 + t394 * t445;
t444 = cos(qJ(5));
t547 = t332 * t444;
t305 = t317 * t439 + t547;
t486 = qJDD(1) * t556;
t507 = qJDD(1) * t446;
t577 = t434 * t576;
t338 = t441 * t507 + t442 * t486 + t577;
t433 = qJDD(2) + qJDD(3);
t470 = t386 * t445 - t434 * t440;
t316 = -qJD(4) * t470 + t338 * t440 - t445 * t433;
t504 = t556 * pkin(1);
t524 = -qJD(3) * t428 + qJDD(2) * t504;
t373 = -pkin(2) * t433 - t524;
t311 = pkin(3) * t316 + t373;
t360 = -t386 * t440 - t445 * t434;
t490 = t556 * qJD(2);
t479 = pkin(1) * t490;
t395 = -t434 * pkin(2) - t479;
t347 = t360 * pkin(3) + t395;
t530 = t440 * t444;
t532 = t439 * t445;
t391 = t530 + t532;
t437 = qJ(4) + qJ(5);
t429 = sin(t437);
t533 = t439 * t440;
t389 = -t444 * t445 + t533;
t558 = qJD(4) + qJD(5);
t354 = t558 * t389;
t525 = -t389 * t576 + t354;
t578 = -t305 * t386 + t311 * t391 - t347 * t525 - t429 * t573;
t515 = qJD(4) * t440;
t575 = (-t440 * t576 + t515) * pkin(3);
t574 = -qJD(5) * t440 - t515;
t541 = t470 * t439;
t327 = t444 * t360 - t541;
t375 = qJD(5) + t381;
t572 = t327 * t375;
t471 = t360 * t439 + t444 * t470;
t571 = t375 * t471;
t570 = t558 * t391;
t569 = -t391 * t576 + t570;
t568 = t373 + t552;
t517 = qJD(3) * t441;
t502 = pkin(1) * t517;
t567 = t502 + t575;
t513 = qJD(5) * t439;
t324 = t332 * t513;
t431 = cos(t437);
t535 = t432 * t443;
t366 = t429 * t447 - t431 * t535;
t534 = t432 * t447;
t368 = t429 * t443 + t431 * t534;
t420 = g(3) * t430;
t566 = g(1) * t368 - g(2) * t366 + t327 * t347 + t431 * t420 + t324;
t508 = qJDD(1) * t442;
t473 = t441 * t508 - t446 * t486;
t392 = t441 * t446 + t442 * t556;
t357 = t434 * t392;
t521 = qJD(1) * t357;
t339 = t473 + t521;
t509 = qJD(1) * qJD(2);
t488 = t442 * t509;
t310 = pkin(2) * t339 - pkin(5) * t338 + (t488 - t507) * pkin(1);
t309 = t445 * t310;
t336 = qJDD(4) + t339;
t489 = t556 * qJD(3);
t506 = qJDD(2) * t441;
t374 = t433 * pkin(5) + (qJD(2) * t489 + t506) * pkin(1);
t293 = pkin(3) * t336 - qJD(4) * t332 - t374 * t440 + t309;
t292 = t444 * t293;
t298 = qJD(4) * t331 + t310 * t440 + t374 * t445;
t365 = t429 * t535 + t431 * t447;
t367 = -t429 * t534 + t431 * t443;
t565 = -g(1) * t367 + g(2) * t365 - t439 * t298 + t347 * t471 + t429 * t420 + t292;
t335 = qJDD(5) + t336;
t564 = t335 * MDP(29) + (-t327 ^ 2 + t471 ^ 2) * MDP(26) - t327 * MDP(25) * t471;
t351 = -pkin(2) * t386 - pkin(5) * t576;
t345 = pkin(1) * t520 + t351;
t424 = pkin(5) + t555;
t563 = t381 * (qJD(4) * t424 + t345);
t561 = t574 * t439;
t514 = qJD(4) * t445;
t560 = -qJD(5) * t445 - t514;
t559 = qJD(1) * t392;
t315 = t445 * t338 + t386 * t515 + t440 * t433 + t434 * t514;
t485 = t315 * t439 + t444 * t316;
t297 = -qJD(5) * t471 + t485;
t465 = -t441 * t442 + t446 * t556;
t356 = t434 * t465;
t518 = qJD(2) * t442;
t322 = pkin(1) * t518 + pkin(2) * t357 - pkin(5) * t356;
t352 = -pkin(1) * t446 - pkin(2) * t465 - pkin(5) * t392;
t337 = -pkin(3) * t465 + t352 * t445;
t480 = qJD(5) * t317 + t298;
t557 = -t375 * (qJD(5) * t337 + t322 * t440 + t352 * t514) + t465 * t480 - t352 * t440 * t335;
t551 = t445 * pkin(3);
t550 = pkin(1) * qJD(3);
t549 = t315 * t440;
t548 = t317 * t444;
t546 = t335 * t389;
t545 = t335 * t391;
t544 = t356 * t440;
t543 = t360 * t381;
t542 = t470 * t381;
t540 = t392 * t395;
t539 = t392 * t440;
t538 = t392 * t445;
t536 = t395 * t576;
t531 = t440 * t443;
t529 = t440 * t447;
t528 = t443 * t445;
t527 = t445 * t447;
t435 = t442 ^ 2;
t523 = -t446 ^ 2 + t435;
t516 = qJD(4) * t381;
t511 = qJD(5) * t444;
t500 = pkin(5) * t516;
t498 = t444 * t315 - t439 * t316 - t360 * t511;
t497 = t440 * t556;
t496 = t445 * t556;
t495 = t556 * t434;
t483 = t381 * t445;
t481 = -qJD(4) * t346 - t374;
t478 = g(1) * t534 + g(2) * t535 + t503 * t576 + t420;
t425 = -t504 - pkin(2);
t477 = -t428 + t575;
t474 = -pkin(5) * t336 - t536;
t472 = t322 * t381 + t352 * t336;
t469 = t331 * t386 + t395 * t515 + t445 * t466;
t468 = -t332 * t386 + t395 * t514 + t440 * t568;
t467 = t481 + t420;
t464 = t392 * t514 + t544;
t463 = t356 * t445 - t392 * t515;
t296 = t470 * t513 + t498;
t459 = -t386 * t503 + t524 + t573;
t457 = -t352 * t516 + t395 * t356 + t373 * t392;
t454 = -t337 * t335 - (pkin(3) * t357 + t322 * t445 + t352 * t574) * t375;
t453 = -pkin(1) * t381 * t489 - t336 * t424 - t536;
t304 = -t332 * t439 + t548;
t452 = t304 * t386 + t311 * t389 + t569 * t347 + t431 * t573;
t450 = (-t296 * t389 - t297 * t391 + t327 * t525 + t471 * t569) * MDP(26) + (t296 * t391 + t471 * t525) * MDP(25) + ((t315 - t543) * t445 + (-t316 + t542) * t440) * MDP(19) + (-t375 * t525 - t386 * t471 + t545) * MDP(27) + (-t327 * t386 - t375 * t569 - t546) * MDP(28) + (-t470 * t483 + t549) * MDP(18) + (-t381 ^ 2 * t440 + t336 * t445 - t360 * t386) * MDP(21) + (t336 * t440 + t381 * t483 - t386 * t470) * MDP(20) + (t338 - t577) * MDP(13) + (-t473 + (-t386 - t559) * t434) * MDP(14) + (t386 ^ 2 - t576 ^ 2) * MDP(12) + t433 * MDP(15) + (MDP(11) * t576 + MDP(22) * t381 + MDP(29) * t375) * t386;
t448 = qJD(2) ^ 2;
t426 = -pkin(2) - t551;
t401 = t425 - t551;
t380 = t386 * pkin(3);
t379 = t432 * t527 + t531;
t378 = -t432 * t529 + t528;
t377 = -t432 * t528 + t529;
t376 = t432 * t531 + t527;
t350 = t389 * t392;
t349 = t391 * t392;
t348 = t445 * t351;
t341 = t440 * t351 + t445 * t479;
t325 = t345 * t445 - t380;
t323 = -t440 * t479 + t348 - t380;
t303 = t356 * t532 + (t538 * t558 + t544) * t444 + t561 * t392;
t302 = -t356 * t389 - t392 * t570;
t1 = [t476 * MDP(3) + (qJDD(1) * t435 + 0.2e1 * t446 * t488) * MDP(4) + (-MDP(10) * t442 + MDP(16) * t432 - MDP(17) * t430 + MDP(9) * t446 + MDP(2)) * (g(1) * t443 - g(2) * t447) + (t338 * t465 - t339 * t392 + t356 * t576 + t357 * t386) * MDP(12) + (((-t386 + t559) * t518 + (-qJD(1) * t356 - qJDD(1) * t392 - t338) * t446) * MDP(17) + ((-qJD(1) * t465 - t576) * t518 + (qJDD(1) * t465 - t339 - t521) * t446) * MDP(16)) * pkin(1) + qJDD(1) * MDP(1) + (qJDD(2) * t442 + t446 * t448) * MDP(6) + (qJDD(2) * t446 - t442 * t448) * MDP(7) + (t315 * t538 - t463 * t470) * MDP(18) + ((-t360 * t445 + t440 * t470) * t356 + (-t549 - t316 * t445 + (t360 * t440 + t445 * t470) * qJD(4)) * t392) * MDP(19) + (-g(1) * t366 - g(2) * t368 - t292 * t465 + t347 * t303 + t304 * t357 + t311 * t349 + (qJD(5) * t332 * t465 - t454) * t444 + t557 * t439 + (t297 * t539 + t327 * t464) * pkin(3)) * MDP(30) + (-t357 * t434 + t433 * t465) * MDP(14) + (t297 * t465 - t303 * t375 - t327 * t357 - t335 * t349) * MDP(28) + (-t335 * t465 + t357 * t375) * MDP(29) + (-t336 * t465 + t357 * t381) * MDP(22) + (t316 * t465 - t336 * t539 - t357 * t360 - t381 * t464) * MDP(21) + (-g(1) * t376 - g(2) * t378 + t298 * t465 - t332 * t357 + t457 * t445 + (-qJD(4) * t540 - t472) * t440) * MDP(24) + (-g(1) * t377 - g(2) * t379 - t309 * t465 + t331 * t357 + ((t394 * t465 + t540) * qJD(4) + t472) * t445 + (-t465 * t481 + t457) * t440) * MDP(23) + (-t315 * t465 + t336 * t538 - t357 * t470 + t381 * t463) * MDP(20) + (-g(1) * t365 - g(2) * t367 + t347 * t302 - t305 * t357 - t311 * t350 - t324 * t465 + (t293 * t465 + t454) * t439 + t557 * t444 + (t296 * t539 - t464 * t471) * pkin(3)) * MDP(31) + (-t296 * t465 + t302 * t375 - t335 * t350 - t357 * t471) * MDP(27) + (-t296 * t349 + t297 * t350 - t302 * t327 + t303 * t471) * MDP(26) + (-t296 * t350 - t302 * t471) * MDP(25) + (t356 * t434 + t392 * t433) * MDP(13) + (t338 * t392 - t356 * t386) * MDP(11) + 0.2e1 * (t442 * t507 - t509 * t523) * MDP(5); (t360 * t502 + t425 * t316 + t453 * t440 + (-t568 - t563) * t445 + t469) * MDP(23) + (-t470 * t502 + t425 * t315 + t453 * t445 + (-t466 + t563) * t440 + t468) * MDP(24) + ((t386 * t520 + (-qJDD(2) - t433) * t441 + (-t490 - t495) * qJD(3)) * pkin(1) + t478) * MDP(17) + ((t433 * t556 - t434 * t517 + t520 * t576) * pkin(1) + t459) * MDP(16) + (t401 * t297 - t424 * t545 + t452 + ((-t439 * t496 - t444 * t497) * t550 + t424 * t354 - t325 * t444 + t345 * t533) * t375 + t567 * t327) * MDP(30) + (t401 * t296 + t424 * t546 + (-(-t439 * t497 + t444 * t496) * t550 + t570 * t424 + t325 * t439 + t345 * t530) * t375 - t567 * t471 + t578) * MDP(31) + (-g(3) * t446 + t442 * t476) * MDP(9) + (g(3) * t442 + t446 * t476) * MDP(10) + qJDD(2) * MDP(8) + t450 + MDP(7) * t507 + MDP(6) * t508 + (-MDP(4) * t442 * t446 + MDP(5) * t523) * qJD(1) ^ 2; (t426 * t297 - (t323 * t444 - t341 * t439) * t375 + t477 * t327 + ((t444 * t560 - t561) * t375 - t545) * pkin(5) + t452) * MDP(30) + (t426 * t296 + (t323 * t439 + t341 * t444) * t375 - t477 * t471 + (-(t439 * t560 - t440 * t511 - t444 * t515) * t375 + t546) * pkin(5) + t578) * MDP(31) + (-t360 * t428 - pkin(2) * t316 - t348 * t381 + (t381 * t479 + t474) * t440 + (-t568 - t500) * t445 + t469) * MDP(23) + (t470 * t428 - pkin(2) * t315 + t341 * t381 + t474 * t445 + (-t466 + t500) * t440 + t468) * MDP(24) + ((-t506 + (-t489 + t495) * qJD(2)) * pkin(1) + t478) * MDP(17) + (t428 * t434 + t459) * MDP(16) + t450; -t470 * t360 * MDP(18) + (-t360 ^ 2 + t470 ^ 2) * MDP(19) + (t315 + t543) * MDP(20) + (-t316 - t542) * MDP(21) + t336 * MDP(22) + (-g(1) * t378 + g(2) * t376 + t332 * t381 - t394 * t514 + t395 * t470 + t440 * t467 + t309) * MDP(23) + (g(1) * t379 - g(2) * t377 + t331 * t381 + t360 * t395 + (qJD(4) * t394 - t310) * t440 + t467 * t445) * MDP(24) + (t296 + t572) * MDP(27) + (-t297 - t571) * MDP(28) + (-(-t331 * t439 - t547) * t375 - t305 * qJD(5) + (t327 * t470 + t444 * t335 - t375 * t513) * pkin(3) + t565) * MDP(30) + ((-t332 * t375 - t293) * t439 + (t331 * t375 - t480) * t444 + (-t439 * t335 - t375 * t511 - t470 * t471) * pkin(3) + t566) * MDP(31) + t564; (t498 + t572) * MDP(27) + (-t485 - t571) * MDP(28) + (t305 * t375 + t565) * MDP(30) + (-t439 * t293 - t444 * t298 + t304 * t375 + t566) * MDP(31) + (MDP(27) * t541 + MDP(28) * t471 - MDP(30) * t305 - MDP(31) * t548) * qJD(5) + t564;];
tau  = t1;
