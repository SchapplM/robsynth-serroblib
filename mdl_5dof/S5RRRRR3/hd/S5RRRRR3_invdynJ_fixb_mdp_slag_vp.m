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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:56:44
% EndTime: 2019-12-05 18:56:55
% DurationCPUTime: 7.52s
% Computational Cost: add. (3966->423), mult. (8791->567), div. (0->0), fcn. (7088->14), ass. (0->189)
t437 = qJ(2) + qJ(3);
t429 = sin(t437);
t442 = sin(qJ(1));
t446 = cos(qJ(1));
t475 = g(1) * t446 + g(2) * t442;
t465 = t475 * t429;
t431 = cos(t437);
t551 = g(3) * t431;
t578 = t551 - t465;
t440 = sin(qJ(3));
t441 = sin(qJ(2));
t555 = cos(qJ(3));
t490 = qJD(1) * t555;
t445 = cos(qJ(2));
t519 = qJD(1) * t445;
t385 = -t440 * t519 - t441 * t490;
t502 = pkin(1) * t519;
t520 = qJD(1) * t441;
t575 = -t440 * t520 + t445 * t490;
t345 = -pkin(2) * t575 + pkin(5) * t385 - t502;
t554 = pkin(1) * t440;
t427 = qJD(2) * t554;
t433 = qJD(2) + qJD(3);
t393 = pkin(5) * t433 + t427;
t439 = sin(qJ(4));
t444 = cos(qJ(4));
t330 = t444 * t345 - t393 * t439;
t380 = qJD(4) - t575;
t316 = pkin(3) * t380 + t330;
t438 = sin(qJ(5));
t331 = t345 * t439 + t393 * t444;
t443 = cos(qJ(5));
t546 = t331 * t443;
t304 = t316 * t438 + t546;
t485 = qJDD(1) * t555;
t506 = qJDD(1) * t445;
t576 = t575 * t433;
t337 = t440 * t506 + t441 * t485 + t576;
t432 = qJDD(2) + qJDD(3);
t469 = t385 * t444 - t433 * t439;
t315 = -qJD(4) * t469 + t337 * t439 - t444 * t432;
t503 = t555 * pkin(1);
t523 = -qJD(3) * t427 + qJDD(2) * t503;
t372 = -pkin(2) * t432 - t523;
t310 = pkin(3) * t315 + t372;
t359 = -t385 * t439 - t444 * t433;
t489 = t555 * qJD(2);
t478 = pkin(1) * t489;
t394 = -t433 * pkin(2) - t478;
t346 = t359 * pkin(3) + t394;
t529 = t439 * t443;
t531 = t438 * t444;
t390 = t529 + t531;
t436 = qJ(4) + qJ(5);
t428 = sin(t436);
t532 = t438 * t439;
t388 = -t443 * t444 + t532;
t557 = qJD(4) + qJD(5);
t353 = t557 * t388;
t524 = -t388 * t575 + t353;
t577 = -t304 * t385 + t310 * t390 - t346 * t524 + t578 * t428;
t515 = qJD(4) * t439;
t574 = (-t439 * t575 + t515) * pkin(3);
t573 = -qJD(5) * t439 - t515;
t540 = t469 * t438;
t326 = t443 * t359 - t540;
t374 = qJD(5) + t380;
t571 = t326 * t374;
t470 = t359 * t438 + t443 * t469;
t570 = t374 * t470;
t569 = t557 * t390;
t568 = -t390 * t575 + t569;
t567 = t372 + t551;
t517 = qJD(3) * t440;
t501 = pkin(1) * t517;
t566 = t501 + t574;
t513 = qJD(5) * t438;
t323 = t331 * t513;
t430 = cos(t436);
t534 = t431 * t442;
t365 = t428 * t446 - t430 * t534;
t533 = t431 * t446;
t367 = t428 * t442 + t430 * t533;
t419 = g(3) * t429;
t565 = g(1) * t367 - g(2) * t365 + t326 * t346 + t430 * t419 + t323;
t507 = qJDD(1) * t441;
t472 = t440 * t507 - t445 * t485;
t391 = t440 * t445 + t441 * t555;
t356 = t433 * t391;
t521 = qJD(1) * t356;
t338 = t472 + t521;
t508 = qJD(1) * qJD(2);
t487 = t441 * t508;
t309 = pkin(2) * t338 - pkin(5) * t337 + (t487 - t506) * pkin(1);
t308 = t444 * t309;
t335 = qJDD(4) + t338;
t488 = t555 * qJD(3);
t505 = qJDD(2) * t440;
t373 = t432 * pkin(5) + (qJD(2) * t488 + t505) * pkin(1);
t292 = pkin(3) * t335 - qJD(4) * t331 - t373 * t439 + t308;
t291 = t443 * t292;
t297 = qJD(4) * t330 + t309 * t439 + t373 * t444;
t364 = t428 * t534 + t430 * t446;
t366 = -t428 * t533 + t430 * t442;
t564 = -g(1) * t366 + g(2) * t364 - t438 * t297 + t346 * t470 + t428 * t419 + t291;
t334 = qJDD(5) + t335;
t563 = t334 * MDP(29) + (-t326 ^ 2 + t470 ^ 2) * MDP(26) - t326 * t470 * MDP(25);
t350 = -pkin(2) * t385 - pkin(5) * t575;
t344 = pkin(1) * t520 + t350;
t423 = pkin(5) + t554;
t562 = t380 * (qJD(4) * t423 + t344);
t560 = t573 * t438;
t514 = qJD(4) * t444;
t559 = -qJD(5) * t444 - t514;
t558 = qJD(1) * t391;
t314 = t444 * t337 + t385 * t515 + t439 * t432 + t433 * t514;
t484 = t314 * t438 + t443 * t315;
t296 = -qJD(5) * t470 + t484;
t464 = -t440 * t441 + t445 * t555;
t355 = t433 * t464;
t518 = qJD(2) * t441;
t321 = pkin(1) * t518 + pkin(2) * t356 - pkin(5) * t355;
t351 = -pkin(1) * t445 - pkin(2) * t464 - pkin(5) * t391;
t336 = -pkin(3) * t464 + t351 * t444;
t479 = qJD(5) * t316 + t297;
t556 = -t374 * (qJD(5) * t336 + t321 * t439 + t351 * t514) + t464 * t479 - t351 * t439 * t334;
t550 = t444 * pkin(3);
t549 = pkin(1) * qJD(3);
t548 = t314 * t439;
t547 = t316 * t443;
t545 = t334 * t388;
t544 = t334 * t390;
t543 = t355 * t439;
t542 = t359 * t380;
t541 = t469 * t380;
t539 = t391 * t394;
t538 = t391 * t439;
t537 = t391 * t444;
t535 = t394 * t575;
t530 = t439 * t442;
t528 = t439 * t446;
t527 = t442 * t444;
t526 = t444 * t446;
t434 = t441 ^ 2;
t522 = -t445 ^ 2 + t434;
t516 = qJD(4) * t380;
t511 = qJD(5) * t443;
t499 = pkin(5) * t516;
t497 = t443 * t314 - t438 * t315 - t359 * t511;
t496 = t439 * t555;
t495 = t444 * t555;
t494 = t555 * t433;
t482 = t380 * t444;
t480 = -qJD(4) * t345 - t373;
t477 = g(1) * t533 + g(2) * t534 + t502 * t575 + t419;
t424 = -t503 - pkin(2);
t476 = -t427 + t574;
t473 = -pkin(5) * t335 - t535;
t471 = t321 * t380 + t351 * t335;
t468 = t330 * t385 + t394 * t515 + t444 * t465;
t467 = -t331 * t385 + t394 * t514 + t439 * t567;
t466 = t480 + t419;
t463 = t391 * t514 + t543;
t462 = t355 * t444 - t391 * t515;
t295 = t469 * t513 + t497;
t458 = -t385 * t502 + t523 - t578;
t456 = -t351 * t516 + t394 * t355 + t372 * t391;
t453 = -t336 * t334 - (pkin(3) * t356 + t321 * t444 + t573 * t351) * t374;
t452 = -pkin(1) * t380 * t488 - t335 * t423 - t535;
t303 = -t331 * t438 + t547;
t451 = t303 * t385 + t310 * t388 + t568 * t346 - t430 * t578;
t449 = (-t295 * t388 - t296 * t390 + t326 * t524 + t470 * t568) * MDP(26) + (t295 * t390 + t470 * t524) * MDP(25) + ((t314 - t542) * t444 + (-t315 + t541) * t439) * MDP(19) + (-t374 * t524 - t385 * t470 + t544) * MDP(27) + (-t326 * t385 - t374 * t568 - t545) * MDP(28) + (-t469 * t482 + t548) * MDP(18) + (-t380 ^ 2 * t439 + t335 * t444 - t359 * t385) * MDP(21) + (t335 * t439 + t380 * t482 - t385 * t469) * MDP(20) + (t337 - t576) * MDP(13) + (-t472 + (-t385 - t558) * t433) * MDP(14) + (t385 ^ 2 - t575 ^ 2) * MDP(12) + t432 * MDP(15) + (MDP(11) * t575 + MDP(22) * t380 + MDP(29) * t374) * t385;
t447 = qJD(2) ^ 2;
t425 = -pkin(2) - t550;
t400 = t424 - t550;
t379 = t385 * pkin(3);
t378 = t431 * t526 + t530;
t377 = -t431 * t528 + t527;
t376 = -t431 * t527 + t528;
t375 = t431 * t530 + t526;
t349 = t388 * t391;
t348 = t390 * t391;
t347 = t444 * t350;
t340 = t439 * t350 + t444 * t478;
t324 = t344 * t444 - t379;
t322 = -t439 * t478 + t347 - t379;
t302 = t355 * t531 + (t537 * t557 + t543) * t443 + t560 * t391;
t301 = -t355 * t388 - t391 * t569;
t1 = [(((-t385 + t558) * t518 + (-qJD(1) * t355 - qJDD(1) * t391 - t337) * t445) * MDP(17) + ((-qJD(1) * t464 - t575) * t518 + (qJDD(1) * t464 - t338 - t521) * t445) * MDP(16)) * pkin(1) + (t337 * t464 - t338 * t391 + t355 * t575 + t356 * t385) * MDP(12) + 0.2e1 * (t441 * t506 - t508 * t522) * MDP(5) + (qJDD(2) * t441 + t445 * t447) * MDP(6) + (qJDD(2) * t445 - t441 * t447) * MDP(7) + ((-t359 * t444 + t439 * t469) * t355 + (-t548 - t315 * t444 + (t359 * t439 + t444 * t469) * qJD(4)) * t391) * MDP(19) + (t314 * t537 - t462 * t469) * MDP(18) + (-t314 * t464 + t335 * t537 - t356 * t469 + t380 * t462) * MDP(20) + (-g(1) * t364 - g(2) * t366 + t346 * t301 - t304 * t356 - t310 * t349 - t323 * t464 + (t292 * t464 + t453) * t438 + t556 * t443 + (t295 * t538 - t463 * t470) * pkin(3)) * MDP(31) + (-t295 * t464 + t301 * t374 - t334 * t349 - t356 * t470) * MDP(27) + (-g(1) * t365 - g(2) * t367 - t291 * t464 + t346 * t302 + t303 * t356 + t310 * t348 + (qJD(5) * t331 * t464 - t453) * t443 + t556 * t438 + (t296 * t538 + t326 * t463) * pkin(3)) * MDP(30) + (-g(1) * t375 - g(2) * t377 + t297 * t464 - t331 * t356 + t456 * t444 + (-qJD(4) * t539 - t471) * t439) * MDP(24) + (-g(1) * t376 - g(2) * t378 - t308 * t464 + t330 * t356 + ((t393 * t464 + t539) * qJD(4) + t471) * t444 + (-t464 * t480 + t456) * t439) * MDP(23) + (t315 * t464 - t335 * t538 - t356 * t359 - t380 * t463) * MDP(21) + (-t356 * t433 + t432 * t464) * MDP(14) + (-t335 * t464 + t356 * t380) * MDP(22) + (t296 * t464 - t302 * t374 - t326 * t356 - t334 * t348) * MDP(28) + (-t334 * t464 + t356 * t374) * MDP(29) + (-t295 * t348 + t296 * t349 - t301 * t326 + t302 * t470) * MDP(26) + (-t295 * t349 - t301 * t470) * MDP(25) + (t355 * t433 + t391 * t432) * MDP(13) + (t337 * t391 - t355 * t385) * MDP(11) + qJDD(1) * MDP(1) + (qJDD(1) * t434 + 0.2e1 * t445 * t487) * MDP(4) + t475 * MDP(3) + (-MDP(10) * t441 + MDP(16) * t431 - MDP(17) * t429 + MDP(9) * t445 + MDP(2)) * (g(1) * t442 - g(2) * t446); t449 + (-g(3) * t445 + t441 * t475) * MDP(9) + (g(3) * t441 + t445 * t475) * MDP(10) + (t400 * t296 - t423 * t544 + t451 + ((-t438 * t495 - t443 * t496) * t549 + t423 * t353 - t324 * t443 + t344 * t532) * t374 + t566 * t326) * MDP(30) + (t400 * t295 + t423 * t545 + (-(-t438 * t496 + t443 * t495) * t549 + t569 * t423 + t324 * t438 + t344 * t529) * t374 - t566 * t470 + t577) * MDP(31) + ((t385 * t520 + (-qJDD(2) - t432) * t440 + (-t489 - t494) * qJD(3)) * pkin(1) + t477) * MDP(17) + ((t432 * t555 - t433 * t517 + t520 * t575) * pkin(1) + t458) * MDP(16) + qJDD(2) * MDP(8) + MDP(6) * t507 + (t359 * t501 + t424 * t315 + t452 * t439 + (-t567 - t562) * t444 + t468) * MDP(23) + MDP(7) * t506 + (-t469 * t501 + t424 * t314 + t452 * t444 + (-t465 + t562) * t439 + t467) * MDP(24) + (-MDP(4) * t441 * t445 + MDP(5) * t522) * qJD(1) ^ 2; t449 + (-t359 * t427 - pkin(2) * t315 - t347 * t380 + (t380 * t478 + t473) * t439 + (-t567 - t499) * t444 + t468) * MDP(23) + (t469 * t427 - pkin(2) * t314 + t340 * t380 + t473 * t444 + (-t465 + t499) * t439 + t467) * MDP(24) + ((-t505 + (-t488 + t494) * qJD(2)) * pkin(1) + t477) * MDP(17) + (t427 * t433 + t458) * MDP(16) + (t425 * t295 + (t322 * t438 + t340 * t443) * t374 - t476 * t470 + (-(t438 * t559 - t439 * t511 - t443 * t515) * t374 + t545) * pkin(5) + t577) * MDP(31) + (t425 * t296 - (t322 * t443 - t340 * t438) * t374 + t476 * t326 + ((t443 * t559 - t560) * t374 - t544) * pkin(5) + t451) * MDP(30); -t469 * t359 * MDP(18) + (-t359 ^ 2 + t469 ^ 2) * MDP(19) + (t314 + t542) * MDP(20) + (-t315 - t541) * MDP(21) + t335 * MDP(22) + (-g(1) * t377 + g(2) * t375 + t331 * t380 - t393 * t514 + t394 * t469 + t439 * t466 + t308) * MDP(23) + (g(1) * t378 - g(2) * t376 + t330 * t380 + t359 * t394 + (qJD(4) * t393 - t309) * t439 + t466 * t444) * MDP(24) + (t295 + t571) * MDP(27) + (-t296 - t570) * MDP(28) + (-(-t330 * t438 - t546) * t374 - t304 * qJD(5) + (t326 * t469 + t443 * t334 - t374 * t513) * pkin(3) + t564) * MDP(30) + ((-t331 * t374 - t292) * t438 + (t330 * t374 - t479) * t443 + (-t438 * t334 - t374 * t511 - t469 * t470) * pkin(3) + t565) * MDP(31) + t563; (t497 + t571) * MDP(27) + (-t484 - t570) * MDP(28) + (t304 * t374 + t564) * MDP(30) + (-t438 * t292 - t443 * t297 + t303 * t374 + t565) * MDP(31) + (MDP(27) * t540 + MDP(28) * t470 - MDP(30) * t304 - MDP(31) * t547) * qJD(5) + t563;];
tau = t1;
