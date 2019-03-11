% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:57
% EndTime: 2019-03-08 20:08:04
% DurationCPUTime: 6.12s
% Computational Cost: add. (3913->448), mult. (9140->595), div. (0->0), fcn. (7430->14), ass. (0->211)
t437 = sin(pkin(11));
t439 = cos(pkin(11));
t443 = sin(qJ(4));
t446 = cos(qJ(4));
t401 = t437 * t443 - t439 * t446;
t559 = pkin(8) + qJ(3);
t407 = t559 * t437;
t408 = t559 * t439;
t478 = -t407 * t446 - t408 * t443;
t330 = -qJD(3) * t401 + qJD(4) * t478;
t438 = sin(pkin(6));
t447 = cos(qJ(2));
t537 = t438 * t447;
t460 = t401 * t537;
t373 = qJD(1) * t460;
t584 = t330 + t373;
t398 = t401 * qJD(4);
t402 = t437 * t446 + t439 * t443;
t399 = t402 * qJD(4);
t444 = sin(qJ(2));
t520 = qJD(1) * t438;
t502 = t444 * t520;
t583 = pkin(4) * t399 + pkin(9) * t398 - t502;
t404 = qJD(2) * qJ(3) + t502;
t558 = cos(pkin(6));
t491 = qJD(1) * t558;
t419 = t439 * t491;
t365 = t419 + (-pkin(8) * qJD(2) - t404) * t437;
t381 = t404 * t439 + t437 * t491;
t519 = qJD(2) * t439;
t366 = pkin(8) * t519 + t381;
t322 = t365 * t443 + t366 * t446;
t582 = qJD(4) * t322;
t397 = t402 * qJD(2);
t442 = sin(qJ(5));
t445 = cos(qJ(5));
t518 = qJD(2) * t443;
t500 = t437 * t518;
t516 = qJD(2) * t446;
t423 = t439 * t516;
t508 = qJDD(2) * t446;
t509 = qJDD(2) * t443;
t503 = qJD(4) * t423 + t437 * t508 + t439 * t509;
t458 = qJD(4) * t500 - t503;
t512 = t445 * qJD(4);
t514 = qJD(5) * t442;
t315 = -qJD(5) * t512 - qJDD(4) * t442 + t397 * t514 + t445 * t458;
t422 = t439 * t508;
t483 = -t437 * t509 + t422;
t356 = qJD(2) * t399 - t483;
t350 = qJDD(5) + t356;
t377 = qJD(4) * t442 + t397 * t445;
t501 = t447 * t520;
t510 = qJDD(2) * qJ(3);
t511 = qJDD(1) * t438;
t379 = t444 * t511 + t510 + (qJD(3) + t501) * qJD(2);
t490 = qJDD(1) * t558;
t417 = t439 * t490;
t555 = pkin(8) * qJDD(2);
t332 = t417 + (-t379 - t555) * t437;
t348 = t379 * t439 + t437 * t490;
t333 = t439 * t555 + t348;
t482 = t332 * t443 + t333 * t446;
t573 = t365 * t446 - t366 * t443;
t297 = qJDD(4) * pkin(9) + qJD(4) * t573 + t482;
t318 = qJD(4) * pkin(9) + t322;
t426 = pkin(3) * t439 + pkin(2);
t484 = qJD(3) - t501;
t385 = -qJD(2) * t426 + t484;
t395 = -t423 + t500;
t325 = pkin(4) * t395 - pkin(9) * t397 + t385;
t305 = t318 * t445 + t325 * t442;
t496 = t447 * t511;
t517 = qJD(2) * t444;
t497 = qJD(1) * t517;
t570 = t438 * t497 + qJDD(3);
t470 = -t496 + t570;
t368 = -qJDD(2) * t426 + t470;
t310 = pkin(4) * t356 + pkin(9) * t458 + t368;
t308 = t445 * t310;
t452 = -qJD(5) * t305 - t442 * t297 + t308;
t290 = pkin(5) * t350 + qJ(6) * t315 - qJD(6) * t377 + t452;
t375 = t397 * t442 - t512;
t300 = -qJ(6) * t375 + t305;
t386 = qJD(5) + t395;
t581 = t300 * t386 + t290;
t316 = qJD(5) * t377 - qJDD(4) * t445 - t442 * t458;
t513 = qJD(5) * t445;
t465 = -t297 * t445 - t310 * t442 + t318 * t514 - t325 * t513;
t291 = -qJ(6) * t316 - qJD(6) * t375 - t465;
t304 = -t318 * t442 + t325 * t445;
t299 = -qJ(6) * t377 + t304;
t295 = pkin(5) * t386 + t299;
t580 = -t295 * t386 + t291;
t536 = t439 * MDP(5);
t579 = -t437 * MDP(6) + t536;
t521 = t437 ^ 2 + t439 ^ 2;
t578 = MDP(7) * t521;
t436 = pkin(11) + qJ(4);
t429 = sin(t436);
t557 = cos(pkin(10));
t486 = t558 * t557;
t556 = sin(pkin(10));
t391 = t444 * t556 - t447 * t486;
t485 = t558 * t556;
t393 = t444 * t557 + t447 * t485;
t488 = g(1) * t393 + g(2) * t391;
t572 = g(3) * t537 - t488;
t575 = t572 * t429;
t574 = -t373 * t442 + t445 * t583;
t353 = pkin(4) * t401 - pkin(9) * t402 - t426;
t571 = t353 * t513 + t442 * t583 + t445 * t584;
t392 = t444 * t486 + t447 * t556;
t394 = -t444 * t485 + t447 * t557;
t487 = g(1) * t394 + g(2) * t392;
t569 = pkin(5) * t316 + qJDD(6);
t546 = (-t404 * t437 + t419) * t437;
t480 = -t381 * t439 + t546;
t567 = t447 * t480 - (-qJD(2) * pkin(2) + t484) * t444;
t566 = t377 ^ 2;
t560 = g(3) * t438;
t440 = -qJ(6) - pkin(9);
t554 = qJDD(2) * pkin(2);
t553 = qJDD(4) * pkin(4);
t552 = t315 * t442;
t550 = t375 * t395;
t549 = t375 * t397;
t548 = t377 * t386;
t547 = t377 * t397;
t545 = t386 * t442;
t544 = t395 * t442;
t543 = t402 * t442;
t542 = t402 * t445;
t430 = cos(t436);
t540 = t430 * t442;
t538 = t438 * t444;
t535 = t442 * t350;
t534 = t442 * t447;
t339 = t445 * t350;
t371 = -t407 * t443 + t408 * t446;
t364 = t445 * t371;
t448 = qJD(2) ^ 2;
t533 = t447 * t448;
t532 = qJDD(1) - g(3);
t477 = qJ(6) * t398 - qJD(6) * t402;
t531 = pkin(5) * t399 - t330 * t442 + t477 * t445 + (-t364 + (qJ(6) * t402 - t353) * t442) * qJD(5) + t574;
t498 = t402 * t513;
t530 = -qJ(6) * t498 + (-qJD(5) * t371 + t477) * t442 + t571;
t529 = -t299 + t295;
t528 = -t316 * t442 - t375 * t513;
t351 = pkin(4) * t397 + pkin(9) * t395;
t527 = t351 * t442 + t445 * t573;
t461 = t402 * t537;
t526 = -qJD(1) * t461 + qJD(3) * t402 + qJD(4) * t371;
t525 = -t386 * t544 + t339;
t524 = t353 * t442 + t364;
t494 = qJD(5) * t440;
t523 = -qJ(6) * t544 + qJD(6) * t445 + t442 * t494 - t527;
t342 = t445 * t351;
t522 = -pkin(5) * t397 - t342 + (-qJ(6) * t395 + t494) * t445 + (-qJD(6) + t573) * t442;
t515 = qJD(5) * t386;
t507 = g(3) * t538;
t506 = t438 * t534;
t505 = t445 * t537;
t499 = t438 * t517;
t493 = t438 * t557;
t492 = t438 * t556;
t489 = t386 * t445;
t389 = -t437 * t538 + t439 * t558;
t390 = t437 * t558 + t439 * t538;
t479 = t389 * t446 - t390 * t443;
t337 = t389 * t443 + t390 * t446;
t382 = t470 - t554;
t475 = -t382 + t488;
t473 = -t332 * t446 + t443 * t333 + t582;
t317 = -qJD(4) * pkin(4) - t573;
t326 = -t337 * t442 - t505;
t471 = -t337 * t445 + t506;
t468 = -t398 * t442 + t498;
t467 = -t398 * t445 - t402 * t514;
t466 = -pkin(9) * t350 + t317 * t386;
t357 = t392 * t429 + t430 * t493;
t359 = t394 * t429 - t430 * t492;
t383 = t429 * t538 - t430 * t558;
t464 = g(1) * t359 + g(2) * t357 + g(3) * t383;
t358 = t392 * t430 - t429 * t493;
t360 = t394 * t430 + t429 * t492;
t384 = t429 * t558 + t430 * t538;
t463 = g(1) * t360 + g(2) * t358 + g(3) * t384;
t298 = t473 - t553;
t456 = t572 * t430;
t455 = -t572 + t496;
t347 = -t379 * t437 + t417;
t454 = -t347 * t437 + t348 * t439 - t487;
t451 = pkin(9) * t515 + t298 - t464;
t450 = t464 - t473;
t449 = -g(1) * (-t360 * t442 + t393 * t445) - g(2) * (-t358 * t442 + t391 * t445) - g(3) * (-t384 * t442 - t505);
t428 = pkin(5) * t445 + pkin(4);
t410 = t440 * t445;
t409 = t440 * t442;
t374 = t375 ^ 2;
t345 = t445 * t353;
t320 = qJD(2) * t461 + qJD(4) * t337;
t319 = -qJD(2) * t460 + qJD(4) * t479;
t312 = -qJ(6) * t543 + t524;
t311 = pkin(5) * t375 + qJD(6) + t317;
t309 = pkin(5) * t401 - qJ(6) * t542 - t371 * t442 + t345;
t303 = qJD(5) * t471 - t319 * t442 + t445 * t499;
t302 = qJD(5) * t326 + t319 * t445 + t442 * t499;
t292 = t298 + t569;
t1 = [t532 * MDP(1) + (-t389 * t437 + t390 * t439) * qJDD(2) * MDP(7) + (t347 * t389 + t348 * t390 - g(3)) * MDP(8) + (-qJD(4) * t320 + qJDD(4) * t479) * MDP(14) + (-t319 * qJD(4) - t337 * qJDD(4)) * MDP(15) + (t303 * t386 - t316 * t479 + t320 * t375 + t326 * t350) * MDP(21) + (-t302 * t386 + t315 * t479 + t320 * t377 + t350 * t471) * MDP(22) + (-t302 * t375 - t303 * t377 + t315 * t326 + t316 * t471) * MDP(23) + (t290 * t326 - t291 * t471 - t292 * t479 + t295 * t303 + t300 * t302 + t311 * t320 - g(3)) * MDP(24) + ((-qJDD(2) * t444 - t533) * MDP(4) + t533 * t578 + (-qJD(2) * t567 - t382 * t447) * MDP(8) + (-t356 * t447 + t395 * t517) * MDP(14) + (t397 * t517 + t447 * t458) * MDP(15) + (MDP(3) + t579) * (qJDD(2) * t447 - t444 * t448)) * t438; qJDD(2) * MDP(2) + t455 * MDP(3) + (-t532 * t538 + t487) * MDP(4) + (-t507 + t454 + (qJD(2) * t484 + t510) * t521) * MDP(7) + (-t480 * qJD(3) + t475 * pkin(2) + t454 * qJ(3) + (-g(3) * (pkin(2) * t447 + qJ(3) * t444) + t567 * qJD(1)) * t438) * MDP(8) + (-t397 * t398 - t402 * t458) * MDP(9) + (-t402 * t356 + t398 * t395 - t397 * t399 + t401 * t458) * MDP(10) + (-qJD(4) * t398 + qJDD(4) * t402) * MDP(11) + (-qJD(4) * t399 - qJDD(4) * t401) * MDP(12) + (-qJD(4) * t526 + qJDD(4) * t478 - t356 * t426 + t368 * t401 + t385 * t399 - t395 * t502 - t456) * MDP(14) + (-t371 * qJDD(4) - t426 * t503 + t368 * t402 - t385 * t398 - t397 * t502 + (t426 * t500 - t584) * qJD(4) + t575) * MDP(15) + (-t315 * t542 + t377 * t467) * MDP(16) + (-(-t375 * t445 - t377 * t442) * t398 + (t552 - t316 * t445 + (t375 * t442 - t377 * t445) * qJD(5)) * t402) * MDP(17) + (-t315 * t401 + t339 * t402 + t377 * t399 + t386 * t467) * MDP(18) + (-t316 * t401 - t375 * t399 - t386 * t468 - t402 * t535) * MDP(19) + (t350 * t401 + t386 * t399) * MDP(20) + (t304 * t399 + t308 * t401 - t478 * t316 + t345 * t350 + t574 * t386 + t526 * t375 + (-t456 + (t317 * t402 - t318 * t401 - t371 * t386) * qJD(5)) * t445 + ((-qJD(5) * t353 - t330) * t386 - t371 * t350 + (-qJD(5) * t325 - t297) * t401 + t298 * t402 - t317 * t398 - t507 - t487) * t442) * MDP(21) + (-t524 * t350 + t465 * t401 - t305 * t399 + t478 * t315 + t298 * t542 - g(1) * (t393 * t540 + t394 * t445) - g(2) * (t391 * t540 + t392 * t445) - (-t430 * t534 + t444 * t445) * t560 + (t371 * t514 - t571) * t386 + t526 * t377 + t467 * t317) * MDP(22) + (t309 * t315 - t312 * t316 - (-t295 * t445 - t300 * t442) * t398 - t531 * t377 - t530 * t375 - t575 + (-t290 * t445 - t291 * t442 + (t295 * t442 - t300 * t445) * qJD(5)) * t402) * MDP(23) + (t291 * t312 + t290 * t309 + t292 * (pkin(5) * t543 - t478) + (pkin(5) * t468 + t526) * t311 + t530 * t300 + t531 * t295 + (-t444 * t560 - t487) * (pkin(5) * t442 + t559) + (-t447 * t560 + t488) * (t428 * t430 - t429 * t440 + t426)) * MDP(24) + t579 * (t438 * (-g(3) * t447 + t497) + t475 + t554); (qJD(2) * t546 - t381 * t519 - t455 + t570) * MDP(8) - t422 * MDP(14) + t503 * MDP(15) + (t525 - t549) * MDP(21) - MDP(22) * t547 + t528 * MDP(23) + (-t311 * t397 + t572) * MDP(24) - t448 * t578 + (-t536 - pkin(2) * MDP(8) + (MDP(14) * t443 + MDP(6)) * t437) * qJDD(2) + (-MDP(21) * t515 - t350 * MDP(22) + MDP(23) * t548 + MDP(24) * t580) * t442 + ((t437 * t516 + t439 * t518 + t397) * MDP(14) + (-t395 - t500) * MDP(15)) * qJD(4) + ((t315 - t550) * MDP(23) + t581 * MDP(24) - t386 ^ 2 * MDP(22)) * t445; -t395 ^ 2 * MDP(10) + ((t395 - t500) * qJD(4) + t503) * MDP(11) + t483 * MDP(12) + qJDD(4) * MDP(13) + (t450 + t582) * MDP(14) + (t385 * t395 + t463 - t482) * MDP(15) + (t377 * t489 - t552) * MDP(16) + ((-t315 - t550) * t445 - t377 * t545 + t528) * MDP(17) + (t386 * t489 + t535 - t547) * MDP(18) + (-t386 * t514 + t525 + t549) * MDP(19) + (-pkin(4) * t316 - t322 * t375 - t342 * t386 + (t386 * t573 + t466) * t442 - t451 * t445) * MDP(21) + (pkin(4) * t315 - t322 * t377 + t386 * t527 + t442 * t451 + t445 * t466) * MDP(22) + (t315 * t409 + t316 * t410 - t523 * t375 - t522 * t377 - t442 * t581 + t445 * t580 - t463) * MDP(23) + (-t291 * t410 + t290 * t409 - t292 * t428 - g(1) * (-t359 * t428 - t360 * t440) - g(2) * (-t357 * t428 - t358 * t440) - g(3) * (-t383 * t428 - t384 * t440) + (pkin(5) * t545 - t322) * t311 + t523 * t300 + t522 * t295) * MDP(24) + (MDP(10) * t397 - MDP(14) * t385 - MDP(20) * t386 - MDP(21) * t304 + MDP(22) * t305 + MDP(9) * t395) * t397; t377 * t375 * MDP(16) + (-t374 + t566) * MDP(17) + (t375 * t386 - t315) * MDP(18) + (-t316 + t548) * MDP(19) + t350 * MDP(20) + (t305 * t386 - t317 * t377 + t449 + t452) * MDP(21) + (t304 * t386 + t317 * t375 - g(1) * (-t360 * t445 - t393 * t442) - g(2) * (-t358 * t445 - t391 * t442) - g(3) * (-t384 * t445 + t506) + t465) * MDP(22) + (pkin(5) * t315 - t375 * t529) * MDP(23) + (t529 * t300 + (-t311 * t377 + t290 + t449) * pkin(5)) * MDP(24); (-t374 - t566) * MDP(23) + (t295 * t377 + t300 * t375 - t450 - t553 + t569) * MDP(24);];
tau  = t1;
