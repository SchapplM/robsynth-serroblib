% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:31
% EndTime: 2019-03-09 03:35:38
% DurationCPUTime: 5.51s
% Computational Cost: add. (5203->409), mult. (11577->529), div. (0->0), fcn. (8855->16), ass. (0->198)
t477 = sin(pkin(11));
t479 = cos(pkin(11));
t484 = sin(qJ(3));
t488 = cos(qJ(3));
t436 = -t477 * t484 + t479 * t488;
t427 = t436 * qJD(1);
t487 = cos(qJ(5));
t417 = t487 * t427;
t437 = t477 * t488 + t479 * t484;
t429 = t437 * qJD(1);
t483 = sin(qJ(5));
t381 = -t429 * t483 + t417;
t486 = cos(qJ(6));
t545 = qJD(6) * t486;
t603 = -t381 * t486 + t545;
t470 = qJ(3) + pkin(11) + qJ(5);
t458 = sin(t470);
t474 = qJ(1) + pkin(10);
t467 = sin(t474);
t468 = cos(t474);
t518 = g(1) * t468 + g(2) * t467;
t602 = t518 * t458;
t511 = t427 * t483 + t487 * t429;
t428 = t437 * qJD(3);
t388 = -qJD(1) * t428 + qJDD(1) * t436;
t542 = qJD(1) * qJD(3);
t532 = t488 * t542;
t533 = t484 * t542;
t389 = qJDD(1) * t437 - t477 * t533 + t479 * t532;
t547 = qJD(5) * t483;
t339 = qJD(5) * t417 + t483 * t388 + t487 * t389 - t429 * t547;
t472 = qJDD(3) + qJDD(5);
t473 = qJD(3) + qJD(5);
t482 = sin(qJ(6));
t534 = t486 * t339 + t482 * t472 + t473 * t545;
t546 = qJD(6) * t482;
t317 = -t511 * t546 + t534;
t316 = t317 * t486;
t375 = t473 * t482 + t486 * t511;
t455 = t486 * t472;
t318 = qJD(6) * t375 + t339 * t482 - t455;
t373 = -t486 * t473 + t482 * t511;
t601 = -t482 * t318 - t603 * t373 + t316;
t315 = t317 * t482;
t340 = qJD(5) * t511 - t487 * t388 + t389 * t483;
t337 = qJDD(6) + t340;
t334 = t482 * t337;
t564 = t381 * t473;
t566 = t511 * t473;
t568 = t375 * t511;
t594 = qJD(6) - t381;
t600 = t472 * MDP(18) + (-t340 + t566) * MDP(17) - t381 ^ 2 * MDP(15) + (-MDP(14) * t381 + t511 * MDP(15) - MDP(25) * t594) * t511 + (t339 - t564) * MDP(16) + (t603 * t375 + t315) * MDP(21) + (t603 * t594 + t334 - t568) * MDP(23);
t478 = sin(pkin(10));
t460 = pkin(1) * t478 + pkin(7);
t557 = qJ(4) + t460;
t524 = t557 * qJD(1);
t405 = qJD(2) * t484 + t488 * t524;
t395 = t477 * t405;
t404 = t488 * qJD(2) - t484 * t524;
t572 = qJD(3) * pkin(3);
t399 = t404 + t572;
t362 = t479 * t399 - t395;
t577 = pkin(8) * t429;
t349 = qJD(3) * pkin(4) + t362 - t577;
t558 = t479 * t405;
t363 = t477 * t399 + t558;
t578 = pkin(8) * t427;
t350 = t363 + t578;
t327 = t349 * t487 - t350 * t483;
t323 = -pkin(5) * t473 - t327;
t599 = t323 * t381;
t469 = t488 * qJDD(2);
t445 = t460 * qJDD(1);
t497 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t445;
t508 = t524 * qJD(3);
t360 = qJDD(3) * pkin(3) - t484 * t497 - t488 * t508 + t469;
t364 = (qJDD(2) - t508) * t484 + t497 * t488;
t331 = t479 * t360 - t364 * t477;
t325 = qJDD(3) * pkin(4) - pkin(8) * t389 + t331;
t332 = t477 * t360 + t479 * t364;
t326 = pkin(8) * t388 + t332;
t328 = t349 * t483 + t350 * t487;
t582 = qJD(5) * t328 - t487 * t325 + t483 * t326;
t306 = -pkin(5) * t472 + t582;
t459 = cos(t470);
t574 = g(3) * t459;
t597 = t306 + t574;
t464 = pkin(3) * t488 + pkin(2);
t480 = cos(pkin(10));
t580 = pkin(1) * t480;
t509 = -t464 - t580;
t425 = qJD(1) * t509 + qJD(4);
t390 = -pkin(4) * t427 + t425;
t452 = g(3) * t458;
t583 = (qJD(5) * t349 + t326) * t487 + t483 * t325 - t350 * t547;
t593 = -t390 * t381 + t459 * t518 + t452 - t583;
t590 = pkin(5) * t511;
t569 = t373 * t511;
t589 = t594 * (pkin(9) * t594 + t590);
t461 = pkin(3) * t479 + pkin(4);
t579 = pkin(3) * t477;
t550 = t483 * t461 + t487 * t579;
t588 = pkin(3) * t533 + t509 * qJDD(1) + qJDD(4);
t324 = pkin(9) * t473 + t328;
t341 = -pkin(5) * t381 - pkin(9) * t511 + t390;
t309 = -t324 * t482 + t341 * t486;
t587 = -t309 * t511 + t323 * t546 + t486 * t602;
t310 = t324 * t486 + t341 * t482;
t586 = t310 * t511 + t323 * t545 + t482 * t597;
t585 = -t390 * t511 - t574 - t582 + t602;
t526 = qJD(3) * t557;
t412 = qJD(4) * t488 - t484 * t526;
t413 = -qJD(4) * t484 - t488 * t526;
t369 = -t412 * t477 + t479 * t413;
t431 = t436 * qJD(3);
t355 = -pkin(8) * t431 + t369;
t370 = t479 * t412 + t477 * t413;
t356 = -pkin(8) * t428 + t370;
t434 = t557 * t484;
t435 = t557 * t488;
t386 = -t479 * t434 - t435 * t477;
t371 = -pkin(8) * t437 + t386;
t387 = -t477 * t434 + t479 * t435;
t372 = pkin(8) * t436 + t387;
t512 = t371 * t487 - t372 * t483;
t311 = qJD(5) * t512 + t355 * t483 + t356 * t487;
t343 = t371 * t483 + t372 * t487;
t392 = t436 * t483 + t437 * t487;
t403 = -pkin(4) * t436 + t509;
t510 = t487 * t436 - t437 * t483;
t345 = -pkin(5) * t510 - pkin(9) * t392 + t403;
t351 = qJD(5) * t510 - t428 * t483 + t431 * t487;
t305 = pkin(9) * t472 + t583;
t522 = qJD(6) * t341 + t305;
t581 = t306 * t392 + t323 * t351 - t343 * t337 - (qJD(6) * t345 + t311) * t594 + t510 * t522;
t573 = g(3) * t488;
t571 = t323 * t392;
t570 = t345 * t337;
t567 = t375 * t482;
t562 = t467 * t482;
t561 = t467 * t486;
t560 = t468 * t482;
t559 = t468 * t486;
t335 = t486 * t337;
t556 = qJDD(2) - g(3);
t352 = qJD(5) * t392 + t487 * t428 + t431 * t483;
t555 = -t317 * t510 + t375 * t352;
t366 = t479 * t404 - t395;
t365 = -t404 * t477 - t558;
t353 = t365 - t578;
t354 = t366 - t577;
t506 = t461 * t487 - t483 * t579;
t552 = -t506 * qJD(5) + t353 * t483 + t354 * t487;
t551 = qJD(5) * t550 + t353 * t487 - t354 * t483;
t475 = t484 ^ 2;
t549 = -t488 ^ 2 + t475;
t462 = -pkin(2) - t580;
t448 = qJD(1) * t462;
t541 = qJDD(1) * t488;
t466 = t484 * t572;
t539 = t392 * t334;
t538 = t392 * t335;
t407 = pkin(4) * t428 + t466;
t406 = t484 * qJD(1) * pkin(3) + pkin(4) * t429;
t525 = t594 * t482;
t361 = -pkin(4) * t388 + t588;
t308 = pkin(5) * t340 - pkin(9) * t339 + t361;
t521 = qJD(6) * t324 - t308;
t424 = pkin(9) + t550;
t519 = -pkin(9) * t381 + qJD(6) * t424 + t406 + t590;
t517 = g(1) * t467 - g(2) * t468;
t485 = sin(qJ(1));
t489 = cos(qJ(1));
t516 = g(1) * t485 - g(2) * t489;
t515 = t318 * t510 - t352 * t373;
t514 = -t337 * t424 - t599;
t513 = t351 * t473 + t392 * t472;
t507 = t335 - (-t381 * t482 + t546) * t594;
t504 = -t351 * t482 - t392 * t545;
t503 = -t351 * t486 + t392 * t546;
t501 = -pkin(9) * t337 + t327 * t594 - t599;
t499 = -qJD(1) * t448 - t445 + t518;
t498 = 0.2e1 * qJD(3) * t448 - qJDD(3) * t460;
t490 = qJD(3) ^ 2;
t495 = -0.2e1 * qJDD(1) * t462 - t460 * t490 + t517;
t481 = -qJ(4) - pkin(7);
t443 = qJDD(3) * t488 - t484 * t490;
t442 = qJDD(3) * t484 + t488 * t490;
t423 = -pkin(5) - t506;
t411 = t459 * t559 + t562;
t410 = -t459 * t560 + t561;
t409 = -t459 * t561 + t560;
t408 = t459 * t562 + t559;
t338 = -t352 * t473 + t472 * t510;
t319 = pkin(5) * t352 - pkin(9) * t351 + t407;
t312 = qJD(5) * t343 - t355 * t487 + t356 * t483;
t307 = t486 * t308;
t1 = [(t484 * t498 + t488 * t495) * MDP(10) + (-t484 * t495 + t488 * t498) * MDP(11) + (t504 * t594 + t515 - t539) * MDP(24) + (-t503 * t594 + t538 + t555) * MDP(23) + (-t337 * t510 + t352 * t594) * MDP(25) + (-g(1) * t409 - g(2) * t411 - t307 * t510 + t309 * t352 + t312 * t373 - t512 * t318 + (t319 * t594 + t570 + (t324 * t510 - t343 * t594 + t571) * qJD(6)) * t486 + t581 * t482) * MDP(26) + (-g(1) * t408 - g(2) * t410 - t310 * t352 + t312 * t375 - t512 * t317 + (-(-qJD(6) * t343 + t319) * t594 - t570 - t521 * t510 - qJD(6) * t571) * t482 + t581 * t486) * MDP(27) + (-t311 * t473 + t339 * t403 - t343 * t472 + t351 * t390 + t361 * t392 + t407 * t511 - t458 * t517) * MDP(20) + (t339 * t392 + t351 * t511) * MDP(14) + (-t331 * t437 + t332 * t436 - t362 * t431 - t363 * t428 - t369 * t429 + t370 * t427 - t386 * t389 + t387 * t388 - t518) * MDP(12) + (g(1) * t489 + g(2) * t485) * MDP(3) + (t516 + (t478 ^ 2 + t480 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t332 * t387 + t363 * t370 + t331 * t386 + t362 * t369 + t425 * t466 - g(1) * (-pkin(1) * t485 - t464 * t467 - t468 * t481) - g(2) * (pkin(1) * t489 + t464 * t468 - t467 * t481) + t588 * t509) * MDP(13) + t513 * MDP(16) + t516 * MDP(2) + (qJDD(1) * t475 + 0.2e1 * t484 * t532) * MDP(5) + qJDD(1) * MDP(1) + 0.2e1 * (t484 * t541 - t542 * t549) * MDP(6) + (t316 * t392 - t375 * t503) * MDP(21) + ((-t373 * t486 - t567) * t351 + (-t315 - t318 * t486 + (t373 * t482 - t375 * t486) * qJD(6)) * t392) * MDP(22) + t338 * MDP(17) + t442 * MDP(7) + t443 * MDP(8) + (t339 * t510 - t340 * t392 + t351 * t381 - t352 * t511) * MDP(15) + (-t312 * t473 + t340 * t403 + t352 * t390 - t361 * t510 - t381 * t407 + t459 * t517 + t472 * t512) * MDP(19); t556 * MDP(4) + t443 * MDP(10) - t442 * MDP(11) + (t388 * t437 - t389 * t436 + t427 * t431 + t428 * t429) * MDP(12) + (t331 * t436 + t332 * t437 - t362 * t428 + t363 * t431 - g(3)) * MDP(13) + t338 * MDP(19) - t513 * MDP(20) + (-t515 - t539) * MDP(26) + (-t538 + t555) * MDP(27) + (MDP(26) * t504 + MDP(27) * t503) * t594; (-MDP(5) * t484 * t488 + MDP(6) * t549) * qJD(1) ^ 2 + (-t406 * t511 - t472 * t550 + t473 * t552 + t593) * MDP(20) + ((t363 + t365) * t429 + (t362 - t366) * t427 + (t388 * t477 - t389 * t479) * pkin(3)) * MDP(12) + t484 * qJDD(1) * MDP(7) + MDP(8) * t541 + qJDD(3) * MDP(9) + (-t484 * t556 + t488 * t499) * MDP(11) + (-t567 * t594 + t601) * MDP(22) + (t507 + t569) * MDP(24) + (t484 * t499 + t469 - t573) * MDP(10) + (-t362 * t365 - t363 * t366 + (-t573 + t331 * t479 + t332 * t477 + (-qJD(1) * t425 + t518) * t484) * pkin(3)) * MDP(13) + (t423 * t318 - t597 * t486 + t514 * t482 + t551 * t373 + (t482 * t552 - t486 * t519) * t594 + t587) * MDP(26) + (t423 * t317 + t514 * t486 - t482 * t602 + t551 * t375 + (t482 * t519 + t486 * t552) * t594 + t586) * MDP(27) + (t381 * t406 + t472 * t506 - t473 * t551 + t585) * MDP(19) + t600; (-t427 ^ 2 - t429 ^ 2) * MDP(12) + (t340 + t566) * MDP(19) + (t339 + t564) * MDP(20) + (t507 - t569) * MDP(26) + (-t486 * t594 ^ 2 - t334 - t568) * MDP(27) + (t362 * t429 - t363 * t427 - t517 + t588) * MDP(13); (t328 * t473 + t585) * MDP(19) + (t327 * t473 + t593) * MDP(20) + (-t375 * t525 + t601) * MDP(22) + (-t525 * t594 + t335 + t569) * MDP(24) + (-pkin(5) * t318 - t328 * t373 + t501 * t482 + (-t597 - t589) * t486 + t587) * MDP(26) + (-pkin(5) * t317 - t328 * t375 + t501 * t486 + (-t602 + t589) * t482 + t586) * MDP(27) + t600; t375 * t373 * MDP(21) + (-t373 ^ 2 + t375 ^ 2) * MDP(22) + (t373 * t594 + t534) * MDP(23) + (t375 * t594 + t455) * MDP(24) + t337 * MDP(25) + (-g(1) * t410 + g(2) * t408 + t310 * t594 - t323 * t375 + t307) * MDP(26) + (g(1) * t411 - g(2) * t409 + t309 * t594 + t323 * t373) * MDP(27) + ((-t305 + t452) * MDP(27) + (-MDP(24) * t511 - MDP(26) * t324 - MDP(27) * t341) * qJD(6)) * t486 + (-qJD(6) * t511 * MDP(23) + (-qJD(6) * t473 - t339) * MDP(24) + (-t522 + t452) * MDP(26) + t521 * MDP(27)) * t482;];
tau  = t1;
