% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:24
% EndTime: 2019-03-09 06:01:32
% DurationCPUTime: 4.41s
% Computational Cost: add. (4239->396), mult. (9834->528), div. (0->0), fcn. (6536->8), ass. (0->183)
t456 = cos(qJ(4));
t455 = sin(qJ(3));
t457 = cos(qJ(3));
t539 = t456 * t457;
t473 = pkin(4) * t455 - pkin(9) * t539;
t478 = pkin(3) * t455 - pkin(8) * t457;
t414 = t478 * qJD(1);
t454 = sin(qJ(4));
t440 = sin(pkin(10)) * pkin(1) + pkin(7);
t422 = t440 * qJD(1);
t566 = qJD(2) * t457 - t455 * t422;
t484 = t456 * t414 - t454 * t566;
t557 = pkin(8) + pkin(9);
t504 = qJD(4) * t557;
t580 = qJD(1) * t473 + t456 * t504 + t484;
t519 = qJD(1) * t457;
t503 = t454 * t519;
t527 = t454 * t414 + t456 * t566;
t579 = pkin(9) * t503 - t454 * t504 - t527;
t512 = qJD(4) * t455;
t495 = qJD(1) * t512;
t509 = qJD(1) * qJD(3);
t496 = t457 * t509;
t576 = qJD(3) * qJD(4) + t496;
t367 = -t454 * t495 + t456 * t576;
t515 = qJD(3) * t456;
t520 = qJD(1) * t455;
t406 = -t454 * t520 + t515;
t517 = qJD(3) * t454;
t407 = t456 * t520 + t517;
t453 = sin(qJ(5));
t556 = cos(qJ(5));
t498 = t556 * qJD(5);
t505 = t454 * t576 + t456 * t495;
t510 = qJD(5) * t453;
t313 = -t556 * t367 - t406 * t498 + t407 * t510 + t453 * t505;
t472 = t453 * t406 + t407 * t556;
t314 = qJD(5) * t472 + t453 * t367 + t556 * t505;
t353 = -t556 * t406 + t407 * t453;
t351 = t353 ^ 2;
t438 = -qJD(4) + t519;
t430 = -qJD(5) + t438;
t516 = qJD(3) * t455;
t493 = MDP(23) * t516;
t558 = t472 ^ 2;
t578 = qJD(1) * t493 + (-t430 * t472 - t314) * MDP(22) + t353 * MDP(19) * t472 + (-t353 * t430 - t313) * MDP(21) + (-t351 + t558) * MDP(20);
t577 = t353 * qJ(6);
t545 = t453 * t454;
t471 = t456 * t556 - t545;
t562 = qJD(4) + qJD(5);
t563 = t556 * qJD(4) + t498;
t529 = -t456 * t563 + t471 * t519 + t545 * t562;
t409 = t453 * t456 + t454 * t556;
t359 = t562 * t409;
t528 = -t409 * t519 + t359;
t448 = t455 * qJD(2);
t380 = t457 * t422 + t448;
t513 = qJD(4) * t454;
t567 = -t380 + (-t503 + t513) * pkin(4);
t511 = qJD(4) * t456;
t500 = t455 * t511;
t514 = qJD(3) * t457;
t502 = t454 * t514;
t575 = t500 + t502;
t371 = -qJD(3) * pkin(3) - t566;
t350 = -pkin(4) * t406 + t371;
t372 = qJD(3) * pkin(8) + t380;
t441 = -cos(pkin(10)) * pkin(1) - pkin(2);
t400 = -pkin(3) * t457 - pkin(8) * t455 + t441;
t375 = t400 * qJD(1);
t544 = t454 * t375;
t337 = t372 * t456 + t544;
t373 = t566 * qJD(3);
t417 = t478 * qJD(3);
t399 = qJD(1) * t417;
t485 = t454 * t373 - t456 * t399;
t464 = -qJD(4) * t337 - t485;
t497 = t455 * t509;
t303 = pkin(4) * t497 - pkin(9) * t367 + t464;
t506 = t456 * t373 + t375 * t511 + t454 * t399;
t469 = -t372 * t513 + t506;
t306 = -pkin(9) * t505 + t469;
t336 = -t372 * t454 + t456 * t375;
t329 = -pkin(9) * t407 + t336;
t322 = -pkin(4) * t438 + t329;
t330 = pkin(9) * t406 + t337;
t482 = -t453 * t303 - t556 * t306 - t322 * t498 + t330 * t510;
t574 = t350 * t353 + t482;
t572 = MDP(5) * t455;
t449 = t455 ^ 2;
t571 = MDP(6) * (-t457 ^ 2 + t449);
t570 = qJ(6) * t472;
t323 = pkin(5) * t353 + qJD(6) + t350;
t569 = t323 * t472;
t568 = t457 * t505;
t386 = t456 * t400;
t541 = t455 * t456;
t546 = t440 * t454;
t344 = -pkin(9) * t541 + t386 + (-pkin(4) - t546) * t457;
t413 = t440 * t539;
t523 = t454 * t400 + t413;
t543 = t454 * t455;
t349 = -pkin(9) * t543 + t523;
t530 = t453 * t344 + t556 * t349;
t425 = t557 * t454;
t426 = t557 * t456;
t524 = -t453 * t425 + t556 * t426;
t565 = qJD(5) * t524 + t579 * t453 + t580 * t556;
t564 = -t425 * t498 - t426 * t510 - t580 * t453 + t579 * t556;
t328 = t556 * t330;
t298 = t453 * t322 + t328;
t462 = -qJD(5) * t298 + t556 * t303 - t453 * t306;
t561 = -t350 * t472 + t462;
t476 = qJD(1) * t449 - t438 * t457;
t501 = t454 * t512;
t560 = -t438 * t501 - t476 * t515;
t555 = t367 * t454;
t554 = t371 * t454;
t553 = t371 * t456;
t374 = qJD(3) * t448 + t422 * t514;
t552 = t374 * t454;
t551 = t374 * t456;
t550 = t406 * t438;
t549 = t406 * t455;
t548 = t407 * t438;
t547 = t438 * t456;
t326 = t453 * t330;
t542 = t454 * t457;
t458 = qJD(3) ^ 2;
t540 = t455 * t458;
t538 = t457 * t458;
t297 = t556 * t322 - t326;
t293 = t297 - t570;
t292 = -pkin(5) * t430 + t293;
t537 = t292 - t293;
t536 = -qJ(6) * t528 + qJD(6) * t471 + t564;
t535 = -pkin(5) * t520 + qJ(6) * t529 - t409 * qJD(6) - t565;
t479 = t556 * t514;
t331 = t359 * t455 + t453 * t502 - t456 * t479;
t384 = t471 * t455;
t534 = -t384 * t314 + t331 * t353;
t332 = t454 * t479 - t453 * t501 - t510 * t543 + (t453 * t514 + t455 * t563) * t456;
t383 = t409 * t455;
t533 = t332 * t430 - t383 * t497;
t532 = t556 * t329 - t326;
t526 = t400 * t511 + t454 * t417;
t525 = t456 * t417 + t516 * t546;
t389 = pkin(4) * t543 + t455 * t440;
t423 = qJD(1) * t441;
t360 = pkin(4) * t575 + t440 * t514;
t447 = -pkin(4) * t456 - pkin(3);
t494 = MDP(16) * t516;
t492 = t313 * t457 + t472 * t516;
t491 = -t329 * t453 - t328;
t489 = t556 * t344 - t349 * t453;
t487 = -t367 * t457 + t407 * t516;
t486 = t438 * t440 + t372;
t483 = -t556 * t425 - t426 * t453;
t480 = t438 * t500;
t477 = -t313 * t383 + t332 * t472;
t474 = 0.2e1 * qJD(3) * t423;
t470 = t314 * t457 - t353 * t516;
t317 = t473 * qJD(3) + (-t413 + (pkin(9) * t455 - t400) * t454) * qJD(4) + t525;
t319 = (-t455 * t515 - t457 * t513) * t440 - t575 * pkin(9) + t526;
t468 = t453 * t317 + t556 * t319 + t344 * t498 - t349 * t510;
t466 = t476 * t454;
t339 = pkin(4) * t505 + t374;
t465 = -t331 * t430 - t384 * t497;
t299 = t314 * pkin(5) + t339;
t461 = -qJD(5) * t530 + t556 * t317 - t453 * t319;
t446 = pkin(4) * t556 + pkin(5);
t347 = qJ(6) * t471 + t524;
t346 = -qJ(6) * t409 + t483;
t309 = -qJ(6) * t383 + t530;
t307 = -pkin(5) * t457 - qJ(6) * t384 + t489;
t296 = t532 - t570;
t295 = t491 + t577;
t294 = t298 - t577;
t291 = -qJ(6) * t332 - qJD(6) * t383 + t468;
t290 = pkin(5) * t516 + t331 * qJ(6) - t384 * qJD(6) + t461;
t289 = -qJ(6) * t314 - qJD(6) * t353 - t482;
t288 = pkin(5) * t497 + t313 * qJ(6) - qJD(6) * t472 + t462;
t1 = [0.2e1 * t496 * t572 - 0.2e1 * t509 * t571 + MDP(7) * t538 - MDP(8) * t540 + (-t440 * t538 + t455 * t474) * MDP(10) + (t440 * t540 + t457 * t474) * MDP(11) + (t367 * t541 + (t456 * t514 - t501) * t407) * MDP(12) + ((t406 * t456 - t407 * t454) * t514 + (-t456 * t505 - t555 + (-t406 * t454 - t407 * t456) * qJD(4)) * t455) * MDP(13) + (t487 - t560) * MDP(14) + (t480 + t568 + (-t466 + t549) * qJD(3)) * MDP(15) + (-t438 - t519) * t494 + (-(-t400 * t513 + t525) * t438 + ((-t406 * t440 + t554) * qJD(3) + (t456 * t486 + t544) * qJD(4) + t485) * t457 + (t440 * t505 + t552 + t371 * t511 + ((-t440 * t542 + t386) * qJD(1) + t336) * qJD(3)) * t455) * MDP(17) + (t526 * t438 + (-t486 * t513 + (t407 * t440 + t553) * qJD(3) + t506) * t457 + (-t371 * t513 + t440 * t367 + t551 + (-qJD(1) * t523 - t440 * t547 - t337) * qJD(3)) * t455) * MDP(18) + (-t313 * t384 - t331 * t472) * MDP(19) + (-t477 + t534) * MDP(20) + (-t465 + t492) * MDP(21) + (t470 + t533) * MDP(22) + (-t430 - t519) * t493 + (t297 * t516 + t389 * t314 + t350 * t332 + t339 * t383 + t360 * t353 - t430 * t461 - t457 * t462 + t489 * t497) * MDP(24) + (t468 * t430 - t482 * t457 + t360 * t472 - t389 * t313 + t339 * t384 - t350 * t331 + (-qJD(1) * t530 - t298) * t516) * MDP(25) + (-t288 * t384 - t289 * t383 - t290 * t472 - t291 * t353 + t292 * t331 - t294 * t332 + t307 * t313 - t309 * t314) * MDP(26) + (t289 * t309 + t294 * t291 + t288 * t307 + t292 * t290 + t299 * (pkin(5) * t383 + t389) + t323 * (pkin(5) * t332 + t360)) * MDP(27); -MDP(10) * t540 - MDP(11) * t538 + (t480 - t568 + (-t466 - t549) * qJD(3)) * MDP(17) + (t487 + t560) * MDP(18) + (-t470 + t533) * MDP(24) + (t465 + t492) * MDP(25) + (t477 + t534) * MDP(26) + (-t288 * t383 + t289 * t384 - t292 * t332 - t294 * t331 - t299 * t457 + t323 * t516) * MDP(27); (qJD(3) * t380 - t374) * MDP(10) - t423 * t519 * MDP(11) + (-t407 * t547 + t555) * MDP(12) + ((t367 - t550) * t456 + (-t505 + t548) * t454) * MDP(13) + (-t438 * t511 + (t438 * t539 + (-t407 + t517) * t455) * qJD(1)) * MDP(14) + (t438 * t513 + (-t438 * t542 + (-t406 + t515) * t455) * qJD(1)) * MDP(15) + (-pkin(3) * t505 - t551 + t484 * t438 + t380 * t406 + (pkin(8) * t547 + t554) * qJD(4) + (-t336 * t455 + (-pkin(8) * t516 - t371 * t457) * t454) * qJD(1)) * MDP(17) + (-pkin(3) * t367 + t552 - t527 * t438 - t380 * t407 + (-pkin(8) * t438 * t454 + t553) * qJD(4) + (-t371 * t539 + (-pkin(8) * t515 + t337) * t455) * qJD(1)) * MDP(18) + (-t313 * t409 - t472 * t529) * MDP(19) + (-t313 * t471 - t314 * t409 + t353 * t529 - t472 * t528) * MDP(20) + (t447 * t314 - t339 * t471 + t350 * t528 + t353 * t567 + t483 * t497) * MDP(24) + (-t447 * t313 + t339 * t409 - t529 * t350 + t567 * t472) * MDP(25) + (-t288 * t409 + t289 * t471 + t292 * t529 - t294 * t528 + t313 * t346 - t314 * t347 - t353 * t536 - t472 * t535) * MDP(26) + (t289 * t347 + t288 * t346 + t299 * (-pkin(5) * t471 + t447) + (pkin(5) * t528 + t567) * t323 + t536 * t294 + t535 * t292) * MDP(27) + (-t457 * t572 + t571) * qJD(1) ^ 2 + (-t423 * MDP(10) + t438 * MDP(16) + (qJD(3) * t409 - t472) * MDP(21) + (qJD(3) * t471 + t353) * MDP(22) - t297 * MDP(24) + (-qJD(3) * t524 + t298) * MDP(25)) * t520 + (t529 * MDP(21) + t528 * MDP(22) + MDP(23) * t520 + MDP(24) * t565 + MDP(25) * t564) * t430; -t407 * t406 * MDP(12) + (-t406 ^ 2 + t407 ^ 2) * MDP(13) + (t367 + t550) * MDP(14) + (-t505 - t548) * MDP(15) + qJD(1) * t494 + (-t337 * t438 - t371 * t407 + t464) * MDP(17) + (-t336 * t438 - t371 * t406 - t469) * MDP(18) + (t491 * t430 + (-t407 * t353 + t430 * t510 + t497 * t556) * pkin(4) + t561) * MDP(24) + (-t532 * t430 + (-t407 * t472 + t430 * t498 - t453 * t497) * pkin(4) + t574) * MDP(25) + (-t292 * t353 + t294 * t472 + t295 * t472 + t296 * t353 + t446 * t313 + (-t314 * t453 + (-t353 * t556 + t453 * t472) * qJD(5)) * pkin(4)) * MDP(26) + (-pkin(5) * t569 + t288 * t446 - t292 * t295 - t294 * t296 + (t289 * t453 - t323 * t407 + (-t292 * t453 + t294 * t556) * qJD(5)) * pkin(4)) * MDP(27) + t578; (-t298 * t430 + t561) * MDP(24) + (-t297 * t430 + t574) * MDP(25) + (pkin(5) * t313 - t353 * t537) * MDP(26) + (t537 * t294 + (t288 - t569) * pkin(5)) * MDP(27) + t578; (-t351 - t558) * MDP(26) + (t292 * t472 + t294 * t353 + t299) * MDP(27);];
tauc  = t1;
