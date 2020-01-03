% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:50
% EndTime: 2019-12-31 21:24:59
% DurationCPUTime: 5.80s
% Computational Cost: add. (3642->397), mult. (9240->565), div. (0->0), fcn. (6532->8), ass. (0->191)
t473 = cos(qJ(5));
t471 = sin(qJ(3));
t472 = sin(qJ(2));
t529 = qJD(1) * t472;
t512 = t471 * t529;
t474 = cos(qJ(3));
t519 = t474 * qJD(2);
t435 = t512 - t519;
t527 = qJD(2) * t471;
t437 = t474 * t529 + t527;
t468 = sin(pkin(9));
t469 = cos(pkin(9));
t489 = -t435 * t469 - t437 * t468;
t547 = t473 * t489;
t384 = t435 * t468 - t437 * t469;
t470 = sin(qJ(5));
t561 = t384 * t470;
t336 = t547 + t561;
t475 = cos(qJ(2));
t528 = qJD(1) * t475;
t456 = -qJD(3) + t528;
t448 = -qJD(5) + t456;
t562 = t336 * t448;
t518 = qJD(1) * qJD(2);
t506 = t475 * t518;
t523 = qJD(3) * t472;
t509 = t471 * t523;
t517 = qJD(2) * qJD(3);
t399 = -qJD(1) * t509 + (t506 + t517) * t474;
t462 = pkin(6) * t528;
t446 = qJD(2) * pkin(7) + t462;
t441 = -pkin(2) * t475 - pkin(7) * t472 - pkin(1);
t425 = t441 * qJD(1);
t553 = t471 * t425;
t393 = t446 * t474 + t553;
t495 = pkin(2) * t472 - pkin(7) * t475;
t439 = t495 * qJD(2);
t426 = qJD(1) * t439;
t507 = t472 * t518;
t496 = pkin(6) * t507;
t538 = -t474 * t426 - t471 * t496;
t479 = -t393 * qJD(3) - t538;
t313 = pkin(3) * t507 - qJ(4) * t399 - qJD(4) * t437 + t479;
t522 = qJD(3) * t474;
t508 = t472 * t522;
t525 = qJD(2) * t475;
t511 = t471 * t525;
t574 = t508 + t511;
t400 = qJD(1) * t574 + t471 * t517;
t524 = qJD(3) * t471;
t486 = t425 * t522 + t471 * t426 - t446 * t524;
t478 = -t474 * t496 + t486;
t318 = -qJ(4) * t400 - qJD(4) * t435 + t478;
t296 = t469 * t313 - t318 * t468;
t351 = t399 * t469 - t400 * t468;
t294 = pkin(4) * t507 - pkin(8) * t351 + t296;
t297 = t468 * t313 + t469 * t318;
t350 = -t399 * t468 - t400 * t469;
t295 = pkin(8) * t350 + t297;
t392 = t474 * t425 - t446 * t471;
t360 = -qJ(4) * t437 + t392;
t352 = -pkin(3) * t456 + t360;
t361 = -qJ(4) * t435 + t393;
t554 = t469 * t361;
t315 = t468 * t352 + t554;
t571 = pkin(8) * t489;
t306 = t315 + t571;
t520 = qJD(5) * t470;
t305 = t306 * t520;
t445 = -qJD(2) * pkin(2) + pkin(6) * t529;
t397 = pkin(3) * t435 + qJD(4) + t445;
t342 = -pkin(4) * t489 + t397;
t583 = -t470 * t294 - t473 * t295 - t342 * t336 + t305;
t504 = MDP(24) * t529;
t573 = -t473 * t384 + t470 * t489;
t582 = qJD(2) * t504 + (-t336 ^ 2 + t573 ^ 2) * MDP(21) - t336 * MDP(20) * t573;
t563 = t573 * t448;
t546 = t474 * t475;
t485 = pkin(3) * t472 - qJ(4) * t546;
t564 = -qJ(4) - pkin(7);
t502 = qJD(3) * t564;
t438 = t495 * qJD(1);
t535 = pkin(6) * t512 + t474 * t438;
t580 = -qJD(1) * t485 - qJD(4) * t471 + t474 * t502 - t535;
t421 = t471 * t438;
t521 = qJD(4) * t474;
t550 = t472 * t474;
t551 = t471 * t475;
t579 = t421 + (-pkin(6) * t550 - qJ(4) * t551) * qJD(1) - t471 * t502 - t521;
t500 = t473 * t294 - t470 * t295;
t578 = -t342 * t573 + t500;
t577 = pkin(8) * t384;
t429 = t468 * t474 + t469 * t471;
t481 = t429 * t475;
t576 = qJD(1) * t481 - t429 * qJD(3);
t488 = t468 * t471 - t469 * t474;
t575 = t456 * t488;
t572 = -0.2e1 * t518;
t570 = MDP(4) * t472;
t466 = t472 ^ 2;
t568 = (-t475 ^ 2 + t466) * MDP(5);
t541 = t579 * t468 + t580 * t469;
t540 = t580 * t468 - t579 * t469;
t567 = -t462 + (-t471 * t528 + t524) * pkin(3);
t499 = -t473 * t350 + t351 * t470;
t301 = qJD(5) * t573 + t499;
t566 = pkin(3) * t468;
t565 = pkin(6) * t471;
t560 = t399 * t471;
t559 = t435 * t456;
t558 = t437 * t456;
t557 = t445 * t471;
t556 = t445 * t474;
t555 = t456 * t474;
t356 = t468 * t361;
t552 = t471 * t472;
t476 = qJD(2) ^ 2;
t549 = t472 * t476;
t314 = t469 * t352 - t356;
t304 = -pkin(4) * t456 + t314 + t577;
t548 = t473 * t304;
t545 = t475 * t476;
t477 = qJD(1) ^ 2;
t544 = t475 * t477;
t490 = -t429 * t470 - t473 * t488;
t543 = qJD(5) * t490 + t576 * t470 + t575 * t473;
t382 = t429 * t473 - t470 * t488;
t542 = qJD(5) * t382 + t575 * t470 - t576 * t473;
t458 = pkin(6) * t546;
t526 = qJD(2) * t472;
t536 = t474 * t439 + t526 * t565;
t329 = -t472 * t521 + t485 * qJD(2) + (-t458 + (qJ(4) * t472 - t441) * t471) * qJD(3) + t536;
t537 = t471 * t439 + t441 * t522;
t338 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t550 + (-qJD(4) * t472 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t475) * t471 + t537;
t303 = t468 * t329 + t469 * t338;
t321 = t469 * t360 - t356;
t431 = t474 * t441;
t389 = -qJ(4) * t550 + t431 + (-pkin(3) - t565) * t475;
t534 = t471 * t441 + t458;
t394 = -qJ(4) * t552 + t534;
t340 = t468 * t389 + t469 * t394;
t539 = -t576 * pkin(4) + t567;
t443 = t564 * t471;
t444 = t564 * t474;
t396 = t468 * t443 - t469 * t444;
t532 = pkin(3) * t552 + t472 * pkin(6);
t515 = qJD(5) * t547 + t470 * t350 + t473 * t351;
t514 = t574 * pkin(3) + pkin(6) * t525;
t513 = -pkin(3) * t474 - pkin(2);
t510 = t475 * t519;
t503 = MDP(15) * t526;
t383 = pkin(3) * t400 + pkin(6) * t506;
t501 = pkin(1) * t572;
t302 = t469 * t329 - t338 * t468;
t320 = -t360 * t468 - t554;
t339 = t469 * t389 - t394 * t468;
t395 = t469 * t443 + t444 * t468;
t498 = t435 + t519;
t497 = -t437 + t527;
t372 = -pkin(8) * t488 + t396;
t494 = pkin(4) * t529 + t575 * pkin(8) + qJD(5) * t372 - t541;
t371 = -pkin(8) * t429 + t395;
t493 = t576 * pkin(8) + qJD(5) * t371 + t540;
t292 = t470 * t304 + t473 * t306;
t412 = t488 * t472;
t322 = -pkin(4) * t475 + pkin(8) * t412 + t339;
t411 = t429 * t472;
t324 = -pkin(8) * t411 + t340;
t492 = t322 * t470 + t324 * t473;
t491 = -t473 * t411 + t412 * t470;
t364 = -t411 * t470 - t412 * t473;
t487 = qJD(1) * t466 - t456 * t475;
t459 = pkin(3) * t469 + pkin(4);
t484 = t459 * t470 + t473 * t566;
t483 = t459 * t473 - t470 * t566;
t300 = t384 * t520 + t515;
t405 = pkin(4) * t488 + t513;
t390 = pkin(4) * t411 + t532;
t366 = t429 * t523 + t468 * t511 - t469 * t510;
t365 = -qJD(2) * t481 + t488 * t523;
t353 = pkin(3) * t437 - pkin(4) * t384;
t341 = -pkin(4) * t365 + t514;
t323 = -pkin(4) * t350 + t383;
t310 = t321 + t577;
t309 = t320 - t571;
t308 = qJD(5) * t364 - t473 * t365 - t366 * t470;
t307 = qJD(5) * t491 + t365 * t470 - t366 * t473;
t299 = pkin(8) * t365 + t303;
t298 = pkin(4) * t526 + pkin(8) * t366 + t302;
t291 = -t306 * t470 + t548;
t1 = [0.2e1 * t506 * t570 + t568 * t572 + MDP(6) * t545 - MDP(7) * t549 + (-pkin(6) * t545 + t472 * t501) * MDP(9) + (pkin(6) * t549 + t475 * t501) * MDP(10) + (t399 * t550 + (-t509 + t510) * t437) * MDP(11) + ((-t435 * t474 - t437 * t471) * t525 + (-t560 - t400 * t474 + (t435 * t471 - t437 * t474) * qJD(3)) * t472) * MDP(12) + (t456 * t509 - t399 * t475 + (t437 * t472 + t474 * t487) * qJD(2)) * MDP(13) + (t456 * t508 + t400 * t475 + (-t435 * t472 - t471 * t487) * qJD(2)) * MDP(14) + (-t456 - t528) * t503 + (-(-t441 * t524 + t536) * t456 + (t445 * t522 + pkin(6) * t400 + (t431 * qJD(1) + t392) * qJD(2)) * t472 + ((pkin(6) * t435 + t557) * qJD(2) + (t553 + (pkin(6) * t456 + t446) * t474) * qJD(3) + t538) * t475) * MDP(16) + ((-pkin(6) * t475 * t524 + t537) * t456 + t486 * t475 + (pkin(6) * t399 - t445 * t524) * t472 + ((pkin(6) * t437 + t556) * t475 + (-pkin(6) * t555 - t534 * qJD(1) - t393) * t472) * qJD(2)) * MDP(17) + (t296 * t412 - t297 * t411 + t302 * t384 + t303 * t489 + t314 * t366 + t315 * t365 - t339 * t351 + t340 * t350) * MDP(18) + (t296 * t339 + t297 * t340 + t314 * t302 + t315 * t303 + t383 * t532 + t397 * t514) * MDP(19) + (t300 * t364 + t307 * t573) * MDP(20) + (t300 * t491 - t301 * t364 + t307 * t336 - t308 * t573) * MDP(21) + (-t300 * t475 - t307 * t448 + (qJD(1) * t364 + t573) * t526) * MDP(22) + (t301 * t475 + t308 * t448 + (qJD(1) * t491 + t336) * t526) * MDP(23) + (-t448 - t528) * MDP(24) * t526 + (-(t298 * t473 - t299 * t470) * t448 - t500 * t475 - t341 * t336 + t390 * t301 - t323 * t491 + t342 * t308 + (t292 * t475 + t448 * t492) * qJD(5) + ((t322 * t473 - t324 * t470) * qJD(1) + t291) * t526) * MDP(25) + (t390 * t300 - t305 * t475 + t342 * t307 + t323 * t364 + t341 * t573 + ((-qJD(5) * t324 + t298) * t448 + t294 * t475) * t470 + ((qJD(5) * t322 + t299) * t448 + (qJD(5) * t304 + t295) * t475) * t473 + (-qJD(1) * t492 - t292) * t526) * MDP(26); -t544 * t570 + t477 * t568 + (-t437 * t555 + t560) * MDP(11) + ((t399 + t559) * t474 + (-t400 + t558) * t471) * MDP(12) + (-t456 * t522 + (t456 * t546 + t472 * t497) * qJD(1)) * MDP(13) + (t456 * t524 + (-t456 * t551 + t472 * t498) * qJD(1)) * MDP(14) + t456 * MDP(15) * t529 + (-pkin(2) * t400 + t535 * t456 + (pkin(7) * t555 + t557) * qJD(3) + ((-pkin(7) * t527 - t392) * t472 + (-pkin(6) * t498 - t557) * t475) * qJD(1)) * MDP(16) + (-pkin(2) * t399 - t421 * t456 + (-pkin(7) * t456 * t471 + t556) * qJD(3) + (-t445 * t546 + (-pkin(7) * t519 + t393) * t472 + (t456 * t550 + t475 * t497) * pkin(6)) * qJD(1)) * MDP(17) + (-t296 * t429 - t297 * t488 - t575 * t314 + t576 * t315 + t350 * t396 - t351 * t395 + t541 * t384 + t540 * t489) * MDP(18) + (t296 * t395 + t297 * t396 + t541 * t314 + t540 * t315 + t383 * t513 + t567 * t397) * MDP(19) + (t300 * t382 + t543 * t573) * MDP(20) + (t300 * t490 - t301 * t382 + t336 * t543 - t542 * t573) * MDP(21) + (-t543 * t448 + (qJD(2) * t382 - t573) * t529) * MDP(22) + (t542 * t448 + (qJD(2) * t490 - t336) * t529) * MDP(23) + t448 * t504 + (t405 * t301 - t323 * t490 + (t470 * t493 + t473 * t494) * t448 + t542 * t342 - t539 * t336 + ((t371 * t473 - t372 * t470) * qJD(2) - t291) * t529) * MDP(25) + (t405 * t300 + t323 * t382 + (-t470 * t494 + t473 * t493) * t448 + t543 * t342 + t539 * t573 + (-(t371 * t470 + t372 * t473) * qJD(2) + t292) * t529) * MDP(26) + (t477 * t472 * MDP(9) + MDP(10) * t544) * pkin(1); t437 * t435 * MDP(11) + (-t435 ^ 2 + t437 ^ 2) * MDP(12) + (t399 - t559) * MDP(13) + (-t400 - t558) * MDP(14) + qJD(1) * t503 + (-t393 * t456 - t437 * t445 + t479) * MDP(16) + (-t392 * t456 + t435 * t445 - t478) * MDP(17) + ((t350 * t468 - t351 * t469) * pkin(3) + (t314 - t321) * t489 + (-t315 - t320) * t384) * MDP(18) + (-t314 * t320 - t315 * t321 + (t296 * t469 + t297 * t468 - t397 * t437) * pkin(3)) * MDP(19) + (t300 + t562) * MDP(22) + (-t301 - t563) * MDP(23) + (t483 * t507 + (t309 * t473 - t310 * t470) * t448 + t353 * t336 + (t448 * t484 - t292) * qJD(5) + t578) * MDP(25) + (-t484 * t507 - (t309 * t470 + t310 * t473) * t448 - t353 * t573 + (t448 * t483 - t548) * qJD(5) + t583) * MDP(26) + t582; (-t384 ^ 2 - t489 ^ 2) * MDP(18) + (-t314 * t384 - t315 * t489 + t383) * MDP(19) + (t301 - t563) * MDP(25) + (t300 - t562) * MDP(26); (t515 + t562) * MDP(22) + (-t499 - t563) * MDP(23) + (-t292 * t448 + t578) * MDP(25) + (-t291 * t448 + t583) * MDP(26) + (MDP(22) * t561 - MDP(23) * t573 - MDP(25) * t292 - MDP(26) * t548) * qJD(5) + t582;];
tauc = t1;
