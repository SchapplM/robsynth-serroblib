% Calculate kinetic energy for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:18
% EndTime: 2019-03-09 14:18:21
% DurationCPUTime: 2.76s
% Computational Cost: add. (3421->387), mult. (5799->592), div. (0->0), fcn. (7000->14), ass. (0->178)
t529 = sin(qJ(2));
t530 = sin(qJ(1));
t532 = cos(qJ(2));
t533 = cos(qJ(1));
t570 = cos(pkin(6));
t551 = t533 * t570;
t502 = t529 * t530 - t532 * t551;
t503 = t529 * t551 + t530 * t532;
t525 = sin(pkin(6));
t564 = t525 * t533;
t459 = Icges(3,5) * t503 - Icges(3,6) * t502 - Icges(3,3) * t564;
t552 = t530 * t570;
t504 = t533 * t529 + t532 * t552;
t505 = -t529 * t552 + t533 * t532;
t566 = t525 * t530;
t460 = Icges(3,5) * t505 - Icges(3,6) * t504 + Icges(3,3) * t566;
t579 = t525 * (t459 * t533 - t460 * t530);
t524 = sin(pkin(12));
t526 = cos(pkin(12));
t483 = -t503 * t524 - t526 * t564;
t556 = t524 * t564;
t484 = t503 * t526 - t556;
t427 = Icges(4,5) * t484 + Icges(4,6) * t483 + Icges(4,3) * t502;
t461 = Icges(3,4) * t503 - Icges(3,2) * t502 - Icges(3,6) * t564;
t578 = -t427 + t461;
t485 = -t505 * t524 + t526 * t566;
t557 = t524 * t566;
t486 = t505 * t526 + t557;
t428 = Icges(4,5) * t486 + Icges(4,6) * t485 + Icges(4,3) * t504;
t462 = Icges(3,4) * t505 - Icges(3,2) * t504 + Icges(3,6) * t566;
t577 = t428 - t462;
t567 = t525 * t529;
t500 = -t524 * t567 + t526 * t570;
t550 = t570 * t524;
t501 = t526 * t567 + t550;
t565 = t525 * t532;
t453 = Icges(4,5) * t501 + Icges(4,6) * t500 - Icges(4,3) * t565;
t491 = Icges(3,6) * t570 + (Icges(3,4) * t529 + Icges(3,2) * t532) * t525;
t576 = t453 - t491;
t573 = pkin(3) * t526;
t531 = cos(qJ(5));
t572 = pkin(5) * t531;
t528 = sin(qJ(5));
t569 = t502 * t528;
t568 = t504 * t528;
t473 = pkin(2) * t503 + qJ(3) * t502;
t474 = pkin(2) * t505 + qJ(3) * t504;
t560 = qJD(2) * t525;
t515 = t530 * t560;
t553 = t533 * t560;
t562 = t473 * t515 + t474 * t553;
t487 = qJD(4) * t504 + t515;
t561 = qJD(1) * (pkin(1) * t530 - pkin(8) * t564);
t559 = qJD(3) * t532;
t516 = qJD(2) * t570 + qJD(1);
t558 = pkin(12) + qJ(4);
t555 = t528 * t565;
t508 = qJD(1) * (pkin(1) * t533 + pkin(8) * t566);
t554 = qJD(3) * t502 + t516 * t474 + t508;
t520 = sin(t558);
t549 = cos(t558);
t543 = t525 * t549;
t481 = t505 * t520 - t530 * t543;
t438 = qJD(5) * t481 + t487;
t548 = qJD(3) * t504 - t561;
t506 = (pkin(2) * t529 - qJ(3) * t532) * t525;
t545 = (-rSges(4,1) * t501 - rSges(4,2) * t500 + rSges(4,3) * t565 - t506) * t560;
t544 = (-pkin(3) * t550 - (-pkin(9) * t532 + t529 * t573) * t525 - t506) * t560;
t488 = qJD(4) * t502 - t553;
t479 = t503 * t520 + t533 * t543;
t439 = qJD(5) * t479 + t488;
t507 = -qJD(4) * t565 + t516;
t493 = t520 * t567 - t549 * t570;
t470 = qJD(5) * t493 + t507;
t424 = -pkin(3) * t556 + pkin(9) * t502 + t503 * t573;
t425 = pkin(3) * t557 + pkin(9) * t504 + t505 * t573;
t541 = t424 * t515 + t425 * t553 - t525 * t559 + t562;
t540 = t516 * t425 + t530 * t544 + t554;
t480 = t503 * t549 - t520 * t564;
t436 = pkin(4) * t480 + pkin(10) * t479;
t482 = t505 * t549 + t520 * t566;
t437 = pkin(4) * t482 + pkin(10) * t481;
t539 = t487 * t436 - t437 * t488 + t541;
t494 = t520 * t570 + t529 * t543;
t457 = pkin(4) * t494 + pkin(10) * t493;
t538 = t507 * t437 - t457 * t487 + t540;
t537 = (-t424 - t473) * t516 + t533 * t544 + t548;
t536 = -t436 * t507 + t488 * t457 + t537;
t523 = qJ(5) + qJ(6);
t522 = cos(t523);
t521 = sin(t523);
t512 = rSges(2,1) * t533 - rSges(2,2) * t530;
t511 = rSges(2,1) * t530 + rSges(2,2) * t533;
t495 = t570 * rSges(3,3) + (rSges(3,1) * t529 + rSges(3,2) * t532) * t525;
t492 = Icges(3,5) * t570 + (Icges(3,1) * t529 + Icges(3,4) * t532) * t525;
t490 = Icges(3,3) * t570 + (Icges(3,5) * t529 + Icges(3,6) * t532) * t525;
t478 = t494 * t531 - t555;
t477 = -t494 * t528 - t531 * t565;
t472 = t494 * t522 - t521 * t565;
t471 = -t494 * t521 - t522 * t565;
t469 = rSges(3,1) * t505 - rSges(3,2) * t504 + rSges(3,3) * t566;
t468 = rSges(3,1) * t503 - rSges(3,2) * t502 - rSges(3,3) * t564;
t464 = Icges(3,1) * t505 - Icges(3,4) * t504 + Icges(3,5) * t566;
t463 = Icges(3,1) * t503 - Icges(3,4) * t502 - Icges(3,5) * t564;
t455 = Icges(4,1) * t501 + Icges(4,4) * t500 - Icges(4,5) * t565;
t454 = Icges(4,4) * t501 + Icges(4,2) * t500 - Icges(4,6) * t565;
t452 = rSges(5,1) * t494 - rSges(5,2) * t493 - rSges(5,3) * t565;
t451 = Icges(5,1) * t494 - Icges(5,4) * t493 - Icges(5,5) * t565;
t450 = Icges(5,4) * t494 - Icges(5,2) * t493 - Icges(5,6) * t565;
t449 = Icges(5,5) * t494 - Icges(5,6) * t493 - Icges(5,3) * t565;
t448 = t482 * t531 + t568;
t447 = -t482 * t528 + t504 * t531;
t446 = t480 * t531 + t569;
t445 = -t480 * t528 + t502 * t531;
t444 = t482 * t522 + t504 * t521;
t443 = -t482 * t521 + t504 * t522;
t442 = t480 * t522 + t502 * t521;
t441 = -t480 * t521 + t502 * t522;
t440 = qJD(6) * t493 + t470;
t434 = rSges(4,1) * t486 + rSges(4,2) * t485 + rSges(4,3) * t504;
t433 = rSges(4,1) * t484 + rSges(4,2) * t483 + rSges(4,3) * t502;
t432 = Icges(4,1) * t486 + Icges(4,4) * t485 + Icges(4,5) * t504;
t431 = Icges(4,1) * t484 + Icges(4,4) * t483 + Icges(4,5) * t502;
t430 = Icges(4,4) * t486 + Icges(4,2) * t485 + Icges(4,6) * t504;
t429 = Icges(4,4) * t484 + Icges(4,2) * t483 + Icges(4,6) * t502;
t423 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t504;
t422 = rSges(5,1) * t480 - rSges(5,2) * t479 + rSges(5,3) * t502;
t421 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t504;
t420 = Icges(5,1) * t480 - Icges(5,4) * t479 + Icges(5,5) * t502;
t419 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t504;
t418 = Icges(5,4) * t480 - Icges(5,2) * t479 + Icges(5,6) * t502;
t417 = Icges(5,5) * t482 - Icges(5,6) * t481 + Icges(5,3) * t504;
t416 = Icges(5,5) * t480 - Icges(5,6) * t479 + Icges(5,3) * t502;
t415 = t469 * t516 - t495 * t515 + t508;
t414 = -t468 * t516 - t495 * t553 - t561;
t413 = rSges(6,1) * t478 + rSges(6,2) * t477 + rSges(6,3) * t493;
t412 = Icges(6,1) * t478 + Icges(6,4) * t477 + Icges(6,5) * t493;
t411 = Icges(6,4) * t478 + Icges(6,2) * t477 + Icges(6,6) * t493;
t410 = Icges(6,5) * t478 + Icges(6,6) * t477 + Icges(6,3) * t493;
t406 = qJD(6) * t479 + t439;
t405 = qJD(6) * t481 + t438;
t403 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t493;
t402 = (t468 * t530 + t469 * t533) * t560;
t401 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t493;
t400 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t493;
t399 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t493;
t398 = -pkin(5) * t555 + pkin(11) * t493 + t494 * t572;
t397 = rSges(6,1) * t448 + rSges(6,2) * t447 + rSges(6,3) * t481;
t396 = rSges(6,1) * t446 + rSges(6,2) * t445 + rSges(6,3) * t479;
t395 = Icges(6,1) * t448 + Icges(6,4) * t447 + Icges(6,5) * t481;
t394 = Icges(6,1) * t446 + Icges(6,4) * t445 + Icges(6,5) * t479;
t393 = Icges(6,4) * t448 + Icges(6,2) * t447 + Icges(6,6) * t481;
t392 = Icges(6,4) * t446 + Icges(6,2) * t445 + Icges(6,6) * t479;
t391 = Icges(6,5) * t448 + Icges(6,6) * t447 + Icges(6,3) * t481;
t390 = Icges(6,5) * t446 + Icges(6,6) * t445 + Icges(6,3) * t479;
t389 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t481;
t388 = rSges(7,1) * t442 + rSges(7,2) * t441 + rSges(7,3) * t479;
t387 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t481;
t386 = Icges(7,1) * t442 + Icges(7,4) * t441 + Icges(7,5) * t479;
t385 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t481;
t384 = Icges(7,4) * t442 + Icges(7,2) * t441 + Icges(7,6) * t479;
t383 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t481;
t382 = Icges(7,5) * t442 + Icges(7,6) * t441 + Icges(7,3) * t479;
t381 = pkin(5) * t568 + pkin(11) * t481 + t482 * t572;
t380 = pkin(5) * t569 + pkin(11) * t479 + t480 * t572;
t379 = t434 * t516 + t530 * t545 + t554;
t378 = (-t433 - t473) * t516 + t533 * t545 + t548;
t377 = (-t559 + (t433 * t530 + t434 * t533) * qJD(2)) * t525 + t562;
t376 = t423 * t507 - t452 * t487 + t540;
t375 = -t422 * t507 + t452 * t488 + t537;
t374 = t422 * t487 - t423 * t488 + t541;
t373 = t397 * t470 - t413 * t438 + t538;
t372 = -t396 * t470 + t413 * t439 + t536;
t371 = t396 * t438 - t397 * t439 + t539;
t370 = t381 * t470 + t389 * t440 - t398 * t438 - t403 * t405 + t538;
t369 = -t380 * t470 - t388 * t440 + t398 * t439 + t403 * t406 + t536;
t368 = t380 * t438 - t381 * t439 + t388 * t405 - t389 * t406 + t539;
t1 = m(3) * (t402 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + t406 * ((t383 * t479 + t385 * t441 + t387 * t442) * t405 + (t479 * t382 + t441 * t384 + t442 * t386) * t406 + (t399 * t479 + t400 * t441 + t401 * t442) * t440) / 0.2e1 + t440 * ((t383 * t493 + t385 * t471 + t387 * t472) * t405 + (t382 * t493 + t384 * t471 + t386 * t472) * t406 + (t493 * t399 + t471 * t400 + t472 * t401) * t440) / 0.2e1 + t405 * ((t481 * t383 + t443 * t385 + t444 * t387) * t405 + (t382 * t481 + t384 * t443 + t386 * t444) * t406 + (t399 * t481 + t400 * t443 + t401 * t444) * t440) / 0.2e1 + t470 * ((t391 * t493 + t393 * t477 + t395 * t478) * t438 + (t390 * t493 + t392 * t477 + t394 * t478) * t439 + (t493 * t410 + t477 * t411 + t478 * t412) * t470) / 0.2e1 + t439 * ((t391 * t479 + t393 * t445 + t395 * t446) * t438 + (t479 * t390 + t445 * t392 + t446 * t394) * t439 + (t410 * t479 + t411 * t445 + t412 * t446) * t470) / 0.2e1 + t438 * ((t481 * t391 + t447 * t393 + t448 * t395) * t438 + (t390 * t481 + t392 * t447 + t394 * t448) * t439 + (t410 * t481 + t411 * t447 + t412 * t448) * t470) / 0.2e1 + t487 * ((t504 * t417 - t481 * t419 + t482 * t421) * t487 + (t416 * t504 - t418 * t481 + t420 * t482) * t488 + (t449 * t504 - t450 * t481 + t451 * t482) * t507) / 0.2e1 + t488 * ((t417 * t502 - t419 * t479 + t421 * t480) * t487 + (t502 * t416 - t479 * t418 + t480 * t420) * t488 + (t449 * t502 - t450 * t479 + t451 * t480) * t507) / 0.2e1 + t507 * ((-t417 * t565 - t419 * t493 + t421 * t494) * t487 + (-t416 * t565 - t418 * t493 + t420 * t494) * t488 + (-t449 * t565 - t493 * t450 + t494 * t451) * t507) / 0.2e1 + m(7) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(6) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(5) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(4) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + (((t430 * t500 + t432 * t501) * t530 - (t429 * t500 + t431 * t501) * t533 + (t427 * t533 - t428 * t530) * t565) * t560 + (t570 * t460 + (t462 * t532 + t464 * t529) * t525) * t515 - (t570 * t459 + (t461 * t532 + t463 * t529) * t525) * t553 + (-t453 * t565 + t500 * t454 + t501 * t455 + t570 * t490 + (t491 * t532 + t492 * t529) * t525) * t516) * t516 / 0.2e1 + (m(2) * (t511 ^ 2 + t512 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t429 * t485 - t431 * t486 - t463 * t505 + t504 * t578) * t533 + (t430 * t485 + t432 * t486 + t505 * t464 + t504 * t577 - t579) * t530) * t560 + (t454 * t485 + t455 * t486 + t490 * t566 + t492 * t505 + t504 * t576) * t516) * t515 / 0.2e1 - (((-t429 * t483 - t431 * t484 - t503 * t463 + t502 * t578 + t579) * t533 + (t430 * t483 + t432 * t484 + t464 * t503 + t502 * t577) * t530) * t560 + (t454 * t483 + t455 * t484 - t490 * t564 + t492 * t503 + t502 * t576) * t516) * t553 / 0.2e1;
T  = t1;
