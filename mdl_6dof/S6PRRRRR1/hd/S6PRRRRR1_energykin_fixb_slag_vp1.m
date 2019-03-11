% Calculate kinetic energy for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:31
% EndTime: 2019-03-09 00:37:33
% DurationCPUTime: 2.44s
% Computational Cost: add. (3454->382), mult. (5503->602), div. (0->0), fcn. (6530->14), ass. (0->174)
t573 = qJD(2) ^ 2;
t537 = cos(qJ(3));
t571 = t537 * pkin(3);
t529 = sin(pkin(12));
t530 = sin(pkin(6));
t569 = t529 * t530;
t531 = cos(pkin(12));
t568 = t530 * t531;
t534 = sin(qJ(3));
t567 = t530 * t534;
t535 = sin(qJ(2));
t566 = t530 * t535;
t565 = t530 * t537;
t538 = cos(qJ(2));
t564 = t530 * t538;
t532 = cos(pkin(6));
t563 = t532 * t534;
t562 = t532 * t535;
t561 = t532 * t538;
t528 = qJ(3) + qJ(4);
t525 = cos(t528);
t560 = pkin(4) * t525;
t509 = t529 * t561 + t531 * t535;
t558 = qJD(2) * t530;
t520 = t529 * t558;
t490 = qJD(3) * t509 + t520;
t523 = qJD(2) * t532;
t556 = -qJD(3) - qJD(4);
t555 = t529 * t567;
t554 = t531 * t567;
t553 = qJ(5) + t528;
t455 = qJD(4) * t509 + t490;
t552 = t531 * t558;
t507 = t529 * t535 - t531 * t561;
t508 = t529 * t538 + t531 * t562;
t473 = t508 * pkin(2) + t507 * pkin(8);
t510 = -t529 * t562 + t531 * t538;
t474 = t510 * pkin(2) + t509 * pkin(8);
t551 = t473 * t520 + t474 * t552 + qJD(1);
t524 = sin(t528);
t550 = pkin(4) * t524;
t439 = qJD(5) * t509 + t455;
t549 = cos(t553);
t491 = qJD(3) * t507 - t552;
t548 = t530 * t549;
t513 = (pkin(2) * t535 - pkin(8) * t538) * t530;
t547 = t474 * t523 - t513 * t520;
t456 = qJD(4) * t507 + t491;
t440 = qJD(5) * t507 + t456;
t423 = -pkin(3) * t554 + pkin(9) * t507 + t508 * t571;
t424 = pkin(3) * t555 + pkin(9) * t509 + t510 * t571;
t546 = t490 * t423 - t424 * t491 + t551;
t489 = t523 + (-qJD(5) + t556) * t564;
t545 = (-t473 * t532 - t513 * t568) * qJD(2);
t471 = pkin(3) * t563 + (-pkin(9) * t538 + t535 * t571) * t530;
t514 = -qJD(3) * t564 + t523;
t544 = t514 * t424 - t471 * t490 + t547;
t395 = pkin(10) * t507 + t508 * t560 - t550 * t568;
t396 = pkin(10) * t509 + t510 * t560 + t550 * t569;
t543 = t455 * t395 - t396 * t456 + t546;
t542 = -t423 * t514 + t491 * t471 + t545;
t437 = t550 * t532 + (-pkin(10) * t538 + t535 * t560) * t530;
t494 = t556 * t564 + t523;
t541 = t494 * t396 - t437 * t455 + t544;
t540 = -t395 * t494 + t456 * t437 + t542;
t536 = cos(qJ(6));
t533 = sin(qJ(6));
t521 = sin(t553);
t512 = t535 * t565 + t563;
t511 = t532 * t537 - t534 * t566;
t500 = t524 * t532 + t525 * t566;
t499 = -t524 * t566 + t525 * t532;
t498 = rSges(3,3) * t532 + (rSges(3,1) * t535 + rSges(3,2) * t538) * t530;
t497 = Icges(3,5) * t532 + (Icges(3,1) * t535 + Icges(3,4) * t538) * t530;
t496 = Icges(3,6) * t532 + (Icges(3,4) * t535 + Icges(3,2) * t538) * t530;
t495 = Icges(3,3) * t532 + (Icges(3,5) * t535 + Icges(3,6) * t538) * t530;
t493 = t532 * t521 + t535 * t548;
t492 = t521 * t566 - t532 * t549;
t488 = t510 * t537 + t555;
t487 = -t510 * t534 + t529 * t565;
t486 = t508 * t537 - t554;
t485 = -t508 * t534 - t531 * t565;
t484 = t510 * t525 + t524 * t569;
t483 = -t510 * t524 + t525 * t569;
t482 = t508 * t525 - t524 * t568;
t481 = -t508 * t524 - t525 * t568;
t480 = t493 * t536 - t533 * t564;
t479 = -t493 * t533 - t536 * t564;
t478 = t510 * t549 + t521 * t569;
t477 = t510 * t521 - t529 * t548;
t476 = t508 * t549 - t521 * t568;
t475 = t508 * t521 + t531 * t548;
t470 = rSges(4,1) * t512 + rSges(4,2) * t511 - rSges(4,3) * t564;
t469 = Icges(4,1) * t512 + Icges(4,4) * t511 - Icges(4,5) * t564;
t468 = Icges(4,4) * t512 + Icges(4,2) * t511 - Icges(4,6) * t564;
t467 = Icges(4,5) * t512 + Icges(4,6) * t511 - Icges(4,3) * t564;
t464 = rSges(3,1) * t510 - rSges(3,2) * t509 + rSges(3,3) * t569;
t463 = rSges(3,1) * t508 - rSges(3,2) * t507 - rSges(3,3) * t568;
t462 = Icges(3,1) * t510 - Icges(3,4) * t509 + Icges(3,5) * t569;
t461 = Icges(3,1) * t508 - Icges(3,4) * t507 - Icges(3,5) * t568;
t460 = Icges(3,4) * t510 - Icges(3,2) * t509 + Icges(3,6) * t569;
t459 = Icges(3,4) * t508 - Icges(3,2) * t507 - Icges(3,6) * t568;
t458 = Icges(3,5) * t510 - Icges(3,6) * t509 + Icges(3,3) * t569;
t457 = Icges(3,5) * t508 - Icges(3,6) * t507 - Icges(3,3) * t568;
t454 = pkin(5) * t493 + pkin(11) * t492;
t453 = rSges(5,1) * t500 + rSges(5,2) * t499 - rSges(5,3) * t564;
t452 = Icges(5,1) * t500 + Icges(5,4) * t499 - Icges(5,5) * t564;
t451 = Icges(5,4) * t500 + Icges(5,2) * t499 - Icges(5,6) * t564;
t450 = Icges(5,5) * t500 + Icges(5,6) * t499 - Icges(5,3) * t564;
t449 = qJD(6) * t492 + t489;
t448 = rSges(6,1) * t493 - rSges(6,2) * t492 - rSges(6,3) * t564;
t447 = Icges(6,1) * t493 - Icges(6,4) * t492 - Icges(6,5) * t564;
t446 = Icges(6,4) * t493 - Icges(6,2) * t492 - Icges(6,6) * t564;
t445 = Icges(6,5) * t493 - Icges(6,6) * t492 - Icges(6,3) * t564;
t444 = t478 * t536 + t509 * t533;
t443 = -t478 * t533 + t509 * t536;
t442 = t476 * t536 + t507 * t533;
t441 = -t476 * t533 + t507 * t536;
t436 = pkin(5) * t478 + pkin(11) * t477;
t435 = pkin(5) * t476 + pkin(11) * t475;
t434 = (-t463 * t532 - t498 * t568) * qJD(2);
t433 = (t464 * t532 - t498 * t569) * qJD(2);
t432 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t509;
t431 = rSges(4,1) * t486 + rSges(4,2) * t485 + rSges(4,3) * t507;
t430 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t509;
t429 = Icges(4,1) * t486 + Icges(4,4) * t485 + Icges(4,5) * t507;
t428 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t509;
t427 = Icges(4,4) * t486 + Icges(4,2) * t485 + Icges(4,6) * t507;
t426 = Icges(4,5) * t488 + Icges(4,6) * t487 + Icges(4,3) * t509;
t425 = Icges(4,5) * t486 + Icges(4,6) * t485 + Icges(4,3) * t507;
t422 = rSges(5,1) * t484 + rSges(5,2) * t483 + rSges(5,3) * t509;
t421 = rSges(5,1) * t482 + rSges(5,2) * t481 + rSges(5,3) * t507;
t420 = Icges(5,1) * t484 + Icges(5,4) * t483 + Icges(5,5) * t509;
t419 = Icges(5,1) * t482 + Icges(5,4) * t481 + Icges(5,5) * t507;
t418 = Icges(5,4) * t484 + Icges(5,2) * t483 + Icges(5,6) * t509;
t417 = Icges(5,4) * t482 + Icges(5,2) * t481 + Icges(5,6) * t507;
t416 = Icges(5,5) * t484 + Icges(5,6) * t483 + Icges(5,3) * t509;
t415 = Icges(5,5) * t482 + Icges(5,6) * t481 + Icges(5,3) * t507;
t414 = rSges(6,1) * t478 - rSges(6,2) * t477 + rSges(6,3) * t509;
t413 = rSges(6,1) * t476 - rSges(6,2) * t475 + rSges(6,3) * t507;
t412 = Icges(6,1) * t478 - Icges(6,4) * t477 + Icges(6,5) * t509;
t411 = Icges(6,1) * t476 - Icges(6,4) * t475 + Icges(6,5) * t507;
t410 = Icges(6,4) * t478 - Icges(6,2) * t477 + Icges(6,6) * t509;
t409 = Icges(6,4) * t476 - Icges(6,2) * t475 + Icges(6,6) * t507;
t408 = Icges(6,5) * t478 - Icges(6,6) * t477 + Icges(6,3) * t509;
t407 = Icges(6,5) * t476 - Icges(6,6) * t475 + Icges(6,3) * t507;
t406 = rSges(7,1) * t480 + rSges(7,2) * t479 + rSges(7,3) * t492;
t405 = Icges(7,1) * t480 + Icges(7,4) * t479 + Icges(7,5) * t492;
t404 = Icges(7,4) * t480 + Icges(7,2) * t479 + Icges(7,6) * t492;
t403 = Icges(7,5) * t480 + Icges(7,6) * t479 + Icges(7,3) * t492;
t402 = qJD(6) * t475 + t440;
t401 = qJD(6) * t477 + t439;
t399 = qJD(1) + (t463 * t529 + t464 * t531) * t558;
t393 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t477;
t392 = rSges(7,1) * t442 + rSges(7,2) * t441 + rSges(7,3) * t475;
t391 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t477;
t390 = Icges(7,1) * t442 + Icges(7,4) * t441 + Icges(7,5) * t475;
t389 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t477;
t388 = Icges(7,4) * t442 + Icges(7,2) * t441 + Icges(7,6) * t475;
t387 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t477;
t386 = Icges(7,5) * t442 + Icges(7,6) * t441 + Icges(7,3) * t475;
t384 = -t431 * t514 + t470 * t491 + t545;
t383 = t432 * t514 - t470 * t490 + t547;
t382 = t431 * t490 - t432 * t491 + t551;
t381 = -t421 * t494 + t453 * t456 + t542;
t380 = t422 * t494 - t453 * t455 + t544;
t379 = t421 * t455 - t422 * t456 + t546;
t378 = -t413 * t489 + t440 * t448 + t540;
t377 = t414 * t489 - t439 * t448 + t541;
t376 = t413 * t439 - t414 * t440 + t543;
t375 = -t392 * t449 + t402 * t406 - t435 * t489 + t440 * t454 + t540;
t374 = t393 * t449 - t401 * t406 + t436 * t489 - t439 * t454 + t541;
t373 = t392 * t401 - t393 * t402 + t435 * t439 - t436 * t440 + t543;
t1 = t402 * ((t387 * t475 + t389 * t441 + t391 * t442) * t401 + (t475 * t386 + t441 * t388 + t442 * t390) * t402 + (t403 * t475 + t404 * t441 + t405 * t442) * t449) / 0.2e1 + t401 * ((t387 * t477 + t389 * t443 + t391 * t444) * t401 + (t386 * t477 + t388 * t443 + t390 * t444) * t402 + (t403 * t477 + t404 * t443 + t405 * t444) * t449) / 0.2e1 + t449 * ((t387 * t492 + t389 * t479 + t391 * t480) * t401 + (t386 * t492 + t388 * t479 + t390 * t480) * t402 + (t492 * t403 + t479 * t404 + t480 * t405) * t449) / 0.2e1 + t440 * ((t408 * t507 - t410 * t475 + t412 * t476) * t439 + (t507 * t407 - t475 * t409 + t476 * t411) * t440 + (t445 * t507 - t446 * t475 + t447 * t476) * t489) / 0.2e1 + t489 * ((-t408 * t564 - t410 * t492 + t412 * t493) * t439 + (-t407 * t564 - t409 * t492 + t411 * t493) * t440 + (-t445 * t564 - t446 * t492 + t447 * t493) * t489) / 0.2e1 + t439 * ((t509 * t408 - t477 * t410 + t478 * t412) * t439 + (t407 * t509 - t409 * t477 + t411 * t478) * t440 + (t445 * t509 - t446 * t477 + t447 * t478) * t489) / 0.2e1 + t455 * ((t509 * t416 + t483 * t418 + t484 * t420) * t455 + (t415 * t509 + t417 * t483 + t419 * t484) * t456 + (t450 * t509 + t451 * t483 + t452 * t484) * t494) / 0.2e1 + t456 * ((t416 * t507 + t418 * t481 + t420 * t482) * t455 + (t507 * t415 + t481 * t417 + t482 * t419) * t456 + (t450 * t507 + t451 * t481 + t452 * t482) * t494) / 0.2e1 + t494 * ((-t416 * t564 + t418 * t499 + t420 * t500) * t455 + (-t415 * t564 + t417 * t499 + t419 * t500) * t456 + (-t450 * t564 + t451 * t499 + t452 * t500) * t494) / 0.2e1 + t491 * ((t426 * t507 + t428 * t485 + t430 * t486) * t490 + (t425 * t507 + t427 * t485 + t429 * t486) * t491 + (t467 * t507 + t468 * t485 + t469 * t486) * t514) / 0.2e1 + t490 * ((t426 * t509 + t428 * t487 + t430 * t488) * t490 + (t425 * t509 + t427 * t487 + t429 * t488) * t491 + (t467 * t509 + t468 * t487 + t469 * t488) * t514) / 0.2e1 + t514 * ((-t426 * t564 + t428 * t511 + t430 * t512) * t490 + (-t425 * t564 + t427 * t511 + t429 * t512) * t491 + (-t467 * t564 + t468 * t511 + t469 * t512) * t514) / 0.2e1 + m(7) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(6) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(5) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(4) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 + m(3) * (t399 ^ 2 + t433 ^ 2 + t434 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 - t573 * ((-t458 * t568 - t460 * t507 + t462 * t508) * t569 - (-t457 * t568 - t459 * t507 + t461 * t508) * t568 + (-t495 * t568 - t496 * t507 + t497 * t508) * t532) * t568 / 0.2e1 + (t532 * (t532 ^ 2 * t495 + (((t460 * t538 + t462 * t535) * t529 - (t459 * t538 + t461 * t535) * t531) * t530 + (-t457 * t531 + t458 * t529 + t496 * t538 + t497 * t535) * t532) * t530) + ((t458 * t569 - t460 * t509 + t462 * t510) * t569 - (t457 * t569 - t459 * t509 + t461 * t510) * t568 + (t495 * t569 - t496 * t509 + t497 * t510) * t532) * t569) * t573 / 0.2e1;
T  = t1;
