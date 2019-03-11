% Calculate kinetic energy for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:00:56
% EndTime: 2019-03-09 23:00:59
% DurationCPUTime: 2.98s
% Computational Cost: add. (2883->343), mult. (5307->521), div. (0->0), fcn. (6300->12), ass. (0->159)
t570 = Icges(5,1) + Icges(6,2);
t569 = Icges(6,1) + Icges(5,3);
t568 = -Icges(5,4) - Icges(6,6);
t567 = Icges(6,4) - Icges(5,5);
t566 = Icges(6,5) - Icges(5,6);
t565 = Icges(5,2) + Icges(6,3);
t513 = sin(qJ(2));
t514 = sin(qJ(1));
t517 = cos(qJ(2));
t518 = cos(qJ(1));
t548 = cos(pkin(6));
t533 = t518 * t548;
t492 = t513 * t514 - t517 * t533;
t493 = t513 * t533 + t514 * t517;
t510 = sin(pkin(6));
t544 = t510 * t518;
t449 = Icges(3,5) * t493 - Icges(3,6) * t492 - Icges(3,3) * t544;
t534 = t514 * t548;
t494 = t518 * t513 + t517 * t534;
t495 = -t513 * t534 + t518 * t517;
t547 = t510 * t514;
t450 = Icges(3,5) * t495 - Icges(3,6) * t494 + Icges(3,3) * t547;
t564 = (t449 * t518 - t450 * t514) * t510;
t543 = qJ(3) + qJ(4);
t536 = cos(t543);
t531 = t510 * t536;
t535 = sin(t543);
t467 = t493 * t535 + t518 * t531;
t530 = t510 * t535;
t468 = t493 * t536 - t518 * t530;
t563 = t565 * t467 + t568 * t468 + t566 * t492;
t469 = t495 * t535 - t514 * t531;
t470 = t495 * t536 + t514 * t530;
t562 = t565 * t469 + t568 * t470 + t566 * t494;
t561 = t566 * t467 - t567 * t468 + t569 * t492;
t560 = t566 * t469 - t567 * t470 + t569 * t494;
t559 = t568 * t467 + t570 * t468 - t567 * t492;
t558 = t568 * t469 + t570 * t470 - t567 * t494;
t484 = t513 * t530 - t536 * t548;
t485 = t513 * t531 + t535 * t548;
t545 = t510 * t517;
t557 = t565 * t484 + t568 * t485 - t566 * t545;
t556 = t568 * t484 + t570 * t485 + t567 * t545;
t555 = t566 * t484 - t567 * t485 - t569 * t545;
t516 = cos(qJ(3));
t550 = pkin(3) * t516;
t546 = t510 * t516;
t461 = pkin(2) * t493 + pkin(9) * t492;
t462 = pkin(2) * t495 + pkin(9) * t494;
t540 = qJD(2) * t510;
t506 = t514 * t540;
t537 = t518 * t540;
t542 = t461 * t506 + t462 * t537;
t476 = qJD(3) * t494 + t506;
t541 = qJD(1) * (pkin(1) * t514 - pkin(8) * t544);
t507 = qJD(2) * t548 + qJD(1);
t512 = sin(qJ(3));
t539 = t512 * t547;
t538 = t512 * t544;
t445 = qJD(4) * t494 + t476;
t532 = t548 * t512;
t477 = qJD(3) * t492 - t537;
t412 = -pkin(3) * t538 + pkin(10) * t492 + t493 * t550;
t413 = pkin(3) * t539 + pkin(10) * t494 + t495 * t550;
t528 = t476 * t412 - t413 * t477 + t542;
t446 = qJD(4) * t492 + t477;
t496 = (pkin(2) * t513 - pkin(9) * t517) * t510;
t498 = qJD(1) * (pkin(1) * t518 + pkin(8) * t547);
t527 = t507 * t462 - t496 * t506 + t498;
t478 = (-qJD(3) - qJD(4)) * t545 + t507;
t425 = pkin(4) * t468 + qJ(5) * t467;
t526 = qJD(5) * t484 + t445 * t425 + t528;
t525 = -t461 * t507 - t496 * t537 - t541;
t455 = pkin(3) * t532 + (-pkin(10) * t517 + t513 * t550) * t510;
t497 = -qJD(3) * t545 + t507;
t524 = t497 * t413 - t455 * t476 + t527;
t426 = pkin(4) * t470 + qJ(5) * t469;
t523 = qJD(5) * t467 + t478 * t426 + t524;
t522 = -t412 * t497 + t477 * t455 + t525;
t447 = pkin(4) * t485 + qJ(5) * t484;
t521 = qJD(5) * t469 + t446 * t447 + t522;
t515 = cos(qJ(6));
t511 = sin(qJ(6));
t503 = rSges(2,1) * t518 - rSges(2,2) * t514;
t502 = rSges(2,1) * t514 + rSges(2,2) * t518;
t491 = t513 * t546 + t532;
t490 = -t510 * t513 * t512 + t516 * t548;
t483 = t548 * rSges(3,3) + (rSges(3,1) * t513 + rSges(3,2) * t517) * t510;
t482 = Icges(3,5) * t548 + (Icges(3,1) * t513 + Icges(3,4) * t517) * t510;
t481 = Icges(3,6) * t548 + (Icges(3,4) * t513 + Icges(3,2) * t517) * t510;
t480 = Icges(3,3) * t548 + (Icges(3,5) * t513 + Icges(3,6) * t517) * t510;
t475 = t495 * t516 + t539;
t474 = -t495 * t512 + t514 * t546;
t473 = t493 * t516 - t538;
t472 = -t493 * t512 - t516 * t544;
t471 = -pkin(5) * t545 + pkin(11) * t485;
t466 = t484 * t511 - t515 * t545;
t465 = t484 * t515 + t511 * t545;
t458 = rSges(3,1) * t495 - rSges(3,2) * t494 + rSges(3,3) * t547;
t457 = rSges(3,1) * t493 - rSges(3,2) * t492 - rSges(3,3) * t544;
t454 = Icges(3,1) * t495 - Icges(3,4) * t494 + Icges(3,5) * t547;
t453 = Icges(3,1) * t493 - Icges(3,4) * t492 - Icges(3,5) * t544;
t452 = Icges(3,4) * t495 - Icges(3,2) * t494 + Icges(3,6) * t547;
t451 = Icges(3,4) * t493 - Icges(3,2) * t492 - Icges(3,6) * t544;
t448 = rSges(4,1) * t491 + rSges(4,2) * t490 - rSges(4,3) * t545;
t444 = Icges(4,1) * t491 + Icges(4,4) * t490 - Icges(4,5) * t545;
t443 = Icges(4,4) * t491 + Icges(4,2) * t490 - Icges(4,6) * t545;
t442 = Icges(4,5) * t491 + Icges(4,6) * t490 - Icges(4,3) * t545;
t441 = qJD(6) * t485 + t478;
t440 = rSges(5,1) * t485 - rSges(5,2) * t484 - rSges(5,3) * t545;
t439 = -rSges(6,1) * t545 - rSges(6,2) * t485 + rSges(6,3) * t484;
t432 = pkin(5) * t494 + pkin(11) * t470;
t431 = pkin(5) * t492 + pkin(11) * t468;
t430 = t469 * t511 + t494 * t515;
t429 = t469 * t515 - t494 * t511;
t428 = t467 * t511 + t492 * t515;
t427 = t467 * t515 - t492 * t511;
t423 = qJD(6) * t468 + t446;
t422 = qJD(6) * t470 + t445;
t421 = rSges(4,1) * t475 + rSges(4,2) * t474 + rSges(4,3) * t494;
t420 = rSges(4,1) * t473 + rSges(4,2) * t472 + rSges(4,3) * t492;
t419 = Icges(4,1) * t475 + Icges(4,4) * t474 + Icges(4,5) * t494;
t418 = Icges(4,1) * t473 + Icges(4,4) * t472 + Icges(4,5) * t492;
t417 = Icges(4,4) * t475 + Icges(4,2) * t474 + Icges(4,6) * t494;
t416 = Icges(4,4) * t473 + Icges(4,2) * t472 + Icges(4,6) * t492;
t415 = Icges(4,5) * t475 + Icges(4,6) * t474 + Icges(4,3) * t494;
t414 = Icges(4,5) * t473 + Icges(4,6) * t472 + Icges(4,3) * t492;
t411 = rSges(5,1) * t470 - rSges(5,2) * t469 + rSges(5,3) * t494;
t410 = rSges(5,1) * t468 - rSges(5,2) * t467 + rSges(5,3) * t492;
t409 = rSges(6,1) * t494 - rSges(6,2) * t470 + rSges(6,3) * t469;
t408 = rSges(6,1) * t492 - rSges(6,2) * t468 + rSges(6,3) * t467;
t394 = rSges(7,1) * t466 + rSges(7,2) * t465 + rSges(7,3) * t485;
t393 = Icges(7,1) * t466 + Icges(7,4) * t465 + Icges(7,5) * t485;
t392 = Icges(7,4) * t466 + Icges(7,2) * t465 + Icges(7,6) * t485;
t391 = Icges(7,5) * t466 + Icges(7,6) * t465 + Icges(7,3) * t485;
t389 = t458 * t507 - t483 * t506 + t498;
t388 = -t457 * t507 - t483 * t537 - t541;
t386 = (t457 * t514 + t458 * t518) * t540;
t383 = rSges(7,1) * t430 + rSges(7,2) * t429 + rSges(7,3) * t470;
t382 = rSges(7,1) * t428 + rSges(7,2) * t427 + rSges(7,3) * t468;
t381 = Icges(7,1) * t430 + Icges(7,4) * t429 + Icges(7,5) * t470;
t380 = Icges(7,1) * t428 + Icges(7,4) * t427 + Icges(7,5) * t468;
t379 = Icges(7,4) * t430 + Icges(7,2) * t429 + Icges(7,6) * t470;
t378 = Icges(7,4) * t428 + Icges(7,2) * t427 + Icges(7,6) * t468;
t377 = Icges(7,5) * t430 + Icges(7,6) * t429 + Icges(7,3) * t470;
t376 = Icges(7,5) * t428 + Icges(7,6) * t427 + Icges(7,3) * t468;
t375 = t421 * t497 - t448 * t476 + t527;
t374 = -t420 * t497 + t448 * t477 + t525;
t373 = t420 * t476 - t421 * t477 + t542;
t372 = t411 * t478 - t440 * t445 + t524;
t371 = -t410 * t478 + t440 * t446 + t522;
t370 = t410 * t445 - t411 * t446 + t528;
t369 = t409 * t478 + (-t439 - t447) * t445 + t523;
t368 = t439 * t446 + (-t408 - t425) * t478 + t521;
t367 = t408 * t445 + (-t409 - t426) * t446 + t526;
t366 = t383 * t441 - t394 * t422 + t432 * t478 + (-t447 - t471) * t445 + t523;
t365 = -t382 * t441 + t394 * t423 + t446 * t471 + (-t425 - t431) * t478 + t521;
t364 = t382 * t422 - t383 * t423 + t431 * t445 + (-t426 - t432) * t446 + t526;
t1 = ((t480 * t547 - t481 * t494 + t482 * t495) * t507 + (-(-t451 * t494 + t453 * t495) * t518 + (-t494 * t452 + t495 * t454 - t564) * t514) * t540) * t506 / 0.2e1 - ((-t480 * t544 - t481 * t492 + t482 * t493) * t507 + ((-t452 * t492 + t454 * t493) * t514 + (t492 * t451 - t493 * t453 + t564) * t518) * t540) * t537 / 0.2e1 + t497 * ((-t415 * t545 + t417 * t490 + t419 * t491) * t476 + (-t414 * t545 + t416 * t490 + t418 * t491) * t477 + (-t442 * t545 + t443 * t490 + t444 * t491) * t497) / 0.2e1 + t422 * ((t470 * t377 + t429 * t379 + t430 * t381) * t422 + (t376 * t470 + t378 * t429 + t380 * t430) * t423 + (t391 * t470 + t392 * t429 + t393 * t430) * t441) / 0.2e1 + t423 * ((t377 * t468 + t379 * t427 + t381 * t428) * t422 + (t468 * t376 + t427 * t378 + t428 * t380) * t423 + (t391 * t468 + t392 * t427 + t393 * t428) * t441) / 0.2e1 + t441 * ((t377 * t485 + t379 * t465 + t381 * t466) * t422 + (t376 * t485 + t378 * t465 + t380 * t466) * t423 + (t485 * t391 + t465 * t392 + t466 * t393) * t441) / 0.2e1 + m(3) * (t386 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(4) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(5) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(6) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(7) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + t507 * ((t548 * t450 + (t452 * t517 + t454 * t513) * t510) * t506 - (t548 * t449 + (t451 * t517 + t453 * t513) * t510) * t537 + (t548 * t480 + (t481 * t517 + t482 * t513) * t510) * t507) / 0.2e1 + t476 * ((t494 * t415 + t474 * t417 + t475 * t419) * t476 + (t414 * t494 + t416 * t474 + t418 * t475) * t477 + (t442 * t494 + t443 * t474 + t444 * t475) * t497) / 0.2e1 + t477 * ((t415 * t492 + t417 * t472 + t419 * t473) * t476 + (t414 * t492 + t416 * t472 + t418 * t473) * t477 + (t442 * t492 + t443 * t472 + t444 * t473) * t497) / 0.2e1 + ((t469 * t557 + t470 * t556 + t494 * t555) * t478 + (t469 * t563 + t470 * t559 + t494 * t561) * t446 + (t562 * t469 + t558 * t470 + t560 * t494) * t445) * t445 / 0.2e1 + ((t467 * t557 + t468 * t556 + t492 * t555) * t478 + (t563 * t467 + t559 * t468 + t561 * t492) * t446 + (t467 * t562 + t468 * t558 + t492 * t560) * t445) * t446 / 0.2e1 + ((t557 * t484 + t556 * t485 - t555 * t545) * t478 + (t484 * t563 + t485 * t559 - t545 * t561) * t446 + (t484 * t562 + t485 * t558 - t545 * t560) * t445) * t478 / 0.2e1 + (m(2) * (t502 ^ 2 + t503 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
