% Calculate kinetic energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:49
% EndTime: 2019-03-09 19:24:52
% DurationCPUTime: 3.06s
% Computational Cost: add. (2701->338), mult. (6799->514), div. (0->0), fcn. (8423->12), ass. (0->152)
t572 = Icges(4,1) + Icges(5,1);
t571 = -Icges(4,4) + Icges(5,5);
t570 = Icges(5,4) + Icges(4,5);
t569 = Icges(4,2) + Icges(5,3);
t568 = Icges(5,2) + Icges(4,3);
t567 = Icges(4,6) - Icges(5,6);
t525 = sin(qJ(2));
t526 = sin(qJ(1));
t528 = cos(qJ(2));
t529 = cos(qJ(1));
t550 = cos(pkin(6));
t540 = t529 * t550;
t504 = t525 * t526 - t528 * t540;
t505 = t525 * t540 + t526 * t528;
t521 = sin(pkin(6));
t547 = t521 * t529;
t464 = Icges(3,5) * t505 - Icges(3,6) * t504 - Icges(3,3) * t547;
t541 = t526 * t550;
t506 = t529 * t525 + t528 * t541;
t507 = -t525 * t541 + t529 * t528;
t549 = t521 * t526;
t465 = Icges(3,5) * t507 - Icges(3,6) * t506 + Icges(3,3) * t549;
t566 = (t464 * t529 - t465 * t526) * t521;
t524 = sin(qJ(3));
t552 = cos(qJ(3));
t543 = t521 * t552;
t484 = t505 * t524 + t529 * t543;
t485 = t505 * t552 - t524 * t547;
t565 = t569 * t484 + t571 * t485 - t567 * t504;
t486 = t507 * t524 - t526 * t543;
t487 = t507 * t552 + t524 * t549;
t564 = t569 * t486 + t571 * t487 - t567 * t506;
t563 = -t567 * t484 + t570 * t485 + t568 * t504;
t562 = -t567 * t486 + t570 * t487 + t568 * t506;
t561 = t571 * t484 + t572 * t485 + t570 * t504;
t560 = t571 * t486 + t572 * t487 + t570 * t506;
t502 = t521 * t524 * t525 - t550 * t552;
t503 = t524 * t550 + t525 * t543;
t548 = t521 * t528;
t559 = t569 * t502 + t571 * t503 + t567 * t548;
t558 = -t567 * t502 + t570 * t503 - t568 * t548;
t557 = t571 * t502 + t572 * t503 - t570 * t548;
t551 = cos(qJ(5));
t478 = pkin(2) * t505 + pkin(9) * t504;
t479 = pkin(2) * t507 + pkin(9) * t506;
t544 = qJD(2) * t521;
t517 = t526 * t544;
t542 = t529 * t544;
t546 = t478 * t517 + t479 * t542;
t488 = qJD(3) * t506 + t517;
t545 = qJD(1) * (pkin(1) * t526 - pkin(8) * t547);
t518 = qJD(2) * t550 + qJD(1);
t460 = -qJD(5) * t506 + t488;
t447 = pkin(3) * t485 + qJ(4) * t484;
t539 = qJD(4) * t502 + t488 * t447 + t546;
t489 = qJD(3) * t504 - t542;
t461 = -qJD(5) * t504 + t489;
t509 = -qJD(3) * t548 + t518;
t508 = (pkin(2) * t525 - pkin(9) * t528) * t521;
t510 = qJD(1) * (pkin(1) * t529 + pkin(8) * t549);
t537 = t518 * t479 - t508 * t517 + t510;
t491 = qJD(5) * t548 + t509;
t448 = pkin(3) * t487 + qJ(4) * t486;
t536 = qJD(4) * t484 + t509 * t448 + t537;
t452 = pkin(4) * t485 - pkin(10) * t504;
t453 = pkin(4) * t487 - pkin(10) * t506;
t535 = t488 * t452 + (-t448 - t453) * t489 + t539;
t534 = -t478 * t518 - t508 * t542 - t545;
t477 = pkin(3) * t503 + qJ(4) * t502;
t533 = qJD(4) * t486 + t489 * t477 + t534;
t490 = pkin(4) * t503 + pkin(10) * t548;
t532 = t509 * t453 + (-t477 - t490) * t488 + t536;
t531 = t489 * t490 + (-t447 - t452) * t509 + t533;
t527 = cos(qJ(6));
t523 = sin(qJ(5));
t522 = sin(qJ(6));
t513 = rSges(2,1) * t529 - rSges(2,2) * t526;
t512 = rSges(2,1) * t526 + rSges(2,2) * t529;
t495 = t550 * rSges(3,3) + (rSges(3,1) * t525 + rSges(3,2) * t528) * t521;
t494 = Icges(3,5) * t550 + (Icges(3,1) * t525 + Icges(3,4) * t528) * t521;
t493 = Icges(3,6) * t550 + (Icges(3,4) * t525 + Icges(3,2) * t528) * t521;
t492 = Icges(3,3) * t550 + (Icges(3,5) * t525 + Icges(3,6) * t528) * t521;
t476 = t502 * t523 + t503 * t551;
t475 = -t502 * t551 + t503 * t523;
t472 = rSges(3,1) * t507 - rSges(3,2) * t506 + rSges(3,3) * t549;
t471 = rSges(3,1) * t505 - rSges(3,2) * t504 - rSges(3,3) * t547;
t469 = Icges(3,1) * t507 - Icges(3,4) * t506 + Icges(3,5) * t549;
t468 = Icges(3,1) * t505 - Icges(3,4) * t504 - Icges(3,5) * t547;
t467 = Icges(3,4) * t507 - Icges(3,2) * t506 + Icges(3,6) * t549;
t466 = Icges(3,4) * t505 - Icges(3,2) * t504 - Icges(3,6) * t547;
t463 = rSges(4,1) * t503 - rSges(4,2) * t502 - rSges(4,3) * t548;
t462 = rSges(5,1) * t503 - rSges(5,2) * t548 + rSges(5,3) * t502;
t450 = t476 * t527 + t522 * t548;
t449 = -t476 * t522 + t527 * t548;
t444 = qJD(6) * t475 + t491;
t443 = t486 * t523 + t487 * t551;
t442 = -t486 * t551 + t487 * t523;
t441 = t484 * t523 + t485 * t551;
t440 = -t484 * t551 + t485 * t523;
t437 = rSges(4,1) * t487 - rSges(4,2) * t486 + rSges(4,3) * t506;
t436 = rSges(5,1) * t487 + rSges(5,2) * t506 + rSges(5,3) * t486;
t435 = rSges(4,1) * t485 - rSges(4,2) * t484 + rSges(4,3) * t504;
t434 = rSges(5,1) * t485 + rSges(5,2) * t504 + rSges(5,3) * t484;
t421 = pkin(5) * t476 + pkin(11) * t475;
t420 = t443 * t527 - t506 * t522;
t419 = -t443 * t522 - t506 * t527;
t418 = t441 * t527 - t504 * t522;
t417 = -t441 * t522 - t504 * t527;
t415 = rSges(6,1) * t476 - rSges(6,2) * t475 + rSges(6,3) * t548;
t414 = Icges(6,1) * t476 - Icges(6,4) * t475 + Icges(6,5) * t548;
t413 = Icges(6,4) * t476 - Icges(6,2) * t475 + Icges(6,6) * t548;
t412 = Icges(6,5) * t476 - Icges(6,6) * t475 + Icges(6,3) * t548;
t411 = t472 * t518 - t495 * t517 + t510;
t410 = -t471 * t518 - t495 * t542 - t545;
t409 = (t471 * t526 + t472 * t529) * t544;
t408 = qJD(6) * t440 + t461;
t407 = qJD(6) * t442 + t460;
t406 = pkin(5) * t443 + pkin(11) * t442;
t405 = pkin(5) * t441 + pkin(11) * t440;
t404 = rSges(7,1) * t450 + rSges(7,2) * t449 + rSges(7,3) * t475;
t403 = Icges(7,1) * t450 + Icges(7,4) * t449 + Icges(7,5) * t475;
t402 = Icges(7,4) * t450 + Icges(7,2) * t449 + Icges(7,6) * t475;
t401 = Icges(7,5) * t450 + Icges(7,6) * t449 + Icges(7,3) * t475;
t400 = rSges(6,1) * t443 - rSges(6,2) * t442 - rSges(6,3) * t506;
t399 = rSges(6,1) * t441 - rSges(6,2) * t440 - rSges(6,3) * t504;
t398 = Icges(6,1) * t443 - Icges(6,4) * t442 - Icges(6,5) * t506;
t397 = Icges(6,1) * t441 - Icges(6,4) * t440 - Icges(6,5) * t504;
t396 = Icges(6,4) * t443 - Icges(6,2) * t442 - Icges(6,6) * t506;
t395 = Icges(6,4) * t441 - Icges(6,2) * t440 - Icges(6,6) * t504;
t394 = Icges(6,5) * t443 - Icges(6,6) * t442 - Icges(6,3) * t506;
t393 = Icges(6,5) * t441 - Icges(6,6) * t440 - Icges(6,3) * t504;
t392 = t437 * t509 - t463 * t488 + t537;
t391 = -t435 * t509 + t463 * t489 + t534;
t390 = rSges(7,1) * t420 + rSges(7,2) * t419 + rSges(7,3) * t442;
t389 = rSges(7,1) * t418 + rSges(7,2) * t417 + rSges(7,3) * t440;
t388 = Icges(7,1) * t420 + Icges(7,4) * t419 + Icges(7,5) * t442;
t387 = Icges(7,1) * t418 + Icges(7,4) * t417 + Icges(7,5) * t440;
t386 = Icges(7,4) * t420 + Icges(7,2) * t419 + Icges(7,6) * t442;
t385 = Icges(7,4) * t418 + Icges(7,2) * t417 + Icges(7,6) * t440;
t384 = Icges(7,5) * t420 + Icges(7,6) * t419 + Icges(7,3) * t442;
t383 = Icges(7,5) * t418 + Icges(7,6) * t417 + Icges(7,3) * t440;
t382 = t435 * t488 - t437 * t489 + t546;
t381 = t436 * t509 + (-t462 - t477) * t488 + t536;
t380 = t462 * t489 + (-t434 - t447) * t509 + t533;
t379 = t434 * t488 + (-t436 - t448) * t489 + t539;
t378 = t400 * t491 - t415 * t460 + t532;
t377 = -t399 * t491 + t415 * t461 + t531;
t376 = t399 * t460 - t400 * t461 + t535;
t375 = t390 * t444 - t404 * t407 + t406 * t491 - t421 * t460 + t532;
t374 = -t389 * t444 + t404 * t408 - t405 * t491 + t421 * t461 + t531;
t373 = t389 * t407 - t390 * t408 + t405 * t460 - t406 * t461 + t535;
t1 = t518 * ((t550 * t465 + (t467 * t528 + t469 * t525) * t521) * t517 - (t550 * t464 + (t466 * t528 + t468 * t525) * t521) * t542 + (t550 * t492 + (t493 * t528 + t494 * t525) * t521) * t518) / 0.2e1 + t460 * ((-t394 * t506 - t396 * t442 + t398 * t443) * t460 + (-t393 * t506 - t395 * t442 + t397 * t443) * t461 + (-t412 * t506 - t413 * t442 + t414 * t443) * t491) / 0.2e1 + t461 * ((-t394 * t504 - t396 * t440 + t398 * t441) * t460 + (-t393 * t504 - t395 * t440 + t397 * t441) * t461 + (-t412 * t504 - t413 * t440 + t414 * t441) * t491) / 0.2e1 + t491 * ((t394 * t548 - t396 * t475 + t398 * t476) * t460 + (t393 * t548 - t395 * t475 + t397 * t476) * t461 + (t412 * t548 - t413 * t475 + t414 * t476) * t491) / 0.2e1 + t407 * ((t442 * t384 + t419 * t386 + t420 * t388) * t407 + (t383 * t442 + t385 * t419 + t387 * t420) * t408 + (t401 * t442 + t402 * t419 + t403 * t420) * t444) / 0.2e1 + t408 * ((t384 * t440 + t386 * t417 + t388 * t418) * t407 + (t440 * t383 + t417 * t385 + t418 * t387) * t408 + (t401 * t440 + t402 * t417 + t403 * t418) * t444) / 0.2e1 + t444 * ((t384 * t475 + t386 * t449 + t388 * t450) * t407 + (t383 * t475 + t385 * t449 + t387 * t450) * t408 + (t475 * t401 + t449 * t402 + t450 * t403) * t444) / 0.2e1 + m(3) * (t409 ^ 2 + t410 ^ 2 + t411 ^ 2) / 0.2e1 + m(4) * (t382 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(5) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(6) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(7) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + ((t492 * t549 - t493 * t506 + t494 * t507) * t518 + (-(-t466 * t506 + t468 * t507) * t529 + (-t506 * t467 + t507 * t469 - t566) * t526) * t544) * t517 / 0.2e1 - ((-t492 * t547 - t493 * t504 + t494 * t505) * t518 + ((-t467 * t504 + t469 * t505) * t526 + (t504 * t466 - t505 * t468 + t566) * t529) * t544) * t542 / 0.2e1 + ((t486 * t559 + t487 * t557 + t506 * t558) * t509 + (t486 * t565 + t487 * t561 + t506 * t563) * t489 + (t564 * t486 + t560 * t487 + t562 * t506) * t488) * t488 / 0.2e1 + ((t484 * t559 + t485 * t557 + t504 * t558) * t509 + (t565 * t484 + t561 * t485 + t563 * t504) * t489 + (t484 * t564 + t485 * t560 + t504 * t562) * t488) * t489 / 0.2e1 + ((t559 * t502 + t557 * t503 - t558 * t548) * t509 + (t502 * t565 + t561 * t503 - t563 * t548) * t489 + (t502 * t564 + t503 * t560 - t548 * t562) * t488) * t509 / 0.2e1 + (m(2) * (t512 ^ 2 + t513 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
