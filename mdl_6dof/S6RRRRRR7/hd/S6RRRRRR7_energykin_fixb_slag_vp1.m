% Calculate kinetic energy for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:28:59
% EndTime: 2019-03-10 04:29:02
% DurationCPUTime: 2.46s
% Computational Cost: add. (3365->388), mult. (6949->607), div. (0->0), fcn. (8532->14), ass. (0->174)
t530 = sin(qJ(2));
t531 = sin(qJ(1));
t533 = cos(qJ(2));
t534 = cos(qJ(1));
t562 = cos(pkin(6));
t546 = t534 * t562;
t500 = t530 * t531 - t533 * t546;
t501 = t530 * t546 + t531 * t533;
t527 = sin(pkin(6));
t557 = t527 * t534;
t459 = Icges(3,5) * t501 - Icges(3,6) * t500 - Icges(3,3) * t557;
t547 = t531 * t562;
t502 = t534 * t530 + t533 * t547;
t503 = -t530 * t547 + t534 * t533;
t559 = t527 * t531;
t460 = Icges(3,5) * t503 - Icges(3,6) * t502 + Icges(3,3) * t559;
t568 = (t459 * t534 - t460 * t531) * t527;
t566 = cos(qJ(3));
t532 = cos(qJ(4));
t564 = t532 * pkin(4);
t528 = sin(qJ(4));
t561 = t500 * t528;
t560 = t502 * t528;
t558 = t527 * t533;
t526 = qJ(4) + qJ(5);
t473 = pkin(2) * t501 + pkin(9) * t500;
t474 = pkin(2) * t503 + pkin(9) * t502;
t552 = qJD(2) * t527;
t514 = t531 * t552;
t549 = t534 * t552;
t556 = t473 * t514 + t474 * t549;
t522 = cos(t526);
t555 = pkin(5) * t522;
t488 = qJD(3) * t502 + t514;
t553 = qJD(1) * (pkin(1) * t531 - pkin(8) * t557);
t515 = qJD(2) * t562 + qJD(1);
t551 = t528 * t558;
t529 = sin(qJ(3));
t550 = t527 * t566;
t486 = t503 * t529 - t531 * t550;
t444 = qJD(4) * t486 + t488;
t521 = sin(t526);
t548 = pkin(5) * t521;
t419 = qJD(5) * t486 + t444;
t489 = qJD(3) * t500 - t549;
t484 = t501 * t529 + t534 * t550;
t485 = t501 * t566 - t529 * t557;
t438 = pkin(3) * t485 + pkin(10) * t484;
t487 = t503 * t566 + t529 * t559;
t439 = pkin(3) * t487 + pkin(10) * t486;
t544 = t488 * t438 - t439 * t489 + t556;
t445 = qJD(4) * t484 + t489;
t505 = -qJD(3) * t558 + t515;
t504 = (pkin(2) * t530 - pkin(9) * t533) * t527;
t506 = qJD(1) * (pkin(1) * t534 + pkin(8) * t559);
t543 = t515 * t474 - t504 * t514 + t506;
t420 = qJD(5) * t484 + t445;
t498 = t527 * t529 * t530 - t562 * t566;
t475 = qJD(4) * t498 + t505;
t454 = qJD(5) * t498 + t475;
t386 = pkin(4) * t561 + pkin(11) * t484 + t485 * t564;
t387 = pkin(4) * t560 + pkin(11) * t486 + t487 * t564;
t542 = t444 * t386 - t387 * t445 + t544;
t541 = -t473 * t515 - t504 * t549 - t553;
t499 = t529 * t562 + t530 * t550;
t470 = pkin(3) * t499 + pkin(10) * t498;
t540 = t505 * t439 - t470 * t488 + t543;
t539 = -t438 * t505 + t489 * t470 + t541;
t422 = -pkin(4) * t551 + pkin(11) * t498 + t499 * t564;
t538 = t475 * t387 - t422 * t444 + t540;
t537 = -t386 * t475 + t445 * t422 + t539;
t523 = qJ(6) + t526;
t518 = cos(t523);
t517 = sin(t523);
t511 = rSges(2,1) * t534 - rSges(2,2) * t531;
t510 = rSges(2,1) * t531 + rSges(2,2) * t534;
t493 = t562 * rSges(3,3) + (rSges(3,1) * t530 + rSges(3,2) * t533) * t527;
t492 = Icges(3,5) * t562 + (Icges(3,1) * t530 + Icges(3,4) * t533) * t527;
t491 = Icges(3,6) * t562 + (Icges(3,4) * t530 + Icges(3,2) * t533) * t527;
t490 = Icges(3,3) * t562 + (Icges(3,5) * t530 + Icges(3,6) * t533) * t527;
t483 = t499 * t532 - t551;
t482 = -t499 * t528 - t532 * t558;
t477 = t499 * t522 - t521 * t558;
t476 = -t499 * t521 - t522 * t558;
t472 = t499 * t518 - t517 * t558;
t471 = -t499 * t517 - t518 * t558;
t467 = rSges(3,1) * t503 - rSges(3,2) * t502 + rSges(3,3) * t559;
t466 = rSges(3,1) * t501 - rSges(3,2) * t500 - rSges(3,3) * t557;
t464 = Icges(3,1) * t503 - Icges(3,4) * t502 + Icges(3,5) * t559;
t463 = Icges(3,1) * t501 - Icges(3,4) * t500 - Icges(3,5) * t557;
t462 = Icges(3,4) * t503 - Icges(3,2) * t502 + Icges(3,6) * t559;
t461 = Icges(3,4) * t501 - Icges(3,2) * t500 - Icges(3,6) * t557;
t458 = rSges(4,1) * t499 - rSges(4,2) * t498 - rSges(4,3) * t558;
t457 = Icges(4,1) * t499 - Icges(4,4) * t498 - Icges(4,5) * t558;
t456 = Icges(4,4) * t499 - Icges(4,2) * t498 - Icges(4,6) * t558;
t455 = Icges(4,5) * t499 - Icges(4,6) * t498 - Icges(4,3) * t558;
t453 = t487 * t532 + t560;
t452 = -t487 * t528 + t502 * t532;
t451 = t485 * t532 + t561;
t450 = -t485 * t528 + t500 * t532;
t449 = t487 * t522 + t502 * t521;
t448 = -t487 * t521 + t502 * t522;
t447 = t485 * t522 + t500 * t521;
t446 = -t485 * t521 + t500 * t522;
t443 = t487 * t518 + t502 * t517;
t442 = -t487 * t517 + t502 * t518;
t441 = t485 * t518 + t500 * t517;
t440 = -t485 * t517 + t500 * t518;
t436 = qJD(6) * t498 + t454;
t434 = rSges(4,1) * t487 - rSges(4,2) * t486 + rSges(4,3) * t502;
t433 = rSges(4,1) * t485 - rSges(4,2) * t484 + rSges(4,3) * t500;
t432 = Icges(4,1) * t487 - Icges(4,4) * t486 + Icges(4,5) * t502;
t431 = Icges(4,1) * t485 - Icges(4,4) * t484 + Icges(4,5) * t500;
t430 = Icges(4,4) * t487 - Icges(4,2) * t486 + Icges(4,6) * t502;
t429 = Icges(4,4) * t485 - Icges(4,2) * t484 + Icges(4,6) * t500;
t428 = Icges(4,5) * t487 - Icges(4,6) * t486 + Icges(4,3) * t502;
t427 = Icges(4,5) * t485 - Icges(4,6) * t484 + Icges(4,3) * t500;
t426 = rSges(5,1) * t483 + rSges(5,2) * t482 + rSges(5,3) * t498;
t425 = Icges(5,1) * t483 + Icges(5,4) * t482 + Icges(5,5) * t498;
t424 = Icges(5,4) * t483 + Icges(5,2) * t482 + Icges(5,6) * t498;
t423 = Icges(5,5) * t483 + Icges(5,6) * t482 + Icges(5,3) * t498;
t421 = rSges(6,1) * t477 + rSges(6,2) * t476 + rSges(6,3) * t498;
t417 = Icges(6,1) * t477 + Icges(6,4) * t476 + Icges(6,5) * t498;
t416 = Icges(6,4) * t477 + Icges(6,2) * t476 + Icges(6,6) * t498;
t415 = Icges(6,5) * t477 + Icges(6,6) * t476 + Icges(6,3) * t498;
t414 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t498;
t413 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t498;
t412 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t498;
t411 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t498;
t410 = t467 * t515 - t493 * t514 + t506;
t409 = -t466 * t515 - t493 * t549 - t553;
t408 = (t466 * t531 + t467 * t534) * t552;
t407 = pkin(12) * t498 + t499 * t555 - t548 * t558;
t406 = qJD(6) * t484 + t420;
t405 = qJD(6) * t486 + t419;
t404 = rSges(5,1) * t453 + rSges(5,2) * t452 + rSges(5,3) * t486;
t403 = rSges(5,1) * t451 + rSges(5,2) * t450 + rSges(5,3) * t484;
t402 = Icges(5,1) * t453 + Icges(5,4) * t452 + Icges(5,5) * t486;
t401 = Icges(5,1) * t451 + Icges(5,4) * t450 + Icges(5,5) * t484;
t400 = Icges(5,4) * t453 + Icges(5,2) * t452 + Icges(5,6) * t486;
t399 = Icges(5,4) * t451 + Icges(5,2) * t450 + Icges(5,6) * t484;
t398 = Icges(5,5) * t453 + Icges(5,6) * t452 + Icges(5,3) * t486;
t397 = Icges(5,5) * t451 + Icges(5,6) * t450 + Icges(5,3) * t484;
t395 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t486;
t394 = rSges(6,1) * t447 + rSges(6,2) * t446 + rSges(6,3) * t484;
t393 = Icges(6,1) * t449 + Icges(6,4) * t448 + Icges(6,5) * t486;
t392 = Icges(6,1) * t447 + Icges(6,4) * t446 + Icges(6,5) * t484;
t391 = Icges(6,4) * t449 + Icges(6,2) * t448 + Icges(6,6) * t486;
t390 = Icges(6,4) * t447 + Icges(6,2) * t446 + Icges(6,6) * t484;
t389 = Icges(6,5) * t449 + Icges(6,6) * t448 + Icges(6,3) * t486;
t388 = Icges(6,5) * t447 + Icges(6,6) * t446 + Icges(6,3) * t484;
t385 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t486;
t384 = rSges(7,1) * t441 + rSges(7,2) * t440 + rSges(7,3) * t484;
t383 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t486;
t382 = Icges(7,1) * t441 + Icges(7,4) * t440 + Icges(7,5) * t484;
t381 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t486;
t380 = Icges(7,4) * t441 + Icges(7,2) * t440 + Icges(7,6) * t484;
t379 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t486;
t378 = Icges(7,5) * t441 + Icges(7,6) * t440 + Icges(7,3) * t484;
t376 = pkin(12) * t486 + t487 * t555 + t502 * t548;
t375 = pkin(12) * t484 + t485 * t555 + t500 * t548;
t373 = t434 * t505 - t458 * t488 + t543;
t372 = -t433 * t505 + t458 * t489 + t541;
t371 = t433 * t488 - t434 * t489 + t556;
t370 = t404 * t475 - t426 * t444 + t540;
t369 = -t403 * t475 + t426 * t445 + t539;
t368 = t403 * t444 - t404 * t445 + t544;
t367 = t395 * t454 - t419 * t421 + t538;
t366 = -t394 * t454 + t420 * t421 + t537;
t365 = t394 * t419 - t395 * t420 + t542;
t364 = t376 * t454 + t385 * t436 - t405 * t414 - t407 * t419 + t538;
t363 = -t375 * t454 - t384 * t436 + t406 * t414 + t407 * t420 + t537;
t362 = t375 * t419 - t376 * t420 + t384 * t405 - t385 * t406 + t542;
t1 = m(5) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(6) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(7) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(3) * (t408 ^ 2 + t409 ^ 2 + t410 ^ 2) / 0.2e1 + t515 * ((t562 * t460 + (t462 * t533 + t464 * t530) * t527) * t514 - (t562 * t459 + (t461 * t533 + t463 * t530) * t527) * t549 + (t562 * t490 + (t491 * t533 + t492 * t530) * t527) * t515) / 0.2e1 + t488 * ((t428 * t502 - t430 * t486 + t432 * t487) * t488 + (t427 * t502 - t429 * t486 + t431 * t487) * t489 + (t455 * t502 - t456 * t486 + t457 * t487) * t505) / 0.2e1 + t489 * ((t428 * t500 - t430 * t484 + t432 * t485) * t488 + (t427 * t500 - t429 * t484 + t431 * t485) * t489 + (t455 * t500 - t456 * t484 + t457 * t485) * t505) / 0.2e1 + t505 * ((-t428 * t558 - t430 * t498 + t432 * t499) * t488 + (-t427 * t558 - t429 * t498 + t431 * t499) * t489 + (-t455 * t558 - t456 * t498 + t457 * t499) * t505) / 0.2e1 + t444 * ((t398 * t486 + t400 * t452 + t402 * t453) * t444 + (t397 * t486 + t399 * t452 + t401 * t453) * t445 + (t423 * t486 + t424 * t452 + t425 * t453) * t475) / 0.2e1 + t445 * ((t398 * t484 + t400 * t450 + t402 * t451) * t444 + (t397 * t484 + t399 * t450 + t401 * t451) * t445 + (t423 * t484 + t424 * t450 + t425 * t451) * t475) / 0.2e1 + t475 * ((t398 * t498 + t400 * t482 + t402 * t483) * t444 + (t397 * t498 + t399 * t482 + t401 * t483) * t445 + (t423 * t498 + t424 * t482 + t425 * t483) * t475) / 0.2e1 + t419 * ((t486 * t389 + t448 * t391 + t449 * t393) * t419 + (t388 * t486 + t390 * t448 + t392 * t449) * t420 + (t415 * t486 + t416 * t448 + t417 * t449) * t454) / 0.2e1 + t420 * ((t389 * t484 + t391 * t446 + t393 * t447) * t419 + (t484 * t388 + t446 * t390 + t447 * t392) * t420 + (t415 * t484 + t416 * t446 + t417 * t447) * t454) / 0.2e1 + t454 * ((t389 * t498 + t391 * t476 + t393 * t477) * t419 + (t388 * t498 + t390 * t476 + t392 * t477) * t420 + (t415 * t498 + t416 * t476 + t417 * t477) * t454) / 0.2e1 + t405 * ((t486 * t379 + t442 * t381 + t443 * t383) * t405 + (t378 * t486 + t380 * t442 + t382 * t443) * t406 + (t411 * t486 + t412 * t442 + t413 * t443) * t436) / 0.2e1 + t406 * ((t379 * t484 + t381 * t440 + t383 * t441) * t405 + (t484 * t378 + t440 * t380 + t441 * t382) * t406 + (t411 * t484 + t412 * t440 + t413 * t441) * t436) / 0.2e1 + t436 * ((t379 * t498 + t381 * t471 + t383 * t472) * t405 + (t378 * t498 + t380 * t471 + t382 * t472) * t406 + (t498 * t411 + t471 * t412 + t472 * t413) * t436) / 0.2e1 + m(4) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 - ((-t490 * t557 - t491 * t500 + t492 * t501) * t515 + ((-t462 * t500 + t464 * t501) * t531 + (t500 * t461 - t501 * t463 + t568) * t534) * t552) * t549 / 0.2e1 + ((t490 * t559 - t491 * t502 + t492 * t503) * t515 + (-(-t461 * t502 + t463 * t503) * t534 + (-t502 * t462 + t503 * t464 - t568) * t531) * t552) * t514 / 0.2e1 + (m(2) * (t510 ^ 2 + t511 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
