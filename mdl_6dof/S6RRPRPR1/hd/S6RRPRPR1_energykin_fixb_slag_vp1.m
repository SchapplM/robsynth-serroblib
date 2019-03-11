% Calculate kinetic energy for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:52
% EndTime: 2019-03-09 10:07:54
% DurationCPUTime: 2.43s
% Computational Cost: add. (1988->310), mult. (1871->476), div. (0->0), fcn. (1772->12), ass. (0->168)
t555 = Icges(3,3) + Icges(4,3);
t469 = qJ(2) + pkin(10);
t459 = sin(t469);
t461 = cos(t469);
t474 = sin(qJ(2));
t476 = cos(qJ(2));
t554 = Icges(3,5) * t476 + Icges(4,5) * t461 - Icges(3,6) * t474 - Icges(4,6) * t459;
t475 = sin(qJ(1));
t477 = cos(qJ(1));
t553 = t554 * t475 - t477 * t555;
t552 = t475 * t555 + t554 * t477;
t551 = Icges(3,5) * t474 + Icges(4,5) * t459 + Icges(3,6) * t476 + Icges(4,6) * t461;
t535 = Icges(4,4) * t459;
t437 = Icges(4,2) * t461 + t535;
t534 = Icges(4,4) * t461;
t438 = Icges(4,1) * t459 + t534;
t537 = Icges(3,4) * t474;
t448 = Icges(3,2) * t476 + t537;
t536 = Icges(3,4) * t476;
t449 = Icges(3,1) * t474 + t536;
t550 = -t437 * t459 + t438 * t461 - t448 * t474 + t449 * t476;
t498 = -Icges(4,2) * t459 + t534;
t405 = Icges(4,6) * t475 + t477 * t498;
t501 = Icges(4,1) * t461 - t535;
t407 = Icges(4,5) * t475 + t477 * t501;
t499 = -Icges(3,2) * t474 + t536;
t422 = Icges(3,6) * t475 + t477 * t499;
t502 = Icges(3,1) * t476 - t537;
t424 = Icges(3,5) * t475 + t477 * t502;
t549 = -t405 * t459 + t407 * t461 - t422 * t474 + t424 * t476;
t404 = -Icges(4,6) * t477 + t475 * t498;
t406 = -Icges(4,5) * t477 + t475 * t501;
t421 = -Icges(3,6) * t477 + t475 * t499;
t423 = -Icges(3,5) * t477 + t475 * t502;
t548 = t404 * t459 - t406 * t461 + t421 * t474 - t423 * t476;
t542 = pkin(2) * t474;
t540 = t476 * pkin(2);
t471 = cos(pkin(11));
t539 = pkin(5) * t471;
t462 = qJ(4) + t469;
t454 = sin(t462);
t533 = Icges(5,4) * t454;
t455 = cos(t462);
t532 = Icges(5,4) * t455;
t531 = t454 * t475;
t530 = t454 * t477;
t468 = pkin(11) + qJ(6);
t458 = sin(t468);
t529 = t458 * t475;
t528 = t458 * t477;
t460 = cos(t468);
t527 = t460 * t475;
t526 = t460 * t477;
t470 = sin(pkin(11));
t525 = t470 * t475;
t524 = t470 * t477;
t523 = t471 * t475;
t522 = t471 * t477;
t400 = -qJ(3) * t477 + t475 * t540;
t401 = qJ(3) * t475 + t477 * t540;
t465 = qJD(2) * t475;
t516 = qJD(2) * t477;
t520 = t400 * t465 + t401 * t516;
t453 = pkin(1) * t475 - pkin(7) * t477;
t519 = -t400 - t453;
t518 = pkin(3) * t461;
t445 = qJD(4) * t475 + t465;
t515 = qJD(5) * t454;
t514 = qJD(6) * t454;
t377 = -pkin(8) * t477 + t475 * t518;
t513 = -t377 + t519;
t446 = (-qJD(2) - qJD(4)) * t477;
t378 = pkin(8) * t475 + t477 * t518;
t510 = t377 * t465 + t378 * t516 + t520;
t504 = pkin(4) * t455 + qJ(5) * t454;
t415 = t504 * t475;
t509 = -t415 + t513;
t441 = qJD(1) * (pkin(1) * t477 + pkin(7) * t475);
t508 = qJD(1) * t401 - qJD(3) * t477 + t441;
t507 = rSges(3,1) * t476 - rSges(3,2) * t474;
t506 = rSges(4,1) * t461 - rSges(4,2) * t459;
t505 = rSges(5,1) * t455 - rSges(5,2) * t454;
t503 = qJD(2) * (-rSges(4,1) * t459 - rSges(4,2) * t461 - t542);
t500 = Icges(5,1) * t455 - t533;
t497 = -Icges(5,2) * t454 + t532;
t494 = Icges(5,5) * t455 - Icges(5,6) * t454;
t487 = qJD(2) * (-pkin(3) * t459 - t542);
t486 = -qJD(5) * t455 + t445 * t415 + t510;
t485 = (Icges(5,5) * t454 + Icges(5,6) * t455) * qJD(1) + (-Icges(5,3) * t477 + t475 * t494) * t446 + (Icges(5,3) * t475 + t477 * t494) * t445;
t484 = pkin(9) * t454 + t455 * t539;
t464 = qJD(3) * t475;
t483 = t477 * t487 + t464;
t434 = pkin(4) * t454 - qJ(5) * t455;
t482 = t446 * t434 + t477 * t515 + t483;
t481 = qJD(1) * t378 + t475 * t487 + t508;
t416 = t504 * t477;
t480 = qJD(1) * t416 + t475 * t515 + t481;
t393 = -Icges(5,6) * t477 + t475 * t497;
t394 = Icges(5,6) * t475 + t477 * t497;
t395 = -Icges(5,5) * t477 + t475 * t500;
t396 = Icges(5,5) * t475 + t477 * t500;
t432 = Icges(5,2) * t455 + t533;
t433 = Icges(5,1) * t454 + t532;
t479 = (-t394 * t454 + t396 * t455) * t445 + (-t393 * t454 + t395 * t455) * t446 + (-t432 * t454 + t433 * t455) * qJD(1);
t452 = rSges(2,1) * t477 - rSges(2,2) * t475;
t451 = rSges(2,1) * t475 + rSges(2,2) * t477;
t450 = rSges(3,1) * t474 + rSges(3,2) * t476;
t444 = -qJD(6) * t455 + qJD(1);
t435 = rSges(5,1) * t454 + rSges(5,2) * t455;
t430 = rSges(3,3) * t475 + t477 * t507;
t429 = -rSges(3,3) * t477 + t475 * t507;
t428 = t455 * t522 + t525;
t427 = -t455 * t524 + t523;
t426 = t455 * t523 - t524;
t425 = -t455 * t525 - t522;
t418 = t475 * t514 + t446;
t417 = t477 * t514 + t445;
t414 = t455 * t526 + t529;
t413 = -t455 * t528 + t527;
t412 = t455 * t527 - t528;
t411 = -t455 * t529 - t526;
t410 = rSges(4,3) * t475 + t477 * t506;
t409 = -rSges(4,3) * t477 + t475 * t506;
t399 = rSges(5,3) * t475 + t477 * t505;
t398 = -rSges(5,3) * t477 + t475 * t505;
t387 = -rSges(6,3) * t455 + (rSges(6,1) * t471 - rSges(6,2) * t470) * t454;
t386 = -Icges(6,5) * t455 + (Icges(6,1) * t471 - Icges(6,4) * t470) * t454;
t385 = -Icges(6,6) * t455 + (Icges(6,4) * t471 - Icges(6,2) * t470) * t454;
t384 = -Icges(6,3) * t455 + (Icges(6,5) * t471 - Icges(6,6) * t470) * t454;
t382 = -rSges(7,3) * t455 + (rSges(7,1) * t460 - rSges(7,2) * t458) * t454;
t381 = -Icges(7,5) * t455 + (Icges(7,1) * t460 - Icges(7,4) * t458) * t454;
t380 = -Icges(7,6) * t455 + (Icges(7,4) * t460 - Icges(7,2) * t458) * t454;
t379 = -Icges(7,3) * t455 + (Icges(7,5) * t460 - Icges(7,6) * t458) * t454;
t376 = -pkin(9) * t455 + t454 * t539;
t372 = qJD(1) * t430 - t450 * t465 + t441;
t371 = -t450 * t516 + (-t429 - t453) * qJD(1);
t370 = (t429 * t475 + t430 * t477) * qJD(2);
t369 = rSges(6,1) * t428 + rSges(6,2) * t427 + rSges(6,3) * t530;
t368 = rSges(6,1) * t426 + rSges(6,2) * t425 + rSges(6,3) * t531;
t367 = Icges(6,1) * t428 + Icges(6,4) * t427 + Icges(6,5) * t530;
t366 = Icges(6,1) * t426 + Icges(6,4) * t425 + Icges(6,5) * t531;
t365 = Icges(6,4) * t428 + Icges(6,2) * t427 + Icges(6,6) * t530;
t364 = Icges(6,4) * t426 + Icges(6,2) * t425 + Icges(6,6) * t531;
t363 = Icges(6,5) * t428 + Icges(6,6) * t427 + Icges(6,3) * t530;
t362 = Icges(6,5) * t426 + Icges(6,6) * t425 + Icges(6,3) * t531;
t361 = pkin(5) * t525 + t477 * t484;
t360 = -pkin(5) * t524 + t475 * t484;
t359 = rSges(7,1) * t414 + rSges(7,2) * t413 + rSges(7,3) * t530;
t358 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t531;
t357 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t530;
t356 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t531;
t355 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t530;
t354 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t531;
t353 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t530;
t352 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t531;
t351 = qJD(1) * t410 + t475 * t503 + t508;
t350 = t464 + t477 * t503 + (-t409 + t519) * qJD(1);
t349 = (t409 * t475 + t410 * t477) * qJD(2) + t520;
t348 = qJD(1) * t399 - t435 * t445 + t481;
t347 = t435 * t446 + (-t398 + t513) * qJD(1) + t483;
t346 = t398 * t445 - t399 * t446 + t510;
t345 = qJD(1) * t369 + (-t387 - t434) * t445 + t480;
t344 = t387 * t446 + (-t368 + t509) * qJD(1) + t482;
t343 = t368 * t445 + (-t369 - t416) * t446 + t486;
t342 = qJD(1) * t361 + t359 * t444 - t382 * t417 + (-t376 - t434) * t445 + t480;
t341 = -t358 * t444 + t376 * t446 + t382 * t418 + (-t360 + t509) * qJD(1) + t482;
t340 = t358 * t417 - t359 * t418 + t360 * t445 + (-t361 - t416) * t446 + t486;
t1 = t417 * ((t353 * t530 + t413 * t355 + t414 * t357) * t417 + (t352 * t530 + t354 * t413 + t356 * t414) * t418 + (t379 * t530 + t380 * t413 + t381 * t414) * t444) / 0.2e1 + t418 * ((t353 * t531 + t355 * t411 + t357 * t412) * t417 + (t352 * t531 + t411 * t354 + t412 * t356) * t418 + (t379 * t531 + t380 * t411 + t381 * t412) * t444) / 0.2e1 + t444 * ((-t352 * t418 - t353 * t417 - t379 * t444) * t455 + ((-t355 * t458 + t357 * t460) * t417 + (-t354 * t458 + t356 * t460) * t418 + (-t380 * t458 + t381 * t460) * t444) * t454) / 0.2e1 + m(3) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(7) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(6) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(5) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(4) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + ((t363 * t530 + t427 * t365 + t428 * t367) * t445 + (t362 * t530 + t364 * t427 + t366 * t428) * t446 + (t384 * t530 + t385 * t427 + t428 * t386) * qJD(1) + t485 * t475 + t479 * t477) * t445 / 0.2e1 + (t479 * t475 - t485 * t477 + (t363 * t531 + t365 * t425 + t367 * t426) * t445 + (t362 * t531 + t425 * t364 + t426 * t366) * t446 + (t384 * t531 + t385 * t425 + t386 * t426) * qJD(1)) * t446 / 0.2e1 + (m(2) * (t451 ^ 2 + t452 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t552 * t475 ^ 2 + (t548 * t477 + (t549 - t553) * t475) * t477) * qJD(2) + (t551 * t475 + t550 * t477) * qJD(1)) * t465 / 0.2e1 - ((t553 * t477 ^ 2 + (t549 * t475 + (t548 - t552) * t477) * t475) * qJD(2) + (t550 * t475 - t551 * t477) * qJD(1)) * t516 / 0.2e1 + ((t394 * t455 + t396 * t454) * t445 + (t393 * t455 + t395 * t454) * t446 + (-t362 * t446 - t363 * t445) * t455 + ((-t365 * t470 + t367 * t471) * t445 + (-t364 * t470 + t366 * t471) * t446) * t454 + ((-t404 * t461 - t406 * t459 - t421 * t476 - t423 * t474) * t477 + (t405 * t461 + t407 * t459 + t422 * t476 + t424 * t474) * t475) * qJD(2) + (t461 * t437 + t459 * t438 + t476 * t448 + t474 * t449 + (t432 - t384) * t455 + (-t385 * t470 + t386 * t471 + t433) * t454) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
