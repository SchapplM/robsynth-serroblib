% Calculate kinetic energy for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:11
% EndTime: 2019-03-09 15:21:13
% DurationCPUTime: 2.62s
% Computational Cost: add. (2010->312), mult. (1893->467), div. (0->0), fcn. (1794->12), ass. (0->170)
t548 = Icges(4,3) + Icges(5,3);
t469 = qJ(2) + qJ(3);
t460 = pkin(10) + t469;
t454 = sin(t460);
t455 = cos(t460);
t464 = sin(t469);
t465 = cos(t469);
t547 = -Icges(4,5) * t465 - Icges(5,5) * t455 + Icges(4,6) * t464 + Icges(5,6) * t454;
t474 = sin(qJ(1));
t476 = cos(qJ(1));
t531 = Icges(5,4) * t455;
t496 = -Icges(5,2) * t454 + t531;
t392 = -Icges(5,6) * t476 + t474 * t496;
t393 = Icges(5,6) * t474 + t476 * t496;
t532 = Icges(5,4) * t454;
t499 = Icges(5,1) * t455 - t532;
t394 = -Icges(5,5) * t476 + t474 * t499;
t395 = Icges(5,5) * t474 + t476 * t499;
t533 = Icges(4,4) * t465;
t497 = -Icges(4,2) * t464 + t533;
t404 = -Icges(4,6) * t476 + t474 * t497;
t405 = Icges(4,6) * t474 + t476 * t497;
t534 = Icges(4,4) * t464;
t500 = Icges(4,1) * t465 - t534;
t406 = -Icges(4,5) * t476 + t474 * t500;
t407 = Icges(4,5) * t474 + t476 * t500;
t431 = Icges(5,2) * t455 + t532;
t432 = Icges(5,1) * t454 + t531;
t437 = Icges(4,2) * t465 + t534;
t438 = Icges(4,1) * t464 + t533;
t463 = qJD(2) * t474;
t445 = qJD(3) * t474 + t463;
t446 = (-qJD(2) - qJD(3)) * t476;
t546 = (-t392 * t454 + t394 * t455 - t404 * t464 + t406 * t465) * t446 + (-t393 * t454 + t395 * t455 - t405 * t464 + t407 * t465) * t445 + (-t431 * t454 + t432 * t455 - t437 * t464 + t438 * t465) * qJD(1);
t545 = (t547 * t474 + t548 * t476) * t446 + (-t548 * t474 + t547 * t476) * t445 + (-Icges(4,5) * t464 - Icges(5,5) * t454 - Icges(4,6) * t465 - Icges(5,6) * t455) * qJD(1);
t541 = pkin(3) * t464;
t475 = cos(qJ(2));
t539 = t475 * pkin(2);
t471 = cos(pkin(11));
t538 = pkin(5) * t471;
t473 = sin(qJ(2));
t536 = Icges(3,4) * t473;
t535 = Icges(3,4) * t475;
t530 = t454 * t474;
t529 = t454 * t476;
t468 = pkin(11) + qJ(6);
t458 = sin(t468);
t528 = t458 * t474;
t527 = t458 * t476;
t459 = cos(t468);
t526 = t459 * t474;
t525 = t459 * t476;
t470 = sin(pkin(11));
t524 = t470 * t474;
t523 = t470 * t476;
t522 = t471 * t474;
t521 = t471 * t476;
t516 = pkin(3) * t465;
t377 = qJ(4) * t474 + t476 * t516;
t502 = pkin(4) * t455 + qJ(5) * t454;
t415 = t502 * t476;
t519 = -t377 - t415;
t399 = -pkin(8) * t476 + t474 * t539;
t400 = pkin(8) * t474 + t476 * t539;
t514 = qJD(2) * t476;
t518 = t399 * t463 + t400 * t514;
t453 = pkin(1) * t474 - pkin(7) * t476;
t517 = -t399 - t453;
t513 = qJD(5) * t454;
t512 = qJD(6) * t454;
t511 = pkin(2) * qJD(2) * t473;
t376 = -qJ(4) * t476 + t474 * t516;
t510 = t445 * t376 + t518;
t509 = -t376 + t517;
t433 = pkin(4) * t454 - qJ(5) * t455;
t508 = -t433 - t541;
t507 = t476 * t511;
t414 = t502 * t474;
t506 = -t414 + t509;
t505 = rSges(3,1) * t475 - rSges(3,2) * t473;
t504 = rSges(4,1) * t465 - rSges(4,2) * t464;
t503 = rSges(5,1) * t455 - rSges(5,2) * t454;
t501 = Icges(3,1) * t475 - t536;
t498 = -Icges(3,2) * t473 + t535;
t495 = Icges(3,5) * t475 - Icges(3,6) * t473;
t420 = -Icges(3,6) * t476 + t474 * t498;
t422 = -Icges(3,5) * t476 + t474 * t501;
t492 = t420 * t473 - t422 * t475;
t421 = Icges(3,6) * t474 + t476 * t498;
t423 = Icges(3,5) * t474 + t476 * t501;
t491 = -t421 * t473 + t423 * t475;
t448 = Icges(3,2) * t475 + t536;
t449 = Icges(3,1) * t473 + t535;
t490 = -t448 * t473 + t449 * t475;
t441 = qJD(1) * (pkin(1) * t476 + pkin(7) * t474);
t489 = qJD(1) * t400 - t474 * t511 + t441;
t488 = qJD(4) * t474 + t446 * t541 - t507;
t487 = -qJD(5) * t455 + t445 * t414 + t510;
t484 = pkin(9) * t454 + t455 * t538;
t483 = t446 * t433 + t476 * t513 + t488;
t482 = qJD(1) * t377 - qJD(4) * t476 + t489;
t481 = qJD(1) * t415 + t474 * t513 + t482;
t452 = rSges(2,1) * t476 - rSges(2,2) * t474;
t451 = rSges(2,1) * t474 + rSges(2,2) * t476;
t450 = rSges(3,1) * t473 + rSges(3,2) * t475;
t447 = Icges(3,5) * t473 + Icges(3,6) * t475;
t444 = -qJD(6) * t455 + qJD(1);
t439 = rSges(4,1) * t464 + rSges(4,2) * t465;
t434 = rSges(5,1) * t454 + rSges(5,2) * t455;
t429 = rSges(3,3) * t474 + t476 * t505;
t428 = -rSges(3,3) * t476 + t474 * t505;
t427 = t455 * t521 + t524;
t426 = -t455 * t523 + t522;
t425 = t455 * t522 - t523;
t424 = -t455 * t524 - t521;
t419 = Icges(3,3) * t474 + t476 * t495;
t418 = -Icges(3,3) * t476 + t474 * t495;
t417 = t474 * t512 + t446;
t416 = t476 * t512 + t445;
t413 = rSges(4,3) * t474 + t476 * t504;
t412 = -rSges(4,3) * t476 + t474 * t504;
t411 = t455 * t525 + t528;
t410 = -t455 * t527 + t526;
t409 = t455 * t526 - t527;
t408 = -t455 * t528 - t525;
t397 = rSges(5,3) * t474 + t476 * t503;
t396 = -rSges(5,3) * t476 + t474 * t503;
t386 = -rSges(6,3) * t455 + (rSges(6,1) * t471 - rSges(6,2) * t470) * t454;
t385 = -Icges(6,5) * t455 + (Icges(6,1) * t471 - Icges(6,4) * t470) * t454;
t384 = -Icges(6,6) * t455 + (Icges(6,4) * t471 - Icges(6,2) * t470) * t454;
t383 = -Icges(6,3) * t455 + (Icges(6,5) * t471 - Icges(6,6) * t470) * t454;
t381 = -rSges(7,3) * t455 + (rSges(7,1) * t459 - rSges(7,2) * t458) * t454;
t380 = -Icges(7,5) * t455 + (Icges(7,1) * t459 - Icges(7,4) * t458) * t454;
t379 = -Icges(7,6) * t455 + (Icges(7,4) * t459 - Icges(7,2) * t458) * t454;
t378 = -Icges(7,3) * t455 + (Icges(7,5) * t459 - Icges(7,6) * t458) * t454;
t375 = -pkin(9) * t455 + t454 * t538;
t373 = qJD(1) * t429 - t450 * t463 + t441;
t372 = -t450 * t514 + (-t428 - t453) * qJD(1);
t371 = (t428 * t474 + t429 * t476) * qJD(2);
t369 = rSges(6,1) * t427 + rSges(6,2) * t426 + rSges(6,3) * t529;
t368 = rSges(6,1) * t425 + rSges(6,2) * t424 + rSges(6,3) * t530;
t367 = Icges(6,1) * t427 + Icges(6,4) * t426 + Icges(6,5) * t529;
t366 = Icges(6,1) * t425 + Icges(6,4) * t424 + Icges(6,5) * t530;
t365 = Icges(6,4) * t427 + Icges(6,2) * t426 + Icges(6,6) * t529;
t364 = Icges(6,4) * t425 + Icges(6,2) * t424 + Icges(6,6) * t530;
t363 = Icges(6,5) * t427 + Icges(6,6) * t426 + Icges(6,3) * t529;
t362 = Icges(6,5) * t425 + Icges(6,6) * t424 + Icges(6,3) * t530;
t361 = pkin(5) * t524 + t476 * t484;
t360 = -pkin(5) * t523 + t474 * t484;
t359 = rSges(7,1) * t411 + rSges(7,2) * t410 + rSges(7,3) * t529;
t358 = rSges(7,1) * t409 + rSges(7,2) * t408 + rSges(7,3) * t530;
t357 = Icges(7,1) * t411 + Icges(7,4) * t410 + Icges(7,5) * t529;
t356 = Icges(7,1) * t409 + Icges(7,4) * t408 + Icges(7,5) * t530;
t355 = Icges(7,4) * t411 + Icges(7,2) * t410 + Icges(7,6) * t529;
t354 = Icges(7,4) * t409 + Icges(7,2) * t408 + Icges(7,6) * t530;
t353 = Icges(7,5) * t411 + Icges(7,6) * t410 + Icges(7,3) * t529;
t352 = Icges(7,5) * t409 + Icges(7,6) * t408 + Icges(7,3) * t530;
t351 = qJD(1) * t413 - t439 * t445 + t489;
t350 = -t507 + t439 * t446 + (-t412 + t517) * qJD(1);
t349 = t412 * t445 - t413 * t446 + t518;
t348 = qJD(1) * t397 + (-t434 - t541) * t445 + t482;
t347 = t434 * t446 + (-t396 + t509) * qJD(1) + t488;
t346 = t396 * t445 + (-t377 - t397) * t446 + t510;
t345 = qJD(1) * t369 + (-t386 + t508) * t445 + t481;
t344 = t386 * t446 + (-t368 + t506) * qJD(1) + t483;
t343 = t368 * t445 + (-t369 + t519) * t446 + t487;
t342 = qJD(1) * t361 + t359 * t444 - t381 * t416 + (-t375 + t508) * t445 + t481;
t341 = -t358 * t444 + t375 * t446 + t381 * t417 + (-t360 + t506) * qJD(1) + t483;
t340 = t358 * t416 - t359 * t417 + t360 * t445 + (-t361 + t519) * t446 + t487;
t1 = m(7) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(5) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(6) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(4) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(3) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + t416 * ((t353 * t529 + t410 * t355 + t411 * t357) * t416 + (t352 * t529 + t354 * t410 + t356 * t411) * t417 + (t378 * t529 + t379 * t410 + t380 * t411) * t444) / 0.2e1 + t417 * ((t353 * t530 + t355 * t408 + t357 * t409) * t416 + (t352 * t530 + t408 * t354 + t409 * t356) * t417 + (t378 * t530 + t379 * t408 + t380 * t409) * t444) / 0.2e1 + t444 * ((-t352 * t417 - t353 * t416 - t378 * t444) * t455 + ((-t355 * t458 + t357 * t459) * t416 + (-t354 * t458 + t356 * t459) * t417 + (-t379 * t458 + t380 * t459) * t444) * t454) / 0.2e1 - ((-t476 * t447 + t474 * t490) * qJD(1) + (t476 ^ 2 * t418 + (t491 * t474 + (-t419 + t492) * t476) * t474) * qJD(2)) * t514 / 0.2e1 + ((t474 * t447 + t476 * t490) * qJD(1) + (t474 ^ 2 * t419 + (t492 * t476 + (-t418 + t491) * t474) * t476) * qJD(2)) * t463 / 0.2e1 + (m(2) * (t451 ^ 2 + t452 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t363 * t529 + t426 * t365 + t427 * t367) * t445 + (t362 * t529 + t364 * t426 + t366 * t427) * t446 + (t383 * t529 + t384 * t426 + t385 * t427) * qJD(1) + t546 * t476 - t545 * t474) * t445 / 0.2e1 + ((t363 * t530 + t365 * t424 + t367 * t425) * t445 + (t362 * t530 + t424 * t364 + t425 * t366) * t446 + (t383 * t530 + t384 * t424 + t385 * t425) * qJD(1) + t545 * t476 + t546 * t474) * t446 / 0.2e1 + (((t421 * t475 + t423 * t473) * t474 - (t420 * t475 + t422 * t473) * t476) * qJD(2) + (t404 * t465 + t406 * t464 + (-t362 + t392) * t455 + (-t364 * t470 + t366 * t471 + t394) * t454) * t446 + (t405 * t465 + t407 * t464 + (-t363 + t393) * t455 + (-t365 * t470 + t367 * t471 + t395) * t454) * t445 + (t465 * t437 + t464 * t438 + t448 * t475 + t473 * t449 + (-t383 + t431) * t455 + (-t384 * t470 + t385 * t471 + t432) * t454) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
