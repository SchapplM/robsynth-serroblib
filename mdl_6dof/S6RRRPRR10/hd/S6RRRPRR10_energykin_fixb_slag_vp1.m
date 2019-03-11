% Calculate kinetic energy for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:36
% EndTime: 2019-03-09 19:15:39
% DurationCPUTime: 3.17s
% Computational Cost: add. (1549->324), mult. (3281->502), div. (0->0), fcn. (3624->10), ass. (0->159)
t530 = Icges(4,1) + Icges(5,1);
t529 = -Icges(4,4) + Icges(5,5);
t528 = Icges(5,4) + Icges(4,5);
t527 = Icges(4,2) + Icges(5,3);
t526 = -Icges(5,6) + Icges(4,6);
t525 = -Icges(4,3) - Icges(5,2);
t470 = sin(qJ(3));
t474 = cos(qJ(3));
t476 = cos(qJ(1));
t472 = sin(qJ(1));
t475 = cos(qJ(2));
t502 = t472 * t475;
t438 = t470 * t502 + t474 * t476;
t439 = -t470 * t476 + t474 * t502;
t471 = sin(qJ(2));
t504 = t471 * t472;
t524 = t527 * t438 + t529 * t439 - t526 * t504;
t501 = t475 * t476;
t440 = t470 * t501 - t472 * t474;
t441 = t470 * t472 + t474 * t501;
t503 = t471 * t476;
t523 = t527 * t440 + t529 * t441 - t526 * t503;
t522 = -t526 * t438 + t528 * t439 - t525 * t504;
t521 = -t526 * t440 + t528 * t441 - t525 * t503;
t520 = t529 * t438 + t530 * t439 + t528 * t504;
t519 = t529 * t440 + t530 * t441 + t528 * t503;
t518 = t526 * t475 + (t527 * t470 + t529 * t474) * t471;
t517 = t525 * t475 + (-t526 * t470 + t528 * t474) * t471;
t516 = -t528 * t475 + (t529 * t470 + t530 * t474) * t471;
t473 = cos(qJ(5));
t511 = pkin(5) * t473;
t509 = Icges(3,4) * t471;
t508 = Icges(3,4) * t475;
t469 = sin(qJ(5));
t507 = t438 * t469;
t506 = t440 * t469;
t505 = t469 * t470;
t493 = pkin(2) * t475 + pkin(8) * t471;
t443 = t493 * t472;
t444 = t493 * t476;
t465 = qJD(2) * t472;
t499 = qJD(2) * t476;
t500 = t443 * t465 + t444 * t499;
t498 = qJD(3) * t471;
t445 = t476 * t498 + t465;
t497 = qJD(5) * t471;
t496 = t471 * pkin(10);
t446 = t472 * t498 - t499;
t495 = t471 * (-qJD(5) - qJD(6));
t461 = -qJD(3) * t475 + qJD(1);
t401 = pkin(3) * t439 + qJ(4) * t438;
t494 = qJD(4) * t471 * t470 + t445 * t401 + t500;
t492 = rSges(3,1) * t475 - rSges(3,2) * t471;
t491 = Icges(3,1) * t475 - t509;
t490 = -Icges(3,2) * t471 + t508;
t489 = Icges(3,5) * t475 - Icges(3,6) * t471;
t420 = -Icges(3,6) * t476 + t472 * t490;
t424 = -Icges(3,5) * t476 + t472 * t491;
t488 = t420 * t471 - t424 * t475;
t421 = Icges(3,6) * t472 + t476 * t490;
t425 = Icges(3,5) * t472 + t476 * t491;
t487 = -t421 * t471 + t425 * t475;
t451 = Icges(3,2) * t475 + t509;
t452 = Icges(3,1) * t471 + t508;
t486 = -t451 * t471 + t452 * t475;
t449 = qJD(1) * (pkin(1) * t476 + pkin(7) * t472);
t456 = pkin(2) * t471 - pkin(8) * t475;
t485 = qJD(1) * t444 - t456 * t465 + t449;
t402 = pkin(3) * t441 + qJ(4) * t440;
t484 = qJD(4) * t438 + t461 * t402 + t485;
t457 = pkin(1) * t472 - pkin(7) * t476;
t483 = (-t443 - t457) * qJD(1) - t456 * t499;
t408 = pkin(4) * t439 - pkin(9) * t504;
t409 = pkin(4) * t441 - pkin(9) * t503;
t482 = t445 * t408 + (-t402 - t409) * t446 + t494;
t442 = (pkin(3) * t474 + qJ(4) * t470) * t471;
t481 = qJD(4) * t440 + t446 * t442 + t483;
t447 = pkin(4) * t471 * t474 + pkin(9) * t475;
t480 = t461 * t409 + (-t442 - t447) * t445 + t484;
t479 = t446 * t447 + (-t401 - t408) * t461 + t481;
t468 = qJ(5) + qJ(6);
t467 = cos(t468);
t466 = sin(t468);
t464 = qJD(5) * t475;
t455 = rSges(2,1) * t476 - rSges(2,2) * t472;
t454 = rSges(2,1) * t472 + rSges(2,2) * t476;
t453 = rSges(3,1) * t471 + rSges(3,2) * t475;
t450 = Icges(3,5) * t471 + Icges(3,6) * t475;
t448 = t461 + t464;
t437 = qJD(1) + t464 + (-qJD(3) + qJD(6)) * t475;
t433 = (t473 * t474 + t505) * t471;
t432 = (-t469 * t474 + t470 * t473) * t471;
t429 = rSges(3,3) * t472 + t476 * t492;
t428 = -rSges(3,3) * t476 + t472 * t492;
t427 = -rSges(4,3) * t475 + (rSges(4,1) * t474 - rSges(4,2) * t470) * t471;
t426 = -rSges(5,2) * t475 + (rSges(5,1) * t474 + rSges(5,3) * t470) * t471;
t417 = Icges(3,3) * t472 + t476 * t489;
t416 = -Icges(3,3) * t476 + t472 * t489;
t413 = -t472 * t497 + t446;
t412 = -t476 * t497 + t445;
t411 = (t466 * t470 + t467 * t474) * t471;
t410 = (-t466 * t474 + t467 * t470) * t471;
t406 = t472 * t495 + t446;
t405 = t476 * t495 + t445;
t399 = t441 * t473 + t506;
t398 = t440 * t473 - t441 * t469;
t397 = t439 * t473 + t507;
t396 = t438 * t473 - t439 * t469;
t395 = pkin(10) * t475 + (pkin(5) * t505 + t474 * t511) * t471;
t394 = t440 * t466 + t441 * t467;
t393 = t440 * t467 - t441 * t466;
t392 = t438 * t466 + t439 * t467;
t391 = t438 * t467 - t439 * t466;
t390 = rSges(4,1) * t441 - rSges(4,2) * t440 + rSges(4,3) * t503;
t389 = rSges(5,1) * t441 + rSges(5,2) * t503 + rSges(5,3) * t440;
t388 = rSges(4,1) * t439 - rSges(4,2) * t438 + rSges(4,3) * t504;
t387 = rSges(5,1) * t439 + rSges(5,2) * t504 + rSges(5,3) * t438;
t373 = rSges(6,1) * t433 + rSges(6,2) * t432 + rSges(6,3) * t475;
t372 = Icges(6,1) * t433 + Icges(6,4) * t432 + Icges(6,5) * t475;
t371 = Icges(6,4) * t433 + Icges(6,2) * t432 + Icges(6,6) * t475;
t370 = Icges(6,5) * t433 + Icges(6,6) * t432 + Icges(6,3) * t475;
t369 = qJD(1) * t429 - t453 * t465 + t449;
t368 = -t453 * t499 + (-t428 - t457) * qJD(1);
t366 = (t428 * t472 + t429 * t476) * qJD(2);
t365 = rSges(7,1) * t411 + rSges(7,2) * t410 + rSges(7,3) * t475;
t364 = Icges(7,1) * t411 + Icges(7,4) * t410 + Icges(7,5) * t475;
t363 = Icges(7,4) * t411 + Icges(7,2) * t410 + Icges(7,6) * t475;
t362 = Icges(7,5) * t411 + Icges(7,6) * t410 + Icges(7,3) * t475;
t361 = pkin(5) * t506 + t441 * t511 - t476 * t496;
t360 = pkin(5) * t507 + t439 * t511 - t472 * t496;
t359 = rSges(6,1) * t399 + rSges(6,2) * t398 - rSges(6,3) * t503;
t358 = rSges(6,1) * t397 + rSges(6,2) * t396 - rSges(6,3) * t504;
t357 = Icges(6,1) * t399 + Icges(6,4) * t398 - Icges(6,5) * t503;
t356 = Icges(6,1) * t397 + Icges(6,4) * t396 - Icges(6,5) * t504;
t355 = Icges(6,4) * t399 + Icges(6,2) * t398 - Icges(6,6) * t503;
t354 = Icges(6,4) * t397 + Icges(6,2) * t396 - Icges(6,6) * t504;
t353 = Icges(6,5) * t399 + Icges(6,6) * t398 - Icges(6,3) * t503;
t352 = Icges(6,5) * t397 + Icges(6,6) * t396 - Icges(6,3) * t504;
t351 = rSges(7,1) * t394 + rSges(7,2) * t393 - rSges(7,3) * t503;
t350 = rSges(7,1) * t392 + rSges(7,2) * t391 - rSges(7,3) * t504;
t349 = Icges(7,1) * t394 + Icges(7,4) * t393 - Icges(7,5) * t503;
t348 = Icges(7,1) * t392 + Icges(7,4) * t391 - Icges(7,5) * t504;
t347 = Icges(7,4) * t394 + Icges(7,2) * t393 - Icges(7,6) * t503;
t346 = Icges(7,4) * t392 + Icges(7,2) * t391 - Icges(7,6) * t504;
t345 = Icges(7,5) * t394 + Icges(7,6) * t393 - Icges(7,3) * t503;
t344 = Icges(7,5) * t392 + Icges(7,6) * t391 - Icges(7,3) * t504;
t343 = t390 * t461 - t427 * t445 + t485;
t342 = -t388 * t461 + t427 * t446 + t483;
t341 = t388 * t445 - t390 * t446 + t500;
t340 = t389 * t461 + (-t426 - t442) * t445 + t484;
t339 = t426 * t446 + (-t387 - t401) * t461 + t481;
t338 = t387 * t445 + (-t389 - t402) * t446 + t494;
t337 = t359 * t448 - t373 * t412 + t480;
t336 = -t358 * t448 + t373 * t413 + t479;
t335 = t358 * t412 - t359 * t413 + t482;
t334 = t351 * t437 + t361 * t448 - t365 * t405 - t395 * t412 + t480;
t333 = -t350 * t437 - t360 * t448 + t365 * t406 + t395 * t413 + t479;
t332 = t350 * t405 - t351 * t406 + t360 * t412 - t361 * t413 + t482;
t1 = ((t472 * t450 + t476 * t486) * qJD(1) + (t472 ^ 2 * t417 + (t488 * t476 + (-t416 + t487) * t472) * t476) * qJD(2)) * t465 / 0.2e1 - ((-t476 * t450 + t472 * t486) * qJD(1) + (t476 ^ 2 * t416 + (t487 * t472 + (-t417 + t488) * t476) * t472) * qJD(2)) * t499 / 0.2e1 + m(3) * (t366 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + t405 * ((-t345 * t503 + t393 * t347 + t394 * t349) * t405 + (-t344 * t503 + t346 * t393 + t348 * t394) * t406 + (-t362 * t503 + t363 * t393 + t364 * t394) * t437) / 0.2e1 + t406 * ((-t345 * t504 + t347 * t391 + t349 * t392) * t405 + (-t344 * t504 + t391 * t346 + t392 * t348) * t406 + (-t362 * t504 + t363 * t391 + t364 * t392) * t437) / 0.2e1 + t437 * ((t345 * t475 + t347 * t410 + t349 * t411) * t405 + (t344 * t475 + t346 * t410 + t348 * t411) * t406 + (t475 * t362 + t410 * t363 + t411 * t364) * t437) / 0.2e1 + qJD(1) * ((t475 * t451 + t471 * t452) * qJD(1) + ((t421 * t475 + t425 * t471) * t472 - (t420 * t475 + t424 * t471) * t476) * qJD(2)) / 0.2e1 + t448 * ((t353 * t475 + t355 * t432 + t357 * t433) * t412 + (t352 * t475 + t354 * t432 + t356 * t433) * t413 + (t475 * t370 + t432 * t371 + t433 * t372) * t448) / 0.2e1 + t412 * ((-t353 * t503 + t398 * t355 + t399 * t357) * t412 + (-t352 * t503 + t354 * t398 + t356 * t399) * t413 + (-t370 * t503 + t371 * t398 + t372 * t399) * t448) / 0.2e1 + t413 * ((-t353 * t504 + t355 * t396 + t357 * t397) * t412 + (-t352 * t504 + t396 * t354 + t397 * t356) * t413 + (-t370 * t504 + t371 * t396 + t372 * t397) * t448) / 0.2e1 + m(7) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(5) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(6) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(4) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + ((t518 * t440 + t516 * t441 + t517 * t503) * t461 + (t524 * t440 + t520 * t441 + t522 * t503) * t446 + (t523 * t440 + t519 * t441 + t521 * t503) * t445) * t445 / 0.2e1 + ((t518 * t438 + t516 * t439 + t517 * t504) * t461 + (t524 * t438 + t520 * t439 + t522 * t504) * t446 + (t523 * t438 + t519 * t439 + t521 * t504) * t445) * t446 / 0.2e1 + ((-t521 * t445 - t522 * t446 - t517 * t461) * t475 + ((t518 * t470 + t516 * t474) * t461 + (t524 * t470 + t520 * t474) * t446 + (t523 * t470 + t519 * t474) * t445) * t471) * t461 / 0.2e1 + (Icges(2,3) + m(2) * (t454 ^ 2 + t455 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
