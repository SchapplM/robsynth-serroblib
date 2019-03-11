% Calculate kinetic energy for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:10
% EndTime: 2019-03-10 03:39:13
% DurationCPUTime: 2.65s
% Computational Cost: add. (2271->337), mult. (2323->549), div. (0->0), fcn. (2280->12), ass. (0->175)
t477 = cos(qJ(2));
t532 = pkin(2) * t477;
t476 = cos(qJ(4));
t531 = t476 * pkin(4);
t474 = sin(qJ(2));
t528 = Icges(3,4) * t474;
t527 = Icges(3,4) * t477;
t472 = qJ(2) + qJ(3);
t465 = sin(t472);
t526 = Icges(4,4) * t465;
t467 = cos(t472);
t525 = Icges(4,4) * t467;
t475 = sin(qJ(1));
t524 = t465 * t475;
t478 = cos(qJ(1));
t523 = t465 * t478;
t522 = t467 * t475;
t521 = t467 * t478;
t473 = sin(qJ(4));
t520 = t473 * t475;
t519 = t473 * t478;
t518 = t475 * t476;
t517 = t476 * t478;
t471 = qJ(4) + qJ(5);
t398 = -pkin(8) * t478 + t475 * t532;
t399 = pkin(8) * t475 + t478 * t532;
t463 = qJD(2) * t475;
t512 = qJD(2) * t478;
t516 = t398 * t463 + t399 * t512;
t456 = pkin(1) * t475 - pkin(7) * t478;
t515 = -t398 - t456;
t466 = cos(t471);
t514 = pkin(5) * t466;
t444 = qJD(3) * t475 + t463;
t511 = qJD(4) * t465;
t510 = qJD(5) * t465;
t509 = qJD(6) * t465;
t508 = -qJD(4) - qJD(5);
t507 = pkin(2) * qJD(2) * t474;
t425 = t478 * t511 + t444;
t464 = sin(t471);
t506 = pkin(5) * t464;
t445 = (-qJD(2) - qJD(3)) * t478;
t505 = t478 * t507;
t388 = t478 * t510 + t425;
t504 = pkin(3) * t467 + pkin(9) * t465;
t503 = rSges(3,1) * t477 - rSges(3,2) * t474;
t502 = rSges(4,1) * t467 - rSges(4,2) * t465;
t501 = Icges(3,1) * t477 - t528;
t500 = Icges(4,1) * t467 - t526;
t499 = -Icges(3,2) * t474 + t527;
t498 = -Icges(4,2) * t465 + t525;
t497 = Icges(3,5) * t477 - Icges(3,6) * t474;
t496 = Icges(4,5) * t467 - Icges(4,6) * t465;
t421 = -Icges(3,6) * t478 + t475 * t499;
t423 = -Icges(3,5) * t478 + t475 * t501;
t495 = t421 * t474 - t423 * t477;
t422 = Icges(3,6) * t475 + t478 * t499;
t424 = Icges(3,5) * t475 + t478 * t501;
t494 = -t422 * t474 + t424 * t477;
t447 = Icges(3,2) * t477 + t528;
t448 = Icges(3,1) * t474 + t527;
t493 = -t447 * t474 + t448 * t477;
t426 = t475 * t511 + t445;
t429 = t504 * t475;
t430 = t504 * t478;
t492 = t444 * t429 - t430 * t445 + t516;
t389 = t475 * t510 + t426;
t442 = qJD(1) * (pkin(1) * t478 + pkin(7) * t475);
t491 = qJD(1) * t399 - t475 * t507 + t442;
t490 = pkin(10) * t465 + t467 * t531;
t489 = (Icges(4,5) * t465 + Icges(4,6) * t467) * qJD(1) + (-Icges(4,3) * t478 + t475 * t496) * t445 + (Icges(4,3) * t475 + t478 * t496) * t444;
t488 = pkin(11) * t465 + t467 * t514;
t363 = -pkin(4) * t519 + t475 * t490;
t364 = pkin(4) * t520 + t478 * t490;
t487 = t425 * t363 - t364 * t426 + t492;
t440 = pkin(3) * t465 - pkin(9) * t467;
t486 = qJD(1) * t430 - t440 * t444 + t491;
t485 = t445 * t440 + (-t429 + t515) * qJD(1) - t505;
t379 = -pkin(10) * t467 + t465 * t531;
t457 = -qJD(4) * t467 + qJD(1);
t484 = t457 * t364 - t379 * t425 + t486;
t483 = -t363 * t457 + t426 * t379 + t485;
t403 = -Icges(4,6) * t478 + t475 * t498;
t404 = Icges(4,6) * t475 + t478 * t498;
t405 = -Icges(4,5) * t478 + t475 * t500;
t406 = Icges(4,5) * t475 + t478 * t500;
t437 = Icges(4,2) * t467 + t526;
t438 = Icges(4,1) * t465 + t525;
t482 = (-t404 * t465 + t406 * t467) * t444 + (-t403 * t465 + t405 * t467) * t445 + (-t437 * t465 + t438 * t467) * qJD(1);
t468 = qJ(6) + t471;
t459 = cos(t468);
t458 = sin(t468);
t451 = rSges(2,1) * t478 - rSges(2,2) * t475;
t450 = rSges(2,1) * t475 + rSges(2,2) * t478;
t449 = rSges(3,1) * t474 + rSges(3,2) * t477;
t446 = Icges(3,5) * t474 + Icges(3,6) * t477;
t439 = rSges(4,1) * t465 + rSges(4,2) * t467;
t435 = t467 * t508 + qJD(1);
t434 = t467 * t517 + t520;
t433 = -t467 * t519 + t518;
t432 = t467 * t518 - t519;
t431 = -t467 * t520 - t517;
t428 = rSges(3,3) * t475 + t478 * t503;
t427 = -rSges(3,3) * t478 + t475 * t503;
t420 = Icges(3,3) * t475 + t478 * t497;
t419 = -Icges(3,3) * t478 + t475 * t497;
t417 = t464 * t475 + t466 * t521;
t416 = -t464 * t521 + t466 * t475;
t415 = -t464 * t478 + t466 * t522;
t414 = -t464 * t522 - t466 * t478;
t413 = qJD(1) + (-qJD(6) + t508) * t467;
t412 = rSges(4,3) * t475 + t478 * t502;
t411 = -rSges(4,3) * t478 + t475 * t502;
t410 = t458 * t475 + t459 * t521;
t409 = -t458 * t521 + t459 * t475;
t408 = -t458 * t478 + t459 * t522;
t407 = -t458 * t522 - t459 * t478;
t397 = -rSges(5,3) * t467 + (rSges(5,1) * t476 - rSges(5,2) * t473) * t465;
t396 = -Icges(5,5) * t467 + (Icges(5,1) * t476 - Icges(5,4) * t473) * t465;
t395 = -Icges(5,6) * t467 + (Icges(5,4) * t476 - Icges(5,2) * t473) * t465;
t394 = -Icges(5,3) * t467 + (Icges(5,5) * t476 - Icges(5,6) * t473) * t465;
t390 = -rSges(6,3) * t467 + (rSges(6,1) * t466 - rSges(6,2) * t464) * t465;
t386 = -Icges(6,5) * t467 + (Icges(6,1) * t466 - Icges(6,4) * t464) * t465;
t385 = -Icges(6,6) * t467 + (Icges(6,4) * t466 - Icges(6,2) * t464) * t465;
t384 = -Icges(6,3) * t467 + (Icges(6,5) * t466 - Icges(6,6) * t464) * t465;
t383 = -rSges(7,3) * t467 + (rSges(7,1) * t459 - rSges(7,2) * t458) * t465;
t382 = -Icges(7,5) * t467 + (Icges(7,1) * t459 - Icges(7,4) * t458) * t465;
t381 = -Icges(7,6) * t467 + (Icges(7,4) * t459 - Icges(7,2) * t458) * t465;
t380 = -Icges(7,3) * t467 + (Icges(7,5) * t459 - Icges(7,6) * t458) * t465;
t378 = t475 * t509 + t389;
t377 = t478 * t509 + t388;
t376 = -pkin(11) * t467 + t465 * t514;
t375 = qJD(1) * t428 - t449 * t463 + t442;
t374 = -t449 * t512 + (-t427 - t456) * qJD(1);
t373 = rSges(5,1) * t434 + rSges(5,2) * t433 + rSges(5,3) * t523;
t372 = rSges(5,1) * t432 + rSges(5,2) * t431 + rSges(5,3) * t524;
t371 = Icges(5,1) * t434 + Icges(5,4) * t433 + Icges(5,5) * t523;
t370 = Icges(5,1) * t432 + Icges(5,4) * t431 + Icges(5,5) * t524;
t369 = Icges(5,4) * t434 + Icges(5,2) * t433 + Icges(5,6) * t523;
t368 = Icges(5,4) * t432 + Icges(5,2) * t431 + Icges(5,6) * t524;
t367 = Icges(5,5) * t434 + Icges(5,6) * t433 + Icges(5,3) * t523;
t366 = Icges(5,5) * t432 + Icges(5,6) * t431 + Icges(5,3) * t524;
t365 = (t427 * t475 + t428 * t478) * qJD(2);
t361 = rSges(6,1) * t417 + rSges(6,2) * t416 + rSges(6,3) * t523;
t360 = rSges(6,1) * t415 + rSges(6,2) * t414 + rSges(6,3) * t524;
t359 = Icges(6,1) * t417 + Icges(6,4) * t416 + Icges(6,5) * t523;
t358 = Icges(6,1) * t415 + Icges(6,4) * t414 + Icges(6,5) * t524;
t357 = Icges(6,4) * t417 + Icges(6,2) * t416 + Icges(6,6) * t523;
t356 = Icges(6,4) * t415 + Icges(6,2) * t414 + Icges(6,6) * t524;
t355 = Icges(6,5) * t417 + Icges(6,6) * t416 + Icges(6,3) * t523;
t354 = Icges(6,5) * t415 + Icges(6,6) * t414 + Icges(6,3) * t524;
t352 = rSges(7,1) * t410 + rSges(7,2) * t409 + rSges(7,3) * t523;
t351 = rSges(7,1) * t408 + rSges(7,2) * t407 + rSges(7,3) * t524;
t350 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t523;
t349 = Icges(7,1) * t408 + Icges(7,4) * t407 + Icges(7,5) * t524;
t348 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t523;
t347 = Icges(7,4) * t408 + Icges(7,2) * t407 + Icges(7,6) * t524;
t346 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t523;
t345 = Icges(7,5) * t408 + Icges(7,6) * t407 + Icges(7,3) * t524;
t343 = t475 * t506 + t478 * t488;
t342 = t475 * t488 - t478 * t506;
t341 = qJD(1) * t412 - t439 * t444 + t491;
t340 = -t505 + t439 * t445 + (-t411 + t515) * qJD(1);
t339 = t411 * t444 - t412 * t445 + t516;
t338 = t373 * t457 - t397 * t425 + t486;
t337 = -t372 * t457 + t397 * t426 + t485;
t336 = t372 * t425 - t373 * t426 + t492;
t335 = t361 * t435 - t388 * t390 + t484;
t334 = -t360 * t435 + t389 * t390 + t483;
t333 = t360 * t388 - t361 * t389 + t487;
t332 = t343 * t435 + t352 * t413 - t376 * t388 - t377 * t383 + t484;
t331 = -t342 * t435 - t351 * t413 + t376 * t389 + t378 * t383 + t483;
t330 = t342 * t388 - t343 * t389 + t351 * t377 - t352 * t378 + t487;
t1 = -((-t478 * t446 + t475 * t493) * qJD(1) + (t478 ^ 2 * t419 + (t494 * t475 + (-t420 + t495) * t478) * t475) * qJD(2)) * t512 / 0.2e1 + ((t475 * t446 + t478 * t493) * qJD(1) + (t475 ^ 2 * t420 + (t495 * t478 + (-t419 + t494) * t475) * t478) * qJD(2)) * t463 / 0.2e1 + t413 * ((-t345 * t378 - t346 * t377 - t380 * t413) * t467 + ((-t348 * t458 + t350 * t459) * t377 + (-t347 * t458 + t349 * t459) * t378 + (-t381 * t458 + t382 * t459) * t413) * t465) / 0.2e1 + t388 * ((t355 * t523 + t416 * t357 + t417 * t359) * t388 + (t354 * t523 + t356 * t416 + t358 * t417) * t389 + (t384 * t523 + t385 * t416 + t386 * t417) * t435) / 0.2e1 + t389 * ((t355 * t524 + t357 * t414 + t359 * t415) * t388 + (t354 * t524 + t414 * t356 + t415 * t358) * t389 + (t384 * t524 + t385 * t414 + t386 * t415) * t435) / 0.2e1 + t435 * ((-t354 * t389 - t355 * t388 - t384 * t435) * t467 + ((-t357 * t464 + t359 * t466) * t388 + (-t356 * t464 + t358 * t466) * t389 + (-t385 * t464 + t386 * t466) * t435) * t465) / 0.2e1 + t377 * ((t346 * t523 + t409 * t348 + t410 * t350) * t377 + (t345 * t523 + t347 * t409 + t349 * t410) * t378 + (t380 * t523 + t381 * t409 + t382 * t410) * t413) / 0.2e1 + t378 * ((t346 * t524 + t348 * t407 + t350 * t408) * t377 + (t345 * t524 + t407 * t347 + t408 * t349) * t378 + (t380 * t524 + t381 * t407 + t382 * t408) * t413) / 0.2e1 + t457 * ((-t366 * t426 - t367 * t425 - t394 * t457) * t467 + ((-t369 * t473 + t371 * t476) * t425 + (-t368 * t473 + t370 * t476) * t426 + (-t395 * t473 + t396 * t476) * t457) * t465) / 0.2e1 + t425 * ((t367 * t523 + t433 * t369 + t434 * t371) * t425 + (t366 * t523 + t368 * t433 + t370 * t434) * t426 + (t394 * t523 + t395 * t433 + t396 * t434) * t457) / 0.2e1 + t426 * ((t367 * t524 + t369 * t431 + t371 * t432) * t425 + (t366 * t524 + t431 * t368 + t432 * t370) * t426 + (t394 * t524 + t395 * t431 + t396 * t432) * t457) / 0.2e1 + t444 * (t489 * t475 + t482 * t478) / 0.2e1 + t445 * (t482 * t475 - t489 * t478) / 0.2e1 + m(6) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(5) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(4) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(3) * (t365 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(7) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t450 ^ 2 + t451 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t404 * t467 + t406 * t465) * t444 + (t403 * t467 + t405 * t465) * t445 + ((t422 * t477 + t424 * t474) * t475 - (t421 * t477 + t423 * t474) * t478) * qJD(2) + (t467 * t437 + t465 * t438 + t477 * t447 + t474 * t448) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
