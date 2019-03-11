% Calculate kinetic energy for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:10
% EndTime: 2019-03-09 22:06:13
% DurationCPUTime: 3.19s
% Computational Cost: add. (2199->331), mult. (2263->516), div. (0->0), fcn. (2220->12), ass. (0->169)
t541 = -Icges(6,3) - Icges(5,3);
t471 = qJ(4) + pkin(11);
t462 = sin(t471);
t463 = cos(t471);
t479 = cos(qJ(1));
t472 = qJ(2) + qJ(3);
t468 = cos(t472);
t476 = sin(qJ(1));
t521 = t468 * t476;
t414 = -t462 * t521 - t463 * t479;
t415 = -t462 * t479 + t463 * t521;
t477 = cos(qJ(4));
t516 = t477 * t479;
t474 = sin(qJ(4));
t519 = t474 * t476;
t431 = -t468 * t519 - t516;
t517 = t476 * t477;
t518 = t474 * t479;
t432 = t468 * t517 - t518;
t467 = sin(t472);
t523 = t467 * t476;
t540 = Icges(5,5) * t432 + Icges(6,5) * t415 + Icges(5,6) * t431 + Icges(6,6) * t414 - t541 * t523;
t520 = t468 * t479;
t416 = -t462 * t520 + t463 * t476;
t417 = t462 * t476 + t463 * t520;
t433 = -t468 * t518 + t517;
t434 = t468 * t516 + t519;
t522 = t467 * t479;
t539 = Icges(5,5) * t434 + Icges(6,5) * t417 + Icges(5,6) * t433 + Icges(6,6) * t416 - t541 * t522;
t538 = t541 * t468 + (Icges(5,5) * t477 + Icges(6,5) * t463 - Icges(5,6) * t474 - Icges(6,6) * t462) * t467;
t478 = cos(qJ(2));
t531 = pkin(2) * t478;
t530 = t477 * pkin(4);
t475 = sin(qJ(2));
t527 = Icges(3,4) * t475;
t526 = Icges(3,4) * t478;
t525 = Icges(4,4) * t467;
t524 = Icges(4,4) * t468;
t403 = -pkin(8) * t479 + t476 * t531;
t404 = pkin(8) * t476 + t479 * t531;
t466 = qJD(2) * t476;
t511 = qJD(2) * t479;
t515 = t403 * t466 + t404 * t511;
t456 = pkin(1) * t476 - pkin(7) * t479;
t514 = -t403 - t456;
t513 = pkin(5) * t463;
t444 = qJD(3) * t476 + t466;
t510 = qJD(4) * t467;
t509 = qJD(5) * t467;
t508 = qJD(6) * t467;
t507 = pkin(2) * qJD(2) * t475;
t425 = t479 * t510 + t444;
t506 = pkin(5) * t462;
t445 = (-qJD(2) - qJD(3)) * t479;
t505 = t479 * t507;
t504 = pkin(3) * t468 + pkin(9) * t467;
t503 = rSges(3,1) * t478 - rSges(3,2) * t475;
t502 = rSges(4,1) * t468 - rSges(4,2) * t467;
t501 = Icges(3,1) * t478 - t527;
t500 = Icges(4,1) * t468 - t525;
t499 = -Icges(3,2) * t475 + t526;
t498 = -Icges(4,2) * t467 + t524;
t497 = Icges(3,5) * t478 - Icges(3,6) * t475;
t496 = Icges(4,5) * t468 - Icges(4,6) * t467;
t421 = -Icges(3,6) * t479 + t476 * t499;
t423 = -Icges(3,5) * t479 + t476 * t501;
t495 = t421 * t475 - t423 * t478;
t422 = Icges(3,6) * t476 + t479 * t499;
t424 = Icges(3,5) * t476 + t479 * t501;
t494 = -t422 * t475 + t424 * t478;
t447 = Icges(3,2) * t478 + t527;
t448 = Icges(3,1) * t475 + t526;
t493 = -t447 * t475 + t448 * t478;
t426 = t476 * t510 + t445;
t429 = t504 * t476;
t430 = t504 * t479;
t492 = t444 * t429 - t430 * t445 + t515;
t443 = qJD(1) * (pkin(1) * t479 + pkin(7) * t476);
t491 = qJD(1) * t404 - t476 * t507 + t443;
t490 = qJ(5) * t467 + t468 * t530;
t489 = (Icges(4,5) * t467 + Icges(4,6) * t468) * qJD(1) + (-Icges(4,3) * t479 + t476 * t496) * t445 + (Icges(4,3) * t476 + t479 * t496) * t444;
t488 = pkin(10) * t467 + t468 * t513;
t366 = -pkin(4) * t518 + t476 * t490;
t487 = -qJD(5) * t468 + t425 * t366 + t492;
t440 = pkin(3) * t467 - pkin(9) * t468;
t486 = qJD(1) * t430 - t440 * t444 + t491;
t485 = t445 * t440 + (-t429 + t514) * qJD(1) - t505;
t367 = pkin(4) * t519 + t479 * t490;
t457 = -qJD(4) * t468 + qJD(1);
t484 = t457 * t367 + t476 * t509 + t486;
t383 = -qJ(5) * t468 + t467 * t530;
t483 = t426 * t383 + t479 * t509 + t485;
t408 = -Icges(4,6) * t479 + t476 * t498;
t409 = Icges(4,6) * t476 + t479 * t498;
t410 = -Icges(4,5) * t479 + t476 * t500;
t411 = Icges(4,5) * t476 + t479 * t500;
t437 = Icges(4,2) * t468 + t525;
t438 = Icges(4,1) * t467 + t524;
t482 = (-t409 * t467 + t411 * t468) * t444 + (-t408 * t467 + t410 * t468) * t445 + (-t437 * t467 + t438 * t468) * qJD(1);
t464 = qJ(6) + t471;
t459 = cos(t464);
t458 = sin(t464);
t451 = rSges(2,1) * t479 - rSges(2,2) * t476;
t450 = rSges(2,1) * t476 + rSges(2,2) * t479;
t449 = rSges(3,1) * t475 + rSges(3,2) * t478;
t446 = Icges(3,5) * t475 + Icges(3,6) * t478;
t439 = rSges(4,1) * t467 + rSges(4,2) * t468;
t435 = qJD(1) + (-qJD(4) - qJD(6)) * t468;
t428 = rSges(3,3) * t476 + t479 * t503;
t427 = -rSges(3,3) * t479 + t476 * t503;
t420 = Icges(3,3) * t476 + t479 * t497;
t419 = -Icges(3,3) * t479 + t476 * t497;
t413 = rSges(4,3) * t476 + t479 * t502;
t412 = -rSges(4,3) * t479 + t476 * t502;
t402 = t458 * t476 + t459 * t520;
t401 = -t458 * t520 + t459 * t476;
t400 = -t458 * t479 + t459 * t521;
t399 = -t458 * t521 - t459 * t479;
t398 = -rSges(5,3) * t468 + (rSges(5,1) * t477 - rSges(5,2) * t474) * t467;
t397 = -Icges(5,5) * t468 + (Icges(5,1) * t477 - Icges(5,4) * t474) * t467;
t396 = -Icges(5,6) * t468 + (Icges(5,4) * t477 - Icges(5,2) * t474) * t467;
t391 = t476 * t508 + t426;
t390 = t479 * t508 + t425;
t388 = -rSges(6,3) * t468 + (rSges(6,1) * t463 - rSges(6,2) * t462) * t467;
t387 = -Icges(6,5) * t468 + (Icges(6,1) * t463 - Icges(6,4) * t462) * t467;
t386 = -Icges(6,6) * t468 + (Icges(6,4) * t463 - Icges(6,2) * t462) * t467;
t384 = -rSges(7,3) * t468 + (rSges(7,1) * t459 - rSges(7,2) * t458) * t467;
t382 = -Icges(7,5) * t468 + (Icges(7,1) * t459 - Icges(7,4) * t458) * t467;
t381 = -Icges(7,6) * t468 + (Icges(7,4) * t459 - Icges(7,2) * t458) * t467;
t380 = -Icges(7,3) * t468 + (Icges(7,5) * t459 - Icges(7,6) * t458) * t467;
t379 = -pkin(10) * t468 + t467 * t513;
t378 = qJD(1) * t428 - t449 * t466 + t443;
t377 = -t449 * t511 + (-t427 - t456) * qJD(1);
t376 = rSges(5,1) * t434 + rSges(5,2) * t433 + rSges(5,3) * t522;
t375 = rSges(5,1) * t432 + rSges(5,2) * t431 + rSges(5,3) * t523;
t374 = Icges(5,1) * t434 + Icges(5,4) * t433 + Icges(5,5) * t522;
t373 = Icges(5,1) * t432 + Icges(5,4) * t431 + Icges(5,5) * t523;
t372 = Icges(5,4) * t434 + Icges(5,2) * t433 + Icges(5,6) * t522;
t371 = Icges(5,4) * t432 + Icges(5,2) * t431 + Icges(5,6) * t523;
t368 = (t427 * t476 + t428 * t479) * qJD(2);
t364 = rSges(6,1) * t417 + rSges(6,2) * t416 + rSges(6,3) * t522;
t363 = rSges(6,1) * t415 + rSges(6,2) * t414 + rSges(6,3) * t523;
t362 = Icges(6,1) * t417 + Icges(6,4) * t416 + Icges(6,5) * t522;
t361 = Icges(6,1) * t415 + Icges(6,4) * t414 + Icges(6,5) * t523;
t360 = Icges(6,4) * t417 + Icges(6,2) * t416 + Icges(6,6) * t522;
t359 = Icges(6,4) * t415 + Icges(6,2) * t414 + Icges(6,6) * t523;
t355 = rSges(7,1) * t402 + rSges(7,2) * t401 + rSges(7,3) * t522;
t354 = rSges(7,1) * t400 + rSges(7,2) * t399 + rSges(7,3) * t523;
t353 = Icges(7,1) * t402 + Icges(7,4) * t401 + Icges(7,5) * t522;
t352 = Icges(7,1) * t400 + Icges(7,4) * t399 + Icges(7,5) * t523;
t351 = Icges(7,4) * t402 + Icges(7,2) * t401 + Icges(7,6) * t522;
t350 = Icges(7,4) * t400 + Icges(7,2) * t399 + Icges(7,6) * t523;
t349 = Icges(7,5) * t402 + Icges(7,6) * t401 + Icges(7,3) * t522;
t348 = Icges(7,5) * t400 + Icges(7,6) * t399 + Icges(7,3) * t523;
t346 = t476 * t506 + t479 * t488;
t345 = t476 * t488 - t479 * t506;
t344 = qJD(1) * t413 - t439 * t444 + t491;
t343 = -t505 + t439 * t445 + (-t412 + t514) * qJD(1);
t342 = t412 * t444 - t413 * t445 + t515;
t341 = t376 * t457 - t398 * t425 + t486;
t340 = -t375 * t457 + t398 * t426 + t485;
t339 = t375 * t425 - t376 * t426 + t492;
t338 = t364 * t457 + (-t383 - t388) * t425 + t484;
t337 = t388 * t426 + (-t363 - t366) * t457 + t483;
t336 = t363 * t425 + (-t364 - t367) * t426 + t487;
t335 = t346 * t457 + t355 * t435 - t384 * t390 + (-t379 - t383) * t425 + t484;
t334 = -t354 * t435 + t379 * t426 + t384 * t391 + (-t345 - t366) * t457 + t483;
t333 = t345 * t425 + t354 * t390 - t355 * t391 + (-t346 - t367) * t426 + t487;
t1 = t390 * ((t349 * t522 + t401 * t351 + t402 * t353) * t390 + (t348 * t522 + t350 * t401 + t352 * t402) * t391 + (t380 * t522 + t381 * t401 + t382 * t402) * t435) / 0.2e1 + t391 * ((t349 * t523 + t351 * t399 + t353 * t400) * t390 + (t348 * t523 + t399 * t350 + t400 * t352) * t391 + (t380 * t523 + t381 * t399 + t382 * t400) * t435) / 0.2e1 + t435 * ((-t348 * t391 - t349 * t390 - t380 * t435) * t468 + ((-t351 * t458 + t353 * t459) * t390 + (-t350 * t458 + t352 * t459) * t391 + (-t381 * t458 + t382 * t459) * t435) * t467) / 0.2e1 + t444 * (t476 * t489 + t479 * t482) / 0.2e1 + t445 * (t476 * t482 - t479 * t489) / 0.2e1 + m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(5) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(6) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(4) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(3) * (t368 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 - ((-t479 * t446 + t476 * t493) * qJD(1) + (t479 ^ 2 * t419 + (t494 * t476 + (-t420 + t495) * t479) * t476) * qJD(2)) * t511 / 0.2e1 + ((t476 * t446 + t479 * t493) * qJD(1) + (t476 ^ 2 * t420 + (t495 * t479 + (-t419 + t494) * t476) * t479) * qJD(2)) * t466 / 0.2e1 + ((t386 * t416 + t387 * t417 + t396 * t433 + t397 * t434 + t538 * t522) * t457 + (t359 * t416 + t361 * t417 + t371 * t433 + t373 * t434 + t540 * t522) * t426 + (t360 * t416 + t362 * t417 + t372 * t433 + t374 * t434 + t539 * t522) * t425) * t425 / 0.2e1 + ((t386 * t414 + t387 * t415 + t396 * t431 + t397 * t432 + t538 * t523) * t457 + (t359 * t414 + t361 * t415 + t371 * t431 + t373 * t432 + t540 * t523) * t426 + (t360 * t414 + t362 * t415 + t372 * t431 + t374 * t432 + t539 * t523) * t425) * t426 / 0.2e1 + ((-t539 * t425 - t540 * t426 - t538 * t457) * t468 + ((-t386 * t462 + t387 * t463 - t396 * t474 + t397 * t477) * t457 + (-t359 * t462 + t361 * t463 - t371 * t474 + t373 * t477) * t426 + (-t360 * t462 + t362 * t463 - t372 * t474 + t374 * t477) * t425) * t467) * t457 / 0.2e1 + (m(2) * (t450 ^ 2 + t451 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t409 * t468 + t411 * t467) * t444 + (t408 * t468 + t410 * t467) * t445 + ((t422 * t478 + t424 * t475) * t476 - (t421 * t478 + t423 * t475) * t479) * qJD(2) + (t468 * t437 + t467 * t438 + t478 * t447 + t475 * t448) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
