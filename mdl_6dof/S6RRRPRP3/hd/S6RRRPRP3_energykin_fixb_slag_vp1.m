% Calculate kinetic energy for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:05
% EndTime: 2019-03-09 16:39:08
% DurationCPUTime: 2.68s
% Computational Cost: add. (1928->286), mult. (2160->439), div. (0->0), fcn. (2127->10), ass. (0->155)
t551 = Icges(6,1) + Icges(7,1);
t550 = -Icges(6,4) + Icges(7,5);
t549 = Icges(7,4) + Icges(6,5);
t548 = Icges(6,2) + Icges(7,3);
t547 = -Icges(7,6) + Icges(6,6);
t546 = -Icges(6,3) - Icges(7,2);
t545 = rSges(7,1) + pkin(5);
t544 = rSges(7,3) + qJ(6);
t467 = pkin(10) + qJ(5);
t462 = cos(t467);
t468 = qJ(2) + qJ(3);
t466 = cos(t468);
t475 = cos(qJ(1));
t461 = sin(t467);
t473 = sin(qJ(1));
t520 = t461 * t473;
t419 = t462 * t475 + t466 * t520;
t512 = t473 * t462;
t420 = -t461 * t475 + t466 * t512;
t465 = sin(t468);
t519 = t465 * t473;
t543 = t548 * t419 + t550 * t420 - t547 * t519;
t517 = t466 * t475;
t421 = t461 * t517 - t512;
t422 = t462 * t517 + t520;
t518 = t465 * t475;
t542 = t548 * t421 + t550 * t422 - t547 * t518;
t541 = -t547 * t419 + t549 * t420 - t546 * t519;
t540 = -t547 * t421 + t549 * t422 - t546 * t518;
t539 = t550 * t419 + t551 * t420 + t549 * t519;
t538 = t550 * t421 + t551 * t422 + t549 * t518;
t537 = t547 * t466 + (t548 * t461 + t550 * t462) * t465;
t536 = t546 * t466 + (-t547 * t461 + t549 * t462) * t465;
t535 = -t549 * t466 + (t550 * t461 + t551 * t462) * t465;
t474 = cos(qJ(2));
t527 = pkin(2) * t474;
t470 = cos(pkin(10));
t526 = pkin(4) * t470;
t472 = sin(qJ(2));
t524 = Icges(3,4) * t472;
t523 = Icges(3,4) * t474;
t522 = Icges(4,4) * t465;
t521 = Icges(4,4) * t466;
t469 = sin(pkin(10));
t516 = t469 * t473;
t515 = t469 * t475;
t514 = t470 * t473;
t513 = t470 * t475;
t510 = rSges(7,2) * t519 + t544 * t419 + t545 * t420;
t509 = rSges(7,2) * t518 + t544 * t421 + t545 * t422;
t508 = -rSges(7,2) * t466 + (t544 * t461 + t545 * t462) * t465;
t408 = -pkin(8) * t475 + t473 * t527;
t409 = pkin(8) * t473 + t475 * t527;
t464 = qJD(2) * t473;
t505 = qJD(2) * t475;
t507 = t408 * t464 + t409 * t505;
t456 = pkin(1) * t473 - pkin(7) * t475;
t506 = -t408 - t456;
t446 = qJD(3) * t473 + t464;
t504 = qJD(4) * t465;
t503 = qJD(5) * t465;
t502 = pkin(2) * qJD(2) * t472;
t497 = pkin(3) * t466 + qJ(4) * t465;
t434 = t497 * t473;
t501 = -t434 + t506;
t447 = (-qJD(2) - qJD(3)) * t475;
t500 = t475 * t502;
t499 = rSges(3,1) * t474 - rSges(3,2) * t472;
t498 = rSges(4,1) * t466 - rSges(4,2) * t465;
t496 = Icges(3,1) * t474 - t524;
t495 = Icges(4,1) * t466 - t522;
t494 = -Icges(3,2) * t472 + t523;
t493 = -Icges(4,2) * t465 + t521;
t492 = Icges(3,5) * t474 - Icges(3,6) * t472;
t491 = Icges(4,5) * t466 - Icges(4,6) * t465;
t426 = -Icges(3,6) * t475 + t473 * t494;
t428 = -Icges(3,5) * t475 + t473 * t496;
t490 = t426 * t472 - t428 * t474;
t427 = Icges(3,6) * t473 + t475 * t494;
t429 = Icges(3,5) * t473 + t475 * t496;
t489 = -t427 * t472 + t429 * t474;
t449 = Icges(3,2) * t474 + t524;
t450 = Icges(3,1) * t472 + t523;
t488 = -t449 * t472 + t450 * t474;
t487 = -qJD(4) * t466 + t446 * t434 + t507;
t445 = qJD(1) * (pkin(1) * t475 + pkin(7) * t473);
t486 = qJD(1) * t409 - t473 * t502 + t445;
t443 = pkin(3) * t465 - qJ(4) * t466;
t485 = t447 * t443 + t475 * t504 - t500;
t484 = (Icges(4,5) * t465 + Icges(4,6) * t466) * qJD(1) + (-Icges(4,3) * t475 + t473 * t491) * t447 + (Icges(4,3) * t473 + t475 * t491) * t446;
t483 = pkin(9) * t465 + t466 * t526;
t435 = t497 * t475;
t482 = qJD(1) * t435 + t473 * t504 + t486;
t374 = -pkin(4) * t515 + t473 * t483;
t375 = pkin(4) * t516 + t475 * t483;
t481 = t446 * t374 + (-t375 - t435) * t447 + t487;
t390 = -pkin(9) * t466 + t465 * t526;
t480 = qJD(1) * t375 + (-t390 - t443) * t446 + t482;
t479 = t447 * t390 + (-t374 + t501) * qJD(1) + t485;
t412 = -Icges(4,6) * t475 + t473 * t493;
t413 = Icges(4,6) * t473 + t475 * t493;
t414 = -Icges(4,5) * t475 + t473 * t495;
t415 = Icges(4,5) * t473 + t475 * t495;
t441 = Icges(4,2) * t466 + t522;
t442 = Icges(4,1) * t465 + t521;
t478 = (-t413 * t465 + t415 * t466) * t446 + (-t412 * t465 + t414 * t466) * t447 + (-t441 * t465 + t442 * t466) * qJD(1);
t457 = -qJD(5) * t466 + qJD(1);
t453 = rSges(2,1) * t475 - rSges(2,2) * t473;
t452 = rSges(2,1) * t473 + rSges(2,2) * t475;
t451 = rSges(3,1) * t472 + rSges(3,2) * t474;
t448 = Icges(3,5) * t472 + Icges(3,6) * t474;
t444 = rSges(4,1) * t465 + rSges(4,2) * t466;
t439 = t466 * t513 + t516;
t438 = -t466 * t515 + t514;
t437 = t466 * t514 - t515;
t436 = -t466 * t516 - t513;
t433 = rSges(3,3) * t473 + t475 * t499;
t432 = -rSges(3,3) * t475 + t473 * t499;
t431 = t473 * t503 + t447;
t430 = t475 * t503 + t446;
t425 = Icges(3,3) * t473 + t475 * t492;
t424 = -Icges(3,3) * t475 + t473 * t492;
t418 = rSges(4,3) * t473 + t475 * t498;
t417 = -rSges(4,3) * t475 + t473 * t498;
t406 = -rSges(5,3) * t466 + (rSges(5,1) * t470 - rSges(5,2) * t469) * t465;
t404 = -Icges(5,5) * t466 + (Icges(5,1) * t470 - Icges(5,4) * t469) * t465;
t403 = -Icges(5,6) * t466 + (Icges(5,4) * t470 - Icges(5,2) * t469) * t465;
t402 = -Icges(5,3) * t466 + (Icges(5,5) * t470 - Icges(5,6) * t469) * t465;
t398 = -rSges(6,3) * t466 + (rSges(6,1) * t462 - rSges(6,2) * t461) * t465;
t388 = qJD(1) * t433 - t451 * t464 + t445;
t387 = -t451 * t505 + (-t432 - t456) * qJD(1);
t384 = (t432 * t473 + t433 * t475) * qJD(2);
t383 = rSges(5,1) * t439 + rSges(5,2) * t438 + rSges(5,3) * t518;
t382 = rSges(5,1) * t437 + rSges(5,2) * t436 + rSges(5,3) * t519;
t381 = Icges(5,1) * t439 + Icges(5,4) * t438 + Icges(5,5) * t518;
t380 = Icges(5,1) * t437 + Icges(5,4) * t436 + Icges(5,5) * t519;
t379 = Icges(5,4) * t439 + Icges(5,2) * t438 + Icges(5,6) * t518;
t378 = Icges(5,4) * t437 + Icges(5,2) * t436 + Icges(5,6) * t519;
t377 = Icges(5,5) * t439 + Icges(5,6) * t438 + Icges(5,3) * t518;
t376 = Icges(5,5) * t437 + Icges(5,6) * t436 + Icges(5,3) * t519;
t372 = rSges(6,1) * t422 - rSges(6,2) * t421 + rSges(6,3) * t518;
t370 = rSges(6,1) * t420 - rSges(6,2) * t419 + rSges(6,3) * t519;
t355 = qJD(1) * t418 - t444 * t446 + t486;
t354 = -t500 + t444 * t447 + (-t417 + t506) * qJD(1);
t353 = t417 * t446 - t418 * t447 + t507;
t352 = qJD(1) * t383 + (-t406 - t443) * t446 + t482;
t351 = t406 * t447 + (-t382 + t501) * qJD(1) + t485;
t350 = t382 * t446 + (-t383 - t435) * t447 + t487;
t349 = t372 * t457 - t398 * t430 + t480;
t348 = -t370 * t457 + t398 * t431 + t479;
t347 = t370 * t430 - t372 * t431 + t481;
t346 = qJD(6) * t419 - t430 * t508 + t457 * t509 + t480;
t345 = qJD(6) * t421 + t431 * t508 - t457 * t510 + t479;
t344 = qJD(6) * t461 * t465 + t430 * t510 - t431 * t509 + t481;
t1 = m(3) * (t384 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + m(4) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(5) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(6) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(7) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 - ((-t475 * t448 + t473 * t488) * qJD(1) + (t475 ^ 2 * t424 + (t489 * t473 + (-t425 + t490) * t475) * t473) * qJD(2)) * t505 / 0.2e1 + ((t473 * t448 + t475 * t488) * qJD(1) + (t473 ^ 2 * t425 + (t490 * t475 + (-t424 + t489) * t473) * t475) * qJD(2)) * t464 / 0.2e1 + ((t537 * t421 + t535 * t422 + t536 * t518) * t457 + (t543 * t421 + t539 * t422 + t541 * t518) * t431 + (t542 * t421 + t538 * t422 + t540 * t518) * t430) * t430 / 0.2e1 + ((t537 * t419 + t535 * t420 + t536 * t519) * t457 + (t543 * t419 + t539 * t420 + t541 * t519) * t431 + (t542 * t419 + t538 * t420 + t540 * t519) * t430) * t431 / 0.2e1 + (t473 * t484 + t475 * t478 + (t377 * t518 + t438 * t379 + t439 * t381) * t446 + (t376 * t518 + t378 * t438 + t380 * t439) * t447 + (t402 * t518 + t403 * t438 + t404 * t439) * qJD(1)) * t446 / 0.2e1 + (t473 * t478 - t475 * t484 + (t377 * t519 + t379 * t436 + t381 * t437) * t446 + (t376 * t519 + t436 * t378 + t437 * t380) * t447 + (t402 * t519 + t403 * t436 + t404 * t437) * qJD(1)) * t447 / 0.2e1 + ((-t540 * t430 - t541 * t431 - t536 * t457) * t466 + ((t537 * t461 + t535 * t462) * t457 + (t543 * t461 + t539 * t462) * t431 + (t542 * t461 + t538 * t462) * t430) * t465) * t457 / 0.2e1 + (m(2) * (t452 ^ 2 + t453 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t427 * t474 + t429 * t472) * t473 - (t426 * t474 + t472 * t428) * t475) * qJD(2) + (t413 * t466 + t415 * t465) * t446 + (t412 * t466 + t414 * t465) * t447 + (-t376 * t447 - t377 * t446) * t466 + ((-t379 * t469 + t381 * t470) * t446 + (-t378 * t469 + t380 * t470) * t447) * t465 + (t474 * t449 + t472 * t450 + (t441 - t402) * t466 + (-t403 * t469 + t404 * t470 + t442) * t465) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
