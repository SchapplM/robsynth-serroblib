% Calculate kinetic energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:44
% EndTime: 2019-03-09 15:47:46
% DurationCPUTime: 2.62s
% Computational Cost: add. (2778->336), mult. (5097->491), div. (0->0), fcn. (6048->12), ass. (0->156)
t570 = Icges(5,1) + Icges(6,2);
t569 = -Icges(5,4) - Icges(6,6);
t568 = Icges(6,4) - Icges(5,5);
t567 = Icges(6,5) - Icges(5,6);
t566 = Icges(5,2) + Icges(6,3);
t565 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t511 = sin(qJ(2));
t512 = sin(qJ(1));
t515 = cos(qJ(2));
t516 = cos(qJ(1));
t548 = cos(pkin(6));
t532 = t516 * t548;
t489 = t511 * t512 - t515 * t532;
t490 = t511 * t532 + t512 * t515;
t507 = sin(pkin(6));
t544 = t507 * t516;
t447 = Icges(3,5) * t490 - Icges(3,6) * t489 - Icges(3,3) * t544;
t533 = t512 * t548;
t491 = t516 * t511 + t515 * t533;
t492 = -t511 * t533 + t516 * t515;
t547 = t507 * t512;
t448 = Icges(3,5) * t492 - Icges(3,6) * t491 + Icges(3,3) * t547;
t564 = t507 * (t447 * t516 - t448 * t512);
t537 = qJ(3) + pkin(11);
t530 = cos(t537);
t528 = t507 * t530;
t529 = sin(t537);
t465 = t490 * t529 + t516 * t528;
t527 = t507 * t529;
t466 = t490 * t530 - t516 * t527;
t563 = t566 * t465 + t569 * t466 + t567 * t489;
t467 = t492 * t529 - t512 * t528;
t468 = t492 * t530 + t512 * t527;
t562 = t566 * t467 + t569 * t468 + t567 * t491;
t561 = t569 * t465 + t570 * t466 - t568 * t489;
t560 = t569 * t467 + t570 * t468 - t568 * t491;
t480 = t511 * t527 - t530 * t548;
t481 = t511 * t528 + t529 * t548;
t545 = t507 * t515;
t559 = t566 * t480 + t569 * t481 - t567 * t545;
t558 = t569 * t480 + t570 * t481 + t568 * t545;
t510 = sin(qJ(3));
t514 = cos(qJ(3));
t470 = -t490 * t510 - t514 * t544;
t535 = t510 * t544;
t471 = t490 * t514 - t535;
t557 = Icges(4,5) * t471 + Icges(4,6) * t470 + t567 * t465 - t568 * t466 + t565 * t489;
t546 = t507 * t514;
t472 = -t492 * t510 + t512 * t546;
t536 = t510 * t547;
t473 = t492 * t514 + t536;
t556 = Icges(4,5) * t473 + Icges(4,6) * t472 + t567 * t467 - t568 * t468 + t565 * t491;
t487 = -t507 * t511 * t510 + t514 * t548;
t531 = t548 * t510;
t488 = t511 * t546 + t531;
t555 = Icges(4,5) * t488 + Icges(4,6) * t487 + t567 * t480 - t568 * t481 - t565 * t545;
t550 = pkin(3) * t514;
t410 = -pkin(3) * t535 + qJ(4) * t489 + t490 * t550;
t422 = pkin(4) * t466 + qJ(5) * t465;
t543 = -t410 - t422;
t411 = pkin(3) * t536 + qJ(4) * t491 + t492 * t550;
t423 = pkin(4) * t468 + qJ(5) * t467;
t542 = -t411 - t423;
t441 = pkin(4) * t481 + qJ(5) * t480;
t445 = pkin(3) * t531 + (-qJ(4) * t515 + t511 * t550) * t507;
t541 = -t441 - t445;
t461 = pkin(2) * t490 + pkin(9) * t489;
t462 = pkin(2) * t492 + pkin(9) * t491;
t538 = qJD(2) * t507;
t503 = t512 * t538;
t534 = t516 * t538;
t540 = t461 * t503 + t462 * t534;
t474 = qJD(3) * t491 + t503;
t539 = qJD(1) * (pkin(1) * t512 - pkin(8) * t544);
t504 = qJD(2) * t548 + qJD(1);
t475 = qJD(3) * t489 - t534;
t494 = -qJD(3) * t545 + t504;
t493 = (pkin(2) * t511 - pkin(9) * t515) * t507;
t495 = qJD(1) * (pkin(1) * t516 + pkin(8) * t547);
t525 = t504 * t462 - t493 * t503 + t495;
t524 = -qJD(4) * t545 + t474 * t410 + t540;
t523 = qJD(4) * t489 + t494 * t411 + t525;
t522 = qJD(5) * t480 + t474 * t422 + t524;
t521 = -t461 * t504 - t493 * t534 - t539;
t520 = qJD(5) * t465 + t494 * t423 + t523;
t519 = qJD(4) * t491 + t475 * t445 + t521;
t518 = qJD(5) * t467 + t475 * t441 + t519;
t513 = cos(qJ(6));
t509 = sin(qJ(6));
t500 = rSges(2,1) * t516 - rSges(2,2) * t512;
t499 = rSges(2,1) * t512 + rSges(2,2) * t516;
t482 = t548 * rSges(3,3) + (rSges(3,1) * t511 + rSges(3,2) * t515) * t507;
t479 = Icges(3,5) * t548 + (Icges(3,1) * t511 + Icges(3,4) * t515) * t507;
t478 = Icges(3,6) * t548 + (Icges(3,4) * t511 + Icges(3,2) * t515) * t507;
t477 = Icges(3,3) * t548 + (Icges(3,5) * t511 + Icges(3,6) * t515) * t507;
t469 = -pkin(5) * t545 + pkin(10) * t481;
t464 = t480 * t509 - t513 * t545;
t463 = t480 * t513 + t509 * t545;
t458 = qJD(6) * t481 + t494;
t455 = rSges(3,1) * t492 - rSges(3,2) * t491 + rSges(3,3) * t547;
t454 = rSges(3,1) * t490 - rSges(3,2) * t489 - rSges(3,3) * t544;
t452 = Icges(3,1) * t492 - Icges(3,4) * t491 + Icges(3,5) * t547;
t451 = Icges(3,1) * t490 - Icges(3,4) * t489 - Icges(3,5) * t544;
t450 = Icges(3,4) * t492 - Icges(3,2) * t491 + Icges(3,6) * t547;
t449 = Icges(3,4) * t490 - Icges(3,2) * t489 - Icges(3,6) * t544;
t446 = rSges(4,1) * t488 + rSges(4,2) * t487 - rSges(4,3) * t545;
t444 = Icges(4,1) * t488 + Icges(4,4) * t487 - Icges(4,5) * t545;
t443 = Icges(4,4) * t488 + Icges(4,2) * t487 - Icges(4,6) * t545;
t440 = pkin(5) * t491 + pkin(10) * t468;
t439 = pkin(5) * t489 + pkin(10) * t466;
t438 = rSges(5,1) * t481 - rSges(5,2) * t480 - rSges(5,3) * t545;
t437 = -rSges(6,1) * t545 - rSges(6,2) * t481 + rSges(6,3) * t480;
t430 = t467 * t509 + t491 * t513;
t429 = t467 * t513 - t491 * t509;
t428 = t465 * t509 + t489 * t513;
t427 = t465 * t513 - t489 * t509;
t426 = qJD(6) * t466 + t475;
t425 = qJD(6) * t468 + t474;
t420 = rSges(4,1) * t473 + rSges(4,2) * t472 + rSges(4,3) * t491;
t419 = rSges(4,1) * t471 + rSges(4,2) * t470 + rSges(4,3) * t489;
t418 = Icges(4,1) * t473 + Icges(4,4) * t472 + Icges(4,5) * t491;
t417 = Icges(4,1) * t471 + Icges(4,4) * t470 + Icges(4,5) * t489;
t416 = Icges(4,4) * t473 + Icges(4,2) * t472 + Icges(4,6) * t491;
t415 = Icges(4,4) * t471 + Icges(4,2) * t470 + Icges(4,6) * t489;
t409 = rSges(5,1) * t468 - rSges(5,2) * t467 + rSges(5,3) * t491;
t408 = rSges(5,1) * t466 - rSges(5,2) * t465 + rSges(5,3) * t489;
t407 = rSges(6,1) * t491 - rSges(6,2) * t468 + rSges(6,3) * t467;
t406 = rSges(6,1) * t489 - rSges(6,2) * t466 + rSges(6,3) * t465;
t393 = t455 * t504 - t482 * t503 + t495;
t392 = -t454 * t504 - t482 * t534 - t539;
t391 = rSges(7,1) * t464 + rSges(7,2) * t463 + rSges(7,3) * t481;
t390 = Icges(7,1) * t464 + Icges(7,4) * t463 + Icges(7,5) * t481;
t389 = Icges(7,4) * t464 + Icges(7,2) * t463 + Icges(7,6) * t481;
t388 = Icges(7,5) * t464 + Icges(7,6) * t463 + Icges(7,3) * t481;
t385 = (t454 * t512 + t455 * t516) * t538;
t383 = rSges(7,1) * t430 + rSges(7,2) * t429 + rSges(7,3) * t468;
t382 = rSges(7,1) * t428 + rSges(7,2) * t427 + rSges(7,3) * t466;
t381 = Icges(7,1) * t430 + Icges(7,4) * t429 + Icges(7,5) * t468;
t380 = Icges(7,1) * t428 + Icges(7,4) * t427 + Icges(7,5) * t466;
t379 = Icges(7,4) * t430 + Icges(7,2) * t429 + Icges(7,6) * t468;
t378 = Icges(7,4) * t428 + Icges(7,2) * t427 + Icges(7,6) * t466;
t377 = Icges(7,5) * t430 + Icges(7,6) * t429 + Icges(7,3) * t468;
t376 = Icges(7,5) * t428 + Icges(7,6) * t427 + Icges(7,3) * t466;
t375 = t420 * t494 - t446 * t474 + t525;
t374 = -t419 * t494 + t446 * t475 + t521;
t373 = t419 * t474 - t420 * t475 + t540;
t372 = t409 * t494 + (-t438 - t445) * t474 + t523;
t371 = t438 * t475 + (-t408 - t410) * t494 + t519;
t370 = t408 * t474 + (-t409 - t411) * t475 + t524;
t369 = t407 * t494 + (-t437 + t541) * t474 + t520;
t368 = t437 * t475 + (-t406 + t543) * t494 + t518;
t367 = t406 * t474 + (-t407 + t542) * t475 + t522;
t366 = t383 * t458 - t391 * t425 + t440 * t494 + (-t469 + t541) * t474 + t520;
t365 = -t382 * t458 + t391 * t426 + t469 * t475 + (-t439 + t543) * t494 + t518;
t364 = t382 * t425 - t383 * t426 + t439 * t474 + (-t440 + t542) * t475 + t522;
t1 = m(6) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(7) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + t504 * ((t548 * t448 + (t450 * t515 + t452 * t511) * t507) * t503 - (t548 * t447 + (t449 * t515 + t451 * t511) * t507) * t534 + (t548 * t477 + (t478 * t515 + t479 * t511) * t507) * t504) / 0.2e1 + t425 * ((t468 * t377 + t429 * t379 + t430 * t381) * t425 + (t376 * t468 + t378 * t429 + t380 * t430) * t426 + (t388 * t468 + t389 * t429 + t390 * t430) * t458) / 0.2e1 + t426 * ((t377 * t466 + t379 * t427 + t381 * t428) * t425 + (t466 * t376 + t427 * t378 + t428 * t380) * t426 + (t388 * t466 + t389 * t427 + t390 * t428) * t458) / 0.2e1 + t458 * ((t377 * t481 + t379 * t463 + t381 * t464) * t425 + (t376 * t481 + t378 * t463 + t380 * t464) * t426 + (t481 * t388 + t463 * t389 + t464 * t390) * t458) / 0.2e1 + m(3) * (t385 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + m(4) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(5) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 - ((-t477 * t544 - t478 * t489 + t479 * t490) * t504 + ((-t450 * t489 + t452 * t490) * t512 + (t489 * t449 - t490 * t451 + t564) * t516) * t538) * t534 / 0.2e1 + ((t477 * t547 - t478 * t491 + t479 * t492) * t504 + (-(-t449 * t491 + t451 * t492) * t516 + (-t491 * t450 + t492 * t452 - t564) * t512) * t538) * t503 / 0.2e1 + (m(2) * (t499 ^ 2 + t500 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t443 * t472 + t444 * t473 + t467 * t559 + t468 * t558 + t491 * t555) * t494 + (t415 * t472 + t417 * t473 + t467 * t563 + t468 * t561 + t491 * t557) * t475 + (t472 * t416 + t473 * t418 + t467 * t562 + t468 * t560 + t491 * t556) * t474) * t474 / 0.2e1 + ((t443 * t470 + t444 * t471 + t465 * t559 + t466 * t558 + t489 * t555) * t494 + (t470 * t415 + t471 * t417 + t465 * t563 + t466 * t561 + t489 * t557) * t475 + (t416 * t470 + t418 * t471 + t465 * t562 + t466 * t560 + t489 * t556) * t474) * t475 / 0.2e1 + ((t487 * t443 + t488 * t444 + t480 * t559 + t481 * t558 - t545 * t555) * t494 + (t415 * t487 + t417 * t488 + t480 * t563 + t481 * t561 - t545 * t557) * t475 + (t416 * t487 + t418 * t488 + t480 * t562 + t481 * t560 - t545 * t556) * t474) * t494 / 0.2e1;
T  = t1;
