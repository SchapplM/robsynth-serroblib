% Calculate kinetic energy for
% S6RRRPPR9
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:36
% EndTime: 2019-03-09 16:11:39
% DurationCPUTime: 2.65s
% Computational Cost: add. (2807->332), mult. (7107->484), div. (0->0), fcn. (8873->12), ass. (0->149)
t560 = Icges(5,1) + Icges(6,1);
t559 = -Icges(5,4) + Icges(6,5);
t558 = Icges(6,4) + Icges(5,5);
t557 = Icges(5,2) + Icges(6,3);
t556 = -Icges(5,6) + Icges(6,6);
t555 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t510 = sin(qJ(2));
t511 = sin(qJ(1));
t513 = cos(qJ(2));
t514 = cos(qJ(1));
t539 = cos(pkin(6));
t525 = t514 * t539;
t490 = t510 * t511 - t513 * t525;
t491 = t510 * t525 + t511 * t513;
t507 = sin(pkin(6));
t535 = t507 * t514;
t453 = Icges(3,5) * t491 - Icges(3,6) * t490 - Icges(3,3) * t535;
t526 = t511 * t539;
t492 = t514 * t510 + t513 * t526;
t493 = -t510 * t526 + t514 * t513;
t537 = t507 * t511;
t454 = Icges(3,5) * t493 - Icges(3,6) * t492 + Icges(3,3) * t537;
t554 = t507 * (t453 * t514 - t454 * t511);
t509 = sin(qJ(3));
t540 = cos(qJ(3));
t474 = t491 * t540 - t509 * t535;
t506 = sin(pkin(11));
t538 = cos(pkin(11));
t444 = t474 * t506 - t490 * t538;
t445 = t474 * t538 + t490 * t506;
t528 = t507 * t540;
t473 = t491 * t509 + t514 * t528;
t553 = t557 * t444 + t559 * t445 + t556 * t473;
t476 = t493 * t540 + t509 * t537;
t446 = t476 * t506 - t492 * t538;
t447 = t476 * t538 + t492 * t506;
t475 = t493 * t509 - t511 * t528;
t552 = t557 * t446 + t559 * t447 + t556 * t475;
t551 = t559 * t444 + t560 * t445 + t558 * t473;
t550 = t559 * t446 + t560 * t447 + t558 * t475;
t489 = t509 * t539 + t510 * t528;
t536 = t507 * t513;
t471 = t489 * t506 + t536 * t538;
t472 = t489 * t538 - t506 * t536;
t488 = t507 * t509 * t510 - t539 * t540;
t549 = t557 * t471 + t559 * t472 + t556 * t488;
t548 = t559 * t471 + t560 * t472 + t558 * t488;
t547 = -Icges(4,4) * t474 - Icges(4,6) * t490 + t556 * t444 + t558 * t445 + t555 * t473;
t546 = -Icges(4,4) * t476 - Icges(4,6) * t492 + t556 * t446 + t558 * t447 + t555 * t475;
t545 = -Icges(4,4) * t489 + Icges(4,6) * t536 + t556 * t471 + t558 * t472 + t555 * t488;
t408 = pkin(4) * t445 + qJ(5) * t444;
t438 = pkin(3) * t474 + qJ(4) * t473;
t534 = -t408 - t438;
t409 = pkin(4) * t447 + qJ(5) * t446;
t439 = pkin(3) * t476 + qJ(4) * t475;
t533 = -t409 - t439;
t436 = pkin(4) * t472 + qJ(5) * t471;
t464 = pkin(3) * t489 + qJ(4) * t488;
t532 = -t436 - t464;
t465 = pkin(2) * t491 + pkin(9) * t490;
t466 = pkin(2) * t493 + pkin(9) * t492;
t529 = qJD(2) * t507;
t502 = t511 * t529;
t527 = t514 * t529;
t531 = t465 * t502 + t466 * t527;
t477 = qJD(3) * t492 + t502;
t530 = qJD(1) * (pkin(1) * t511 - pkin(8) * t535);
t503 = qJD(2) * t539 + qJD(1);
t524 = qJD(4) * t488 + t477 * t438 + t531;
t478 = qJD(3) * t490 - t527;
t495 = -qJD(3) * t536 + t503;
t522 = qJD(5) * t471 + t477 * t408 + t524;
t494 = (pkin(2) * t510 - pkin(9) * t513) * t507;
t496 = qJD(1) * (pkin(1) * t514 + pkin(8) * t537);
t521 = t503 * t466 - t494 * t502 + t496;
t520 = qJD(4) * t473 + t495 * t439 + t521;
t519 = -t465 * t503 - t494 * t527 - t530;
t518 = qJD(5) * t444 + t495 * t409 + t520;
t517 = qJD(4) * t475 + t478 * t464 + t519;
t516 = qJD(5) * t446 + t478 * t436 + t517;
t512 = cos(qJ(6));
t508 = sin(qJ(6));
t499 = rSges(2,1) * t514 - rSges(2,2) * t511;
t498 = rSges(2,1) * t511 + rSges(2,2) * t514;
t482 = t539 * rSges(3,3) + (rSges(3,1) * t510 + rSges(3,2) * t513) * t507;
t481 = Icges(3,5) * t539 + (Icges(3,1) * t510 + Icges(3,4) * t513) * t507;
t480 = Icges(3,6) * t539 + (Icges(3,4) * t510 + Icges(3,2) * t513) * t507;
t479 = Icges(3,3) * t539 + (Icges(3,5) * t510 + Icges(3,6) * t513) * t507;
t467 = -qJD(6) * t488 + t495;
t461 = rSges(3,1) * t493 - rSges(3,2) * t492 + rSges(3,3) * t537;
t460 = rSges(3,1) * t491 - rSges(3,2) * t490 - rSges(3,3) * t535;
t458 = Icges(3,1) * t493 - Icges(3,4) * t492 + Icges(3,5) * t537;
t457 = Icges(3,1) * t491 - Icges(3,4) * t490 - Icges(3,5) * t535;
t456 = Icges(3,4) * t493 - Icges(3,2) * t492 + Icges(3,6) * t537;
t455 = Icges(3,4) * t491 - Icges(3,2) * t490 - Icges(3,6) * t535;
t452 = rSges(4,1) * t489 - rSges(4,2) * t488 - rSges(4,3) * t536;
t451 = Icges(4,1) * t489 - Icges(4,4) * t488 - Icges(4,5) * t536;
t449 = Icges(4,5) * t489 - Icges(4,6) * t488 - Icges(4,3) * t536;
t448 = pkin(5) * t472 - pkin(10) * t488;
t441 = -qJD(6) * t473 + t478;
t440 = -qJD(6) * t475 + t477;
t435 = t471 * t508 + t472 * t512;
t434 = t471 * t512 - t472 * t508;
t432 = rSges(4,1) * t476 - rSges(4,2) * t475 + rSges(4,3) * t492;
t431 = rSges(4,1) * t474 - rSges(4,2) * t473 + rSges(4,3) * t490;
t430 = Icges(4,1) * t476 - Icges(4,4) * t475 + Icges(4,5) * t492;
t429 = Icges(4,1) * t474 - Icges(4,4) * t473 + Icges(4,5) * t490;
t426 = Icges(4,5) * t476 - Icges(4,6) * t475 + Icges(4,3) * t492;
t425 = Icges(4,5) * t474 - Icges(4,6) * t473 + Icges(4,3) * t490;
t424 = rSges(5,1) * t472 - rSges(5,2) * t471 + rSges(5,3) * t488;
t423 = rSges(6,1) * t472 + rSges(6,2) * t488 + rSges(6,3) * t471;
t416 = pkin(5) * t447 - pkin(10) * t475;
t415 = pkin(5) * t445 - pkin(10) * t473;
t412 = t461 * t503 - t482 * t502 + t496;
t411 = -t460 * t503 - t482 * t527 - t530;
t410 = (t460 * t511 + t461 * t514) * t529;
t407 = t446 * t508 + t447 * t512;
t406 = t446 * t512 - t447 * t508;
t405 = t444 * t508 + t445 * t512;
t404 = t444 * t512 - t445 * t508;
t401 = rSges(5,1) * t447 - rSges(5,2) * t446 + rSges(5,3) * t475;
t400 = rSges(6,1) * t447 + rSges(6,2) * t475 + rSges(6,3) * t446;
t399 = rSges(5,1) * t445 - rSges(5,2) * t444 + rSges(5,3) * t473;
t398 = rSges(6,1) * t445 + rSges(6,2) * t473 + rSges(6,3) * t444;
t385 = rSges(7,1) * t435 + rSges(7,2) * t434 - rSges(7,3) * t488;
t384 = Icges(7,1) * t435 + Icges(7,4) * t434 - Icges(7,5) * t488;
t383 = Icges(7,4) * t435 + Icges(7,2) * t434 - Icges(7,6) * t488;
t382 = Icges(7,5) * t435 + Icges(7,6) * t434 - Icges(7,3) * t488;
t381 = t432 * t495 - t452 * t477 + t521;
t380 = -t431 * t495 + t452 * t478 + t519;
t379 = rSges(7,1) * t407 + rSges(7,2) * t406 - rSges(7,3) * t475;
t378 = rSges(7,1) * t405 + rSges(7,2) * t404 - rSges(7,3) * t473;
t377 = Icges(7,1) * t407 + Icges(7,4) * t406 - Icges(7,5) * t475;
t376 = Icges(7,1) * t405 + Icges(7,4) * t404 - Icges(7,5) * t473;
t375 = Icges(7,4) * t407 + Icges(7,2) * t406 - Icges(7,6) * t475;
t374 = Icges(7,4) * t405 + Icges(7,2) * t404 - Icges(7,6) * t473;
t373 = Icges(7,5) * t407 + Icges(7,6) * t406 - Icges(7,3) * t475;
t372 = Icges(7,5) * t405 + Icges(7,6) * t404 - Icges(7,3) * t473;
t371 = t431 * t477 - t432 * t478 + t531;
t370 = t401 * t495 + (-t424 - t464) * t477 + t520;
t369 = t424 * t478 + (-t399 - t438) * t495 + t517;
t368 = t399 * t477 + (-t401 - t439) * t478 + t524;
t367 = t400 * t495 + (-t423 + t532) * t477 + t518;
t366 = t423 * t478 + (-t398 + t534) * t495 + t516;
t365 = t398 * t477 + (-t400 + t533) * t478 + t522;
t364 = t379 * t467 - t385 * t440 + t416 * t495 + (-t448 + t532) * t477 + t518;
t363 = -t378 * t467 + t385 * t441 + t448 * t478 + (-t415 + t534) * t495 + t516;
t362 = t378 * t440 - t379 * t441 + t415 * t477 + (-t416 + t533) * t478 + t522;
t1 = t441 * ((-t373 * t473 + t375 * t404 + t377 * t405) * t440 + (-t473 * t372 + t404 * t374 + t405 * t376) * t441 + (-t382 * t473 + t383 * t404 + t384 * t405) * t467) / 0.2e1 + t467 * ((-t373 * t488 + t375 * t434 + t377 * t435) * t440 + (-t372 * t488 + t374 * t434 + t376 * t435) * t441 + (-t382 * t488 + t383 * t434 + t384 * t435) * t467) / 0.2e1 + m(3) * (t410 ^ 2 + t411 ^ 2 + t412 ^ 2) / 0.2e1 + m(4) * (t371 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(5) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(6) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(7) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + t503 * ((t539 * t454 + (t456 * t513 + t458 * t510) * t507) * t502 - (t539 * t453 + (t455 * t513 + t457 * t510) * t507) * t527 + (t539 * t479 + (t480 * t513 + t481 * t510) * t507) * t503) / 0.2e1 + t440 * ((-t475 * t373 + t406 * t375 + t407 * t377) * t440 + (-t372 * t475 + t374 * t406 + t376 * t407) * t441 + (-t382 * t475 + t383 * t406 + t384 * t407) * t467) / 0.2e1 - ((-t479 * t535 - t480 * t490 + t481 * t491) * t503 + ((-t456 * t490 + t458 * t491) * t511 + (t490 * t455 - t491 * t457 + t554) * t514) * t529) * t527 / 0.2e1 + ((t479 * t537 - t480 * t492 + t481 * t493) * t503 + (-(-t455 * t492 + t457 * t493) * t514 + (-t492 * t456 + t493 * t458 - t554) * t511) * t529) * t502 / 0.2e1 + (Icges(2,3) + m(2) * (t498 ^ 2 + t499 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t446 * t549 + t447 * t548 + t449 * t492 + t451 * t476 + t475 * t545) * t495 + (t425 * t492 + t429 * t476 + t446 * t553 + t447 * t551 + t475 * t547) * t478 + (t426 * t492 + t430 * t476 + t446 * t552 + t447 * t550 + t475 * t546) * t477) * t477 / 0.2e1 + ((t444 * t549 + t445 * t548 + t449 * t490 + t451 * t474 + t473 * t545) * t495 + (t425 * t490 + t429 * t474 + t444 * t553 + t445 * t551 + t473 * t547) * t478 + (t426 * t490 + t430 * t474 + t444 * t552 + t445 * t550 + t473 * t546) * t477) * t478 / 0.2e1 + ((-t449 * t536 + t451 * t489 + t471 * t549 + t472 * t548 + t488 * t545) * t495 + (-t425 * t536 + t429 * t489 + t471 * t553 + t472 * t551 + t547 * t488) * t478 + (-t426 * t536 + t430 * t489 + t471 * t552 + t472 * t550 + t488 * t546) * t477) * t495 / 0.2e1;
T  = t1;
