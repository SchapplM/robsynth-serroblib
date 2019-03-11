% Calculate kinetic energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:07
% EndTime: 2019-03-10 01:18:10
% DurationCPUTime: 3.31s
% Computational Cost: add. (2125->319), mult. (2762->501), div. (0->0), fcn. (2775->10), ass. (0->163)
t560 = Icges(6,1) + Icges(7,1);
t559 = Icges(6,4) + Icges(7,4);
t558 = -Icges(7,5) - Icges(6,5);
t557 = Icges(6,2) + Icges(7,2);
t556 = -Icges(7,6) - Icges(6,6);
t555 = -Icges(7,3) - Icges(6,3);
t483 = sin(qJ(1));
t486 = cos(qJ(1));
t482 = sin(qJ(2));
t485 = cos(qJ(2));
t480 = qJ(3) + qJ(4);
t474 = cos(t480);
t516 = pkin(4) * t474;
t492 = pkin(10) * t482 + t485 * t516;
t473 = sin(t480);
t508 = pkin(4) * t473;
t355 = t483 * t508 + t486 * t492;
t397 = -pkin(10) * t485 + t482 * t516;
t472 = qJD(2) * t483;
t512 = qJD(3) * t482;
t448 = t486 * t512 + t472;
t511 = qJD(4) * t482;
t418 = t486 * t511 + t448;
t509 = -qJD(3) - qJD(4);
t451 = t485 * t509 + qJD(1);
t554 = t451 * t355 - t397 * t418;
t354 = t483 * t492 - t486 * t508;
t513 = qJD(2) * t486;
t449 = t483 * t512 - t513;
t419 = t483 * t511 + t449;
t553 = -t354 * t451 + t419 * t397;
t476 = qJ(5) + t480;
t469 = sin(t476);
t470 = cos(t476);
t525 = t485 * t483;
t414 = -t469 * t525 - t470 * t486;
t415 = -t469 * t486 + t470 * t525;
t527 = t482 * t483;
t552 = -t556 * t414 - t558 * t415 - t555 * t527;
t524 = t485 * t486;
t416 = -t469 * t524 + t470 * t483;
t417 = t469 * t483 + t470 * t524;
t526 = t482 * t486;
t551 = -t556 * t416 - t558 * t417 - t555 * t526;
t550 = t557 * t414 + t559 * t415 - t556 * t527;
t549 = t557 * t416 + t559 * t417 - t556 * t526;
t548 = t559 * t414 + t560 * t415 - t558 * t527;
t547 = t559 * t416 + t560 * t417 - t558 * t526;
t484 = cos(qJ(3));
t537 = t484 * pkin(3);
t495 = pkin(9) * t482 + t485 * t537;
t481 = sin(qJ(3));
t529 = t481 * t483;
t396 = pkin(3) * t529 + t486 * t495;
t409 = -pkin(9) * t485 + t482 * t537;
t466 = -qJD(3) * t485 + qJD(1);
t546 = t466 * t396 - t409 * t448;
t528 = t481 * t486;
t395 = -pkin(3) * t528 + t483 * t495;
t545 = -t395 * t466 + t449 * t409;
t544 = t555 * t485 + (t556 * t469 - t558 * t470) * t482;
t543 = t556 * t485 + (-t557 * t469 + t559 * t470) * t482;
t542 = t558 * t485 + (-t559 * t469 + t560 * t470) * t482;
t535 = Icges(3,4) * t482;
t534 = Icges(3,4) * t485;
t520 = pkin(5) * t470;
t493 = qJ(6) * t482 + t485 * t520;
t518 = pkin(5) * t469;
t523 = rSges(7,1) * t415 + rSges(7,2) * t414 + rSges(7,3) * t527 + t483 * t493 - t486 * t518;
t522 = rSges(7,1) * t417 + rSges(7,2) * t416 + rSges(7,3) * t526 + t483 * t518 + t486 * t493;
t521 = (-qJ(6) - rSges(7,3)) * t485 + (rSges(7,1) * t470 - rSges(7,2) * t469 + t520) * t482;
t506 = pkin(2) * t485 + pkin(8) * t482;
t446 = t506 * t483;
t447 = t506 * t486;
t519 = t446 * t472 + t447 * t513;
t452 = qJD(1) * (pkin(1) * t486 + pkin(7) * t483);
t517 = qJD(1) * t447 + t452;
t510 = qJD(5) * t482;
t461 = pkin(1) * t483 - pkin(7) * t486;
t507 = (-t446 - t461) * qJD(1);
t505 = rSges(3,1) * t485 - rSges(3,2) * t482;
t504 = Icges(3,1) * t485 - t535;
t503 = -Icges(3,2) * t482 + t534;
t502 = Icges(3,5) * t485 - Icges(3,6) * t482;
t424 = -Icges(3,6) * t486 + t483 * t503;
t427 = -Icges(3,5) * t486 + t483 * t504;
t501 = t424 * t482 - t427 * t485;
t425 = Icges(3,6) * t483 + t486 * t503;
t428 = Icges(3,5) * t483 + t486 * t504;
t500 = -t425 * t482 + t428 * t485;
t455 = Icges(3,2) * t485 + t535;
t456 = Icges(3,1) * t482 + t534;
t499 = -t455 * t482 + t456 * t485;
t460 = pkin(2) * t482 - pkin(8) * t485;
t498 = -qJD(2) * t460 + qJD(6) * t482;
t497 = t448 * t395 - t396 * t449 + t519;
t496 = -t460 * t472 + t517;
t494 = -t460 * t513 + t507;
t491 = t418 * t354 - t355 * t419 + t497;
t490 = t496 + t546;
t489 = t494 + t545;
t459 = rSges(2,1) * t486 - rSges(2,2) * t483;
t458 = rSges(2,1) * t483 + rSges(2,2) * t486;
t457 = rSges(3,1) * t482 + rSges(3,2) * t485;
t454 = Icges(3,5) * t482 + Icges(3,6) * t485;
t445 = t484 * t524 + t529;
t444 = -t481 * t524 + t483 * t484;
t443 = t484 * t525 - t528;
t442 = -t481 * t525 - t484 * t486;
t441 = qJD(1) + (-qJD(5) + t509) * t485;
t436 = t473 * t483 + t474 * t524;
t435 = -t473 * t524 + t474 * t483;
t434 = -t473 * t486 + t474 * t525;
t433 = -t473 * t525 - t474 * t486;
t432 = rSges(3,3) * t483 + t486 * t505;
t431 = -rSges(3,3) * t486 + t483 * t505;
t430 = -rSges(4,3) * t485 + (rSges(4,1) * t484 - rSges(4,2) * t481) * t482;
t426 = -Icges(4,5) * t485 + (Icges(4,1) * t484 - Icges(4,4) * t481) * t482;
t423 = -Icges(4,6) * t485 + (Icges(4,4) * t484 - Icges(4,2) * t481) * t482;
t422 = Icges(3,3) * t483 + t486 * t502;
t421 = -Icges(3,3) * t486 + t483 * t502;
t420 = -Icges(4,3) * t485 + (Icges(4,5) * t484 - Icges(4,6) * t481) * t482;
t413 = -rSges(5,3) * t485 + (rSges(5,1) * t474 - rSges(5,2) * t473) * t482;
t412 = -Icges(5,5) * t485 + (Icges(5,1) * t474 - Icges(5,4) * t473) * t482;
t411 = -Icges(5,6) * t485 + (Icges(5,4) * t474 - Icges(5,2) * t473) * t482;
t410 = -Icges(5,3) * t485 + (Icges(5,5) * t474 - Icges(5,6) * t473) * t482;
t408 = -rSges(6,3) * t485 + (rSges(6,1) * t470 - rSges(6,2) * t469) * t482;
t400 = t483 * t510 + t419;
t399 = t486 * t510 + t418;
t394 = rSges(4,1) * t445 + rSges(4,2) * t444 + rSges(4,3) * t526;
t393 = rSges(4,1) * t443 + rSges(4,2) * t442 + rSges(4,3) * t527;
t392 = Icges(4,1) * t445 + Icges(4,4) * t444 + Icges(4,5) * t526;
t391 = Icges(4,1) * t443 + Icges(4,4) * t442 + Icges(4,5) * t527;
t390 = Icges(4,4) * t445 + Icges(4,2) * t444 + Icges(4,6) * t526;
t389 = Icges(4,4) * t443 + Icges(4,2) * t442 + Icges(4,6) * t527;
t388 = Icges(4,5) * t445 + Icges(4,6) * t444 + Icges(4,3) * t526;
t387 = Icges(4,5) * t443 + Icges(4,6) * t442 + Icges(4,3) * t527;
t386 = qJD(1) * t432 - t457 * t472 + t452;
t385 = -t457 * t513 + (-t431 - t461) * qJD(1);
t384 = (t431 * t483 + t432 * t486) * qJD(2);
t382 = rSges(5,1) * t436 + rSges(5,2) * t435 + rSges(5,3) * t526;
t381 = rSges(5,1) * t434 + rSges(5,2) * t433 + rSges(5,3) * t527;
t380 = Icges(5,1) * t436 + Icges(5,4) * t435 + Icges(5,5) * t526;
t379 = Icges(5,1) * t434 + Icges(5,4) * t433 + Icges(5,5) * t527;
t378 = Icges(5,4) * t436 + Icges(5,2) * t435 + Icges(5,6) * t526;
t377 = Icges(5,4) * t434 + Icges(5,2) * t433 + Icges(5,6) * t527;
t376 = Icges(5,5) * t436 + Icges(5,6) * t435 + Icges(5,3) * t526;
t375 = Icges(5,5) * t434 + Icges(5,6) * t433 + Icges(5,3) * t527;
t373 = rSges(6,1) * t417 + rSges(6,2) * t416 + rSges(6,3) * t526;
t371 = rSges(6,1) * t415 + rSges(6,2) * t414 + rSges(6,3) * t527;
t349 = t394 * t466 - t430 * t448 + t496;
t348 = -t393 * t466 + t430 * t449 + t494;
t347 = t393 * t448 - t394 * t449 + t519;
t346 = t382 * t451 - t413 * t418 + t490;
t345 = -t381 * t451 + t413 * t419 + t489;
t344 = t381 * t418 - t382 * t419 + t497;
t343 = t373 * t441 - t399 * t408 + t490 + t554;
t342 = -t371 * t441 + t400 * t408 + t489 + t553;
t341 = t371 * t399 - t373 * t400 + t491;
t340 = -t399 * t521 + t441 * t522 + t483 * t498 + t517 + t546 + t554;
t339 = t400 * t521 - t441 * t523 + t486 * t498 + t507 + t545 + t553;
t338 = -qJD(6) * t485 + t399 * t523 - t400 * t522 + t491;
t1 = -((-t486 * t454 + t483 * t499) * qJD(1) + (t486 ^ 2 * t421 + (t500 * t483 + (-t422 + t501) * t486) * t483) * qJD(2)) * t513 / 0.2e1 + ((t483 * t454 + t486 * t499) * qJD(1) + (t483 ^ 2 * t422 + (t501 * t486 + (-t421 + t500) * t483) * t486) * qJD(2)) * t472 / 0.2e1 + m(3) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(4) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(5) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(6) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(7) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + qJD(1) * ((t485 * t455 + t482 * t456) * qJD(1) + ((t425 * t485 + t428 * t482) * t483 - (t424 * t485 + t427 * t482) * t486) * qJD(2)) / 0.2e1 + t466 * ((-t387 * t449 - t388 * t448 - t420 * t466) * t485 + ((-t390 * t481 + t392 * t484) * t448 + (-t389 * t481 + t391 * t484) * t449 + (-t423 * t481 + t426 * t484) * t466) * t482) / 0.2e1 + t448 * ((t388 * t526 + t444 * t390 + t445 * t392) * t448 + (t387 * t526 + t389 * t444 + t391 * t445) * t449 + (t420 * t526 + t423 * t444 + t426 * t445) * t466) / 0.2e1 + t449 * ((t388 * t527 + t390 * t442 + t392 * t443) * t448 + (t387 * t527 + t442 * t389 + t443 * t391) * t449 + (t420 * t527 + t423 * t442 + t426 * t443) * t466) / 0.2e1 + t418 * ((t376 * t526 + t435 * t378 + t436 * t380) * t418 + (t375 * t526 + t377 * t435 + t379 * t436) * t419 + (t410 * t526 + t411 * t435 + t412 * t436) * t451) / 0.2e1 + t419 * ((t376 * t527 + t378 * t433 + t380 * t434) * t418 + (t375 * t527 + t433 * t377 + t434 * t379) * t419 + (t410 * t527 + t411 * t433 + t412 * t434) * t451) / 0.2e1 + t451 * ((-t375 * t419 - t376 * t418 - t410 * t451) * t485 + ((-t378 * t473 + t380 * t474) * t418 + (-t377 * t473 + t379 * t474) * t419 + (-t411 * t473 + t412 * t474) * t451) * t482) / 0.2e1 + ((t543 * t416 + t542 * t417 + t544 * t526) * t441 + (t550 * t416 + t548 * t417 + t552 * t526) * t400 + (t549 * t416 + t547 * t417 + t551 * t526) * t399) * t399 / 0.2e1 + ((t543 * t414 + t542 * t415 + t544 * t527) * t441 + (t550 * t414 + t548 * t415 + t552 * t527) * t400 + (t549 * t414 + t547 * t415 + t551 * t527) * t399) * t400 / 0.2e1 + ((-t551 * t399 - t552 * t400 - t544 * t441) * t485 + ((-t543 * t469 + t542 * t470) * t441 + (-t550 * t469 + t548 * t470) * t400 + (-t549 * t469 + t547 * t470) * t399) * t482) * t441 / 0.2e1 + (Icges(2,3) + m(2) * (t458 ^ 2 + t459 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
