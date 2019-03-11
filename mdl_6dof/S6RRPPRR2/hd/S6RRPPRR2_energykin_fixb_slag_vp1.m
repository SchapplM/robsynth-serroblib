% Calculate kinetic energy for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:44
% EndTime: 2019-03-09 08:49:47
% DurationCPUTime: 3.06s
% Computational Cost: add. (2081->333), mult. (2154->509), div. (0->0), fcn. (2111->12), ass. (0->168)
t554 = Icges(3,3) + Icges(4,3);
t475 = qJ(2) + pkin(10);
t466 = sin(t475);
t468 = cos(t475);
t480 = sin(qJ(2));
t482 = cos(qJ(2));
t553 = Icges(3,5) * t482 + Icges(4,5) * t468 - Icges(3,6) * t480 - Icges(4,6) * t466;
t481 = sin(qJ(1));
t483 = cos(qJ(1));
t552 = t553 * t481 - t554 * t483;
t551 = t554 * t481 + t553 * t483;
t550 = Icges(3,5) * t480 + Icges(4,5) * t466 + Icges(3,6) * t482 + Icges(4,6) * t468;
t535 = Icges(4,4) * t466;
t442 = Icges(4,2) * t468 + t535;
t534 = Icges(4,4) * t468;
t443 = Icges(4,1) * t466 + t534;
t537 = Icges(3,4) * t480;
t450 = Icges(3,2) * t482 + t537;
t536 = Icges(3,4) * t482;
t451 = Icges(3,1) * t480 + t536;
t549 = -t442 * t466 + t443 * t468 - t450 * t480 + t451 * t482;
t502 = -Icges(4,2) * t466 + t534;
t412 = Icges(4,6) * t481 + t483 * t502;
t504 = Icges(4,1) * t468 - t535;
t414 = Icges(4,5) * t481 + t483 * t504;
t503 = -Icges(3,2) * t480 + t536;
t427 = Icges(3,6) * t481 + t483 * t503;
t505 = Icges(3,1) * t482 - t537;
t429 = Icges(3,5) * t481 + t483 * t505;
t548 = -t412 * t466 + t414 * t468 - t427 * t480 + t429 * t482;
t411 = -Icges(4,6) * t483 + t481 * t502;
t413 = -Icges(4,5) * t483 + t481 * t504;
t426 = -Icges(3,6) * t483 + t481 * t503;
t428 = -Icges(3,5) * t483 + t481 * t505;
t547 = t411 * t466 - t413 * t468 + t426 * t480 - t428 * t482;
t543 = pkin(2) * t480;
t540 = pkin(2) * t482;
t477 = cos(pkin(11));
t539 = t477 * pkin(4);
t533 = t466 * t481;
t532 = t466 * t483;
t531 = t468 * t481;
t530 = t468 * t483;
t476 = sin(pkin(11));
t529 = t476 * t481;
t528 = t476 * t483;
t527 = t477 * t481;
t526 = t477 * t483;
t407 = -qJ(3) * t483 + t481 * t540;
t408 = qJ(3) * t481 + t483 * t540;
t471 = qJD(2) * t481;
t519 = qJD(2) * t483;
t524 = t407 * t471 + t408 * t519;
t459 = pkin(1) * t481 - pkin(7) * t483;
t523 = -t407 - t459;
t474 = pkin(11) + qJ(5);
t467 = cos(t474);
t522 = pkin(5) * t467;
t470 = qJD(3) * t481;
t518 = qJD(4) * t466;
t521 = t483 * t518 + t470;
t517 = qJD(5) * t466;
t439 = t483 * t517 + t471;
t516 = qJD(6) * t466;
t507 = pkin(3) * t468 + qJ(4) * t466;
t430 = t507 * t481;
t515 = -t430 + t523;
t512 = -pkin(3) * t466 + qJ(4) * t468 - t543;
t465 = sin(t474);
t511 = pkin(5) * t465;
t440 = t481 * t517 - t519;
t448 = qJD(1) * (pkin(1) * t483 + pkin(7) * t481);
t510 = qJD(1) * t408 - qJD(3) * t483 + t448;
t509 = rSges(3,1) * t482 - rSges(3,2) * t480;
t508 = rSges(4,1) * t468 - rSges(4,2) * t466;
t506 = qJD(2) * (-rSges(4,1) * t466 - rSges(4,2) * t468 - t543);
t493 = qJD(2) * (pkin(8) * t468 - t466 * t539 + t512);
t492 = qJD(2) * (rSges(5,3) * t468 - (rSges(5,1) * t477 - rSges(5,2) * t476) * t466 + t512);
t431 = t507 * t483;
t491 = qJD(1) * t431 + t481 * t518 + t510;
t490 = -qJD(4) * t468 + t430 * t471 + t431 * t519 + t524;
t489 = pkin(8) * t466 + t468 * t539;
t488 = pkin(9) * t466 + t468 * t522;
t371 = -pkin(4) * t528 + t481 * t489;
t372 = pkin(4) * t529 + t483 * t489;
t487 = t371 * t471 + t372 * t519 + t490;
t486 = qJD(1) * t372 + t481 * t493 + t491;
t485 = (-t371 + t515) * qJD(1) + t483 * t493 + t521;
t469 = qJ(6) + t474;
t462 = cos(t469);
t461 = sin(t469);
t460 = -qJD(5) * t468 + qJD(1);
t458 = rSges(2,1) * t483 - rSges(2,2) * t481;
t457 = rSges(2,1) * t481 + rSges(2,2) * t483;
t456 = rSges(3,1) * t480 + rSges(3,2) * t482;
t438 = qJD(1) + (-qJD(5) - qJD(6)) * t468;
t437 = t468 * t526 + t529;
t436 = -t468 * t528 + t527;
t435 = t468 * t527 - t528;
t434 = -t468 * t529 - t526;
t433 = rSges(3,3) * t481 + t483 * t509;
t432 = -rSges(3,3) * t483 + t481 * t509;
t422 = t465 * t481 + t467 * t530;
t421 = -t465 * t530 + t467 * t481;
t420 = -t465 * t483 + t467 * t531;
t419 = -t465 * t531 - t467 * t483;
t418 = rSges(4,3) * t481 + t483 * t508;
t417 = -rSges(4,3) * t483 + t481 * t508;
t406 = t461 * t481 + t462 * t530;
t405 = -t461 * t530 + t462 * t481;
t404 = -t461 * t483 + t462 * t531;
t403 = -t461 * t531 - t462 * t483;
t402 = t481 * t516 + t440;
t401 = t483 * t516 + t439;
t398 = -Icges(5,5) * t468 + (Icges(5,1) * t477 - Icges(5,4) * t476) * t466;
t397 = -Icges(5,6) * t468 + (Icges(5,4) * t477 - Icges(5,2) * t476) * t466;
t396 = -Icges(5,3) * t468 + (Icges(5,5) * t477 - Icges(5,6) * t476) * t466;
t393 = -rSges(6,3) * t468 + (rSges(6,1) * t467 - rSges(6,2) * t465) * t466;
t392 = -Icges(6,5) * t468 + (Icges(6,1) * t467 - Icges(6,4) * t465) * t466;
t391 = -Icges(6,6) * t468 + (Icges(6,4) * t467 - Icges(6,2) * t465) * t466;
t390 = -Icges(6,3) * t468 + (Icges(6,5) * t467 - Icges(6,6) * t465) * t466;
t389 = -rSges(7,3) * t468 + (rSges(7,1) * t462 - rSges(7,2) * t461) * t466;
t388 = -Icges(7,5) * t468 + (Icges(7,1) * t462 - Icges(7,4) * t461) * t466;
t387 = -Icges(7,6) * t468 + (Icges(7,4) * t462 - Icges(7,2) * t461) * t466;
t386 = -Icges(7,3) * t468 + (Icges(7,5) * t462 - Icges(7,6) * t461) * t466;
t384 = qJD(1) * t433 - t456 * t471 + t448;
t383 = -t456 * t519 + (-t432 - t459) * qJD(1);
t382 = (t432 * t481 + t433 * t483) * qJD(2);
t381 = -pkin(9) * t468 + t466 * t522;
t380 = rSges(5,1) * t437 + rSges(5,2) * t436 + rSges(5,3) * t532;
t379 = rSges(5,1) * t435 + rSges(5,2) * t434 + rSges(5,3) * t533;
t378 = Icges(5,1) * t437 + Icges(5,4) * t436 + Icges(5,5) * t532;
t377 = Icges(5,1) * t435 + Icges(5,4) * t434 + Icges(5,5) * t533;
t376 = Icges(5,4) * t437 + Icges(5,2) * t436 + Icges(5,6) * t532;
t375 = Icges(5,4) * t435 + Icges(5,2) * t434 + Icges(5,6) * t533;
t374 = Icges(5,5) * t437 + Icges(5,6) * t436 + Icges(5,3) * t532;
t373 = Icges(5,5) * t435 + Icges(5,6) * t434 + Icges(5,3) * t533;
t367 = rSges(6,1) * t422 + rSges(6,2) * t421 + rSges(6,3) * t532;
t366 = rSges(6,1) * t420 + rSges(6,2) * t419 + rSges(6,3) * t533;
t365 = Icges(6,1) * t422 + Icges(6,4) * t421 + Icges(6,5) * t532;
t364 = Icges(6,1) * t420 + Icges(6,4) * t419 + Icges(6,5) * t533;
t363 = Icges(6,4) * t422 + Icges(6,2) * t421 + Icges(6,6) * t532;
t362 = Icges(6,4) * t420 + Icges(6,2) * t419 + Icges(6,6) * t533;
t361 = Icges(6,5) * t422 + Icges(6,6) * t421 + Icges(6,3) * t532;
t360 = Icges(6,5) * t420 + Icges(6,6) * t419 + Icges(6,3) * t533;
t359 = rSges(7,1) * t406 + rSges(7,2) * t405 + rSges(7,3) * t532;
t358 = rSges(7,1) * t404 + rSges(7,2) * t403 + rSges(7,3) * t533;
t357 = Icges(7,1) * t406 + Icges(7,4) * t405 + Icges(7,5) * t532;
t356 = Icges(7,1) * t404 + Icges(7,4) * t403 + Icges(7,5) * t533;
t355 = Icges(7,4) * t406 + Icges(7,2) * t405 + Icges(7,6) * t532;
t354 = Icges(7,4) * t404 + Icges(7,2) * t403 + Icges(7,6) * t533;
t353 = Icges(7,5) * t406 + Icges(7,6) * t405 + Icges(7,3) * t532;
t352 = Icges(7,5) * t404 + Icges(7,6) * t403 + Icges(7,3) * t533;
t351 = t481 * t511 + t483 * t488;
t350 = t481 * t488 - t483 * t511;
t349 = qJD(1) * t418 + t481 * t506 + t510;
t348 = t470 + t483 * t506 + (-t417 + t523) * qJD(1);
t347 = (t417 * t481 + t418 * t483) * qJD(2) + t524;
t346 = qJD(1) * t380 + t481 * t492 + t491;
t345 = t483 * t492 + (-t379 + t515) * qJD(1) + t521;
t344 = (t379 * t481 + t380 * t483) * qJD(2) + t490;
t343 = t367 * t460 - t393 * t439 + t486;
t342 = -t366 * t460 + t393 * t440 + t485;
t341 = t366 * t439 - t367 * t440 + t487;
t340 = t351 * t460 + t359 * t438 - t381 * t439 - t389 * t401 + t486;
t339 = -t350 * t460 - t358 * t438 + t381 * t440 + t389 * t402 + t485;
t338 = t350 * t439 - t351 * t440 + t358 * t401 - t359 * t402 + t487;
t1 = t401 * ((t353 * t532 + t405 * t355 + t406 * t357) * t401 + (t352 * t532 + t354 * t405 + t356 * t406) * t402 + (t386 * t532 + t387 * t405 + t388 * t406) * t438) / 0.2e1 + t438 * ((-t352 * t402 - t353 * t401 - t386 * t438) * t468 + ((-t355 * t461 + t357 * t462) * t401 + (-t354 * t461 + t356 * t462) * t402 + (-t387 * t461 + t388 * t462) * t438) * t466) / 0.2e1 + t460 * ((-t360 * t440 - t361 * t439 - t390 * t460) * t468 + ((-t363 * t465 + t365 * t467) * t439 + (-t362 * t465 + t364 * t467) * t440 + (-t391 * t465 + t392 * t467) * t460) * t466) / 0.2e1 + t402 * ((t353 * t533 + t355 * t403 + t357 * t404) * t401 + (t352 * t533 + t403 * t354 + t404 * t356) * t402 + (t386 * t533 + t387 * t403 + t388 * t404) * t438) / 0.2e1 + t439 * ((t361 * t532 + t421 * t363 + t422 * t365) * t439 + (t360 * t532 + t362 * t421 + t364 * t422) * t440 + (t390 * t532 + t391 * t421 + t392 * t422) * t460) / 0.2e1 + t440 * ((t361 * t533 + t363 * t419 + t365 * t420) * t439 + (t360 * t533 + t419 * t362 + t420 * t364) * t440 + (t390 * t533 + t391 * t419 + t392 * t420) * t460) / 0.2e1 + m(7) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(6) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(4) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(5) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 + (m(2) * (t457 ^ 2 + t458 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t426 * t482 - t428 * t480 + (t373 - t411) * t468 + (t375 * t476 - t377 * t477 - t413) * t466) * t483 + (t427 * t482 + t429 * t480 + (-t374 + t412) * t468 + (-t376 * t476 + t378 * t477 + t414) * t466) * t481) * qJD(2) + (t482 * t450 + t480 * t451 + (-t396 + t442) * t468 + (-t397 * t476 + t398 * t477 + t443) * t466) * qJD(1)) * qJD(1) / 0.2e1 + (((-t373 * t532 - t375 * t436 - t377 * t437 + t547 * t483) * t483 + (t374 * t532 + t376 * t436 + t378 * t437 + (t548 - t552) * t483 + t551 * t481) * t481) * qJD(2) + (t396 * t532 + t397 * t436 + t398 * t437 + t481 * t550 + t483 * t549) * qJD(1)) * t471 / 0.2e1 - (((t374 * t533 + t376 * t434 + t378 * t435 + t548 * t481) * t481 + (-t373 * t533 - t375 * t434 - t377 * t435 + (t547 - t551) * t481 + t552 * t483) * t483) * qJD(2) + (t396 * t533 + t397 * t434 + t398 * t435 + t481 * t549 - t483 * t550) * qJD(1)) * t519 / 0.2e1;
T  = t1;
