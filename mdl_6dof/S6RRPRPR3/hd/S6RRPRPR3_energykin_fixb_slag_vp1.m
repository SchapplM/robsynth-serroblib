% Calculate kinetic energy for
% S6RRPRPR3
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:48
% EndTime: 2019-03-09 10:15:51
% DurationCPUTime: 2.95s
% Computational Cost: add. (2135->327), mult. (2199->497), div. (0->0), fcn. (2156->12), ass. (0->165)
t557 = Icges(3,3) + Icges(4,3);
t556 = -Icges(6,3) - Icges(5,3);
t474 = qJ(2) + pkin(10);
t465 = sin(t474);
t467 = cos(t474);
t478 = sin(qJ(2));
t481 = cos(qJ(2));
t555 = Icges(3,5) * t481 + Icges(4,5) * t467 - Icges(3,6) * t478 - Icges(4,6) * t465;
t473 = qJ(4) + pkin(11);
t464 = sin(t473);
t466 = cos(t473);
t482 = cos(qJ(1));
t479 = sin(qJ(1));
t526 = t467 * t479;
t418 = -t464 * t526 - t466 * t482;
t419 = -t464 * t482 + t466 * t526;
t480 = cos(qJ(4));
t521 = t480 * t482;
t477 = sin(qJ(4));
t524 = t477 * t479;
t433 = -t467 * t524 - t521;
t522 = t479 * t480;
t523 = t477 * t482;
t434 = t467 * t522 - t523;
t528 = t465 * t479;
t554 = Icges(5,5) * t434 + Icges(6,5) * t419 + Icges(5,6) * t433 + Icges(6,6) * t418 - t556 * t528;
t525 = t467 * t482;
t420 = -t464 * t525 + t466 * t479;
t421 = t464 * t479 + t466 * t525;
t435 = -t467 * t523 + t522;
t436 = t467 * t521 + t524;
t527 = t465 * t482;
t553 = Icges(5,5) * t436 + Icges(6,5) * t421 + Icges(5,6) * t435 + Icges(6,6) * t420 - t556 * t527;
t552 = t556 * t467 + (Icges(5,5) * t480 + Icges(6,5) * t466 - Icges(5,6) * t477 - Icges(6,6) * t464) * t465;
t551 = t555 * t479 - t557 * t482;
t550 = t557 * t479 + t555 * t482;
t549 = Icges(3,5) * t478 + Icges(4,5) * t465 + Icges(3,6) * t481 + Icges(4,6) * t467;
t530 = Icges(4,4) * t465;
t441 = Icges(4,2) * t467 + t530;
t529 = Icges(4,4) * t467;
t442 = Icges(4,1) * t465 + t529;
t532 = Icges(3,4) * t478;
t449 = Icges(3,2) * t481 + t532;
t531 = Icges(3,4) * t481;
t450 = Icges(3,1) * t478 + t531;
t548 = -t441 * t465 + t442 * t467 - t449 * t478 + t450 * t481;
t499 = -Icges(4,2) * t465 + t529;
t410 = -Icges(4,6) * t482 + t479 * t499;
t501 = Icges(4,1) * t467 - t530;
t412 = -Icges(4,5) * t482 + t479 * t501;
t500 = -Icges(3,2) * t478 + t531;
t425 = -Icges(3,6) * t482 + t479 * t500;
t502 = Icges(3,1) * t481 - t532;
t427 = -Icges(3,5) * t482 + t479 * t502;
t547 = t410 * t465 - t412 * t467 + t425 * t478 - t427 * t481;
t411 = Icges(4,6) * t479 + t482 * t499;
t413 = Icges(4,5) * t479 + t482 * t501;
t426 = Icges(3,6) * t479 + t482 * t500;
t428 = Icges(3,5) * t479 + t482 * t502;
t546 = -t411 * t465 + t413 * t467 - t426 * t478 + t428 * t481;
t539 = pkin(2) * t478;
t536 = pkin(2) * t481;
t535 = t480 * pkin(4);
t406 = -qJ(3) * t482 + t479 * t536;
t407 = qJ(3) * t479 + t482 * t536;
t470 = qJD(2) * t479;
t516 = qJD(2) * t482;
t520 = t406 * t470 + t407 * t516;
t458 = pkin(1) * t479 - pkin(7) * t482;
t519 = -t406 - t458;
t518 = pkin(5) * t466;
t515 = qJD(4) * t465;
t438 = t482 * t515 + t470;
t514 = qJD(5) * t465;
t513 = qJD(6) * t465;
t510 = pkin(5) * t464;
t439 = t479 * t515 - t516;
t508 = pkin(3) * t467 + pkin(8) * t465;
t431 = t508 * t479;
t432 = t508 * t482;
t509 = t431 * t470 + t432 * t516 + t520;
t447 = qJD(1) * (pkin(1) * t482 + pkin(7) * t479);
t507 = qJD(1) * t407 - qJD(3) * t482 + t447;
t506 = rSges(3,1) * t481 - rSges(3,2) * t478;
t505 = rSges(4,1) * t467 - rSges(4,2) * t465;
t504 = qJD(2) * (-rSges(4,1) * t465 - rSges(4,2) * t467 - t539);
t503 = qJD(2) * (-pkin(3) * t465 + pkin(8) * t467 - t539);
t490 = qJ(5) * t465 + t467 * t535;
t369 = -pkin(4) * t523 + t479 * t490;
t489 = -qJD(5) * t467 + t438 * t369 + t509;
t488 = pkin(9) * t465 + t467 * t518;
t487 = qJD(1) * t432 + t479 * t503 + t507;
t469 = qJD(3) * t479;
t486 = t469 + (-t431 + t519) * qJD(1) + t482 * t503;
t370 = pkin(4) * t524 + t482 * t490;
t459 = -qJD(4) * t467 + qJD(1);
t485 = t459 * t370 + t479 * t514 + t487;
t384 = -qJ(5) * t467 + t465 * t535;
t484 = t439 * t384 + t482 * t514 + t486;
t468 = qJ(6) + t473;
t461 = cos(t468);
t460 = sin(t468);
t457 = rSges(2,1) * t482 - rSges(2,2) * t479;
t456 = rSges(2,1) * t479 + rSges(2,2) * t482;
t455 = rSges(3,1) * t478 + rSges(3,2) * t481;
t437 = qJD(1) + (-qJD(4) - qJD(6)) * t467;
t430 = rSges(3,3) * t479 + t482 * t506;
t429 = -rSges(3,3) * t482 + t479 * t506;
t415 = rSges(4,3) * t479 + t482 * t505;
t414 = -rSges(4,3) * t482 + t479 * t505;
t405 = t460 * t479 + t461 * t525;
t404 = -t460 * t525 + t461 * t479;
t403 = -t460 * t482 + t461 * t526;
t402 = -t460 * t526 - t461 * t482;
t401 = t479 * t513 + t439;
t400 = t482 * t513 + t438;
t399 = -rSges(5,3) * t467 + (rSges(5,1) * t480 - rSges(5,2) * t477) * t465;
t397 = -Icges(5,5) * t467 + (Icges(5,1) * t480 - Icges(5,4) * t477) * t465;
t396 = -Icges(5,6) * t467 + (Icges(5,4) * t480 - Icges(5,2) * t477) * t465;
t392 = -rSges(6,3) * t467 + (rSges(6,1) * t466 - rSges(6,2) * t464) * t465;
t391 = -Icges(6,5) * t467 + (Icges(6,1) * t466 - Icges(6,4) * t464) * t465;
t390 = -Icges(6,6) * t467 + (Icges(6,4) * t466 - Icges(6,2) * t464) * t465;
t388 = -rSges(7,3) * t467 + (rSges(7,1) * t461 - rSges(7,2) * t460) * t465;
t387 = -Icges(7,5) * t467 + (Icges(7,1) * t461 - Icges(7,4) * t460) * t465;
t386 = -Icges(7,6) * t467 + (Icges(7,4) * t461 - Icges(7,2) * t460) * t465;
t385 = -Icges(7,3) * t467 + (Icges(7,5) * t461 - Icges(7,6) * t460) * t465;
t383 = qJD(1) * t430 - t455 * t470 + t447;
t382 = -t455 * t516 + (-t429 - t458) * qJD(1);
t381 = -pkin(9) * t467 + t465 * t518;
t380 = (t429 * t479 + t430 * t482) * qJD(2);
t379 = rSges(5,1) * t436 + rSges(5,2) * t435 + rSges(5,3) * t527;
t378 = rSges(5,1) * t434 + rSges(5,2) * t433 + rSges(5,3) * t528;
t377 = Icges(5,1) * t436 + Icges(5,4) * t435 + Icges(5,5) * t527;
t376 = Icges(5,1) * t434 + Icges(5,4) * t433 + Icges(5,5) * t528;
t375 = Icges(5,4) * t436 + Icges(5,2) * t435 + Icges(5,6) * t527;
t374 = Icges(5,4) * t434 + Icges(5,2) * t433 + Icges(5,6) * t528;
t368 = rSges(6,1) * t421 + rSges(6,2) * t420 + rSges(6,3) * t527;
t367 = rSges(6,1) * t419 + rSges(6,2) * t418 + rSges(6,3) * t528;
t366 = Icges(6,1) * t421 + Icges(6,4) * t420 + Icges(6,5) * t527;
t365 = Icges(6,1) * t419 + Icges(6,4) * t418 + Icges(6,5) * t528;
t364 = Icges(6,4) * t421 + Icges(6,2) * t420 + Icges(6,6) * t527;
t363 = Icges(6,4) * t419 + Icges(6,2) * t418 + Icges(6,6) * t528;
t359 = rSges(7,1) * t405 + rSges(7,2) * t404 + rSges(7,3) * t527;
t358 = rSges(7,1) * t403 + rSges(7,2) * t402 + rSges(7,3) * t528;
t357 = Icges(7,1) * t405 + Icges(7,4) * t404 + Icges(7,5) * t527;
t356 = Icges(7,1) * t403 + Icges(7,4) * t402 + Icges(7,5) * t528;
t355 = Icges(7,4) * t405 + Icges(7,2) * t404 + Icges(7,6) * t527;
t354 = Icges(7,4) * t403 + Icges(7,2) * t402 + Icges(7,6) * t528;
t353 = Icges(7,5) * t405 + Icges(7,6) * t404 + Icges(7,3) * t527;
t352 = Icges(7,5) * t403 + Icges(7,6) * t402 + Icges(7,3) * t528;
t350 = t479 * t510 + t482 * t488;
t349 = t479 * t488 - t482 * t510;
t348 = qJD(1) * t415 + t479 * t504 + t507;
t347 = t469 + t482 * t504 + (-t414 + t519) * qJD(1);
t346 = (t414 * t479 + t415 * t482) * qJD(2) + t520;
t345 = t379 * t459 - t399 * t438 + t487;
t344 = -t378 * t459 + t399 * t439 + t486;
t343 = t378 * t438 - t379 * t439 + t509;
t342 = t368 * t459 + (-t384 - t392) * t438 + t485;
t341 = t392 * t439 + (-t367 - t369) * t459 + t484;
t340 = t367 * t438 + (-t368 - t370) * t439 + t489;
t339 = t350 * t459 + t359 * t437 - t388 * t400 + (-t381 - t384) * t438 + t485;
t338 = -t358 * t437 + t381 * t439 + t388 * t401 + (-t349 - t369) * t459 + t484;
t337 = t349 * t438 + t358 * t400 - t359 * t401 + (-t350 - t370) * t439 + t489;
t1 = t401 * ((t353 * t528 + t355 * t402 + t357 * t403) * t400 + (t352 * t528 + t402 * t354 + t403 * t356) * t401 + (t385 * t528 + t386 * t402 + t387 * t403) * t437) / 0.2e1 + t437 * ((-t352 * t401 - t353 * t400 - t385 * t437) * t467 + ((-t355 * t460 + t357 * t461) * t400 + (-t354 * t460 + t356 * t461) * t401 + (-t386 * t460 + t387 * t461) * t437) * t465) / 0.2e1 + t400 * ((t353 * t527 + t404 * t355 + t405 * t357) * t400 + (t352 * t527 + t354 * t404 + t356 * t405) * t401 + (t385 * t527 + t386 * t404 + t387 * t405) * t437) / 0.2e1 + m(7) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(6) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(5) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(4) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(3) * (t380 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + ((t390 * t420 + t391 * t421 + t396 * t435 + t397 * t436 + t527 * t552) * t459 + (t363 * t420 + t365 * t421 + t374 * t435 + t376 * t436 + t527 * t554) * t439 + (t420 * t364 + t421 * t366 + t435 * t375 + t436 * t377 + t553 * t527) * t438) * t438 / 0.2e1 + ((t390 * t418 + t391 * t419 + t396 * t433 + t397 * t434 + t528 * t552) * t459 + (t418 * t363 + t419 * t365 + t433 * t374 + t434 * t376 + t554 * t528) * t439 + (t364 * t418 + t366 * t419 + t375 * t433 + t377 * t434 + t528 * t553) * t438) * t439 / 0.2e1 + ((-t553 * t438 - t439 * t554 - t552 * t459) * t467 + ((-t390 * t464 + t391 * t466 - t396 * t477 + t397 * t480) * t459 + (-t363 * t464 + t365 * t466 - t374 * t477 + t376 * t480) * t439 + (-t364 * t464 + t366 * t466 - t375 * t477 + t377 * t480) * t438) * t465) * t459 / 0.2e1 + (m(2) * (t456 ^ 2 + t457 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t410 * t467 - t412 * t465 - t425 * t481 - t427 * t478) * t482 + (t411 * t467 + t413 * t465 + t426 * t481 + t428 * t478) * t479) * qJD(2) + (t467 * t441 + t465 * t442 + t481 * t449 + t478 * t450) * qJD(1)) * qJD(1) / 0.2e1 + ((t550 * t479 ^ 2 + (t547 * t482 + (t546 - t551) * t479) * t482) * qJD(2) + (t479 * t549 + t482 * t548) * qJD(1)) * t470 / 0.2e1 - ((t551 * t482 ^ 2 + (t546 * t479 + (t547 - t550) * t482) * t479) * qJD(2) + (t479 * t548 - t482 * t549) * qJD(1)) * t516 / 0.2e1;
T  = t1;
