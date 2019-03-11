% Calculate kinetic energy for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:13
% EndTime: 2019-03-09 06:14:15
% DurationCPUTime: 2.25s
% Computational Cost: add. (1913->258), mult. (1987->403), div. (0->0), fcn. (1977->10), ass. (0->139)
t535 = Icges(6,1) + Icges(7,1);
t534 = Icges(6,4) + Icges(7,4);
t533 = -Icges(7,5) - Icges(6,5);
t532 = Icges(6,2) + Icges(7,2);
t531 = -Icges(7,6) - Icges(6,6);
t530 = -Icges(7,3) - Icges(6,3);
t464 = cos(qJ(1));
t456 = pkin(10) + qJ(3);
t448 = sin(t456);
t449 = cos(t456);
t463 = cos(qJ(4));
t512 = t463 * pkin(4);
t470 = pkin(9) * t448 + t449 * t512;
t461 = sin(qJ(4));
t462 = sin(qJ(1));
t500 = t461 * t462;
t378 = pkin(4) * t500 + t464 * t470;
t391 = -pkin(9) * t449 + t448 * t512;
t450 = qJD(3) * t462;
t487 = qJD(4) * t448;
t430 = t464 * t487 + t450;
t445 = -qJD(4) * t449 + qJD(1);
t529 = t445 * t378 - t391 * t430;
t457 = qJ(4) + qJ(5);
t453 = cos(t457);
t501 = t453 * t464;
t452 = sin(t457);
t504 = t452 * t462;
t418 = -t449 * t504 - t501;
t502 = t453 * t462;
t503 = t452 * t464;
t419 = t449 * t502 - t503;
t506 = t448 * t462;
t528 = -t531 * t418 - t533 * t419 - t530 * t506;
t420 = -t449 * t503 + t502;
t421 = t449 * t501 + t504;
t505 = t448 * t464;
t527 = -t531 * t420 - t533 * t421 - t530 * t505;
t526 = t532 * t418 + t534 * t419 - t531 * t506;
t525 = t532 * t420 + t534 * t421 - t531 * t505;
t524 = t534 * t418 + t535 * t419 - t533 * t506;
t523 = t534 * t420 + t535 * t421 - t533 * t505;
t499 = t461 * t464;
t377 = -pkin(4) * t499 + t462 * t470;
t488 = qJD(3) * t464;
t431 = t462 * t487 - t488;
t522 = -t377 * t445 + t431 * t391;
t521 = t530 * t449 + (t531 * t452 - t533 * t453) * t448;
t520 = t531 * t449 + (-t532 * t452 + t534 * t453) * t448;
t519 = t533 * t449 + (-t534 * t452 + t535 * t453) * t448;
t459 = cos(pkin(10));
t513 = pkin(2) * t459;
t510 = Icges(4,4) * t448;
t509 = Icges(4,4) * t449;
t498 = t462 * t463;
t497 = t463 * t464;
t490 = pkin(5) * t453;
t469 = qJ(6) * t448 + t449 * t490;
t485 = pkin(5) * t452;
t495 = rSges(7,1) * t419 + rSges(7,2) * t418 + rSges(7,3) * t506 + t462 * t469 - t464 * t485;
t494 = rSges(7,1) * t421 + rSges(7,2) * t420 + rSges(7,3) * t505 + t462 * t485 + t464 * t469;
t493 = (-qJ(6) - rSges(7,3)) * t449 + (rSges(7,1) * t453 - rSges(7,2) * t452 + t490) * t448;
t442 = pkin(1) * t462 - qJ(2) * t464;
t492 = pkin(7) * t464 - t462 * t513 - t442;
t484 = pkin(3) * t449 + pkin(8) * t448;
t423 = t484 * t462;
t424 = t484 * t464;
t491 = t423 * t450 + t424 * t488;
t486 = qJD(5) * t448;
t438 = qJD(1) * (pkin(1) * t464 + qJ(2) * t462);
t483 = -qJD(2) * t464 + qJD(1) * (pkin(7) * t462 + t464 * t513) + t438;
t458 = sin(pkin(10));
t482 = rSges(3,1) * t459 - rSges(3,2) * t458;
t481 = rSges(4,1) * t449 - rSges(4,2) * t448;
t480 = Icges(4,1) * t449 - t510;
t479 = -Icges(4,2) * t448 + t509;
t478 = Icges(4,5) * t449 - Icges(4,6) * t448;
t410 = -Icges(4,6) * t464 + t462 * t479;
t412 = -Icges(4,5) * t464 + t462 * t480;
t477 = t410 * t448 - t412 * t449;
t411 = Icges(4,6) * t462 + t464 * t479;
t413 = Icges(4,5) * t462 + t464 * t480;
t476 = -t411 * t448 + t413 * t449;
t433 = Icges(4,2) * t449 + t510;
t434 = Icges(4,1) * t448 + t509;
t475 = -t433 * t448 + t434 * t449;
t436 = pkin(3) * t448 - pkin(8) * t449;
t474 = -qJD(3) * t436 + qJD(6) * t448;
t473 = t430 * t377 - t378 * t431 + t491;
t472 = qJD(1) * t424 + t483;
t451 = qJD(2) * t462;
t471 = t451 + (-t423 + t492) * qJD(1);
t468 = -t436 * t450 + t472;
t467 = -t436 * t488 + t471;
t444 = rSges(2,1) * t464 - rSges(2,2) * t462;
t443 = rSges(2,1) * t462 + rSges(2,2) * t464;
t435 = rSges(4,1) * t448 + rSges(4,2) * t449;
t432 = Icges(4,5) * t448 + Icges(4,6) * t449;
t429 = qJD(1) + (-qJD(4) - qJD(5)) * t449;
t428 = t449 * t497 + t500;
t427 = -t449 * t499 + t498;
t426 = t449 * t498 - t499;
t425 = -t449 * t500 - t497;
t415 = rSges(4,3) * t462 + t464 * t481;
t414 = -rSges(4,3) * t464 + t462 * t481;
t409 = Icges(4,3) * t462 + t464 * t478;
t408 = -Icges(4,3) * t464 + t462 * t478;
t407 = t462 * t486 + t431;
t406 = t464 * t486 + t430;
t404 = -rSges(5,3) * t449 + (rSges(5,1) * t463 - rSges(5,2) * t461) * t448;
t403 = -Icges(5,5) * t449 + (Icges(5,1) * t463 - Icges(5,4) * t461) * t448;
t402 = -Icges(5,6) * t449 + (Icges(5,4) * t463 - Icges(5,2) * t461) * t448;
t401 = -Icges(5,3) * t449 + (Icges(5,5) * t463 - Icges(5,6) * t461) * t448;
t399 = -rSges(6,3) * t449 + (rSges(6,1) * t453 - rSges(6,2) * t452) * t448;
t390 = qJD(1) * t462 * rSges(3,3) + t438 + (qJD(1) * t482 - qJD(2)) * t464;
t389 = t451 + (t464 * rSges(3,3) - t462 * t482 - t442) * qJD(1);
t387 = rSges(5,1) * t428 + rSges(5,2) * t427 + rSges(5,3) * t505;
t386 = rSges(5,1) * t426 + rSges(5,2) * t425 + rSges(5,3) * t506;
t384 = Icges(5,1) * t428 + Icges(5,4) * t427 + Icges(5,5) * t505;
t383 = Icges(5,1) * t426 + Icges(5,4) * t425 + Icges(5,5) * t506;
t382 = Icges(5,4) * t428 + Icges(5,2) * t427 + Icges(5,6) * t505;
t381 = Icges(5,4) * t426 + Icges(5,2) * t425 + Icges(5,6) * t506;
t380 = Icges(5,5) * t428 + Icges(5,6) * t427 + Icges(5,3) * t505;
t379 = Icges(5,5) * t426 + Icges(5,6) * t425 + Icges(5,3) * t506;
t376 = rSges(6,1) * t421 + rSges(6,2) * t420 + rSges(6,3) * t505;
t374 = rSges(6,1) * t419 + rSges(6,2) * t418 + rSges(6,3) * t506;
t360 = (t414 * t462 + t415 * t464) * qJD(3);
t355 = qJD(1) * t415 - t435 * t450 + t483;
t354 = -t435 * t488 + t451 + (-t414 + t492) * qJD(1);
t353 = t386 * t430 - t387 * t431 + t491;
t352 = t387 * t445 - t404 * t430 + t468;
t351 = -t386 * t445 + t404 * t431 + t467;
t350 = t376 * t429 - t399 * t406 + t468 + t529;
t349 = -t374 * t429 + t399 * t407 + t467 + t522;
t348 = t374 * t406 - t376 * t407 + t473;
t347 = -t406 * t493 + t429 * t494 + t462 * t474 + t472 + t529;
t346 = t407 * t493 - t429 * t495 + t464 * t474 + t471 + t522;
t345 = -qJD(6) * t449 + t406 * t495 - t407 * t494 + t473;
t1 = m(6) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(5) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(4) * (t354 ^ 2 + t355 ^ 2 + t360 ^ 2) / 0.2e1 + m(3) * (t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(7) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + t445 * ((-t379 * t431 - t380 * t430 - t401 * t445) * t449 + ((-t382 * t461 + t384 * t463) * t430 + (-t381 * t461 + t383 * t463) * t431 + (-t402 * t461 + t403 * t463) * t445) * t448) / 0.2e1 + t430 * ((t380 * t505 + t427 * t382 + t428 * t384) * t430 + (t379 * t505 + t381 * t427 + t383 * t428) * t431 + (t401 * t505 + t402 * t427 + t403 * t428) * t445) / 0.2e1 + t431 * ((t380 * t506 + t382 * t425 + t384 * t426) * t430 + (t379 * t506 + t425 * t381 + t426 * t383) * t431 + (t401 * t506 + t402 * t425 + t403 * t426) * t445) / 0.2e1 - ((-t464 * t432 + t462 * t475) * qJD(1) + (t464 ^ 2 * t408 + (t476 * t462 + (-t409 + t477) * t464) * t462) * qJD(3)) * t488 / 0.2e1 + qJD(1) * ((t449 * t433 + t448 * t434) * qJD(1) + ((t411 * t449 + t413 * t448) * t462 - (t410 * t449 + t412 * t448) * t464) * qJD(3)) / 0.2e1 + ((t462 * t432 + t464 * t475) * qJD(1) + (t462 ^ 2 * t409 + (t477 * t464 + (-t408 + t476) * t462) * t464) * qJD(3)) * t450 / 0.2e1 + ((t520 * t420 + t519 * t421 + t521 * t505) * t429 + (t526 * t420 + t524 * t421 + t528 * t505) * t407 + (t525 * t420 + t523 * t421 + t527 * t505) * t406) * t406 / 0.2e1 + ((t520 * t418 + t519 * t419 + t521 * t506) * t429 + (t526 * t418 + t524 * t419 + t528 * t506) * t407 + (t525 * t418 + t523 * t419 + t527 * t506) * t406) * t407 / 0.2e1 + ((-t527 * t406 - t528 * t407 - t521 * t429) * t449 + ((-t520 * t452 + t519 * t453) * t429 + (-t526 * t452 + t524 * t453) * t407 + (-t525 * t452 + t523 * t453) * t406) * t448) * t429 / 0.2e1 + (m(2) * (t443 ^ 2 + t444 ^ 2) + Icges(3,2) * t459 ^ 2 + (Icges(3,1) * t458 + 0.2e1 * Icges(3,4) * t459) * t458 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
