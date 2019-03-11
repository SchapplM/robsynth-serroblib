% Calculate kinetic energy for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:26
% EndTime: 2019-03-09 20:58:29
% DurationCPUTime: 2.88s
% Computational Cost: add. (2017->302), mult. (2649->459), div. (0->0), fcn. (2672->10), ass. (0->150)
t536 = Icges(6,1) + Icges(7,1);
t535 = -Icges(6,4) + Icges(7,5);
t534 = Icges(7,4) + Icges(6,5);
t533 = Icges(6,2) + Icges(7,3);
t532 = -Icges(7,6) + Icges(6,6);
t531 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t530 = rSges(7,1) + pkin(5);
t529 = rSges(7,3) + qJ(6);
t467 = qJ(3) + qJ(4);
t461 = pkin(10) + t467;
t458 = sin(t461);
t459 = cos(t461);
t473 = cos(qJ(1));
t470 = sin(qJ(1));
t472 = cos(qJ(2));
t506 = t472 * t470;
t406 = t458 * t506 + t459 * t473;
t407 = -t458 * t473 + t459 * t506;
t469 = sin(qJ(2));
t508 = t469 * t470;
t528 = t406 * t533 + t407 * t535 - t508 * t532;
t505 = t472 * t473;
t408 = t458 * t505 - t459 * t470;
t409 = t458 * t470 + t459 * t505;
t507 = t469 * t473;
t527 = t408 * t533 + t409 * t535 - t507 * t532;
t526 = t406 * t535 + t407 * t536 + t508 * t534;
t525 = t408 * t535 + t409 * t536 + t507 * t534;
t524 = t532 * t472 + (t458 * t533 + t459 * t535) * t469;
t523 = -t534 * t472 + (t458 * t535 + t459 * t536) * t469;
t463 = sin(t467);
t464 = cos(t467);
t425 = -t463 * t506 - t464 * t473;
t426 = -t463 * t473 + t464 * t506;
t522 = Icges(5,5) * t426 + Icges(5,6) * t425 - t406 * t532 + t407 * t534 - t508 * t531;
t427 = -t463 * t505 + t464 * t470;
t428 = t463 * t470 + t464 * t505;
t521 = Icges(5,5) * t428 + Icges(5,6) * t427 - t408 * t532 + t409 * t534 - t507 * t531;
t520 = t531 * t472 + (Icges(5,5) * t464 - Icges(5,6) * t463 - t458 * t532 + t459 * t534) * t469;
t471 = cos(qJ(3));
t514 = t471 * pkin(3);
t512 = Icges(3,4) * t469;
t511 = Icges(3,4) * t472;
t468 = sin(qJ(3));
t510 = t468 * t470;
t509 = t468 * t473;
t504 = rSges(7,2) * t508 + t406 * t529 + t407 * t530;
t503 = rSges(7,2) * t507 + t408 * t529 + t409 * t530;
t502 = -rSges(7,2) * t472 + (t458 * t529 + t459 * t530) * t469;
t493 = pkin(2) * t472 + pkin(8) * t469;
t436 = t493 * t470;
t437 = t493 * t473;
t462 = qJD(2) * t470;
t498 = qJD(2) * t473;
t501 = t436 * t462 + t437 * t498;
t500 = pkin(4) * t464;
t497 = qJD(3) * t469;
t438 = t473 * t497 + t462;
t496 = qJD(4) * t469;
t495 = qJD(5) * t469;
t494 = pkin(4) * t463;
t439 = t470 * t497 - t498;
t492 = rSges(3,1) * t472 - rSges(3,2) * t469;
t491 = Icges(3,1) * t472 - t512;
t490 = -Icges(3,2) * t469 + t511;
t489 = Icges(3,5) * t472 - Icges(3,6) * t469;
t417 = -Icges(3,6) * t473 + t470 * t490;
t420 = -Icges(3,5) * t473 + t470 * t491;
t488 = t417 * t469 - t420 * t472;
t418 = Icges(3,6) * t470 + t473 * t490;
t421 = Icges(3,5) * t470 + t473 * t491;
t487 = -t418 * t469 + t421 * t472;
t445 = Icges(3,2) * t472 + t512;
t446 = Icges(3,1) * t469 + t511;
t486 = -t445 * t469 + t446 * t472;
t483 = pkin(9) * t469 + t472 * t514;
t389 = -pkin(3) * t509 + t470 * t483;
t390 = pkin(3) * t510 + t473 * t483;
t485 = t389 * t438 - t390 * t439 + t501;
t442 = qJD(1) * (pkin(1) * t473 + pkin(7) * t470);
t451 = pkin(2) * t469 - pkin(8) * t472;
t484 = qJD(1) * t437 - t451 * t462 + t442;
t452 = pkin(1) * t470 - pkin(7) * t473;
t482 = (-t436 - t452) * qJD(1) - t451 * t498;
t481 = qJ(5) * t469 + t472 * t500;
t401 = -pkin(9) * t472 + t469 * t514;
t457 = -qJD(3) * t472 + qJD(1);
t480 = t390 * t457 - t401 * t438 + t484;
t347 = t470 * t481 - t473 * t494;
t411 = t473 * t496 + t438;
t479 = -qJD(5) * t472 + t347 * t411 + t485;
t348 = t470 * t494 + t473 * t481;
t441 = qJD(1) + (-qJD(3) - qJD(4)) * t472;
t478 = t348 * t441 + t470 * t495 + t480;
t477 = -t389 * t457 + t401 * t439 + t482;
t391 = -qJ(5) * t472 + t469 * t500;
t412 = t470 * t496 + t439;
t476 = t391 * t412 + t473 * t495 + t477;
t449 = rSges(2,1) * t473 - rSges(2,2) * t470;
t448 = rSges(2,1) * t470 + rSges(2,2) * t473;
t447 = rSges(3,1) * t469 + rSges(3,2) * t472;
t444 = Icges(3,5) * t469 + Icges(3,6) * t472;
t435 = t471 * t505 + t510;
t434 = -t468 * t505 + t470 * t471;
t433 = t471 * t506 - t509;
t432 = -t468 * t506 - t471 * t473;
t424 = rSges(3,3) * t470 + t473 * t492;
t423 = -rSges(3,3) * t473 + t470 * t492;
t422 = -rSges(4,3) * t472 + (rSges(4,1) * t471 - rSges(4,2) * t468) * t469;
t419 = -Icges(4,5) * t472 + (Icges(4,1) * t471 - Icges(4,4) * t468) * t469;
t416 = -Icges(4,6) * t472 + (Icges(4,4) * t471 - Icges(4,2) * t468) * t469;
t415 = Icges(3,3) * t470 + t473 * t489;
t414 = -Icges(3,3) * t473 + t470 * t489;
t413 = -Icges(4,3) * t472 + (Icges(4,5) * t471 - Icges(4,6) * t468) * t469;
t405 = -rSges(5,3) * t472 + (rSges(5,1) * t464 - rSges(5,2) * t463) * t469;
t404 = -Icges(5,5) * t472 + (Icges(5,1) * t464 - Icges(5,4) * t463) * t469;
t403 = -Icges(5,6) * t472 + (Icges(5,4) * t464 - Icges(5,2) * t463) * t469;
t400 = -rSges(6,3) * t472 + (rSges(6,1) * t459 - rSges(6,2) * t458) * t469;
t388 = rSges(4,1) * t435 + rSges(4,2) * t434 + rSges(4,3) * t507;
t387 = rSges(4,1) * t433 + rSges(4,2) * t432 + rSges(4,3) * t508;
t386 = Icges(4,1) * t435 + Icges(4,4) * t434 + Icges(4,5) * t507;
t385 = Icges(4,1) * t433 + Icges(4,4) * t432 + Icges(4,5) * t508;
t384 = Icges(4,4) * t435 + Icges(4,2) * t434 + Icges(4,6) * t507;
t383 = Icges(4,4) * t433 + Icges(4,2) * t432 + Icges(4,6) * t508;
t382 = Icges(4,5) * t435 + Icges(4,6) * t434 + Icges(4,3) * t507;
t381 = Icges(4,5) * t433 + Icges(4,6) * t432 + Icges(4,3) * t508;
t380 = qJD(1) * t424 - t447 * t462 + t442;
t379 = -t447 * t498 + (-t423 - t452) * qJD(1);
t378 = (t423 * t470 + t424 * t473) * qJD(2);
t374 = rSges(5,1) * t428 + rSges(5,2) * t427 + rSges(5,3) * t507;
t373 = rSges(5,1) * t426 + rSges(5,2) * t425 + rSges(5,3) * t508;
t372 = Icges(5,1) * t428 + Icges(5,4) * t427 + Icges(5,5) * t507;
t371 = Icges(5,1) * t426 + Icges(5,4) * t425 + Icges(5,5) * t508;
t370 = Icges(5,4) * t428 + Icges(5,2) * t427 + Icges(5,6) * t507;
t369 = Icges(5,4) * t426 + Icges(5,2) * t425 + Icges(5,6) * t508;
t365 = rSges(6,1) * t409 - rSges(6,2) * t408 + rSges(6,3) * t507;
t363 = rSges(6,1) * t407 - rSges(6,2) * t406 + rSges(6,3) * t508;
t344 = t388 * t457 - t422 * t438 + t484;
t343 = -t387 * t457 + t422 * t439 + t482;
t342 = t387 * t438 - t388 * t439 + t501;
t341 = t374 * t441 - t405 * t411 + t480;
t340 = -t373 * t441 + t405 * t412 + t477;
t339 = t373 * t411 - t374 * t412 + t485;
t338 = t365 * t441 + (-t391 - t400) * t411 + t478;
t337 = t400 * t412 + (-t347 - t363) * t441 + t476;
t336 = t363 * t411 + (-t348 - t365) * t412 + t479;
t335 = qJD(6) * t406 + t503 * t441 + (-t391 - t502) * t411 + t478;
t334 = qJD(6) * t408 + t502 * t412 + (-t347 - t504) * t441 + t476;
t333 = qJD(6) * t458 * t469 + t504 * t411 + (-t348 - t503) * t412 + t479;
t1 = m(6) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(5) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(4) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(3) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + t438 * ((t382 * t507 + t434 * t384 + t435 * t386) * t438 + (t381 * t507 + t383 * t434 + t385 * t435) * t439 + (t413 * t507 + t416 * t434 + t419 * t435) * t457) / 0.2e1 + t439 * ((t382 * t508 + t384 * t432 + t386 * t433) * t438 + (t381 * t508 + t432 * t383 + t433 * t385) * t439 + (t413 * t508 + t416 * t432 + t419 * t433) * t457) / 0.2e1 + t457 * ((-t381 * t439 - t382 * t438 - t413 * t457) * t472 + ((-t384 * t468 + t386 * t471) * t438 + (-t383 * t468 + t385 * t471) * t439 + (-t416 * t468 + t419 * t471) * t457) * t469) / 0.2e1 + qJD(1) * ((t472 * t445 + t469 * t446) * qJD(1) + ((t418 * t472 + t421 * t469) * t470 - (t417 * t472 + t420 * t469) * t473) * qJD(2)) / 0.2e1 - ((-t444 * t473 + t470 * t486) * qJD(1) + (t473 ^ 2 * t414 + (t487 * t470 + (-t415 + t488) * t473) * t470) * qJD(2)) * t498 / 0.2e1 + ((t470 * t444 + t473 * t486) * qJD(1) + (t470 ^ 2 * t415 + (t488 * t473 + (-t414 + t487) * t470) * t473) * qJD(2)) * t462 / 0.2e1 + (m(2) * (t448 ^ 2 + t449 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t403 * t427 + t404 * t428 + t408 * t524 + t409 * t523 + t507 * t520) * t441 + (t369 * t427 + t371 * t428 + t408 * t528 + t409 * t526 + t507 * t522) * t412 + (t427 * t370 + t428 * t372 + t527 * t408 + t525 * t409 + t521 * t507) * t411) * t411 / 0.2e1 + ((t403 * t425 + t404 * t426 + t406 * t524 + t407 * t523 + t508 * t520) * t441 + (t425 * t369 + t426 * t371 + t528 * t406 + t526 * t407 + t522 * t508) * t412 + (t370 * t425 + t372 * t426 + t527 * t406 + t525 * t407 + t521 * t508) * t411) * t412 / 0.2e1 + ((-t411 * t521 - t412 * t522 - t441 * t520) * t472 + ((-t403 * t463 + t404 * t464 + t458 * t524 + t459 * t523) * t441 + (-t369 * t463 + t371 * t464 + t458 * t528 + t459 * t526) * t412 + (-t370 * t463 + t372 * t464 + t458 * t527 + t459 * t525) * t411) * t469) * t441 / 0.2e1;
T  = t1;
