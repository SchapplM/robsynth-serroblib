% Calculate kinetic energy for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:31
% EndTime: 2019-03-09 10:49:33
% DurationCPUTime: 3.00s
% Computational Cost: add. (1820->319), mult. (2797->485), div. (0->0), fcn. (2974->10), ass. (0->153)
t524 = Icges(5,1) + Icges(6,1);
t523 = -Icges(5,4) + Icges(6,5);
t522 = Icges(6,4) + Icges(5,5);
t521 = Icges(5,2) + Icges(6,3);
t520 = -Icges(6,6) + Icges(5,6);
t519 = -Icges(5,3) - Icges(6,2);
t458 = pkin(10) + qJ(4);
t455 = sin(t458);
t456 = cos(t458);
t467 = cos(qJ(1));
t464 = sin(qJ(1));
t466 = cos(qJ(2));
t497 = t466 * t464;
t412 = t455 * t497 + t456 * t467;
t413 = -t455 * t467 + t456 * t497;
t463 = sin(qJ(2));
t499 = t463 * t464;
t518 = t521 * t412 + t523 * t413 - t520 * t499;
t496 = t466 * t467;
t414 = t455 * t496 - t464 * t456;
t415 = t455 * t464 + t456 * t496;
t498 = t463 * t467;
t517 = t521 * t414 + t523 * t415 - t520 * t498;
t516 = -t520 * t412 + t522 * t413 - t519 * t499;
t515 = -t520 * t414 + t522 * t415 - t519 * t498;
t514 = t523 * t412 + t524 * t413 + t522 * t499;
t513 = t523 * t414 + t524 * t415 + t522 * t498;
t512 = t520 * t466 + (t521 * t455 + t523 * t456) * t463;
t511 = t519 * t466 + (-t520 * t455 + t522 * t456) * t463;
t510 = -t522 * t466 + (t523 * t455 + t524 * t456) * t463;
t460 = cos(pkin(10));
t504 = pkin(3) * t460;
t503 = Icges(3,4) * t463;
t502 = Icges(3,4) * t466;
t459 = sin(pkin(10));
t501 = t459 * t464;
t500 = t459 * t467;
t482 = pkin(2) * t466 + qJ(3) * t463;
t432 = t482 * t464;
t447 = pkin(1) * t464 - pkin(7) * t467;
t494 = -t432 - t447;
t457 = qJD(2) * t464;
t491 = qJD(4) * t463;
t435 = t467 * t491 + t457;
t493 = qJD(2) * t467;
t492 = qJD(3) * t463;
t490 = qJD(6) * t463;
t433 = t482 * t467;
t438 = qJD(1) * (pkin(1) * t467 + pkin(7) * t464);
t489 = qJD(1) * t433 + t464 * t492 + t438;
t443 = pkin(2) * t463 - qJ(3) * t466;
t486 = qJD(2) * (pkin(8) * t466 - t463 * t504 - t443);
t485 = qJD(2) * (rSges(4,3) * t466 - (rSges(4,1) * t460 - rSges(4,2) * t459) * t463 - t443);
t436 = t464 * t491 - t493;
t484 = -qJD(3) * t466 + t432 * t457 + t433 * t493;
t483 = rSges(3,1) * t466 - rSges(3,2) * t463;
t481 = Icges(3,1) * t466 - t503;
t480 = -Icges(3,2) * t463 + t502;
t479 = Icges(3,5) * t466 - Icges(3,6) * t463;
t418 = -Icges(3,6) * t467 + t464 * t480;
t420 = -Icges(3,5) * t467 + t464 * t481;
t478 = t418 * t463 - t420 * t466;
t419 = Icges(3,6) * t464 + t467 * t480;
t421 = Icges(3,5) * t464 + t467 * t481;
t477 = -t419 * t463 + t421 * t466;
t440 = Icges(3,2) * t466 + t503;
t441 = Icges(3,1) * t463 + t502;
t476 = -t440 * t463 + t441 * t466;
t474 = pkin(8) * t463 + t466 * t504;
t388 = -pkin(3) * t500 + t464 * t474;
t389 = pkin(3) * t501 + t467 * t474;
t475 = t388 * t457 + t389 * t493 + t484;
t378 = pkin(4) * t413 + qJ(5) * t412;
t473 = qJD(5) * t463 * t455 + t435 * t378 + t475;
t472 = qJD(1) * t389 + t464 * t486 + t489;
t379 = pkin(4) * t415 + qJ(5) * t414;
t453 = -qJD(4) * t466 + qJD(1);
t471 = qJD(5) * t412 + t453 * t379 + t472;
t452 = t467 * t492;
t470 = t452 + (-t388 + t494) * qJD(1) + t467 * t486;
t424 = (pkin(4) * t456 + qJ(5) * t455) * t463;
t469 = qJD(5) * t414 + t436 * t424 + t470;
t465 = cos(qJ(6));
t462 = sin(qJ(6));
t446 = rSges(2,1) * t467 - rSges(2,2) * t464;
t445 = rSges(2,1) * t464 + rSges(2,2) * t467;
t444 = rSges(3,1) * t463 + rSges(3,2) * t466;
t439 = Icges(3,5) * t463 + Icges(3,6) * t466;
t437 = qJD(1) + (-qJD(4) + qJD(6)) * t466;
t434 = pkin(5) * t456 * t463 + pkin(9) * t466;
t431 = t460 * t496 + t501;
t430 = -t459 * t496 + t460 * t464;
t429 = t460 * t497 - t500;
t428 = -t459 * t497 - t460 * t467;
t426 = rSges(3,3) * t464 + t467 * t483;
t425 = -rSges(3,3) * t467 + t464 * t483;
t417 = Icges(3,3) * t464 + t467 * t479;
t416 = -Icges(3,3) * t467 + t464 * t479;
t410 = -t464 * t490 + t436;
t409 = -t467 * t490 + t435;
t408 = -Icges(4,5) * t466 + (Icges(4,1) * t460 - Icges(4,4) * t459) * t463;
t407 = -Icges(4,6) * t466 + (Icges(4,4) * t460 - Icges(4,2) * t459) * t463;
t406 = -Icges(4,3) * t466 + (Icges(4,5) * t460 - Icges(4,6) * t459) * t463;
t403 = (t455 * t462 + t456 * t465) * t463;
t402 = (t455 * t465 - t456 * t462) * t463;
t401 = -rSges(5,3) * t466 + (rSges(5,1) * t456 - rSges(5,2) * t455) * t463;
t400 = -rSges(6,2) * t466 + (rSges(6,1) * t456 + rSges(6,3) * t455) * t463;
t392 = pkin(5) * t415 - pkin(9) * t498;
t391 = pkin(5) * t413 - pkin(9) * t499;
t387 = rSges(4,1) * t431 + rSges(4,2) * t430 + rSges(4,3) * t498;
t386 = rSges(4,1) * t429 + rSges(4,2) * t428 + rSges(4,3) * t499;
t385 = Icges(4,1) * t431 + Icges(4,4) * t430 + Icges(4,5) * t498;
t384 = Icges(4,1) * t429 + Icges(4,4) * t428 + Icges(4,5) * t499;
t383 = Icges(4,4) * t431 + Icges(4,2) * t430 + Icges(4,6) * t498;
t382 = Icges(4,4) * t429 + Icges(4,2) * t428 + Icges(4,6) * t499;
t381 = Icges(4,5) * t431 + Icges(4,6) * t430 + Icges(4,3) * t498;
t380 = Icges(4,5) * t429 + Icges(4,6) * t428 + Icges(4,3) * t499;
t376 = t414 * t462 + t415 * t465;
t375 = t414 * t465 - t415 * t462;
t374 = t412 * t462 + t413 * t465;
t373 = t412 * t465 - t413 * t462;
t372 = qJD(1) * t426 - t444 * t457 + t438;
t371 = -t444 * t493 + (-t425 - t447) * qJD(1);
t368 = (t425 * t464 + t426 * t467) * qJD(2);
t367 = rSges(5,1) * t415 - rSges(5,2) * t414 + rSges(5,3) * t498;
t366 = rSges(6,1) * t415 + rSges(6,2) * t498 + rSges(6,3) * t414;
t365 = rSges(5,1) * t413 - rSges(5,2) * t412 + rSges(5,3) * t499;
t364 = rSges(6,1) * t413 + rSges(6,2) * t499 + rSges(6,3) * t412;
t350 = rSges(7,1) * t403 + rSges(7,2) * t402 + rSges(7,3) * t466;
t349 = Icges(7,1) * t403 + Icges(7,4) * t402 + Icges(7,5) * t466;
t348 = Icges(7,4) * t403 + Icges(7,2) * t402 + Icges(7,6) * t466;
t347 = Icges(7,5) * t403 + Icges(7,6) * t402 + Icges(7,3) * t466;
t345 = rSges(7,1) * t376 + rSges(7,2) * t375 - rSges(7,3) * t498;
t344 = rSges(7,1) * t374 + rSges(7,2) * t373 - rSges(7,3) * t499;
t343 = Icges(7,1) * t376 + Icges(7,4) * t375 - Icges(7,5) * t498;
t342 = Icges(7,1) * t374 + Icges(7,4) * t373 - Icges(7,5) * t499;
t341 = Icges(7,4) * t376 + Icges(7,2) * t375 - Icges(7,6) * t498;
t340 = Icges(7,4) * t374 + Icges(7,2) * t373 - Icges(7,6) * t499;
t339 = Icges(7,5) * t376 + Icges(7,6) * t375 - Icges(7,3) * t498;
t338 = Icges(7,5) * t374 + Icges(7,6) * t373 - Icges(7,3) * t499;
t337 = qJD(1) * t387 + t464 * t485 + t489;
t336 = t452 + t467 * t485 + (-t386 + t494) * qJD(1);
t335 = (t386 * t464 + t387 * t467) * qJD(2) + t484;
t334 = t367 * t453 - t401 * t435 + t472;
t333 = -t365 * t453 + t401 * t436 + t470;
t332 = t365 * t435 - t367 * t436 + t475;
t331 = t366 * t453 + (-t400 - t424) * t435 + t471;
t330 = t400 * t436 + (-t364 - t378) * t453 + t469;
t329 = t364 * t435 + (-t366 - t379) * t436 + t473;
t328 = t345 * t437 - t350 * t409 + t392 * t453 + (-t424 - t434) * t435 + t471;
t327 = -t344 * t437 + t350 * t410 + t434 * t436 + (-t378 - t391) * t453 + t469;
t326 = t344 * t409 - t345 * t410 + t391 * t435 + (-t379 - t392) * t436 + t473;
t1 = m(7) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + m(5) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(6) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(4) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(3) * (t368 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + t437 * ((t339 * t466 + t341 * t402 + t343 * t403) * t409 + (t338 * t466 + t340 * t402 + t342 * t403) * t410 + (t466 * t347 + t402 * t348 + t403 * t349) * t437) / 0.2e1 + t409 * ((-t339 * t498 + t375 * t341 + t376 * t343) * t409 + (-t338 * t498 + t340 * t375 + t342 * t376) * t410 + (-t347 * t498 + t348 * t375 + t349 * t376) * t437) / 0.2e1 + t410 * ((-t339 * t499 + t341 * t373 + t343 * t374) * t409 + (-t338 * t499 + t373 * t340 + t374 * t342) * t410 + (-t347 * t499 + t348 * t373 + t349 * t374) * t437) / 0.2e1 + ((t512 * t414 + t510 * t415 + t511 * t498) * t453 + (t518 * t414 + t514 * t415 + t516 * t498) * t436 + (t517 * t414 + t513 * t415 + t515 * t498) * t435) * t435 / 0.2e1 + ((t512 * t412 + t510 * t413 + t511 * t499) * t453 + (t518 * t412 + t514 * t413 + t516 * t499) * t436 + (t517 * t412 + t513 * t413 + t515 * t499) * t435) * t436 / 0.2e1 + ((-t515 * t435 - t516 * t436 - t511 * t453) * t466 + ((t512 * t455 + t510 * t456) * t453 + (t518 * t455 + t514 * t456) * t436 + (t517 * t455 + t513 * t456) * t435) * t463) * t453 / 0.2e1 + (m(2) * (t445 ^ 2 + t446 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t419 * t466 + t421 * t463) * t464 - (t418 * t466 + t420 * t463) * t467 + (t380 * t467 - t381 * t464) * t466 + ((-t383 * t459 + t385 * t460) * t464 - (-t382 * t459 + t384 * t460) * t467) * t463) * qJD(2) + ((t440 - t406) * t466 + (-t407 * t459 + t408 * t460 + t441) * t463) * qJD(1)) * qJD(1) / 0.2e1 + (((-t380 * t498 - t382 * t430 - t384 * t431 + t478 * t467) * t467 + (t381 * t498 + t383 * t430 + t385 * t431 + (-t416 + t477) * t467 + t417 * t464) * t464) * qJD(2) + (t406 * t498 + t407 * t430 + t408 * t431 + t464 * t439 + t467 * t476) * qJD(1)) * t457 / 0.2e1 - (((-t380 * t499 - t382 * t428 - t384 * t429 + t416 * t467) * t467 + ((-t417 + t478) * t467 + t381 * t499 + t383 * t428 + t385 * t429 + t477 * t464) * t464) * qJD(2) + (t406 * t499 + t407 * t428 + t408 * t429 - t467 * t439 + t464 * t476) * qJD(1)) * t493 / 0.2e1;
T  = t1;
