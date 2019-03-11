% Calculate kinetic energy for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:13
% EndTime: 2019-03-09 07:12:15
% DurationCPUTime: 2.10s
% Computational Cost: add. (2083->303), mult. (2012->494), div. (0->0), fcn. (2002->12), ass. (0->156)
t463 = cos(pkin(11));
t515 = pkin(2) * t463;
t467 = cos(qJ(4));
t514 = t467 * pkin(4);
t459 = pkin(11) + qJ(3);
t451 = sin(t459);
t512 = Icges(4,4) * t451;
t452 = cos(t459);
t511 = Icges(4,4) * t452;
t466 = sin(qJ(1));
t510 = t451 * t466;
t468 = cos(qJ(1));
t509 = t451 * t468;
t508 = t452 * t466;
t507 = t452 * t468;
t461 = qJ(4) + qJ(5);
t455 = sin(t461);
t506 = t455 * t466;
t505 = t455 * t468;
t456 = cos(t461);
t504 = t456 * t466;
t503 = t456 * t468;
t465 = sin(qJ(4));
t502 = t465 * t466;
t501 = t465 * t468;
t500 = t466 * t467;
t499 = t467 * t468;
t443 = pkin(1) * t466 - qJ(2) * t468;
t497 = pkin(7) * t468 - t466 * t515 - t443;
t487 = pkin(3) * t452 + pkin(8) * t451;
t422 = t487 * t466;
t423 = t487 * t468;
t453 = qJD(3) * t466;
t493 = qJD(3) * t468;
t496 = t422 * t453 + t423 * t493;
t495 = pkin(5) * t456;
t492 = qJD(4) * t451;
t429 = t468 * t492 + t453;
t491 = qJD(5) * t451;
t490 = qJD(6) * t451;
t489 = -qJD(4) - qJD(5);
t400 = t468 * t491 + t429;
t488 = pkin(5) * t455;
t430 = t466 * t492 - t493;
t437 = qJD(1) * (pkin(1) * t468 + qJ(2) * t466);
t486 = -qJD(2) * t468 + qJD(1) * (pkin(7) * t466 + t468 * t515) + t437;
t401 = t466 * t491 + t430;
t462 = sin(pkin(11));
t485 = rSges(3,1) * t463 - rSges(3,2) * t462;
t484 = rSges(4,1) * t452 - rSges(4,2) * t451;
t483 = Icges(4,1) * t452 - t512;
t482 = -Icges(4,2) * t451 + t511;
t481 = Icges(4,5) * t452 - Icges(4,6) * t451;
t404 = -Icges(4,6) * t468 + t466 * t482;
t406 = -Icges(4,5) * t468 + t466 * t483;
t480 = t404 * t451 - t406 * t452;
t405 = Icges(4,6) * t466 + t468 * t482;
t407 = Icges(4,5) * t466 + t468 * t483;
t479 = -t405 * t451 + t407 * t452;
t432 = Icges(4,2) * t452 + t512;
t433 = Icges(4,1) * t451 + t511;
t478 = -t432 * t451 + t433 * t452;
t476 = pkin(9) * t451 + t452 * t514;
t369 = -pkin(4) * t501 + t466 * t476;
t370 = pkin(4) * t502 + t468 * t476;
t477 = t429 * t369 - t370 * t430 + t496;
t475 = pkin(10) * t451 + t452 * t495;
t435 = pkin(3) * t451 - pkin(8) * t452;
t474 = qJD(1) * t423 - t435 * t453 + t486;
t454 = qJD(2) * t466;
t473 = t454 + (-t422 + t497) * qJD(1) - t435 * t493;
t385 = -pkin(9) * t452 + t451 * t514;
t446 = -qJD(4) * t452 + qJD(1);
t472 = t446 * t370 - t385 * t429 + t474;
t471 = -t369 * t446 + t430 * t385 + t473;
t457 = qJ(6) + t461;
t449 = cos(t457);
t448 = sin(t457);
t445 = rSges(2,1) * t468 - rSges(2,2) * t466;
t444 = rSges(2,1) * t466 + rSges(2,2) * t468;
t434 = rSges(4,1) * t451 + rSges(4,2) * t452;
t431 = Icges(4,5) * t451 + Icges(4,6) * t452;
t428 = t452 * t489 + qJD(1);
t427 = t452 * t499 + t502;
t426 = -t452 * t501 + t500;
t425 = t452 * t500 - t501;
t424 = -t452 * t502 - t499;
t420 = t452 * t503 + t506;
t419 = -t452 * t505 + t504;
t418 = t452 * t504 - t505;
t417 = -t452 * t506 - t503;
t416 = qJD(1) + (-qJD(6) + t489) * t452;
t413 = rSges(4,3) * t466 + t468 * t484;
t412 = -rSges(4,3) * t468 + t466 * t484;
t411 = t448 * t466 + t449 * t507;
t410 = -t448 * t507 + t449 * t466;
t409 = -t448 * t468 + t449 * t508;
t408 = -t448 * t508 - t449 * t468;
t403 = Icges(4,3) * t466 + t468 * t481;
t402 = -Icges(4,3) * t468 + t466 * t481;
t398 = -rSges(5,3) * t452 + (rSges(5,1) * t467 - rSges(5,2) * t465) * t451;
t397 = -Icges(5,5) * t452 + (Icges(5,1) * t467 - Icges(5,4) * t465) * t451;
t396 = -Icges(5,6) * t452 + (Icges(5,4) * t467 - Icges(5,2) * t465) * t451;
t395 = -Icges(5,3) * t452 + (Icges(5,5) * t467 - Icges(5,6) * t465) * t451;
t393 = -rSges(6,3) * t452 + (rSges(6,1) * t456 - rSges(6,2) * t455) * t451;
t392 = -Icges(6,5) * t452 + (Icges(6,1) * t456 - Icges(6,4) * t455) * t451;
t391 = -Icges(6,6) * t452 + (Icges(6,4) * t456 - Icges(6,2) * t455) * t451;
t390 = -Icges(6,3) * t452 + (Icges(6,5) * t456 - Icges(6,6) * t455) * t451;
t389 = -rSges(7,3) * t452 + (rSges(7,1) * t449 - rSges(7,2) * t448) * t451;
t388 = -Icges(7,5) * t452 + (Icges(7,1) * t449 - Icges(7,4) * t448) * t451;
t387 = -Icges(7,6) * t452 + (Icges(7,4) * t449 - Icges(7,2) * t448) * t451;
t386 = -Icges(7,3) * t452 + (Icges(7,5) * t449 - Icges(7,6) * t448) * t451;
t384 = t466 * t490 + t401;
t383 = t468 * t490 + t400;
t382 = qJD(1) * t466 * rSges(3,3) + t437 + (qJD(1) * t485 - qJD(2)) * t468;
t381 = t454 + (t468 * rSges(3,3) - t466 * t485 - t443) * qJD(1);
t380 = -pkin(10) * t452 + t451 * t495;
t379 = rSges(5,1) * t427 + rSges(5,2) * t426 + rSges(5,3) * t509;
t378 = rSges(5,1) * t425 + rSges(5,2) * t424 + rSges(5,3) * t510;
t376 = Icges(5,1) * t427 + Icges(5,4) * t426 + Icges(5,5) * t509;
t375 = Icges(5,1) * t425 + Icges(5,4) * t424 + Icges(5,5) * t510;
t374 = Icges(5,4) * t427 + Icges(5,2) * t426 + Icges(5,6) * t509;
t373 = Icges(5,4) * t425 + Icges(5,2) * t424 + Icges(5,6) * t510;
t372 = Icges(5,5) * t427 + Icges(5,6) * t426 + Icges(5,3) * t509;
t371 = Icges(5,5) * t425 + Icges(5,6) * t424 + Icges(5,3) * t510;
t368 = rSges(6,1) * t420 + rSges(6,2) * t419 + rSges(6,3) * t509;
t367 = rSges(6,1) * t418 + rSges(6,2) * t417 + rSges(6,3) * t510;
t366 = Icges(6,1) * t420 + Icges(6,4) * t419 + Icges(6,5) * t509;
t365 = Icges(6,1) * t418 + Icges(6,4) * t417 + Icges(6,5) * t510;
t364 = Icges(6,4) * t420 + Icges(6,2) * t419 + Icges(6,6) * t509;
t363 = Icges(6,4) * t418 + Icges(6,2) * t417 + Icges(6,6) * t510;
t362 = Icges(6,5) * t420 + Icges(6,6) * t419 + Icges(6,3) * t509;
t361 = Icges(6,5) * t418 + Icges(6,6) * t417 + Icges(6,3) * t510;
t360 = (t412 * t466 + t413 * t468) * qJD(3);
t358 = rSges(7,1) * t411 + rSges(7,2) * t410 + rSges(7,3) * t509;
t357 = rSges(7,1) * t409 + rSges(7,2) * t408 + rSges(7,3) * t510;
t356 = Icges(7,1) * t411 + Icges(7,4) * t410 + Icges(7,5) * t509;
t355 = Icges(7,1) * t409 + Icges(7,4) * t408 + Icges(7,5) * t510;
t354 = Icges(7,4) * t411 + Icges(7,2) * t410 + Icges(7,6) * t509;
t353 = Icges(7,4) * t409 + Icges(7,2) * t408 + Icges(7,6) * t510;
t352 = Icges(7,5) * t411 + Icges(7,6) * t410 + Icges(7,3) * t509;
t351 = Icges(7,5) * t409 + Icges(7,6) * t408 + Icges(7,3) * t510;
t349 = t466 * t488 + t468 * t475;
t348 = t466 * t475 - t468 * t488;
t347 = qJD(1) * t413 - t434 * t453 + t486;
t346 = -t434 * t493 + t454 + (-t412 + t497) * qJD(1);
t345 = t378 * t429 - t379 * t430 + t496;
t344 = t379 * t446 - t398 * t429 + t474;
t343 = -t378 * t446 + t398 * t430 + t473;
t342 = t368 * t428 - t393 * t400 + t472;
t341 = -t367 * t428 + t393 * t401 + t471;
t340 = t367 * t400 - t368 * t401 + t477;
t339 = t349 * t428 + t358 * t416 - t380 * t400 - t383 * t389 + t472;
t338 = -t348 * t428 - t357 * t416 + t380 * t401 + t384 * t389 + t471;
t337 = t348 * t400 - t349 * t401 + t357 * t383 - t358 * t384 + t477;
t1 = ((t466 * t431 + t468 * t478) * qJD(1) + (t466 ^ 2 * t403 + (t480 * t468 + (-t402 + t479) * t466) * t468) * qJD(3)) * t453 / 0.2e1 + t400 * ((t362 * t509 + t419 * t364 + t420 * t366) * t400 + (t361 * t509 + t363 * t419 + t365 * t420) * t401 + (t390 * t509 + t391 * t419 + t392 * t420) * t428) / 0.2e1 + t383 * ((t352 * t509 + t410 * t354 + t411 * t356) * t383 + (t351 * t509 + t353 * t410 + t355 * t411) * t384 + (t386 * t509 + t387 * t410 + t388 * t411) * t416) / 0.2e1 + t384 * ((t352 * t510 + t354 * t408 + t356 * t409) * t383 + (t351 * t510 + t408 * t353 + t409 * t355) * t384 + (t386 * t510 + t387 * t408 + t388 * t409) * t416) / 0.2e1 + t416 * ((-t351 * t384 - t352 * t383 - t386 * t416) * t452 + ((-t354 * t448 + t356 * t449) * t383 + (-t353 * t448 + t355 * t449) * t384 + (-t387 * t448 + t388 * t449) * t416) * t451) / 0.2e1 + t401 * ((t362 * t510 + t364 * t417 + t366 * t418) * t400 + (t361 * t510 + t417 * t363 + t418 * t365) * t401 + (t390 * t510 + t391 * t417 + t392 * t418) * t428) / 0.2e1 + t428 * ((-t361 * t401 - t362 * t400 - t390 * t428) * t452 + ((-t364 * t455 + t366 * t456) * t400 + (-t363 * t455 + t365 * t456) * t401 + (-t391 * t455 + t392 * t456) * t428) * t451) / 0.2e1 + t429 * ((t372 * t509 + t426 * t374 + t427 * t376) * t429 + (t371 * t509 + t373 * t426 + t375 * t427) * t430 + (t395 * t509 + t396 * t426 + t397 * t427) * t446) / 0.2e1 + t430 * ((t372 * t510 + t374 * t424 + t376 * t425) * t429 + (t371 * t510 + t424 * t373 + t425 * t375) * t430 + (t395 * t510 + t396 * t424 + t397 * t425) * t446) / 0.2e1 + t446 * ((-t371 * t430 - t372 * t429 - t395 * t446) * t452 + ((-t374 * t465 + t376 * t467) * t429 + (-t373 * t465 + t375 * t467) * t430 + (-t396 * t465 + t397 * t467) * t446) * t451) / 0.2e1 - ((-t468 * t431 + t478 * t466) * qJD(1) + (t468 ^ 2 * t402 + (t479 * t466 + (-t403 + t480) * t468) * t466) * qJD(3)) * t493 / 0.2e1 + qJD(1) * ((t452 * t432 + t451 * t433) * qJD(1) + ((t405 * t452 + t407 * t451) * t466 - (t404 * t452 + t406 * t451) * t468) * qJD(3)) / 0.2e1 + m(7) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(6) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(5) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(4) * (t346 ^ 2 + t347 ^ 2 + t360 ^ 2) / 0.2e1 + m(3) * (t381 ^ 2 + t382 ^ 2) / 0.2e1 + (Icges(2,3) + Icges(3,2) * t463 ^ 2 + (Icges(3,1) * t462 + 0.2e1 * Icges(3,4) * t463) * t462 + m(2) * (t444 ^ 2 + t445 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
