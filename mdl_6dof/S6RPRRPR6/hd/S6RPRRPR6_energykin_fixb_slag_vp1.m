% Calculate kinetic energy for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:19
% EndTime: 2019-03-09 05:15:22
% DurationCPUTime: 2.60s
% Computational Cost: add. (2011->297), mult. (1952->461), div. (0->0), fcn. (1942->12), ass. (0->148)
t517 = -Icges(6,3) - Icges(5,3);
t457 = qJ(4) + pkin(11);
t448 = sin(t457);
t450 = cos(t457);
t465 = cos(qJ(1));
t498 = t450 * t465;
t456 = pkin(10) + qJ(3);
t449 = cos(t456);
t463 = sin(qJ(1));
t501 = t449 * t463;
t413 = -t448 * t501 - t498;
t499 = t450 * t463;
t414 = -t448 * t465 + t449 * t499;
t464 = cos(qJ(4));
t494 = t464 * t465;
t462 = sin(qJ(4));
t497 = t462 * t463;
t420 = -t449 * t497 - t494;
t495 = t463 * t464;
t496 = t462 * t465;
t421 = t449 * t495 - t496;
t447 = sin(t456);
t503 = t447 * t463;
t516 = Icges(5,5) * t421 + Icges(6,5) * t414 + Icges(5,6) * t420 + Icges(6,6) * t413 - t503 * t517;
t500 = t449 * t465;
t415 = -t448 * t500 + t499;
t416 = t448 * t463 + t449 * t498;
t422 = -t449 * t496 + t495;
t423 = t449 * t494 + t497;
t502 = t447 * t465;
t515 = Icges(5,5) * t423 + Icges(6,5) * t416 + Icges(5,6) * t422 + Icges(6,6) * t415 - t502 * t517;
t514 = t517 * t449 + (Icges(5,5) * t464 + Icges(6,5) * t450 - Icges(5,6) * t462 - Icges(6,6) * t448) * t447;
t459 = cos(pkin(10));
t508 = pkin(2) * t459;
t507 = t464 * pkin(4);
t505 = Icges(4,4) * t447;
t504 = Icges(4,4) * t449;
t439 = pkin(1) * t463 - qJ(2) * t465;
t492 = pkin(7) * t465 - t463 * t508 - t439;
t483 = pkin(3) * t449 + pkin(8) * t447;
t418 = t483 * t463;
t419 = t483 * t465;
t452 = qJD(3) * t463;
t488 = qJD(3) * t465;
t491 = t418 * t452 + t419 * t488;
t490 = pkin(5) * t450;
t487 = qJD(4) * t447;
t425 = t465 * t487 + t452;
t486 = qJD(5) * t447;
t485 = qJD(6) * t447;
t484 = pkin(5) * t448;
t426 = t463 * t487 - t488;
t433 = qJD(1) * (pkin(1) * t465 + qJ(2) * t463);
t482 = -qJD(2) * t465 + qJD(1) * (pkin(7) * t463 + t465 * t508) + t433;
t458 = sin(pkin(10));
t481 = rSges(3,1) * t459 - rSges(3,2) * t458;
t480 = rSges(4,1) * t449 - rSges(4,2) * t447;
t479 = Icges(4,1) * t449 - t505;
t478 = -Icges(4,2) * t447 + t504;
t477 = Icges(4,5) * t449 - Icges(4,6) * t447;
t405 = -Icges(4,6) * t465 + t463 * t478;
t407 = -Icges(4,5) * t465 + t463 * t479;
t476 = t405 * t447 - t407 * t449;
t406 = Icges(4,6) * t463 + t465 * t478;
t408 = Icges(4,5) * t463 + t465 * t479;
t475 = -t406 * t447 + t408 * t449;
t428 = Icges(4,2) * t449 + t505;
t429 = Icges(4,1) * t447 + t504;
t474 = -t428 * t447 + t429 * t449;
t472 = qJ(5) * t447 + t449 * t507;
t368 = -pkin(4) * t496 + t463 * t472;
t473 = -qJD(5) * t449 + t425 * t368 + t491;
t471 = pkin(9) * t447 + t449 * t490;
t431 = pkin(3) * t447 - pkin(8) * t449;
t470 = qJD(1) * t419 - t431 * t452 + t482;
t453 = qJD(2) * t463;
t469 = t453 + (-t418 + t492) * qJD(1) - t431 * t488;
t369 = pkin(4) * t497 + t465 * t472;
t442 = -qJD(4) * t449 + qJD(1);
t468 = t442 * t369 + t463 * t486 + t470;
t382 = -qJ(5) * t449 + t447 * t507;
t467 = t426 * t382 + t465 * t486 + t469;
t451 = qJ(6) + t457;
t444 = cos(t451);
t443 = sin(t451);
t441 = rSges(2,1) * t465 - rSges(2,2) * t463;
t440 = rSges(2,1) * t463 + rSges(2,2) * t465;
t430 = rSges(4,1) * t447 + rSges(4,2) * t449;
t427 = Icges(4,5) * t447 + Icges(4,6) * t449;
t424 = qJD(1) + (-qJD(4) - qJD(6)) * t449;
t410 = rSges(4,3) * t463 + t465 * t480;
t409 = -rSges(4,3) * t465 + t463 * t480;
t404 = Icges(4,3) * t463 + t465 * t477;
t403 = -Icges(4,3) * t465 + t463 * t477;
t402 = t443 * t463 + t444 * t500;
t401 = -t443 * t500 + t444 * t463;
t400 = -t443 * t465 + t444 * t501;
t399 = -t443 * t501 - t444 * t465;
t398 = t463 * t485 + t426;
t397 = t465 * t485 + t425;
t395 = -rSges(5,3) * t449 + (rSges(5,1) * t464 - rSges(5,2) * t462) * t447;
t394 = -Icges(5,5) * t449 + (Icges(5,1) * t464 - Icges(5,4) * t462) * t447;
t393 = -Icges(5,6) * t449 + (Icges(5,4) * t464 - Icges(5,2) * t462) * t447;
t390 = -rSges(6,3) * t449 + (rSges(6,1) * t450 - rSges(6,2) * t448) * t447;
t389 = -Icges(6,5) * t449 + (Icges(6,1) * t450 - Icges(6,4) * t448) * t447;
t388 = -Icges(6,6) * t449 + (Icges(6,4) * t450 - Icges(6,2) * t448) * t447;
t386 = -rSges(7,3) * t449 + (rSges(7,1) * t444 - rSges(7,2) * t443) * t447;
t385 = -Icges(7,5) * t449 + (Icges(7,1) * t444 - Icges(7,4) * t443) * t447;
t384 = -Icges(7,6) * t449 + (Icges(7,4) * t444 - Icges(7,2) * t443) * t447;
t383 = -Icges(7,3) * t449 + (Icges(7,5) * t444 - Icges(7,6) * t443) * t447;
t381 = qJD(1) * t463 * rSges(3,3) + t433 + (qJD(1) * t481 - qJD(2)) * t465;
t380 = t453 + (t465 * rSges(3,3) - t463 * t481 - t439) * qJD(1);
t379 = -pkin(9) * t449 + t447 * t490;
t378 = rSges(5,1) * t423 + rSges(5,2) * t422 + rSges(5,3) * t502;
t377 = rSges(5,1) * t421 + rSges(5,2) * t420 + rSges(5,3) * t503;
t376 = Icges(5,1) * t423 + Icges(5,4) * t422 + Icges(5,5) * t502;
t375 = Icges(5,1) * t421 + Icges(5,4) * t420 + Icges(5,5) * t503;
t374 = Icges(5,4) * t423 + Icges(5,2) * t422 + Icges(5,6) * t502;
t373 = Icges(5,4) * t421 + Icges(5,2) * t420 + Icges(5,6) * t503;
t367 = (t409 * t463 + t410 * t465) * qJD(3);
t366 = rSges(6,1) * t416 + rSges(6,2) * t415 + rSges(6,3) * t502;
t365 = rSges(6,1) * t414 + rSges(6,2) * t413 + rSges(6,3) * t503;
t364 = Icges(6,1) * t416 + Icges(6,4) * t415 + Icges(6,5) * t502;
t363 = Icges(6,1) * t414 + Icges(6,4) * t413 + Icges(6,5) * t503;
t362 = Icges(6,4) * t416 + Icges(6,2) * t415 + Icges(6,6) * t502;
t361 = Icges(6,4) * t414 + Icges(6,2) * t413 + Icges(6,6) * t503;
t357 = rSges(7,1) * t402 + rSges(7,2) * t401 + rSges(7,3) * t502;
t356 = rSges(7,1) * t400 + rSges(7,2) * t399 + rSges(7,3) * t503;
t355 = Icges(7,1) * t402 + Icges(7,4) * t401 + Icges(7,5) * t502;
t354 = Icges(7,1) * t400 + Icges(7,4) * t399 + Icges(7,5) * t503;
t353 = Icges(7,4) * t402 + Icges(7,2) * t401 + Icges(7,6) * t502;
t352 = Icges(7,4) * t400 + Icges(7,2) * t399 + Icges(7,6) * t503;
t351 = Icges(7,5) * t402 + Icges(7,6) * t401 + Icges(7,3) * t502;
t350 = Icges(7,5) * t400 + Icges(7,6) * t399 + Icges(7,3) * t503;
t348 = t463 * t484 + t465 * t471;
t347 = t463 * t471 - t465 * t484;
t346 = qJD(1) * t410 - t430 * t452 + t482;
t345 = -t430 * t488 + t453 + (-t409 + t492) * qJD(1);
t344 = t377 * t425 - t378 * t426 + t491;
t343 = t378 * t442 - t395 * t425 + t470;
t342 = -t377 * t442 + t395 * t426 + t469;
t341 = t366 * t442 + (-t382 - t390) * t425 + t468;
t340 = t390 * t426 + (-t365 - t368) * t442 + t467;
t339 = t365 * t425 + (-t366 - t369) * t426 + t473;
t338 = t348 * t442 + t357 * t424 - t386 * t397 + (-t379 - t382) * t425 + t468;
t337 = -t356 * t424 + t379 * t426 + t386 * t398 + (-t347 - t368) * t442 + t467;
t336 = t347 * t425 + t356 * t397 - t357 * t398 + (-t348 - t369) * t426 + t473;
t1 = ((t463 * t427 + t465 * t474) * qJD(1) + (t463 ^ 2 * t404 + (t476 * t465 + (-t403 + t475) * t463) * t465) * qJD(3)) * t452 / 0.2e1 + t397 * ((t351 * t502 + t401 * t353 + t402 * t355) * t397 + (t350 * t502 + t352 * t401 + t354 * t402) * t398 + (t383 * t502 + t384 * t401 + t385 * t402) * t424) / 0.2e1 + t398 * ((t351 * t503 + t353 * t399 + t355 * t400) * t397 + (t350 * t503 + t399 * t352 + t400 * t354) * t398 + (t383 * t503 + t384 * t399 + t385 * t400) * t424) / 0.2e1 + t424 * ((-t350 * t398 - t351 * t397 - t383 * t424) * t449 + ((-t353 * t443 + t355 * t444) * t397 + (-t352 * t443 + t354 * t444) * t398 + (-t384 * t443 + t385 * t444) * t424) * t447) / 0.2e1 - ((-t465 * t427 + t463 * t474) * qJD(1) + (t465 ^ 2 * t403 + (t475 * t463 + (-t404 + t476) * t465) * t463) * qJD(3)) * t488 / 0.2e1 + qJD(1) * ((t449 * t428 + t447 * t429) * qJD(1) + ((t406 * t449 + t408 * t447) * t463 - (t405 * t449 + t447 * t407) * t465) * qJD(3)) / 0.2e1 + m(6) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(7) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(5) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(4) * (t345 ^ 2 + t346 ^ 2 + t367 ^ 2) / 0.2e1 + m(3) * (t380 ^ 2 + t381 ^ 2) / 0.2e1 + ((t388 * t415 + t389 * t416 + t393 * t422 + t394 * t423 + t502 * t514) * t442 + (t361 * t415 + t363 * t416 + t373 * t422 + t375 * t423 + t502 * t516) * t426 + (t415 * t362 + t416 * t364 + t422 * t374 + t423 * t376 + t515 * t502) * t425) * t425 / 0.2e1 + ((t388 * t413 + t389 * t414 + t393 * t420 + t394 * t421 + t503 * t514) * t442 + (t413 * t361 + t414 * t363 + t420 * t373 + t421 * t375 + t516 * t503) * t426 + (t362 * t413 + t364 * t414 + t374 * t420 + t376 * t421 + t503 * t515) * t425) * t426 / 0.2e1 + ((-t515 * t425 - t426 * t516 - t514 * t442) * t449 + ((-t388 * t448 + t389 * t450 - t393 * t462 + t394 * t464) * t442 + (-t361 * t448 + t363 * t450 - t373 * t462 + t375 * t464) * t426 + (-t362 * t448 + t364 * t450 - t374 * t462 + t376 * t464) * t425) * t447) * t442 / 0.2e1 + (Icges(3,2) * t459 ^ 2 + (Icges(3,1) * t458 + 0.2e1 * Icges(3,4) * t459) * t458 + Icges(2,3) + m(2) * (t440 ^ 2 + t441 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
