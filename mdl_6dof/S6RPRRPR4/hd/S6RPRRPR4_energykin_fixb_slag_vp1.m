% Calculate kinetic energy for
% S6RPRRPR4
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:09:01
% EndTime: 2019-03-09 05:09:03
% DurationCPUTime: 1.83s
% Computational Cost: add. (1864->280), mult. (1624->440), div. (0->0), fcn. (1558->12), ass. (0->149)
t452 = cos(pkin(10));
t507 = t452 * pkin(2);
t451 = cos(pkin(11));
t506 = pkin(5) * t451;
t448 = pkin(10) + qJ(3);
t438 = sin(t448);
t505 = Icges(4,4) * t438;
t440 = cos(t448);
t504 = Icges(4,4) * t440;
t441 = qJ(4) + t448;
t433 = sin(t441);
t503 = Icges(5,4) * t433;
t434 = cos(t441);
t502 = Icges(5,4) * t434;
t455 = sin(qJ(1));
t501 = t433 * t455;
t456 = cos(qJ(1));
t500 = t433 * t456;
t447 = pkin(11) + qJ(6);
t437 = sin(t447);
t499 = t437 * t455;
t498 = t437 * t456;
t439 = cos(t447);
t497 = t439 * t455;
t496 = t439 * t456;
t449 = sin(pkin(11));
t495 = t449 * t455;
t494 = t449 * t456;
t493 = t451 * t455;
t492 = t451 * t456;
t487 = pkin(3) * t440;
t368 = -pkin(8) * t456 + t455 * t487;
t369 = pkin(8) * t455 + t456 * t487;
t443 = qJD(3) * t455;
t485 = qJD(3) * t456;
t489 = t368 * t443 + t369 * t485;
t430 = pkin(1) * t455 - qJ(2) * t456;
t488 = pkin(7) * t456 - t455 * t507 - t430;
t428 = qJD(4) * t455 + t443;
t484 = qJD(5) * t433;
t483 = qJD(6) * t433;
t482 = pkin(3) * qJD(3) * t438;
t481 = -t368 + t488;
t429 = (-qJD(3) - qJD(4)) * t456;
t475 = pkin(4) * t434 + qJ(5) * t433;
t406 = t475 * t455;
t480 = -t406 + t481;
t424 = qJD(1) * (pkin(1) * t456 + qJ(2) * t455);
t479 = -qJD(2) * t456 + qJD(1) * (pkin(7) * t455 + t456 * t507) + t424;
t450 = sin(pkin(10));
t478 = rSges(3,1) * t452 - rSges(3,2) * t450;
t477 = rSges(4,1) * t440 - rSges(4,2) * t438;
t476 = rSges(5,1) * t434 - rSges(5,2) * t433;
t474 = Icges(4,1) * t440 - t505;
t473 = Icges(5,1) * t434 - t503;
t472 = -Icges(4,2) * t438 + t504;
t471 = -Icges(5,2) * t433 + t502;
t470 = Icges(4,5) * t440 - Icges(4,6) * t438;
t469 = Icges(5,5) * t434 - Icges(5,6) * t433;
t395 = -Icges(4,6) * t456 + t455 * t472;
t397 = -Icges(4,5) * t456 + t455 * t474;
t468 = t395 * t438 - t397 * t440;
t396 = Icges(4,6) * t455 + t456 * t472;
t398 = Icges(4,5) * t455 + t456 * t474;
t467 = -t396 * t438 + t398 * t440;
t420 = Icges(4,2) * t440 + t505;
t421 = Icges(4,1) * t438 + t504;
t466 = -t420 * t438 + t421 * t440;
t444 = qJD(2) * t455;
t465 = -t456 * t482 + t444;
t464 = -qJD(5) * t434 + t428 * t406 + t489;
t417 = pkin(4) * t433 - qJ(5) * t434;
t463 = t429 * t417 + t456 * t484 + t465;
t462 = (Icges(5,5) * t433 + Icges(5,6) * t434) * qJD(1) + (-Icges(5,3) * t456 + t455 * t469) * t429 + (Icges(5,3) * t455 + t456 * t469) * t428;
t461 = pkin(9) * t433 + t434 * t506;
t460 = qJD(1) * t369 - t455 * t482 + t479;
t407 = t475 * t456;
t459 = qJD(1) * t407 + t455 * t484 + t460;
t386 = -Icges(5,6) * t456 + t455 * t471;
t387 = Icges(5,6) * t455 + t456 * t471;
t388 = -Icges(5,5) * t456 + t455 * t473;
t389 = Icges(5,5) * t455 + t456 * t473;
t415 = Icges(5,2) * t434 + t503;
t416 = Icges(5,1) * t433 + t502;
t458 = (-t387 * t433 + t389 * t434) * t428 + (-t386 * t433 + t388 * t434) * t429 + (-t415 * t433 + t416 * t434) * qJD(1);
t432 = rSges(2,1) * t456 - rSges(2,2) * t455;
t431 = rSges(2,1) * t455 + rSges(2,2) * t456;
t427 = -qJD(6) * t434 + qJD(1);
t422 = rSges(4,1) * t438 + rSges(4,2) * t440;
t419 = Icges(4,5) * t438 + Icges(4,6) * t440;
t418 = rSges(5,1) * t433 + rSges(5,2) * t434;
t413 = t434 * t492 + t495;
t412 = -t434 * t494 + t493;
t411 = t434 * t493 - t494;
t410 = -t434 * t495 - t492;
t409 = t455 * t483 + t429;
t408 = t456 * t483 + t428;
t405 = t434 * t496 + t499;
t404 = -t434 * t498 + t497;
t403 = t434 * t497 - t498;
t402 = -t434 * t499 - t496;
t401 = rSges(4,3) * t455 + t456 * t477;
t400 = -rSges(4,3) * t456 + t455 * t477;
t394 = Icges(4,3) * t455 + t456 * t470;
t393 = -Icges(4,3) * t456 + t455 * t470;
t391 = rSges(5,3) * t455 + t456 * t476;
t390 = -rSges(5,3) * t456 + t455 * t476;
t381 = -rSges(6,3) * t434 + (rSges(6,1) * t451 - rSges(6,2) * t449) * t433;
t380 = -Icges(6,5) * t434 + (Icges(6,1) * t451 - Icges(6,4) * t449) * t433;
t379 = -Icges(6,6) * t434 + (Icges(6,4) * t451 - Icges(6,2) * t449) * t433;
t378 = -Icges(6,3) * t434 + (Icges(6,5) * t451 - Icges(6,6) * t449) * t433;
t376 = -rSges(7,3) * t434 + (rSges(7,1) * t439 - rSges(7,2) * t437) * t433;
t375 = -Icges(7,5) * t434 + (Icges(7,1) * t439 - Icges(7,4) * t437) * t433;
t374 = -Icges(7,6) * t434 + (Icges(7,4) * t439 - Icges(7,2) * t437) * t433;
t373 = -Icges(7,3) * t434 + (Icges(7,5) * t439 - Icges(7,6) * t437) * t433;
t372 = qJD(1) * t455 * rSges(3,3) + t424 + (qJD(1) * t478 - qJD(2)) * t456;
t371 = t444 + (t456 * rSges(3,3) - t455 * t478 - t430) * qJD(1);
t370 = -pkin(9) * t434 + t433 * t506;
t364 = rSges(6,1) * t413 + rSges(6,2) * t412 + rSges(6,3) * t500;
t363 = rSges(6,1) * t411 + rSges(6,2) * t410 + rSges(6,3) * t501;
t362 = Icges(6,1) * t413 + Icges(6,4) * t412 + Icges(6,5) * t500;
t361 = Icges(6,1) * t411 + Icges(6,4) * t410 + Icges(6,5) * t501;
t360 = Icges(6,4) * t413 + Icges(6,2) * t412 + Icges(6,6) * t500;
t359 = Icges(6,4) * t411 + Icges(6,2) * t410 + Icges(6,6) * t501;
t358 = Icges(6,5) * t413 + Icges(6,6) * t412 + Icges(6,3) * t500;
t357 = Icges(6,5) * t411 + Icges(6,6) * t410 + Icges(6,3) * t501;
t356 = pkin(5) * t495 + t456 * t461;
t355 = -pkin(5) * t494 + t455 * t461;
t354 = (t400 * t455 + t401 * t456) * qJD(3);
t353 = rSges(7,1) * t405 + rSges(7,2) * t404 + rSges(7,3) * t500;
t352 = rSges(7,1) * t403 + rSges(7,2) * t402 + rSges(7,3) * t501;
t351 = Icges(7,1) * t405 + Icges(7,4) * t404 + Icges(7,5) * t500;
t350 = Icges(7,1) * t403 + Icges(7,4) * t402 + Icges(7,5) * t501;
t349 = Icges(7,4) * t405 + Icges(7,2) * t404 + Icges(7,6) * t500;
t348 = Icges(7,4) * t403 + Icges(7,2) * t402 + Icges(7,6) * t501;
t347 = Icges(7,5) * t405 + Icges(7,6) * t404 + Icges(7,3) * t500;
t346 = Icges(7,5) * t403 + Icges(7,6) * t402 + Icges(7,3) * t501;
t345 = qJD(1) * t401 - t422 * t443 + t479;
t344 = -t422 * t485 + t444 + (-t400 + t488) * qJD(1);
t343 = qJD(1) * t391 - t418 * t428 + t460;
t342 = t418 * t429 + (-t390 + t481) * qJD(1) + t465;
t341 = t390 * t428 - t391 * t429 + t489;
t340 = qJD(1) * t364 + (-t381 - t417) * t428 + t459;
t339 = t381 * t429 + (-t363 + t480) * qJD(1) + t463;
t338 = t363 * t428 + (-t364 - t407) * t429 + t464;
t337 = qJD(1) * t356 + t353 * t427 - t376 * t408 + (-t370 - t417) * t428 + t459;
t336 = -t352 * t427 + t370 * t429 + t376 * t409 + (-t355 + t480) * qJD(1) + t463;
t335 = t352 * t408 - t353 * t409 + t355 * t428 + (-t356 - t407) * t429 + t464;
t1 = ((t455 * t419 + t456 * t466) * qJD(1) + (t455 ^ 2 * t394 + (t468 * t456 + (-t393 + t467) * t455) * t456) * qJD(3)) * t443 / 0.2e1 + t408 * ((t347 * t500 + t404 * t349 + t405 * t351) * t408 + (t346 * t500 + t348 * t404 + t350 * t405) * t409 + (t373 * t500 + t374 * t404 + t375 * t405) * t427) / 0.2e1 + t409 * ((t347 * t501 + t349 * t402 + t351 * t403) * t408 + (t346 * t501 + t348 * t402 + t350 * t403) * t409 + (t373 * t501 + t374 * t402 + t375 * t403) * t427) / 0.2e1 + t427 * ((-t346 * t409 - t347 * t408 - t373 * t427) * t434 + ((-t349 * t437 + t351 * t439) * t408 + (-t348 * t437 + t350 * t439) * t409 + (-t374 * t437 + t375 * t439) * t427) * t433) / 0.2e1 - ((-t456 * t419 + t455 * t466) * qJD(1) + (t456 ^ 2 * t393 + (t467 * t455 + (-t394 + t468) * t456) * t455) * qJD(3)) * t485 / 0.2e1 + m(5) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(6) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(7) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(4) * (t344 ^ 2 + t345 ^ 2 + t354 ^ 2) / 0.2e1 + m(3) * (t371 ^ 2 + t372 ^ 2) / 0.2e1 + (t462 * t455 + t458 * t456 + (t358 * t500 + t360 * t412 + t362 * t413) * t428 + (t357 * t500 + t359 * t412 + t361 * t413) * t429 + (t378 * t500 + t379 * t412 + t380 * t413) * qJD(1)) * t428 / 0.2e1 + (t458 * t455 - t462 * t456 + (t358 * t501 + t360 * t410 + t362 * t411) * t428 + (t357 * t501 + t359 * t410 + t361 * t411) * t429 + (t378 * t501 + t379 * t410 + t380 * t411) * qJD(1)) * t429 / 0.2e1 + (Icges(2,3) + Icges(3,2) * t452 ^ 2 + (Icges(3,1) * t450 + 0.2e1 * Icges(3,4) * t452) * t450 + m(2) * (t431 ^ 2 + t432 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((-t357 * t429 - t358 * t428) * t434 + ((-t360 * t449 + t362 * t451) * t428 + (-t359 * t449 + t361 * t451) * t429) * t433 + (t387 * t434 + t389 * t433) * t428 + (t386 * t434 + t388 * t433) * t429 + ((t396 * t440 + t398 * t438) * t455 - (t395 * t440 + t397 * t438) * t456) * qJD(3) + (t440 * t420 + t438 * t421 + (-t378 + t415) * t434 + (-t379 * t449 + t380 * t451 + t416) * t433) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
