% Calculate kinetic energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:44
% EndTime: 2019-03-09 02:46:46
% DurationCPUTime: 1.92s
% Computational Cost: add. (1551->258), mult. (2145->385), div. (0->0), fcn. (2299->10), ass. (0->130)
t496 = Icges(5,1) + Icges(6,1);
t495 = Icges(5,4) - Icges(6,5);
t494 = Icges(6,4) + Icges(5,5);
t493 = Icges(5,2) + Icges(6,3);
t492 = -Icges(6,6) + Icges(5,6);
t491 = -Icges(5,3) - Icges(6,2);
t430 = pkin(9) + qJ(3);
t428 = cos(t430);
t433 = cos(pkin(10));
t439 = cos(qJ(1));
t472 = t433 * t439;
t431 = sin(pkin(10));
t437 = sin(qJ(1));
t474 = t431 * t437;
t405 = t428 * t474 + t472;
t471 = t437 * t433;
t473 = t431 * t439;
t406 = t428 * t471 - t473;
t427 = sin(t430);
t476 = t427 * t437;
t490 = -t493 * t405 + t495 * t406 + t492 * t476;
t407 = t428 * t473 - t471;
t408 = t428 * t472 + t474;
t475 = t427 * t439;
t489 = t493 * t407 - t495 * t408 - t492 * t475;
t488 = t492 * t405 - t494 * t406 + t491 * t476;
t487 = -t492 * t407 + t494 * t408 - t491 * t475;
t486 = t495 * t405 - t496 * t406 - t494 * t476;
t485 = -t495 * t407 + t496 * t408 + t494 * t475;
t484 = t492 * t428 + (t493 * t431 - t495 * t433) * t427;
t483 = t491 * t428 + (-t492 * t431 + t494 * t433) * t427;
t482 = -t494 * t428 + (-t495 * t431 + t496 * t433) * t427;
t434 = cos(pkin(9));
t479 = pkin(2) * t434;
t478 = Icges(4,4) * t427;
t477 = Icges(4,4) * t428;
t421 = pkin(1) * t437 - qJ(2) * t439;
t469 = pkin(7) * t439 - t437 * t479 - t421;
t415 = pkin(3) * t427 - qJ(4) * t428;
t468 = -(pkin(4) * t433 + qJ(5) * t431) * t427 - t415;
t429 = qJD(2) * t437;
t464 = qJD(4) * t427;
t467 = t439 * t464 + t429;
t466 = qJD(3) * t437;
t465 = qJD(3) * t439;
t463 = qJD(6) * t427;
t450 = pkin(3) * t428 + qJ(4) * t427;
t403 = t450 * t437;
t462 = -t403 + t469;
t461 = qJD(5) * t407 + t467;
t458 = qJD(3) * (rSges(5,3) * t428 - (rSges(5,1) * t433 - rSges(5,2) * t431) * t427 - t415);
t371 = pkin(4) * t406 + qJ(5) * t405;
t457 = -t371 + t462;
t456 = qJD(3) * (rSges(6,2) * t428 - (rSges(6,1) * t433 + rSges(6,3) * t431) * t427 + t468);
t455 = qJD(3) * (-pkin(5) * t427 * t433 - pkin(8) * t428 + t468);
t417 = qJD(1) * (pkin(1) * t439 + qJ(2) * t437);
t454 = -qJD(2) * t439 + qJD(1) * (pkin(7) * t437 + t439 * t479) + t417;
t404 = t450 * t439;
t453 = -qJD(4) * t428 + t403 * t466 + t404 * t465;
t432 = sin(pkin(9));
t452 = rSges(3,1) * t434 - rSges(3,2) * t432;
t451 = rSges(4,1) * t428 - rSges(4,2) * t427;
t449 = Icges(4,1) * t428 - t478;
t448 = -Icges(4,2) * t427 + t477;
t447 = Icges(4,5) * t428 - Icges(4,6) * t427;
t389 = -Icges(4,6) * t439 + t437 * t448;
t391 = -Icges(4,5) * t439 + t437 * t449;
t446 = t389 * t427 - t391 * t428;
t390 = Icges(4,6) * t437 + t439 * t448;
t392 = Icges(4,5) * t437 + t439 * t449;
t445 = -t390 * t427 + t392 * t428;
t413 = Icges(4,2) * t428 + t478;
t414 = Icges(4,1) * t427 + t477;
t444 = -t413 * t427 + t414 * t428;
t443 = qJD(1) * t404 + t437 * t464 + t454;
t372 = pkin(4) * t408 + qJ(5) * t407;
t442 = qJD(5) * t427 * t431 + t371 * t466 + t372 * t465 + t453;
t441 = qJD(1) * t372 + qJD(5) * t405 + t443;
t438 = cos(qJ(6));
t436 = sin(qJ(6));
t424 = qJD(6) * t428 + qJD(1);
t423 = rSges(2,1) * t439 - rSges(2,2) * t437;
t422 = rSges(2,1) * t437 + rSges(2,2) * t439;
t416 = rSges(4,1) * t427 + rSges(4,2) * t428;
t412 = Icges(4,5) * t427 + Icges(4,6) * t428;
t411 = -t437 * t463 - t465;
t410 = -t439 * t463 + t466;
t398 = (t431 * t436 + t433 * t438) * t427;
t397 = (t431 * t438 - t433 * t436) * t427;
t396 = rSges(4,3) * t437 + t439 * t451;
t395 = -rSges(4,3) * t439 + t437 * t451;
t388 = Icges(4,3) * t437 + t439 * t447;
t387 = -Icges(4,3) * t439 + t437 * t447;
t376 = pkin(5) * t408 - pkin(8) * t475;
t375 = pkin(5) * t406 - pkin(8) * t476;
t374 = qJD(1) * t437 * rSges(3,3) + t417 + (qJD(1) * t452 - qJD(2)) * t439;
t373 = t429 + (t439 * rSges(3,3) - t437 * t452 - t421) * qJD(1);
t369 = t407 * t436 + t408 * t438;
t368 = t407 * t438 - t408 * t436;
t367 = t405 * t436 + t406 * t438;
t366 = t405 * t438 - t406 * t436;
t363 = rSges(5,1) * t408 - rSges(5,2) * t407 + rSges(5,3) * t475;
t362 = rSges(6,1) * t408 + rSges(6,2) * t475 + rSges(6,3) * t407;
t361 = rSges(5,1) * t406 - rSges(5,2) * t405 + rSges(5,3) * t476;
t360 = rSges(6,1) * t406 + rSges(6,2) * t476 + rSges(6,3) * t405;
t347 = rSges(7,1) * t398 + rSges(7,2) * t397 + rSges(7,3) * t428;
t346 = (t395 * t437 + t396 * t439) * qJD(3);
t345 = Icges(7,1) * t398 + Icges(7,4) * t397 + Icges(7,5) * t428;
t344 = Icges(7,4) * t398 + Icges(7,2) * t397 + Icges(7,6) * t428;
t343 = Icges(7,5) * t398 + Icges(7,6) * t397 + Icges(7,3) * t428;
t342 = qJD(1) * t396 - t416 * t466 + t454;
t341 = -t416 * t465 + t429 + (-t395 + t469) * qJD(1);
t340 = rSges(7,1) * t369 + rSges(7,2) * t368 - rSges(7,3) * t475;
t339 = rSges(7,1) * t367 + rSges(7,2) * t366 - rSges(7,3) * t476;
t338 = Icges(7,1) * t369 + Icges(7,4) * t368 - Icges(7,5) * t475;
t337 = Icges(7,1) * t367 + Icges(7,4) * t366 - Icges(7,5) * t476;
t336 = Icges(7,4) * t369 + Icges(7,2) * t368 - Icges(7,6) * t475;
t335 = Icges(7,4) * t367 + Icges(7,2) * t366 - Icges(7,6) * t476;
t334 = Icges(7,5) * t369 + Icges(7,6) * t368 - Icges(7,3) * t475;
t333 = Icges(7,5) * t367 + Icges(7,6) * t366 - Icges(7,3) * t476;
t332 = (t361 * t437 + t363 * t439) * qJD(3) + t453;
t331 = qJD(1) * t363 + t437 * t458 + t443;
t330 = t439 * t458 + (-t361 + t462) * qJD(1) + t467;
t329 = qJD(1) * t362 + t437 * t456 + t441;
t328 = t439 * t456 + (-t360 + t457) * qJD(1) + t461;
t327 = (t360 * t437 + t362 * t439) * qJD(3) + t442;
t326 = qJD(1) * t376 + t340 * t424 - t347 * t410 + t437 * t455 + t441;
t325 = -t339 * t424 + t347 * t411 + t439 * t455 + (-t375 + t457) * qJD(1) + t461;
t324 = t339 * t410 - t340 * t411 + (t375 * t437 + t376 * t439) * qJD(3) + t442;
t1 = t410 * ((-t334 * t475 + t368 * t336 + t369 * t338) * t410 + (-t333 * t475 + t335 * t368 + t337 * t369) * t411 + (-t343 * t475 + t344 * t368 + t345 * t369) * t424) / 0.2e1 + t411 * ((-t334 * t476 + t336 * t366 + t338 * t367) * t410 + (-t333 * t476 + t366 * t335 + t367 * t337) * t411 + (-t343 * t476 + t344 * t366 + t345 * t367) * t424) / 0.2e1 + t424 * ((t334 * t428 + t336 * t397 + t338 * t398) * t410 + (t333 * t428 + t335 * t397 + t337 * t398) * t411 + (t428 * t343 + t397 * t344 + t398 * t345) * t424) / 0.2e1 + m(5) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(6) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(7) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t341 ^ 2 + t342 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t373 ^ 2 + t374 ^ 2) / 0.2e1 + (Icges(2,3) + Icges(3,2) * t434 ^ 2 + (Icges(3,1) * t432 + 0.2e1 * Icges(3,4) * t434) * t432 + m(2) * (t422 ^ 2 + t423 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((((-t389 - t488) * t439 + (t390 - t487) * t437) * t428 + ((t490 * t431 + t486 * t433 - t391) * t439 + (t489 * t431 + t485 * t433 + t392) * t437) * t427) * qJD(3) + ((t413 - t483) * t428 + (t484 * t431 + t482 * t433 + t414) * t427) * qJD(1)) * qJD(1) / 0.2e1 + (((t490 * t407 + t486 * t408 + t446 * t439 + t488 * t475) * t439 + ((-t387 + t445) * t439 + t437 * t388 + t487 * t475 + t485 * t408 + t489 * t407) * t437) * qJD(3) + (t484 * t407 + t482 * t408 + t437 * t412 + t439 * t444 + t483 * t475) * qJD(1)) * t466 / 0.2e1 - (((t439 * t387 + t490 * t405 + t486 * t406 + t488 * t476) * t439 + (t445 * t437 + (-t388 + t446) * t439 + t487 * t476 + t485 * t406 + t489 * t405) * t437) * qJD(3) + (t484 * t405 + t482 * t406 - t439 * t412 + t437 * t444 + t483 * t476) * qJD(1)) * t465 / 0.2e1;
T  = t1;
