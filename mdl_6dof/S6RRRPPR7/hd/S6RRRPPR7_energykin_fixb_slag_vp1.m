% Calculate kinetic energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:31
% EndTime: 2019-03-09 15:56:33
% DurationCPUTime: 2.91s
% Computational Cost: add. (1513->319), mult. (3221->476), div. (0->0), fcn. (3564->10), ass. (0->157)
t532 = Icges(4,1) + Icges(5,1);
t531 = -Icges(4,4) + Icges(5,5);
t530 = Icges(5,4) + Icges(4,5);
t529 = Icges(4,2) + Icges(5,3);
t528 = -Icges(5,6) + Icges(4,6);
t527 = Icges(6,3) + Icges(4,3) + Icges(5,2);
t469 = sin(qJ(3));
t472 = cos(qJ(3));
t474 = cos(qJ(1));
t471 = sin(qJ(1));
t473 = cos(qJ(2));
t503 = t473 * t471;
t435 = t469 * t503 + t472 * t474;
t436 = -t469 * t474 + t472 * t503;
t470 = sin(qJ(2));
t505 = t470 * t471;
t526 = t529 * t435 + t531 * t436 - t528 * t505;
t502 = t473 * t474;
t437 = t469 * t502 - t471 * t472;
t438 = t469 * t471 + t472 * t502;
t504 = t470 * t474;
t525 = t529 * t437 + t531 * t438 - t528 * t504;
t524 = t531 * t435 + t532 * t436 + t530 * t505;
t523 = t531 * t437 + t532 * t438 + t530 * t504;
t402 = pkin(3) * t438 + qJ(4) * t437;
t458 = -qJD(3) * t473 + qJD(1);
t522 = qJD(4) * t435 + t458 * t402;
t439 = (pkin(3) * t472 + qJ(4) * t469) * t470;
t494 = qJD(3) * t470;
t495 = qJD(2) * t474;
t443 = t471 * t494 - t495;
t521 = qJD(4) * t437 + t443 * t439;
t520 = t528 * t473 + (t529 * t469 + t531 * t472) * t470;
t519 = -t530 * t473 + (t531 * t469 + t532 * t472) * t470;
t466 = sin(pkin(10));
t467 = cos(pkin(10));
t396 = t435 * t467 - t436 * t466;
t508 = t435 * t466;
t397 = t436 * t467 + t508;
t518 = Icges(6,5) * t397 + Icges(6,6) * t396 + t528 * t435 - t530 * t436 - t527 * t505;
t398 = t437 * t467 - t438 * t466;
t507 = t437 * t466;
t399 = t438 * t467 + t507;
t517 = Icges(6,5) * t399 + Icges(6,6) * t398 + t528 * t437 - t530 * t438 - t527 * t504;
t430 = (-t466 * t472 + t467 * t469) * t470;
t506 = t466 * t469;
t431 = (t467 * t472 + t506) * t470;
t516 = Icges(6,5) * t431 + Icges(6,6) * t430 + (t528 * t469 - t530 * t472) * t470 + t527 * t473;
t511 = pkin(5) * t467;
t510 = Icges(3,4) * t470;
t509 = Icges(3,4) * t473;
t401 = pkin(3) * t436 + qJ(4) * t435;
t406 = pkin(4) * t436 - qJ(5) * t505;
t500 = -t401 - t406;
t407 = pkin(4) * t438 - qJ(5) * t504;
t499 = -t402 - t407;
t489 = pkin(2) * t473 + pkin(8) * t470;
t440 = t489 * t471;
t441 = t489 * t474;
t464 = qJD(2) * t471;
t498 = t440 * t464 + t441 * t495;
t446 = qJD(1) * (pkin(1) * t474 + pkin(7) * t471);
t497 = qJD(1) * t441 + t446;
t444 = pkin(4) * t470 * t472 + qJ(5) * t473;
t496 = -t439 - t444;
t442 = t474 * t494 + t464;
t493 = qJD(6) * t470;
t492 = t470 * pkin(9);
t454 = pkin(1) * t471 - pkin(7) * t474;
t491 = (-t440 - t454) * qJD(1);
t490 = qJD(4) * t470 * t469 + t442 * t401 + t498;
t488 = rSges(3,1) * t473 - rSges(3,2) * t470;
t487 = Icges(3,1) * t473 - t510;
t486 = -Icges(3,2) * t470 + t509;
t485 = Icges(3,5) * t473 - Icges(3,6) * t470;
t418 = -Icges(3,6) * t474 + t471 * t486;
t422 = -Icges(3,5) * t474 + t471 * t487;
t484 = t418 * t470 - t422 * t473;
t419 = Icges(3,6) * t471 + t474 * t486;
t423 = Icges(3,5) * t471 + t474 * t487;
t483 = -t419 * t470 + t423 * t473;
t448 = Icges(3,2) * t473 + t510;
t449 = Icges(3,1) * t470 + t509;
t482 = -t448 * t470 + t449 * t473;
t453 = pkin(2) * t470 - pkin(8) * t473;
t481 = -qJD(2) * t453 - qJD(5) * t470;
t480 = -t453 * t464 + t497;
t479 = qJD(5) * t473 + t442 * t406 + t490;
t478 = -t453 * t495 + t491;
t477 = t458 * t407 + t471 * t481 + t497 + t522;
t476 = t443 * t444 + t474 * t481 + t491 + t521;
t465 = pkin(10) + qJ(6);
t462 = cos(t465);
t461 = sin(t465);
t452 = rSges(2,1) * t474 - rSges(2,2) * t471;
t451 = rSges(2,1) * t471 + rSges(2,2) * t474;
t450 = rSges(3,1) * t470 + rSges(3,2) * t473;
t447 = Icges(3,5) * t470 + Icges(3,6) * t473;
t445 = qJD(1) + (-qJD(3) + qJD(6)) * t473;
t427 = rSges(3,3) * t471 + t474 * t488;
t426 = -rSges(3,3) * t474 + t471 * t488;
t425 = -rSges(4,3) * t473 + (rSges(4,1) * t472 - rSges(4,2) * t469) * t470;
t424 = -rSges(5,2) * t473 + (rSges(5,1) * t472 + rSges(5,3) * t469) * t470;
t415 = Icges(3,3) * t471 + t474 * t485;
t414 = -Icges(3,3) * t474 + t471 * t485;
t411 = -t471 * t493 + t443;
t410 = -t474 * t493 + t442;
t409 = (t461 * t469 + t462 * t472) * t470;
t408 = (-t461 * t472 + t462 * t469) * t470;
t395 = pkin(9) * t473 + (pkin(5) * t506 + t472 * t511) * t470;
t394 = rSges(4,1) * t438 - rSges(4,2) * t437 + rSges(4,3) * t504;
t393 = rSges(5,1) * t438 + rSges(5,2) * t504 + rSges(5,3) * t437;
t392 = rSges(4,1) * t436 - rSges(4,2) * t435 + rSges(4,3) * t505;
t391 = rSges(5,1) * t436 + rSges(5,2) * t505 + rSges(5,3) * t435;
t377 = t437 * t461 + t438 * t462;
t376 = t437 * t462 - t438 * t461;
t375 = t435 * t461 + t436 * t462;
t374 = t435 * t462 - t436 * t461;
t373 = rSges(6,1) * t431 + rSges(6,2) * t430 + rSges(6,3) * t473;
t372 = Icges(6,1) * t431 + Icges(6,4) * t430 + Icges(6,5) * t473;
t371 = Icges(6,4) * t431 + Icges(6,2) * t430 + Icges(6,6) * t473;
t369 = qJD(1) * t427 - t450 * t464 + t446;
t368 = -t450 * t495 + (-t426 - t454) * qJD(1);
t366 = (t426 * t471 + t427 * t474) * qJD(2);
t365 = rSges(7,1) * t409 + rSges(7,2) * t408 + rSges(7,3) * t473;
t364 = Icges(7,1) * t409 + Icges(7,4) * t408 + Icges(7,5) * t473;
t363 = Icges(7,4) * t409 + Icges(7,2) * t408 + Icges(7,6) * t473;
t362 = Icges(7,5) * t409 + Icges(7,6) * t408 + Icges(7,3) * t473;
t361 = pkin(5) * t507 + t438 * t511 - t474 * t492;
t360 = pkin(5) * t508 + t436 * t511 - t471 * t492;
t359 = rSges(6,1) * t399 + rSges(6,2) * t398 - rSges(6,3) * t504;
t358 = rSges(6,1) * t397 + rSges(6,2) * t396 - rSges(6,3) * t505;
t357 = Icges(6,1) * t399 + Icges(6,4) * t398 - Icges(6,5) * t504;
t356 = Icges(6,1) * t397 + Icges(6,4) * t396 - Icges(6,5) * t505;
t355 = Icges(6,4) * t399 + Icges(6,2) * t398 - Icges(6,6) * t504;
t354 = Icges(6,4) * t397 + Icges(6,2) * t396 - Icges(6,6) * t505;
t351 = rSges(7,1) * t377 + rSges(7,2) * t376 - rSges(7,3) * t504;
t350 = rSges(7,1) * t375 + rSges(7,2) * t374 - rSges(7,3) * t505;
t349 = Icges(7,1) * t377 + Icges(7,4) * t376 - Icges(7,5) * t504;
t348 = Icges(7,1) * t375 + Icges(7,4) * t374 - Icges(7,5) * t505;
t347 = Icges(7,4) * t377 + Icges(7,2) * t376 - Icges(7,6) * t504;
t346 = Icges(7,4) * t375 + Icges(7,2) * t374 - Icges(7,6) * t505;
t345 = Icges(7,5) * t377 + Icges(7,6) * t376 - Icges(7,3) * t504;
t344 = Icges(7,5) * t375 + Icges(7,6) * t374 - Icges(7,3) * t505;
t343 = t394 * t458 - t425 * t442 + t480;
t342 = -t392 * t458 + t425 * t443 + t478;
t341 = t392 * t442 - t394 * t443 + t498;
t340 = t393 * t458 + (-t424 - t439) * t442 + t480 + t522;
t339 = t424 * t443 + (-t391 - t401) * t458 + t478 + t521;
t338 = t391 * t442 + (-t393 - t402) * t443 + t490;
t337 = t359 * t458 + (-t373 + t496) * t442 + t477;
t336 = t373 * t443 + (-t358 + t500) * t458 + t476;
t335 = t358 * t442 + (-t359 + t499) * t443 + t479;
t334 = t351 * t445 + t361 * t458 - t365 * t410 + (-t395 + t496) * t442 + t477;
t333 = -t350 * t445 + t365 * t411 + t395 * t443 + (-t360 + t500) * t458 + t476;
t332 = t350 * t410 - t351 * t411 + t360 * t442 + (-t361 + t499) * t443 + t479;
t1 = ((t471 * t447 + t474 * t482) * qJD(1) + (t471 ^ 2 * t415 + (t484 * t474 + (-t414 + t483) * t471) * t474) * qJD(2)) * t464 / 0.2e1 - ((-t474 * t447 + t471 * t482) * qJD(1) + (t474 ^ 2 * t414 + (t483 * t471 + (-t415 + t484) * t474) * t471) * qJD(2)) * t495 / 0.2e1 + m(4) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(5) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(6) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(7) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(3) * (t366 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + t411 * ((-t345 * t505 + t347 * t374 + t349 * t375) * t410 + (-t344 * t505 + t374 * t346 + t375 * t348) * t411 + (-t362 * t505 + t363 * t374 + t364 * t375) * t445) / 0.2e1 + t445 * ((t345 * t473 + t347 * t408 + t349 * t409) * t410 + (t344 * t473 + t346 * t408 + t348 * t409) * t411 + (t473 * t362 + t408 * t363 + t409 * t364) * t445) / 0.2e1 + t410 * ((-t345 * t504 + t376 * t347 + t377 * t349) * t410 + (-t344 * t504 + t346 * t376 + t348 * t377) * t411 + (-t362 * t504 + t363 * t376 + t364 * t377) * t445) / 0.2e1 + qJD(1) * ((t473 * t448 + t470 * t449) * qJD(1) + ((t419 * t473 + t423 * t470) * t471 - (t418 * t473 + t422 * t470) * t474) * qJD(2)) / 0.2e1 + (Icges(2,3) + m(2) * (t451 ^ 2 + t452 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t371 * t398 + t372 * t399 + t520 * t437 + t519 * t438 - t516 * t504) * t458 + (t354 * t398 + t356 * t399 + t526 * t437 + t524 * t438 - t518 * t504) * t443 + (t398 * t355 + t399 * t357 + t525 * t437 + t523 * t438 - t517 * t504) * t442) * t442 / 0.2e1 + ((t371 * t396 + t372 * t397 + t520 * t435 + t519 * t436 - t516 * t505) * t458 + (t396 * t354 + t397 * t356 + t526 * t435 + t524 * t436 - t518 * t505) * t443 + (t355 * t396 + t357 * t397 + t525 * t435 + t523 * t436 - t517 * t505) * t442) * t443 / 0.2e1 + ((t355 * t430 + t357 * t431) * t442 + (t354 * t430 + t356 * t431) * t443 + (t430 * t371 + t431 * t372) * t458 + (t517 * t442 + t518 * t443 + t516 * t458) * t473 + ((t520 * t469 + t519 * t472) * t458 + (t526 * t469 + t524 * t472) * t443 + (t525 * t469 + t523 * t472) * t442) * t470) * t458 / 0.2e1;
T  = t1;
