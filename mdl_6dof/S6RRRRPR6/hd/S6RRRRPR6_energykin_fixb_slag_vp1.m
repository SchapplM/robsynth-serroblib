% Calculate kinetic energy for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:38
% EndTime: 2019-03-09 22:18:41
% DurationCPUTime: 3.53s
% Computational Cost: add. (2244->356), mult. (2727->557), div. (0->0), fcn. (2740->12), ass. (0->170)
t547 = -Icges(6,3) - Icges(5,3);
t489 = qJ(3) + qJ(4);
t481 = pkin(11) + t489;
t475 = sin(t481);
t476 = cos(t481);
t495 = cos(qJ(1));
t492 = sin(qJ(1));
t494 = cos(qJ(2));
t531 = t492 * t494;
t418 = -t475 * t531 - t476 * t495;
t419 = -t475 * t495 + t476 * t531;
t483 = sin(t489);
t484 = cos(t489);
t437 = -t483 * t531 - t484 * t495;
t438 = -t483 * t495 + t484 * t531;
t491 = sin(qJ(2));
t533 = t491 * t492;
t546 = Icges(5,5) * t438 + Icges(6,5) * t419 + Icges(5,6) * t437 + Icges(6,6) * t418 - t547 * t533;
t530 = t494 * t495;
t420 = -t475 * t530 + t476 * t492;
t421 = t475 * t492 + t476 * t530;
t439 = -t483 * t530 + t484 * t492;
t440 = t483 * t492 + t484 * t530;
t532 = t491 * t495;
t545 = Icges(5,5) * t440 + Icges(6,5) * t421 + Icges(5,6) * t439 + Icges(6,6) * t420 - t547 * t532;
t544 = t547 * t494 + (Icges(5,5) * t484 + Icges(6,5) * t476 - Icges(5,6) * t483 - Icges(6,6) * t475) * t491;
t493 = cos(qJ(3));
t539 = t493 * pkin(3);
t537 = Icges(3,4) * t491;
t536 = Icges(3,4) * t494;
t490 = sin(qJ(3));
t535 = t490 * t492;
t534 = t490 * t495;
t529 = pkin(5) * t476;
t516 = pkin(2) * t494 + pkin(8) * t491;
t450 = t516 * t492;
t451 = t516 * t495;
t482 = qJD(2) * t492;
t523 = qJD(2) * t495;
t528 = t450 * t482 + t451 * t523;
t527 = pkin(5) * t475;
t526 = pkin(4) * t484;
t522 = qJD(3) * t491;
t452 = t495 * t522 + t482;
t521 = qJD(4) * t491;
t520 = qJD(5) * t491;
t519 = qJD(6) * t491;
t518 = -qJD(3) - qJD(4);
t422 = t495 * t521 + t452;
t517 = pkin(4) * t483;
t453 = t492 * t522 - t523;
t423 = t492 * t521 + t453;
t515 = rSges(3,1) * t494 - rSges(3,2) * t491;
t514 = Icges(3,1) * t494 - t537;
t513 = -Icges(3,2) * t491 + t536;
t512 = Icges(3,5) * t494 - Icges(3,6) * t491;
t429 = -Icges(3,6) * t495 + t492 * t513;
t432 = -Icges(3,5) * t495 + t492 * t514;
t511 = t429 * t491 - t432 * t494;
t430 = Icges(3,6) * t492 + t495 * t513;
t433 = Icges(3,5) * t492 + t495 * t514;
t510 = -t430 * t491 + t433 * t494;
t459 = Icges(3,2) * t494 + t537;
t460 = Icges(3,1) * t491 + t536;
t509 = -t459 * t491 + t460 * t494;
t506 = pkin(9) * t491 + t494 * t539;
t395 = -pkin(3) * t534 + t492 * t506;
t396 = pkin(3) * t535 + t495 * t506;
t508 = t452 * t395 - t396 * t453 + t528;
t456 = qJD(1) * (pkin(1) * t495 + pkin(7) * t492);
t464 = pkin(2) * t491 - pkin(8) * t494;
t507 = qJD(1) * t451 - t464 * t482 + t456;
t465 = pkin(1) * t492 - pkin(7) * t495;
t505 = (-t450 - t465) * qJD(1) - t464 * t523;
t504 = pkin(10) * t491 + t494 * t529;
t503 = qJ(5) * t491 + t494 * t526;
t413 = -pkin(9) * t494 + t491 * t539;
t474 = -qJD(3) * t494 + qJD(1);
t502 = t474 * t396 - t413 * t452 + t507;
t354 = t492 * t503 - t495 * t517;
t501 = -qJD(5) * t494 + t422 * t354 + t508;
t355 = t492 * t517 + t495 * t503;
t455 = t494 * t518 + qJD(1);
t500 = t455 * t355 + t492 * t520 + t502;
t499 = -t395 * t474 + t453 * t413 + t505;
t397 = -qJ(5) * t494 + t491 * t526;
t498 = t423 * t397 + t495 * t520 + t499;
t479 = qJ(6) + t481;
t467 = cos(t479);
t466 = sin(t479);
t463 = rSges(2,1) * t495 - rSges(2,2) * t492;
t462 = rSges(2,1) * t492 + rSges(2,2) * t495;
t461 = rSges(3,1) * t491 + rSges(3,2) * t494;
t458 = Icges(3,5) * t491 + Icges(3,6) * t494;
t449 = t493 * t530 + t535;
t448 = -t490 * t530 + t492 * t493;
t447 = t493 * t531 - t534;
t446 = -t490 * t531 - t493 * t495;
t445 = qJD(1) + (-qJD(6) + t518) * t494;
t436 = rSges(3,3) * t492 + t495 * t515;
t435 = -rSges(3,3) * t495 + t492 * t515;
t434 = -rSges(4,3) * t494 + (rSges(4,1) * t493 - rSges(4,2) * t490) * t491;
t431 = -Icges(4,5) * t494 + (Icges(4,1) * t493 - Icges(4,4) * t490) * t491;
t428 = -Icges(4,6) * t494 + (Icges(4,4) * t493 - Icges(4,2) * t490) * t491;
t427 = Icges(3,3) * t492 + t495 * t512;
t426 = -Icges(3,3) * t495 + t492 * t512;
t425 = -Icges(4,3) * t494 + (Icges(4,5) * t493 - Icges(4,6) * t490) * t491;
t417 = -rSges(5,3) * t494 + (rSges(5,1) * t484 - rSges(5,2) * t483) * t491;
t416 = -Icges(5,5) * t494 + (Icges(5,1) * t484 - Icges(5,4) * t483) * t491;
t415 = -Icges(5,6) * t494 + (Icges(5,4) * t484 - Icges(5,2) * t483) * t491;
t412 = t466 * t492 + t467 * t530;
t411 = -t466 * t530 + t467 * t492;
t410 = -t466 * t495 + t467 * t531;
t409 = -t466 * t531 - t467 * t495;
t408 = -rSges(6,3) * t494 + (rSges(6,1) * t476 - rSges(6,2) * t475) * t491;
t407 = -Icges(6,5) * t494 + (Icges(6,1) * t476 - Icges(6,4) * t475) * t491;
t406 = -Icges(6,6) * t494 + (Icges(6,4) * t476 - Icges(6,2) * t475) * t491;
t404 = t492 * t519 + t423;
t403 = t495 * t519 + t422;
t402 = -rSges(7,3) * t494 + (rSges(7,1) * t467 - rSges(7,2) * t466) * t491;
t401 = -Icges(7,5) * t494 + (Icges(7,1) * t467 - Icges(7,4) * t466) * t491;
t400 = -Icges(7,6) * t494 + (Icges(7,4) * t467 - Icges(7,2) * t466) * t491;
t399 = -Icges(7,3) * t494 + (Icges(7,5) * t467 - Icges(7,6) * t466) * t491;
t394 = rSges(4,1) * t449 + rSges(4,2) * t448 + rSges(4,3) * t532;
t393 = rSges(4,1) * t447 + rSges(4,2) * t446 + rSges(4,3) * t533;
t392 = Icges(4,1) * t449 + Icges(4,4) * t448 + Icges(4,5) * t532;
t391 = Icges(4,1) * t447 + Icges(4,4) * t446 + Icges(4,5) * t533;
t390 = Icges(4,4) * t449 + Icges(4,2) * t448 + Icges(4,6) * t532;
t389 = Icges(4,4) * t447 + Icges(4,2) * t446 + Icges(4,6) * t533;
t388 = Icges(4,5) * t449 + Icges(4,6) * t448 + Icges(4,3) * t532;
t387 = Icges(4,5) * t447 + Icges(4,6) * t446 + Icges(4,3) * t533;
t386 = qJD(1) * t436 - t461 * t482 + t456;
t385 = -t461 * t523 + (-t435 - t465) * qJD(1);
t384 = (t435 * t492 + t436 * t495) * qJD(2);
t382 = rSges(5,1) * t440 + rSges(5,2) * t439 + rSges(5,3) * t532;
t381 = rSges(5,1) * t438 + rSges(5,2) * t437 + rSges(5,3) * t533;
t380 = Icges(5,1) * t440 + Icges(5,4) * t439 + Icges(5,5) * t532;
t379 = Icges(5,1) * t438 + Icges(5,4) * t437 + Icges(5,5) * t533;
t378 = Icges(5,4) * t440 + Icges(5,2) * t439 + Icges(5,6) * t532;
t377 = Icges(5,4) * t438 + Icges(5,2) * t437 + Icges(5,6) * t533;
t374 = -pkin(10) * t494 + t491 * t529;
t372 = rSges(6,1) * t421 + rSges(6,2) * t420 + rSges(6,3) * t532;
t371 = rSges(6,1) * t419 + rSges(6,2) * t418 + rSges(6,3) * t533;
t370 = Icges(6,1) * t421 + Icges(6,4) * t420 + Icges(6,5) * t532;
t369 = Icges(6,1) * t419 + Icges(6,4) * t418 + Icges(6,5) * t533;
t368 = Icges(6,4) * t421 + Icges(6,2) * t420 + Icges(6,6) * t532;
t367 = Icges(6,4) * t419 + Icges(6,2) * t418 + Icges(6,6) * t533;
t363 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t532;
t362 = rSges(7,1) * t410 + rSges(7,2) * t409 + rSges(7,3) * t533;
t361 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t532;
t360 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t533;
t359 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t532;
t358 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t533;
t357 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t532;
t356 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t533;
t351 = t492 * t527 + t495 * t504;
t350 = t492 * t504 - t495 * t527;
t349 = t394 * t474 - t434 * t452 + t507;
t348 = -t393 * t474 + t434 * t453 + t505;
t347 = t393 * t452 - t394 * t453 + t528;
t346 = t382 * t455 - t417 * t422 + t502;
t345 = -t381 * t455 + t417 * t423 + t499;
t344 = t381 * t422 - t382 * t423 + t508;
t343 = t372 * t455 + (-t397 - t408) * t422 + t500;
t342 = t408 * t423 + (-t354 - t371) * t455 + t498;
t341 = t371 * t422 + (-t355 - t372) * t423 + t501;
t340 = t351 * t455 + t363 * t445 - t402 * t403 + (-t374 - t397) * t422 + t500;
t339 = -t362 * t445 + t374 * t423 + t402 * t404 + (-t350 - t354) * t455 + t498;
t338 = t350 * t422 + t362 * t403 - t363 * t404 + (-t351 - t355) * t423 + t501;
t1 = m(6) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(7) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(5) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(4) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(3) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + t404 * ((t357 * t533 + t359 * t409 + t361 * t410) * t403 + (t356 * t533 + t409 * t358 + t410 * t360) * t404 + (t399 * t533 + t400 * t409 + t401 * t410) * t445) / 0.2e1 + t403 * ((t357 * t532 + t411 * t359 + t412 * t361) * t403 + (t356 * t532 + t358 * t411 + t360 * t412) * t404 + (t399 * t532 + t400 * t411 + t401 * t412) * t445) / 0.2e1 + t445 * ((-t356 * t404 - t357 * t403 - t399 * t445) * t494 + ((-t359 * t466 + t361 * t467) * t403 + (-t358 * t466 + t360 * t467) * t404 + (-t400 * t466 + t401 * t467) * t445) * t491) / 0.2e1 + t453 * ((t388 * t533 + t390 * t446 + t392 * t447) * t452 + (t387 * t533 + t446 * t389 + t447 * t391) * t453 + (t425 * t533 + t428 * t446 + t431 * t447) * t474) / 0.2e1 + t474 * ((-t387 * t453 - t388 * t452 - t425 * t474) * t494 + ((-t390 * t490 + t392 * t493) * t452 + (-t389 * t490 + t391 * t493) * t453 + (-t428 * t490 + t431 * t493) * t474) * t491) / 0.2e1 + t452 * ((t388 * t532 + t448 * t390 + t449 * t392) * t452 + (t387 * t532 + t389 * t448 + t391 * t449) * t453 + (t425 * t532 + t428 * t448 + t431 * t449) * t474) / 0.2e1 + qJD(1) * ((t459 * t494 + t460 * t491) * qJD(1) + ((t430 * t494 + t433 * t491) * t492 - (t429 * t494 + t432 * t491) * t495) * qJD(2)) / 0.2e1 - ((-t495 * t458 + t492 * t509) * qJD(1) + (t495 ^ 2 * t426 + (t510 * t492 + (-t427 + t511) * t495) * t492) * qJD(2)) * t523 / 0.2e1 + ((t492 * t458 + t495 * t509) * qJD(1) + (t492 ^ 2 * t427 + (t511 * t495 + (-t426 + t510) * t492) * t495) * qJD(2)) * t482 / 0.2e1 + ((t406 * t420 + t407 * t421 + t415 * t439 + t416 * t440 + t532 * t544) * t455 + (t367 * t420 + t369 * t421 + t377 * t439 + t379 * t440 + t532 * t546) * t423 + (t420 * t368 + t421 * t370 + t439 * t378 + t440 * t380 + t532 * t545) * t422) * t422 / 0.2e1 + ((t406 * t418 + t407 * t419 + t415 * t437 + t416 * t438 + t533 * t544) * t455 + (t418 * t367 + t419 * t369 + t437 * t377 + t438 * t379 + t533 * t546) * t423 + (t368 * t418 + t370 * t419 + t378 * t437 + t380 * t438 + t533 * t545) * t422) * t423 / 0.2e1 + ((-t545 * t422 - t423 * t546 - t544 * t455) * t494 + ((-t406 * t475 + t407 * t476 - t415 * t483 + t416 * t484) * t455 + (-t367 * t475 + t369 * t476 - t377 * t483 + t379 * t484) * t423 + (-t368 * t475 + t370 * t476 - t378 * t483 + t380 * t484) * t422) * t491) * t455 / 0.2e1 + (m(2) * (t462 ^ 2 + t463 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
