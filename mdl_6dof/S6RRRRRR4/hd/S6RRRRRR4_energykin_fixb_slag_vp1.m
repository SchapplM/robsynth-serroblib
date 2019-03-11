% Calculate kinetic energy for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:41
% EndTime: 2019-03-10 03:46:44
% DurationCPUTime: 2.81s
% Computational Cost: add. (2280->362), mult. (2787->590), div. (0->0), fcn. (2800->12), ass. (0->176)
t493 = cos(qJ(3));
t540 = t493 * pkin(3);
t491 = sin(qJ(2));
t538 = Icges(3,4) * t491;
t494 = cos(qJ(2));
t537 = Icges(3,4) * t494;
t490 = sin(qJ(3));
t492 = sin(qJ(1));
t536 = t490 * t492;
t495 = cos(qJ(1));
t535 = t490 * t495;
t534 = t491 * t492;
t533 = t491 * t495;
t532 = t494 * t492;
t531 = t494 * t495;
t489 = qJ(3) + qJ(4);
t484 = qJ(5) + t489;
t478 = cos(t484);
t530 = pkin(5) * t478;
t516 = pkin(2) * t494 + pkin(8) * t491;
t450 = t516 * t492;
t451 = t516 * t495;
t481 = qJD(2) * t492;
t524 = qJD(2) * t495;
t529 = t450 * t481 + t451 * t524;
t477 = sin(t484);
t528 = pkin(5) * t477;
t483 = cos(t489);
t527 = pkin(4) * t483;
t523 = qJD(3) * t491;
t452 = t495 * t523 + t481;
t522 = qJD(4) * t491;
t521 = qJD(5) * t491;
t520 = qJD(6) * t491;
t519 = -qJD(3) - qJD(4);
t422 = t495 * t522 + t452;
t482 = sin(t489);
t518 = pkin(4) * t482;
t517 = -qJD(5) + t519;
t453 = t492 * t523 - t524;
t402 = t495 * t521 + t422;
t423 = t492 * t522 + t453;
t515 = rSges(3,1) * t494 - rSges(3,2) * t491;
t514 = Icges(3,1) * t494 - t538;
t513 = -Icges(3,2) * t491 + t537;
t512 = Icges(3,5) * t494 - Icges(3,6) * t491;
t428 = -Icges(3,6) * t495 + t492 * t513;
t431 = -Icges(3,5) * t495 + t492 * t514;
t511 = t428 * t491 - t431 * t494;
t429 = Icges(3,6) * t492 + t495 * t513;
t432 = Icges(3,5) * t492 + t495 * t514;
t510 = -t429 * t491 + t432 * t494;
t459 = Icges(3,2) * t494 + t538;
t460 = Icges(3,1) * t491 + t537;
t509 = -t459 * t491 + t460 * t494;
t506 = pkin(9) * t491 + t494 * t540;
t392 = -pkin(3) * t535 + t492 * t506;
t393 = pkin(3) * t536 + t495 * t506;
t508 = t452 * t392 - t393 * t453 + t529;
t456 = qJD(1) * (pkin(1) * t495 + pkin(7) * t492);
t464 = pkin(2) * t491 - pkin(8) * t494;
t507 = qJD(1) * t451 - t464 * t481 + t456;
t403 = t492 * t521 + t423;
t465 = pkin(1) * t492 - pkin(7) * t495;
t505 = (-t450 - t465) * qJD(1) - t464 * t524;
t504 = pkin(11) * t491 + t494 * t530;
t503 = pkin(10) * t491 + t494 * t527;
t351 = t492 * t503 - t495 * t518;
t352 = t492 * t518 + t495 * t503;
t502 = t422 * t351 - t352 * t423 + t508;
t408 = -pkin(9) * t494 + t491 * t540;
t474 = -qJD(3) * t494 + qJD(1);
t501 = t474 * t393 - t408 * t452 + t507;
t500 = -t392 * t474 + t453 * t408 + t505;
t394 = -pkin(10) * t494 + t491 * t527;
t455 = t494 * t519 + qJD(1);
t499 = t455 * t352 - t394 * t422 + t501;
t498 = -t351 * t455 + t423 * t394 + t500;
t480 = qJ(6) + t484;
t473 = cos(t480);
t472 = sin(t480);
t463 = rSges(2,1) * t495 - rSges(2,2) * t492;
t462 = rSges(2,1) * t492 + rSges(2,2) * t495;
t461 = rSges(3,1) * t491 + rSges(3,2) * t494;
t458 = Icges(3,5) * t491 + Icges(3,6) * t494;
t449 = t493 * t531 + t536;
t448 = -t490 * t531 + t492 * t493;
t447 = t493 * t532 - t535;
t446 = -t490 * t532 - t493 * t495;
t445 = t494 * t517 + qJD(1);
t440 = t482 * t492 + t483 * t531;
t439 = -t482 * t531 + t483 * t492;
t438 = -t482 * t495 + t483 * t532;
t437 = -t482 * t532 - t483 * t495;
t436 = rSges(3,3) * t492 + t495 * t515;
t435 = -rSges(3,3) * t495 + t492 * t515;
t434 = -rSges(4,3) * t494 + (rSges(4,1) * t493 - rSges(4,2) * t490) * t491;
t430 = -Icges(4,5) * t494 + (Icges(4,1) * t493 - Icges(4,4) * t490) * t491;
t427 = -Icges(4,6) * t494 + (Icges(4,4) * t493 - Icges(4,2) * t490) * t491;
t426 = Icges(3,3) * t492 + t495 * t512;
t425 = -Icges(3,3) * t495 + t492 * t512;
t424 = -Icges(4,3) * t494 + (Icges(4,5) * t493 - Icges(4,6) * t490) * t491;
t421 = qJD(1) + (-qJD(6) + t517) * t494;
t420 = t477 * t492 + t478 * t531;
t419 = -t477 * t531 + t478 * t492;
t418 = -t477 * t495 + t478 * t532;
t417 = -t477 * t532 - t478 * t495;
t416 = -rSges(5,3) * t494 + (rSges(5,1) * t483 - rSges(5,2) * t482) * t491;
t415 = -Icges(5,5) * t494 + (Icges(5,1) * t483 - Icges(5,4) * t482) * t491;
t414 = -Icges(5,6) * t494 + (Icges(5,4) * t483 - Icges(5,2) * t482) * t491;
t413 = -Icges(5,3) * t494 + (Icges(5,5) * t483 - Icges(5,6) * t482) * t491;
t412 = t472 * t492 + t473 * t531;
t411 = -t472 * t531 + t473 * t492;
t410 = -t472 * t495 + t473 * t532;
t409 = -t472 * t532 - t473 * t495;
t407 = -rSges(6,3) * t494 + (rSges(6,1) * t478 - rSges(6,2) * t477) * t491;
t406 = -Icges(6,5) * t494 + (Icges(6,1) * t478 - Icges(6,4) * t477) * t491;
t405 = -Icges(6,6) * t494 + (Icges(6,4) * t478 - Icges(6,2) * t477) * t491;
t404 = -Icges(6,3) * t494 + (Icges(6,5) * t478 - Icges(6,6) * t477) * t491;
t401 = -rSges(7,3) * t494 + (rSges(7,1) * t473 - rSges(7,2) * t472) * t491;
t400 = -Icges(7,5) * t494 + (Icges(7,1) * t473 - Icges(7,4) * t472) * t491;
t399 = -Icges(7,6) * t494 + (Icges(7,4) * t473 - Icges(7,2) * t472) * t491;
t398 = -Icges(7,3) * t494 + (Icges(7,5) * t473 - Icges(7,6) * t472) * t491;
t397 = t492 * t520 + t403;
t396 = t495 * t520 + t402;
t391 = rSges(4,1) * t449 + rSges(4,2) * t448 + rSges(4,3) * t533;
t390 = rSges(4,1) * t447 + rSges(4,2) * t446 + rSges(4,3) * t534;
t389 = Icges(4,1) * t449 + Icges(4,4) * t448 + Icges(4,5) * t533;
t388 = Icges(4,1) * t447 + Icges(4,4) * t446 + Icges(4,5) * t534;
t387 = Icges(4,4) * t449 + Icges(4,2) * t448 + Icges(4,6) * t533;
t386 = Icges(4,4) * t447 + Icges(4,2) * t446 + Icges(4,6) * t534;
t385 = Icges(4,5) * t449 + Icges(4,6) * t448 + Icges(4,3) * t533;
t384 = Icges(4,5) * t447 + Icges(4,6) * t446 + Icges(4,3) * t534;
t383 = qJD(1) * t436 - t461 * t481 + t456;
t382 = -t461 * t524 + (-t435 - t465) * qJD(1);
t381 = (t435 * t492 + t436 * t495) * qJD(2);
t379 = rSges(5,1) * t440 + rSges(5,2) * t439 + rSges(5,3) * t533;
t378 = rSges(5,1) * t438 + rSges(5,2) * t437 + rSges(5,3) * t534;
t377 = Icges(5,1) * t440 + Icges(5,4) * t439 + Icges(5,5) * t533;
t376 = Icges(5,1) * t438 + Icges(5,4) * t437 + Icges(5,5) * t534;
t375 = Icges(5,4) * t440 + Icges(5,2) * t439 + Icges(5,6) * t533;
t374 = Icges(5,4) * t438 + Icges(5,2) * t437 + Icges(5,6) * t534;
t373 = Icges(5,5) * t440 + Icges(5,6) * t439 + Icges(5,3) * t533;
t372 = Icges(5,5) * t438 + Icges(5,6) * t437 + Icges(5,3) * t534;
t371 = -pkin(11) * t494 + t491 * t530;
t370 = rSges(6,1) * t420 + rSges(6,2) * t419 + rSges(6,3) * t533;
t369 = rSges(6,1) * t418 + rSges(6,2) * t417 + rSges(6,3) * t534;
t367 = Icges(6,1) * t420 + Icges(6,4) * t419 + Icges(6,5) * t533;
t366 = Icges(6,1) * t418 + Icges(6,4) * t417 + Icges(6,5) * t534;
t365 = Icges(6,4) * t420 + Icges(6,2) * t419 + Icges(6,6) * t533;
t364 = Icges(6,4) * t418 + Icges(6,2) * t417 + Icges(6,6) * t534;
t363 = Icges(6,5) * t420 + Icges(6,6) * t419 + Icges(6,3) * t533;
t362 = Icges(6,5) * t418 + Icges(6,6) * t417 + Icges(6,3) * t534;
t360 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t533;
t359 = rSges(7,1) * t410 + rSges(7,2) * t409 + rSges(7,3) * t534;
t358 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t533;
t357 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t534;
t356 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t533;
t355 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t534;
t354 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t533;
t353 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t534;
t348 = t492 * t528 + t495 * t504;
t347 = t492 * t504 - t495 * t528;
t346 = t391 * t474 - t434 * t452 + t507;
t345 = -t390 * t474 + t434 * t453 + t505;
t344 = t390 * t452 - t391 * t453 + t529;
t343 = t379 * t455 - t416 * t422 + t501;
t342 = -t378 * t455 + t416 * t423 + t500;
t341 = t378 * t422 - t379 * t423 + t508;
t340 = t370 * t445 - t402 * t407 + t499;
t339 = -t369 * t445 + t403 * t407 + t498;
t338 = t369 * t402 - t370 * t403 + t502;
t337 = t348 * t445 + t360 * t421 - t371 * t402 - t396 * t401 + t499;
t336 = -t347 * t445 - t359 * t421 + t371 * t403 + t397 * t401 + t498;
t335 = t347 * t402 - t348 * t403 + t359 * t396 - t360 * t397 + t502;
t1 = m(4) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + t397 * ((t354 * t534 + t356 * t409 + t358 * t410) * t396 + (t353 * t534 + t409 * t355 + t410 * t357) * t397 + (t398 * t534 + t399 * t409 + t400 * t410) * t421) / 0.2e1 + t421 * ((-t353 * t397 - t354 * t396 - t398 * t421) * t494 + ((-t356 * t472 + t358 * t473) * t396 + (-t355 * t472 + t357 * t473) * t397 + (-t399 * t472 + t400 * t473) * t421) * t491) / 0.2e1 + t396 * ((t354 * t533 + t411 * t356 + t412 * t358) * t396 + (t353 * t533 + t355 * t411 + t357 * t412) * t397 + (t398 * t533 + t399 * t411 + t400 * t412) * t421) / 0.2e1 + t402 * ((t363 * t533 + t419 * t365 + t420 * t367) * t402 + (t362 * t533 + t364 * t419 + t366 * t420) * t403 + (t404 * t533 + t405 * t419 + t406 * t420) * t445) / 0.2e1 + t403 * ((t363 * t534 + t365 * t417 + t367 * t418) * t402 + (t362 * t534 + t417 * t364 + t418 * t366) * t403 + (t404 * t534 + t405 * t417 + t406 * t418) * t445) / 0.2e1 + t445 * ((-t362 * t403 - t363 * t402 - t404 * t445) * t494 + ((-t365 * t477 + t367 * t478) * t402 + (-t364 * t477 + t366 * t478) * t403 + (-t405 * t477 + t406 * t478) * t445) * t491) / 0.2e1 + t422 * ((t373 * t533 + t439 * t375 + t440 * t377) * t422 + (t372 * t533 + t374 * t439 + t376 * t440) * t423 + (t413 * t533 + t414 * t439 + t415 * t440) * t455) / 0.2e1 + t423 * ((t373 * t534 + t375 * t437 + t377 * t438) * t422 + (t372 * t534 + t437 * t374 + t438 * t376) * t423 + (t413 * t534 + t414 * t437 + t415 * t438) * t455) / 0.2e1 + t455 * ((-t372 * t423 - t373 * t422 - t413 * t455) * t494 + ((-t375 * t482 + t377 * t483) * t422 + (-t374 * t482 + t376 * t483) * t423 + (-t414 * t482 + t415 * t483) * t455) * t491) / 0.2e1 + t453 * ((t385 * t534 + t387 * t446 + t389 * t447) * t452 + (t384 * t534 + t446 * t386 + t447 * t388) * t453 + (t424 * t534 + t427 * t446 + t430 * t447) * t474) / 0.2e1 + t474 * ((-t384 * t453 - t385 * t452 - t424 * t474) * t494 + ((-t387 * t490 + t389 * t493) * t452 + (-t386 * t490 + t388 * t493) * t453 + (-t427 * t490 + t430 * t493) * t474) * t491) / 0.2e1 + t452 * ((t385 * t533 + t448 * t387 + t449 * t389) * t452 + (t384 * t533 + t386 * t448 + t388 * t449) * t453 + (t424 * t533 + t427 * t448 + t430 * t449) * t474) / 0.2e1 + qJD(1) * ((t494 * t459 + t491 * t460) * qJD(1) + ((t429 * t494 + t432 * t491) * t492 - (t428 * t494 + t431 * t491) * t495) * qJD(2)) / 0.2e1 + m(7) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(6) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(5) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 - ((-t495 * t458 + t492 * t509) * qJD(1) + (t495 ^ 2 * t425 + (t510 * t492 + (-t426 + t511) * t495) * t492) * qJD(2)) * t524 / 0.2e1 + ((t492 * t458 + t495 * t509) * qJD(1) + (t492 ^ 2 * t426 + (t511 * t495 + (-t425 + t510) * t492) * t495) * qJD(2)) * t481 / 0.2e1 + (m(2) * (t462 ^ 2 + t463 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
