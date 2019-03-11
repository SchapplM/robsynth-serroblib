% Calculate kinetic energy for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:17
% EndTime: 2019-03-09 14:02:19
% DurationCPUTime: 2.80s
% Computational Cost: add. (2184->362), mult. (2627->576), div. (0->0), fcn. (2640->12), ass. (0->173)
t487 = cos(pkin(11));
t539 = t487 * pkin(3);
t489 = sin(qJ(2));
t538 = Icges(3,4) * t489;
t491 = cos(qJ(2));
t537 = Icges(3,4) * t491;
t486 = sin(pkin(11));
t490 = sin(qJ(1));
t536 = t486 * t490;
t492 = cos(qJ(1));
t535 = t486 * t492;
t534 = t489 * t490;
t533 = t489 * t492;
t532 = t491 * t490;
t531 = t491 * t492;
t485 = pkin(11) + qJ(4);
t479 = qJ(5) + t485;
t474 = cos(t479);
t529 = pkin(5) * t474;
t473 = sin(t479);
t528 = pkin(5) * t473;
t509 = pkin(2) * t491 + qJ(3) * t489;
t446 = t509 * t490;
t461 = pkin(1) * t490 - pkin(7) * t492;
t527 = -t446 - t461;
t478 = cos(t485);
t526 = pkin(4) * t478;
t480 = qJD(2) * t490;
t521 = qJD(4) * t489;
t448 = t492 * t521 + t480;
t523 = qJD(2) * t492;
t522 = qJD(3) * t489;
t520 = qJD(5) * t489;
t519 = qJD(6) * t489;
t518 = -qJD(4) - qJD(5);
t447 = t509 * t492;
t453 = qJD(1) * (pkin(1) * t492 + pkin(7) * t490);
t517 = qJD(1) * t447 + t490 * t522 + t453;
t422 = t492 * t520 + t448;
t477 = sin(t485);
t514 = pkin(4) * t477;
t457 = pkin(2) * t489 - qJ(3) * t491;
t513 = qJD(2) * (pkin(8) * t491 - t489 * t539 - t457);
t512 = qJD(2) * (rSges(4,3) * t491 - (rSges(4,1) * t487 - rSges(4,2) * t486) * t489 - t457);
t449 = t490 * t521 - t523;
t423 = t490 * t520 + t449;
t511 = -qJD(3) * t491 + t446 * t480 + t447 * t523;
t510 = rSges(3,1) * t491 - rSges(3,2) * t489;
t508 = Icges(3,1) * t491 - t538;
t507 = -Icges(3,2) * t489 + t537;
t506 = Icges(3,5) * t491 - Icges(3,6) * t489;
t431 = -Icges(3,6) * t492 + t490 * t507;
t433 = -Icges(3,5) * t492 + t490 * t508;
t505 = t431 * t489 - t433 * t491;
t432 = Icges(3,6) * t490 + t492 * t507;
t434 = Icges(3,5) * t490 + t492 * t508;
t504 = -t432 * t489 + t434 * t491;
t455 = Icges(3,2) * t491 + t538;
t456 = Icges(3,1) * t489 + t537;
t503 = -t455 * t489 + t456 * t491;
t501 = pkin(8) * t489 + t491 * t539;
t392 = -pkin(3) * t535 + t490 * t501;
t393 = pkin(3) * t536 + t492 * t501;
t502 = t392 * t480 + t393 * t523 + t511;
t500 = pkin(10) * t489 + t491 * t529;
t499 = pkin(9) * t489 + t491 * t526;
t498 = qJD(1) * t393 + t490 * t513 + t517;
t350 = t490 * t499 - t492 * t514;
t351 = t490 * t514 + t492 * t499;
t497 = t448 * t350 - t351 * t449 + t502;
t469 = t492 * t522;
t496 = t469 + (-t392 + t527) * qJD(1) + t492 * t513;
t394 = -pkin(9) * t491 + t489 * t526;
t470 = -qJD(4) * t491 + qJD(1);
t495 = t470 * t351 - t394 * t448 + t498;
t494 = -t350 * t470 + t449 * t394 + t496;
t476 = qJ(6) + t479;
t463 = cos(t476);
t462 = sin(t476);
t460 = rSges(2,1) * t492 - rSges(2,2) * t490;
t459 = rSges(2,1) * t490 + rSges(2,2) * t492;
t458 = rSges(3,1) * t489 + rSges(3,2) * t491;
t454 = Icges(3,5) * t489 + Icges(3,6) * t491;
t452 = t491 * t518 + qJD(1);
t445 = t487 * t531 + t536;
t444 = -t486 * t531 + t487 * t490;
t443 = t487 * t532 - t535;
t442 = -t486 * t532 - t487 * t492;
t441 = qJD(1) + (-qJD(6) + t518) * t491;
t439 = rSges(3,3) * t490 + t492 * t510;
t438 = -rSges(3,3) * t492 + t490 * t510;
t430 = Icges(3,3) * t490 + t492 * t506;
t429 = -Icges(3,3) * t492 + t490 * t506;
t428 = t477 * t490 + t478 * t531;
t427 = -t477 * t531 + t478 * t490;
t426 = -t477 * t492 + t478 * t532;
t425 = -t477 * t532 - t478 * t492;
t421 = -Icges(4,5) * t491 + (Icges(4,1) * t487 - Icges(4,4) * t486) * t489;
t420 = -Icges(4,6) * t491 + (Icges(4,4) * t487 - Icges(4,2) * t486) * t489;
t419 = -Icges(4,3) * t491 + (Icges(4,5) * t487 - Icges(4,6) * t486) * t489;
t417 = t473 * t490 + t474 * t531;
t416 = -t473 * t531 + t474 * t490;
t415 = -t473 * t492 + t474 * t532;
t414 = -t473 * t532 - t474 * t492;
t413 = -rSges(5,3) * t491 + (rSges(5,1) * t478 - rSges(5,2) * t477) * t489;
t412 = -Icges(5,5) * t491 + (Icges(5,1) * t478 - Icges(5,4) * t477) * t489;
t411 = -Icges(5,6) * t491 + (Icges(5,4) * t478 - Icges(5,2) * t477) * t489;
t410 = -Icges(5,3) * t491 + (Icges(5,5) * t478 - Icges(5,6) * t477) * t489;
t409 = t462 * t490 + t463 * t531;
t408 = -t462 * t531 + t463 * t490;
t407 = -t462 * t492 + t463 * t532;
t406 = -t462 * t532 - t463 * t492;
t404 = -rSges(6,3) * t491 + (rSges(6,1) * t474 - rSges(6,2) * t473) * t489;
t403 = -Icges(6,5) * t491 + (Icges(6,1) * t474 - Icges(6,4) * t473) * t489;
t402 = -Icges(6,6) * t491 + (Icges(6,4) * t474 - Icges(6,2) * t473) * t489;
t401 = -Icges(6,3) * t491 + (Icges(6,5) * t474 - Icges(6,6) * t473) * t489;
t400 = t490 * t519 + t423;
t399 = t492 * t519 + t422;
t398 = -rSges(7,3) * t491 + (rSges(7,1) * t463 - rSges(7,2) * t462) * t489;
t397 = -Icges(7,5) * t491 + (Icges(7,1) * t463 - Icges(7,4) * t462) * t489;
t396 = -Icges(7,6) * t491 + (Icges(7,4) * t463 - Icges(7,2) * t462) * t489;
t395 = -Icges(7,3) * t491 + (Icges(7,5) * t463 - Icges(7,6) * t462) * t489;
t391 = rSges(4,1) * t445 + rSges(4,2) * t444 + rSges(4,3) * t533;
t390 = rSges(4,1) * t443 + rSges(4,2) * t442 + rSges(4,3) * t534;
t389 = Icges(4,1) * t445 + Icges(4,4) * t444 + Icges(4,5) * t533;
t388 = Icges(4,1) * t443 + Icges(4,4) * t442 + Icges(4,5) * t534;
t387 = Icges(4,4) * t445 + Icges(4,2) * t444 + Icges(4,6) * t533;
t386 = Icges(4,4) * t443 + Icges(4,2) * t442 + Icges(4,6) * t534;
t385 = Icges(4,5) * t445 + Icges(4,6) * t444 + Icges(4,3) * t533;
t384 = Icges(4,5) * t443 + Icges(4,6) * t442 + Icges(4,3) * t534;
t382 = qJD(1) * t439 - t458 * t480 + t453;
t381 = -t458 * t523 + (-t438 - t461) * qJD(1);
t378 = (t438 * t490 + t439 * t492) * qJD(2);
t377 = rSges(5,1) * t428 + rSges(5,2) * t427 + rSges(5,3) * t533;
t376 = rSges(5,1) * t426 + rSges(5,2) * t425 + rSges(5,3) * t534;
t375 = Icges(5,1) * t428 + Icges(5,4) * t427 + Icges(5,5) * t533;
t374 = Icges(5,1) * t426 + Icges(5,4) * t425 + Icges(5,5) * t534;
t373 = Icges(5,4) * t428 + Icges(5,2) * t427 + Icges(5,6) * t533;
t372 = Icges(5,4) * t426 + Icges(5,2) * t425 + Icges(5,6) * t534;
t371 = Icges(5,5) * t428 + Icges(5,6) * t427 + Icges(5,3) * t533;
t370 = Icges(5,5) * t426 + Icges(5,6) * t425 + Icges(5,3) * t534;
t368 = -pkin(10) * t491 + t489 * t529;
t367 = rSges(6,1) * t417 + rSges(6,2) * t416 + rSges(6,3) * t533;
t366 = rSges(6,1) * t415 + rSges(6,2) * t414 + rSges(6,3) * t534;
t365 = Icges(6,1) * t417 + Icges(6,4) * t416 + Icges(6,5) * t533;
t364 = Icges(6,1) * t415 + Icges(6,4) * t414 + Icges(6,5) * t534;
t363 = Icges(6,4) * t417 + Icges(6,2) * t416 + Icges(6,6) * t533;
t362 = Icges(6,4) * t415 + Icges(6,2) * t414 + Icges(6,6) * t534;
t361 = Icges(6,5) * t417 + Icges(6,6) * t416 + Icges(6,3) * t533;
t360 = Icges(6,5) * t415 + Icges(6,6) * t414 + Icges(6,3) * t534;
t359 = rSges(7,1) * t409 + rSges(7,2) * t408 + rSges(7,3) * t533;
t358 = rSges(7,1) * t407 + rSges(7,2) * t406 + rSges(7,3) * t534;
t357 = Icges(7,1) * t409 + Icges(7,4) * t408 + Icges(7,5) * t533;
t356 = Icges(7,1) * t407 + Icges(7,4) * t406 + Icges(7,5) * t534;
t355 = Icges(7,4) * t409 + Icges(7,2) * t408 + Icges(7,6) * t533;
t354 = Icges(7,4) * t407 + Icges(7,2) * t406 + Icges(7,6) * t534;
t353 = Icges(7,5) * t409 + Icges(7,6) * t408 + Icges(7,3) * t533;
t352 = Icges(7,5) * t407 + Icges(7,6) * t406 + Icges(7,3) * t534;
t347 = t490 * t528 + t492 * t500;
t346 = t490 * t500 - t492 * t528;
t345 = qJD(1) * t391 + t490 * t512 + t517;
t344 = t469 + t492 * t512 + (-t390 + t527) * qJD(1);
t343 = (t390 * t490 + t391 * t492) * qJD(2) + t511;
t342 = t377 * t470 - t413 * t448 + t498;
t341 = -t376 * t470 + t413 * t449 + t496;
t340 = t376 * t448 - t377 * t449 + t502;
t339 = t367 * t452 - t404 * t422 + t495;
t338 = -t366 * t452 + t404 * t423 + t494;
t337 = t366 * t422 - t367 * t423 + t497;
t336 = t347 * t452 + t359 * t441 - t368 * t422 - t398 * t399 + t495;
t335 = -t346 * t452 - t358 * t441 + t368 * t423 + t398 * t400 + t494;
t334 = t346 * t422 - t347 * t423 + t358 * t399 - t359 * t400 + t497;
t1 = m(6) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + t422 * ((t361 * t533 + t416 * t363 + t417 * t365) * t422 + (t360 * t533 + t362 * t416 + t364 * t417) * t423 + (t401 * t533 + t402 * t416 + t403 * t417) * t452) / 0.2e1 + t448 * ((t371 * t533 + t427 * t373 + t428 * t375) * t448 + (t370 * t533 + t372 * t427 + t374 * t428) * t449 + (t410 * t533 + t411 * t427 + t412 * t428) * t470) / 0.2e1 + t449 * ((t371 * t534 + t373 * t425 + t375 * t426) * t448 + (t370 * t534 + t425 * t372 + t426 * t374) * t449 + (t410 * t534 + t411 * t425 + t412 * t426) * t470) / 0.2e1 + t470 * ((-t370 * t449 - t371 * t448 - t410 * t470) * t491 + ((-t373 * t477 + t375 * t478) * t448 + (-t372 * t477 + t374 * t478) * t449 + (-t411 * t477 + t412 * t478) * t470) * t489) / 0.2e1 + m(7) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(5) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(4) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(3) * (t378 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + t423 * ((t361 * t534 + t363 * t414 + t365 * t415) * t422 + (t360 * t534 + t414 * t362 + t415 * t364) * t423 + (t401 * t534 + t402 * t414 + t403 * t415) * t452) / 0.2e1 + t400 * ((t353 * t534 + t355 * t406 + t357 * t407) * t399 + (t352 * t534 + t406 * t354 + t407 * t356) * t400 + (t395 * t534 + t396 * t406 + t397 * t407) * t441) / 0.2e1 + t441 * ((-t352 * t400 - t353 * t399 - t395 * t441) * t491 + ((-t355 * t462 + t357 * t463) * t399 + (-t354 * t462 + t356 * t463) * t400 + (-t396 * t462 + t397 * t463) * t441) * t489) / 0.2e1 + t399 * ((t353 * t533 + t408 * t355 + t409 * t357) * t399 + (t352 * t533 + t354 * t408 + t356 * t409) * t400 + (t395 * t533 + t396 * t408 + t397 * t409) * t441) / 0.2e1 + t452 * ((-t360 * t423 - t361 * t422 - t401 * t452) * t491 + ((-t363 * t473 + t365 * t474) * t422 + (-t362 * t473 + t364 * t474) * t423 + (-t402 * t473 + t403 * t474) * t452) * t489) / 0.2e1 + (m(2) * (t459 ^ 2 + t460 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t384 * t492 - t385 * t490) * t491 + ((-t387 * t486 + t389 * t487) * t490 - (-t386 * t486 + t388 * t487) * t492) * t489 + (t432 * t491 + t434 * t489) * t490 - (t431 * t491 + t433 * t489) * t492) * qJD(2) + ((-t419 + t455) * t491 + (-t420 * t486 + t421 * t487 + t456) * t489) * qJD(1)) * qJD(1) / 0.2e1 + (((-t384 * t533 - t386 * t444 - t388 * t445 + t505 * t492) * t492 + ((-t429 + t504) * t492 + t385 * t533 + t387 * t444 + t389 * t445 + t430 * t490) * t490) * qJD(2) + (t419 * t533 + t420 * t444 + t421 * t445 + t490 * t454 + t492 * t503) * qJD(1)) * t480 / 0.2e1 - (((-t384 * t534 - t386 * t442 - t388 * t443 + t429 * t492) * t492 + (t385 * t534 + t387 * t442 + t389 * t443 + (-t430 + t505) * t492 + t504 * t490) * t490) * qJD(2) + (t419 * t534 + t420 * t442 + t421 * t443 - t492 * t454 + t490 * t503) * qJD(1)) * t523 / 0.2e1;
T  = t1;
