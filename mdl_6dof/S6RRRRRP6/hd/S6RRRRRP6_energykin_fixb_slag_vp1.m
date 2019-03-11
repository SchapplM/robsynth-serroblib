% Calculate kinetic energy for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:31
% EndTime: 2019-03-10 01:26:35
% DurationCPUTime: 3.28s
% Computational Cost: add. (2053->309), mult. (2709->492), div. (0->0), fcn. (2732->10), ass. (0->157)
t544 = Icges(6,1) + Icges(7,1);
t543 = -Icges(6,4) + Icges(7,5);
t542 = Icges(7,4) + Icges(6,5);
t541 = Icges(6,2) + Icges(7,3);
t540 = -Icges(7,6) + Icges(6,6);
t539 = -Icges(6,3) - Icges(7,2);
t538 = rSges(7,1) + pkin(5);
t537 = rSges(7,3) + qJ(6);
t474 = qJ(3) + qJ(4);
t471 = qJ(5) + t474;
t465 = sin(t471);
t466 = cos(t471);
t480 = cos(qJ(1));
t477 = sin(qJ(1));
t479 = cos(qJ(2));
t514 = t479 * t477;
t412 = t465 * t514 + t466 * t480;
t413 = -t465 * t480 + t466 * t514;
t476 = sin(qJ(2));
t516 = t476 * t477;
t536 = t541 * t412 + t543 * t413 - t540 * t516;
t513 = t479 * t480;
t414 = t465 * t513 - t477 * t466;
t415 = t465 * t477 + t466 * t513;
t515 = t476 * t480;
t535 = t541 * t414 + t543 * t415 - t540 * t515;
t534 = -t540 * t412 + t542 * t413 - t539 * t516;
t533 = -t540 * t414 + t542 * t415 - t539 * t515;
t532 = t543 * t412 + t544 * t413 + t542 * t516;
t531 = t543 * t414 + t544 * t415 + t542 * t515;
t530 = t540 * t479 + (t541 * t465 + t543 * t466) * t476;
t529 = t539 * t479 + (-t540 * t465 + t542 * t466) * t476;
t528 = -t542 * t479 + (t543 * t465 + t544 * t466) * t476;
t478 = cos(qJ(3));
t522 = t478 * pkin(3);
t520 = Icges(3,4) * t476;
t519 = Icges(3,4) * t479;
t475 = sin(qJ(3));
t518 = t475 * t477;
t517 = t475 * t480;
t512 = rSges(7,2) * t516 + t537 * t412 + t538 * t413;
t511 = rSges(7,2) * t515 + t537 * t414 + t538 * t415;
t510 = -rSges(7,2) * t479 + (t537 * t465 + t538 * t466) * t476;
t500 = pkin(2) * t479 + pkin(8) * t476;
t443 = t500 * t477;
t444 = t500 * t480;
t468 = qJD(2) * t477;
t506 = qJD(2) * t480;
t509 = t443 * t468 + t444 * t506;
t470 = cos(t474);
t508 = pkin(4) * t470;
t505 = qJD(3) * t476;
t445 = t480 * t505 + t468;
t504 = qJD(4) * t476;
t503 = qJD(5) * t476;
t502 = -qJD(3) - qJD(4);
t417 = t480 * t504 + t445;
t469 = sin(t474);
t501 = pkin(4) * t469;
t446 = t477 * t505 - t506;
t418 = t477 * t504 + t446;
t499 = rSges(3,1) * t479 - rSges(3,2) * t476;
t498 = Icges(3,1) * t479 - t520;
t497 = -Icges(3,2) * t476 + t519;
t496 = Icges(3,5) * t479 - Icges(3,6) * t476;
t423 = -Icges(3,6) * t480 + t477 * t497;
t426 = -Icges(3,5) * t480 + t477 * t498;
t495 = t423 * t476 - t426 * t479;
t424 = Icges(3,6) * t477 + t480 * t497;
t427 = Icges(3,5) * t477 + t480 * t498;
t494 = -t424 * t476 + t427 * t479;
t452 = Icges(3,2) * t479 + t520;
t453 = Icges(3,1) * t476 + t519;
t493 = -t452 * t476 + t453 * t479;
t490 = pkin(9) * t476 + t479 * t522;
t393 = -pkin(3) * t517 + t477 * t490;
t394 = pkin(3) * t518 + t480 * t490;
t492 = t445 * t393 - t394 * t446 + t509;
t449 = qJD(1) * (pkin(1) * t480 + pkin(7) * t477);
t457 = pkin(2) * t476 - pkin(8) * t479;
t491 = qJD(1) * t444 - t457 * t468 + t449;
t458 = pkin(1) * t477 - pkin(7) * t480;
t489 = (-t443 - t458) * qJD(1) - t457 * t506;
t488 = pkin(10) * t476 + t479 * t508;
t351 = t477 * t488 - t480 * t501;
t352 = t477 * t501 + t480 * t488;
t487 = t417 * t351 - t352 * t418 + t492;
t407 = -pkin(9) * t479 + t476 * t522;
t464 = -qJD(3) * t479 + qJD(1);
t486 = t464 * t394 - t407 * t445 + t491;
t485 = -t393 * t464 + t446 * t407 + t489;
t395 = -pkin(10) * t479 + t476 * t508;
t448 = t479 * t502 + qJD(1);
t484 = t448 * t352 - t395 * t417 + t486;
t483 = -t351 * t448 + t418 * t395 + t485;
t456 = rSges(2,1) * t480 - rSges(2,2) * t477;
t455 = rSges(2,1) * t477 + rSges(2,2) * t480;
t454 = rSges(3,1) * t476 + rSges(3,2) * t479;
t451 = Icges(3,5) * t476 + Icges(3,6) * t479;
t442 = t478 * t513 + t518;
t441 = -t475 * t513 + t477 * t478;
t440 = t478 * t514 - t517;
t439 = -t475 * t514 - t478 * t480;
t438 = qJD(1) + (-qJD(5) + t502) * t479;
t434 = t469 * t477 + t470 * t513;
t433 = -t469 * t513 + t470 * t477;
t432 = -t469 * t480 + t470 * t514;
t431 = -t469 * t514 - t470 * t480;
t430 = rSges(3,3) * t477 + t480 * t499;
t429 = -rSges(3,3) * t480 + t477 * t499;
t428 = -rSges(4,3) * t479 + (rSges(4,1) * t478 - rSges(4,2) * t475) * t476;
t425 = -Icges(4,5) * t479 + (Icges(4,1) * t478 - Icges(4,4) * t475) * t476;
t422 = -Icges(4,6) * t479 + (Icges(4,4) * t478 - Icges(4,2) * t475) * t476;
t421 = Icges(3,3) * t477 + t480 * t496;
t420 = -Icges(3,3) * t480 + t477 * t496;
t419 = -Icges(4,3) * t479 + (Icges(4,5) * t478 - Icges(4,6) * t475) * t476;
t411 = -rSges(5,3) * t479 + (rSges(5,1) * t470 - rSges(5,2) * t469) * t476;
t410 = -Icges(5,5) * t479 + (Icges(5,1) * t470 - Icges(5,4) * t469) * t476;
t409 = -Icges(5,6) * t479 + (Icges(5,4) * t470 - Icges(5,2) * t469) * t476;
t408 = -Icges(5,3) * t479 + (Icges(5,5) * t470 - Icges(5,6) * t469) * t476;
t406 = -rSges(6,3) * t479 + (rSges(6,1) * t466 - rSges(6,2) * t465) * t476;
t398 = t477 * t503 + t418;
t397 = t480 * t503 + t417;
t392 = rSges(4,1) * t442 + rSges(4,2) * t441 + rSges(4,3) * t515;
t391 = rSges(4,1) * t440 + rSges(4,2) * t439 + rSges(4,3) * t516;
t390 = Icges(4,1) * t442 + Icges(4,4) * t441 + Icges(4,5) * t515;
t389 = Icges(4,1) * t440 + Icges(4,4) * t439 + Icges(4,5) * t516;
t388 = Icges(4,4) * t442 + Icges(4,2) * t441 + Icges(4,6) * t515;
t387 = Icges(4,4) * t440 + Icges(4,2) * t439 + Icges(4,6) * t516;
t386 = Icges(4,5) * t442 + Icges(4,6) * t441 + Icges(4,3) * t515;
t385 = Icges(4,5) * t440 + Icges(4,6) * t439 + Icges(4,3) * t516;
t382 = qJD(1) * t430 - t454 * t468 + t449;
t381 = -t454 * t506 + (-t429 - t458) * qJD(1);
t380 = (t429 * t477 + t430 * t480) * qJD(2);
t378 = rSges(5,1) * t434 + rSges(5,2) * t433 + rSges(5,3) * t515;
t377 = rSges(5,1) * t432 + rSges(5,2) * t431 + rSges(5,3) * t516;
t376 = Icges(5,1) * t434 + Icges(5,4) * t433 + Icges(5,5) * t515;
t375 = Icges(5,1) * t432 + Icges(5,4) * t431 + Icges(5,5) * t516;
t374 = Icges(5,4) * t434 + Icges(5,2) * t433 + Icges(5,6) * t515;
t373 = Icges(5,4) * t432 + Icges(5,2) * t431 + Icges(5,6) * t516;
t372 = Icges(5,5) * t434 + Icges(5,6) * t433 + Icges(5,3) * t515;
t371 = Icges(5,5) * t432 + Icges(5,6) * t431 + Icges(5,3) * t516;
t370 = rSges(6,1) * t415 - rSges(6,2) * t414 + rSges(6,3) * t515;
t368 = rSges(6,1) * t413 - rSges(6,2) * t412 + rSges(6,3) * t516;
t348 = t392 * t464 - t428 * t445 + t491;
t347 = -t391 * t464 + t428 * t446 + t489;
t346 = t391 * t445 - t392 * t446 + t509;
t345 = t378 * t448 - t411 * t417 + t486;
t344 = -t377 * t448 + t411 * t418 + t485;
t343 = t377 * t417 - t378 * t418 + t492;
t342 = t370 * t438 - t397 * t406 + t484;
t341 = -t368 * t438 + t398 * t406 + t483;
t340 = t368 * t397 - t370 * t398 + t487;
t339 = qJD(6) * t412 - t397 * t510 + t438 * t511 + t484;
t338 = qJD(6) * t414 + t398 * t510 - t438 * t512 + t483;
t337 = qJD(6) * t465 * t476 + t397 * t512 - t398 * t511 + t487;
t1 = t418 * ((t372 * t516 + t374 * t431 + t376 * t432) * t417 + (t371 * t516 + t431 * t373 + t432 * t375) * t418 + (t408 * t516 + t409 * t431 + t410 * t432) * t448) / 0.2e1 + t448 * ((-t371 * t418 - t372 * t417 - t408 * t448) * t479 + ((-t374 * t469 + t376 * t470) * t417 + (-t373 * t469 + t375 * t470) * t418 + (-t409 * t469 + t410 * t470) * t448) * t476) / 0.2e1 + t417 * ((t372 * t515 + t433 * t374 + t434 * t376) * t417 + (t371 * t515 + t373 * t433 + t375 * t434) * t418 + (t408 * t515 + t409 * t433 + t410 * t434) * t448) / 0.2e1 + t445 * ((t386 * t515 + t441 * t388 + t442 * t390) * t445 + (t385 * t515 + t387 * t441 + t389 * t442) * t446 + (t419 * t515 + t422 * t441 + t425 * t442) * t464) / 0.2e1 + t446 * ((t386 * t516 + t388 * t439 + t440 * t390) * t445 + (t385 * t516 + t439 * t387 + t440 * t389) * t446 + (t419 * t516 + t422 * t439 + t425 * t440) * t464) / 0.2e1 + t464 * ((-t385 * t446 - t386 * t445 - t419 * t464) * t479 + ((-t388 * t475 + t390 * t478) * t445 + (-t387 * t475 + t389 * t478) * t446 + (-t422 * t475 + t425 * t478) * t464) * t476) / 0.2e1 + qJD(1) * ((t452 * t479 + t453 * t476) * qJD(1) + ((t424 * t479 + t427 * t476) * t477 - (t423 * t479 + t426 * t476) * t480) * qJD(2)) / 0.2e1 + m(7) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(6) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(5) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(4) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(3) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 - ((-t480 * t451 + t477 * t493) * qJD(1) + (t480 ^ 2 * t420 + (t494 * t477 + (-t421 + t495) * t480) * t477) * qJD(2)) * t506 / 0.2e1 + ((t477 * t451 + t480 * t493) * qJD(1) + (t477 ^ 2 * t421 + (t495 * t480 + (-t420 + t494) * t477) * t480) * qJD(2)) * t468 / 0.2e1 + ((t530 * t414 + t528 * t415 + t529 * t515) * t438 + (t536 * t414 + t532 * t415 + t534 * t515) * t398 + (t535 * t414 + t531 * t415 + t533 * t515) * t397) * t397 / 0.2e1 + ((t530 * t412 + t528 * t413 + t529 * t516) * t438 + (t536 * t412 + t532 * t413 + t534 * t516) * t398 + (t535 * t412 + t531 * t413 + t533 * t516) * t397) * t398 / 0.2e1 + ((-t533 * t397 - t534 * t398 - t529 * t438) * t479 + ((t530 * t465 + t528 * t466) * t438 + (t536 * t465 + t532 * t466) * t398 + (t535 * t465 + t531 * t466) * t397) * t476) * t438 / 0.2e1 + (m(2) * (t455 ^ 2 + t456 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
