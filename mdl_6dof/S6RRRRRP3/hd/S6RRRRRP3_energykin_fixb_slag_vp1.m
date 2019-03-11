% Calculate kinetic energy for
% S6RRRRRP3
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:49
% EndTime: 2019-03-10 01:05:52
% DurationCPUTime: 2.85s
% Computational Cost: add. (2101->293), mult. (2298->459), div. (0->0), fcn. (2255->10), ass. (0->161)
t555 = Icges(6,1) + Icges(7,1);
t554 = Icges(6,4) + Icges(7,4);
t553 = -Icges(7,5) - Icges(6,5);
t552 = Icges(6,2) + Icges(7,2);
t551 = -Icges(7,6) - Icges(6,6);
t550 = -Icges(7,3) - Icges(6,3);
t472 = cos(qJ(1));
t466 = qJ(2) + qJ(3);
t460 = sin(t466);
t462 = cos(t466);
t470 = cos(qJ(4));
t530 = t470 * pkin(4);
t482 = pkin(10) * t460 + t462 * t530;
t467 = sin(qJ(4));
t469 = sin(qJ(1));
t516 = t467 * t469;
t370 = pkin(4) * t516 + t472 * t482;
t383 = -pkin(10) * t462 + t460 * t530;
t458 = qJD(2) * t469;
t443 = qJD(3) * t469 + t458;
t503 = qJD(4) * t460;
t424 = t472 * t503 + t443;
t454 = -qJD(4) * t462 + qJD(1);
t549 = t454 * t370 - t383 * t424;
t465 = qJ(4) + qJ(5);
t459 = sin(t465);
t461 = cos(t465);
t518 = t462 * t469;
t413 = -t459 * t518 - t461 * t472;
t414 = -t459 * t472 + t461 * t518;
t520 = t460 * t469;
t548 = -t551 * t413 - t553 * t414 - t550 * t520;
t517 = t462 * t472;
t415 = -t459 * t517 + t461 * t469;
t416 = t459 * t469 + t461 * t517;
t519 = t460 * t472;
t547 = -t551 * t415 - t553 * t416 - t550 * t519;
t546 = t552 * t413 + t554 * t414 - t551 * t520;
t545 = t552 * t415 + t554 * t416 - t551 * t519;
t544 = t554 * t413 + t555 * t414 - t553 * t520;
t543 = t554 * t415 + t555 * t416 - t553 * t519;
t515 = t467 * t472;
t369 = -pkin(4) * t515 + t469 * t482;
t444 = (-qJD(2) - qJD(3)) * t472;
t425 = t469 * t503 + t444;
t542 = -t369 * t454 + t425 * t383;
t541 = t550 * t462 + (t551 * t459 - t553 * t461) * t460;
t540 = t551 * t462 + (-t552 * t459 + t554 * t461) * t460;
t539 = t553 * t462 + (-t554 * t459 + t555 * t461) * t460;
t498 = pkin(3) * t462 + pkin(9) * t460;
t429 = t498 * t472;
t439 = pkin(3) * t460 - pkin(9) * t462;
t538 = qJD(1) * t429 - t439 * t443;
t471 = cos(qJ(2));
t531 = pkin(2) * t471;
t468 = sin(qJ(2));
t527 = Icges(3,4) * t468;
t526 = Icges(3,4) * t471;
t525 = Icges(4,4) * t460;
t524 = Icges(4,4) * t462;
t514 = t469 * t470;
t513 = t470 * t472;
t506 = pkin(5) * t461;
t480 = qJ(6) * t460 + t462 * t506;
t500 = pkin(5) * t459;
t512 = rSges(7,1) * t414 + rSges(7,2) * t413 + rSges(7,3) * t520 + t469 * t480 - t472 * t500;
t511 = rSges(7,1) * t416 + rSges(7,2) * t415 + rSges(7,3) * t519 + t469 * t500 + t472 * t480;
t510 = (-qJ(6) - rSges(7,3)) * t462 + (rSges(7,1) * t461 - rSges(7,2) * t459 + t506) * t460;
t402 = -pkin(8) * t472 + t469 * t531;
t403 = pkin(8) * t469 + t472 * t531;
t504 = qJD(2) * t472;
t509 = t402 * t458 + t403 * t504;
t441 = qJD(1) * (pkin(1) * t472 + pkin(7) * t469);
t508 = qJD(1) * t403 + t441;
t453 = pkin(1) * t469 - pkin(7) * t472;
t507 = -t402 - t453;
t502 = qJD(5) * t460;
t501 = pkin(2) * qJD(2) * t468;
t499 = t472 * t501;
t497 = rSges(3,1) * t471 - rSges(3,2) * t468;
t496 = rSges(4,1) * t462 - rSges(4,2) * t460;
t495 = Icges(3,1) * t471 - t527;
t494 = Icges(4,1) * t462 - t525;
t493 = -Icges(3,2) * t468 + t526;
t492 = -Icges(4,2) * t460 + t524;
t491 = Icges(3,5) * t471 - Icges(3,6) * t468;
t490 = Icges(4,5) * t462 - Icges(4,6) * t460;
t420 = -Icges(3,6) * t472 + t469 * t493;
t422 = -Icges(3,5) * t472 + t469 * t495;
t489 = t420 * t468 - t422 * t471;
t421 = Icges(3,6) * t469 + t472 * t493;
t423 = Icges(3,5) * t469 + t472 * t495;
t488 = -t421 * t468 + t423 * t471;
t446 = Icges(3,2) * t471 + t527;
t447 = Icges(3,1) * t468 + t526;
t487 = -t446 * t468 + t447 * t471;
t428 = t498 * t469;
t486 = t443 * t428 - t429 * t444 + t509;
t485 = qJD(6) * t460 - t501;
t484 = t444 * t439 + (-t428 + t507) * qJD(1);
t483 = -t469 * t501 + t508;
t481 = (Icges(4,5) * t460 + Icges(4,6) * t462) * qJD(1) + (-Icges(4,3) * t472 + t469 * t490) * t444 + (Icges(4,3) * t469 + t472 * t490) * t443;
t479 = t424 * t369 - t370 * t425 + t486;
t478 = t483 + t538;
t477 = t484 - t499;
t407 = -Icges(4,6) * t472 + t469 * t492;
t408 = Icges(4,6) * t469 + t472 * t492;
t409 = -Icges(4,5) * t472 + t469 * t494;
t410 = Icges(4,5) * t469 + t472 * t494;
t436 = Icges(4,2) * t462 + t525;
t437 = Icges(4,1) * t460 + t524;
t476 = (-t408 * t460 + t410 * t462) * t443 + (-t407 * t460 + t409 * t462) * t444 + (-t436 * t460 + t437 * t462) * qJD(1);
t450 = rSges(2,1) * t472 - rSges(2,2) * t469;
t449 = rSges(2,1) * t469 + rSges(2,2) * t472;
t448 = rSges(3,1) * t468 + rSges(3,2) * t471;
t445 = Icges(3,5) * t468 + Icges(3,6) * t471;
t438 = rSges(4,1) * t460 + rSges(4,2) * t462;
t434 = qJD(1) + (-qJD(4) - qJD(5)) * t462;
t433 = t462 * t513 + t516;
t432 = -t462 * t515 + t514;
t431 = t462 * t514 - t515;
t430 = -t462 * t516 - t513;
t427 = rSges(3,3) * t469 + t472 * t497;
t426 = -rSges(3,3) * t472 + t469 * t497;
t419 = Icges(3,3) * t469 + t472 * t491;
t418 = -Icges(3,3) * t472 + t469 * t491;
t412 = rSges(4,3) * t469 + t472 * t496;
t411 = -rSges(4,3) * t472 + t469 * t496;
t401 = -rSges(5,3) * t462 + (rSges(5,1) * t470 - rSges(5,2) * t467) * t460;
t400 = -Icges(5,5) * t462 + (Icges(5,1) * t470 - Icges(5,4) * t467) * t460;
t399 = -Icges(5,6) * t462 + (Icges(5,4) * t470 - Icges(5,2) * t467) * t460;
t398 = -Icges(5,3) * t462 + (Icges(5,5) * t470 - Icges(5,6) * t467) * t460;
t394 = -rSges(6,3) * t462 + (rSges(6,1) * t461 - rSges(6,2) * t459) * t460;
t392 = t469 * t502 + t425;
t391 = t472 * t502 + t424;
t381 = qJD(1) * t427 - t448 * t458 + t441;
t380 = -t448 * t504 + (-t426 - t453) * qJD(1);
t379 = rSges(5,1) * t433 + rSges(5,2) * t432 + rSges(5,3) * t519;
t378 = rSges(5,1) * t431 + rSges(5,2) * t430 + rSges(5,3) * t520;
t377 = Icges(5,1) * t433 + Icges(5,4) * t432 + Icges(5,5) * t519;
t376 = Icges(5,1) * t431 + Icges(5,4) * t430 + Icges(5,5) * t520;
t375 = Icges(5,4) * t433 + Icges(5,2) * t432 + Icges(5,6) * t519;
t374 = Icges(5,4) * t431 + Icges(5,2) * t430 + Icges(5,6) * t520;
t373 = Icges(5,5) * t433 + Icges(5,6) * t432 + Icges(5,3) * t519;
t372 = Icges(5,5) * t431 + Icges(5,6) * t430 + Icges(5,3) * t520;
t371 = (t426 * t469 + t427 * t472) * qJD(2);
t367 = rSges(6,1) * t416 + rSges(6,2) * t415 + rSges(6,3) * t519;
t365 = rSges(6,1) * t414 + rSges(6,2) * t413 + rSges(6,3) * t520;
t347 = qJD(1) * t412 - t438 * t443 + t483;
t346 = -t499 + t438 * t444 + (-t411 + t507) * qJD(1);
t345 = t411 * t443 - t412 * t444 + t509;
t344 = t379 * t454 - t401 * t424 + t478;
t343 = -t378 * t454 + t401 * t425 + t477;
t342 = t378 * t424 - t379 * t425 + t486;
t341 = t367 * t434 - t391 * t394 + t478 + t549;
t340 = -t365 * t434 + t392 * t394 + t477 + t542;
t339 = t365 * t391 - t367 * t392 + t479;
t338 = -t391 * t510 + t434 * t511 + t469 * t485 + t508 + t538 + t549;
t337 = t392 * t510 - t434 * t512 + t472 * t485 + t484 + t542;
t336 = -qJD(6) * t462 + t391 * t512 - t392 * t511 + t479;
t1 = ((t469 * t445 + t472 * t487) * qJD(1) + (t469 ^ 2 * t419 + (t489 * t472 + (-t418 + t488) * t469) * t472) * qJD(2)) * t458 / 0.2e1 - ((-t472 * t445 + t469 * t487) * qJD(1) + (t472 ^ 2 * t418 + (t488 * t469 + (-t419 + t489) * t472) * t469) * qJD(2)) * t504 / 0.2e1 + t425 * ((t373 * t520 + t375 * t430 + t377 * t431) * t424 + (t372 * t520 + t430 * t374 + t431 * t376) * t425 + (t398 * t520 + t399 * t430 + t400 * t431) * t454) / 0.2e1 + t454 * ((-t372 * t425 - t373 * t424 - t398 * t454) * t462 + ((-t375 * t467 + t377 * t470) * t424 + (-t374 * t467 + t376 * t470) * t425 + (-t399 * t467 + t400 * t470) * t454) * t460) / 0.2e1 + t424 * ((t373 * t519 + t432 * t375 + t433 * t377) * t424 + (t372 * t519 + t374 * t432 + t376 * t433) * t425 + (t398 * t519 + t399 * t432 + t400 * t433) * t454) / 0.2e1 + t443 * (t481 * t469 + t476 * t472) / 0.2e1 + t444 * (t469 * t476 - t481 * t472) / 0.2e1 + m(7) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(6) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(5) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(4) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(3) * (t371 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + ((t415 * t540 + t416 * t539 + t519 * t541) * t434 + (t415 * t546 + t416 * t544 + t519 * t548) * t392 + (t545 * t415 + t543 * t416 + t547 * t519) * t391) * t391 / 0.2e1 + ((t413 * t540 + t414 * t539 + t520 * t541) * t434 + (t546 * t413 + t544 * t414 + t548 * t520) * t392 + (t413 * t545 + t414 * t543 + t520 * t547) * t391) * t392 / 0.2e1 + ((-t391 * t547 - t392 * t548 - t434 * t541) * t462 + ((-t459 * t540 + t461 * t539) * t434 + (-t459 * t546 + t461 * t544) * t392 + (-t459 * t545 + t461 * t543) * t391) * t460) * t434 / 0.2e1 + (Icges(2,3) + m(2) * (t449 ^ 2 + t450 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t408 * t462 + t410 * t460) * t443 + (t407 * t462 + t409 * t460) * t444 + ((t421 * t471 + t423 * t468) * t469 - (t420 * t471 + t422 * t468) * t472) * qJD(2) + (t462 * t436 + t460 * t437 + t471 * t446 + t468 * t447) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
