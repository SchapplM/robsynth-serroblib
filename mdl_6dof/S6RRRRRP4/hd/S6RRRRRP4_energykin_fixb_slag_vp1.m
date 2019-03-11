% Calculate kinetic energy for
% S6RRRRRP4
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:11:58
% EndTime: 2019-03-10 01:12:01
% DurationCPUTime: 2.86s
% Computational Cost: add. (2054->284), mult. (2265->450), div. (0->0), fcn. (2232->10), ass. (0->157)
t551 = Icges(6,1) + Icges(7,1);
t550 = -Icges(6,4) + Icges(7,5);
t549 = Icges(7,4) + Icges(6,5);
t548 = Icges(6,2) + Icges(7,3);
t547 = -Icges(7,6) + Icges(6,6);
t546 = -Icges(6,3) - Icges(7,2);
t545 = rSges(7,1) + pkin(5);
t544 = rSges(7,3) + qJ(6);
t470 = qJ(4) + qJ(5);
t468 = cos(t470);
t471 = qJ(2) + qJ(3);
t469 = cos(t471);
t477 = cos(qJ(1));
t466 = sin(t470);
t474 = sin(qJ(1));
t521 = t466 * t474;
t421 = t468 * t477 + t469 * t521;
t515 = t474 * t468;
t422 = -t466 * t477 + t469 * t515;
t467 = sin(t471);
t520 = t467 * t474;
t543 = t548 * t421 + t550 * t422 - t547 * t520;
t518 = t469 * t477;
t423 = t466 * t518 - t515;
t424 = t468 * t518 + t521;
t519 = t467 * t477;
t542 = t548 * t423 + t550 * t424 - t547 * t519;
t541 = -t547 * t421 + t549 * t422 - t546 * t520;
t540 = -t547 * t423 + t549 * t424 - t546 * t519;
t539 = t550 * t421 + t551 * t422 + t549 * t520;
t538 = t550 * t423 + t551 * t424 + t549 * t519;
t537 = t547 * t469 + (t548 * t466 + t550 * t468) * t467;
t536 = t546 * t469 + (-t547 * t466 + t549 * t468) * t467;
t535 = -t549 * t469 + (t550 * t466 + t551 * t468) * t467;
t476 = cos(qJ(2));
t529 = pkin(2) * t476;
t475 = cos(qJ(4));
t528 = pkin(4) * t475;
t473 = sin(qJ(2));
t525 = Icges(3,4) * t473;
t524 = Icges(3,4) * t476;
t523 = Icges(4,4) * t467;
t522 = Icges(4,4) * t469;
t472 = sin(qJ(4));
t517 = t472 * t474;
t516 = t472 * t477;
t514 = t474 * t475;
t513 = t475 * t477;
t512 = rSges(7,2) * t520 + t544 * t421 + t545 * t422;
t511 = rSges(7,2) * t519 + t544 * t423 + t545 * t424;
t510 = -rSges(7,2) * t469 + (t544 * t466 + t545 * t468) * t467;
t409 = -pkin(8) * t477 + t474 * t529;
t410 = pkin(8) * t474 + t477 * t529;
t465 = qJD(2) * t474;
t507 = qJD(2) * t477;
t509 = t409 * t465 + t410 * t507;
t459 = pkin(1) * t474 - pkin(7) * t477;
t508 = -t409 - t459;
t449 = qJD(3) * t474 + t465;
t506 = qJD(4) * t467;
t505 = qJD(5) * t467;
t504 = pkin(2) * qJD(2) * t473;
t432 = t477 * t506 + t449;
t450 = (-qJD(2) - qJD(3)) * t477;
t503 = t477 * t504;
t502 = pkin(3) * t469 + pkin(9) * t467;
t501 = rSges(3,1) * t476 - rSges(3,2) * t473;
t500 = rSges(4,1) * t469 - rSges(4,2) * t467;
t499 = Icges(3,1) * t476 - t525;
t498 = Icges(4,1) * t469 - t523;
t497 = -Icges(3,2) * t473 + t524;
t496 = -Icges(4,2) * t467 + t522;
t495 = Icges(3,5) * t476 - Icges(3,6) * t473;
t494 = Icges(4,5) * t469 - Icges(4,6) * t467;
t428 = -Icges(3,6) * t477 + t474 * t497;
t430 = -Icges(3,5) * t477 + t474 * t499;
t493 = t428 * t473 - t430 * t476;
t429 = Icges(3,6) * t474 + t477 * t497;
t431 = Icges(3,5) * t474 + t477 * t499;
t492 = -t429 * t473 + t431 * t476;
t452 = Icges(3,2) * t476 + t525;
t453 = Icges(3,1) * t473 + t524;
t491 = -t452 * t473 + t453 * t476;
t433 = t474 * t506 + t450;
t436 = t502 * t474;
t437 = t502 * t477;
t490 = t449 * t436 - t437 * t450 + t509;
t448 = qJD(1) * (pkin(1) * t477 + pkin(7) * t474);
t489 = qJD(1) * t410 - t474 * t504 + t448;
t488 = pkin(10) * t467 + t469 * t528;
t487 = qJD(1) * (Icges(4,5) * t467 + Icges(4,6) * t469) + (-Icges(4,3) * t477 + t474 * t494) * t450 + (Icges(4,3) * t474 + t477 * t494) * t449;
t375 = -pkin(4) * t516 + t474 * t488;
t376 = pkin(4) * t517 + t477 * t488;
t486 = t432 * t375 - t376 * t433 + t490;
t447 = pkin(3) * t467 - pkin(9) * t469;
t485 = qJD(1) * t437 - t447 * t449 + t489;
t484 = t450 * t447 + (-t436 + t508) * qJD(1) - t503;
t390 = -pkin(10) * t469 + t467 * t528;
t460 = -qJD(4) * t469 + qJD(1);
t483 = t460 * t376 - t390 * t432 + t485;
t482 = -t375 * t460 + t433 * t390 + t484;
t414 = -Icges(4,6) * t477 + t474 * t496;
t415 = Icges(4,6) * t474 + t477 * t496;
t416 = -Icges(4,5) * t477 + t474 * t498;
t417 = Icges(4,5) * t474 + t477 * t498;
t444 = Icges(4,2) * t469 + t523;
t445 = Icges(4,1) * t467 + t522;
t481 = (-t415 * t467 + t417 * t469) * t449 + (-t414 * t467 + t416 * t469) * t450 + (-t444 * t467 + t445 * t469) * qJD(1);
t456 = rSges(2,1) * t477 - rSges(2,2) * t474;
t455 = rSges(2,1) * t474 + rSges(2,2) * t477;
t454 = rSges(3,1) * t473 + rSges(3,2) * t476;
t451 = Icges(3,5) * t473 + Icges(3,6) * t476;
t446 = rSges(4,1) * t467 + rSges(4,2) * t469;
t442 = qJD(1) + (-qJD(4) - qJD(5)) * t469;
t441 = t469 * t513 + t517;
t440 = -t469 * t516 + t514;
t439 = t469 * t514 - t516;
t438 = -t469 * t517 - t513;
t435 = rSges(3,3) * t474 + t477 * t501;
t434 = -rSges(3,3) * t477 + t474 * t501;
t427 = Icges(3,3) * t474 + t477 * t495;
t426 = -Icges(3,3) * t477 + t474 * t495;
t419 = rSges(4,3) * t474 + t477 * t500;
t418 = -rSges(4,3) * t477 + t474 * t500;
t408 = -rSges(5,3) * t469 + (rSges(5,1) * t475 - rSges(5,2) * t472) * t467;
t407 = -Icges(5,5) * t469 + (Icges(5,1) * t475 - Icges(5,4) * t472) * t467;
t406 = -Icges(5,6) * t469 + (Icges(5,4) * t475 - Icges(5,2) * t472) * t467;
t405 = -Icges(5,3) * t469 + (Icges(5,5) * t475 - Icges(5,6) * t472) * t467;
t401 = -rSges(6,3) * t469 + (rSges(6,1) * t468 - rSges(6,2) * t466) * t467;
t399 = t474 * t505 + t433;
t398 = t477 * t505 + t432;
t387 = qJD(1) * t435 - t454 * t465 + t448;
t386 = -t454 * t507 + (-t434 - t459) * qJD(1);
t385 = rSges(5,1) * t441 + rSges(5,2) * t440 + rSges(5,3) * t519;
t384 = rSges(5,1) * t439 + rSges(5,2) * t438 + rSges(5,3) * t520;
t383 = Icges(5,1) * t441 + Icges(5,4) * t440 + Icges(5,5) * t519;
t382 = Icges(5,1) * t439 + Icges(5,4) * t438 + Icges(5,5) * t520;
t381 = Icges(5,4) * t441 + Icges(5,2) * t440 + Icges(5,6) * t519;
t380 = Icges(5,4) * t439 + Icges(5,2) * t438 + Icges(5,6) * t520;
t379 = Icges(5,5) * t441 + Icges(5,6) * t440 + Icges(5,3) * t519;
t378 = Icges(5,5) * t439 + Icges(5,6) * t438 + Icges(5,3) * t520;
t377 = (t434 * t474 + t435 * t477) * qJD(2);
t373 = rSges(6,1) * t424 - rSges(6,2) * t423 + rSges(6,3) * t519;
t371 = rSges(6,1) * t422 - rSges(6,2) * t421 + rSges(6,3) * t520;
t355 = qJD(1) * t419 - t446 * t449 + t489;
t354 = -t503 + t446 * t450 + (-t418 + t508) * qJD(1);
t353 = t418 * t449 - t419 * t450 + t509;
t352 = t385 * t460 - t408 * t432 + t485;
t351 = -t384 * t460 + t408 * t433 + t484;
t350 = t384 * t432 - t385 * t433 + t490;
t349 = t373 * t442 - t398 * t401 + t483;
t348 = -t371 * t442 + t399 * t401 + t482;
t347 = t371 * t398 - t373 * t399 + t486;
t346 = qJD(6) * t421 - t398 * t510 + t442 * t511 + t483;
t345 = qJD(6) * t423 + t399 * t510 - t442 * t512 + t482;
t344 = qJD(6) * t466 * t467 + t398 * t512 - t399 * t511 + t486;
t1 = m(3) * (t377 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(4) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(5) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(6) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(7) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + t449 * (t487 * t474 + t481 * t477) / 0.2e1 + t450 * (t481 * t474 - t487 * t477) / 0.2e1 + t460 * ((-t378 * t433 - t379 * t432 - t405 * t460) * t469 + ((-t381 * t472 + t383 * t475) * t432 + (-t380 * t472 + t382 * t475) * t433 + (-t406 * t472 + t407 * t475) * t460) * t467) / 0.2e1 + t433 * ((t379 * t520 + t381 * t438 + t383 * t439) * t432 + (t378 * t520 + t438 * t380 + t439 * t382) * t433 + (t405 * t520 + t406 * t438 + t407 * t439) * t460) / 0.2e1 + t432 * ((t379 * t519 + t440 * t381 + t441 * t383) * t432 + (t378 * t519 + t380 * t440 + t382 * t441) * t433 + (t405 * t519 + t406 * t440 + t407 * t441) * t460) / 0.2e1 - ((-t477 * t451 + t491 * t474) * qJD(1) + (t477 ^ 2 * t426 + (t492 * t474 + (-t427 + t493) * t477) * t474) * qJD(2)) * t507 / 0.2e1 + ((t474 * t451 + t477 * t491) * qJD(1) + (t474 ^ 2 * t427 + (t493 * t477 + (-t426 + t492) * t474) * t477) * qJD(2)) * t465 / 0.2e1 + ((t423 * t537 + t424 * t535 + t519 * t536) * t442 + (t423 * t543 + t539 * t424 + t541 * t519) * t399 + (t542 * t423 + t538 * t424 + t540 * t519) * t398) * t398 / 0.2e1 + ((t421 * t537 + t422 * t535 + t520 * t536) * t442 + (t543 * t421 + t539 * t422 + t541 * t520) * t399 + (t421 * t542 + t422 * t538 + t520 * t540) * t398) * t399 / 0.2e1 + ((-t398 * t540 - t399 * t541 - t442 * t536) * t469 + ((t466 * t537 + t468 * t535) * t442 + (t466 * t543 + t539 * t468) * t399 + (t466 * t542 + t468 * t538) * t398) * t467) * t442 / 0.2e1 + (m(2) * (t455 ^ 2 + t456 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t429 * t476 + t431 * t473) * t474 - (t428 * t476 + t430 * t473) * t477) * qJD(2) + (t415 * t469 + t417 * t467) * t449 + (t414 * t469 + t416 * t467) * t450 + (t469 * t444 + t467 * t445 + t452 * t476 + t473 * t453) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
