% Calculate kinetic energy for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:14
% EndTime: 2019-03-09 17:14:16
% DurationCPUTime: 2.86s
% Computational Cost: add. (1366->278), mult. (3256->413), div. (0->0), fcn. (3599->8), ass. (0->144)
t556 = Icges(4,1) + Icges(5,1);
t555 = Icges(6,1) + Icges(7,1);
t554 = -Icges(4,4) + Icges(5,5);
t553 = Icges(5,4) + Icges(4,5);
t552 = Icges(6,4) + Icges(7,4);
t551 = Icges(6,5) + Icges(7,5);
t550 = Icges(4,2) + Icges(5,3);
t549 = Icges(6,2) + Icges(7,2);
t548 = -Icges(5,6) + Icges(4,6);
t547 = Icges(6,6) + Icges(7,6);
t546 = -Icges(4,3) - Icges(5,2);
t545 = Icges(6,3) + Icges(7,3);
t472 = sin(qJ(3));
t476 = cos(qJ(3));
t478 = cos(qJ(1));
t474 = sin(qJ(1));
t477 = cos(qJ(2));
t508 = t477 * t474;
t443 = t472 * t508 + t476 * t478;
t444 = -t472 * t478 + t476 * t508;
t471 = sin(qJ(5));
t475 = cos(qJ(5));
t406 = t443 * t475 - t444 * t471;
t513 = t443 * t471;
t407 = t444 * t475 + t513;
t473 = sin(qJ(2));
t510 = t473 * t474;
t544 = t547 * t406 + t551 * t407 - t545 * t510;
t507 = t477 * t478;
t445 = t472 * t507 - t474 * t476;
t446 = t472 * t474 + t476 * t507;
t408 = t445 * t475 - t446 * t471;
t512 = t445 * t471;
t409 = t446 * t475 + t512;
t509 = t473 * t478;
t543 = t547 * t408 + t551 * t409 - t545 * t509;
t542 = t549 * t406 + t552 * t407 - t547 * t510;
t541 = t549 * t408 + t552 * t409 - t547 * t509;
t540 = t552 * t406 + t555 * t407 - t551 * t510;
t539 = t552 * t408 + t555 * t409 - t551 * t509;
t438 = (-t471 * t476 + t472 * t475) * t473;
t511 = t471 * t472;
t439 = (t475 * t476 + t511) * t473;
t538 = t547 * t438 + t551 * t439 + t545 * t477;
t537 = t549 * t438 + t552 * t439 + t547 * t477;
t536 = t552 * t438 + t555 * t439 + t551 * t477;
t535 = t550 * t443 + t554 * t444 - t548 * t510;
t534 = t550 * t445 + t554 * t446 - t548 * t509;
t533 = -t548 * t443 + t553 * t444 - t546 * t510;
t532 = -t548 * t445 + t553 * t446 - t546 * t509;
t531 = t554 * t443 + t556 * t444 + t553 * t510;
t530 = t554 * t445 + t556 * t446 + t553 * t509;
t412 = pkin(3) * t446 + qJ(4) * t445;
t466 = -qJD(3) * t477 + qJD(1);
t529 = qJD(4) * t443 + t466 * t412;
t447 = (pkin(3) * t476 + qJ(4) * t472) * t473;
t500 = qJD(3) * t473;
t501 = qJD(2) * t478;
t451 = t474 * t500 - t501;
t528 = qJD(4) * t445 + t451 * t447;
t527 = t548 * t477 + (t550 * t472 + t554 * t476) * t473;
t526 = t546 * t477 + (-t548 * t472 + t553 * t476) * t473;
t525 = -t553 * t477 + (t554 * t472 + t556 * t476) * t473;
t517 = pkin(5) * t475;
t515 = Icges(3,4) * t473;
t514 = Icges(3,4) * t477;
t496 = t473 * qJ(6);
t506 = rSges(7,1) * t407 + rSges(7,2) * t406 - rSges(7,3) * t510 + pkin(5) * t513 + t444 * t517 - t474 * t496;
t505 = rSges(7,1) * t409 + rSges(7,2) * t408 - rSges(7,3) * t509 + pkin(5) * t512 + t446 * t517 - t478 * t496;
t504 = rSges(7,1) * t439 + rSges(7,2) * t438 + (pkin(5) * t511 + t476 * t517) * t473 + (rSges(7,3) + qJ(6)) * t477;
t493 = pkin(2) * t477 + pkin(8) * t473;
t448 = t493 * t474;
t449 = t493 * t478;
t469 = qJD(2) * t474;
t503 = t448 * t469 + t449 * t501;
t454 = qJD(1) * (pkin(1) * t478 + pkin(7) * t474);
t502 = qJD(1) * t449 + t454;
t450 = t478 * t500 + t469;
t499 = qJD(5) * t473;
t461 = pkin(2) * t473 - pkin(8) * t477;
t498 = t461 * t469;
t497 = t461 * t501;
t462 = pkin(1) * t474 - pkin(7) * t478;
t495 = (-t448 - t462) * qJD(1);
t411 = pkin(3) * t444 + qJ(4) * t443;
t494 = qJD(4) * t473 * t472 + t450 * t411 + t503;
t492 = rSges(3,1) * t477 - rSges(3,2) * t473;
t491 = Icges(3,1) * t477 - t515;
t490 = -Icges(3,2) * t473 + t514;
t489 = Icges(3,5) * t477 - Icges(3,6) * t473;
t426 = -Icges(3,6) * t478 + t474 * t490;
t430 = -Icges(3,5) * t478 + t474 * t491;
t488 = t426 * t473 - t430 * t477;
t427 = Icges(3,6) * t474 + t478 * t490;
t431 = Icges(3,5) * t474 + t478 * t491;
t487 = -t427 * t473 + t431 * t477;
t456 = Icges(3,2) * t477 + t515;
t457 = Icges(3,1) * t473 + t514;
t486 = -t456 * t473 + t457 * t477;
t485 = -qJD(2) * t461 - qJD(6) * t473;
t484 = -t498 + t502;
t483 = t495 - t497;
t416 = pkin(4) * t444 - pkin(9) * t510;
t417 = pkin(4) * t446 - pkin(9) * t509;
t482 = t450 * t416 + (-t412 - t417) * t451 + t494;
t452 = pkin(4) * t473 * t476 + pkin(9) * t477;
t481 = t466 * t417 + (-t447 - t452) * t450 + t502 + t529;
t480 = t451 * t452 + (-t411 - t416) * t466 + t495 + t528;
t460 = rSges(2,1) * t478 - rSges(2,2) * t474;
t459 = rSges(2,1) * t474 + rSges(2,2) * t478;
t458 = rSges(3,1) * t473 + rSges(3,2) * t477;
t455 = Icges(3,5) * t473 + Icges(3,6) * t477;
t453 = qJD(1) + (-qJD(3) + qJD(5)) * t477;
t435 = rSges(3,3) * t474 + t478 * t492;
t434 = -rSges(3,3) * t478 + t474 * t492;
t433 = -rSges(4,3) * t477 + (rSges(4,1) * t476 - rSges(4,2) * t472) * t473;
t432 = -rSges(5,2) * t477 + (rSges(5,1) * t476 + rSges(5,3) * t472) * t473;
t423 = Icges(3,3) * t474 + t478 * t489;
t422 = -Icges(3,3) * t478 + t474 * t489;
t419 = -t474 * t499 + t451;
t418 = -t478 * t499 + t450;
t404 = rSges(4,1) * t446 - rSges(4,2) * t445 + rSges(4,3) * t509;
t403 = rSges(5,1) * t446 + rSges(5,2) * t509 + rSges(5,3) * t445;
t402 = rSges(4,1) * t444 - rSges(4,2) * t443 + rSges(4,3) * t510;
t401 = rSges(5,1) * t444 + rSges(5,2) * t510 + rSges(5,3) * t443;
t387 = rSges(6,1) * t439 + rSges(6,2) * t438 + rSges(6,3) * t477;
t379 = qJD(1) * t435 - t458 * t469 + t454;
t378 = -t458 * t501 + (-t434 - t462) * qJD(1);
t376 = (t434 * t474 + t435 * t478) * qJD(2);
t373 = rSges(6,1) * t409 + rSges(6,2) * t408 - rSges(6,3) * t509;
t371 = rSges(6,1) * t407 + rSges(6,2) * t406 - rSges(6,3) * t510;
t357 = t404 * t466 - t433 * t450 + t484;
t356 = -t402 * t466 + t433 * t451 + t483;
t355 = t402 * t450 - t404 * t451 + t503;
t354 = t403 * t466 + (-t432 - t447) * t450 + t484 + t529;
t353 = t432 * t451 + (-t401 - t411) * t466 + t483 + t528;
t352 = t401 * t450 + (-t403 - t412) * t451 + t494;
t351 = t373 * t453 - t387 * t418 + t481 - t498;
t350 = -t371 * t453 + t387 * t419 + t480 - t497;
t349 = t371 * t418 - t373 * t419 + t482;
t348 = -t418 * t504 + t453 * t505 + t474 * t485 + t481;
t347 = t419 * t504 - t453 * t506 + t478 * t485 + t480;
t346 = qJD(6) * t477 + t418 * t506 - t419 * t505 + t482;
t1 = m(3) * (t376 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(4) * (t355 ^ 2 + t356 ^ 2 + t357 ^ 2) / 0.2e1 + m(5) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(6) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(7) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + qJD(1) * ((t456 * t477 + t457 * t473) * qJD(1) + ((t427 * t477 + t431 * t473) * t474 - (t426 * t477 + t430 * t473) * t478) * qJD(2)) / 0.2e1 - ((-t455 * t478 + t474 * t486) * qJD(1) + (t422 * t478 ^ 2 + (t487 * t474 + (-t423 + t488) * t478) * t474) * qJD(2)) * t501 / 0.2e1 + ((t455 * t474 + t478 * t486) * qJD(1) + (t423 * t474 ^ 2 + (t488 * t478 + (-t422 + t487) * t474) * t478) * qJD(2)) * t469 / 0.2e1 + ((t537 * t408 + t536 * t409 - t538 * t509) * t453 + (t542 * t408 + t540 * t409 - t544 * t509) * t419 + (t541 * t408 + t539 * t409 - t543 * t509) * t418) * t418 / 0.2e1 + ((t537 * t406 + t536 * t407 - t538 * t510) * t453 + (t542 * t406 + t540 * t407 - t544 * t510) * t419 + (t541 * t406 + t539 * t407 - t543 * t510) * t418) * t419 / 0.2e1 + ((t527 * t445 + t525 * t446 + t526 * t509) * t466 + (t535 * t445 + t531 * t446 + t533 * t509) * t451 + (t534 * t445 + t530 * t446 + t532 * t509) * t450) * t450 / 0.2e1 + ((t527 * t443 + t525 * t444 + t526 * t510) * t466 + (t535 * t443 + t531 * t444 + t533 * t510) * t451 + (t534 * t443 + t530 * t444 + t532 * t510) * t450) * t451 / 0.2e1 + ((t537 * t438 + t536 * t439 + t538 * t477) * t453 + (t542 * t438 + t540 * t439 + t544 * t477) * t419 + (t541 * t438 + t539 * t439 + t543 * t477) * t418) * t453 / 0.2e1 + ((-t532 * t450 - t533 * t451 - t526 * t466) * t477 + ((t527 * t472 + t525 * t476) * t466 + (t535 * t472 + t531 * t476) * t451 + (t534 * t472 + t530 * t476) * t450) * t473) * t466 / 0.2e1 + (m(2) * (t459 ^ 2 + t460 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
