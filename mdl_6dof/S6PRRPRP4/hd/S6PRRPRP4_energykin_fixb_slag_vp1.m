% Calculate kinetic energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:14
% EndTime: 2019-03-08 21:40:16
% DurationCPUTime: 2.42s
% Computational Cost: add. (2224->283), mult. (5653->419), div. (0->0), fcn. (6872->10), ass. (0->137)
t577 = Icges(4,1) + Icges(5,2);
t576 = Icges(5,1) + Icges(4,3);
t575 = Icges(6,1) + Icges(7,1);
t574 = -Icges(4,4) - Icges(5,6);
t573 = Icges(5,4) - Icges(4,5);
t572 = Icges(6,4) + Icges(7,4);
t571 = Icges(5,5) - Icges(4,6);
t570 = Icges(6,5) + Icges(7,5);
t569 = Icges(4,2) + Icges(5,3);
t568 = Icges(6,2) + Icges(7,2);
t567 = Icges(6,6) + Icges(7,6);
t566 = Icges(6,3) + Icges(7,3);
t565 = rSges(7,3) + qJ(6);
t502 = sin(pkin(10));
t504 = cos(pkin(10));
t510 = cos(qJ(2));
t505 = cos(pkin(6));
t508 = sin(qJ(2));
t529 = t505 * t508;
t488 = t502 * t510 + t504 * t529;
t503 = sin(pkin(6));
t539 = cos(qJ(3));
t522 = t503 * t539;
t538 = sin(qJ(3));
t471 = t488 * t538 + t504 * t522;
t528 = t505 * t510;
t487 = t502 * t508 - t504 * t528;
t507 = sin(qJ(5));
t509 = cos(qJ(5));
t439 = t471 * t509 - t487 * t507;
t535 = t471 * t507;
t440 = t487 * t509 + t535;
t521 = t503 * t538;
t472 = t488 * t539 - t504 * t521;
t564 = t567 * t439 + t570 * t440 + t566 * t472;
t490 = -t502 * t529 + t504 * t510;
t473 = t490 * t538 - t502 * t522;
t489 = t502 * t528 + t504 * t508;
t441 = t473 * t509 - t489 * t507;
t534 = t473 * t507;
t442 = t489 * t509 + t534;
t474 = t490 * t539 + t502 * t521;
t563 = t567 * t441 + t570 * t442 + t566 * t474;
t562 = t568 * t439 + t572 * t440 + t567 * t472;
t561 = t568 * t441 + t572 * t442 + t567 * t474;
t560 = t572 * t439 + t575 * t440 + t570 * t472;
t559 = t572 * t441 + t575 * t442 + t570 * t474;
t558 = t569 * t471 + t574 * t472 + t571 * t487;
t557 = t569 * t473 + t574 * t474 + t571 * t489;
t556 = t571 * t471 - t573 * t472 + t576 * t487;
t555 = t571 * t473 - t573 * t474 + t576 * t489;
t554 = t574 * t471 + t577 * t472 - t573 * t487;
t553 = t574 * t473 + t577 * t474 - t573 * t489;
t491 = -t505 * t539 + t508 * t521;
t530 = t503 * t510;
t475 = t491 * t509 + t507 * t530;
t533 = t491 * t507;
t476 = -t509 * t530 + t533;
t492 = t505 * t538 + t508 * t522;
t552 = t567 * t475 + t570 * t476 + t566 * t492;
t551 = t568 * t475 + t572 * t476 + t567 * t492;
t550 = t572 * t475 + t575 * t476 + t570 * t492;
t549 = t569 * t491 + t574 * t492 - t571 * t530;
t548 = t574 * t491 + t577 * t492 + t573 * t530;
t547 = t571 * t491 - t573 * t492 - t576 * t530;
t546 = qJD(2) ^ 2;
t537 = pkin(5) * t509;
t532 = t502 * t503;
t531 = t503 * t504;
t527 = rSges(7,1) * t440 + rSges(7,2) * t439 + pkin(5) * t535 + t565 * t472 + t487 * t537;
t526 = rSges(7,1) * t442 + rSges(7,2) * t441 + pkin(5) * t534 + t565 * t474 + t489 * t537;
t525 = rSges(7,1) * t476 + rSges(7,2) * t475 + pkin(5) * t533 + t565 * t492 - t530 * t537;
t524 = qJD(2) * t503;
t498 = t502 * t524;
t477 = qJD(3) * t489 + t498;
t501 = qJD(2) * t505;
t520 = t504 * t524;
t465 = pkin(2) * t488 + pkin(8) * t487;
t466 = pkin(2) * t490 + pkin(8) * t489;
t519 = t465 * t498 + t466 * t520 + qJD(1);
t478 = qJD(3) * t487 - t520;
t494 = -qJD(3) * t530 + t501;
t435 = pkin(3) * t472 + qJ(4) * t471;
t518 = qJD(4) * t491 + t477 * t435 + t519;
t493 = (pkin(2) * t508 - pkin(8) * t510) * t503;
t517 = t466 * t501 - t493 * t498;
t436 = pkin(3) * t474 + qJ(4) * t473;
t516 = qJD(4) * t471 + t494 * t436 + t517;
t515 = (-t465 * t505 - t493 * t531) * qJD(2);
t444 = pkin(4) * t487 + pkin(9) * t472;
t445 = pkin(4) * t489 + pkin(9) * t474;
t514 = t477 * t444 + (-t436 - t445) * t478 + t518;
t467 = pkin(3) * t492 + qJ(4) * t491;
t513 = qJD(4) * t473 + t478 * t467 + t515;
t479 = -pkin(4) * t530 + t492 * pkin(9);
t512 = t494 * t445 + (-t467 - t479) * t477 + t516;
t511 = t478 * t479 + (-t435 - t444) * t494 + t513;
t483 = t505 * rSges(3,3) + (rSges(3,1) * t508 + rSges(3,2) * t510) * t503;
t482 = Icges(3,5) * t505 + (Icges(3,1) * t508 + Icges(3,4) * t510) * t503;
t481 = Icges(3,6) * t505 + (Icges(3,4) * t508 + Icges(3,2) * t510) * t503;
t480 = Icges(3,3) * t505 + (Icges(3,5) * t508 + Icges(3,6) * t510) * t503;
t470 = qJD(5) * t492 + t494;
t463 = t492 * rSges(4,1) - t491 * rSges(4,2) - rSges(4,3) * t530;
t462 = -rSges(5,1) * t530 - t492 * rSges(5,2) + t491 * rSges(5,3);
t453 = rSges(3,1) * t490 - rSges(3,2) * t489 + rSges(3,3) * t532;
t452 = rSges(3,1) * t488 - rSges(3,2) * t487 - rSges(3,3) * t531;
t451 = Icges(3,1) * t490 - Icges(3,4) * t489 + Icges(3,5) * t532;
t450 = Icges(3,1) * t488 - Icges(3,4) * t487 - Icges(3,5) * t531;
t449 = Icges(3,4) * t490 - Icges(3,2) * t489 + Icges(3,6) * t532;
t448 = Icges(3,4) * t488 - Icges(3,2) * t487 - Icges(3,6) * t531;
t447 = Icges(3,5) * t490 - Icges(3,6) * t489 + Icges(3,3) * t532;
t446 = Icges(3,5) * t488 - Icges(3,6) * t487 - Icges(3,3) * t531;
t438 = qJD(5) * t472 + t478;
t437 = qJD(5) * t474 + t477;
t429 = (-t452 * t505 - t483 * t531) * qJD(2);
t428 = (t453 * t505 - t483 * t532) * qJD(2);
t427 = rSges(6,1) * t476 + rSges(6,2) * t475 + rSges(6,3) * t492;
t419 = rSges(4,1) * t474 - rSges(4,2) * t473 + rSges(4,3) * t489;
t418 = rSges(4,1) * t472 - rSges(4,2) * t471 + rSges(4,3) * t487;
t417 = rSges(5,1) * t489 - rSges(5,2) * t474 + rSges(5,3) * t473;
t416 = rSges(5,1) * t487 - rSges(5,2) * t472 + rSges(5,3) * t471;
t402 = qJD(1) + (t452 * t502 + t453 * t504) * t524;
t399 = rSges(6,1) * t442 + rSges(6,2) * t441 + rSges(6,3) * t474;
t397 = rSges(6,1) * t440 + rSges(6,2) * t439 + rSges(6,3) * t472;
t383 = -t418 * t494 + t463 * t478 + t515;
t382 = t419 * t494 - t463 * t477 + t517;
t381 = t418 * t477 - t419 * t478 + t519;
t380 = t462 * t478 + (-t416 - t435) * t494 + t513;
t379 = t417 * t494 + (-t462 - t467) * t477 + t516;
t378 = t416 * t477 + (-t417 - t436) * t478 + t518;
t377 = -t397 * t470 + t427 * t438 + t511;
t376 = t399 * t470 - t427 * t437 + t512;
t375 = t397 * t437 - t399 * t438 + t514;
t374 = qJD(6) * t474 + t438 * t525 - t470 * t527 + t511;
t373 = qJD(6) * t472 - t437 * t525 + t470 * t526 + t512;
t372 = qJD(6) * t492 + t437 * t527 - t438 * t526 + t514;
t1 = m(3) * (t402 ^ 2 + t428 ^ 2 + t429 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(4) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(5) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(6) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(7) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 - t546 * ((-t447 * t531 - t449 * t487 + t451 * t488) * t532 - (-t446 * t531 - t448 * t487 + t450 * t488) * t531 + (-t480 * t531 - t481 * t487 + t482 * t488) * t505) * t531 / 0.2e1 + ((t551 * t441 + t550 * t442 + t552 * t474) * t470 + (t562 * t441 + t560 * t442 + t564 * t474) * t438 + (t561 * t441 + t559 * t442 + t563 * t474) * t437) * t437 / 0.2e1 + ((t551 * t439 + t550 * t440 + t552 * t472) * t470 + (t562 * t439 + t560 * t440 + t564 * t472) * t438 + (t561 * t439 + t559 * t440 + t563 * t472) * t437) * t438 / 0.2e1 + ((t551 * t475 + t550 * t476 + t552 * t492) * t470 + (t562 * t475 + t560 * t476 + t564 * t492) * t438 + (t561 * t475 + t559 * t476 + t563 * t492) * t437) * t470 / 0.2e1 + ((t549 * t473 + t548 * t474 + t547 * t489) * t494 + (t558 * t473 + t554 * t474 + t556 * t489) * t478 + (t557 * t473 + t553 * t474 + t555 * t489) * t477) * t477 / 0.2e1 + ((t549 * t471 + t548 * t472 + t547 * t487) * t494 + (t558 * t471 + t554 * t472 + t556 * t487) * t478 + (t557 * t471 + t553 * t472 + t555 * t487) * t477) * t478 / 0.2e1 + ((t549 * t491 + t548 * t492 - t547 * t530) * t494 + (t558 * t491 + t554 * t492 - t556 * t530) * t478 + (t557 * t491 + t553 * t492 - t555 * t530) * t477) * t494 / 0.2e1 + (t505 * (t505 ^ 2 * t480 + (((t449 * t510 + t451 * t508) * t502 - (t448 * t510 + t450 * t508) * t504) * t503 + (-t446 * t504 + t447 * t502 + t481 * t510 + t482 * t508) * t505) * t503) + ((t447 * t532 - t449 * t489 + t451 * t490) * t532 - (t446 * t532 - t448 * t489 + t450 * t490) * t531 + (t480 * t532 - t481 * t489 + t482 * t490) * t505) * t532) * t546 / 0.2e1;
T  = t1;
