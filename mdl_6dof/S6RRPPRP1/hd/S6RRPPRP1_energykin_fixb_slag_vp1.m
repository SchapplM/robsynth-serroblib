% Calculate kinetic energy for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:23
% EndTime: 2019-03-09 08:25:26
% DurationCPUTime: 2.52s
% Computational Cost: add. (1864->280), mult. (2096->410), div. (0->0), fcn. (2063->10), ass. (0->151)
t563 = Icges(6,1) + Icges(7,1);
t562 = -Icges(6,4) + Icges(7,5);
t561 = Icges(7,4) + Icges(6,5);
t560 = Icges(6,2) + Icges(7,3);
t559 = -Icges(7,6) + Icges(6,6);
t558 = Icges(3,3) + Icges(4,3);
t557 = -Icges(6,3) - Icges(7,2);
t465 = qJ(2) + pkin(9);
t460 = sin(t465);
t462 = cos(t465);
t470 = sin(qJ(2));
t472 = cos(qJ(2));
t556 = Icges(3,5) * t472 + Icges(4,5) * t462 - Icges(3,6) * t470 - Icges(4,6) * t460;
t555 = rSges(7,1) + pkin(5);
t554 = rSges(7,3) + qJ(6);
t464 = pkin(10) + qJ(5);
t461 = cos(t464);
t473 = cos(qJ(1));
t459 = sin(t464);
t471 = sin(qJ(1));
t523 = t459 * t471;
t419 = t461 * t473 + t462 * t523;
t515 = t471 * t461;
t420 = -t459 * t473 + t462 * t515;
t522 = t460 * t471;
t553 = t560 * t419 + t562 * t420 - t559 * t522;
t520 = t462 * t473;
t421 = t459 * t520 - t515;
t422 = t461 * t520 + t523;
t521 = t460 * t473;
t552 = t560 * t421 + t562 * t422 - t559 * t521;
t551 = -t559 * t419 + t561 * t420 - t557 * t522;
t550 = -t559 * t421 + t561 * t422 - t557 * t521;
t549 = t562 * t419 + t563 * t420 + t561 * t522;
t548 = t562 * t421 + t563 * t422 + t561 * t521;
t547 = t559 * t462 + (t560 * t459 + t562 * t461) * t460;
t546 = t557 * t462 + (-t559 * t459 + t561 * t461) * t460;
t545 = -t561 * t462 + (t562 * t459 + t563 * t461) * t460;
t544 = t556 * t471 - t558 * t473;
t543 = t558 * t471 + t556 * t473;
t542 = Icges(3,5) * t470 + Icges(4,5) * t460 + Icges(3,6) * t472 + Icges(4,6) * t462;
t525 = Icges(4,4) * t460;
t441 = Icges(4,2) * t462 + t525;
t524 = Icges(4,4) * t462;
t442 = Icges(4,1) * t460 + t524;
t527 = Icges(3,4) * t470;
t447 = Icges(3,2) * t472 + t527;
t526 = Icges(3,4) * t472;
t448 = Icges(3,1) * t470 + t526;
t541 = -t441 * t460 + t442 * t462 - t447 * t470 + t448 * t472;
t491 = -Icges(4,2) * t460 + t524;
t411 = Icges(4,6) * t471 + t473 * t491;
t493 = Icges(4,1) * t462 - t525;
t413 = Icges(4,5) * t471 + t473 * t493;
t492 = -Icges(3,2) * t470 + t526;
t427 = Icges(3,6) * t471 + t473 * t492;
t494 = Icges(3,1) * t472 - t527;
t429 = Icges(3,5) * t471 + t473 * t494;
t540 = -t411 * t460 + t413 * t462 - t427 * t470 + t429 * t472;
t410 = -Icges(4,6) * t473 + t471 * t491;
t412 = -Icges(4,5) * t473 + t471 * t493;
t426 = -Icges(3,6) * t473 + t471 * t492;
t428 = -Icges(3,5) * t473 + t471 * t494;
t539 = t410 * t460 - t412 * t462 + t426 * t470 - t428 * t472;
t532 = pkin(2) * t470;
t530 = pkin(2) * t472;
t467 = cos(pkin(10));
t529 = pkin(4) * t467;
t466 = sin(pkin(10));
t519 = t466 * t471;
t518 = t466 * t473;
t517 = t467 * t471;
t516 = t467 * t473;
t513 = rSges(7,2) * t522 + t554 * t419 + t420 * t555;
t512 = rSges(7,2) * t521 + t554 * t421 + t422 * t555;
t511 = -t462 * rSges(7,2) + (t554 * t459 + t461 * t555) * t460;
t406 = -qJ(3) * t473 + t471 * t530;
t407 = qJ(3) * t471 + t473 * t530;
t506 = qJD(2) * t473;
t507 = qJD(2) * t471;
t510 = t406 * t507 + t407 * t506;
t454 = pkin(1) * t471 - pkin(7) * t473;
t509 = -t406 - t454;
t463 = qJD(3) * t471;
t505 = qJD(4) * t460;
t508 = t473 * t505 + t463;
t504 = qJD(5) * t460;
t496 = pkin(3) * t462 + qJ(4) * t460;
t430 = t496 * t471;
t503 = -t430 + t509;
t500 = -pkin(3) * t460 + qJ(4) * t462 - t532;
t445 = qJD(1) * (pkin(1) * t473 + pkin(7) * t471);
t499 = qJD(1) * t407 - qJD(3) * t473 + t445;
t498 = rSges(3,1) * t472 - rSges(3,2) * t470;
t497 = rSges(4,1) * t462 - rSges(4,2) * t460;
t495 = qJD(2) * (-rSges(4,1) * t460 - rSges(4,2) * t462 - t532);
t482 = qJD(2) * (pkin(8) * t462 - t460 * t529 + t500);
t481 = qJD(2) * (t462 * rSges(5,3) - (rSges(5,1) * t467 - rSges(5,2) * t466) * t460 + t500);
t431 = t496 * t473;
t480 = qJD(1) * t431 + t471 * t505 + t499;
t479 = -qJD(4) * t462 + t430 * t507 + t431 * t506 + t510;
t478 = pkin(8) * t460 + t462 * t529;
t375 = -pkin(4) * t518 + t471 * t478;
t376 = pkin(4) * t519 + t473 * t478;
t477 = t375 * t507 + t376 * t506 + t479;
t476 = qJD(1) * t376 + t471 * t482 + t480;
t475 = (-t375 + t503) * qJD(1) + t473 * t482 + t508;
t455 = -qJD(5) * t462 + qJD(1);
t453 = rSges(2,1) * t473 - rSges(2,2) * t471;
t452 = rSges(2,1) * t471 + rSges(2,2) * t473;
t451 = rSges(3,1) * t470 + rSges(3,2) * t472;
t439 = t471 * t504 - t506;
t438 = t473 * t504 + t507;
t437 = t462 * t516 + t519;
t436 = -t462 * t518 + t517;
t435 = t462 * t517 - t518;
t434 = -t462 * t519 - t516;
t433 = t471 * rSges(3,3) + t473 * t498;
t432 = -t473 * rSges(3,3) + t471 * t498;
t418 = t471 * rSges(4,3) + t473 * t497;
t417 = -t473 * rSges(4,3) + t471 * t497;
t403 = -Icges(5,5) * t462 + (Icges(5,1) * t467 - Icges(5,4) * t466) * t460;
t402 = -Icges(5,6) * t462 + (Icges(5,4) * t467 - Icges(5,2) * t466) * t460;
t401 = -Icges(5,3) * t462 + (Icges(5,5) * t467 - Icges(5,6) * t466) * t460;
t398 = -t462 * rSges(6,3) + (rSges(6,1) * t461 - rSges(6,2) * t459) * t460;
t389 = qJD(1) * t433 - t451 * t507 + t445;
t388 = -t451 * t506 + (-t432 - t454) * qJD(1);
t385 = (t432 * t471 + t433 * t473) * qJD(2);
t384 = rSges(5,1) * t437 + rSges(5,2) * t436 + rSges(5,3) * t521;
t383 = rSges(5,1) * t435 + rSges(5,2) * t434 + rSges(5,3) * t522;
t382 = Icges(5,1) * t437 + Icges(5,4) * t436 + Icges(5,5) * t521;
t381 = Icges(5,1) * t435 + Icges(5,4) * t434 + Icges(5,5) * t522;
t380 = Icges(5,4) * t437 + Icges(5,2) * t436 + Icges(5,6) * t521;
t379 = Icges(5,4) * t435 + Icges(5,2) * t434 + Icges(5,6) * t522;
t378 = Icges(5,5) * t437 + Icges(5,6) * t436 + Icges(5,3) * t521;
t377 = Icges(5,5) * t435 + Icges(5,6) * t434 + Icges(5,3) * t522;
t371 = rSges(6,1) * t422 - rSges(6,2) * t421 + rSges(6,3) * t521;
t369 = rSges(6,1) * t420 - rSges(6,2) * t419 + rSges(6,3) * t522;
t355 = qJD(1) * t418 + t471 * t495 + t499;
t354 = t463 + t473 * t495 + (-t417 + t509) * qJD(1);
t353 = (t417 * t471 + t418 * t473) * qJD(2) + t510;
t352 = qJD(1) * t384 + t471 * t481 + t480;
t351 = t473 * t481 + (-t383 + t503) * qJD(1) + t508;
t350 = (t383 * t471 + t384 * t473) * qJD(2) + t479;
t349 = t371 * t455 - t398 * t438 + t476;
t348 = -t369 * t455 + t398 * t439 + t475;
t347 = t369 * t438 - t371 * t439 + t477;
t346 = qJD(6) * t419 - t438 * t511 + t455 * t512 + t476;
t345 = qJD(6) * t421 + t439 * t511 - t455 * t513 + t475;
t344 = qJD(6) * t459 * t460 + t438 * t513 - t439 * t512 + t477;
t1 = m(7) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t385 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(4) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(5) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(6) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + ((t547 * t421 + t545 * t422 + t546 * t521) * t455 + (t553 * t421 + t549 * t422 + t551 * t521) * t439 + (t552 * t421 + t548 * t422 + t550 * t521) * t438) * t438 / 0.2e1 + ((t547 * t419 + t545 * t420 + t546 * t522) * t455 + (t553 * t419 + t549 * t420 + t551 * t522) * t439 + (t552 * t419 + t548 * t420 + t550 * t522) * t438) * t439 / 0.2e1 + ((-t550 * t438 - t551 * t439 - t546 * t455) * t462 + ((t547 * t459 + t545 * t461) * t455 + (t553 * t459 + t549 * t461) * t439 + (t552 * t459 + t548 * t461) * t438) * t460) * t455 / 0.2e1 + (m(2) * (t452 ^ 2 + t453 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t472 * t426 - t470 * t428 + (-t410 + t377) * t462 + (t379 * t466 - t381 * t467 - t412) * t460) * t473 + (t472 * t427 + t470 * t429 + (t411 - t378) * t462 + (-t380 * t466 + t382 * t467 + t413) * t460) * t471) * qJD(2) + (t472 * t447 + t470 * t448 + (t441 - t401) * t462 + (-t402 * t466 + t403 * t467 + t442) * t460) * qJD(1)) * qJD(1) / 0.2e1 + (((-t377 * t521 - t379 * t436 - t381 * t437 + t539 * t473) * t473 + (t378 * t521 + t436 * t380 + t437 * t382 + (t540 - t544) * t473 + t543 * t471) * t471) * qJD(2) + (t401 * t521 + t436 * t402 + t437 * t403 + t542 * t471 + t541 * t473) * qJD(1)) * t507 / 0.2e1 - (((t378 * t522 + t380 * t434 + t382 * t435 + t540 * t471) * t471 + (-t377 * t522 - t434 * t379 - t435 * t381 + (t539 - t543) * t471 + t544 * t473) * t473) * qJD(2) + (t401 * t522 + t434 * t402 + t435 * t403 + t541 * t471 - t542 * t473) * qJD(1)) * t506 / 0.2e1;
T  = t1;
