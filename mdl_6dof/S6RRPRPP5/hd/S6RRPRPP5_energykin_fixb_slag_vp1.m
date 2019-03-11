% Calculate kinetic energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:28
% EndTime: 2019-03-09 10:03:31
% DurationCPUTime: 2.67s
% Computational Cost: add. (880->215), mult. (2094->314), div. (0->0), fcn. (2075->6), ass. (0->123)
t541 = Icges(3,4) + Icges(4,6);
t540 = Icges(3,1) + Icges(4,2);
t539 = -Icges(3,2) - Icges(4,3);
t444 = cos(qJ(2));
t538 = t541 * t444;
t441 = sin(qJ(2));
t537 = t541 * t441;
t536 = -Icges(4,4) + Icges(3,5);
t535 = Icges(4,5) - Icges(3,6);
t534 = t539 * t441 + t538;
t533 = -t540 * t444 + t537;
t532 = Icges(4,1) + Icges(3,3);
t442 = sin(qJ(1));
t445 = cos(qJ(1));
t531 = t534 * t442 + t535 * t445;
t530 = -t535 * t442 + t534 * t445;
t529 = t533 * t442 + t536 * t445;
t528 = t536 * t442 - t533 * t445;
t527 = t539 * t444 - t537;
t526 = t540 * t441 + t538;
t525 = t535 * t441 + t536 * t444;
t524 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t523 = Icges(5,4) - Icges(7,4) - Icges(6,5);
t522 = -Icges(7,5) + Icges(6,4) + Icges(5,5);
t521 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t520 = Icges(6,6) - Icges(7,6) - Icges(5,6);
t519 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t518 = rSges(7,1) + pkin(5);
t517 = -rSges(7,3) - qJ(6);
t443 = cos(qJ(4));
t484 = t443 * t445;
t440 = sin(qJ(4));
t488 = t440 * t442;
t409 = -t441 * t484 + t488;
t486 = t442 * t443;
t487 = t440 * t445;
t410 = t441 * t487 + t486;
t369 = pkin(4) * t410 + qJ(5) * t409;
t411 = t441 * t486 + t487;
t439 = qJD(4) * t441 + qJD(1);
t516 = -qJD(5) * t411 + t439 * t369;
t414 = (-pkin(4) * t440 + qJ(5) * t443) * t444;
t475 = qJD(4) * t444;
t477 = qJD(2) * t445;
t418 = t442 * t475 - t477;
t515 = qJD(5) * t409 + t418 * t414;
t514 = t525 * t442 - t445 * t532;
t513 = t442 * t532 + t525 * t445;
t512 = t536 * t441 - t535 * t444;
t511 = t441 * t527 + t444 * t526;
t510 = -t441 * t530 + t444 * t528;
t509 = t441 * t531 + t444 * t529;
t483 = t444 * t445;
t508 = t409 * t520 + t410 * t522 + t483 * t519;
t412 = t441 * t488 - t484;
t485 = t442 * t444;
t507 = -t411 * t520 + t412 * t522 + t485 * t519;
t506 = t409 * t521 - t410 * t523 + t483 * t520;
t505 = -t411 * t521 - t412 * t523 + t485 * t520;
t504 = -t409 * t523 + t410 * t524 + t483 * t522;
t503 = t411 * t523 + t412 * t524 + t485 * t522;
t502 = (-t440 * t522 + t443 * t520) * t444 + t519 * t441;
t501 = (t440 * t523 + t443 * t521) * t444 + t520 * t441;
t500 = (-t440 * t524 - t443 * t523) * t444 + t522 * t441;
t482 = rSges(7,2) * t409 + t518 * t410 + t517 * t483;
t481 = -rSges(7,2) * t411 + t518 * t412 + t517 * t485;
t480 = (rSges(7,2) * t443 - t518 * t440) * t444 + t517 * t441;
t466 = pkin(2) * t444 + qJ(3) * t441;
t413 = t466 * t442;
t434 = pkin(1) * t442 - pkin(7) * t445;
t479 = -t413 - t434;
t478 = qJD(2) * t442;
t476 = qJD(3) * t441;
t415 = t466 * t445;
t422 = qJD(1) * (pkin(1) * t445 + pkin(7) * t442);
t474 = qJD(1) * t415 + t442 * t476 + t422;
t429 = pkin(2) * t441 - qJ(3) * t444;
t471 = qJD(2) * (rSges(4,2) * t441 + rSges(4,3) * t444 - t429);
t420 = pkin(3) * t442 + pkin(8) * t483;
t470 = qJD(1) * t420 + t474;
t469 = -qJD(3) * t444 + t413 * t478 + t415 * t477;
t468 = rSges(3,1) * t444 - rSges(3,2) * t441;
t467 = -rSges(4,2) * t444 + rSges(4,3) * t441;
t465 = (-pkin(8) * t441 - t429) * qJD(2);
t421 = -pkin(3) * t445 + pkin(8) * t485;
t438 = t445 * t476;
t452 = t438 + (-t421 + t479) * qJD(1);
t451 = t420 * t477 + t421 * t478 + t469;
t450 = -qJD(6) * t444 + t465;
t370 = pkin(4) * t412 - qJ(5) * t411;
t417 = t445 * t475 + t478;
t449 = qJD(5) * t444 * t443 + t417 * t370 + t451;
t448 = t442 * t465 + t470;
t447 = t445 * t465 + t452;
t433 = rSges(2,1) * t445 - rSges(2,2) * t442;
t432 = rSges(2,1) * t442 + rSges(2,2) * t445;
t431 = rSges(3,1) * t441 + rSges(3,2) * t444;
t403 = -rSges(4,1) * t445 + t442 * t467;
t402 = rSges(4,1) * t442 + t445 * t467;
t401 = rSges(3,3) * t442 + t445 * t468;
t400 = rSges(5,3) * t441 + (-rSges(5,1) * t440 - rSges(5,2) * t443) * t444;
t399 = rSges(6,2) * t441 + (-rSges(6,1) * t440 + rSges(6,3) * t443) * t444;
t397 = -rSges(3,3) * t445 + t442 * t468;
t368 = rSges(5,1) * t412 + rSges(5,2) * t411 + rSges(5,3) * t485;
t367 = rSges(6,1) * t412 + rSges(6,2) * t485 - rSges(6,3) * t411;
t365 = rSges(5,1) * t410 - rSges(5,2) * t409 + rSges(5,3) * t483;
t364 = rSges(6,1) * t410 + rSges(6,2) * t483 + rSges(6,3) * t409;
t343 = qJD(1) * t401 - t431 * t478 + t422;
t342 = -t431 * t477 + (-t397 - t434) * qJD(1);
t340 = (t397 * t442 + t401 * t445) * qJD(2);
t339 = qJD(1) * t402 + t442 * t471 + t474;
t338 = t438 + t445 * t471 + (-t403 + t479) * qJD(1);
t337 = (t402 * t445 + t403 * t442) * qJD(2) + t469;
t336 = t365 * t439 - t400 * t417 + t448;
t335 = -t368 * t439 + t400 * t418 + t447;
t334 = -t365 * t418 + t368 * t417 + t451;
t333 = t364 * t439 + (-t399 - t414) * t417 + t448 + t516;
t332 = t399 * t418 + (-t367 - t370) * t439 + t447 + t515;
t331 = t367 * t417 + (-t364 - t369) * t418 + t449;
t330 = t482 * t439 + t450 * t442 + (-t414 - t480) * t417 + t470 + t516;
t329 = t480 * t418 + t450 * t445 + (-t370 - t481) * t439 + t452 + t515;
t328 = -qJD(6) * t441 + t481 * t417 + (-t369 - t482) * t418 + t449;
t1 = m(6) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(7) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(5) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(3) * (t340 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(4) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t432 ^ 2 + t433 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t441 * t529 - t444 * t531) * t445 + (t441 * t528 + t444 * t530) * t442) * qJD(2) + (t526 * t441 - t527 * t444) * qJD(1)) * qJD(1) / 0.2e1 + ((t513 * t442 ^ 2 + (t509 * t445 + (t510 - t514) * t442) * t445) * qJD(2) + (t442 * t512 + t445 * t511) * qJD(1)) * t478 / 0.2e1 - ((t514 * t445 ^ 2 + (t510 * t442 + (t509 - t513) * t445) * t442) * qJD(2) + (t442 * t511 - t445 * t512) * qJD(1)) * t477 / 0.2e1 + ((t409 * t501 + t410 * t500 + t483 * t502) * t439 + (t409 * t505 + t410 * t503 + t483 * t507) * t418 + (t506 * t409 + t504 * t410 + t508 * t483) * t417) * t417 / 0.2e1 + ((-t411 * t501 + t412 * t500 + t485 * t502) * t439 + (-t505 * t411 + t503 * t412 + t507 * t485) * t418 + (-t411 * t506 + t412 * t504 + t485 * t508) * t417) * t418 / 0.2e1 + (((-t440 * t500 + t443 * t501) * t439 + (-t440 * t503 + t443 * t505) * t418 + (-t440 * t504 + t443 * t506) * t417) * t444 + (t417 * t508 + t507 * t418 + t502 * t439) * t441) * t439 / 0.2e1;
T  = t1;
