% Calculate kinetic energy for
% S6RRRPRP9
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:50
% EndTime: 2019-03-09 17:21:53
% DurationCPUTime: 2.90s
% Computational Cost: add. (1359->266), mult. (3274->399), div. (0->0), fcn. (3640->8), ass. (0->138)
t561 = Icges(4,1) + Icges(5,1);
t560 = Icges(6,1) + Icges(7,1);
t559 = -Icges(4,4) + Icges(5,5);
t558 = Icges(5,4) + Icges(4,5);
t557 = -Icges(6,4) + Icges(7,5);
t556 = Icges(7,4) + Icges(6,5);
t555 = Icges(4,2) + Icges(5,3);
t554 = Icges(6,2) + Icges(7,3);
t553 = Icges(7,2) + Icges(6,3);
t552 = -Icges(5,6) + Icges(4,6);
t551 = Icges(6,6) - Icges(7,6);
t550 = -Icges(4,3) - Icges(5,2);
t549 = rSges(7,1) + pkin(5);
t548 = rSges(7,3) + qJ(6);
t484 = sin(qJ(3));
t487 = cos(qJ(3));
t489 = cos(qJ(1));
t486 = sin(qJ(1));
t488 = cos(qJ(2));
t515 = t486 * t488;
t456 = t484 * t515 + t487 * t489;
t457 = -t484 * t489 + t487 * t515;
t483 = sin(qJ(5));
t522 = cos(qJ(5));
t417 = -t456 * t522 + t457 * t483;
t418 = t456 * t483 + t457 * t522;
t485 = sin(qJ(2));
t518 = t485 * t486;
t547 = t554 * t417 + t557 * t418 + t551 * t518;
t514 = t488 * t489;
t458 = t484 * t514 - t486 * t487;
t459 = t484 * t486 + t487 * t514;
t419 = -t458 * t522 + t459 * t483;
t420 = t458 * t483 + t459 * t522;
t516 = t485 * t489;
t546 = t554 * t419 + t557 * t420 + t551 * t516;
t545 = -t551 * t417 + t556 * t418 - t553 * t518;
t544 = -t551 * t419 + t556 * t420 - t553 * t516;
t543 = t557 * t417 + t560 * t418 - t556 * t518;
t542 = t557 * t419 + t560 * t420 - t556 * t516;
t517 = t485 * t487;
t519 = t484 * t485;
t451 = t483 * t517 - t519 * t522;
t452 = (t483 * t484 + t487 * t522) * t485;
t541 = t554 * t451 + t557 * t452 - t551 * t488;
t540 = -t551 * t451 + t556 * t452 + t553 * t488;
t539 = t557 * t451 + t560 * t452 + t556 * t488;
t538 = t555 * t456 + t559 * t457 - t552 * t518;
t537 = t555 * t458 + t559 * t459 - t552 * t516;
t536 = -t552 * t456 + t558 * t457 - t550 * t518;
t535 = -t552 * t458 + t558 * t459 - t550 * t516;
t534 = t559 * t456 + t561 * t457 + t558 * t518;
t533 = t559 * t458 + t561 * t459 + t558 * t516;
t532 = t552 * t488 + (t555 * t484 + t559 * t487) * t485;
t531 = t550 * t488 + (-t552 * t484 + t558 * t487) * t485;
t530 = -t558 * t488 + (t559 * t484 + t561 * t487) * t485;
t521 = Icges(3,4) * t485;
t520 = Icges(3,4) * t488;
t513 = -rSges(7,2) * t518 + t548 * t417 + t549 * t418;
t512 = -rSges(7,2) * t516 + t548 * t419 + t549 * t420;
t511 = rSges(7,2) * t488 + t548 * t451 + t549 * t452;
t505 = pkin(2) * t488 + pkin(8) * t485;
t461 = t505 * t486;
t462 = t505 * t489;
t482 = qJD(2) * t486;
t509 = qJD(2) * t489;
t510 = t461 * t482 + t462 * t509;
t508 = qJD(3) * t485;
t463 = t489 * t508 + t482;
t507 = qJD(5) * t485;
t464 = t486 * t508 - t509;
t422 = pkin(3) * t457 + qJ(4) * t456;
t506 = qJD(4) * t519 + t463 * t422 + t510;
t504 = rSges(3,1) * t488 - rSges(3,2) * t485;
t503 = Icges(3,1) * t488 - t521;
t502 = -Icges(3,2) * t485 + t520;
t501 = Icges(3,5) * t488 - Icges(3,6) * t485;
t437 = -Icges(3,6) * t489 + t486 * t502;
t441 = -Icges(3,5) * t489 + t486 * t503;
t500 = t437 * t485 - t441 * t488;
t438 = Icges(3,6) * t486 + t489 * t502;
t442 = Icges(3,5) * t486 + t489 * t503;
t499 = -t438 * t485 + t442 * t488;
t469 = Icges(3,2) * t488 + t521;
t470 = Icges(3,1) * t485 + t520;
t498 = -t469 * t485 + t470 * t488;
t467 = qJD(1) * (pkin(1) * t489 + pkin(7) * t486);
t474 = pkin(2) * t485 - pkin(8) * t488;
t497 = qJD(1) * t462 - t474 * t482 + t467;
t423 = pkin(3) * t459 + qJ(4) * t458;
t480 = -qJD(3) * t488 + qJD(1);
t496 = qJD(4) * t456 + t480 * t423 + t497;
t475 = pkin(1) * t486 - pkin(7) * t489;
t495 = (-t461 - t475) * qJD(1) - t474 * t509;
t427 = pkin(4) * t457 - pkin(9) * t518;
t428 = pkin(4) * t459 - pkin(9) * t516;
t494 = t463 * t427 + (-t423 - t428) * t464 + t506;
t460 = (pkin(3) * t487 + qJ(4) * t484) * t485;
t493 = qJD(4) * t458 + t464 * t460 + t495;
t465 = pkin(4) * t517 + pkin(9) * t488;
t492 = t480 * t428 + (-t460 - t465) * t463 + t496;
t491 = t464 * t465 + (-t422 - t427) * t480 + t493;
t473 = rSges(2,1) * t489 - rSges(2,2) * t486;
t472 = rSges(2,1) * t486 + rSges(2,2) * t489;
t471 = rSges(3,1) * t485 + rSges(3,2) * t488;
t468 = Icges(3,5) * t485 + Icges(3,6) * t488;
t466 = qJD(1) + (-qJD(3) + qJD(5)) * t488;
t446 = rSges(3,3) * t486 + t489 * t504;
t445 = -rSges(3,3) * t489 + t486 * t504;
t444 = -rSges(4,3) * t488 + (rSges(4,1) * t487 - rSges(4,2) * t484) * t485;
t443 = -rSges(5,2) * t488 + (rSges(5,1) * t487 + rSges(5,3) * t484) * t485;
t434 = Icges(3,3) * t486 + t489 * t501;
t433 = -Icges(3,3) * t489 + t486 * t501;
t430 = -t486 * t507 + t464;
t429 = -t489 * t507 + t463;
t415 = rSges(4,1) * t459 - rSges(4,2) * t458 + rSges(4,3) * t516;
t414 = rSges(5,1) * t459 + rSges(5,2) * t516 + rSges(5,3) * t458;
t413 = rSges(4,1) * t457 - rSges(4,2) * t456 + rSges(4,3) * t518;
t412 = rSges(5,1) * t457 + rSges(5,2) * t518 + rSges(5,3) * t456;
t398 = rSges(6,1) * t452 - rSges(6,2) * t451 + rSges(6,3) * t488;
t390 = qJD(1) * t446 - t471 * t482 + t467;
t389 = -t471 * t509 + (-t445 - t475) * qJD(1);
t387 = (t445 * t486 + t446 * t489) * qJD(2);
t384 = rSges(6,1) * t420 - rSges(6,2) * t419 - rSges(6,3) * t516;
t382 = rSges(6,1) * t418 - rSges(6,2) * t417 - rSges(6,3) * t518;
t368 = t415 * t480 - t444 * t463 + t497;
t367 = -t413 * t480 + t444 * t464 + t495;
t366 = t413 * t463 - t415 * t464 + t510;
t365 = t414 * t480 + (-t443 - t460) * t463 + t496;
t364 = t443 * t464 + (-t412 - t422) * t480 + t493;
t363 = t412 * t463 + (-t414 - t423) * t464 + t506;
t362 = t384 * t466 - t398 * t429 + t492;
t361 = -t382 * t466 + t398 * t430 + t491;
t360 = t382 * t429 - t384 * t430 + t494;
t359 = qJD(6) * t417 - t429 * t511 + t466 * t512 + t492;
t358 = qJD(6) * t419 + t430 * t511 - t466 * t513 + t491;
t357 = qJD(6) * t451 + t429 * t513 - t430 * t512 + t494;
t1 = m(6) * (t360 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 + m(5) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + m(4) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(3) * (t387 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(7) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + qJD(1) * ((t469 * t488 + t470 * t485) * qJD(1) + ((t438 * t488 + t442 * t485) * t486 - (t437 * t488 + t441 * t485) * t489) * qJD(2)) / 0.2e1 - ((-t468 * t489 + t486 * t498) * qJD(1) + (t433 * t489 ^ 2 + (t499 * t486 + (-t434 + t500) * t489) * t486) * qJD(2)) * t509 / 0.2e1 + ((t468 * t486 + t489 * t498) * qJD(1) + (t434 * t486 ^ 2 + (t500 * t489 + (-t433 + t499) * t486) * t489) * qJD(2)) * t482 / 0.2e1 + ((t419 * t541 + t420 * t539 - t516 * t540) * t466 + (t419 * t547 + t543 * t420 - t545 * t516) * t430 + (t546 * t419 + t542 * t420 - t544 * t516) * t429) * t429 / 0.2e1 + ((t417 * t541 + t418 * t539 - t518 * t540) * t466 + (t547 * t417 + t543 * t418 - t545 * t518) * t430 + (t417 * t546 + t418 * t542 - t518 * t544) * t429) * t430 / 0.2e1 + ((t458 * t532 + t459 * t530 + t516 * t531) * t480 + (t458 * t538 + t459 * t534 + t516 * t536) * t464 + (t537 * t458 + t533 * t459 + t535 * t516) * t463) * t463 / 0.2e1 + ((t456 * t532 + t457 * t530 + t518 * t531) * t480 + (t538 * t456 + t534 * t457 + t536 * t518) * t464 + (t456 * t537 + t457 * t533 + t518 * t535) * t463) * t464 / 0.2e1 + ((t541 * t451 + t539 * t452 + t540 * t488) * t466 + (t451 * t547 + t543 * t452 + t545 * t488) * t430 + (t451 * t546 + t452 * t542 + t488 * t544) * t429) * t466 / 0.2e1 + ((-t463 * t535 - t464 * t536 - t480 * t531) * t488 + ((t484 * t532 + t487 * t530) * t480 + (t484 * t538 + t487 * t534) * t464 + (t484 * t537 + t487 * t533) * t463) * t485) * t480 / 0.2e1 + (m(2) * (t472 ^ 2 + t473 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
