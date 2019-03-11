% Calculate kinetic energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:36:59
% EndTime: 2019-03-09 08:37:02
% DurationCPUTime: 2.47s
% Computational Cost: add. (1263->264), mult. (3114->381), div. (0->0), fcn. (3480->8), ass. (0->137)
t562 = Icges(4,1) + Icges(5,1);
t561 = Icges(6,1) + Icges(7,1);
t560 = Icges(4,4) - Icges(5,5);
t559 = Icges(5,4) + Icges(4,5);
t558 = -Icges(6,4) + Icges(7,5);
t557 = Icges(7,4) + Icges(6,5);
t556 = Icges(4,2) + Icges(5,3);
t555 = Icges(6,2) + Icges(7,3);
t554 = Icges(7,2) + Icges(6,3);
t553 = -Icges(5,6) + Icges(4,6);
t552 = Icges(6,6) - Icges(7,6);
t551 = -Icges(4,3) - Icges(5,2);
t550 = rSges(7,1) + pkin(5);
t549 = rSges(7,3) + qJ(6);
t478 = sin(pkin(9));
t479 = cos(pkin(9));
t484 = cos(qJ(1));
t482 = sin(qJ(1));
t483 = cos(qJ(2));
t518 = t482 * t483;
t453 = t478 * t518 + t479 * t484;
t454 = -t478 * t484 + t479 * t518;
t480 = sin(qJ(5));
t526 = cos(qJ(5));
t416 = -t453 * t526 + t454 * t480;
t417 = t453 * t480 + t454 * t526;
t481 = sin(qJ(2));
t520 = t481 * t482;
t548 = t416 * t555 + t417 * t558 + t520 * t552;
t517 = t483 * t484;
t455 = t478 * t517 - t482 * t479;
t456 = t478 * t482 + t479 * t517;
t418 = -t455 * t526 + t456 * t480;
t419 = t455 * t480 + t456 * t526;
t519 = t481 * t484;
t547 = t418 * t555 + t419 * t558 + t519 * t552;
t546 = -t416 * t552 + t417 * t557 - t520 * t554;
t545 = -t418 * t552 + t419 * t557 - t519 * t554;
t544 = t416 * t558 + t417 * t561 - t520 * t557;
t543 = t418 * t558 + t419 * t561 - t519 * t557;
t521 = t479 * t481;
t522 = t478 * t481;
t448 = t480 * t521 - t522 * t526;
t449 = (t478 * t480 + t479 * t526) * t481;
t542 = t448 * t555 + t449 * t558 - t483 * t552;
t541 = -t448 * t552 + t449 * t557 + t483 * t554;
t540 = t448 * t558 + t449 * t561 + t483 * t557;
t539 = -t453 * t556 + t454 * t560 + t520 * t553;
t538 = t455 * t556 - t456 * t560 - t519 * t553;
t537 = t453 * t553 - t454 * t559 + t520 * t551;
t536 = -t455 * t553 + t456 * t559 - t519 * t551;
t535 = t560 * t453 - t454 * t562 - t559 * t520;
t534 = -t560 * t455 + t456 * t562 + t559 * t519;
t533 = t553 * t483 + (t478 * t556 - t479 * t560) * t481;
t532 = t551 * t483 + (-t478 * t553 + t479 * t559) * t481;
t531 = -t559 * t483 + (-t560 * t478 + t479 * t562) * t481;
t524 = Icges(3,4) * t481;
t523 = Icges(3,4) * t483;
t516 = -rSges(7,2) * t520 + t549 * t416 + t550 * t417;
t515 = -rSges(7,2) * t519 + t549 * t418 + t550 * t419;
t514 = rSges(7,2) * t483 + t549 * t448 + t550 * t449;
t508 = qJD(3) * t481;
t475 = t484 * t508;
t513 = qJD(4) * t455 + t475;
t467 = pkin(2) * t481 - qJ(3) * t483;
t512 = -(pkin(3) * t479 + qJ(4) * t478) * t481 - t467;
t497 = pkin(2) * t483 + qJ(3) * t481;
t458 = t497 * t482;
t471 = pkin(1) * t482 - pkin(7) * t484;
t511 = -t458 - t471;
t510 = qJD(2) * t482;
t509 = qJD(2) * t484;
t507 = qJD(5) * t481;
t421 = pkin(3) * t454 + qJ(4) * t453;
t506 = -t421 + t511;
t459 = t497 * t484;
t463 = qJD(1) * (pkin(1) * t484 + pkin(7) * t482);
t505 = qJD(1) * t459 + t482 * t508 + t463;
t502 = qJD(2) * (rSges(4,3) * t483 - (rSges(4,1) * t479 - rSges(4,2) * t478) * t481 - t467);
t501 = qJD(2) * (rSges(5,2) * t483 - (rSges(5,1) * t479 + rSges(5,3) * t478) * t481 + t512);
t500 = qJD(2) * (-pkin(4) * t521 - pkin(8) * t483 + t512);
t499 = -qJD(3) * t483 + t458 * t510 + t459 * t509;
t498 = rSges(3,1) * t483 - rSges(3,2) * t481;
t422 = pkin(3) * t456 + qJ(4) * t455;
t496 = qJD(1) * t422 + qJD(4) * t453 + t505;
t495 = Icges(3,1) * t483 - t524;
t494 = -Icges(3,2) * t481 + t523;
t493 = Icges(3,5) * t483 - Icges(3,6) * t481;
t438 = -Icges(3,6) * t484 + t482 * t494;
t440 = -Icges(3,5) * t484 + t482 * t495;
t492 = t438 * t481 - t440 * t483;
t439 = Icges(3,6) * t482 + t484 * t494;
t441 = Icges(3,5) * t482 + t484 * t495;
t491 = -t439 * t481 + t441 * t483;
t465 = Icges(3,2) * t483 + t524;
t466 = Icges(3,1) * t481 + t523;
t490 = -t465 * t481 + t466 * t483;
t489 = qJD(4) * t522 + t421 * t510 + t422 * t509 + t499;
t426 = pkin(4) * t454 - pkin(8) * t520;
t427 = pkin(4) * t456 - pkin(8) * t519;
t488 = t426 * t510 + t427 * t509 + t489;
t487 = qJD(1) * t427 + t482 * t500 + t496;
t486 = (-t426 + t506) * qJD(1) + t484 * t500 + t513;
t476 = qJD(5) * t483 + qJD(1);
t470 = rSges(2,1) * t484 - rSges(2,2) * t482;
t469 = rSges(2,1) * t482 + rSges(2,2) * t484;
t468 = rSges(3,1) * t481 + rSges(3,2) * t483;
t464 = Icges(3,5) * t481 + Icges(3,6) * t483;
t461 = -t482 * t507 - t509;
t460 = -t484 * t507 + t510;
t447 = rSges(3,3) * t482 + t484 * t498;
t446 = -rSges(3,3) * t484 + t482 * t498;
t437 = Icges(3,3) * t482 + t484 * t493;
t436 = -Icges(3,3) * t484 + t482 * t493;
t412 = rSges(4,1) * t456 - rSges(4,2) * t455 + rSges(4,3) * t519;
t411 = rSges(5,1) * t456 + rSges(5,2) * t519 + rSges(5,3) * t455;
t410 = rSges(4,1) * t454 - rSges(4,2) * t453 + rSges(4,3) * t520;
t409 = rSges(5,1) * t454 + rSges(5,2) * t520 + rSges(5,3) * t453;
t396 = rSges(6,1) * t449 - rSges(6,2) * t448 + rSges(6,3) * t483;
t388 = qJD(1) * t447 - t468 * t510 + t463;
t387 = -t468 * t509 + (-t446 - t471) * qJD(1);
t386 = (t446 * t482 + t447 * t484) * qJD(2);
t383 = rSges(6,1) * t419 - rSges(6,2) * t418 - rSges(6,3) * t519;
t381 = rSges(6,1) * t417 - rSges(6,2) * t416 - rSges(6,3) * t520;
t367 = qJD(1) * t412 + t482 * t502 + t505;
t366 = t475 + t484 * t502 + (-t410 + t511) * qJD(1);
t365 = (t410 * t482 + t412 * t484) * qJD(2) + t499;
t364 = qJD(1) * t411 + t482 * t501 + t496;
t363 = t484 * t501 + (-t409 + t506) * qJD(1) + t513;
t362 = (t409 * t482 + t411 * t484) * qJD(2) + t489;
t361 = t383 * t476 - t396 * t460 + t487;
t360 = -t381 * t476 + t396 * t461 + t486;
t359 = t381 * t460 - t383 * t461 + t488;
t358 = qJD(6) * t416 - t460 * t514 + t476 * t515 + t487;
t357 = qJD(6) * t418 + t461 * t514 - t476 * t516 + t486;
t356 = qJD(6) * t448 + t460 * t516 - t461 * t515 + t488;
t1 = m(4) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(5) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(6) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(7) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(3) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + ((t542 * t418 + t540 * t419 - t541 * t519) * t476 + (t548 * t418 + t544 * t419 - t546 * t519) * t461 + (t547 * t418 + t543 * t419 - t545 * t519) * t460) * t460 / 0.2e1 + ((t542 * t416 + t540 * t417 - t541 * t520) * t476 + (t548 * t416 + t544 * t417 - t546 * t520) * t461 + (t547 * t416 + t543 * t417 - t545 * t520) * t460) * t461 / 0.2e1 + ((t542 * t448 + t540 * t449 + t541 * t483) * t476 + (t548 * t448 + t544 * t449 + t546 * t483) * t461 + (t547 * t448 + t543 * t449 + t545 * t483) * t460) * t476 / 0.2e1 + (m(2) * (t469 ^ 2 + t470 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((((-t438 - t537) * t484 + (t439 - t536) * t482) * t483 + ((t539 * t478 + t535 * t479 - t440) * t484 + (t538 * t478 + t534 * t479 + t441) * t482) * t481) * qJD(2) + ((t465 - t532) * t483 + (t533 * t478 + t531 * t479 + t466) * t481) * qJD(1)) * qJD(1) / 0.2e1 + (((t539 * t455 + t535 * t456 + t492 * t484 + t537 * t519) * t484 + ((-t436 + t491) * t484 + t437 * t482 + t536 * t519 + t534 * t456 + t538 * t455) * t482) * qJD(2) + (t533 * t455 + t531 * t456 + t482 * t464 + t484 * t490 + t532 * t519) * qJD(1)) * t510 / 0.2e1 - (((t436 * t484 + t539 * t453 + t535 * t454 + t537 * t520) * t484 + (t491 * t482 + (-t437 + t492) * t484 + t536 * t520 + t534 * t454 + t538 * t453) * t482) * qJD(2) + (t533 * t453 + t531 * t454 - t484 * t464 + t482 * t490 + t532 * t520) * qJD(1)) * t509 / 0.2e1;
T  = t1;
