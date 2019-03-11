% Calculate kinetic energy for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:49
% EndTime: 2019-03-09 09:58:52
% DurationCPUTime: 3.05s
% Computational Cost: add. (1219->257), mult. (2115->371), div. (0->0), fcn. (2082->8), ass. (0->139)
t582 = Icges(3,4) + Icges(4,6);
t581 = Icges(3,1) + Icges(4,2);
t580 = -Icges(3,2) - Icges(4,3);
t481 = cos(qJ(2));
t579 = t582 * t481;
t478 = sin(qJ(2));
t578 = t582 * t478;
t577 = -Icges(4,4) + Icges(3,5);
t576 = Icges(4,5) - Icges(3,6);
t575 = t580 * t478 + t579;
t574 = -t581 * t481 + t578;
t573 = Icges(4,1) + Icges(3,3);
t572 = Icges(6,1) + Icges(7,1);
t571 = Icges(6,4) - Icges(7,5);
t570 = Icges(7,4) + Icges(6,5);
t569 = Icges(6,2) + Icges(7,3);
t568 = Icges(7,6) - Icges(6,6);
t479 = sin(qJ(1));
t482 = cos(qJ(1));
t567 = t575 * t479 + t576 * t482;
t566 = -t576 * t479 + t575 * t482;
t565 = t574 * t479 + t577 * t482;
t564 = t577 * t479 - t574 * t482;
t563 = t580 * t481 - t578;
t562 = t581 * t478 + t579;
t561 = t576 * t478 + t577 * t481;
t560 = Icges(6,3) + Icges(7,2) + Icges(5,3);
t559 = rSges(7,1) + pkin(5);
t558 = rSges(7,3) + qJ(6);
t477 = sin(qJ(4));
t534 = pkin(4) * t477;
t488 = qJ(5) * t481 + t478 * t534;
t480 = cos(qJ(4));
t532 = pkin(4) * t480;
t400 = t479 * t488 - t482 * t532;
t512 = qJD(4) * t481;
t515 = qJD(2) * t479;
t448 = t482 * t512 + t515;
t557 = qJD(5) * t478 + t448 * t400;
t475 = qJ(4) + pkin(9);
t472 = sin(t475);
t473 = cos(t475);
t525 = t478 * t482;
t410 = t472 * t479 - t473 * t525;
t411 = t472 * t525 + t473 * t479;
t521 = t481 * t482;
t556 = t569 * t410 - t571 * t411 + t568 * t521;
t526 = t478 * t479;
t412 = t472 * t482 + t473 * t526;
t413 = t472 * t526 - t473 * t482;
t522 = t481 * t479;
t555 = -t569 * t412 - t571 * t413 + t568 * t522;
t554 = -t571 * t410 + t572 * t411 + t570 * t521;
t553 = t571 * t412 + t572 * t413 + t570 * t522;
t552 = (t571 * t472 + t569 * t473) * t481 + t568 * t478;
t551 = (-t572 * t472 - t571 * t473) * t481 + t570 * t478;
t550 = t561 * t479 - t573 * t482;
t549 = t573 * t479 + t561 * t482;
t450 = pkin(3) * t479 + pkin(8) * t521;
t451 = -pkin(3) * t482 + pkin(8) * t522;
t514 = qJD(2) * t482;
t548 = t450 * t514 + t451 * t515;
t547 = t577 * t478 - t576 * t481;
t546 = t563 * t478 + t562 * t481;
t545 = -t566 * t478 + t564 * t481;
t544 = t567 * t478 + t565 * t481;
t523 = t480 * t482;
t441 = -t477 * t479 + t478 * t523;
t524 = t479 * t480;
t442 = t477 * t525 + t524;
t543 = Icges(5,5) * t442 + Icges(5,6) * t441 + t568 * t410 + t570 * t411 + t560 * t521;
t443 = t477 * t482 + t478 * t524;
t444 = t477 * t526 - t523;
t542 = Icges(5,5) * t444 + Icges(5,6) * t443 - t568 * t412 + t570 * t413 + t560 * t522;
t541 = (-Icges(5,5) * t477 - Icges(5,6) * t480 - t570 * t472 + t568 * t473) * t481 + t560 * t478;
t520 = rSges(7,2) * t521 + t558 * t410 + t559 * t411;
t519 = rSges(7,2) * t522 - t558 * t412 + t559 * t413;
t518 = rSges(7,2) * t478 + (-t559 * t472 + t558 * t473) * t481;
t503 = pkin(2) * t481 + qJ(3) * t478;
t445 = t503 * t479;
t446 = t503 * t482;
t517 = t445 * t515 + t446 * t514;
t465 = pkin(1) * t479 - pkin(7) * t482;
t516 = -t445 - t465;
t513 = qJD(3) * t478;
t511 = qJD(5) * t481;
t452 = qJD(1) * (pkin(1) * t482 + pkin(7) * t479);
t510 = qJD(1) * t446 + t479 * t513 + t452;
t460 = pkin(2) * t478 - qJ(3) * t481;
t507 = qJD(2) * (rSges(4,2) * t478 + rSges(4,3) * t481 - t460);
t506 = -qJD(3) * t481 + t517;
t505 = rSges(3,1) * t481 - rSges(3,2) * t478;
t504 = -rSges(4,2) * t481 + rSges(4,3) * t478;
t502 = qJD(2) * (-pkin(8) * t478 - t460);
t489 = t506 + t548;
t487 = qJD(1) * t450 + t479 * t502 + t510;
t399 = t479 * t532 + t482 * t488;
t470 = qJD(4) * t478 + qJD(1);
t486 = t470 * t399 + t479 * t511 + t487;
t469 = t482 * t513;
t485 = t469 + (-t451 + t516) * qJD(1) + t482 * t502;
t437 = qJ(5) * t478 - t481 * t534;
t449 = t479 * t512 - t514;
t484 = t449 * t437 + t482 * t511 + t485;
t464 = rSges(2,1) * t482 - rSges(2,2) * t479;
t463 = rSges(2,1) * t479 + rSges(2,2) * t482;
t462 = rSges(3,1) * t478 + rSges(3,2) * t481;
t436 = -rSges(4,1) * t482 + t479 * t504;
t435 = rSges(4,1) * t479 + t482 * t504;
t434 = rSges(3,3) * t479 + t482 * t505;
t433 = rSges(5,3) * t478 + (-rSges(5,1) * t477 - rSges(5,2) * t480) * t481;
t432 = -rSges(3,3) * t482 + t479 * t505;
t420 = Icges(5,5) * t478 + (-Icges(5,1) * t477 - Icges(5,4) * t480) * t481;
t417 = Icges(5,6) * t478 + (-Icges(5,4) * t477 - Icges(5,2) * t480) * t481;
t409 = rSges(6,3) * t478 + (-rSges(6,1) * t472 - rSges(6,2) * t473) * t481;
t398 = rSges(5,1) * t444 + rSges(5,2) * t443 + rSges(5,3) * t522;
t397 = rSges(5,1) * t442 + rSges(5,2) * t441 + rSges(5,3) * t521;
t396 = Icges(5,1) * t444 + Icges(5,4) * t443 + Icges(5,5) * t522;
t395 = Icges(5,1) * t442 + Icges(5,4) * t441 + Icges(5,5) * t521;
t394 = Icges(5,4) * t444 + Icges(5,2) * t443 + Icges(5,6) * t522;
t393 = Icges(5,4) * t442 + Icges(5,2) * t441 + Icges(5,6) * t521;
t387 = qJD(1) * t434 - t462 * t515 + t452;
t386 = -t462 * t514 + (-t432 - t465) * qJD(1);
t385 = (t432 * t479 + t434 * t482) * qJD(2);
t384 = rSges(6,1) * t413 + rSges(6,2) * t412 + rSges(6,3) * t522;
t382 = rSges(6,1) * t411 - rSges(6,2) * t410 + rSges(6,3) * t521;
t367 = qJD(1) * t435 + t479 * t507 + t510;
t366 = t469 + t482 * t507 + (-t436 + t516) * qJD(1);
t365 = (t435 * t482 + t436 * t479) * qJD(2) + t506;
t364 = t397 * t470 - t433 * t448 + t487;
t363 = -t398 * t470 + t433 * t449 + t485;
t362 = -t397 * t449 + t398 * t448 + t489;
t361 = t382 * t470 + (-t409 - t437) * t448 + t486;
t360 = t409 * t449 + (-t384 - t400) * t470 + t484;
t359 = t384 * t448 + (-t382 - t399) * t449 + t489 + t557;
t358 = -qJD(6) * t412 + t520 * t470 + (-t437 - t518) * t448 + t486;
t357 = qJD(6) * t410 + t518 * t449 + (-t400 - t519) * t470 + t484;
t356 = (qJD(6) * t473 - qJD(3)) * t481 + t519 * t448 + (-t399 - t520) * t449 + t517 + t548 + t557;
t1 = m(7) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(3) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(4) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(5) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(6) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t463 ^ 2 + t464 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t565 * t478 - t567 * t481) * t482 + (t564 * t478 + t566 * t481) * t479) * qJD(2) + (t562 * t478 - t563 * t481) * qJD(1)) * qJD(1) / 0.2e1 + ((t549 * t479 ^ 2 + (t544 * t482 + (t545 - t550) * t479) * t482) * qJD(2) + (t479 * t547 + t482 * t546) * qJD(1)) * t515 / 0.2e1 - ((t550 * t482 ^ 2 + (t545 * t479 + (t544 - t549) * t482) * t479) * qJD(2) + (t479 * t546 - t482 * t547) * qJD(1)) * t514 / 0.2e1 + ((t410 * t552 + t411 * t551 + t417 * t441 + t420 * t442 + t521 * t541) * t470 + (t394 * t441 + t396 * t442 + t410 * t555 + t411 * t553 + t521 * t542) * t449 + (t393 * t441 + t395 * t442 + t556 * t410 + t554 * t411 + t543 * t521) * t448) * t448 / 0.2e1 + ((-t412 * t552 + t413 * t551 + t417 * t443 + t420 * t444 + t522 * t541) * t470 + (t394 * t443 + t396 * t444 - t555 * t412 + t553 * t413 + t542 * t522) * t449 + (t393 * t443 + t395 * t444 - t412 * t556 + t413 * t554 + t522 * t543) * t448) * t449 / 0.2e1 + (((-t417 * t480 - t420 * t477 - t472 * t551 + t473 * t552) * t470 + (-t394 * t480 - t396 * t477 - t472 * t553 + t473 * t555) * t449 + (-t393 * t480 - t395 * t477 - t472 * t554 + t473 * t556) * t448) * t481 + (t448 * t543 + t449 * t542 + t470 * t541) * t478) * t470 / 0.2e1;
T  = t1;
