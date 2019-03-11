% Calculate kinetic energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:24
% EndTime: 2019-03-09 08:41:27
% DurationCPUTime: 3.48s
% Computational Cost: add. (1192->264), mult. (2070->383), div. (0->0), fcn. (2037->8), ass. (0->141)
t583 = Icges(3,4) + Icges(4,6);
t582 = Icges(3,1) + Icges(4,2);
t581 = Icges(3,2) + Icges(4,3);
t481 = cos(qJ(2));
t580 = t583 * t481;
t479 = sin(qJ(2));
t579 = t583 * t479;
t578 = -Icges(4,4) + Icges(3,5);
t577 = Icges(4,5) - Icges(3,6);
t576 = t479 * t581 - t580;
t575 = -t481 * t582 + t579;
t574 = Icges(4,1) + Icges(3,3);
t573 = Icges(6,1) + Icges(7,1);
t572 = Icges(6,4) - Icges(7,5);
t571 = Icges(7,4) + Icges(6,5);
t570 = Icges(6,2) + Icges(7,3);
t569 = Icges(7,6) - Icges(6,6);
t568 = Icges(6,3) + Icges(7,2);
t480 = sin(qJ(1));
t482 = cos(qJ(1));
t567 = -t576 * t480 + t577 * t482;
t566 = t577 * t480 + t576 * t482;
t565 = t575 * t480 + t578 * t482;
t564 = t578 * t480 - t575 * t482;
t563 = -t481 * t581 - t579;
t562 = t479 * t582 + t580;
t561 = t577 * t479 + t578 * t481;
t560 = rSges(7,1) + pkin(5);
t559 = rSges(7,3) + qJ(6);
t475 = pkin(9) + qJ(5);
t472 = sin(t475);
t473 = cos(t475);
t528 = t479 * t482;
t414 = t472 * t480 - t473 * t528;
t415 = t472 * t528 + t473 * t480;
t526 = t481 * t482;
t558 = t570 * t414 - t572 * t415 + t569 * t526;
t529 = t479 * t480;
t416 = t472 * t482 + t473 * t529;
t417 = t472 * t529 - t473 * t482;
t527 = t480 * t481;
t557 = -t570 * t416 - t572 * t417 + t569 * t527;
t556 = t569 * t414 + t571 * t415 + t568 * t526;
t555 = -t569 * t416 + t571 * t417 + t568 * t527;
t554 = -t572 * t414 + t573 * t415 + t571 * t526;
t553 = t572 * t416 + t573 * t417 + t571 * t527;
t476 = sin(pkin(9));
t536 = pkin(4) * t476;
t487 = pkin(8) * t481 + t479 * t536;
t477 = cos(pkin(9));
t534 = pkin(4) * t477;
t400 = t480 * t534 + t482 * t487;
t401 = t480 * t487 - t482 * t534;
t517 = qJD(2) * t482;
t518 = qJD(2) * t480;
t552 = t400 * t517 + t401 * t518;
t551 = (t572 * t472 + t570 * t473) * t481 + t569 * t479;
t550 = (-t571 * t472 + t569 * t473) * t481 + t568 * t479;
t549 = (-t573 * t472 - t572 * t473) * t481 + t571 * t479;
t548 = t561 * t480 - t574 * t482;
t547 = t574 * t480 + t561 * t482;
t546 = t578 * t479 - t577 * t481;
t545 = t563 * t479 + t562 * t481;
t544 = t566 * t479 + t564 * t481;
t543 = t567 * t479 + t565 * t481;
t524 = rSges(7,2) * t526 + t559 * t414 + t560 * t415;
t523 = rSges(7,2) * t527 - t559 * t416 + t560 * t417;
t522 = rSges(7,2) * t479 + (-t560 * t472 + t559 * t473) * t481;
t504 = pkin(2) * t481 + qJ(3) * t479;
t445 = t504 * t480;
t446 = t504 * t482;
t521 = t445 * t518 + t446 * t517;
t465 = pkin(1) * t480 - pkin(7) * t482;
t520 = -t445 - t465;
t516 = qJD(3) * t479;
t469 = t482 * t516;
t514 = qJD(4) * t481;
t519 = t482 * t514 + t469;
t515 = qJD(3) * t481;
t513 = qJD(5) * t481;
t452 = qJD(1) * (pkin(1) * t482 + pkin(7) * t480);
t512 = qJD(1) * t446 + t480 * t516 + t452;
t451 = -pkin(3) * t482 + qJ(4) * t527;
t511 = -t451 + t520;
t460 = pkin(2) * t479 - qJ(3) * t481;
t508 = -qJ(4) * t479 - t460;
t507 = qJD(2) * (rSges(4,2) * t479 + rSges(4,3) * t481 - t460);
t506 = rSges(3,1) * t481 - rSges(3,2) * t479;
t505 = -rSges(4,2) * t481 + rSges(4,3) * t479;
t450 = pkin(3) * t480 + qJ(4) * t526;
t503 = qJD(4) * t479 + t450 * t517 + t451 * t518 + t521;
t502 = qJD(1) * t450 + t480 * t514 + t512;
t489 = qJD(2) * (-rSges(5,3) * t479 - (-rSges(5,1) * t476 - rSges(5,2) * t477) * t481 + t508);
t488 = qJD(2) * (-pkin(8) * t479 + t481 * t536 + t508);
t486 = t503 - t515;
t485 = qJD(1) * t400 + t480 * t488 + t502;
t484 = (-t401 + t511) * qJD(1) + t482 * t488 + t519;
t470 = qJD(5) * t479 + qJD(1);
t464 = rSges(2,1) * t482 - rSges(2,2) * t480;
t463 = rSges(2,1) * t480 + rSges(2,2) * t482;
t462 = rSges(3,1) * t479 + rSges(3,2) * t481;
t449 = t480 * t513 - t517;
t448 = t482 * t513 + t518;
t444 = t476 * t529 - t477 * t482;
t443 = t476 * t482 + t477 * t529;
t442 = t476 * t528 + t477 * t480;
t441 = -t476 * t480 + t477 * t528;
t436 = -rSges(4,1) * t482 + t480 * t505;
t435 = rSges(4,1) * t480 + t482 * t505;
t434 = rSges(3,3) * t480 + t482 * t506;
t433 = -rSges(3,3) * t482 + t480 * t506;
t412 = Icges(5,5) * t479 + (-Icges(5,1) * t476 - Icges(5,4) * t477) * t481;
t411 = Icges(5,6) * t479 + (-Icges(5,4) * t476 - Icges(5,2) * t477) * t481;
t410 = Icges(5,3) * t479 + (-Icges(5,5) * t476 - Icges(5,6) * t477) * t481;
t409 = rSges(6,3) * t479 + (-rSges(6,1) * t472 - rSges(6,2) * t473) * t481;
t398 = rSges(5,1) * t444 + rSges(5,2) * t443 + rSges(5,3) * t527;
t397 = rSges(5,1) * t442 + rSges(5,2) * t441 + rSges(5,3) * t526;
t396 = Icges(5,1) * t444 + Icges(5,4) * t443 + Icges(5,5) * t527;
t395 = Icges(5,1) * t442 + Icges(5,4) * t441 + Icges(5,5) * t526;
t394 = Icges(5,4) * t444 + Icges(5,2) * t443 + Icges(5,6) * t527;
t393 = Icges(5,4) * t442 + Icges(5,2) * t441 + Icges(5,6) * t526;
t392 = Icges(5,5) * t444 + Icges(5,6) * t443 + Icges(5,3) * t527;
t391 = Icges(5,5) * t442 + Icges(5,6) * t441 + Icges(5,3) * t526;
t386 = qJD(1) * t434 - t462 * t518 + t452;
t385 = -t462 * t517 + (-t433 - t465) * qJD(1);
t384 = (t433 * t480 + t434 * t482) * qJD(2);
t383 = rSges(6,1) * t417 + rSges(6,2) * t416 + rSges(6,3) * t527;
t381 = rSges(6,1) * t415 - rSges(6,2) * t414 + rSges(6,3) * t526;
t367 = qJD(1) * t435 + t480 * t507 + t512;
t366 = t469 + t482 * t507 + (-t436 + t520) * qJD(1);
t365 = -t515 + (t435 * t482 + t436 * t480) * qJD(2) + t521;
t364 = qJD(1) * t397 + t480 * t489 + t502;
t363 = t482 * t489 + (-t398 + t511) * qJD(1) + t519;
t362 = (t397 * t482 + t398 * t480) * qJD(2) + t486;
t361 = t381 * t470 - t409 * t448 + t485;
t360 = -t383 * t470 + t409 * t449 + t484;
t359 = -t381 * t449 + t383 * t448 + t486 + t552;
t358 = -qJD(6) * t416 - t448 * t522 + t470 * t524 + t485;
t357 = qJD(6) * t414 + t449 * t522 - t470 * t523 + t484;
t356 = (qJD(6) * t473 - qJD(3)) * t481 - t524 * t449 + t523 * t448 + t503 + t552;
t1 = m(6) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(4) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(5) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(3) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(7) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + ((t551 * t414 + t549 * t415 + t550 * t526) * t470 + (t557 * t414 + t553 * t415 + t555 * t526) * t449 + (t558 * t414 + t554 * t415 + t556 * t526) * t448) * t448 / 0.2e1 + ((-t551 * t416 + t549 * t417 + t550 * t527) * t470 + (-t557 * t416 + t553 * t417 + t555 * t527) * t449 + (-t558 * t416 + t554 * t417 + t556 * t527) * t448) * t449 / 0.2e1 + (((-t549 * t472 + t551 * t473) * t470 + (-t553 * t472 + t557 * t473) * t449 + (-t554 * t472 + t558 * t473) * t448) * t481 + (t556 * t448 + t555 * t449 + t550 * t470) * t479) * t470 / 0.2e1 + (Icges(2,3) + m(2) * (t463 ^ 2 + t464 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((((t394 * t477 + t396 * t476 - t567) * t481 + (-t392 + t565) * t479) * t482 + ((-t393 * t477 - t395 * t476 - t566) * t481 + (t391 + t564) * t479) * t480) * qJD(2) + ((-t411 * t477 - t412 * t476 - t563) * t481 + (t410 + t562) * t479) * qJD(1)) * qJD(1) / 0.2e1 + (((-t392 * t526 - t394 * t441 - t396 * t442 + t543 * t482) * t482 + (t391 * t526 + t393 * t441 + t395 * t442 + (t544 - t548) * t482 + t547 * t480) * t480) * qJD(2) + (t410 * t526 + t411 * t441 + t412 * t442 + t546 * t480 + t545 * t482) * qJD(1)) * t518 / 0.2e1 - (((t391 * t527 + t393 * t443 + t395 * t444 + t544 * t480) * t480 + (-t392 * t527 - t394 * t443 - t396 * t444 + (t543 - t547) * t480 + t548 * t482) * t482) * qJD(2) + (t410 * t527 + t411 * t443 + t412 * t444 + t545 * t480 - t546 * t482) * qJD(1)) * t517 / 0.2e1;
T  = t1;
