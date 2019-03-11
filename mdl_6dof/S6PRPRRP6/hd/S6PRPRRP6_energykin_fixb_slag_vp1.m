% Calculate kinetic energy for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:11
% EndTime: 2019-03-08 20:18:14
% DurationCPUTime: 2.71s
% Computational Cost: add. (1989->275), mult. (5108->418), div. (0->0), fcn. (6155->10), ass. (0->136)
t586 = Icges(3,1) + Icges(4,2);
t585 = Icges(6,1) + Icges(7,1);
t584 = Icges(3,4) + Icges(4,6);
t583 = -Icges(6,4) + Icges(7,5);
t582 = Icges(7,4) + Icges(6,5);
t581 = Icges(3,5) - Icges(4,4);
t580 = Icges(3,2) + Icges(4,3);
t579 = Icges(6,2) + Icges(7,3);
t578 = Icges(7,2) + Icges(6,3);
t577 = Icges(3,6) - Icges(4,5);
t576 = Icges(6,6) - Icges(7,6);
t575 = Icges(3,3) + Icges(4,1);
t574 = rSges(7,1) + pkin(5);
t573 = rSges(7,3) + qJ(6);
t515 = cos(pkin(6));
t513 = sin(pkin(6));
t514 = cos(pkin(10));
t543 = t514 * t513;
t512 = sin(pkin(10));
t546 = t512 * t513;
t518 = sin(qJ(2));
t519 = cos(qJ(2));
t556 = t575 * t515 + (t581 * t518 + t577 * t519) * t513;
t541 = t515 * t519;
t497 = t512 * t518 - t514 * t541;
t542 = t515 * t518;
t498 = t512 * t519 + t514 * t542;
t557 = -t577 * t497 + t581 * t498 - t575 * t543;
t499 = t512 * t541 + t514 * t518;
t500 = -t512 * t542 + t514 * t519;
t558 = -t577 * t499 + t581 * t500 + t575 * t546;
t572 = -t556 * t515 + t543 * t557 - t546 * t558;
t517 = sin(qJ(4));
t548 = cos(qJ(4));
t533 = t513 * t548;
t473 = t499 * t517 + t512 * t533;
t516 = sin(qJ(5));
t547 = cos(qJ(5));
t438 = t473 * t516 - t500 * t547;
t439 = t473 * t547 + t500 * t516;
t545 = t513 * t517;
t472 = -t499 * t548 + t512 * t545;
t571 = t579 * t438 + t583 * t439 - t576 * t472;
t475 = t497 * t517 - t514 * t533;
t440 = t475 * t516 - t498 * t547;
t441 = t475 * t547 + t498 * t516;
t474 = t497 * t548 + t517 * t543;
t570 = t579 * t440 + t583 * t441 + t576 * t474;
t569 = -t576 * t438 + t582 * t439 + t578 * t472;
t568 = -t576 * t440 + t582 * t441 - t578 * t474;
t567 = t583 * t438 + t585 * t439 + t582 * t472;
t566 = t583 * t440 + t585 * t441 - t582 * t474;
t502 = t515 * t548 - t519 * t545;
t544 = t513 * t518;
t476 = t502 * t516 - t544 * t547;
t477 = t502 * t547 + t516 * t544;
t501 = t515 * t517 + t519 * t533;
t565 = t579 * t476 + t583 * t477 - t576 * t501;
t564 = -t576 * t476 + t582 * t477 + t578 * t501;
t563 = t583 * t476 + t585 * t477 + t582 * t501;
t562 = t580 * t499 - t584 * t500 - t577 * t546;
t561 = t580 * t497 - t584 * t498 + t577 * t543;
t560 = -t584 * t499 + t586 * t500 + t581 * t546;
t559 = t584 * t497 - t586 * t498 + t581 * t543;
t555 = t577 * t515 + (t584 * t518 + t580 * t519) * t513;
t554 = t581 * t515 + (t586 * t518 + t584 * t519) * t513;
t552 = qJD(2) ^ 2;
t540 = rSges(7,2) * t472 + t573 * t438 + t439 * t574;
t539 = -rSges(7,2) * t474 + t573 * t440 + t441 * t574;
t538 = rSges(7,2) * t501 + t573 * t476 + t477 * t574;
t466 = pkin(2) * t500 + qJ(3) * t499;
t511 = qJD(2) * t515;
t537 = qJD(3) * t497 + t466 * t511;
t536 = qJD(2) * t513;
t509 = t512 * t536;
t478 = qJD(4) * t500 + t509;
t504 = qJD(4) * t544 + t511;
t535 = qJD(3) * t519;
t532 = t514 * t536;
t465 = pkin(2) * t498 + qJ(3) * t497;
t531 = t465 * t509 + t466 * t532 + qJD(1);
t503 = (pkin(2) * t518 - qJ(3) * t519) * t513;
t529 = (-t515 * rSges(4,1) - (-rSges(4,2) * t518 - rSges(4,3) * t519) * t513 - t503) * t513;
t528 = (-pkin(3) * t515 - pkin(8) * t544 - t503) * t513;
t479 = qJD(4) * t498 - t532;
t480 = pkin(3) * t546 + pkin(8) * t500;
t481 = -pkin(3) * t543 + pkin(8) * t498;
t525 = t480 * t532 + t481 * t509 - t513 * t535 + t531;
t524 = qJD(2) * t512 * t528 + t480 * t511 + t537;
t433 = pkin(4) * t473 + pkin(9) * t472;
t434 = pkin(4) * t475 - pkin(9) * t474;
t523 = -t479 * t433 + t478 * t434 + t525;
t496 = qJD(3) * t499;
t522 = t496 + ((-t465 - t481) * t515 + t514 * t528) * qJD(2);
t467 = pkin(4) * t502 + pkin(9) * t501;
t521 = t504 * t433 - t467 * t478 + t524;
t520 = -t434 * t504 + t479 * t467 + t522;
t488 = t515 * rSges(3,3) + (rSges(3,1) * t518 + rSges(3,2) * t519) * t513;
t470 = qJD(5) * t501 + t504;
t463 = rSges(5,1) * t502 - rSges(5,2) * t501 + rSges(5,3) * t544;
t462 = Icges(5,1) * t502 - Icges(5,4) * t501 + Icges(5,5) * t544;
t461 = Icges(5,4) * t502 - Icges(5,2) * t501 + Icges(5,6) * t544;
t460 = Icges(5,5) * t502 - Icges(5,6) * t501 + Icges(5,3) * t544;
t459 = rSges(3,1) * t500 - rSges(3,2) * t499 + rSges(3,3) * t546;
t458 = rSges(3,1) * t498 - rSges(3,2) * t497 - rSges(3,3) * t543;
t457 = -rSges(4,1) * t543 - rSges(4,2) * t498 + rSges(4,3) * t497;
t456 = rSges(4,1) * t546 - rSges(4,2) * t500 + rSges(4,3) * t499;
t437 = -qJD(5) * t474 + t479;
t436 = qJD(5) * t472 + t478;
t430 = (-t458 * t515 - t488 * t543) * qJD(2);
t429 = (t459 * t515 - t488 * t546) * qJD(2);
t428 = rSges(6,1) * t477 - rSges(6,2) * t476 + rSges(6,3) * t501;
t420 = rSges(5,1) * t475 + rSges(5,2) * t474 + rSges(5,3) * t498;
t419 = rSges(5,1) * t473 - rSges(5,2) * t472 + rSges(5,3) * t500;
t418 = Icges(5,1) * t475 + Icges(5,4) * t474 + Icges(5,5) * t498;
t417 = Icges(5,1) * t473 - Icges(5,4) * t472 + Icges(5,5) * t500;
t416 = Icges(5,4) * t475 + Icges(5,2) * t474 + Icges(5,6) * t498;
t415 = Icges(5,4) * t473 - Icges(5,2) * t472 + Icges(5,6) * t500;
t414 = Icges(5,5) * t475 + Icges(5,6) * t474 + Icges(5,3) * t498;
t413 = Icges(5,5) * t473 - Icges(5,6) * t472 + Icges(5,3) * t500;
t411 = qJD(1) + (t458 * t512 + t459 * t514) * t536;
t408 = rSges(6,1) * t441 - rSges(6,2) * t440 - rSges(6,3) * t474;
t406 = rSges(6,1) * t439 - rSges(6,2) * t438 + rSges(6,3) * t472;
t392 = t496 + ((-t457 - t465) * t515 + t514 * t529) * qJD(2);
t391 = (t456 * t515 + t512 * t529) * qJD(2) + t537;
t390 = (-t535 + (t456 * t514 + t457 * t512) * qJD(2)) * t513 + t531;
t389 = -t420 * t504 + t463 * t479 + t522;
t388 = t419 * t504 - t463 * t478 + t524;
t387 = -t479 * t419 + t478 * t420 + t525;
t386 = -t408 * t470 + t428 * t437 + t520;
t385 = t406 * t470 - t428 * t436 + t521;
t384 = -t437 * t406 + t436 * t408 + t523;
t383 = qJD(6) * t438 + t437 * t538 - t470 * t539 + t520;
t382 = qJD(6) * t440 - t436 * t538 + t470 * t540 + t521;
t381 = qJD(6) * t476 + t436 * t539 - t437 * t540 + t523;
t1 = m(6) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(7) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(5) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(4) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(3) * (t411 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t478 * ((t500 * t413 - t472 * t415 + t473 * t417) * t478 + (t414 * t500 - t416 * t472 + t418 * t473) * t479 + (t460 * t500 - t461 * t472 + t462 * t473) * t504) / 0.2e1 + t479 * ((t413 * t498 + t415 * t474 + t417 * t475) * t478 + (t498 * t414 + t474 * t416 + t475 * t418) * t479 + (t460 * t498 + t461 * t474 + t462 * t475) * t504) / 0.2e1 + t504 * ((t413 * t544 - t415 * t501 + t417 * t502) * t478 + (t414 * t544 - t416 * t501 + t418 * t502) * t479 + (t460 * t544 - t461 * t501 + t462 * t502) * t504) / 0.2e1 + ((t438 * t565 + t439 * t563 + t472 * t564) * t470 + (t438 * t570 + t439 * t566 + t472 * t568) * t437 + (t571 * t438 + t567 * t439 + t569 * t472) * t436) * t436 / 0.2e1 + ((t440 * t565 + t441 * t563 - t474 * t564) * t470 + (t570 * t440 + t566 * t441 - t568 * t474) * t437 + (t440 * t571 + t567 * t441 - t569 * t474) * t436) * t437 / 0.2e1 + ((t565 * t476 + t563 * t477 + t564 * t501) * t470 + (t476 * t570 + t477 * t566 + t501 * t568) * t437 + (t476 * t571 + t567 * t477 + t569 * t501) * t436) * t470 / 0.2e1 - ((t497 * t562 + t498 * t560) * t546 + (-t497 * t555 + t498 * t554) * t515 + (-t497 * t561 + t498 * t559 + t572) * t543) * t552 * t543 / 0.2e1 + ((t556 * t515 ^ 2 + (((t518 * t559 + t519 * t561) * t514 + (t518 * t560 - t562 * t519) * t512) * t513 + (t512 * t558 - t557 * t514 + t518 * t554 + t519 * t555) * t515) * t513) * t515 + ((-t499 * t561 + t500 * t559) * t543 + (-t499 * t555 + t500 * t554) * t515 + (t499 * t562 + t500 * t560 - t572) * t546) * t546) * t552 / 0.2e1;
T  = t1;
