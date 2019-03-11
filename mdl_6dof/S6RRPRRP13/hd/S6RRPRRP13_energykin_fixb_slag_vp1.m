% Calculate kinetic energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:33
% EndTime: 2019-03-09 12:55:35
% DurationCPUTime: 2.89s
% Computational Cost: add. (2087->289), mult. (5199->431), div. (0->0), fcn. (6227->10), ass. (0->142)
t588 = Icges(3,1) + Icges(4,2);
t587 = Icges(6,1) + Icges(7,1);
t586 = Icges(3,4) + Icges(4,6);
t585 = Icges(6,4) + Icges(7,4);
t584 = Icges(3,5) - Icges(4,4);
t583 = Icges(6,5) + Icges(7,5);
t582 = Icges(3,2) + Icges(4,3);
t581 = Icges(6,2) + Icges(7,2);
t580 = Icges(3,6) - Icges(4,5);
t579 = Icges(6,6) + Icges(7,6);
t578 = Icges(3,3) + Icges(4,1);
t577 = Icges(6,3) + Icges(7,3);
t576 = rSges(7,3) + qJ(6);
t506 = sin(pkin(6));
t507 = cos(pkin(6));
t514 = cos(qJ(2));
t515 = cos(qJ(1));
t541 = t514 * t515;
t511 = sin(qJ(2));
t512 = sin(qJ(1));
t544 = t511 * t512;
t488 = -t507 * t541 + t544;
t542 = t512 * t514;
t543 = t511 * t515;
t489 = t507 * t543 + t542;
t490 = t507 * t542 + t543;
t491 = -t507 * t544 + t541;
t545 = t506 * t515;
t546 = t506 * t512;
t558 = (-t488 * t580 + t489 * t584 - t545 * t578) * t515 + (t490 * t580 - t491 * t584 - t546 * t578) * t512;
t575 = t558 * t506;
t510 = sin(qJ(4));
t552 = cos(qJ(4));
t531 = t506 * t552;
t466 = t490 * t510 + t512 * t531;
t509 = sin(qJ(5));
t513 = cos(qJ(5));
t429 = -t466 * t509 + t491 * t513;
t548 = t491 * t509;
t430 = t466 * t513 + t548;
t465 = -t490 * t552 + t510 * t546;
t574 = t429 * t579 + t430 * t583 + t465 * t577;
t468 = t488 * t510 - t515 * t531;
t431 = -t468 * t509 + t489 * t513;
t549 = t489 * t509;
t432 = t468 * t513 + t549;
t467 = t488 * t552 + t510 * t545;
t573 = t431 * t579 + t432 * t583 - t467 * t577;
t572 = t429 * t581 + t430 * t585 + t465 * t579;
t571 = t431 * t581 + t432 * t585 - t467 * t579;
t570 = t429 * t585 + t430 * t587 + t465 * t583;
t569 = t431 * t585 + t432 * t587 - t467 * t583;
t487 = -t506 * t514 * t510 + t507 * t552;
t547 = t506 * t511;
t463 = -t487 * t509 + t513 * t547;
t533 = t509 * t547;
t464 = t487 * t513 + t533;
t486 = t507 * t510 + t514 * t531;
t568 = t463 * t579 + t464 * t583 + t486 * t577;
t567 = t463 * t581 + t464 * t585 + t486 * t579;
t566 = t463 * t585 + t464 * t587 + t486 * t583;
t565 = t490 * t582 - t491 * t586 - t546 * t580;
t564 = t488 * t582 - t489 * t586 + t545 * t580;
t563 = -t586 * t490 + t491 * t588 + t584 * t546;
t562 = t586 * t488 - t489 * t588 + t584 * t545;
t561 = t578 * t507 + (t511 * t584 + t514 * t580) * t506;
t560 = t580 * t507 + (t511 * t586 + t514 * t582) * t506;
t559 = t584 * t507 + (t511 * t588 + t586 * t514) * t506;
t551 = pkin(5) * t513;
t540 = rSges(7,1) * t430 + rSges(7,2) * t429 + pkin(5) * t548 + t576 * t465 + t466 * t551;
t539 = rSges(7,1) * t432 + rSges(7,2) * t431 + pkin(5) * t549 - t576 * t467 + t468 * t551;
t538 = rSges(7,1) * t464 + rSges(7,2) * t463 + pkin(5) * t533 + t576 * t486 + t487 * t551;
t457 = pkin(2) * t489 + qJ(3) * t488;
t458 = pkin(2) * t491 + qJ(3) * t490;
t535 = qJD(2) * t506;
t502 = t512 * t535;
t530 = t515 * t535;
t537 = t457 * t502 + t458 * t530;
t469 = qJD(4) * t491 + t502;
t536 = qJD(1) * (pkin(1) * t512 - pkin(8) * t545);
t534 = qJD(3) * t514;
t503 = qJD(2) * t507 + qJD(1);
t494 = qJD(1) * (pkin(1) * t515 + pkin(8) * t546);
t532 = qJD(3) * t488 + t503 * t458 + t494;
t493 = qJD(4) * t547 + t503;
t529 = qJD(3) * t490 - t536;
t492 = (pkin(2) * t511 - qJ(3) * t514) * t506;
t526 = (-rSges(4,1) * t507 - (-rSges(4,2) * t511 - rSges(4,3) * t514) * t506 - t492) * t535;
t525 = (-pkin(3) * t507 - pkin(9) * t547 - t492) * t535;
t470 = qJD(4) * t489 - t530;
t471 = pkin(3) * t546 + pkin(9) * t491;
t472 = -pkin(3) * t545 + pkin(9) * t489;
t522 = t471 * t530 + t472 * t502 - t506 * t534 + t537;
t521 = t503 * t471 + t512 * t525 + t532;
t425 = pkin(4) * t466 + pkin(10) * t465;
t426 = pkin(4) * t468 - pkin(10) * t467;
t520 = -t425 * t470 + t469 * t426 + t522;
t456 = pkin(4) * t487 + pkin(10) * t486;
t519 = t493 * t425 - t456 * t469 + t521;
t518 = (-t457 - t472) * t503 + t515 * t525 + t529;
t517 = -t426 * t493 + t470 * t456 + t518;
t498 = rSges(2,1) * t515 - rSges(2,2) * t512;
t497 = rSges(2,1) * t512 + rSges(2,2) * t515;
t479 = rSges(3,3) * t507 + (rSges(3,1) * t511 + rSges(3,2) * t514) * t506;
t459 = qJD(5) * t486 + t493;
t455 = rSges(3,1) * t491 - rSges(3,2) * t490 + rSges(3,3) * t546;
t454 = rSges(3,1) * t489 - rSges(3,2) * t488 - rSges(3,3) * t545;
t453 = -rSges(4,1) * t545 - rSges(4,2) * t489 + rSges(4,3) * t488;
t452 = rSges(4,1) * t546 - rSges(4,2) * t491 + rSges(4,3) * t490;
t436 = rSges(5,1) * t487 - rSges(5,2) * t486 + rSges(5,3) * t547;
t435 = Icges(5,1) * t487 - Icges(5,4) * t486 + Icges(5,5) * t547;
t434 = Icges(5,4) * t487 - Icges(5,2) * t486 + Icges(5,6) * t547;
t433 = Icges(5,5) * t487 - Icges(5,6) * t486 + Icges(5,3) * t547;
t428 = -qJD(5) * t467 + t470;
t427 = qJD(5) * t465 + t469;
t422 = rSges(5,1) * t468 + rSges(5,2) * t467 + rSges(5,3) * t489;
t421 = rSges(5,1) * t466 - rSges(5,2) * t465 + rSges(5,3) * t491;
t420 = Icges(5,1) * t468 + Icges(5,4) * t467 + Icges(5,5) * t489;
t419 = Icges(5,1) * t466 - Icges(5,4) * t465 + Icges(5,5) * t491;
t418 = Icges(5,4) * t468 + Icges(5,2) * t467 + Icges(5,6) * t489;
t417 = Icges(5,4) * t466 - Icges(5,2) * t465 + Icges(5,6) * t491;
t416 = Icges(5,5) * t468 + Icges(5,6) * t467 + Icges(5,3) * t489;
t415 = Icges(5,5) * t466 - Icges(5,6) * t465 + Icges(5,3) * t491;
t414 = rSges(6,1) * t464 + rSges(6,2) * t463 + rSges(6,3) * t486;
t404 = t455 * t503 - t479 * t502 + t494;
t403 = -t454 * t503 - t479 * t530 - t536;
t402 = (t454 * t512 + t455 * t515) * t535;
t401 = rSges(6,1) * t432 + rSges(6,2) * t431 - rSges(6,3) * t467;
t399 = rSges(6,1) * t430 + rSges(6,2) * t429 + rSges(6,3) * t465;
t383 = t452 * t503 + t512 * t526 + t532;
t382 = (-t453 - t457) * t503 + t515 * t526 + t529;
t381 = (-t534 + (t452 * t515 + t453 * t512) * qJD(2)) * t506 + t537;
t380 = t421 * t493 - t436 * t469 + t521;
t379 = -t422 * t493 + t436 * t470 + t518;
t378 = -t421 * t470 + t422 * t469 + t522;
t377 = t399 * t459 - t414 * t427 + t519;
t376 = -t401 * t459 + t414 * t428 + t517;
t375 = -t399 * t428 + t401 * t427 + t520;
t374 = -qJD(6) * t467 - t427 * t538 + t459 * t540 + t519;
t373 = qJD(6) * t465 + t428 * t538 - t459 * t539 + t517;
t372 = qJD(6) * t486 + t427 * t539 - t428 * t540 + t520;
t1 = t469 * ((t491 * t415 - t465 * t417 + t466 * t419) * t469 + (t416 * t491 - t465 * t418 + t466 * t420) * t470 + (t433 * t491 - t434 * t465 + t435 * t466) * t493) / 0.2e1 + m(4) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(3) * (t402 ^ 2 + t403 ^ 2 + t404 ^ 2) / 0.2e1 + m(7) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(5) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + t470 * ((t415 * t489 + t417 * t467 + t419 * t468) * t469 + (t489 * t416 + t467 * t418 + t468 * t420) * t470 + (t433 * t489 + t434 * t467 + t435 * t468) * t493) / 0.2e1 + t493 * ((t415 * t547 - t417 * t486 + t419 * t487) * t469 + (t416 * t547 - t418 * t486 + t420 * t487) * t470 + (t433 * t547 - t486 * t434 + t487 * t435) * t493) / 0.2e1 + ((t567 * t429 + t566 * t430 + t568 * t465) * t459 + (t571 * t429 + t569 * t430 + t573 * t465) * t428 + (t572 * t429 + t570 * t430 + t574 * t465) * t427) * t427 / 0.2e1 + ((t567 * t431 + t566 * t432 - t568 * t467) * t459 + (t571 * t431 + t569 * t432 - t573 * t467) * t428 + (t572 * t431 + t570 * t432 - t574 * t467) * t427) * t428 / 0.2e1 + ((t567 * t463 + t566 * t464 + t568 * t486) * t459 + (t571 * t463 + t569 * t464 + t573 * t486) * t428 + (t572 * t463 + t570 * t464 + t574 * t486) * t427) * t459 / 0.2e1 + ((-t558 * t507 + ((t562 * t511 + t564 * t514) * t515 + (t563 * t511 - t565 * t514) * t512) * t506) * t535 + (t561 * t507 + (t559 * t511 + t560 * t514) * t506) * t503) * t503 / 0.2e1 + (m(2) * (t497 ^ 2 + t498 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t564 * t490 + t562 * t491) * t515 + (t565 * t490 + t563 * t491 - t575) * t512) * t535 + (-t560 * t490 + t559 * t491 + t561 * t546) * t503) * t502 / 0.2e1 - (((-t564 * t488 + t562 * t489 + t575) * t515 + (t565 * t488 + t563 * t489) * t512) * t535 + (-t560 * t488 + t559 * t489 - t561 * t545) * t503) * t530 / 0.2e1;
T  = t1;
