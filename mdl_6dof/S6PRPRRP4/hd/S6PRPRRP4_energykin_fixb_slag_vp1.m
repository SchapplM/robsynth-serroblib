% Calculate kinetic energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:13
% EndTime: 2019-03-08 20:09:16
% DurationCPUTime: 3.04s
% Computational Cost: add. (3061->325), mult. (5618->497), div. (0->0), fcn. (6813->12), ass. (0->152)
t587 = Icges(6,1) + Icges(7,1);
t586 = -Icges(6,4) + Icges(7,5);
t585 = Icges(7,4) + Icges(6,5);
t584 = Icges(6,2) + Icges(7,3);
t583 = Icges(7,2) + Icges(6,3);
t582 = -Icges(6,6) + Icges(7,6);
t581 = rSges(7,1) + pkin(5);
t580 = rSges(7,3) + qJ(6);
t522 = sin(pkin(10));
t525 = cos(pkin(10));
t530 = cos(qJ(2));
t526 = cos(pkin(6));
t529 = sin(qJ(2));
t558 = t526 * t529;
t509 = t522 * t530 + t525 * t558;
t548 = pkin(11) + qJ(4);
t520 = sin(t548);
t542 = cos(t548);
t523 = sin(pkin(6));
t561 = t523 * t525;
t483 = t509 * t542 - t520 * t561;
t557 = t526 * t530;
t508 = t522 * t529 - t525 * t557;
t528 = sin(qJ(5));
t565 = cos(qJ(5));
t454 = t483 * t528 - t508 * t565;
t455 = t483 * t565 + t508 * t528;
t537 = t523 * t542;
t482 = t509 * t520 + t525 * t537;
t579 = t584 * t454 + t586 * t455 + t582 * t482;
t511 = -t522 * t558 + t525 * t530;
t562 = t522 * t523;
t485 = t511 * t542 + t520 * t562;
t510 = t522 * t557 + t525 * t529;
t456 = t485 * t528 - t510 * t565;
t457 = t485 * t565 + t510 * t528;
t484 = t511 * t520 - t522 * t537;
t578 = t584 * t456 + t586 * t457 + t582 * t484;
t577 = t582 * t454 + t585 * t455 + t583 * t482;
t576 = t582 * t456 + t585 * t457 + t583 * t484;
t575 = t586 * t454 + t587 * t455 + t585 * t482;
t574 = t586 * t456 + t587 * t457 + t585 * t484;
t498 = t526 * t520 + t529 * t537;
t559 = t523 * t530;
t486 = t498 * t528 + t559 * t565;
t487 = t498 * t565 - t528 * t559;
t560 = t523 * t529;
t497 = t520 * t560 - t526 * t542;
t573 = t584 * t486 + t586 * t487 + t582 * t497;
t572 = t582 * t486 + t585 * t487 + t583 * t497;
t571 = t586 * t486 + t587 * t487 + t585 * t497;
t521 = sin(pkin(11));
t524 = cos(pkin(11));
t506 = -t521 * t560 + t524 * t526;
t563 = t521 * t526;
t507 = t524 * t560 + t563;
t462 = Icges(4,5) * t507 + Icges(4,6) * t506 - Icges(4,3) * t559;
t495 = Icges(3,6) * t526 + (Icges(3,4) * t529 + Icges(3,2) * t530) * t523;
t570 = t462 - t495;
t569 = qJD(2) ^ 2;
t564 = pkin(3) * t524;
t555 = rSges(7,2) * t482 + t580 * t454 + t455 * t581;
t554 = rSges(7,2) * t484 + t580 * t456 + t457 * t581;
t553 = rSges(7,2) * t497 + t580 * t486 + t487 * t581;
t480 = pkin(2) * t511 + qJ(3) * t510;
t519 = qJD(2) * t526;
t552 = qJD(3) * t508 + t480 * t519;
t551 = qJD(2) * t523;
t517 = t522 * t551;
t492 = qJD(4) * t510 + t517;
t550 = qJD(3) * t530;
t547 = t521 * t562;
t546 = t521 * t561;
t545 = t525 * t551;
t479 = pkin(2) * t509 + qJ(3) * t508;
t544 = t479 * t517 + t480 * t545 + qJD(1);
t512 = (pkin(2) * t529 - qJ(3) * t530) * t523;
t541 = (-t507 * rSges(4,1) - t506 * rSges(4,2) + rSges(4,3) * t559 - t512) * t523;
t540 = (-pkin(3) * t563 - (-pkin(8) * t530 + t529 * t564) * t523 - t512) * t523;
t493 = qJD(4) * t508 - t545;
t513 = -qJD(4) * t559 + t519;
t435 = -pkin(3) * t546 + pkin(8) * t508 + t509 * t564;
t436 = pkin(3) * t547 + pkin(8) * t510 + t511 * t564;
t536 = t435 * t517 + t436 * t545 - t523 * t550 + t544;
t535 = qJD(2) * t522 * t540 + t436 * t519 + t552;
t449 = pkin(4) * t483 + pkin(9) * t482;
t450 = pkin(4) * t485 + pkin(9) * t484;
t534 = t492 * t449 - t493 * t450 + t536;
t505 = qJD(3) * t510;
t533 = t505 + ((-t435 - t479) * t526 + t525 * t540) * qJD(2);
t472 = pkin(4) * t498 + pkin(9) * t497;
t532 = t513 * t450 - t472 * t492 + t535;
t531 = -t449 * t513 + t493 * t472 + t533;
t499 = t526 * rSges(3,3) + (rSges(3,1) * t529 + rSges(3,2) * t530) * t523;
t496 = Icges(3,5) * t526 + (Icges(3,1) * t529 + Icges(3,4) * t530) * t523;
t494 = Icges(3,3) * t526 + (Icges(3,5) * t529 + Icges(3,6) * t530) * t523;
t491 = t511 * t524 + t547;
t490 = -t511 * t521 + t524 * t562;
t489 = t509 * t524 - t546;
t488 = -t509 * t521 - t524 * t561;
t481 = qJD(5) * t497 + t513;
t477 = rSges(3,1) * t511 - rSges(3,2) * t510 + rSges(3,3) * t562;
t476 = rSges(3,1) * t509 - rSges(3,2) * t508 - rSges(3,3) * t561;
t470 = Icges(3,1) * t511 - Icges(3,4) * t510 + Icges(3,5) * t562;
t469 = Icges(3,1) * t509 - Icges(3,4) * t508 - Icges(3,5) * t561;
t468 = Icges(3,4) * t511 - Icges(3,2) * t510 + Icges(3,6) * t562;
t467 = Icges(3,4) * t509 - Icges(3,2) * t508 - Icges(3,6) * t561;
t466 = Icges(3,5) * t511 - Icges(3,6) * t510 + Icges(3,3) * t562;
t465 = Icges(3,5) * t509 - Icges(3,6) * t508 - Icges(3,3) * t561;
t464 = Icges(4,1) * t507 + Icges(4,4) * t506 - Icges(4,5) * t559;
t463 = Icges(4,4) * t507 + Icges(4,2) * t506 - Icges(4,6) * t559;
t461 = t498 * rSges(5,1) - t497 * rSges(5,2) - rSges(5,3) * t559;
t460 = Icges(5,1) * t498 - Icges(5,4) * t497 - Icges(5,5) * t559;
t459 = Icges(5,4) * t498 - Icges(5,2) * t497 - Icges(5,6) * t559;
t458 = Icges(5,5) * t498 - Icges(5,6) * t497 - Icges(5,3) * t559;
t453 = qJD(5) * t482 + t493;
t452 = qJD(5) * t484 + t492;
t447 = (-t476 * t526 - t499 * t561) * qJD(2);
t446 = (t477 * t526 - t499 * t562) * qJD(2);
t445 = rSges(4,1) * t491 + rSges(4,2) * t490 + rSges(4,3) * t510;
t444 = rSges(4,1) * t489 + rSges(4,2) * t488 + rSges(4,3) * t508;
t443 = Icges(4,1) * t491 + Icges(4,4) * t490 + Icges(4,5) * t510;
t442 = Icges(4,1) * t489 + Icges(4,4) * t488 + Icges(4,5) * t508;
t441 = Icges(4,4) * t491 + Icges(4,2) * t490 + Icges(4,6) * t510;
t440 = Icges(4,4) * t489 + Icges(4,2) * t488 + Icges(4,6) * t508;
t439 = Icges(4,5) * t491 + Icges(4,6) * t490 + Icges(4,3) * t510;
t438 = Icges(4,5) * t489 + Icges(4,6) * t488 + Icges(4,3) * t508;
t434 = rSges(5,1) * t485 - rSges(5,2) * t484 + rSges(5,3) * t510;
t433 = rSges(5,1) * t483 - rSges(5,2) * t482 + rSges(5,3) * t508;
t432 = Icges(5,1) * t485 - Icges(5,4) * t484 + Icges(5,5) * t510;
t431 = Icges(5,1) * t483 - Icges(5,4) * t482 + Icges(5,5) * t508;
t430 = Icges(5,4) * t485 - Icges(5,2) * t484 + Icges(5,6) * t510;
t429 = Icges(5,4) * t483 - Icges(5,2) * t482 + Icges(5,6) * t508;
t428 = Icges(5,5) * t485 - Icges(5,6) * t484 + Icges(5,3) * t510;
t427 = Icges(5,5) * t483 - Icges(5,6) * t482 + Icges(5,3) * t508;
t426 = rSges(6,1) * t487 - rSges(6,2) * t486 + rSges(6,3) * t497;
t414 = qJD(1) + (t476 * t522 + t477 * t525) * t551;
t411 = rSges(6,1) * t457 - rSges(6,2) * t456 + rSges(6,3) * t484;
t409 = rSges(6,1) * t455 - rSges(6,2) * t454 + rSges(6,3) * t482;
t395 = t505 + ((-t444 - t479) * t526 + t525 * t541) * qJD(2);
t394 = (t445 * t526 + t522 * t541) * qJD(2) + t552;
t393 = (-t550 + (t444 * t522 + t445 * t525) * qJD(2)) * t523 + t544;
t392 = -t433 * t513 + t461 * t493 + t533;
t391 = t434 * t513 - t461 * t492 + t535;
t390 = t492 * t433 - t493 * t434 + t536;
t389 = -t409 * t481 + t426 * t453 + t531;
t388 = t411 * t481 - t426 * t452 + t532;
t387 = t452 * t409 - t453 * t411 + t534;
t386 = qJD(6) * t456 + t453 * t553 - t481 * t555 + t531;
t385 = qJD(6) * t454 - t452 * t553 + t481 * t554 + t532;
t384 = qJD(6) * t486 + t452 * t555 - t453 * t554 + t534;
t1 = t493 * ((t428 * t508 - t430 * t482 + t432 * t483) * t492 + (t427 * t508 - t429 * t482 + t431 * t483) * t493 + (t458 * t508 - t459 * t482 + t460 * t483) * t513) / 0.2e1 + t513 * ((-t428 * t559 - t497 * t430 + t498 * t432) * t492 + (-t427 * t559 - t497 * t429 + t498 * t431) * t493 + (-t458 * t559 - t497 * t459 + t498 * t460) * t513) / 0.2e1 + t492 * ((t428 * t510 - t430 * t484 + t432 * t485) * t492 + (t427 * t510 - t429 * t484 + t431 * t485) * t493 + (t458 * t510 - t459 * t484 + t460 * t485) * t513) / 0.2e1 + m(5) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(4) * (t393 ^ 2 + t394 ^ 2 + t395 ^ 2) / 0.2e1 + m(3) * (t414 ^ 2 + t446 ^ 2 + t447 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(6) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(7) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + ((t456 * t573 + t457 * t571 + t484 * t572) * t481 + (t456 * t579 + t575 * t457 + t577 * t484) * t453 + (t578 * t456 + t574 * t457 + t576 * t484) * t452) * t452 / 0.2e1 + ((t454 * t573 + t455 * t571 + t482 * t572) * t481 + (t579 * t454 + t575 * t455 + t577 * t482) * t453 + (t454 * t578 + t455 * t574 + t482 * t576) * t452) * t453 / 0.2e1 + ((t573 * t486 + t571 * t487 + t572 * t497) * t481 + (t486 * t579 + t575 * t487 + t577 * t497) * t453 + (t486 * t578 + t487 * t574 + t497 * t576) * t452) * t481 / 0.2e1 - (((t439 * t508 + t441 * t488 + t443 * t489) * t522 - (t438 * t508 + t440 * t488 + t442 * t489) * t525) * t523 + (-t466 * t561 - t468 * t508 + t470 * t509) * t562 - (-t465 * t561 - t467 * t508 + t469 * t509) * t561 + (t463 * t488 + t464 * t489 - t494 * t561 + t496 * t509 + t508 * t570) * t526) * t569 * t561 / 0.2e1 + (((-t439 * t559 + t506 * t441 + t507 * t443) * t562 - (-t438 * t559 + t506 * t440 + t507 * t442) * t561 + ((t468 * t530 + t470 * t529) * t522 - (t467 * t530 + t469 * t529) * t525) * t523 ^ 2 + (-t462 * t559 + t506 * t463 + t507 * t464 + (-t465 * t525 + t466 * t522 + t495 * t530 + t496 * t529) * t523 + t526 * t494) * t526) * t526 + ((t466 * t562 - t468 * t510 + t470 * t511) * t562 - (t465 * t562 - t467 * t510 + t469 * t511) * t561 + ((t439 * t510 + t441 * t490 + t443 * t491) * t522 - (t438 * t510 + t440 * t490 + t442 * t491) * t525) * t523 + (t463 * t490 + t464 * t491 + t494 * t562 + t496 * t511 + t510 * t570) * t526) * t562) * t569 / 0.2e1;
T  = t1;
