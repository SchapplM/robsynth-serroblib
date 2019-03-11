% Calculate kinetic energy for
% S6RRPRRP14
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:43
% EndTime: 2019-03-09 13:03:46
% DurationCPUTime: 2.89s
% Computational Cost: add. (2058->284), mult. (5160->424), div. (0->0), fcn. (6189->10), ass. (0->139)
t594 = Icges(3,1) + Icges(4,2);
t593 = Icges(6,1) + Icges(7,1);
t592 = Icges(3,4) + Icges(4,6);
t591 = -Icges(6,4) + Icges(7,5);
t590 = Icges(7,4) + Icges(6,5);
t589 = Icges(3,5) - Icges(4,4);
t588 = Icges(3,2) + Icges(4,3);
t587 = Icges(6,2) + Icges(7,3);
t586 = Icges(7,2) + Icges(6,3);
t585 = Icges(3,6) - Icges(4,5);
t584 = Icges(6,6) - Icges(7,6);
t583 = Icges(3,3) + Icges(4,1);
t582 = rSges(7,1) + pkin(5);
t581 = rSges(7,3) + qJ(6);
t517 = sin(pkin(6));
t518 = cos(pkin(6));
t523 = cos(qJ(2));
t524 = cos(qJ(1));
t549 = t523 * t524;
t521 = sin(qJ(2));
t522 = sin(qJ(1));
t553 = t521 * t522;
t499 = -t518 * t549 + t553;
t550 = t522 * t523;
t552 = t521 * t524;
t500 = t518 * t552 + t550;
t501 = t518 * t550 + t552;
t502 = -t518 * t553 + t549;
t551 = t522 * t517;
t554 = t517 * t524;
t563 = (-t499 * t585 + t500 * t589 - t554 * t583) * t524 + (t501 * t585 - t502 * t589 - t551 * t583) * t522;
t580 = t563 * t517;
t520 = sin(qJ(4));
t557 = cos(qJ(4));
t540 = t517 * t557;
t475 = t501 * t520 + t522 * t540;
t519 = sin(qJ(5));
t556 = cos(qJ(5));
t438 = t475 * t519 - t502 * t556;
t439 = t475 * t556 + t502 * t519;
t474 = -t501 * t557 + t520 * t551;
t579 = t438 * t587 + t439 * t591 - t474 * t584;
t477 = t499 * t520 - t524 * t540;
t440 = t477 * t519 - t500 * t556;
t441 = t477 * t556 + t500 * t519;
t476 = t499 * t557 + t520 * t554;
t578 = t440 * t587 + t441 * t591 + t476 * t584;
t577 = -t438 * t584 + t439 * t590 + t474 * t586;
t576 = -t440 * t584 + t441 * t590 - t476 * t586;
t575 = t438 * t591 + t439 * t593 + t474 * t590;
t574 = t440 * t591 + t441 * t593 - t476 * t590;
t498 = -t517 * t523 * t520 + t518 * t557;
t555 = t517 * t521;
t472 = t498 * t519 - t555 * t556;
t473 = t498 * t556 + t519 * t555;
t497 = t518 * t520 + t523 * t540;
t573 = t472 * t587 + t473 * t591 - t497 * t584;
t572 = -t472 * t584 + t473 * t590 + t497 * t586;
t571 = t472 * t591 + t473 * t593 + t497 * t590;
t570 = t501 * t588 - t502 * t592 - t551 * t585;
t569 = t499 * t588 - t500 * t592 + t554 * t585;
t568 = -t592 * t501 + t502 * t594 + t589 * t551;
t567 = t592 * t499 - t500 * t594 + t589 * t554;
t566 = t583 * t518 + (t521 * t589 + t523 * t585) * t517;
t565 = t585 * t518 + (t521 * t592 + t523 * t588) * t517;
t564 = t589 * t518 + (t521 * t594 + t592 * t523) * t517;
t548 = rSges(7,2) * t474 + t581 * t438 + t582 * t439;
t547 = -rSges(7,2) * t476 + t581 * t440 + t582 * t441;
t546 = rSges(7,2) * t497 + t581 * t472 + t582 * t473;
t466 = pkin(2) * t500 + qJ(3) * t499;
t467 = pkin(2) * t502 + qJ(3) * t501;
t543 = qJD(2) * t517;
t514 = t522 * t543;
t539 = t524 * t543;
t545 = t466 * t514 + t467 * t539;
t478 = qJD(4) * t502 + t514;
t544 = qJD(1) * (pkin(1) * t522 - pkin(8) * t554);
t542 = qJD(3) * t523;
t515 = qJD(2) * t518 + qJD(1);
t505 = qJD(1) * (pkin(1) * t524 + pkin(8) * t551);
t541 = qJD(3) * t499 + t515 * t467 + t505;
t504 = qJD(4) * t555 + t515;
t538 = qJD(3) * t501 - t544;
t503 = (pkin(2) * t521 - qJ(3) * t523) * t517;
t535 = (-rSges(4,1) * t518 - (-rSges(4,2) * t521 - rSges(4,3) * t523) * t517 - t503) * t543;
t534 = (-pkin(3) * t518 - pkin(9) * t555 - t503) * t543;
t479 = qJD(4) * t500 - t539;
t480 = pkin(3) * t551 + pkin(9) * t502;
t481 = -pkin(3) * t554 + pkin(9) * t500;
t531 = t480 * t539 + t481 * t514 - t517 * t542 + t545;
t530 = t515 * t480 + t522 * t534 + t541;
t434 = pkin(4) * t475 + pkin(10) * t474;
t435 = pkin(4) * t477 - pkin(10) * t476;
t529 = -t434 * t479 + t478 * t435 + t531;
t465 = pkin(4) * t498 + pkin(10) * t497;
t528 = t504 * t434 - t465 * t478 + t530;
t527 = (-t466 - t481) * t515 + t524 * t534 + t538;
t526 = -t435 * t504 + t479 * t465 + t527;
t509 = rSges(2,1) * t524 - rSges(2,2) * t522;
t508 = rSges(2,1) * t522 + rSges(2,2) * t524;
t488 = rSges(3,3) * t518 + (rSges(3,1) * t521 + rSges(3,2) * t523) * t517;
t468 = qJD(5) * t497 + t504;
t464 = rSges(3,1) * t502 - rSges(3,2) * t501 + rSges(3,3) * t551;
t463 = rSges(3,1) * t500 - rSges(3,2) * t499 - rSges(3,3) * t554;
t462 = -rSges(4,1) * t554 - rSges(4,2) * t500 + rSges(4,3) * t499;
t461 = rSges(4,1) * t551 - rSges(4,2) * t502 + rSges(4,3) * t501;
t445 = rSges(5,1) * t498 - rSges(5,2) * t497 + rSges(5,3) * t555;
t444 = Icges(5,1) * t498 - Icges(5,4) * t497 + Icges(5,5) * t555;
t443 = Icges(5,4) * t498 - Icges(5,2) * t497 + Icges(5,6) * t555;
t442 = Icges(5,5) * t498 - Icges(5,6) * t497 + Icges(5,3) * t555;
t437 = -qJD(5) * t476 + t479;
t436 = qJD(5) * t474 + t478;
t430 = rSges(5,1) * t477 + rSges(5,2) * t476 + rSges(5,3) * t500;
t429 = rSges(5,1) * t475 - rSges(5,2) * t474 + rSges(5,3) * t502;
t428 = Icges(5,1) * t477 + Icges(5,4) * t476 + Icges(5,5) * t500;
t427 = Icges(5,1) * t475 - Icges(5,4) * t474 + Icges(5,5) * t502;
t426 = Icges(5,4) * t477 + Icges(5,2) * t476 + Icges(5,6) * t500;
t425 = Icges(5,4) * t475 - Icges(5,2) * t474 + Icges(5,6) * t502;
t424 = Icges(5,5) * t477 + Icges(5,6) * t476 + Icges(5,3) * t500;
t423 = Icges(5,5) * t475 - Icges(5,6) * t474 + Icges(5,3) * t502;
t422 = rSges(6,1) * t473 - rSges(6,2) * t472 + rSges(6,3) * t497;
t413 = t464 * t515 - t488 * t514 + t505;
t412 = -t463 * t515 - t488 * t539 - t544;
t411 = (t463 * t522 + t464 * t524) * t543;
t408 = rSges(6,1) * t441 - rSges(6,2) * t440 - rSges(6,3) * t476;
t406 = rSges(6,1) * t439 - rSges(6,2) * t438 + rSges(6,3) * t474;
t392 = t461 * t515 + t522 * t535 + t541;
t391 = (-t462 - t466) * t515 + t524 * t535 + t538;
t390 = (-t542 + (t461 * t524 + t462 * t522) * qJD(2)) * t517 + t545;
t389 = t429 * t504 - t445 * t478 + t530;
t388 = -t430 * t504 + t445 * t479 + t527;
t387 = -t429 * t479 + t430 * t478 + t531;
t386 = t406 * t468 - t422 * t436 + t528;
t385 = -t408 * t468 + t422 * t437 + t526;
t384 = -t406 * t437 + t408 * t436 + t529;
t383 = qJD(6) * t440 - t436 * t546 + t468 * t548 + t528;
t382 = qJD(6) * t438 + t437 * t546 - t468 * t547 + t526;
t381 = qJD(6) * t472 + t436 * t547 - t437 * t548 + t529;
t1 = m(7) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(5) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(4) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(3) * (t411 ^ 2 + t412 ^ 2 + t413 ^ 2) / 0.2e1 + m(6) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + t478 * ((t502 * t423 - t474 * t425 + t475 * t427) * t478 + (t424 * t502 - t426 * t474 + t428 * t475) * t479 + (t442 * t502 - t443 * t474 + t444 * t475) * t504) / 0.2e1 + t479 * ((t423 * t500 + t425 * t476 + t427 * t477) * t478 + (t500 * t424 + t476 * t426 + t477 * t428) * t479 + (t442 * t500 + t443 * t476 + t444 * t477) * t504) / 0.2e1 + t504 * ((t423 * t555 - t425 * t497 + t427 * t498) * t478 + (t424 * t555 - t426 * t497 + t428 * t498) * t479 + (t442 * t555 - t497 * t443 + t498 * t444) * t504) / 0.2e1 + ((t573 * t438 + t571 * t439 + t572 * t474) * t468 + (t578 * t438 + t574 * t439 + t576 * t474) * t437 + (t579 * t438 + t575 * t439 + t577 * t474) * t436) * t436 / 0.2e1 + ((t573 * t440 + t571 * t441 - t572 * t476) * t468 + (t578 * t440 + t574 * t441 - t576 * t476) * t437 + (t579 * t440 + t575 * t441 - t577 * t476) * t436) * t437 / 0.2e1 + ((t573 * t472 + t571 * t473 + t572 * t497) * t468 + (t578 * t472 + t574 * t473 + t576 * t497) * t437 + (t579 * t472 + t575 * t473 + t577 * t497) * t436) * t468 / 0.2e1 + ((-t563 * t518 + ((t567 * t521 + t569 * t523) * t524 + (t568 * t521 - t570 * t523) * t522) * t517) * t543 + (t566 * t518 + (t564 * t521 + t565 * t523) * t517) * t515) * t515 / 0.2e1 + (m(2) * (t508 ^ 2 + t509 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t569 * t501 + t567 * t502) * t524 + (t570 * t501 + t568 * t502 - t580) * t522) * t543 + (-t565 * t501 + t564 * t502 + t566 * t551) * t515) * t514 / 0.2e1 - (((-t569 * t499 + t567 * t500 + t580) * t524 + (t570 * t499 + t568 * t500) * t522) * t543 + (-t565 * t499 + t564 * t500 - t566 * t554) * t515) * t539 / 0.2e1;
T  = t1;
