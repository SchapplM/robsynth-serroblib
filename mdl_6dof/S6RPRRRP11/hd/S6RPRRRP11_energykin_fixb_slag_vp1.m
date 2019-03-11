% Calculate kinetic energy for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:27
% EndTime: 2019-03-09 06:34:29
% DurationCPUTime: 2.60s
% Computational Cost: add. (5101->297), mult. (13944->456), div. (0->0), fcn. (18136->14), ass. (0->146)
t620 = Icges(6,1) + Icges(7,1);
t619 = Icges(6,4) + Icges(7,4);
t618 = Icges(6,5) + Icges(7,5);
t617 = Icges(6,2) + Icges(7,2);
t616 = Icges(6,6) + Icges(7,6);
t615 = Icges(6,3) + Icges(7,3);
t614 = rSges(7,3) + qJ(6);
t552 = cos(pkin(6));
t593 = cos(pkin(12));
t597 = sin(qJ(1));
t576 = t597 * t593;
t550 = sin(pkin(12));
t558 = cos(qJ(1));
t586 = t558 * t550;
t538 = t552 * t586 + t576;
t556 = sin(qJ(3));
t578 = t558 * t593;
t581 = t597 * t550;
t568 = -t552 * t578 + t581;
t594 = cos(pkin(7));
t562 = t568 * t594;
t551 = sin(pkin(6));
t592 = sin(pkin(7));
t579 = t551 * t592;
t599 = cos(qJ(3));
t518 = t538 * t599 + (-t558 * t579 - t562) * t556;
t580 = t551 * t594;
t531 = -t558 * t580 + t568 * t592;
t555 = sin(qJ(4));
t598 = cos(qJ(4));
t501 = t518 * t598 + t531 * t555;
t577 = t599 * t592;
t587 = t551 * t558;
t517 = t538 * t556 + t562 * t599 + t577 * t587;
t554 = sin(qJ(5));
t557 = cos(qJ(5));
t475 = -t501 * t554 + t517 * t557;
t591 = t517 * t554;
t476 = t501 * t557 + t591;
t500 = t518 * t555 - t531 * t598;
t613 = t616 * t475 + t618 * t476 + t615 * t500;
t539 = -t552 * t581 + t578;
t569 = t552 * t576 + t586;
t604 = t569 * t594 - t597 * t579;
t520 = t539 * t599 - t604 * t556;
t532 = t569 * t592 + t580 * t597;
t503 = t520 * t598 + t532 * t555;
t519 = t539 * t556 + t604 * t599;
t477 = -t503 * t554 + t519 * t557;
t590 = t519 * t554;
t478 = t503 * t557 + t590;
t502 = t520 * t555 - t532 * t598;
t612 = t616 * t477 + t618 * t478 + t615 * t502;
t611 = t617 * t475 + t619 * t476 + t616 * t500;
t610 = t617 * t477 + t619 * t478 + t616 * t502;
t609 = t619 * t475 + t620 * t476 + t618 * t500;
t608 = t619 * t477 + t620 * t478 + t618 * t502;
t574 = t594 * t593;
t530 = t552 * t592 * t556 + (t550 * t599 + t556 * t574) * t551;
t565 = t552 * t594 - t579 * t593;
t516 = t530 * t598 + t555 * t565;
t588 = t551 * t550;
t529 = -t551 * t574 * t599 - t552 * t577 + t556 * t588;
t498 = -t516 * t554 + t529 * t557;
t589 = t529 * t554;
t499 = t516 * t557 + t589;
t515 = t530 * t555 - t565 * t598;
t607 = t616 * t498 + t618 * t499 + t615 * t515;
t606 = t617 * t498 + t619 * t499 + t616 * t515;
t605 = t619 * t498 + t620 * t499 + t618 * t515;
t596 = pkin(5) * t557;
t585 = rSges(7,1) * t476 + rSges(7,2) * t475 + pkin(5) * t591 + t500 * t614 + t501 * t596;
t584 = rSges(7,1) * t478 + rSges(7,2) * t477 + pkin(5) * t590 + t502 * t614 + t503 * t596;
t583 = rSges(7,1) * t499 + rSges(7,2) * t498 + pkin(5) * t589 + t515 * t614 + t516 * t596;
t527 = qJD(3) * t531;
t504 = qJD(4) * t517 + t527;
t528 = qJD(3) * t532;
t505 = qJD(4) * t519 + t528;
t535 = qJD(3) * t565 + qJD(1);
t582 = t551 * t597;
t521 = qJD(4) * t529 + t535;
t575 = -qJD(2) * t587 + qJD(1) * (t558 * pkin(1) + qJ(2) * t582);
t541 = pkin(1) * t597 - qJ(2) * t587;
t547 = qJD(2) * t582;
t573 = t547 + (-t538 * pkin(2) - pkin(9) * t531 - t541) * qJD(1);
t571 = qJD(1) * (t539 * pkin(2) + pkin(9) * t532) + t575;
t494 = pkin(3) * t518 + pkin(10) * t517;
t495 = pkin(3) * t520 + pkin(10) * t519;
t549 = qJD(2) * t552;
t570 = t494 * t528 - t495 * t527 + t549;
t512 = pkin(3) * t530 + pkin(10) * t529;
t567 = -t494 * t535 + t512 * t527 + t573;
t471 = pkin(4) * t501 + pkin(11) * t500;
t472 = pkin(4) * t503 + pkin(11) * t502;
t566 = t505 * t471 - t472 * t504 + t570;
t564 = t535 * t495 - t512 * t528 + t571;
t493 = pkin(4) * t516 + pkin(11) * t515;
t561 = -t471 * t521 + t504 * t493 + t567;
t560 = t521 * t472 - t493 * t505 + t564;
t545 = t558 * rSges(2,1) - rSges(2,2) * t597;
t544 = rSges(2,1) * t597 + t558 * rSges(2,2);
t511 = qJD(1) * (t539 * rSges(3,1) - rSges(3,2) * t569 + rSges(3,3) * t582) + t575;
t510 = t547 + (-t538 * rSges(3,1) + rSges(3,2) * t568 + rSges(3,3) * t587 - t541) * qJD(1);
t509 = t530 * rSges(4,1) - t529 * rSges(4,2) + rSges(4,3) * t565;
t508 = Icges(4,1) * t530 - Icges(4,4) * t529 + Icges(4,5) * t565;
t507 = Icges(4,4) * t530 - Icges(4,2) * t529 + Icges(4,6) * t565;
t506 = Icges(4,5) * t530 - Icges(4,6) * t529 + Icges(4,3) * t565;
t496 = qJD(5) * t515 + t521;
t490 = rSges(4,1) * t520 - rSges(4,2) * t519 + rSges(4,3) * t532;
t489 = rSges(4,1) * t518 - rSges(4,2) * t517 + rSges(4,3) * t531;
t488 = Icges(4,1) * t520 - Icges(4,4) * t519 + Icges(4,5) * t532;
t487 = Icges(4,1) * t518 - Icges(4,4) * t517 + Icges(4,5) * t531;
t486 = Icges(4,4) * t520 - Icges(4,2) * t519 + Icges(4,6) * t532;
t485 = Icges(4,4) * t518 - Icges(4,2) * t517 + Icges(4,6) * t531;
t484 = Icges(4,5) * t520 - Icges(4,6) * t519 + Icges(4,3) * t532;
t483 = Icges(4,5) * t518 - Icges(4,6) * t517 + Icges(4,3) * t531;
t482 = rSges(5,1) * t516 - rSges(5,2) * t515 + rSges(5,3) * t529;
t481 = Icges(5,1) * t516 - Icges(5,4) * t515 + Icges(5,5) * t529;
t480 = Icges(5,4) * t516 - Icges(5,2) * t515 + Icges(5,6) * t529;
t479 = Icges(5,5) * t516 - Icges(5,6) * t515 + Icges(5,3) * t529;
t474 = qJD(5) * t502 + t505;
t473 = qJD(5) * t500 + t504;
t468 = rSges(5,1) * t503 - rSges(5,2) * t502 + rSges(5,3) * t519;
t467 = rSges(5,1) * t501 - rSges(5,2) * t500 + rSges(5,3) * t517;
t466 = Icges(5,1) * t503 - Icges(5,4) * t502 + Icges(5,5) * t519;
t465 = Icges(5,1) * t501 - Icges(5,4) * t500 + Icges(5,5) * t517;
t464 = Icges(5,4) * t503 - Icges(5,2) * t502 + Icges(5,6) * t519;
t463 = Icges(5,4) * t501 - Icges(5,2) * t500 + Icges(5,6) * t517;
t462 = Icges(5,5) * t503 - Icges(5,6) * t502 + Icges(5,3) * t519;
t461 = Icges(5,5) * t501 - Icges(5,6) * t500 + Icges(5,3) * t517;
t460 = rSges(6,1) * t499 + rSges(6,2) * t498 + rSges(6,3) * t515;
t450 = t490 * t535 - t509 * t528 + t571;
t449 = -t489 * t535 + t509 * t527 + t573;
t448 = t549 + (t489 * t532 - t490 * t531) * qJD(3);
t447 = rSges(6,1) * t478 + rSges(6,2) * t477 + rSges(6,3) * t502;
t445 = rSges(6,1) * t476 + rSges(6,2) * t475 + rSges(6,3) * t500;
t429 = t468 * t521 - t482 * t505 + t564;
t428 = -t467 * t521 + t482 * t504 + t567;
t427 = t467 * t505 - t468 * t504 + t570;
t426 = t447 * t496 - t460 * t474 + t560;
t425 = -t445 * t496 + t460 * t473 + t561;
t424 = t445 * t474 - t447 * t473 + t566;
t423 = qJD(6) * t500 - t474 * t583 + t496 * t584 + t560;
t422 = qJD(6) * t502 + t473 * t583 - t496 * t585 + t561;
t421 = qJD(6) * t515 - t473 * t584 + t474 * t585 + t566;
t1 = t505 * ((t519 * t462 - t502 * t464 + t503 * t466) * t505 + (t461 * t519 - t463 * t502 + t465 * t503) * t504 + (t479 * t519 - t480 * t502 + t481 * t503) * t521) / 0.2e1 + t504 * ((t462 * t517 - t464 * t500 + t466 * t501) * t505 + (t517 * t461 - t500 * t463 + t501 * t465) * t504 + (t479 * t517 - t480 * t500 + t481 * t501) * t521) / 0.2e1 + t521 * ((t462 * t529 - t464 * t515 + t466 * t516) * t505 + (t461 * t529 - t463 * t515 + t465 * t516) * t504 + (t529 * t479 - t515 * t480 + t516 * t481) * t521) / 0.2e1 + ((t506 * t531 - t507 * t517 + t508 * t518) * t535 + ((t484 * t531 - t486 * t517 + t488 * t518) * t532 + (t483 * t531 - t485 * t517 + t487 * t518) * t531) * qJD(3)) * t527 / 0.2e1 + t535 * ((t506 * t565 - t529 * t507 + t530 * t508) * t535 + ((t484 * t565 - t529 * t486 + t530 * t488) * t532 + (t483 * t565 - t529 * t485 + t530 * t487) * t531) * qJD(3)) / 0.2e1 + m(7) * (t421 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + m(6) * (t424 ^ 2 + t425 ^ 2 + t426 ^ 2) / 0.2e1 + m(5) * (t427 ^ 2 + t428 ^ 2 + t429 ^ 2) / 0.2e1 + m(4) * (t448 ^ 2 + t449 ^ 2 + t450 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t552 ^ 2 + t510 ^ 2 + t511 ^ 2) / 0.2e1 + ((t506 * t532 - t507 * t519 + t508 * t520) * t535 + ((t484 * t532 - t486 * t519 + t488 * t520) * t532 + (t483 * t532 - t485 * t519 + t487 * t520) * t531) * qJD(3)) * t528 / 0.2e1 + ((t606 * t475 + t605 * t476 + t607 * t500) * t496 + (t610 * t475 + t608 * t476 + t612 * t500) * t474 + (t611 * t475 + t609 * t476 + t613 * t500) * t473) * t473 / 0.2e1 + ((t606 * t477 + t605 * t478 + t607 * t502) * t496 + (t610 * t477 + t608 * t478 + t612 * t502) * t474 + (t611 * t477 + t609 * t478 + t613 * t502) * t473) * t474 / 0.2e1 + ((t606 * t498 + t605 * t499 + t607 * t515) * t496 + (t610 * t498 + t608 * t499 + t612 * t515) * t474 + (t611 * t498 + t609 * t499 + t613 * t515) * t473) * t496 / 0.2e1 + (Icges(2,3) + (Icges(3,5) * t552 + (Icges(3,1) * t550 + Icges(3,4) * t593) * t551) * t588 + t551 * t593 * (Icges(3,6) * t552 + (Icges(3,4) * t550 + Icges(3,2) * t593) * t551) + t552 * (Icges(3,3) * t552 + (Icges(3,5) * t550 + Icges(3,6) * t593) * t551) + m(2) * (t544 ^ 2 + t545 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
