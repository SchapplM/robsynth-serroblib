% Calculate kinetic energy for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:16
% EndTime: 2019-03-10 02:16:19
% DurationCPUTime: 2.94s
% Computational Cost: add. (3106->332), mult. (6776->509), div. (0->0), fcn. (8331->12), ass. (0->156)
t596 = Icges(6,1) + Icges(7,1);
t595 = -Icges(6,4) + Icges(7,5);
t594 = Icges(7,4) + Icges(6,5);
t593 = Icges(6,2) + Icges(7,3);
t592 = Icges(7,2) + Icges(6,3);
t591 = -Icges(6,6) + Icges(7,6);
t590 = rSges(7,1) + pkin(5);
t589 = rSges(7,3) + qJ(6);
t537 = sin(qJ(2));
t538 = sin(qJ(1));
t540 = cos(qJ(2));
t541 = cos(qJ(1));
t571 = cos(pkin(6));
t553 = t541 * t571;
t516 = t537 * t538 - t540 * t553;
t517 = t537 * t553 + t538 * t540;
t534 = sin(pkin(6));
t566 = t534 * t541;
t478 = Icges(3,5) * t517 - Icges(3,6) * t516 - Icges(3,3) * t566;
t554 = t538 * t571;
t518 = t541 * t537 + t540 * t554;
t519 = -t537 * t554 + t541 * t540;
t568 = t534 * t538;
t479 = Icges(3,5) * t519 - Icges(3,6) * t518 + Icges(3,3) * t568;
t588 = (t478 * t541 - t479 * t538) * t534;
t536 = sin(qJ(3));
t574 = cos(qJ(3));
t500 = t517 * t574 - t536 * t566;
t565 = qJ(4) + qJ(5);
t533 = sin(t565);
t555 = cos(t565);
t465 = t500 * t533 - t516 * t555;
t466 = t500 * t555 + t516 * t533;
t557 = t534 * t574;
t499 = t517 * t536 + t541 * t557;
t587 = t593 * t465 + t595 * t466 + t591 * t499;
t502 = t519 * t574 + t536 * t568;
t467 = t502 * t533 - t518 * t555;
t468 = t502 * t555 + t518 * t533;
t501 = t519 * t536 - t538 * t557;
t586 = t593 * t467 + t595 * t468 + t591 * t501;
t585 = t591 * t465 + t594 * t466 + t592 * t499;
t584 = t591 * t467 + t594 * t468 + t592 * t501;
t583 = t595 * t465 + t596 * t466 + t594 * t499;
t582 = t595 * t467 + t596 * t468 + t594 * t501;
t515 = t536 * t571 + t537 * t557;
t567 = t534 * t540;
t493 = t515 * t533 + t555 * t567;
t494 = t515 * t555 - t533 * t567;
t514 = t534 * t536 * t537 - t571 * t574;
t581 = t593 * t493 + t595 * t494 + t591 * t514;
t580 = t591 * t493 + t594 * t494 + t592 * t514;
t579 = t595 * t493 + t596 * t494 + t594 * t514;
t539 = cos(qJ(4));
t573 = pkin(4) * t539;
t535 = sin(qJ(4));
t570 = t516 * t535;
t569 = t518 * t535;
t564 = rSges(7,2) * t499 + t589 * t465 + t590 * t466;
t563 = rSges(7,2) * t501 + t589 * t467 + t590 * t468;
t562 = rSges(7,2) * t514 + t589 * t493 + t590 * t494;
t490 = pkin(2) * t517 + pkin(9) * t516;
t491 = pkin(2) * t519 + pkin(9) * t518;
t559 = qJD(2) * t534;
t528 = t538 * t559;
t556 = t541 * t559;
t561 = t490 * t528 + t491 * t556;
t503 = qJD(3) * t518 + t528;
t560 = qJD(1) * (pkin(1) * t538 - pkin(8) * t566);
t529 = qJD(2) * t571 + qJD(1);
t558 = t535 * t567;
t463 = qJD(4) * t501 + t503;
t504 = qJD(3) * t516 - t556;
t461 = pkin(3) * t500 + pkin(10) * t499;
t462 = pkin(3) * t502 + pkin(10) * t501;
t551 = t503 * t461 - t462 * t504 + t561;
t464 = qJD(4) * t499 + t504;
t521 = -qJD(3) * t567 + t529;
t520 = (pkin(2) * t537 - pkin(9) * t540) * t534;
t522 = qJD(1) * (pkin(1) * t541 + pkin(8) * t568);
t550 = t529 * t491 - t520 * t528 + t522;
t492 = qJD(4) * t514 + t521;
t402 = pkin(4) * t570 + pkin(11) * t499 + t500 * t573;
t403 = pkin(4) * t569 + pkin(11) * t501 + t502 * t573;
t549 = t463 * t402 - t403 * t464 + t551;
t548 = -t490 * t529 - t520 * t556 - t560;
t489 = pkin(3) * t515 + pkin(10) * t514;
t547 = t521 * t462 - t489 * t503 + t550;
t546 = -t461 * t521 + t504 * t489 + t548;
t445 = -pkin(4) * t558 + pkin(11) * t514 + t515 * t573;
t545 = t492 * t403 - t445 * t463 + t547;
t544 = -t402 * t492 + t464 * t445 + t546;
t525 = rSges(2,1) * t541 - rSges(2,2) * t538;
t524 = rSges(2,1) * t538 + rSges(2,2) * t541;
t510 = t571 * rSges(3,3) + (rSges(3,1) * t537 + rSges(3,2) * t540) * t534;
t509 = Icges(3,5) * t571 + (Icges(3,1) * t537 + Icges(3,4) * t540) * t534;
t508 = Icges(3,6) * t571 + (Icges(3,4) * t537 + Icges(3,2) * t540) * t534;
t507 = Icges(3,3) * t571 + (Icges(3,5) * t537 + Icges(3,6) * t540) * t534;
t498 = t515 * t539 - t558;
t497 = -t515 * t535 - t539 * t567;
t486 = rSges(3,1) * t519 - rSges(3,2) * t518 + rSges(3,3) * t568;
t485 = rSges(3,1) * t517 - rSges(3,2) * t516 - rSges(3,3) * t566;
t483 = Icges(3,1) * t519 - Icges(3,4) * t518 + Icges(3,5) * t568;
t482 = Icges(3,1) * t517 - Icges(3,4) * t516 - Icges(3,5) * t566;
t481 = Icges(3,4) * t519 - Icges(3,2) * t518 + Icges(3,6) * t568;
t480 = Icges(3,4) * t517 - Icges(3,2) * t516 - Icges(3,6) * t566;
t477 = rSges(4,1) * t515 - rSges(4,2) * t514 - rSges(4,3) * t567;
t476 = Icges(4,1) * t515 - Icges(4,4) * t514 - Icges(4,5) * t567;
t475 = Icges(4,4) * t515 - Icges(4,2) * t514 - Icges(4,6) * t567;
t474 = Icges(4,5) * t515 - Icges(4,6) * t514 - Icges(4,3) * t567;
t473 = qJD(5) * t514 + t492;
t472 = t502 * t539 + t569;
t471 = -t502 * t535 + t518 * t539;
t470 = t500 * t539 + t570;
t469 = -t500 * t535 + t516 * t539;
t457 = rSges(4,1) * t502 - rSges(4,2) * t501 + rSges(4,3) * t518;
t456 = rSges(4,1) * t500 - rSges(4,2) * t499 + rSges(4,3) * t516;
t455 = Icges(4,1) * t502 - Icges(4,4) * t501 + Icges(4,5) * t518;
t454 = Icges(4,1) * t500 - Icges(4,4) * t499 + Icges(4,5) * t516;
t453 = Icges(4,4) * t502 - Icges(4,2) * t501 + Icges(4,6) * t518;
t452 = Icges(4,4) * t500 - Icges(4,2) * t499 + Icges(4,6) * t516;
t451 = Icges(4,5) * t502 - Icges(4,6) * t501 + Icges(4,3) * t518;
t450 = Icges(4,5) * t500 - Icges(4,6) * t499 + Icges(4,3) * t516;
t449 = rSges(5,1) * t498 + rSges(5,2) * t497 + rSges(5,3) * t514;
t448 = Icges(5,1) * t498 + Icges(5,4) * t497 + Icges(5,5) * t514;
t447 = Icges(5,4) * t498 + Icges(5,2) * t497 + Icges(5,6) * t514;
t446 = Icges(5,5) * t498 + Icges(5,6) * t497 + Icges(5,3) * t514;
t444 = rSges(6,1) * t494 - rSges(6,2) * t493 + rSges(6,3) * t514;
t442 = qJD(5) * t499 + t464;
t441 = qJD(5) * t501 + t463;
t433 = t486 * t529 - t510 * t528 + t522;
t432 = -t485 * t529 - t510 * t556 - t560;
t431 = (t485 * t538 + t486 * t541) * t559;
t428 = rSges(5,1) * t472 + rSges(5,2) * t471 + rSges(5,3) * t501;
t427 = rSges(5,1) * t470 + rSges(5,2) * t469 + rSges(5,3) * t499;
t426 = Icges(5,1) * t472 + Icges(5,4) * t471 + Icges(5,5) * t501;
t425 = Icges(5,1) * t470 + Icges(5,4) * t469 + Icges(5,5) * t499;
t424 = Icges(5,4) * t472 + Icges(5,2) * t471 + Icges(5,6) * t501;
t423 = Icges(5,4) * t470 + Icges(5,2) * t469 + Icges(5,6) * t499;
t422 = Icges(5,5) * t472 + Icges(5,6) * t471 + Icges(5,3) * t501;
t421 = Icges(5,5) * t470 + Icges(5,6) * t469 + Icges(5,3) * t499;
t419 = rSges(6,1) * t468 - rSges(6,2) * t467 + rSges(6,3) * t501;
t417 = rSges(6,1) * t466 - rSges(6,2) * t465 + rSges(6,3) * t499;
t399 = t457 * t521 - t477 * t503 + t550;
t398 = -t456 * t521 + t477 * t504 + t548;
t397 = t456 * t503 - t457 * t504 + t561;
t396 = t428 * t492 - t449 * t463 + t547;
t395 = -t427 * t492 + t449 * t464 + t546;
t394 = t427 * t463 - t428 * t464 + t551;
t393 = t419 * t473 - t441 * t444 + t545;
t392 = -t417 * t473 + t442 * t444 + t544;
t391 = t417 * t441 - t419 * t442 + t549;
t390 = qJD(6) * t465 - t441 * t562 + t473 * t563 + t545;
t389 = qJD(6) * t467 + t442 * t562 - t473 * t564 + t544;
t388 = qJD(6) * t493 + t441 * t564 - t442 * t563 + t549;
t1 = m(7) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + t529 * ((t571 * t479 + (t481 * t540 + t483 * t537) * t534) * t528 - (t571 * t478 + (t480 * t540 + t482 * t537) * t534) * t556 + (t571 * t507 + (t508 * t540 + t509 * t537) * t534) * t529) / 0.2e1 + t503 * ((t451 * t518 - t453 * t501 + t455 * t502) * t503 + (t450 * t518 - t452 * t501 + t454 * t502) * t504 + (t474 * t518 - t475 * t501 + t476 * t502) * t521) / 0.2e1 + t504 * ((t451 * t516 - t453 * t499 + t455 * t500) * t503 + (t450 * t516 - t452 * t499 + t454 * t500) * t504 + (t474 * t516 - t475 * t499 + t476 * t500) * t521) / 0.2e1 + t521 * ((-t451 * t567 - t453 * t514 + t455 * t515) * t503 + (-t450 * t567 - t452 * t514 + t454 * t515) * t504 + (-t474 * t567 - t514 * t475 + t515 * t476) * t521) / 0.2e1 + t463 * ((t501 * t422 + t471 * t424 + t472 * t426) * t463 + (t421 * t501 + t423 * t471 + t425 * t472) * t464 + (t446 * t501 + t447 * t471 + t448 * t472) * t492) / 0.2e1 + t464 * ((t422 * t499 + t424 * t469 + t426 * t470) * t463 + (t499 * t421 + t469 * t423 + t470 * t425) * t464 + (t446 * t499 + t447 * t469 + t448 * t470) * t492) / 0.2e1 + t492 * ((t422 * t514 + t424 * t497 + t426 * t498) * t463 + (t421 * t514 + t423 * t497 + t425 * t498) * t464 + (t514 * t446 + t497 * t447 + t498 * t448) * t492) / 0.2e1 + m(3) * (t431 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + m(4) * (t397 ^ 2 + t398 ^ 2 + t399 ^ 2) / 0.2e1 + m(5) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(6) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 - ((-t507 * t566 - t508 * t516 + t509 * t517) * t529 + ((-t481 * t516 + t483 * t517) * t538 + (t516 * t480 - t517 * t482 + t588) * t541) * t559) * t556 / 0.2e1 + ((t507 * t568 - t508 * t518 + t509 * t519) * t529 + (-(-t480 * t518 + t482 * t519) * t541 + (-t518 * t481 + t519 * t483 - t588) * t538) * t559) * t528 / 0.2e1 + ((t467 * t581 + t468 * t579 + t501 * t580) * t473 + (t467 * t587 + t468 * t583 + t501 * t585) * t442 + (t586 * t467 + t582 * t468 + t584 * t501) * t441) * t441 / 0.2e1 + ((t465 * t581 + t466 * t579 + t499 * t580) * t473 + (t587 * t465 + t583 * t466 + t585 * t499) * t442 + (t465 * t586 + t466 * t582 + t499 * t584) * t441) * t442 / 0.2e1 + ((t581 * t493 + t579 * t494 + t580 * t514) * t473 + (t493 * t587 + t494 * t583 + t514 * t585) * t442 + (t493 * t586 + t494 * t582 + t514 * t584) * t441) * t473 / 0.2e1 + (Icges(2,3) + m(2) * (t524 ^ 2 + t525 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
