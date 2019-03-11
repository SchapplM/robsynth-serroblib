% Calculate kinetic energy for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:31
% EndTime: 2019-03-09 12:36:33
% DurationCPUTime: 2.97s
% Computational Cost: add. (3130->331), mult. (5670->495), div. (0->0), fcn. (6847->12), ass. (0->158)
t601 = Icges(6,1) + Icges(7,1);
t600 = -Icges(6,4) + Icges(7,5);
t599 = Icges(7,4) + Icges(6,5);
t598 = Icges(6,2) + Icges(7,3);
t597 = Icges(7,2) + Icges(6,3);
t596 = -Icges(6,6) + Icges(7,6);
t595 = rSges(7,1) + pkin(5);
t594 = rSges(7,3) + qJ(6);
t534 = sin(qJ(2));
t535 = sin(qJ(1));
t536 = cos(qJ(2));
t537 = cos(qJ(1));
t573 = cos(pkin(6));
t554 = t537 * t573;
t511 = t534 * t535 - t536 * t554;
t512 = t534 * t554 + t535 * t536;
t530 = sin(pkin(6));
t569 = t530 * t537;
t471 = Icges(3,5) * t512 - Icges(3,6) * t511 - Icges(3,3) * t569;
t555 = t535 * t573;
t513 = t537 * t534 + t536 * t555;
t514 = -t534 * t555 + t537 * t536;
t571 = t530 * t535;
t472 = Icges(3,5) * t514 - Icges(3,6) * t513 + Icges(3,3) * t571;
t593 = t530 * (t471 * t537 - t472 * t535);
t560 = pkin(11) + qJ(4);
t528 = sin(t560);
t552 = cos(t560);
t488 = t512 * t552 - t528 * t569;
t533 = sin(qJ(5));
t575 = cos(qJ(5));
t457 = t488 * t533 - t511 * t575;
t458 = t488 * t575 + t511 * t533;
t546 = t530 * t552;
t487 = t512 * t528 + t537 * t546;
t592 = t598 * t457 + t600 * t458 + t596 * t487;
t490 = t514 * t552 + t528 * t571;
t459 = t490 * t533 - t513 * t575;
t460 = t490 * t575 + t513 * t533;
t489 = t514 * t528 - t535 * t546;
t591 = t598 * t459 + t600 * t460 + t596 * t489;
t590 = t596 * t457 + t599 * t458 + t597 * t487;
t589 = t596 * t459 + t599 * t460 + t597 * t489;
t588 = t600 * t457 + t601 * t458 + t599 * t487;
t587 = t600 * t459 + t601 * t460 + t599 * t489;
t501 = t528 * t573 + t534 * t546;
t570 = t530 * t536;
t485 = t501 * t533 + t570 * t575;
t486 = t501 * t575 - t533 * t570;
t572 = t530 * t534;
t500 = t528 * t572 - t552 * t573;
t586 = t598 * t485 + t600 * t486 + t596 * t500;
t585 = t596 * t485 + t599 * t486 + t597 * t500;
t584 = t600 * t485 + t601 * t486 + t599 * t500;
t529 = sin(pkin(11));
t531 = cos(pkin(11));
t491 = -t512 * t529 - t531 * t569;
t558 = t529 * t569;
t492 = t512 * t531 - t558;
t443 = Icges(4,5) * t492 + Icges(4,6) * t491 + Icges(4,3) * t511;
t473 = Icges(3,4) * t512 - Icges(3,2) * t511 - Icges(3,6) * t569;
t583 = -t443 + t473;
t493 = -t514 * t529 + t531 * t571;
t559 = t529 * t571;
t494 = t514 * t531 + t559;
t444 = Icges(4,5) * t494 + Icges(4,6) * t493 + Icges(4,3) * t513;
t474 = Icges(3,4) * t514 - Icges(3,2) * t513 + Icges(3,6) * t571;
t582 = t444 - t474;
t509 = -t529 * t572 + t531 * t573;
t553 = t573 * t529;
t510 = t531 * t572 + t553;
t465 = Icges(4,5) * t510 + Icges(4,6) * t509 - Icges(4,3) * t570;
t498 = Icges(3,6) * t573 + (Icges(3,4) * t534 + Icges(3,2) * t536) * t530;
t581 = t465 - t498;
t574 = pkin(3) * t531;
t567 = rSges(7,2) * t487 + t594 * t457 + t595 * t458;
t566 = rSges(7,2) * t489 + t594 * t459 + t595 * t460;
t565 = rSges(7,2) * t500 + t594 * t485 + t595 * t486;
t483 = pkin(2) * t512 + qJ(3) * t511;
t484 = pkin(2) * t514 + qJ(3) * t513;
t562 = qJD(2) * t530;
t524 = t535 * t562;
t556 = t537 * t562;
t564 = t483 * t524 + t484 * t556;
t495 = qJD(4) * t513 + t524;
t563 = qJD(1) * (pkin(1) * t535 - pkin(8) * t569);
t561 = qJD(3) * t536;
t525 = qJD(2) * t573 + qJD(1);
t517 = qJD(1) * (pkin(1) * t537 + pkin(8) * t571);
t557 = qJD(3) * t511 + t525 * t484 + t517;
t551 = qJD(3) * t513 - t563;
t515 = (pkin(2) * t534 - qJ(3) * t536) * t530;
t548 = (-rSges(4,1) * t510 - rSges(4,2) * t509 + rSges(4,3) * t570 - t515) * t562;
t547 = (-pkin(3) * t553 - (-pkin(9) * t536 + t534 * t574) * t530 - t515) * t562;
t496 = qJD(4) * t511 - t556;
t516 = -qJD(4) * t570 + t525;
t440 = -pkin(3) * t558 + pkin(9) * t511 + t512 * t574;
t441 = pkin(3) * t559 + pkin(9) * t513 + t514 * t574;
t544 = t440 * t524 + t441 * t556 - t530 * t561 + t564;
t543 = t525 * t441 + t535 * t547 + t557;
t453 = pkin(4) * t488 + pkin(10) * t487;
t454 = pkin(4) * t490 + pkin(10) * t489;
t542 = t495 * t453 - t454 * t496 + t544;
t469 = pkin(4) * t501 + pkin(10) * t500;
t541 = t516 * t454 - t469 * t495 + t543;
t540 = (-t440 - t483) * t525 + t537 * t547 + t551;
t539 = -t453 * t516 + t496 * t469 + t540;
t521 = rSges(2,1) * t537 - rSges(2,2) * t535;
t520 = rSges(2,1) * t535 + rSges(2,2) * t537;
t502 = t573 * rSges(3,3) + (rSges(3,1) * t534 + rSges(3,2) * t536) * t530;
t499 = Icges(3,5) * t573 + (Icges(3,1) * t534 + Icges(3,4) * t536) * t530;
t497 = Icges(3,3) * t573 + (Icges(3,5) * t534 + Icges(3,6) * t536) * t530;
t482 = qJD(5) * t500 + t516;
t481 = rSges(3,1) * t514 - rSges(3,2) * t513 + rSges(3,3) * t571;
t480 = rSges(3,1) * t512 - rSges(3,2) * t511 - rSges(3,3) * t569;
t476 = Icges(3,1) * t514 - Icges(3,4) * t513 + Icges(3,5) * t571;
t475 = Icges(3,1) * t512 - Icges(3,4) * t511 - Icges(3,5) * t569;
t467 = Icges(4,1) * t510 + Icges(4,4) * t509 - Icges(4,5) * t570;
t466 = Icges(4,4) * t510 + Icges(4,2) * t509 - Icges(4,6) * t570;
t464 = rSges(5,1) * t501 - rSges(5,2) * t500 - rSges(5,3) * t570;
t463 = Icges(5,1) * t501 - Icges(5,4) * t500 - Icges(5,5) * t570;
t462 = Icges(5,4) * t501 - Icges(5,2) * t500 - Icges(5,6) * t570;
t461 = Icges(5,5) * t501 - Icges(5,6) * t500 - Icges(5,3) * t570;
t456 = qJD(5) * t487 + t496;
t455 = qJD(5) * t489 + t495;
t450 = rSges(4,1) * t494 + rSges(4,2) * t493 + rSges(4,3) * t513;
t449 = rSges(4,1) * t492 + rSges(4,2) * t491 + rSges(4,3) * t511;
t448 = Icges(4,1) * t494 + Icges(4,4) * t493 + Icges(4,5) * t513;
t447 = Icges(4,1) * t492 + Icges(4,4) * t491 + Icges(4,5) * t511;
t446 = Icges(4,4) * t494 + Icges(4,2) * t493 + Icges(4,6) * t513;
t445 = Icges(4,4) * t492 + Icges(4,2) * t491 + Icges(4,6) * t511;
t439 = rSges(5,1) * t490 - rSges(5,2) * t489 + rSges(5,3) * t513;
t438 = rSges(5,1) * t488 - rSges(5,2) * t487 + rSges(5,3) * t511;
t437 = Icges(5,1) * t490 - Icges(5,4) * t489 + Icges(5,5) * t513;
t436 = Icges(5,1) * t488 - Icges(5,4) * t487 + Icges(5,5) * t511;
t435 = Icges(5,4) * t490 - Icges(5,2) * t489 + Icges(5,6) * t513;
t434 = Icges(5,4) * t488 - Icges(5,2) * t487 + Icges(5,6) * t511;
t433 = Icges(5,5) * t490 - Icges(5,6) * t489 + Icges(5,3) * t513;
t432 = Icges(5,5) * t488 - Icges(5,6) * t487 + Icges(5,3) * t511;
t431 = t481 * t525 - t502 * t524 + t517;
t430 = -t480 * t525 - t502 * t556 - t563;
t429 = rSges(6,1) * t486 - rSges(6,2) * t485 + rSges(6,3) * t500;
t417 = (t480 * t535 + t481 * t537) * t562;
t414 = rSges(6,1) * t460 - rSges(6,2) * t459 + rSges(6,3) * t489;
t412 = rSges(6,1) * t458 - rSges(6,2) * t457 + rSges(6,3) * t487;
t398 = t450 * t525 + t535 * t548 + t557;
t397 = (-t449 - t483) * t525 + t537 * t548 + t551;
t396 = (-t561 + (t449 * t535 + t450 * t537) * qJD(2)) * t530 + t564;
t395 = t439 * t516 - t464 * t495 + t543;
t394 = -t438 * t516 + t464 * t496 + t540;
t393 = t438 * t495 - t439 * t496 + t544;
t392 = t414 * t482 - t429 * t455 + t541;
t391 = -t412 * t482 + t429 * t456 + t539;
t390 = t412 * t455 - t414 * t456 + t542;
t389 = qJD(6) * t457 - t455 * t565 + t482 * t566 + t541;
t388 = qJD(6) * t459 + t456 * t565 - t482 * t567 + t539;
t387 = qJD(6) * t485 + t455 * t567 - t456 * t566 + t542;
t1 = m(6) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(5) * (t393 ^ 2 + t394 ^ 2 + t395 ^ 2) / 0.2e1 + m(4) * (t396 ^ 2 + t397 ^ 2 + t398 ^ 2) / 0.2e1 + m(3) * (t417 ^ 2 + t430 ^ 2 + t431 ^ 2) / 0.2e1 + t495 * ((t433 * t513 - t435 * t489 + t437 * t490) * t495 + (t432 * t513 - t434 * t489 + t436 * t490) * t496 + (t461 * t513 - t462 * t489 + t463 * t490) * t516) / 0.2e1 + t516 * ((-t433 * t570 - t435 * t500 + t437 * t501) * t495 + (-t432 * t570 - t434 * t500 + t436 * t501) * t496 + (-t461 * t570 - t462 * t500 + t463 * t501) * t516) / 0.2e1 + t496 * ((t433 * t511 - t435 * t487 + t437 * t488) * t495 + (t432 * t511 - t434 * t487 + t436 * t488) * t496 + (t461 * t511 - t462 * t487 + t463 * t488) * t516) / 0.2e1 + m(7) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + ((t459 * t586 + t460 * t584 + t489 * t585) * t482 + (t459 * t592 + t460 * t588 + t489 * t590) * t456 + (t591 * t459 + t587 * t460 + t589 * t489) * t455) * t455 / 0.2e1 + ((t457 * t586 + t458 * t584 + t487 * t585) * t482 + (t592 * t457 + t588 * t458 + t590 * t487) * t456 + (t457 * t591 + t458 * t587 + t487 * t589) * t455) * t456 / 0.2e1 + ((t586 * t485 + t584 * t486 + t585 * t500) * t482 + (t485 * t592 + t486 * t588 + t500 * t590) * t456 + (t485 * t591 + t486 * t587 + t500 * t589) * t455) * t482 / 0.2e1 + ((t573 * t472 + (t474 * t536 + t476 * t534) * t530) * t524 - (t573 * t471 + (t473 * t536 + t475 * t534) * t530) * t556 + ((t446 * t509 + t448 * t510) * t535 - (t445 * t509 + t447 * t510) * t537 + (t443 * t537 - t444 * t535) * t570) * t562 + (t573 * t497 + (t498 * t536 + t499 * t534) * t530 - t465 * t570 + t466 * t509 + t467 * t510) * t525) * t525 / 0.2e1 + (m(2) * (t520 ^ 2 + t521 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 - (((-t445 * t491 - t447 * t492 - t512 * t475 + t511 * t583 + t593) * t537 + (t446 * t491 + t448 * t492 + t476 * t512 + t511 * t582) * t535) * t562 + (t466 * t491 + t467 * t492 - t497 * t569 + t499 * t512 + t511 * t581) * t525) * t556 / 0.2e1 + (((-t445 * t493 - t447 * t494 - t475 * t514 + t513 * t583) * t537 + (t446 * t493 + t448 * t494 + t514 * t476 + t513 * t582 - t593) * t535) * t562 + (t466 * t493 + t467 * t494 + t497 * t571 + t499 * t514 + t513 * t581) * t525) * t524 / 0.2e1;
T  = t1;
