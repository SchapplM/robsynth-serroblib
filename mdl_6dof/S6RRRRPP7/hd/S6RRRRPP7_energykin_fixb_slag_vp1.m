% Calculate kinetic energy for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:06
% EndTime: 2019-03-09 21:17:09
% DurationCPUTime: 2.73s
% Computational Cost: add. (3010->326), mult. (6560->479), div. (0->0), fcn. (8055->12), ass. (0->150)
t590 = Icges(6,1) + Icges(7,1);
t589 = -Icges(6,4) + Icges(7,5);
t588 = Icges(7,4) + Icges(6,5);
t587 = Icges(6,2) + Icges(7,3);
t586 = -Icges(6,6) + Icges(7,6);
t585 = Icges(7,2) + Icges(5,3) + Icges(6,3);
t584 = rSges(7,1) + pkin(5);
t583 = rSges(7,3) + qJ(6);
t532 = sin(qJ(2));
t533 = sin(qJ(1));
t535 = cos(qJ(2));
t536 = cos(qJ(1));
t565 = cos(pkin(6));
t548 = t536 * t565;
t510 = t532 * t533 - t535 * t548;
t511 = t532 * t548 + t533 * t535;
t528 = sin(pkin(6));
t560 = t528 * t536;
t472 = Icges(3,5) * t511 - Icges(3,6) * t510 - Icges(3,3) * t560;
t549 = t533 * t565;
t512 = t536 * t532 + t535 * t549;
t513 = -t532 * t549 + t536 * t535;
t562 = t528 * t533;
t473 = Icges(3,5) * t513 - Icges(3,6) * t512 + Icges(3,3) * t562;
t582 = (t472 * t536 - t473 * t533) * t528;
t531 = sin(qJ(3));
t568 = cos(qJ(3));
t494 = t511 * t568 - t531 * t560;
t553 = qJ(4) + pkin(11);
t527 = sin(t553);
t547 = cos(t553);
t460 = t494 * t527 - t510 * t547;
t461 = t494 * t547 + t510 * t527;
t551 = t528 * t568;
t493 = t511 * t531 + t536 * t551;
t581 = t587 * t460 + t589 * t461 + t586 * t493;
t496 = t513 * t568 + t531 * t562;
t462 = t496 * t527 - t512 * t547;
t463 = t496 * t547 + t512 * t527;
t495 = t513 * t531 - t533 * t551;
t580 = t587 * t462 + t589 * t463 + t586 * t495;
t579 = t589 * t460 + t590 * t461 + t588 * t493;
t578 = t589 * t462 + t590 * t463 + t588 * t495;
t509 = t531 * t565 + t532 * t551;
t561 = t528 * t535;
t486 = t509 * t527 + t547 * t561;
t487 = t509 * t547 - t527 * t561;
t508 = t528 * t531 * t532 - t565 * t568;
t577 = t587 * t486 + t589 * t487 + t586 * t508;
t576 = t589 * t486 + t590 * t487 + t588 * t508;
t530 = sin(qJ(4));
t534 = cos(qJ(4));
t464 = -t494 * t530 + t510 * t534;
t564 = t510 * t530;
t465 = t494 * t534 + t564;
t575 = Icges(5,5) * t465 + Icges(5,6) * t464 + t586 * t460 + t588 * t461 + t585 * t493;
t466 = -t496 * t530 + t512 * t534;
t563 = t512 * t530;
t467 = t496 * t534 + t563;
t574 = Icges(5,5) * t467 + Icges(5,6) * t466 + t586 * t462 + t588 * t463 + t585 * t495;
t491 = -t509 * t530 - t534 * t561;
t552 = t530 * t561;
t492 = t509 * t534 - t552;
t573 = Icges(5,5) * t492 + Icges(5,6) * t491 + t586 * t486 + t588 * t487 + t585 * t508;
t567 = pkin(4) * t534;
t559 = rSges(7,2) * t493 + t583 * t460 + t584 * t461;
t558 = rSges(7,2) * t495 + t583 * t462 + t584 * t463;
t557 = rSges(7,2) * t508 + t583 * t486 + t584 * t487;
t484 = pkin(2) * t511 + pkin(9) * t510;
t485 = pkin(2) * t513 + pkin(9) * t512;
t554 = qJD(2) * t528;
t522 = t533 * t554;
t550 = t536 * t554;
t556 = t484 * t522 + t485 * t550;
t497 = qJD(3) * t512 + t522;
t555 = qJD(1) * (pkin(1) * t533 - pkin(8) * t560);
t523 = qJD(2) * t565 + qJD(1);
t498 = qJD(3) * t510 - t550;
t456 = pkin(3) * t494 + pkin(10) * t493;
t457 = pkin(3) * t496 + pkin(10) * t495;
t545 = t497 * t456 - t457 * t498 + t556;
t515 = -qJD(3) * t561 + t523;
t514 = (pkin(2) * t532 - pkin(9) * t535) * t528;
t516 = qJD(1) * (pkin(1) * t536 + pkin(8) * t562);
t544 = t523 * t485 - t514 * t522 + t516;
t399 = pkin(4) * t564 + qJ(5) * t493 + t494 * t567;
t458 = qJD(4) * t495 + t497;
t543 = qJD(5) * t508 + t458 * t399 + t545;
t542 = -t484 * t523 - t514 * t550 - t555;
t483 = pkin(3) * t509 + pkin(10) * t508;
t541 = t515 * t457 - t483 * t497 + t544;
t400 = pkin(4) * t563 + qJ(5) * t495 + t496 * t567;
t488 = qJD(4) * t508 + t515;
t540 = qJD(5) * t493 + t488 * t400 + t541;
t539 = -t456 * t515 + t498 * t483 + t542;
t440 = -pkin(4) * t552 + qJ(5) * t508 + t509 * t567;
t459 = qJD(4) * t493 + t498;
t538 = qJD(5) * t495 + t459 * t440 + t539;
t519 = rSges(2,1) * t536 - rSges(2,2) * t533;
t518 = rSges(2,1) * t533 + rSges(2,2) * t536;
t504 = t565 * rSges(3,3) + (rSges(3,1) * t532 + rSges(3,2) * t535) * t528;
t503 = Icges(3,5) * t565 + (Icges(3,1) * t532 + Icges(3,4) * t535) * t528;
t502 = Icges(3,6) * t565 + (Icges(3,4) * t532 + Icges(3,2) * t535) * t528;
t501 = Icges(3,3) * t565 + (Icges(3,5) * t532 + Icges(3,6) * t535) * t528;
t480 = rSges(3,1) * t513 - rSges(3,2) * t512 + rSges(3,3) * t562;
t479 = rSges(3,1) * t511 - rSges(3,2) * t510 - rSges(3,3) * t560;
t477 = Icges(3,1) * t513 - Icges(3,4) * t512 + Icges(3,5) * t562;
t476 = Icges(3,1) * t511 - Icges(3,4) * t510 - Icges(3,5) * t560;
t475 = Icges(3,4) * t513 - Icges(3,2) * t512 + Icges(3,6) * t562;
t474 = Icges(3,4) * t511 - Icges(3,2) * t510 - Icges(3,6) * t560;
t471 = rSges(4,1) * t509 - rSges(4,2) * t508 - rSges(4,3) * t561;
t470 = Icges(4,1) * t509 - Icges(4,4) * t508 - Icges(4,5) * t561;
t469 = Icges(4,4) * t509 - Icges(4,2) * t508 - Icges(4,6) * t561;
t468 = Icges(4,5) * t509 - Icges(4,6) * t508 - Icges(4,3) * t561;
t452 = rSges(4,1) * t496 - rSges(4,2) * t495 + rSges(4,3) * t512;
t451 = rSges(4,1) * t494 - rSges(4,2) * t493 + rSges(4,3) * t510;
t450 = Icges(4,1) * t496 - Icges(4,4) * t495 + Icges(4,5) * t512;
t449 = Icges(4,1) * t494 - Icges(4,4) * t493 + Icges(4,5) * t510;
t448 = Icges(4,4) * t496 - Icges(4,2) * t495 + Icges(4,6) * t512;
t447 = Icges(4,4) * t494 - Icges(4,2) * t493 + Icges(4,6) * t510;
t446 = Icges(4,5) * t496 - Icges(4,6) * t495 + Icges(4,3) * t512;
t445 = Icges(4,5) * t494 - Icges(4,6) * t493 + Icges(4,3) * t510;
t444 = rSges(5,1) * t492 + rSges(5,2) * t491 + rSges(5,3) * t508;
t443 = Icges(5,1) * t492 + Icges(5,4) * t491 + Icges(5,5) * t508;
t442 = Icges(5,4) * t492 + Icges(5,2) * t491 + Icges(5,6) * t508;
t438 = rSges(6,1) * t487 - rSges(6,2) * t486 + rSges(6,3) * t508;
t430 = t480 * t523 - t504 * t522 + t516;
t429 = -t479 * t523 - t504 * t550 - t555;
t428 = (t479 * t533 + t480 * t536) * t554;
t425 = rSges(5,1) * t467 + rSges(5,2) * t466 + rSges(5,3) * t495;
t424 = rSges(5,1) * t465 + rSges(5,2) * t464 + rSges(5,3) * t493;
t423 = Icges(5,1) * t467 + Icges(5,4) * t466 + Icges(5,5) * t495;
t422 = Icges(5,1) * t465 + Icges(5,4) * t464 + Icges(5,5) * t493;
t421 = Icges(5,4) * t467 + Icges(5,2) * t466 + Icges(5,6) * t495;
t420 = Icges(5,4) * t465 + Icges(5,2) * t464 + Icges(5,6) * t493;
t416 = rSges(6,1) * t463 - rSges(6,2) * t462 + rSges(6,3) * t495;
t414 = rSges(6,1) * t461 - rSges(6,2) * t460 + rSges(6,3) * t493;
t396 = t452 * t515 - t471 * t497 + t544;
t395 = -t451 * t515 + t471 * t498 + t542;
t394 = t451 * t497 - t452 * t498 + t556;
t393 = t425 * t488 - t444 * t458 + t541;
t392 = -t424 * t488 + t444 * t459 + t539;
t391 = t424 * t458 - t425 * t459 + t545;
t390 = t416 * t488 + (-t438 - t440) * t458 + t540;
t389 = t438 * t459 + (-t399 - t414) * t488 + t538;
t388 = t414 * t458 + (-t400 - t416) * t459 + t543;
t387 = qJD(6) * t460 + t558 * t488 + (-t440 - t557) * t458 + t540;
t386 = qJD(6) * t462 + t557 * t459 + (-t399 - t559) * t488 + t538;
t385 = qJD(6) * t486 + t559 * t458 + (-t400 - t558) * t459 + t543;
t1 = ((t501 * t562 - t502 * t512 + t503 * t513) * t523 + (-(-t474 * t512 + t476 * t513) * t536 + (-t512 * t475 + t513 * t477 - t582) * t533) * t554) * t522 / 0.2e1 - ((-t501 * t560 - t502 * t510 + t503 * t511) * t523 + ((-t475 * t510 + t477 * t511) * t533 + (t510 * t474 - t511 * t476 + t582) * t536) * t554) * t550 / 0.2e1 + t497 * ((t512 * t446 - t495 * t448 + t496 * t450) * t497 + (t445 * t512 - t447 * t495 + t449 * t496) * t498 + (t468 * t512 - t469 * t495 + t470 * t496) * t515) / 0.2e1 + t515 * ((-t446 * t561 - t448 * t508 + t450 * t509) * t497 + (-t445 * t561 - t447 * t508 + t449 * t509) * t498 + (-t468 * t561 - t508 * t469 + t509 * t470) * t515) / 0.2e1 + t498 * ((t446 * t510 - t448 * t493 + t450 * t494) * t497 + (t510 * t445 - t493 * t447 + t494 * t449) * t498 + (t468 * t510 - t469 * t493 + t470 * t494) * t515) / 0.2e1 + t523 * ((t565 * t473 + (t475 * t535 + t477 * t532) * t528) * t522 - (t565 * t472 + (t474 * t535 + t476 * t532) * t528) * t550 + (t565 * t501 + (t502 * t535 + t503 * t532) * t528) * t523) / 0.2e1 + m(5) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + m(4) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(3) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(6) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(7) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t518 ^ 2 + t519 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t442 * t466 + t443 * t467 + t462 * t577 + t463 * t576 + t495 * t573) * t488 + (t420 * t466 + t422 * t467 + t462 * t581 + t463 * t579 + t495 * t575) * t459 + (t466 * t421 + t467 * t423 + t580 * t462 + t578 * t463 + t574 * t495) * t458) * t458 / 0.2e1 + ((t442 * t464 + t443 * t465 + t460 * t577 + t461 * t576 + t493 * t573) * t488 + (t464 * t420 + t465 * t422 + t581 * t460 + t579 * t461 + t575 * t493) * t459 + (t421 * t464 + t423 * t465 + t460 * t580 + t461 * t578 + t493 * t574) * t458) * t459 / 0.2e1 + ((t491 * t442 + t492 * t443 + t577 * t486 + t576 * t487 + t573 * t508) * t488 + (t420 * t491 + t422 * t492 + t486 * t581 + t487 * t579 + t508 * t575) * t459 + (t421 * t491 + t423 * t492 + t486 * t580 + t487 * t578 + t508 * t574) * t458) * t488 / 0.2e1;
T  = t1;
