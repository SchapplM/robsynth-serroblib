% Calculate kinetic energy for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:38
% EndTime: 2019-03-08 19:46:41
% DurationCPUTime: 2.92s
% Computational Cost: add. (2130->326), mult. (5021->485), div. (0->0), fcn. (6032->12), ass. (0->151)
t572 = Icges(3,1) + Icges(4,2);
t571 = Icges(3,4) + Icges(4,6);
t570 = Icges(3,5) - Icges(4,4);
t569 = Icges(3,2) + Icges(4,3);
t568 = Icges(5,2) + Icges(6,3);
t567 = Icges(3,6) - Icges(4,5);
t566 = Icges(3,3) + Icges(4,1);
t513 = cos(pkin(6));
t510 = sin(pkin(6));
t512 = cos(pkin(10));
t542 = t510 * t512;
t509 = sin(pkin(10));
t543 = t509 * t510;
t516 = sin(qJ(2));
t517 = cos(qJ(2));
t553 = t566 * t513 + (t570 * t516 + t567 * t517) * t510;
t538 = t513 * t517;
t490 = t509 * t516 - t512 * t538;
t539 = t513 * t516;
t491 = t509 * t517 + t512 * t539;
t556 = -t567 * t490 + t570 * t491 - t566 * t542;
t492 = t509 * t538 + t512 * t516;
t493 = -t509 * t539 + t512 * t517;
t557 = -t567 * t492 + t570 * t493 + t566 * t543;
t565 = -t553 * t513 + t542 * t556 - t543 * t557;
t515 = sin(qJ(4));
t547 = cos(qJ(4));
t531 = t510 * t547;
t467 = t492 * t515 + t509 * t531;
t508 = sin(pkin(11));
t511 = cos(pkin(11));
t428 = -t467 * t508 + t493 * t511;
t544 = t493 * t508;
t429 = t467 * t511 + t544;
t541 = t510 * t515;
t466 = -t492 * t547 + t509 * t541;
t564 = -Icges(5,4) * t467 + Icges(6,5) * t429 - Icges(5,6) * t493 + Icges(6,6) * t428 + t568 * t466;
t469 = t490 * t515 - t512 * t531;
t430 = -t469 * t508 + t491 * t511;
t545 = t491 * t508;
t431 = t469 * t511 + t545;
t468 = t490 * t547 + t512 * t541;
t563 = -Icges(5,4) * t469 + Icges(6,5) * t431 - Icges(5,6) * t491 + Icges(6,6) * t430 - t568 * t468;
t495 = t513 * t547 - t517 * t541;
t540 = t510 * t516;
t470 = -t495 * t508 + t511 * t540;
t532 = t508 * t540;
t471 = t495 * t511 + t532;
t494 = t513 * t515 + t517 * t531;
t562 = -Icges(5,4) * t495 + Icges(6,5) * t471 - Icges(5,6) * t540 + Icges(6,6) * t470 + t568 * t494;
t561 = t569 * t492 - t571 * t493 - t567 * t543;
t560 = t569 * t490 - t571 * t491 + t567 * t542;
t559 = -t571 * t492 + t572 * t493 + t570 * t543;
t558 = t571 * t490 - t572 * t491 + t570 * t542;
t555 = t567 * t513 + (t571 * t516 + t569 * t517) * t510;
t554 = t570 * t513 + (t572 * t516 + t571 * t517) * t510;
t551 = qJD(2) ^ 2;
t546 = pkin(5) * t511;
t456 = pkin(2) * t493 + qJ(3) * t492;
t504 = qJD(2) * t513;
t536 = qJD(3) * t490 + t456 * t504;
t535 = qJD(2) * t510;
t501 = t509 * t535;
t472 = qJD(4) * t493 + t501;
t497 = qJD(4) * t540 + t504;
t534 = qJD(3) * t517;
t530 = t512 * t535;
t455 = pkin(2) * t491 + qJ(3) * t490;
t529 = t455 * t501 + t456 * t530 + qJD(1);
t496 = (pkin(2) * t516 - qJ(3) * t517) * t510;
t527 = (-t513 * rSges(4,1) - (-rSges(4,2) * t516 - rSges(4,3) * t517) * t510 - t496) * t510;
t526 = (-pkin(3) * t513 - pkin(8) * t540 - t496) * t510;
t473 = qJD(4) * t491 - t530;
t474 = pkin(3) * t543 + pkin(8) * t493;
t475 = -pkin(3) * t542 + pkin(8) * t491;
t523 = t474 * t530 + t475 * t501 - t510 * t534 + t529;
t522 = qJD(2) * t509 * t526 + t474 * t504 + t536;
t421 = pkin(4) * t469 - qJ(5) * t468;
t521 = qJD(5) * t494 + t472 * t421 + t523;
t420 = pkin(4) * t467 + qJ(5) * t466;
t520 = -qJD(5) * t468 + t497 * t420 + t522;
t488 = qJD(3) * t492;
t519 = t488 + ((-t455 - t475) * t513 + t512 * t526) * qJD(2);
t457 = pkin(4) * t495 + qJ(5) * t494;
t518 = qJD(5) * t466 + t473 * t457 + t519;
t507 = pkin(11) + qJ(6);
t506 = cos(t507);
t505 = sin(t507);
t482 = t513 * rSges(3,3) + (rSges(3,1) * t516 + rSges(3,2) * t517) * t510;
t464 = qJD(6) * t494 + t497;
t461 = t495 * t506 + t505 * t540;
t460 = -t495 * t505 + t506 * t540;
t453 = rSges(5,1) * t495 - rSges(5,2) * t494 + rSges(5,3) * t540;
t452 = Icges(5,1) * t495 - Icges(5,4) * t494 + Icges(5,5) * t540;
t450 = Icges(5,5) * t495 - Icges(5,6) * t494 + Icges(5,3) * t540;
t449 = rSges(3,1) * t493 - rSges(3,2) * t492 + rSges(3,3) * t543;
t448 = rSges(3,1) * t491 - rSges(3,2) * t490 - rSges(3,3) * t542;
t447 = -rSges(4,1) * t542 - rSges(4,2) * t491 + rSges(4,3) * t490;
t446 = rSges(4,1) * t543 - rSges(4,2) * t493 + rSges(4,3) * t492;
t427 = t469 * t506 + t491 * t505;
t426 = -t469 * t505 + t491 * t506;
t425 = t467 * t506 + t493 * t505;
t424 = -t467 * t505 + t493 * t506;
t423 = -qJD(6) * t468 + t473;
t422 = qJD(6) * t466 + t472;
t417 = (-t448 * t513 - t482 * t542) * qJD(2);
t416 = (t449 * t513 - t482 * t543) * qJD(2);
t415 = rSges(6,1) * t471 + rSges(6,2) * t470 + rSges(6,3) * t494;
t414 = rSges(5,1) * t469 + rSges(5,2) * t468 + rSges(5,3) * t491;
t413 = rSges(5,1) * t467 - rSges(5,2) * t466 + rSges(5,3) * t493;
t412 = Icges(6,1) * t471 + Icges(6,4) * t470 + Icges(6,5) * t494;
t411 = Icges(6,4) * t471 + Icges(6,2) * t470 + Icges(6,6) * t494;
t409 = Icges(5,1) * t469 + Icges(5,4) * t468 + Icges(5,5) * t491;
t408 = Icges(5,1) * t467 - Icges(5,4) * t466 + Icges(5,5) * t493;
t405 = Icges(5,5) * t469 + Icges(5,6) * t468 + Icges(5,3) * t491;
t404 = Icges(5,5) * t467 - Icges(5,6) * t466 + Icges(5,3) * t493;
t403 = pkin(5) * t532 + pkin(9) * t494 + t495 * t546;
t402 = rSges(7,1) * t461 + rSges(7,2) * t460 + rSges(7,3) * t494;
t401 = Icges(7,1) * t461 + Icges(7,4) * t460 + Icges(7,5) * t494;
t400 = Icges(7,4) * t461 + Icges(7,2) * t460 + Icges(7,6) * t494;
t399 = Icges(7,5) * t461 + Icges(7,6) * t460 + Icges(7,3) * t494;
t397 = qJD(1) + (t448 * t509 + t449 * t512) * t535;
t396 = rSges(6,1) * t431 + rSges(6,2) * t430 - rSges(6,3) * t468;
t395 = rSges(6,1) * t429 + rSges(6,2) * t428 + rSges(6,3) * t466;
t394 = Icges(6,1) * t431 + Icges(6,4) * t430 - Icges(6,5) * t468;
t393 = Icges(6,1) * t429 + Icges(6,4) * t428 + Icges(6,5) * t466;
t392 = Icges(6,4) * t431 + Icges(6,2) * t430 - Icges(6,6) * t468;
t391 = Icges(6,4) * t429 + Icges(6,2) * t428 + Icges(6,6) * t466;
t388 = rSges(7,1) * t427 + rSges(7,2) * t426 - rSges(7,3) * t468;
t387 = rSges(7,1) * t425 + rSges(7,2) * t424 + rSges(7,3) * t466;
t386 = Icges(7,1) * t427 + Icges(7,4) * t426 - Icges(7,5) * t468;
t385 = Icges(7,1) * t425 + Icges(7,4) * t424 + Icges(7,5) * t466;
t384 = Icges(7,4) * t427 + Icges(7,2) * t426 - Icges(7,6) * t468;
t383 = Icges(7,4) * t425 + Icges(7,2) * t424 + Icges(7,6) * t466;
t382 = Icges(7,5) * t427 + Icges(7,6) * t426 - Icges(7,3) * t468;
t381 = Icges(7,5) * t425 + Icges(7,6) * t424 + Icges(7,3) * t466;
t380 = pkin(5) * t545 - pkin(9) * t468 + t469 * t546;
t379 = pkin(5) * t544 + pkin(9) * t466 + t467 * t546;
t378 = t488 + ((-t447 - t455) * t513 + t512 * t527) * qJD(2);
t377 = (t446 * t513 + t509 * t527) * qJD(2) + t536;
t376 = (-t534 + (t446 * t512 + t447 * t509) * qJD(2)) * t510 + t529;
t375 = -t414 * t497 + t453 * t473 + t519;
t374 = t413 * t497 - t453 * t472 + t522;
t373 = -t473 * t413 + t472 * t414 + t523;
t372 = t415 * t473 + (-t396 - t421) * t497 + t518;
t371 = t395 * t497 + (-t415 - t457) * t472 + t520;
t370 = t472 * t396 + (-t395 - t420) * t473 + t521;
t369 = -t388 * t464 + t402 * t423 + t403 * t473 + (-t380 - t421) * t497 + t518;
t368 = t379 * t497 + t387 * t464 - t402 * t422 + (-t403 - t457) * t472 + t520;
t367 = t472 * t380 - t423 * t387 + t422 * t388 + (-t379 - t420) * t473 + t521;
t1 = m(7) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(5) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(6) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(4) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(3) * (t397 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t422 * ((t466 * t381 + t424 * t383 + t425 * t385) * t422 + (t382 * t466 + t384 * t424 + t386 * t425) * t423 + (t399 * t466 + t400 * t424 + t401 * t425) * t464) / 0.2e1 + t423 * ((-t381 * t468 + t383 * t426 + t385 * t427) * t422 + (-t382 * t468 + t384 * t426 + t386 * t427) * t423 + (-t399 * t468 + t400 * t426 + t401 * t427) * t464) / 0.2e1 + t464 * ((t381 * t494 + t383 * t460 + t385 * t461) * t422 + (t382 * t494 + t384 * t460 + t386 * t461) * t423 + (t399 * t494 + t400 * t460 + t401 * t461) * t464) / 0.2e1 + ((t411 * t428 + t412 * t429 + t450 * t493 + t452 * t467 + t466 * t562) * t497 + (t392 * t428 + t394 * t429 + t405 * t493 + t409 * t467 + t466 * t563) * t473 + (t391 * t428 + t393 * t429 + t404 * t493 + t408 * t467 + t564 * t466) * t472) * t472 / 0.2e1 + ((t411 * t430 + t412 * t431 + t450 * t491 + t452 * t469 - t468 * t562) * t497 + (t392 * t430 + t394 * t431 + t405 * t491 + t409 * t469 - t563 * t468) * t473 + (t391 * t430 + t393 * t431 + t404 * t491 + t408 * t469 - t468 * t564) * t472) * t473 / 0.2e1 + ((t411 * t470 + t412 * t471 + t450 * t540 + t452 * t495 + t562 * t494) * t497 + (t392 * t470 + t394 * t471 + t405 * t540 + t409 * t495 + t494 * t563) * t473 + (t391 * t470 + t393 * t471 + t404 * t540 + t408 * t495 + t494 * t564) * t472) * t497 / 0.2e1 - ((t490 * t561 + t491 * t559) * t543 + (-t490 * t555 + t491 * t554) * t513 + (-t490 * t560 + t491 * t558 + t565) * t542) * t551 * t542 / 0.2e1 + ((t553 * t513 ^ 2 + (((t516 * t558 + t517 * t560) * t512 + (t516 * t559 - t517 * t561) * t509) * t510 + (t509 * t557 - t512 * t556 + t516 * t554 + t517 * t555) * t513) * t510) * t513 + ((-t492 * t560 + t493 * t558) * t542 + (-t492 * t555 + t493 * t554) * t513 + (t492 * t561 + t493 * t559 - t565) * t543) * t543) * t551 / 0.2e1;
T  = t1;
