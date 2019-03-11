% Calculate kinetic energy for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:12
% EndTime: 2019-03-09 14:35:16
% DurationCPUTime: 3.17s
% Computational Cost: add. (2373->340), mult. (4637->523), div. (0->0), fcn. (5426->12), ass. (0->160)
t577 = Icges(3,1) + Icges(4,2);
t576 = Icges(3,4) + Icges(4,6);
t575 = Icges(3,5) - Icges(4,4);
t574 = Icges(3,2) + Icges(4,3);
t573 = Icges(3,6) - Icges(4,5);
t572 = Icges(3,3) + Icges(4,1);
t514 = sin(pkin(6));
t515 = cos(pkin(6));
t522 = cos(qJ(2));
t523 = cos(qJ(1));
t548 = t522 * t523;
t518 = sin(qJ(2));
t519 = sin(qJ(1));
t551 = t518 * t519;
t494 = -t515 * t548 + t551;
t549 = t519 * t522;
t550 = t518 * t523;
t495 = t515 * t550 + t549;
t496 = t515 * t549 + t550;
t497 = -t515 * t551 + t548;
t553 = t514 * t523;
t555 = t514 * t519;
t563 = (t573 * t494 - t575 * t495 + t572 * t553) * t523 + (-t573 * t496 + t575 * t497 + t572 * t555) * t519;
t571 = t563 * t514;
t570 = t574 * t496 - t576 * t497 - t573 * t555;
t569 = t574 * t494 - t576 * t495 + t573 * t553;
t568 = -t576 * t496 + t577 * t497 + t575 * t555;
t567 = t576 * t494 - t577 * t495 + t575 * t553;
t566 = t573 * t515 + (t576 * t518 + t574 * t522) * t514;
t565 = t575 * t515 + (t577 * t518 + t576 * t522) * t514;
t564 = t572 * t515 + (t575 * t518 + t573 * t522) * t514;
t521 = cos(qJ(4));
t560 = pkin(4) * t521;
t517 = sin(qJ(4));
t558 = t494 * t517;
t557 = t496 * t517;
t556 = t514 * t518;
t554 = t514 * t522;
t552 = t517 * t522;
t547 = qJ(4) + qJ(5);
t455 = pkin(2) * t495 + qJ(3) * t494;
t456 = pkin(2) * t497 + qJ(3) * t496;
t544 = qJD(2) * t514;
t509 = t519 * t544;
t541 = t523 * t544;
t546 = t455 * t509 + t456 * t541;
t470 = qJD(4) * t497 + t509;
t545 = qJD(1) * (pkin(1) * t519 - pkin(8) * t553);
t543 = qJD(3) * t522;
t510 = qJD(2) * t515 + qJD(1);
t500 = qJD(1) * (pkin(1) * t523 + pkin(8) * t555);
t542 = qJD(3) * t494 + t510 * t456 + t500;
t431 = qJD(5) * t497 + t470;
t499 = qJD(4) * t556 + t510;
t540 = cos(t547);
t539 = qJD(3) * t496 - t545;
t474 = qJD(5) * t556 + t499;
t536 = t514 * t540;
t498 = (pkin(2) * t518 - qJ(3) * t522) * t514;
t535 = (-rSges(4,1) * t515 - (-rSges(4,2) * t518 - rSges(4,3) * t522) * t514 - t498) * t544;
t534 = (-pkin(3) * t515 - pkin(9) * t556 - t498) * t544;
t471 = qJD(4) * t495 - t541;
t432 = qJD(5) * t495 + t471;
t472 = pkin(3) * t555 + pkin(9) * t497;
t473 = -pkin(3) * t553 + pkin(9) * t495;
t531 = t472 * t541 + t473 * t509 - t514 * t543 + t546;
t530 = t510 * t472 + t519 * t534 + t542;
t414 = pkin(4) * t557 + pkin(10) * t497 + t555 * t560;
t415 = pkin(4) * t558 + pkin(10) * t495 - t553 * t560;
t529 = -t414 * t471 + t470 * t415 + t531;
t454 = t560 * t515 + (-pkin(4) * t552 + pkin(10) * t518) * t514;
t528 = t499 * t414 - t454 * t470 + t530;
t527 = (-t455 - t473) * t510 + t523 * t534 + t539;
t526 = -t415 * t499 + t471 * t454 + t527;
t520 = cos(qJ(6));
t516 = sin(qJ(6));
t513 = sin(t547);
t504 = rSges(2,1) * t523 - rSges(2,2) * t519;
t503 = rSges(2,1) * t519 + rSges(2,2) * t523;
t493 = -t514 * t552 + t515 * t521;
t492 = -t515 * t517 - t521 * t554;
t485 = -t513 * t554 + t515 * t540;
t484 = t515 * t513 + t522 * t536;
t482 = rSges(3,3) * t515 + (rSges(3,1) * t518 + rSges(3,2) * t522) * t514;
t469 = -t521 * t553 + t558;
t468 = t494 * t521 + t517 * t553;
t467 = t521 * t555 + t557;
t466 = t496 * t521 - t517 * t555;
t465 = t494 * t513 - t523 * t536;
t464 = t494 * t540 + t513 * t553;
t463 = t496 * t513 + t519 * t536;
t462 = -t496 * t540 + t513 * t555;
t459 = t485 * t520 + t516 * t556;
t458 = -t485 * t516 + t520 * t556;
t453 = rSges(3,1) * t497 - rSges(3,2) * t496 + rSges(3,3) * t555;
t452 = rSges(3,1) * t495 - rSges(3,2) * t494 - rSges(3,3) * t553;
t451 = -rSges(4,1) * t553 - rSges(4,2) * t495 + rSges(4,3) * t494;
t450 = rSges(4,1) * t555 - rSges(4,2) * t497 + rSges(4,3) * t496;
t434 = pkin(5) * t485 + pkin(11) * t484;
t433 = rSges(5,1) * t493 + rSges(5,2) * t492 + rSges(5,3) * t556;
t430 = Icges(5,1) * t493 + Icges(5,4) * t492 + Icges(5,5) * t556;
t429 = Icges(5,4) * t493 + Icges(5,2) * t492 + Icges(5,6) * t556;
t428 = Icges(5,5) * t493 + Icges(5,6) * t492 + Icges(5,3) * t556;
t427 = qJD(6) * t484 + t474;
t426 = rSges(6,1) * t485 - rSges(6,2) * t484 + rSges(6,3) * t556;
t425 = Icges(6,1) * t485 - Icges(6,4) * t484 + Icges(6,5) * t556;
t424 = Icges(6,4) * t485 - Icges(6,2) * t484 + Icges(6,6) * t556;
t423 = Icges(6,5) * t485 - Icges(6,6) * t484 + Icges(6,3) * t556;
t422 = t465 * t520 + t495 * t516;
t421 = -t465 * t516 + t495 * t520;
t420 = t463 * t520 + t497 * t516;
t419 = -t463 * t516 + t497 * t520;
t417 = pkin(5) * t465 - pkin(11) * t464;
t416 = pkin(5) * t463 + pkin(11) * t462;
t413 = -qJD(6) * t464 + t432;
t412 = qJD(6) * t462 + t431;
t411 = rSges(5,1) * t469 + rSges(5,2) * t468 + rSges(5,3) * t495;
t410 = rSges(5,1) * t467 + rSges(5,2) * t466 + rSges(5,3) * t497;
t409 = Icges(5,1) * t469 + Icges(5,4) * t468 + Icges(5,5) * t495;
t408 = Icges(5,1) * t467 + Icges(5,4) * t466 + Icges(5,5) * t497;
t407 = Icges(5,4) * t469 + Icges(5,2) * t468 + Icges(5,6) * t495;
t406 = Icges(5,4) * t467 + Icges(5,2) * t466 + Icges(5,6) * t497;
t405 = Icges(5,5) * t469 + Icges(5,6) * t468 + Icges(5,3) * t495;
t404 = Icges(5,5) * t467 + Icges(5,6) * t466 + Icges(5,3) * t497;
t403 = rSges(6,1) * t465 + rSges(6,2) * t464 + rSges(6,3) * t495;
t402 = rSges(6,1) * t463 - rSges(6,2) * t462 + rSges(6,3) * t497;
t401 = Icges(6,1) * t465 + Icges(6,4) * t464 + Icges(6,5) * t495;
t400 = Icges(6,1) * t463 - Icges(6,4) * t462 + Icges(6,5) * t497;
t399 = Icges(6,4) * t465 + Icges(6,2) * t464 + Icges(6,6) * t495;
t398 = Icges(6,4) * t463 - Icges(6,2) * t462 + Icges(6,6) * t497;
t397 = Icges(6,5) * t465 + Icges(6,6) * t464 + Icges(6,3) * t495;
t396 = Icges(6,5) * t463 - Icges(6,6) * t462 + Icges(6,3) * t497;
t394 = rSges(7,1) * t459 + rSges(7,2) * t458 + rSges(7,3) * t484;
t393 = Icges(7,1) * t459 + Icges(7,4) * t458 + Icges(7,5) * t484;
t392 = Icges(7,4) * t459 + Icges(7,2) * t458 + Icges(7,6) * t484;
t391 = Icges(7,5) * t459 + Icges(7,6) * t458 + Icges(7,3) * t484;
t390 = t453 * t510 - t482 * t509 + t500;
t389 = -t452 * t510 - t482 * t541 - t545;
t388 = (t452 * t519 + t453 * t523) * t544;
t386 = rSges(7,1) * t422 + rSges(7,2) * t421 - rSges(7,3) * t464;
t385 = rSges(7,1) * t420 + rSges(7,2) * t419 + rSges(7,3) * t462;
t384 = Icges(7,1) * t422 + Icges(7,4) * t421 - Icges(7,5) * t464;
t383 = Icges(7,1) * t420 + Icges(7,4) * t419 + Icges(7,5) * t462;
t382 = Icges(7,4) * t422 + Icges(7,2) * t421 - Icges(7,6) * t464;
t381 = Icges(7,4) * t420 + Icges(7,2) * t419 + Icges(7,6) * t462;
t380 = Icges(7,5) * t422 + Icges(7,6) * t421 - Icges(7,3) * t464;
t379 = Icges(7,5) * t420 + Icges(7,6) * t419 + Icges(7,3) * t462;
t378 = t450 * t510 + t519 * t535 + t542;
t377 = (-t451 - t455) * t510 + t523 * t535 + t539;
t376 = (-t543 + (t450 * t523 + t451 * t519) * qJD(2)) * t514 + t546;
t375 = t410 * t499 - t433 * t470 + t530;
t374 = -t411 * t499 + t433 * t471 + t527;
t373 = -t410 * t471 + t411 * t470 + t531;
t372 = t402 * t474 - t426 * t431 + t528;
t371 = -t403 * t474 + t426 * t432 + t526;
t370 = -t402 * t432 + t403 * t431 + t529;
t369 = t385 * t427 - t394 * t412 + t416 * t474 - t431 * t434 + t528;
t368 = -t386 * t427 + t394 * t413 - t417 * t474 + t432 * t434 + t526;
t367 = -t385 * t413 + t386 * t412 - t416 * t432 + t417 * t431 + t529;
t1 = m(4) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + t412 * ((t462 * t379 + t419 * t381 + t420 * t383) * t412 + (t380 * t462 + t382 * t419 + t384 * t420) * t413 + (t391 * t462 + t392 * t419 + t393 * t420) * t427) / 0.2e1 + t427 * ((t379 * t484 + t381 * t458 + t383 * t459) * t412 + (t380 * t484 + t382 * t458 + t384 * t459) * t413 + (t484 * t391 + t458 * t392 + t459 * t393) * t427) / 0.2e1 + t413 * ((-t379 * t464 + t381 * t421 + t383 * t422) * t412 + (-t464 * t380 + t421 * t382 + t422 * t384) * t413 + (-t391 * t464 + t392 * t421 + t393 * t422) * t427) / 0.2e1 + t432 * ((t396 * t495 + t398 * t464 + t400 * t465) * t431 + (t495 * t397 + t464 * t399 + t465 * t401) * t432 + (t423 * t495 + t424 * t464 + t425 * t465) * t474) / 0.2e1 + t474 * ((t396 * t556 - t398 * t484 + t400 * t485) * t431 + (t397 * t556 - t399 * t484 + t401 * t485) * t432 + (t423 * t556 - t424 * t484 + t425 * t485) * t474) / 0.2e1 + t431 * ((t497 * t396 - t462 * t398 + t463 * t400) * t431 + (t397 * t497 - t399 * t462 + t401 * t463) * t432 + (t423 * t497 - t424 * t462 + t425 * t463) * t474) / 0.2e1 + t470 * ((t497 * t404 + t466 * t406 + t467 * t408) * t470 + (t405 * t497 + t407 * t466 + t409 * t467) * t471 + (t428 * t497 + t429 * t466 + t430 * t467) * t499) / 0.2e1 + t471 * ((t404 * t495 + t406 * t468 + t408 * t469) * t470 + (t495 * t405 + t468 * t407 + t469 * t409) * t471 + (t428 * t495 + t429 * t468 + t430 * t469) * t499) / 0.2e1 + t499 * ((t404 * t556 + t406 * t492 + t408 * t493) * t470 + (t405 * t556 + t407 * t492 + t409 * t493) * t471 + (t428 * t556 + t429 * t492 + t430 * t493) * t499) / 0.2e1 + m(7) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(6) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(5) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(3) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + ((t563 * t515 + ((t518 * t567 + t522 * t569) * t523 + (t568 * t518 - t522 * t570) * t519) * t514) * t544 + (t564 * t515 + (t518 * t565 + t522 * t566) * t514) * t510) * t510 / 0.2e1 + (m(2) * (t503 ^ 2 + t504 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 - (((-t494 * t569 + t495 * t567 - t571) * t523 + (t494 * t570 + t568 * t495) * t519) * t544 + (-t494 * t566 + t495 * t565 - t553 * t564) * t510) * t541 / 0.2e1 + (((-t496 * t569 + t497 * t567) * t523 + (t496 * t570 + t568 * t497 + t571) * t519) * t544 + (-t496 * t566 + t497 * t565 + t555 * t564) * t510) * t509 / 0.2e1;
T  = t1;
