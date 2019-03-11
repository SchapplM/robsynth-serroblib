% Calculate kinetic energy for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:53
% EndTime: 2019-03-08 20:40:55
% DurationCPUTime: 3.03s
% Computational Cost: add. (2304->331), mult. (4585->517), div. (0->0), fcn. (5392->12), ass. (0->156)
t567 = Icges(3,1) + Icges(4,2);
t566 = Icges(3,4) + Icges(4,6);
t565 = Icges(3,5) - Icges(4,4);
t564 = Icges(3,2) + Icges(4,3);
t563 = Icges(3,6) - Icges(4,5);
t562 = Icges(3,3) + Icges(4,1);
t511 = cos(pkin(6));
t509 = sin(pkin(6));
t510 = cos(pkin(11));
t542 = t510 * t509;
t508 = sin(pkin(11));
t545 = t508 * t509;
t514 = sin(qJ(2));
t517 = cos(qJ(2));
t554 = t562 * t511 + (t565 * t514 + t563 * t517) * t509;
t540 = t511 * t517;
t491 = t508 * t514 - t510 * t540;
t541 = t511 * t514;
t492 = t508 * t517 + t510 * t541;
t555 = -t563 * t491 + t565 * t492 - t562 * t542;
t493 = t508 * t540 + t510 * t514;
t494 = -t508 * t541 + t510 * t517;
t556 = -t563 * t493 + t565 * t494 + t562 * t545;
t561 = -t554 * t511 + t555 * t542 - t556 * t545;
t560 = t564 * t493 - t566 * t494 - t563 * t545;
t559 = t564 * t491 - t566 * t492 + t563 * t542;
t558 = -t566 * t493 + t567 * t494 + t565 * t545;
t557 = t566 * t491 - t567 * t492 + t565 * t542;
t553 = t563 * t511 + (t566 * t514 + t564 * t517) * t509;
t552 = t565 * t511 + (t567 * t514 + t566 * t517) * t509;
t550 = qJD(2) ^ 2;
t516 = cos(qJ(4));
t549 = pkin(4) * t516;
t513 = sin(qJ(4));
t547 = t491 * t513;
t546 = t493 * t513;
t544 = t509 * t514;
t543 = t509 * t517;
t539 = t513 * t517;
t538 = qJ(4) + qJ(5);
t455 = pkin(2) * t494 + qJ(3) * t493;
t506 = qJD(2) * t511;
t537 = qJD(3) * t491 + t455 * t506;
t536 = qJD(2) * t509;
t502 = t508 * t536;
t469 = qJD(4) * t494 + t502;
t498 = qJD(4) * t544 + t506;
t535 = qJD(3) * t517;
t426 = qJD(5) * t494 + t469;
t474 = qJD(5) * t544 + t498;
t533 = t510 * t536;
t454 = pkin(2) * t492 + qJ(3) * t491;
t532 = t454 * t502 + t455 * t533 + qJD(1);
t531 = cos(t538);
t497 = (pkin(2) * t514 - qJ(3) * t517) * t509;
t529 = (-rSges(4,1) * t511 - (-rSges(4,2) * t514 - rSges(4,3) * t517) * t509 - t497) * t509;
t528 = (-t511 * pkin(3) - pkin(8) * t544 - t497) * t509;
t525 = t509 * t531;
t470 = qJD(4) * t492 - t533;
t427 = qJD(5) * t492 + t470;
t471 = pkin(3) * t545 + t494 * pkin(8);
t472 = -pkin(3) * t542 + t492 * pkin(8);
t524 = t471 * t533 + t472 * t502 - t509 * t535 + t532;
t523 = qJD(2) * t508 * t528 + t471 * t506 + t537;
t411 = pkin(4) * t546 + pkin(9) * t494 + t545 * t549;
t412 = pkin(4) * t547 + pkin(9) * t492 - t542 * t549;
t522 = -t411 * t470 + t469 * t412 + t524;
t490 = qJD(3) * t493;
t521 = t490 + ((-t454 - t472) * t511 + t510 * t528) * qJD(2);
t453 = t549 * t511 + (-pkin(4) * t539 + pkin(9) * t514) * t509;
t520 = t498 * t411 - t453 * t469 + t523;
t519 = -t412 * t498 + t470 * t453 + t521;
t515 = cos(qJ(6));
t512 = sin(qJ(6));
t507 = sin(t538);
t496 = -t509 * t539 + t511 * t516;
t495 = -t511 * t513 - t516 * t543;
t484 = -t507 * t543 + t511 * t531;
t483 = t511 * t507 + t517 * t525;
t481 = rSges(3,3) * t511 + (rSges(3,1) * t514 + rSges(3,2) * t517) * t509;
t468 = -t516 * t542 + t547;
t467 = t491 * t516 + t513 * t542;
t466 = t516 * t545 + t546;
t465 = t493 * t516 - t513 * t545;
t463 = t484 * t515 + t512 * t544;
t462 = -t484 * t512 + t515 * t544;
t461 = t491 * t507 - t510 * t525;
t460 = t491 * t531 + t507 * t542;
t459 = t493 * t507 + t508 * t525;
t458 = -t493 * t531 + t507 * t545;
t451 = qJD(6) * t483 + t474;
t450 = pkin(5) * t484 + pkin(10) * t483;
t449 = rSges(5,1) * t496 + rSges(5,2) * t495 + rSges(5,3) * t544;
t448 = Icges(5,1) * t496 + Icges(5,4) * t495 + Icges(5,5) * t544;
t447 = Icges(5,4) * t496 + Icges(5,2) * t495 + Icges(5,6) * t544;
t446 = Icges(5,5) * t496 + Icges(5,6) * t495 + Icges(5,3) * t544;
t445 = rSges(3,1) * t494 - rSges(3,2) * t493 + rSges(3,3) * t545;
t444 = rSges(3,1) * t492 - rSges(3,2) * t491 - rSges(3,3) * t542;
t443 = -rSges(4,1) * t542 - rSges(4,2) * t492 + rSges(4,3) * t491;
t442 = rSges(4,1) * t545 - rSges(4,2) * t494 + rSges(4,3) * t493;
t425 = rSges(6,1) * t484 - rSges(6,2) * t483 + rSges(6,3) * t544;
t424 = Icges(6,1) * t484 - Icges(6,4) * t483 + Icges(6,5) * t544;
t423 = Icges(6,4) * t484 - Icges(6,2) * t483 + Icges(6,6) * t544;
t422 = Icges(6,5) * t484 - Icges(6,6) * t483 + Icges(6,3) * t544;
t421 = t461 * t515 + t492 * t512;
t420 = -t461 * t512 + t492 * t515;
t419 = t459 * t515 + t494 * t512;
t418 = -t459 * t512 + t494 * t515;
t416 = pkin(5) * t461 - pkin(10) * t460;
t415 = pkin(5) * t459 + pkin(10) * t458;
t414 = (-t444 * t511 - t481 * t542) * qJD(2);
t413 = (t445 * t511 - t481 * t545) * qJD(2);
t410 = -qJD(6) * t460 + t427;
t409 = qJD(6) * t458 + t426;
t408 = rSges(5,1) * t468 + rSges(5,2) * t467 + rSges(5,3) * t492;
t407 = rSges(5,1) * t466 + rSges(5,2) * t465 + rSges(5,3) * t494;
t406 = Icges(5,1) * t468 + Icges(5,4) * t467 + Icges(5,5) * t492;
t405 = Icges(5,1) * t466 + Icges(5,4) * t465 + Icges(5,5) * t494;
t404 = Icges(5,4) * t468 + Icges(5,2) * t467 + Icges(5,6) * t492;
t403 = Icges(5,4) * t466 + Icges(5,2) * t465 + Icges(5,6) * t494;
t402 = Icges(5,5) * t468 + Icges(5,6) * t467 + Icges(5,3) * t492;
t401 = Icges(5,5) * t466 + Icges(5,6) * t465 + Icges(5,3) * t494;
t400 = rSges(6,1) * t461 + rSges(6,2) * t460 + rSges(6,3) * t492;
t399 = rSges(6,1) * t459 - rSges(6,2) * t458 + rSges(6,3) * t494;
t398 = Icges(6,1) * t461 + Icges(6,4) * t460 + Icges(6,5) * t492;
t397 = Icges(6,1) * t459 - Icges(6,4) * t458 + Icges(6,5) * t494;
t396 = Icges(6,4) * t461 + Icges(6,2) * t460 + Icges(6,6) * t492;
t395 = Icges(6,4) * t459 - Icges(6,2) * t458 + Icges(6,6) * t494;
t394 = Icges(6,5) * t461 + Icges(6,6) * t460 + Icges(6,3) * t492;
t393 = Icges(6,5) * t459 - Icges(6,6) * t458 + Icges(6,3) * t494;
t392 = rSges(7,1) * t463 + rSges(7,2) * t462 + rSges(7,3) * t483;
t391 = Icges(7,1) * t463 + Icges(7,4) * t462 + Icges(7,5) * t483;
t390 = Icges(7,4) * t463 + Icges(7,2) * t462 + Icges(7,6) * t483;
t389 = Icges(7,5) * t463 + Icges(7,6) * t462 + Icges(7,3) * t483;
t386 = qJD(1) + (t444 * t508 + t445 * t510) * t536;
t385 = rSges(7,1) * t421 + rSges(7,2) * t420 - rSges(7,3) * t460;
t384 = rSges(7,1) * t419 + rSges(7,2) * t418 + rSges(7,3) * t458;
t383 = Icges(7,1) * t421 + Icges(7,4) * t420 - Icges(7,5) * t460;
t382 = Icges(7,1) * t419 + Icges(7,4) * t418 + Icges(7,5) * t458;
t381 = Icges(7,4) * t421 + Icges(7,2) * t420 - Icges(7,6) * t460;
t380 = Icges(7,4) * t419 + Icges(7,2) * t418 + Icges(7,6) * t458;
t379 = Icges(7,5) * t421 + Icges(7,6) * t420 - Icges(7,3) * t460;
t378 = Icges(7,5) * t419 + Icges(7,6) * t418 + Icges(7,3) * t458;
t377 = t490 + ((-t443 - t454) * t511 + t510 * t529) * qJD(2);
t376 = (t442 * t511 + t508 * t529) * qJD(2) + t537;
t375 = (-t535 + (t442 * t510 + t443 * t508) * qJD(2)) * t509 + t532;
t374 = -t408 * t498 + t449 * t470 + t521;
t373 = t407 * t498 - t449 * t469 + t523;
t372 = -t407 * t470 + t408 * t469 + t524;
t371 = -t400 * t474 + t425 * t427 + t519;
t370 = t399 * t474 - t425 * t426 + t520;
t369 = -t399 * t427 + t400 * t426 + t522;
t368 = -t385 * t451 + t392 * t410 - t416 * t474 + t427 * t450 + t519;
t367 = t384 * t451 - t392 * t409 + t415 * t474 - t426 * t450 + t520;
t366 = -t384 * t410 + t385 * t409 - t415 * t427 + t416 * t426 + t522;
t1 = m(3) * (t386 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t451 * ((t378 * t483 + t380 * t462 + t382 * t463) * t409 + (t379 * t483 + t381 * t462 + t383 * t463) * t410 + (t483 * t389 + t462 * t390 + t463 * t391) * t451) / 0.2e1 + t409 * ((t458 * t378 + t418 * t380 + t419 * t382) * t409 + (t379 * t458 + t381 * t418 + t383 * t419) * t410 + (t389 * t458 + t390 * t418 + t391 * t419) * t451) / 0.2e1 + t410 * ((-t378 * t460 + t380 * t420 + t382 * t421) * t409 + (-t460 * t379 + t420 * t381 + t421 * t383) * t410 + (-t389 * t460 + t390 * t420 + t391 * t421) * t451) / 0.2e1 + t426 * ((t494 * t393 - t458 * t395 + t459 * t397) * t426 + (t394 * t494 - t396 * t458 + t398 * t459) * t427 + (t422 * t494 - t423 * t458 + t424 * t459) * t474) / 0.2e1 + t427 * ((t393 * t492 + t395 * t460 + t397 * t461) * t426 + (t492 * t394 + t460 * t396 + t461 * t398) * t427 + (t422 * t492 + t423 * t460 + t424 * t461) * t474) / 0.2e1 + t474 * ((t393 * t544 - t395 * t483 + t397 * t484) * t426 + (t394 * t544 - t396 * t483 + t398 * t484) * t427 + (t422 * t544 - t483 * t423 + t484 * t424) * t474) / 0.2e1 + t498 * ((t401 * t544 + t403 * t495 + t405 * t496) * t469 + (t402 * t544 + t404 * t495 + t406 * t496) * t470 + (t446 * t544 + t495 * t447 + t496 * t448) * t498) / 0.2e1 + t470 * ((t401 * t492 + t403 * t467 + t405 * t468) * t469 + (t492 * t402 + t467 * t404 + t468 * t406) * t470 + (t446 * t492 + t447 * t467 + t448 * t468) * t498) / 0.2e1 + t469 * ((t494 * t401 + t465 * t403 + t466 * t405) * t469 + (t402 * t494 + t404 * t465 + t406 * t466) * t470 + (t446 * t494 + t447 * t465 + t448 * t466) * t498) / 0.2e1 + m(7) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(6) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(5) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(4) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 - ((t491 * t560 + t558 * t492) * t545 + (-t491 * t553 + t492 * t552) * t511 + (-t491 * t559 + t492 * t557 + t561) * t542) * t550 * t542 / 0.2e1 + ((t554 * t511 ^ 2 + (((t514 * t557 + t517 * t559) * t510 + (t558 * t514 - t517 * t560) * t508) * t509 + (t508 * t556 - t510 * t555 + t514 * t552 + t517 * t553) * t511) * t509) * t511 + ((-t493 * t559 + t494 * t557) * t542 + (-t493 * t553 + t494 * t552) * t511 + (t493 * t560 + t558 * t494 - t561) * t545) * t545) * t550 / 0.2e1;
T  = t1;
