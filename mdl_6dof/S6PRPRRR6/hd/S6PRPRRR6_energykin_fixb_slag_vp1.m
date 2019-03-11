% Calculate kinetic energy for
% S6PRPRRR6
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:10
% EndTime: 2019-03-08 20:45:13
% DurationCPUTime: 2.98s
% Computational Cost: add. (2226->331), mult. (5237->515), div. (0->0), fcn. (6308->12), ass. (0->156)
t563 = Icges(3,1) + Icges(4,2);
t562 = Icges(3,4) + Icges(4,6);
t561 = Icges(3,5) - Icges(4,4);
t560 = Icges(3,2) + Icges(4,3);
t559 = Icges(3,6) - Icges(4,5);
t558 = Icges(3,3) + Icges(4,1);
t509 = cos(pkin(6));
t507 = sin(pkin(6));
t508 = cos(pkin(11));
t539 = t507 * t508;
t506 = sin(pkin(11));
t540 = t506 * t507;
t512 = sin(qJ(2));
t514 = cos(qJ(2));
t550 = t558 * t509 + (t561 * t512 + t559 * t514) * t507;
t535 = t509 * t514;
t488 = t506 * t512 - t508 * t535;
t536 = t509 * t512;
t489 = t506 * t514 + t508 * t536;
t551 = -t559 * t488 + t561 * t489 - t558 * t539;
t490 = t506 * t535 + t508 * t512;
t491 = -t506 * t536 + t508 * t514;
t552 = -t559 * t490 + t561 * t491 + t558 * t540;
t557 = -t550 * t509 + t551 * t539 - t552 * t540;
t556 = t560 * t490 - t562 * t491 - t559 * t540;
t555 = t560 * t488 - t562 * t489 + t559 * t539;
t554 = -t562 * t490 + t563 * t491 + t561 * t540;
t553 = t562 * t488 - t563 * t489 + t561 * t539;
t549 = t559 * t509 + (t562 * t512 + t560 * t514) * t507;
t548 = t561 * t509 + (t563 * t512 + t562 * t514) * t507;
t546 = qJD(2) ^ 2;
t545 = cos(qJ(4));
t513 = cos(qJ(5));
t544 = pkin(5) * t513;
t510 = sin(qJ(5));
t542 = t489 * t510;
t541 = t491 * t510;
t511 = sin(qJ(4));
t538 = t507 * t511;
t537 = t507 * t512;
t454 = pkin(2) * t491 + qJ(3) * t490;
t502 = qJD(2) * t509;
t534 = qJD(3) * t488 + t454 * t502;
t533 = qJD(2) * t507;
t499 = t506 * t533;
t470 = qJD(4) * t491 + t499;
t495 = qJD(4) * t537 + t502;
t532 = qJD(3) * t514;
t530 = t510 * t537;
t464 = -t490 * t545 + t506 * t538;
t419 = qJD(5) * t464 + t470;
t529 = t507 * t545;
t492 = t509 * t511 + t514 * t529;
t462 = qJD(5) * t492 + t495;
t528 = t508 * t533;
t453 = pkin(2) * t489 + qJ(3) * t488;
t527 = t453 * t499 + t454 * t528 + qJD(1);
t494 = (pkin(2) * t512 - qJ(3) * t514) * t507;
t525 = (-rSges(4,1) * t509 - (-rSges(4,2) * t512 - rSges(4,3) * t514) * t507 - t494) * t507;
t524 = (-pkin(3) * t509 - pkin(8) * t537 - t494) * t507;
t471 = qJD(4) * t489 - t528;
t466 = t488 * t545 + t508 * t538;
t420 = -qJD(5) * t466 + t471;
t472 = pkin(3) * t540 + pkin(8) * t491;
t473 = -pkin(3) * t539 + pkin(8) * t489;
t521 = t472 * t528 + t473 * t499 - t507 * t532 + t527;
t520 = qJD(2) * t506 * t524 + t472 * t502 + t534;
t465 = t490 * t511 + t506 * t529;
t417 = t465 * pkin(4) + t464 * pkin(9);
t467 = t488 * t511 - t508 * t529;
t418 = t467 * pkin(4) - t466 * pkin(9);
t519 = -t417 * t471 + t470 * t418 + t521;
t486 = qJD(3) * t490;
t518 = t486 + ((-t453 - t473) * t509 + t508 * t524) * qJD(2);
t493 = t509 * t545 - t514 * t538;
t455 = t493 * pkin(4) + t492 * pkin(9);
t517 = t495 * t417 - t455 * t470 + t520;
t516 = -t418 * t495 + t471 * t455 + t518;
t505 = qJ(5) + qJ(6);
t504 = cos(t505);
t503 = sin(t505);
t480 = rSges(3,3) * t509 + (rSges(3,1) * t512 + rSges(3,2) * t514) * t507;
t469 = t493 * t513 + t530;
t468 = -t493 * t510 + t513 * t537;
t459 = t493 * t504 + t503 * t537;
t458 = -t493 * t503 + t504 * t537;
t451 = rSges(5,1) * t493 - rSges(5,2) * t492 + rSges(5,3) * t537;
t450 = Icges(5,1) * t493 - Icges(5,4) * t492 + Icges(5,5) * t537;
t449 = Icges(5,4) * t493 - Icges(5,2) * t492 + Icges(5,6) * t537;
t448 = Icges(5,5) * t493 - Icges(5,6) * t492 + Icges(5,3) * t537;
t447 = rSges(3,1) * t491 - rSges(3,2) * t490 + rSges(3,3) * t540;
t446 = rSges(3,1) * t489 - rSges(3,2) * t488 - rSges(3,3) * t539;
t445 = -rSges(4,1) * t539 - rSges(4,2) * t489 + rSges(4,3) * t488;
t444 = rSges(4,1) * t540 - rSges(4,2) * t491 + rSges(4,3) * t490;
t429 = qJD(6) * t492 + t462;
t428 = t467 * t513 + t542;
t427 = -t467 * t510 + t489 * t513;
t426 = t465 * t513 + t541;
t425 = -t465 * t510 + t491 * t513;
t424 = t467 * t504 + t489 * t503;
t423 = -t467 * t503 + t489 * t504;
t422 = t465 * t504 + t491 * t503;
t421 = -t465 * t503 + t491 * t504;
t414 = (-t446 * t509 - t480 * t539) * qJD(2);
t413 = (t447 * t509 - t480 * t540) * qJD(2);
t412 = rSges(6,1) * t469 + rSges(6,2) * t468 + rSges(6,3) * t492;
t411 = Icges(6,1) * t469 + Icges(6,4) * t468 + Icges(6,5) * t492;
t410 = Icges(6,4) * t469 + Icges(6,2) * t468 + Icges(6,6) * t492;
t409 = Icges(6,5) * t469 + Icges(6,6) * t468 + Icges(6,3) * t492;
t408 = rSges(5,1) * t467 + rSges(5,2) * t466 + rSges(5,3) * t489;
t407 = rSges(5,1) * t465 - rSges(5,2) * t464 + rSges(5,3) * t491;
t406 = Icges(5,1) * t467 + Icges(5,4) * t466 + Icges(5,5) * t489;
t405 = Icges(5,1) * t465 - Icges(5,4) * t464 + Icges(5,5) * t491;
t404 = Icges(5,4) * t467 + Icges(5,2) * t466 + Icges(5,6) * t489;
t403 = Icges(5,4) * t465 - Icges(5,2) * t464 + Icges(5,6) * t491;
t402 = Icges(5,5) * t467 + Icges(5,6) * t466 + Icges(5,3) * t489;
t401 = Icges(5,5) * t465 - Icges(5,6) * t464 + Icges(5,3) * t491;
t400 = pkin(5) * t530 + pkin(10) * t492 + t493 * t544;
t399 = rSges(7,1) * t459 + rSges(7,2) * t458 + rSges(7,3) * t492;
t398 = Icges(7,1) * t459 + Icges(7,4) * t458 + Icges(7,5) * t492;
t397 = Icges(7,4) * t459 + Icges(7,2) * t458 + Icges(7,6) * t492;
t396 = Icges(7,5) * t459 + Icges(7,6) * t458 + Icges(7,3) * t492;
t395 = -qJD(6) * t466 + t420;
t394 = qJD(6) * t464 + t419;
t392 = qJD(1) + (t446 * t506 + t447 * t508) * t533;
t391 = rSges(6,1) * t428 + rSges(6,2) * t427 - rSges(6,3) * t466;
t390 = rSges(6,1) * t426 + rSges(6,2) * t425 + rSges(6,3) * t464;
t389 = Icges(6,1) * t428 + Icges(6,4) * t427 - Icges(6,5) * t466;
t388 = Icges(6,1) * t426 + Icges(6,4) * t425 + Icges(6,5) * t464;
t387 = Icges(6,4) * t428 + Icges(6,2) * t427 - Icges(6,6) * t466;
t386 = Icges(6,4) * t426 + Icges(6,2) * t425 + Icges(6,6) * t464;
t385 = Icges(6,5) * t428 + Icges(6,6) * t427 - Icges(6,3) * t466;
t384 = Icges(6,5) * t426 + Icges(6,6) * t425 + Icges(6,3) * t464;
t383 = rSges(7,1) * t424 + rSges(7,2) * t423 - rSges(7,3) * t466;
t382 = rSges(7,1) * t422 + rSges(7,2) * t421 + rSges(7,3) * t464;
t381 = Icges(7,1) * t424 + Icges(7,4) * t423 - Icges(7,5) * t466;
t380 = Icges(7,1) * t422 + Icges(7,4) * t421 + Icges(7,5) * t464;
t379 = Icges(7,4) * t424 + Icges(7,2) * t423 - Icges(7,6) * t466;
t378 = Icges(7,4) * t422 + Icges(7,2) * t421 + Icges(7,6) * t464;
t377 = Icges(7,5) * t424 + Icges(7,6) * t423 - Icges(7,3) * t466;
t376 = Icges(7,5) * t422 + Icges(7,6) * t421 + Icges(7,3) * t464;
t375 = pkin(5) * t542 - pkin(10) * t466 + t467 * t544;
t374 = pkin(5) * t541 + pkin(10) * t464 + t465 * t544;
t373 = t486 + ((-t445 - t453) * t509 + t508 * t525) * qJD(2);
t372 = (t444 * t509 + t506 * t525) * qJD(2) + t534;
t371 = (-t532 + (t444 * t508 + t445 * t506) * qJD(2)) * t507 + t527;
t370 = -t408 * t495 + t451 * t471 + t518;
t369 = t407 * t495 - t451 * t470 + t520;
t368 = -t407 * t471 + t408 * t470 + t521;
t367 = -t391 * t462 + t412 * t420 + t516;
t366 = t390 * t462 - t412 * t419 + t517;
t365 = -t390 * t420 + t391 * t419 + t519;
t364 = -t375 * t462 - t383 * t429 + t395 * t399 + t400 * t420 + t516;
t363 = t374 * t462 + t382 * t429 - t394 * t399 - t400 * t419 + t517;
t362 = -t374 * t420 + t375 * t419 - t382 * t395 + t383 * t394 + t519;
t1 = t462 * ((t384 * t492 + t386 * t468 + t388 * t469) * t419 + (t385 * t492 + t387 * t468 + t389 * t469) * t420 + (t409 * t492 + t410 * t468 + t411 * t469) * t462) / 0.2e1 + t471 * ((t401 * t489 + t403 * t466 + t405 * t467) * t470 + (t402 * t489 + t404 * t466 + t406 * t467) * t471 + (t448 * t489 + t449 * t466 + t450 * t467) * t495) / 0.2e1 + t495 * ((t401 * t537 - t403 * t492 + t405 * t493) * t470 + (t402 * t537 - t404 * t492 + t406 * t493) * t471 + (t448 * t537 - t449 * t492 + t450 * t493) * t495) / 0.2e1 + t470 * ((t401 * t491 - t403 * t464 + t405 * t465) * t470 + (t402 * t491 - t404 * t464 + t406 * t465) * t471 + (t448 * t491 - t449 * t464 + t450 * t465) * t495) / 0.2e1 + m(7) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(6) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(5) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(3) * (t392 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t429 * ((t376 * t492 + t378 * t458 + t380 * t459) * t394 + (t377 * t492 + t379 * t458 + t381 * t459) * t395 + (t396 * t492 + t397 * t458 + t398 * t459) * t429) / 0.2e1 + t394 * ((t464 * t376 + t421 * t378 + t422 * t380) * t394 + (t377 * t464 + t379 * t421 + t381 * t422) * t395 + (t396 * t464 + t397 * t421 + t398 * t422) * t429) / 0.2e1 + t395 * ((-t376 * t466 + t378 * t423 + t380 * t424) * t394 + (-t466 * t377 + t423 * t379 + t424 * t381) * t395 + (-t396 * t466 + t397 * t423 + t398 * t424) * t429) / 0.2e1 + t419 * ((t384 * t464 + t386 * t425 + t388 * t426) * t419 + (t385 * t464 + t387 * t425 + t389 * t426) * t420 + (t409 * t464 + t410 * t425 + t411 * t426) * t462) / 0.2e1 + t420 * ((-t384 * t466 + t386 * t427 + t388 * t428) * t419 + (-t385 * t466 + t387 * t427 + t389 * t428) * t420 + (-t409 * t466 + t410 * t427 + t411 * t428) * t462) / 0.2e1 - ((t488 * t556 + t554 * t489) * t540 + (-t488 * t549 + t489 * t548) * t509 + (-t488 * t555 + t489 * t553 + t557) * t539) * t546 * t539 / 0.2e1 + ((t550 * t509 ^ 2 + (((t512 * t553 + t514 * t555) * t508 + (t554 * t512 - t514 * t556) * t506) * t507 + (t506 * t552 - t508 * t551 + t512 * t548 + t514 * t549) * t509) * t507) * t509 + ((-t490 * t555 + t491 * t553) * t539 + (-t490 * t549 + t491 * t548) * t509 + (t490 * t556 + t554 * t491 - t557) * t540) * t540) * t546 / 0.2e1;
T  = t1;
