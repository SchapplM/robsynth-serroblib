% Calculate kinetic energy for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR13_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:07
% EndTime: 2019-03-09 11:24:10
% DurationCPUTime: 3.11s
% Computational Cost: add. (2199->335), mult. (5073->491), div. (0->0), fcn. (6066->12), ass. (0->154)
t578 = Icges(3,1) + Icges(4,2);
t577 = Icges(3,4) + Icges(4,6);
t576 = Icges(3,5) - Icges(4,4);
t575 = Icges(3,2) + Icges(4,3);
t574 = Icges(5,2) + Icges(6,3);
t573 = Icges(3,6) - Icges(4,5);
t572 = Icges(3,3) + Icges(4,1);
t512 = sin(pkin(6));
t514 = cos(pkin(6));
t519 = cos(qJ(2));
t520 = cos(qJ(1));
t544 = t519 * t520;
t517 = sin(qJ(2));
t518 = sin(qJ(1));
t547 = t517 * t518;
t490 = -t514 * t544 + t547;
t545 = t518 * t519;
t546 = t517 * t520;
t491 = t514 * t546 + t545;
t492 = t514 * t545 + t546;
t493 = -t514 * t547 + t544;
t548 = t512 * t520;
t549 = t512 * t518;
t560 = (t490 * t573 - t491 * t576 + t548 * t572) * t520 + (-t492 * t573 + t493 * t576 + t549 * t572) * t518;
t571 = t560 * t512;
t516 = sin(qJ(4));
t554 = cos(qJ(4));
t536 = t512 * t554;
t467 = t492 * t516 + t518 * t536;
t511 = sin(pkin(11));
t513 = cos(pkin(11));
t426 = -t467 * t511 + t493 * t513;
t551 = t493 * t511;
t427 = t467 * t513 + t551;
t466 = -t492 * t554 + t516 * t549;
t570 = -Icges(5,4) * t467 + Icges(6,5) * t427 - Icges(5,6) * t493 + Icges(6,6) * t426 + t466 * t574;
t469 = t490 * t516 - t520 * t536;
t428 = -t469 * t511 + t491 * t513;
t552 = t491 * t511;
t429 = t469 * t513 + t552;
t468 = t490 * t554 + t516 * t548;
t569 = -Icges(5,4) * t469 + Icges(6,5) * t429 - Icges(5,6) * t491 + Icges(6,6) * t428 - t468 * t574;
t489 = -t512 * t516 * t519 + t514 * t554;
t550 = t512 * t517;
t464 = -t489 * t511 + t513 * t550;
t538 = t511 * t550;
t465 = t489 * t513 + t538;
t488 = t514 * t516 + t519 * t536;
t568 = -Icges(5,4) * t489 + Icges(6,5) * t465 - Icges(5,6) * t550 + Icges(6,6) * t464 + t488 * t574;
t567 = t492 * t575 - t493 * t577 - t549 * t573;
t566 = t490 * t575 - t491 * t577 + t548 * t573;
t565 = -t492 * t577 + t493 * t578 + t549 * t576;
t564 = t490 * t577 - t491 * t578 + t548 * t576;
t563 = t573 * t514 + (t517 * t577 + t519 * t575) * t512;
t562 = t576 * t514 + (t517 * t578 + t519 * t577) * t512;
t561 = t572 * t514 + (t517 * t576 + t519 * t573) * t512;
t553 = pkin(5) * t513;
t454 = pkin(2) * t491 + qJ(3) * t490;
t455 = pkin(2) * t493 + qJ(3) * t492;
t540 = qJD(2) * t512;
t504 = t518 * t540;
t535 = t520 * t540;
t542 = t454 * t504 + t455 * t535;
t470 = qJD(4) * t493 + t504;
t541 = qJD(1) * (pkin(1) * t518 - pkin(8) * t548);
t539 = qJD(3) * t519;
t505 = qJD(2) * t514 + qJD(1);
t496 = qJD(1) * (pkin(1) * t520 + pkin(8) * t549);
t537 = qJD(3) * t490 + t455 * t505 + t496;
t495 = qJD(4) * t550 + t505;
t534 = qJD(3) * t492 - t541;
t494 = (pkin(2) * t517 - qJ(3) * t519) * t512;
t531 = (-rSges(4,1) * t514 - (-rSges(4,2) * t517 - rSges(4,3) * t519) * t512 - t494) * t540;
t530 = (-pkin(3) * t514 - pkin(9) * t550 - t494) * t540;
t471 = qJD(4) * t491 - t535;
t472 = pkin(3) * t549 + pkin(9) * t493;
t473 = -pkin(3) * t548 + pkin(9) * t491;
t527 = t472 * t535 + t473 * t504 - t512 * t539 + t542;
t419 = pkin(4) * t469 - qJ(5) * t468;
t526 = qJD(5) * t488 + t419 * t470 + t527;
t525 = t472 * t505 + t518 * t530 + t537;
t418 = pkin(4) * t467 + qJ(5) * t466;
t524 = -qJD(5) * t468 + t418 * t495 + t525;
t523 = (-t454 - t473) * t505 + t520 * t530 + t534;
t453 = pkin(4) * t489 + qJ(5) * t488;
t522 = qJD(5) * t466 + t453 * t471 + t523;
t510 = pkin(11) + qJ(6);
t509 = cos(t510);
t508 = sin(t510);
t500 = rSges(2,1) * t520 - rSges(2,2) * t518;
t499 = rSges(2,1) * t518 + rSges(2,2) * t520;
t480 = rSges(3,3) * t514 + (rSges(3,1) * t517 + rSges(3,2) * t519) * t512;
t458 = qJD(6) * t488 + t495;
t457 = t489 * t509 + t508 * t550;
t456 = -t489 * t508 + t509 * t550;
t452 = rSges(3,1) * t493 - rSges(3,2) * t492 + rSges(3,3) * t549;
t451 = rSges(3,1) * t491 - rSges(3,2) * t490 - rSges(3,3) * t548;
t450 = -rSges(4,1) * t548 - rSges(4,2) * t491 + rSges(4,3) * t490;
t449 = rSges(4,1) * t549 - rSges(4,2) * t493 + rSges(4,3) * t492;
t433 = rSges(5,1) * t489 - rSges(5,2) * t488 + rSges(5,3) * t550;
t432 = Icges(5,1) * t489 - Icges(5,4) * t488 + Icges(5,5) * t550;
t430 = Icges(5,5) * t489 - Icges(5,6) * t488 + Icges(5,3) * t550;
t425 = t469 * t509 + t491 * t508;
t424 = -t469 * t508 + t491 * t509;
t423 = t467 * t509 + t493 * t508;
t422 = -t467 * t508 + t493 * t509;
t421 = -qJD(6) * t468 + t471;
t420 = qJD(6) * t466 + t470;
t415 = rSges(5,1) * t469 + rSges(5,2) * t468 + rSges(5,3) * t491;
t414 = rSges(5,1) * t467 - rSges(5,2) * t466 + rSges(5,3) * t493;
t413 = Icges(5,1) * t469 + Icges(5,4) * t468 + Icges(5,5) * t491;
t412 = Icges(5,1) * t467 - Icges(5,4) * t466 + Icges(5,5) * t493;
t409 = Icges(5,5) * t469 + Icges(5,6) * t468 + Icges(5,3) * t491;
t408 = Icges(5,5) * t467 - Icges(5,6) * t466 + Icges(5,3) * t493;
t407 = rSges(6,1) * t465 + rSges(6,2) * t464 + rSges(6,3) * t488;
t406 = Icges(6,1) * t465 + Icges(6,4) * t464 + Icges(6,5) * t488;
t405 = Icges(6,4) * t465 + Icges(6,2) * t464 + Icges(6,6) * t488;
t402 = pkin(5) * t538 + pkin(10) * t488 + t489 * t553;
t401 = rSges(7,1) * t457 + rSges(7,2) * t456 + rSges(7,3) * t488;
t400 = Icges(7,1) * t457 + Icges(7,4) * t456 + Icges(7,5) * t488;
t399 = Icges(7,4) * t457 + Icges(7,2) * t456 + Icges(7,6) * t488;
t398 = Icges(7,5) * t457 + Icges(7,6) * t456 + Icges(7,3) * t488;
t397 = t452 * t505 - t480 * t504 + t496;
t396 = -t451 * t505 - t480 * t535 - t541;
t395 = (t451 * t518 + t452 * t520) * t540;
t394 = rSges(6,1) * t429 + rSges(6,2) * t428 - rSges(6,3) * t468;
t393 = rSges(6,1) * t427 + rSges(6,2) * t426 + rSges(6,3) * t466;
t392 = Icges(6,1) * t429 + Icges(6,4) * t428 - Icges(6,5) * t468;
t391 = Icges(6,1) * t427 + Icges(6,4) * t426 + Icges(6,5) * t466;
t390 = Icges(6,4) * t429 + Icges(6,2) * t428 - Icges(6,6) * t468;
t389 = Icges(6,4) * t427 + Icges(6,2) * t426 + Icges(6,6) * t466;
t386 = rSges(7,1) * t425 + rSges(7,2) * t424 - rSges(7,3) * t468;
t385 = rSges(7,1) * t423 + rSges(7,2) * t422 + rSges(7,3) * t466;
t384 = Icges(7,1) * t425 + Icges(7,4) * t424 - Icges(7,5) * t468;
t383 = Icges(7,1) * t423 + Icges(7,4) * t422 + Icges(7,5) * t466;
t382 = Icges(7,4) * t425 + Icges(7,2) * t424 - Icges(7,6) * t468;
t381 = Icges(7,4) * t423 + Icges(7,2) * t422 + Icges(7,6) * t466;
t380 = Icges(7,5) * t425 + Icges(7,6) * t424 - Icges(7,3) * t468;
t379 = Icges(7,5) * t423 + Icges(7,6) * t422 + Icges(7,3) * t466;
t378 = pkin(5) * t552 - pkin(10) * t468 + t469 * t553;
t377 = pkin(5) * t551 + pkin(10) * t466 + t467 * t553;
t376 = t449 * t505 + t518 * t531 + t537;
t375 = (-t450 - t454) * t505 + t520 * t531 + t534;
t374 = (-t539 + (t449 * t520 + t450 * t518) * qJD(2)) * t512 + t542;
t373 = t414 * t495 - t433 * t470 + t525;
t372 = -t415 * t495 + t433 * t471 + t523;
t371 = -t414 * t471 + t415 * t470 + t527;
t370 = t393 * t495 + (-t407 - t453) * t470 + t524;
t369 = t407 * t471 + (-t394 - t419) * t495 + t522;
t368 = t394 * t470 + (-t393 - t418) * t471 + t526;
t367 = t377 * t495 + t385 * t458 - t401 * t420 + (-t402 - t453) * t470 + t524;
t366 = -t386 * t458 + t401 * t421 + t402 * t471 + (-t378 - t419) * t495 + t522;
t365 = t378 * t470 - t385 * t421 + t386 * t420 + (-t377 - t418) * t471 + t526;
t1 = m(3) * (t395 ^ 2 + t396 ^ 2 + t397 ^ 2) / 0.2e1 + m(7) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(5) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(4) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(6) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + t420 * ((t466 * t379 + t422 * t381 + t423 * t383) * t420 + (t380 * t466 + t382 * t422 + t384 * t423) * t421 + (t398 * t466 + t399 * t422 + t400 * t423) * t458) / 0.2e1 + t421 * ((-t379 * t468 + t381 * t424 + t383 * t425) * t420 + (-t468 * t380 + t424 * t382 + t425 * t384) * t421 + (-t398 * t468 + t399 * t424 + t400 * t425) * t458) / 0.2e1 + t458 * ((t379 * t488 + t381 * t456 + t383 * t457) * t420 + (t380 * t488 + t382 * t456 + t384 * t457) * t421 + (t398 * t488 + t399 * t456 + t400 * t457) * t458) / 0.2e1 + ((t405 * t426 + t406 * t427 + t430 * t493 + t432 * t467 + t466 * t568) * t495 + (t390 * t426 + t392 * t427 + t409 * t493 + t413 * t467 + t466 * t569) * t471 + (t389 * t426 + t391 * t427 + t408 * t493 + t412 * t467 + t570 * t466) * t470) * t470 / 0.2e1 + ((t405 * t428 + t406 * t429 + t430 * t491 + t432 * t469 - t468 * t568) * t495 + (t390 * t428 + t392 * t429 + t409 * t491 + t413 * t469 - t569 * t468) * t471 + (t389 * t428 + t391 * t429 + t408 * t491 + t412 * t469 - t468 * t570) * t470) * t471 / 0.2e1 + ((t405 * t464 + t406 * t465 + t430 * t550 + t432 * t489 + t568 * t488) * t495 + (t390 * t464 + t392 * t465 + t409 * t550 + t413 * t489 + t488 * t569) * t471 + (t389 * t464 + t391 * t465 + t408 * t550 + t412 * t489 + t488 * t570) * t470) * t495 / 0.2e1 + ((t560 * t514 + ((t517 * t564 + t519 * t566) * t520 + (t517 * t565 - t519 * t567) * t518) * t512) * t540 + (t561 * t514 + (t517 * t562 + t519 * t563) * t512) * t505) * t505 / 0.2e1 + (Icges(2,3) + m(2) * (t499 ^ 2 + t500 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t492 * t566 + t493 * t564) * t520 + (t492 * t567 + t493 * t565 + t571) * t518) * t540 + (-t492 * t563 + t493 * t562 + t549 * t561) * t505) * t504 / 0.2e1 - (((-t490 * t566 + t491 * t564 - t571) * t520 + (t490 * t567 + t491 * t565) * t518) * t540 + (-t490 * t563 + t491 * t562 - t548 * t561) * t505) * t535 / 0.2e1;
T  = t1;
