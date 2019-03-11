% Calculate kinetic energy for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:01
% EndTime: 2019-03-09 09:39:04
% DurationCPUTime: 3.08s
% Computational Cost: add. (2268->339), mult. (4427->505), div. (0->0), fcn. (5174->12), ass. (0->160)
t581 = Icges(3,1) + Icges(4,2);
t580 = Icges(3,4) + Icges(4,6);
t579 = Icges(3,5) - Icges(4,4);
t578 = Icges(3,2) + Icges(4,3);
t577 = Icges(3,6) - Icges(4,5);
t576 = Icges(3,3) + Icges(4,1);
t514 = sin(pkin(6));
t516 = cos(pkin(6));
t522 = cos(qJ(2));
t523 = cos(qJ(1));
t551 = t522 * t523;
t519 = sin(qJ(2));
t520 = sin(qJ(1));
t554 = t519 * t520;
t493 = -t516 * t551 + t554;
t552 = t520 * t522;
t553 = t519 * t523;
t494 = t516 * t553 + t552;
t495 = t516 * t552 + t553;
t496 = -t516 * t554 + t551;
t555 = t514 * t523;
t557 = t514 * t520;
t567 = (t493 * t577 - t494 * t579 + t555 * t576) * t523 + (-t495 * t577 + t496 * t579 + t557 * t576) * t520;
t575 = t567 * t514;
t574 = t495 * t578 - t496 * t580 - t557 * t577;
t573 = t493 * t578 - t494 * t580 + t555 * t577;
t572 = -t580 * t495 + t496 * t581 + t579 * t557;
t571 = t580 * t493 - t494 * t581 + t579 * t555;
t570 = t577 * t516 + (t519 * t580 + t522 * t578) * t514;
t569 = t579 * t516 + (t519 * t581 + t580 * t522) * t514;
t568 = t576 * t516 + (t519 * t579 + t522 * t577) * t514;
t513 = sin(pkin(11));
t515 = cos(pkin(11));
t466 = t495 * t515 - t513 * t557;
t559 = t495 * t513;
t467 = t515 * t557 + t559;
t407 = Icges(5,5) * t467 + Icges(5,6) * t466 + Icges(5,3) * t496;
t566 = t407 + t572;
t468 = t493 * t515 + t513 * t555;
t560 = t493 * t513;
t469 = -t515 * t555 + t560;
t408 = Icges(5,5) * t469 + Icges(5,6) * t468 + Icges(5,3) * t494;
t565 = -t408 + t571;
t556 = t514 * t522;
t491 = -t513 * t516 - t515 * t556;
t492 = -t513 * t556 + t515 * t516;
t558 = t514 * t519;
t429 = Icges(5,5) * t492 + Icges(5,6) * t491 + Icges(5,3) * t558;
t564 = t429 + t569;
t561 = pkin(4) * t515;
t455 = pkin(2) * t494 + qJ(3) * t493;
t456 = pkin(2) * t496 + qJ(3) * t495;
t545 = qJD(2) * t514;
t508 = t520 * t545;
t541 = t523 * t545;
t549 = t455 * t508 + t456 * t541;
t473 = -pkin(3) * t555 + qJ(4) * t494;
t548 = -t455 - t473;
t497 = (pkin(2) * t519 - qJ(3) * t522) * t514;
t547 = -pkin(3) * t516 - qJ(4) * t558 - t497;
t470 = qJD(5) * t496 + t508;
t546 = qJD(1) * (pkin(1) * t520 - pkin(8) * t555);
t544 = qJD(3) * t522;
t509 = qJD(2) * t516 + qJD(1);
t543 = pkin(11) + qJ(5);
t499 = qJD(1) * (pkin(1) * t523 + pkin(8) * t557);
t542 = qJD(3) * t493 + t509 * t456 + t499;
t498 = qJD(5) * t558 + t509;
t540 = cos(t543);
t539 = qJD(3) * t495 - t546;
t536 = (-rSges(4,1) * t516 - (-rSges(4,2) * t519 - rSges(4,3) * t522) * t514 - t497) * t545;
t471 = qJD(5) * t494 - t541;
t535 = qJD(4) * t496 + t539;
t534 = t514 * t540;
t472 = pkin(3) * t557 + qJ(4) * t496;
t533 = qJD(4) * t494 + t509 * t472 + t542;
t532 = qJD(4) * t558 + t472 * t541 + t473 * t508 + t549;
t529 = (-rSges(5,1) * t492 - rSges(5,2) * t491 - rSges(5,3) * t558 + t547) * t545;
t528 = (-t561 * t516 - (-pkin(4) * t513 * t522 + pkin(9) * t519) * t514 + t547) * t545;
t415 = pkin(4) * t559 + pkin(9) * t496 + t557 * t561;
t416 = pkin(4) * t560 + pkin(9) * t494 - t555 * t561;
t527 = t415 * t541 + t416 * t508 - t514 * t544 + t532;
t526 = t509 * t415 + t520 * t528 + t533;
t525 = (-t416 + t548) * t509 + t523 * t528 + t535;
t521 = cos(qJ(6));
t518 = sin(qJ(6));
t512 = sin(t543);
t503 = rSges(2,1) * t523 - rSges(2,2) * t520;
t502 = rSges(2,1) * t520 + rSges(2,2) * t523;
t483 = rSges(3,3) * t516 + (rSges(3,1) * t519 + rSges(3,2) * t522) * t514;
t482 = -t512 * t556 + t516 * t540;
t481 = t516 * t512 + t522 * t534;
t463 = t493 * t512 - t523 * t534;
t462 = t493 * t540 + t512 * t555;
t461 = t495 * t512 + t520 * t534;
t460 = -t495 * t540 + t512 * t557;
t458 = t482 * t521 + t518 * t558;
t457 = -t482 * t518 + t521 * t558;
t454 = qJD(6) * t481 + t498;
t452 = rSges(3,1) * t496 - rSges(3,2) * t495 + rSges(3,3) * t557;
t451 = rSges(3,1) * t494 - rSges(3,2) * t493 - rSges(3,3) * t555;
t450 = -rSges(4,1) * t555 - rSges(4,2) * t494 + rSges(4,3) * t493;
t449 = rSges(4,1) * t557 - rSges(4,2) * t496 + rSges(4,3) * t495;
t433 = pkin(5) * t482 + pkin(10) * t481;
t431 = Icges(5,1) * t492 + Icges(5,4) * t491 + Icges(5,5) * t558;
t430 = Icges(5,4) * t492 + Icges(5,2) * t491 + Icges(5,6) * t558;
t428 = rSges(6,1) * t482 - rSges(6,2) * t481 + rSges(6,3) * t558;
t427 = Icges(6,1) * t482 - Icges(6,4) * t481 + Icges(6,5) * t558;
t426 = Icges(6,4) * t482 - Icges(6,2) * t481 + Icges(6,6) * t558;
t425 = Icges(6,5) * t482 - Icges(6,6) * t481 + Icges(6,3) * t558;
t424 = t463 * t521 + t494 * t518;
t423 = -t463 * t518 + t494 * t521;
t422 = t461 * t521 + t496 * t518;
t421 = -t461 * t518 + t496 * t521;
t420 = -qJD(6) * t462 + t471;
t419 = qJD(6) * t460 + t470;
t418 = pkin(5) * t463 - pkin(10) * t462;
t417 = pkin(5) * t461 + pkin(10) * t460;
t414 = rSges(5,1) * t469 + rSges(5,2) * t468 + rSges(5,3) * t494;
t413 = rSges(5,1) * t467 + rSges(5,2) * t466 + rSges(5,3) * t496;
t412 = Icges(5,1) * t469 + Icges(5,4) * t468 + Icges(5,5) * t494;
t411 = Icges(5,1) * t467 + Icges(5,4) * t466 + Icges(5,5) * t496;
t410 = Icges(5,4) * t469 + Icges(5,2) * t468 + Icges(5,6) * t494;
t409 = Icges(5,4) * t467 + Icges(5,2) * t466 + Icges(5,6) * t496;
t403 = rSges(6,1) * t463 + rSges(6,2) * t462 + rSges(6,3) * t494;
t402 = rSges(6,1) * t461 - rSges(6,2) * t460 + rSges(6,3) * t496;
t401 = Icges(6,1) * t463 + Icges(6,4) * t462 + Icges(6,5) * t494;
t400 = Icges(6,1) * t461 - Icges(6,4) * t460 + Icges(6,5) * t496;
t399 = Icges(6,4) * t463 + Icges(6,2) * t462 + Icges(6,6) * t494;
t398 = Icges(6,4) * t461 - Icges(6,2) * t460 + Icges(6,6) * t496;
t397 = Icges(6,5) * t463 + Icges(6,6) * t462 + Icges(6,3) * t494;
t396 = Icges(6,5) * t461 - Icges(6,6) * t460 + Icges(6,3) * t496;
t395 = t452 * t509 - t483 * t508 + t499;
t394 = -t451 * t509 - t483 * t541 - t546;
t393 = rSges(7,1) * t458 + rSges(7,2) * t457 + rSges(7,3) * t481;
t392 = Icges(7,1) * t458 + Icges(7,4) * t457 + Icges(7,5) * t481;
t391 = Icges(7,4) * t458 + Icges(7,2) * t457 + Icges(7,6) * t481;
t390 = Icges(7,5) * t458 + Icges(7,6) * t457 + Icges(7,3) * t481;
t389 = (t451 * t520 + t452 * t523) * t545;
t388 = rSges(7,1) * t424 + rSges(7,2) * t423 - rSges(7,3) * t462;
t387 = rSges(7,1) * t422 + rSges(7,2) * t421 + rSges(7,3) * t460;
t386 = Icges(7,1) * t424 + Icges(7,4) * t423 - Icges(7,5) * t462;
t385 = Icges(7,1) * t422 + Icges(7,4) * t421 + Icges(7,5) * t460;
t384 = Icges(7,4) * t424 + Icges(7,2) * t423 - Icges(7,6) * t462;
t383 = Icges(7,4) * t422 + Icges(7,2) * t421 + Icges(7,6) * t460;
t382 = Icges(7,5) * t424 + Icges(7,6) * t423 - Icges(7,3) * t462;
t381 = Icges(7,5) * t422 + Icges(7,6) * t421 + Icges(7,3) * t460;
t380 = t449 * t509 + t520 * t536 + t542;
t379 = (-t450 - t455) * t509 + t523 * t536 + t539;
t378 = (-t544 + (t449 * t523 + t450 * t520) * qJD(2)) * t514 + t549;
t377 = t413 * t509 + t520 * t529 + t533;
t376 = (-t414 + t548) * t509 + t523 * t529 + t535;
t375 = (-t544 + (t413 * t523 + t414 * t520) * qJD(2)) * t514 + t532;
t374 = t402 * t498 - t428 * t470 + t526;
t373 = -t403 * t498 + t428 * t471 + t525;
t372 = -t402 * t471 + t403 * t470 + t527;
t371 = t387 * t454 - t393 * t419 + t417 * t498 - t433 * t470 + t526;
t370 = -t388 * t454 + t393 * t420 - t418 * t498 + t433 * t471 + t525;
t369 = -t387 * t420 + t388 * t419 - t417 * t471 + t418 * t470 + t527;
t1 = m(7) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(3) * (t389 ^ 2 + t394 ^ 2 + t395 ^ 2) / 0.2e1 + m(6) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(4) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(5) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + t419 * ((t460 * t381 + t421 * t383 + t422 * t385) * t419 + (t382 * t460 + t384 * t421 + t386 * t422) * t420 + (t390 * t460 + t391 * t421 + t392 * t422) * t454) / 0.2e1 + t420 * ((-t381 * t462 + t383 * t423 + t385 * t424) * t419 + (-t462 * t382 + t423 * t384 + t424 * t386) * t420 + (-t390 * t462 + t391 * t423 + t392 * t424) * t454) / 0.2e1 + t454 * ((t381 * t481 + t383 * t457 + t385 * t458) * t419 + (t382 * t481 + t384 * t457 + t386 * t458) * t420 + (t481 * t390 + t457 * t391 + t458 * t392) * t454) / 0.2e1 + t471 * ((t396 * t494 + t398 * t462 + t400 * t463) * t470 + (t494 * t397 + t462 * t399 + t463 * t401) * t471 + (t425 * t494 + t426 * t462 + t427 * t463) * t498) / 0.2e1 + t498 * ((t396 * t558 - t398 * t481 + t400 * t482) * t470 + (t397 * t558 - t399 * t481 + t401 * t482) * t471 + (t425 * t558 - t481 * t426 + t482 * t427) * t498) / 0.2e1 + t470 * ((t496 * t396 - t460 * t398 + t461 * t400) * t470 + (t397 * t496 - t399 * t460 + t401 * t461) * t471 + (t425 * t496 - t426 * t460 + t427 * t461) * t498) / 0.2e1 + (m(2) * (t502 ^ 2 + t503 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t409 * t491 + t411 * t492) * t520 - (t410 * t491 + t412 * t492) * t523 + (t407 * t520 - t408 * t523) * t558 + t567 * t516 + ((t571 * t519 + t573 * t522) * t523 + (t572 * t519 - t574 * t522) * t520) * t514) * t545 + (t429 * t558 + t491 * t430 + t492 * t431 + t568 * t516 + (t569 * t519 + t570 * t522) * t514) * t509) * t509 / 0.2e1 + (((-t410 * t466 - t412 * t467 - t573 * t495 + t565 * t496) * t523 + (t409 * t466 + t411 * t467 + t574 * t495 + t566 * t496 + t575) * t520) * t545 + (t430 * t466 + t431 * t467 - t570 * t495 + t564 * t496 + t568 * t557) * t509) * t508 / 0.2e1 - (((-t410 * t468 - t412 * t469 - t573 * t493 + t565 * t494 - t575) * t523 + (t409 * t468 + t411 * t469 + t574 * t493 + t566 * t494) * t520) * t545 + (t430 * t468 + t431 * t469 - t570 * t493 + t564 * t494 - t568 * t555) * t509) * t541 / 0.2e1;
T  = t1;
