% Calculate kinetic energy for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:37
% EndTime: 2019-03-10 03:56:39
% DurationCPUTime: 2.54s
% Computational Cost: add. (3523->389), mult. (5555->607), div. (0->0), fcn. (6564->14), ass. (0->178)
t532 = sin(pkin(6));
t533 = cos(pkin(6));
t540 = cos(qJ(2));
t541 = cos(qJ(1));
t566 = t540 * t541;
t536 = sin(qJ(2));
t537 = sin(qJ(1));
t569 = t536 * t537;
t507 = -t533 * t566 + t569;
t567 = t537 * t540;
t568 = t536 * t541;
t508 = t533 * t568 + t567;
t509 = t533 * t567 + t568;
t510 = -t533 * t569 + t566;
t571 = t532 * t541;
t574 = t532 * t537;
t553 = (Icges(3,5) * t508 - Icges(3,6) * t507 - Icges(3,3) * t571) * t541 - (Icges(3,5) * t510 - Icges(3,6) * t509 + Icges(3,3) * t574) * t537;
t580 = t532 * t553;
t539 = cos(qJ(3));
t577 = t539 * pkin(3);
t575 = t532 * t536;
t573 = t532 * t539;
t572 = t532 * t540;
t535 = sin(qJ(3));
t570 = t533 * t535;
t531 = qJ(3) + qJ(4);
t473 = pkin(2) * t508 + pkin(9) * t507;
t474 = pkin(2) * t510 + pkin(9) * t509;
t561 = qJD(2) * t532;
t522 = t537 * t561;
t556 = t541 * t561;
t565 = t473 * t522 + t474 * t556;
t528 = cos(t531);
t564 = pkin(4) * t528;
t488 = qJD(3) * t509 + t522;
t562 = qJD(1) * (pkin(1) * t537 - pkin(8) * t571);
t523 = qJD(2) * t533 + qJD(1);
t560 = -qJD(3) - qJD(4);
t559 = t535 * t574;
t558 = t535 * t571;
t557 = qJ(5) + t531;
t456 = qJD(4) * t509 + t488;
t527 = sin(t531);
t555 = pkin(4) * t527;
t437 = qJD(5) * t509 + t456;
t554 = cos(t557);
t489 = qJD(3) * t507 - t556;
t423 = -pkin(3) * t558 + pkin(10) * t507 + t508 * t577;
t424 = pkin(3) * t559 + pkin(10) * t509 + t510 * t577;
t552 = t488 * t423 - t424 * t489 + t565;
t551 = t532 * t554;
t457 = qJD(4) * t507 + t489;
t511 = (pkin(2) * t536 - pkin(9) * t540) * t532;
t513 = qJD(1) * (pkin(1) * t541 + pkin(8) * t574);
t550 = t523 * t474 - t511 * t522 + t513;
t438 = qJD(5) * t507 + t457;
t483 = (-qJD(5) + t560) * t572 + t523;
t393 = pkin(11) * t507 + t508 * t564 - t555 * t571;
t394 = pkin(11) * t509 + t510 * t564 + t555 * t574;
t549 = t456 * t393 - t394 * t457 + t552;
t548 = -t473 * t523 - t511 * t556 - t562;
t465 = pkin(3) * t570 + (-pkin(10) * t540 + t536 * t577) * t532;
t512 = -qJD(3) * t572 + t523;
t547 = t512 * t424 - t465 * t488 + t550;
t546 = -t423 * t512 + t489 * t465 + t548;
t435 = t555 * t533 + (-pkin(11) * t540 + t536 * t564) * t532;
t490 = t560 * t572 + t523;
t545 = t490 * t394 - t435 * t456 + t547;
t544 = -t393 * t490 + t457 * t435 + t546;
t538 = cos(qJ(6));
t534 = sin(qJ(6));
t524 = sin(t557);
t519 = rSges(2,1) * t541 - rSges(2,2) * t537;
t518 = rSges(2,1) * t537 + rSges(2,2) * t541;
t506 = t536 * t573 + t570;
t505 = t533 * t539 - t535 * t575;
t498 = t527 * t533 + t528 * t575;
t497 = -t527 * t575 + t528 * t533;
t496 = rSges(3,3) * t533 + (rSges(3,1) * t536 + rSges(3,2) * t540) * t532;
t495 = Icges(3,5) * t533 + (Icges(3,1) * t536 + Icges(3,4) * t540) * t532;
t494 = Icges(3,6) * t533 + (Icges(3,4) * t536 + Icges(3,2) * t540) * t532;
t493 = Icges(3,3) * t533 + (Icges(3,5) * t536 + Icges(3,6) * t540) * t532;
t492 = t533 * t524 + t536 * t551;
t491 = t524 * t575 - t533 * t554;
t487 = t510 * t539 + t559;
t486 = -t510 * t535 + t537 * t573;
t485 = t508 * t539 - t558;
t484 = -t508 * t535 - t539 * t571;
t482 = t510 * t528 + t527 * t574;
t481 = -t510 * t527 + t528 * t574;
t480 = t508 * t528 - t527 * t571;
t479 = -t508 * t527 - t528 * t571;
t478 = t510 * t554 + t524 * t574;
t477 = t510 * t524 - t537 * t551;
t476 = t508 * t554 - t524 * t571;
t475 = t508 * t524 + t541 * t551;
t472 = t492 * t538 - t534 * t572;
t471 = -t492 * t534 - t538 * t572;
t468 = rSges(3,1) * t510 - rSges(3,2) * t509 + rSges(3,3) * t574;
t467 = rSges(3,1) * t508 - rSges(3,2) * t507 - rSges(3,3) * t571;
t464 = Icges(3,1) * t510 - Icges(3,4) * t509 + Icges(3,5) * t574;
t463 = Icges(3,1) * t508 - Icges(3,4) * t507 - Icges(3,5) * t571;
t462 = Icges(3,4) * t510 - Icges(3,2) * t509 + Icges(3,6) * t574;
t461 = Icges(3,4) * t508 - Icges(3,2) * t507 - Icges(3,6) * t571;
t458 = rSges(4,1) * t506 + rSges(4,2) * t505 - rSges(4,3) * t572;
t455 = Icges(4,1) * t506 + Icges(4,4) * t505 - Icges(4,5) * t572;
t454 = Icges(4,4) * t506 + Icges(4,2) * t505 - Icges(4,6) * t572;
t453 = Icges(4,5) * t506 + Icges(4,6) * t505 - Icges(4,3) * t572;
t452 = pkin(5) * t492 + pkin(12) * t491;
t451 = rSges(5,1) * t498 + rSges(5,2) * t497 - rSges(5,3) * t572;
t450 = Icges(5,1) * t498 + Icges(5,4) * t497 - Icges(5,5) * t572;
t449 = Icges(5,4) * t498 + Icges(5,2) * t497 - Icges(5,6) * t572;
t448 = Icges(5,5) * t498 + Icges(5,6) * t497 - Icges(5,3) * t572;
t447 = qJD(6) * t491 + t483;
t446 = rSges(6,1) * t492 - rSges(6,2) * t491 - rSges(6,3) * t572;
t445 = Icges(6,1) * t492 - Icges(6,4) * t491 - Icges(6,5) * t572;
t444 = Icges(6,4) * t492 - Icges(6,2) * t491 - Icges(6,6) * t572;
t443 = Icges(6,5) * t492 - Icges(6,6) * t491 - Icges(6,3) * t572;
t442 = t478 * t538 + t509 * t534;
t441 = -t478 * t534 + t509 * t538;
t440 = t476 * t538 + t507 * t534;
t439 = -t476 * t534 + t507 * t538;
t434 = pkin(5) * t478 + pkin(12) * t477;
t433 = pkin(5) * t476 + pkin(12) * t475;
t432 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t509;
t431 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t507;
t430 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t509;
t429 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t507;
t428 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t509;
t427 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t507;
t426 = Icges(4,5) * t487 + Icges(4,6) * t486 + Icges(4,3) * t509;
t425 = Icges(4,5) * t485 + Icges(4,6) * t484 + Icges(4,3) * t507;
t422 = rSges(5,1) * t482 + rSges(5,2) * t481 + rSges(5,3) * t509;
t421 = rSges(5,1) * t480 + rSges(5,2) * t479 + rSges(5,3) * t507;
t420 = Icges(5,1) * t482 + Icges(5,4) * t481 + Icges(5,5) * t509;
t419 = Icges(5,1) * t480 + Icges(5,4) * t479 + Icges(5,5) * t507;
t418 = Icges(5,4) * t482 + Icges(5,2) * t481 + Icges(5,6) * t509;
t417 = Icges(5,4) * t480 + Icges(5,2) * t479 + Icges(5,6) * t507;
t416 = Icges(5,5) * t482 + Icges(5,6) * t481 + Icges(5,3) * t509;
t415 = Icges(5,5) * t480 + Icges(5,6) * t479 + Icges(5,3) * t507;
t414 = rSges(6,1) * t478 - rSges(6,2) * t477 + rSges(6,3) * t509;
t413 = rSges(6,1) * t476 - rSges(6,2) * t475 + rSges(6,3) * t507;
t412 = Icges(6,1) * t478 - Icges(6,4) * t477 + Icges(6,5) * t509;
t411 = Icges(6,1) * t476 - Icges(6,4) * t475 + Icges(6,5) * t507;
t410 = Icges(6,4) * t478 - Icges(6,2) * t477 + Icges(6,6) * t509;
t409 = Icges(6,4) * t476 - Icges(6,2) * t475 + Icges(6,6) * t507;
t408 = Icges(6,5) * t478 - Icges(6,6) * t477 + Icges(6,3) * t509;
t407 = Icges(6,5) * t476 - Icges(6,6) * t475 + Icges(6,3) * t507;
t406 = t468 * t523 - t496 * t522 + t513;
t405 = -t467 * t523 - t496 * t556 - t562;
t404 = qJD(6) * t475 + t438;
t403 = qJD(6) * t477 + t437;
t402 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t491;
t401 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t491;
t400 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t491;
t399 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t491;
t397 = (t467 * t537 + t468 * t541) * t561;
t391 = rSges(7,1) * t442 + rSges(7,2) * t441 + rSges(7,3) * t477;
t390 = rSges(7,1) * t440 + rSges(7,2) * t439 + rSges(7,3) * t475;
t389 = Icges(7,1) * t442 + Icges(7,4) * t441 + Icges(7,5) * t477;
t388 = Icges(7,1) * t440 + Icges(7,4) * t439 + Icges(7,5) * t475;
t387 = Icges(7,4) * t442 + Icges(7,2) * t441 + Icges(7,6) * t477;
t386 = Icges(7,4) * t440 + Icges(7,2) * t439 + Icges(7,6) * t475;
t385 = Icges(7,5) * t442 + Icges(7,6) * t441 + Icges(7,3) * t477;
t384 = Icges(7,5) * t440 + Icges(7,6) * t439 + Icges(7,3) * t475;
t382 = t432 * t512 - t458 * t488 + t550;
t381 = -t431 * t512 + t458 * t489 + t548;
t380 = t431 * t488 - t432 * t489 + t565;
t379 = t422 * t490 - t451 * t456 + t547;
t378 = -t421 * t490 + t451 * t457 + t546;
t377 = t421 * t456 - t422 * t457 + t552;
t376 = t414 * t483 - t437 * t446 + t545;
t375 = -t413 * t483 + t438 * t446 + t544;
t374 = t413 * t437 - t414 * t438 + t549;
t373 = t391 * t447 - t402 * t403 + t434 * t483 - t437 * t452 + t545;
t372 = -t390 * t447 + t402 * t404 - t433 * t483 + t438 * t452 + t544;
t371 = t390 * t403 - t391 * t404 + t433 * t437 - t434 * t438 + t549;
t1 = m(7) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(6) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(5) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(4) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(3) * (t397 ^ 2 + t405 ^ 2 + t406 ^ 2) / 0.2e1 + t447 * ((t385 * t491 + t387 * t471 + t389 * t472) * t403 + (t384 * t491 + t386 * t471 + t388 * t472) * t404 + (t491 * t399 + t471 * t400 + t472 * t401) * t447) / 0.2e1 + t403 * ((t477 * t385 + t441 * t387 + t442 * t389) * t403 + (t384 * t477 + t386 * t441 + t388 * t442) * t404 + (t399 * t477 + t400 * t441 + t401 * t442) * t447) / 0.2e1 + t404 * ((t385 * t475 + t387 * t439 + t389 * t440) * t403 + (t475 * t384 + t439 * t386 + t440 * t388) * t404 + (t399 * t475 + t400 * t439 + t401 * t440) * t447) / 0.2e1 + t438 * ((t408 * t507 - t410 * t475 + t412 * t476) * t437 + (t507 * t407 - t475 * t409 + t476 * t411) * t438 + (t443 * t507 - t444 * t475 + t445 * t476) * t483) / 0.2e1 + t437 * ((t509 * t408 - t477 * t410 + t478 * t412) * t437 + (t407 * t509 - t409 * t477 + t411 * t478) * t438 + (t443 * t509 - t444 * t477 + t445 * t478) * t483) / 0.2e1 + t483 * ((-t408 * t572 - t410 * t491 + t412 * t492) * t437 + (-t407 * t572 - t409 * t491 + t411 * t492) * t438 + (-t443 * t572 - t444 * t491 + t445 * t492) * t483) / 0.2e1 + t490 * ((-t416 * t572 + t418 * t497 + t420 * t498) * t456 + (-t415 * t572 + t417 * t497 + t419 * t498) * t457 + (-t448 * t572 + t449 * t497 + t450 * t498) * t490) / 0.2e1 + t456 * ((t509 * t416 + t481 * t418 + t482 * t420) * t456 + (t415 * t509 + t417 * t481 + t419 * t482) * t457 + (t448 * t509 + t449 * t481 + t450 * t482) * t490) / 0.2e1 + t457 * ((t416 * t507 + t418 * t479 + t420 * t480) * t456 + (t507 * t415 + t479 * t417 + t480 * t419) * t457 + (t448 * t507 + t449 * t479 + t450 * t480) * t490) / 0.2e1 + t512 * ((-t426 * t572 + t428 * t505 + t430 * t506) * t488 + (-t425 * t572 + t427 * t505 + t429 * t506) * t489 + (-t453 * t572 + t454 * t505 + t455 * t506) * t512) / 0.2e1 + t488 * ((t426 * t509 + t428 * t486 + t430 * t487) * t488 + (t425 * t509 + t427 * t486 + t429 * t487) * t489 + (t453 * t509 + t454 * t486 + t455 * t487) * t512) / 0.2e1 + t489 * ((t426 * t507 + t428 * t484 + t430 * t485) * t488 + (t425 * t507 + t427 * t484 + t429 * t485) * t489 + (t453 * t507 + t454 * t484 + t455 * t485) * t512) / 0.2e1 + t523 * ((t533 * t493 + (t494 * t540 + t495 * t536) * t532) * t523 + (((t462 * t540 + t464 * t536) * t537 - (t461 * t540 + t463 * t536) * t541) * t532 - t553 * t533) * t561) / 0.2e1 - ((-t493 * t571 - t494 * t507 + t495 * t508) * t523 + ((-t462 * t507 + t464 * t508) * t537 + (t507 * t461 - t508 * t463 + t580) * t541) * t561) * t556 / 0.2e1 + ((t493 * t574 - t494 * t509 + t495 * t510) * t523 + (-(-t461 * t509 + t463 * t510) * t541 + (-t509 * t462 + t510 * t464 - t580) * t537) * t561) * t522 / 0.2e1 + (Icges(2,3) + m(2) * (t518 ^ 2 + t519 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
