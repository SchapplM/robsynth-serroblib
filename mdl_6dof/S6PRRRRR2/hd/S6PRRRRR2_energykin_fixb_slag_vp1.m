% Calculate kinetic energy for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:24
% EndTime: 2019-03-09 00:42:26
% DurationCPUTime: 2.35s
% Computational Cost: add. (3512->381), mult. (6067->600), div. (0->0), fcn. (7350->14), ass. (0->175)
t567 = qJD(2) ^ 2;
t530 = cos(qJ(3));
t566 = pkin(3) * t530;
t529 = cos(qJ(5));
t565 = pkin(5) * t529;
t522 = sin(pkin(12));
t524 = cos(pkin(12));
t528 = sin(qJ(2));
t525 = cos(pkin(6));
t531 = cos(qJ(2));
t552 = t525 * t531;
t503 = t522 * t528 - t524 * t552;
t526 = sin(qJ(5));
t562 = t503 * t526;
t505 = t522 * t552 + t524 * t528;
t561 = t505 * t526;
t523 = sin(pkin(6));
t560 = t522 * t523;
t559 = t523 * t524;
t527 = sin(qJ(3));
t558 = t523 * t527;
t557 = t523 * t528;
t556 = t523 * t530;
t555 = t523 * t531;
t554 = t525 * t527;
t553 = t525 * t528;
t551 = qJ(3) + qJ(4);
t550 = qJD(2) * t523;
t514 = t522 * t550;
t489 = qJD(3) * t505 + t514;
t517 = qJD(2) * t525;
t548 = t522 * t558;
t547 = t524 * t558;
t546 = t526 * t555;
t453 = qJD(4) * t505 + t489;
t545 = t524 * t550;
t504 = t522 * t531 + t524 * t553;
t475 = t504 * pkin(2) + t503 * pkin(8);
t506 = -t522 * t553 + t524 * t531;
t476 = t506 * pkin(2) + t505 * pkin(8);
t544 = t475 * t514 + t476 * t545 + qJD(1);
t543 = cos(t551);
t519 = sin(t551);
t542 = t523 * t543;
t481 = t506 * t519 - t522 * t542;
t433 = qJD(5) * t481 + t453;
t490 = qJD(3) * t503 - t545;
t509 = (pkin(2) * t528 - pkin(8) * t531) * t523;
t541 = t476 * t517 - t509 * t514;
t454 = qJD(4) * t503 + t490;
t492 = t517 + (-qJD(3) - qJD(4)) * t555;
t479 = t504 * t519 + t524 * t542;
t434 = qJD(5) * t479 + t454;
t423 = -pkin(3) * t547 + pkin(9) * t503 + t504 * t566;
t424 = pkin(3) * t548 + pkin(9) * t505 + t506 * t566;
t540 = t489 * t423 - t424 * t490 + t544;
t497 = t519 * t557 - t525 * t543;
t471 = qJD(5) * t497 + t492;
t539 = (-t475 * t525 - t509 * t559) * qJD(2);
t470 = pkin(3) * t554 + (-pkin(9) * t531 + t528 * t566) * t523;
t510 = -qJD(3) * t555 + t517;
t538 = t510 * t424 - t470 * t489 + t541;
t480 = t504 * t543 - t519 * t559;
t438 = pkin(4) * t480 + pkin(10) * t479;
t482 = t506 * t543 + t519 * t560;
t439 = pkin(4) * t482 + pkin(10) * t481;
t537 = t453 * t438 - t439 * t454 + t540;
t536 = -t423 * t510 + t490 * t470 + t539;
t498 = t525 * t519 + t528 * t542;
t469 = pkin(4) * t498 + pkin(10) * t497;
t535 = t492 * t439 - t453 * t469 + t538;
t534 = -t438 * t492 + t454 * t469 + t536;
t521 = qJ(5) + qJ(6);
t520 = cos(t521);
t518 = sin(t521);
t508 = t528 * t556 + t554;
t507 = t525 * t530 - t527 * t557;
t496 = rSges(3,3) * t525 + (rSges(3,1) * t528 + rSges(3,2) * t531) * t523;
t495 = Icges(3,5) * t525 + (Icges(3,1) * t528 + Icges(3,4) * t531) * t523;
t494 = Icges(3,6) * t525 + (Icges(3,4) * t528 + Icges(3,2) * t531) * t523;
t493 = Icges(3,3) * t525 + (Icges(3,5) * t528 + Icges(3,6) * t531) * t523;
t488 = t506 * t530 + t548;
t487 = -t506 * t527 + t522 * t556;
t486 = t504 * t530 - t547;
t485 = -t504 * t527 - t524 * t556;
t484 = t498 * t529 - t546;
t483 = -t498 * t526 - t529 * t555;
t474 = t498 * t520 - t518 * t555;
t473 = -t498 * t518 - t520 * t555;
t468 = rSges(4,1) * t508 + rSges(4,2) * t507 - rSges(4,3) * t555;
t467 = Icges(4,1) * t508 + Icges(4,4) * t507 - Icges(4,5) * t555;
t466 = Icges(4,4) * t508 + Icges(4,2) * t507 - Icges(4,6) * t555;
t465 = Icges(4,5) * t508 + Icges(4,6) * t507 - Icges(4,3) * t555;
t462 = rSges(3,1) * t506 - rSges(3,2) * t505 + rSges(3,3) * t560;
t461 = rSges(3,1) * t504 - rSges(3,2) * t503 - rSges(3,3) * t559;
t460 = Icges(3,1) * t506 - Icges(3,4) * t505 + Icges(3,5) * t560;
t459 = Icges(3,1) * t504 - Icges(3,4) * t503 - Icges(3,5) * t559;
t458 = Icges(3,4) * t506 - Icges(3,2) * t505 + Icges(3,6) * t560;
t457 = Icges(3,4) * t504 - Icges(3,2) * t503 - Icges(3,6) * t559;
t456 = Icges(3,5) * t506 - Icges(3,6) * t505 + Icges(3,3) * t560;
t455 = Icges(3,5) * t504 - Icges(3,6) * t503 - Icges(3,3) * t559;
t452 = rSges(5,1) * t498 - rSges(5,2) * t497 - rSges(5,3) * t555;
t451 = Icges(5,1) * t498 - Icges(5,4) * t497 - Icges(5,5) * t555;
t450 = Icges(5,4) * t498 - Icges(5,2) * t497 - Icges(5,6) * t555;
t449 = Icges(5,5) * t498 - Icges(5,6) * t497 - Icges(5,3) * t555;
t448 = t482 * t529 + t561;
t447 = -t482 * t526 + t505 * t529;
t446 = t480 * t529 + t562;
t445 = -t480 * t526 + t503 * t529;
t444 = t482 * t520 + t505 * t518;
t443 = -t482 * t518 + t505 * t520;
t442 = t480 * t520 + t503 * t518;
t441 = -t480 * t518 + t503 * t520;
t440 = qJD(6) * t497 + t471;
t436 = (-t461 * t525 - t496 * t559) * qJD(2);
t435 = (t462 * t525 - t496 * t560) * qJD(2);
t432 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t505;
t431 = rSges(4,1) * t486 + rSges(4,2) * t485 + rSges(4,3) * t503;
t430 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t505;
t429 = Icges(4,1) * t486 + Icges(4,4) * t485 + Icges(4,5) * t503;
t428 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t505;
t427 = Icges(4,4) * t486 + Icges(4,2) * t485 + Icges(4,6) * t503;
t426 = Icges(4,5) * t488 + Icges(4,6) * t487 + Icges(4,3) * t505;
t425 = Icges(4,5) * t486 + Icges(4,6) * t485 + Icges(4,3) * t503;
t422 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t505;
t421 = rSges(5,1) * t480 - rSges(5,2) * t479 + rSges(5,3) * t503;
t419 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t505;
t418 = Icges(5,1) * t480 - Icges(5,4) * t479 + Icges(5,5) * t503;
t417 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t505;
t416 = Icges(5,4) * t480 - Icges(5,2) * t479 + Icges(5,6) * t503;
t415 = Icges(5,5) * t482 - Icges(5,6) * t481 + Icges(5,3) * t505;
t414 = Icges(5,5) * t480 - Icges(5,6) * t479 + Icges(5,3) * t503;
t413 = rSges(6,1) * t484 + rSges(6,2) * t483 + rSges(6,3) * t497;
t412 = Icges(6,1) * t484 + Icges(6,4) * t483 + Icges(6,5) * t497;
t411 = Icges(6,4) * t484 + Icges(6,2) * t483 + Icges(6,6) * t497;
t410 = Icges(6,5) * t484 + Icges(6,6) * t483 + Icges(6,3) * t497;
t408 = rSges(7,1) * t474 + rSges(7,2) * t473 + rSges(7,3) * t497;
t407 = Icges(7,1) * t474 + Icges(7,4) * t473 + Icges(7,5) * t497;
t406 = Icges(7,4) * t474 + Icges(7,2) * t473 + Icges(7,6) * t497;
t405 = Icges(7,5) * t474 + Icges(7,6) * t473 + Icges(7,3) * t497;
t404 = -pkin(5) * t546 + pkin(11) * t497 + t498 * t565;
t402 = qJD(1) + (t461 * t522 + t462 * t524) * t550;
t400 = qJD(6) * t479 + t434;
t399 = qJD(6) * t481 + t433;
t397 = rSges(6,1) * t448 + rSges(6,2) * t447 + rSges(6,3) * t481;
t396 = rSges(6,1) * t446 + rSges(6,2) * t445 + rSges(6,3) * t479;
t395 = Icges(6,1) * t448 + Icges(6,4) * t447 + Icges(6,5) * t481;
t394 = Icges(6,1) * t446 + Icges(6,4) * t445 + Icges(6,5) * t479;
t393 = Icges(6,4) * t448 + Icges(6,2) * t447 + Icges(6,6) * t481;
t392 = Icges(6,4) * t446 + Icges(6,2) * t445 + Icges(6,6) * t479;
t391 = Icges(6,5) * t448 + Icges(6,6) * t447 + Icges(6,3) * t481;
t390 = Icges(6,5) * t446 + Icges(6,6) * t445 + Icges(6,3) * t479;
t389 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t481;
t388 = rSges(7,1) * t442 + rSges(7,2) * t441 + rSges(7,3) * t479;
t387 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t481;
t386 = Icges(7,1) * t442 + Icges(7,4) * t441 + Icges(7,5) * t479;
t385 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t481;
t384 = Icges(7,4) * t442 + Icges(7,2) * t441 + Icges(7,6) * t479;
t383 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t481;
t382 = Icges(7,5) * t442 + Icges(7,6) * t441 + Icges(7,3) * t479;
t381 = pkin(5) * t561 + pkin(11) * t481 + t482 * t565;
t380 = pkin(5) * t562 + pkin(11) * t479 + t480 * t565;
t379 = -t431 * t510 + t468 * t490 + t539;
t378 = t432 * t510 - t468 * t489 + t541;
t377 = t431 * t489 - t432 * t490 + t544;
t376 = -t421 * t492 + t452 * t454 + t536;
t375 = t422 * t492 - t452 * t453 + t538;
t374 = t421 * t453 - t422 * t454 + t540;
t373 = -t396 * t471 + t413 * t434 + t534;
t372 = t397 * t471 - t413 * t433 + t535;
t371 = t396 * t433 - t397 * t434 + t537;
t370 = -t380 * t471 - t388 * t440 + t400 * t408 + t404 * t434 + t534;
t369 = t381 * t471 + t389 * t440 - t399 * t408 - t404 * t433 + t535;
t368 = t380 * t433 - t381 * t434 + t399 * t388 - t400 * t389 + t537;
t1 = t489 * ((t426 * t505 + t428 * t487 + t430 * t488) * t489 + (t425 * t505 + t427 * t487 + t429 * t488) * t490 + (t465 * t505 + t466 * t487 + t467 * t488) * t510) / 0.2e1 + t490 * ((t426 * t503 + t428 * t485 + t430 * t486) * t489 + (t425 * t503 + t427 * t485 + t429 * t486) * t490 + (t465 * t503 + t466 * t485 + t467 * t486) * t510) / 0.2e1 + t510 * ((-t426 * t555 + t428 * t507 + t430 * t508) * t489 + (-t425 * t555 + t427 * t507 + t429 * t508) * t490 + (-t465 * t555 + t466 * t507 + t467 * t508) * t510) / 0.2e1 + m(4) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(3) * (t402 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(6) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(5) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + t440 * ((t383 * t497 + t385 * t473 + t387 * t474) * t399 + (t382 * t497 + t384 * t473 + t386 * t474) * t400 + (t497 * t405 + t473 * t406 + t474 * t407) * t440) / 0.2e1 + t400 * ((t383 * t479 + t385 * t441 + t387 * t442) * t399 + (t479 * t382 + t441 * t384 + t442 * t386) * t400 + (t405 * t479 + t406 * t441 + t407 * t442) * t440) / 0.2e1 + t399 * ((t481 * t383 + t443 * t385 + t444 * t387) * t399 + (t382 * t481 + t384 * t443 + t386 * t444) * t400 + (t405 * t481 + t406 * t443 + t407 * t444) * t440) / 0.2e1 + t471 * ((t391 * t497 + t393 * t483 + t395 * t484) * t433 + (t390 * t497 + t392 * t483 + t394 * t484) * t434 + (t410 * t497 + t411 * t483 + t412 * t484) * t471) / 0.2e1 + t433 * ((t481 * t391 + t447 * t393 + t448 * t395) * t433 + (t390 * t481 + t392 * t447 + t394 * t448) * t434 + (t410 * t481 + t411 * t447 + t412 * t448) * t471) / 0.2e1 + t434 * ((t391 * t479 + t393 * t445 + t395 * t446) * t433 + (t479 * t390 + t445 * t392 + t446 * t394) * t434 + (t410 * t479 + t411 * t445 + t412 * t446) * t471) / 0.2e1 + t454 * ((t415 * t503 - t417 * t479 + t419 * t480) * t453 + (t503 * t414 - t479 * t416 + t480 * t418) * t454 + (t449 * t503 - t450 * t479 + t451 * t480) * t492) / 0.2e1 + t492 * ((-t415 * t555 - t417 * t497 + t419 * t498) * t453 + (-t414 * t555 - t416 * t497 + t418 * t498) * t454 + (-t449 * t555 - t450 * t497 + t451 * t498) * t492) / 0.2e1 + t453 * ((t505 * t415 - t481 * t417 + t482 * t419) * t453 + (t414 * t505 - t416 * t481 + t418 * t482) * t454 + (t449 * t505 - t481 * t450 + t482 * t451) * t492) / 0.2e1 - t567 * ((-t456 * t559 - t458 * t503 + t460 * t504) * t560 - (-t455 * t559 - t457 * t503 + t459 * t504) * t559 + (-t493 * t559 - t494 * t503 + t495 * t504) * t525) * t559 / 0.2e1 + (t525 * (t525 ^ 2 * t493 + (((t458 * t531 + t460 * t528) * t522 - (t457 * t531 + t459 * t528) * t524) * t523 + (-t455 * t524 + t456 * t522 + t494 * t531 + t495 * t528) * t525) * t523) + ((t456 * t560 - t458 * t505 + t460 * t506) * t560 - (t455 * t560 - t457 * t505 + t459 * t506) * t559 + (t493 * t560 - t494 * t505 + t495 * t506) * t525) * t560) * t567 / 0.2e1;
T  = t1;
