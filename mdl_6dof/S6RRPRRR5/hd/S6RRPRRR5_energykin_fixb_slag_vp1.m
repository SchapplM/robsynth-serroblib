% Calculate kinetic energy for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:01
% EndTime: 2019-03-09 13:38:04
% DurationCPUTime: 3.34s
% Computational Cost: add. (3767->385), mult. (9015->593), div. (0->0), fcn. (11460->14), ass. (0->175)
t544 = sin(qJ(1));
t546 = cos(qJ(1));
t543 = sin(qJ(2));
t583 = sin(pkin(12));
t585 = cos(pkin(6));
t559 = t585 * t583;
t584 = cos(pkin(12));
t560 = t585 * t584;
t590 = cos(qJ(2));
t552 = -t543 * t559 + t560 * t590;
t554 = t543 * t584 + t583 * t590;
t498 = -t544 * t554 + t546 * t552;
t514 = t543 * t560 + t559 * t590;
t522 = -t543 * t583 + t590 * t584;
t499 = t514 * t546 + t522 * t544;
t540 = sin(pkin(6));
t578 = t540 * t546;
t447 = Icges(4,5) * t499 + Icges(4,6) * t498 - Icges(4,3) * t578;
t500 = -t544 * t552 - t546 * t554;
t501 = -t514 * t544 + t522 * t546;
t579 = t540 * t544;
t448 = Icges(4,5) * t501 + Icges(4,6) * t500 + Icges(4,3) * t579;
t562 = t585 * t590;
t516 = -t544 * t543 + t546 * t562;
t569 = t543 * t585;
t517 = t544 * t590 + t546 * t569;
t485 = Icges(3,5) * t517 + Icges(3,6) * t516 - Icges(3,3) * t578;
t518 = -t546 * t543 - t544 * t562;
t519 = -t544 * t569 + t546 * t590;
t486 = Icges(3,5) * t519 + Icges(3,6) * t518 + Icges(3,3) * t579;
t595 = t540 * ((t447 + t485) * t546 + (-t448 - t486) * t544);
t512 = t522 * t540;
t513 = t554 * t540;
t594 = Icges(4,5) * t513 + Icges(4,6) * t512 + (Icges(3,5) * t543 + Icges(3,6) * t590) * t540 + (Icges(4,3) + Icges(3,3)) * t585;
t589 = cos(qJ(4));
t588 = pkin(2) * t590;
t545 = cos(qJ(5));
t587 = pkin(5) * t545;
t541 = sin(qJ(5));
t582 = t498 * t541;
t581 = t500 * t541;
t580 = t512 * t541;
t570 = pkin(2) * t569 - qJ(3) * t540;
t497 = -t544 * t570 + t546 * t588;
t520 = qJD(1) * (pkin(1) * t546 + pkin(8) * t579);
t531 = qJD(2) * t585 + qJD(1);
t577 = t531 * t497 + t520;
t575 = qJD(2) * t540;
t530 = t544 * t575;
t476 = -qJD(4) * t500 + t530;
t576 = qJD(1) * (pkin(1) * t544 - pkin(8) * t578);
t574 = qJD(3) * t546;
t496 = t544 * t588 + t546 * t570;
t571 = t546 * t575;
t573 = qJD(3) * t585 + t496 * t530 + t497 * t571;
t542 = sin(qJ(4));
t572 = t540 * t589;
t474 = t501 * t542 - t544 * t572;
t433 = qJD(5) * t474 + t476;
t505 = -qJD(4) * t512 + t531;
t523 = t540 * t543 * pkin(2) + qJ(3) * t585;
t567 = qJD(2) * (-t513 * rSges(4,1) - t512 * rSges(4,2) - rSges(4,3) * t585 - t523);
t566 = qJD(2) * (-pkin(3) * t513 + pkin(9) * t512 - t523);
t565 = qJD(3) * t579 - t576;
t503 = t513 * t542 - t585 * t589;
t465 = qJD(5) * t503 + t505;
t477 = -qJD(4) * t498 - t571;
t462 = pkin(3) * t499 - pkin(9) * t498;
t463 = pkin(3) * t501 - pkin(9) * t500;
t558 = t462 * t530 + t463 * t571 + t573;
t472 = t499 * t542 + t546 * t572;
t434 = qJD(5) * t472 + t477;
t473 = t499 * t589 - t542 * t578;
t431 = pkin(4) * t473 + pkin(10) * t472;
t475 = t501 * t589 + t542 * t579;
t432 = pkin(4) * t475 + pkin(10) * t474;
t555 = t476 * t431 - t432 * t477 + t558;
t553 = t531 * t463 + (t544 * t566 - t574) * t540 + t577;
t551 = (-t462 - t496) * t531 + t566 * t578 + t565;
t504 = t513 * t589 + t542 * t585;
t464 = pkin(4) * t504 + pkin(10) * t503;
t550 = t505 * t432 - t464 * t476 + t553;
t549 = -t431 * t505 + t477 * t464 + t551;
t539 = qJ(5) + qJ(6);
t538 = cos(t539);
t537 = sin(t539);
t526 = rSges(2,1) * t546 - rSges(2,2) * t544;
t525 = rSges(2,1) * t544 + rSges(2,2) * t546;
t511 = t585 * rSges(3,3) + (rSges(3,1) * t543 + rSges(3,2) * t590) * t540;
t510 = Icges(3,5) * t585 + (Icges(3,1) * t543 + Icges(3,4) * t590) * t540;
t509 = Icges(3,6) * t585 + (Icges(3,4) * t543 + Icges(3,2) * t590) * t540;
t492 = rSges(3,1) * t519 + rSges(3,2) * t518 + rSges(3,3) * t579;
t491 = rSges(3,1) * t517 + rSges(3,2) * t516 - rSges(3,3) * t578;
t490 = Icges(3,1) * t519 + Icges(3,4) * t518 + Icges(3,5) * t579;
t489 = Icges(3,1) * t517 + Icges(3,4) * t516 - Icges(3,5) * t578;
t488 = Icges(3,4) * t519 + Icges(3,2) * t518 + Icges(3,6) * t579;
t487 = Icges(3,4) * t517 + Icges(3,2) * t516 - Icges(3,6) * t578;
t483 = Icges(4,1) * t513 + Icges(4,4) * t512 + Icges(4,5) * t585;
t482 = Icges(4,4) * t513 + Icges(4,2) * t512 + Icges(4,6) * t585;
t469 = t504 * t545 - t580;
t468 = -t504 * t541 - t512 * t545;
t467 = t504 * t538 - t512 * t537;
t466 = -t504 * t537 - t512 * t538;
t461 = rSges(5,1) * t504 - rSges(5,2) * t503 - rSges(5,3) * t512;
t460 = Icges(5,1) * t504 - Icges(5,4) * t503 - Icges(5,5) * t512;
t459 = Icges(5,4) * t504 - Icges(5,2) * t503 - Icges(5,6) * t512;
t458 = Icges(5,5) * t504 - Icges(5,6) * t503 - Icges(5,3) * t512;
t455 = rSges(4,1) * t501 + rSges(4,2) * t500 + rSges(4,3) * t579;
t454 = rSges(4,1) * t499 + rSges(4,2) * t498 - rSges(4,3) * t578;
t452 = Icges(4,1) * t501 + Icges(4,4) * t500 + Icges(4,5) * t579;
t451 = Icges(4,1) * t499 + Icges(4,4) * t498 - Icges(4,5) * t578;
t450 = Icges(4,4) * t501 + Icges(4,2) * t500 + Icges(4,6) * t579;
t449 = Icges(4,4) * t499 + Icges(4,2) * t498 - Icges(4,6) * t578;
t446 = qJD(6) * t503 + t465;
t445 = t492 * t531 - t511 * t530 + t520;
t444 = -t491 * t531 - t511 * t571 - t576;
t443 = t475 * t545 - t581;
t442 = -t475 * t541 - t500 * t545;
t441 = t473 * t545 - t582;
t440 = -t473 * t541 - t498 * t545;
t439 = (t491 * t544 + t492 * t546) * t575;
t438 = t475 * t538 - t500 * t537;
t437 = -t475 * t537 - t500 * t538;
t436 = t473 * t538 - t498 * t537;
t435 = -t473 * t537 - t498 * t538;
t428 = rSges(6,1) * t469 + rSges(6,2) * t468 + rSges(6,3) * t503;
t427 = Icges(6,1) * t469 + Icges(6,4) * t468 + Icges(6,5) * t503;
t426 = Icges(6,4) * t469 + Icges(6,2) * t468 + Icges(6,6) * t503;
t425 = Icges(6,5) * t469 + Icges(6,6) * t468 + Icges(6,3) * t503;
t424 = rSges(5,1) * t475 - rSges(5,2) * t474 - rSges(5,3) * t500;
t423 = rSges(5,1) * t473 - rSges(5,2) * t472 - rSges(5,3) * t498;
t422 = Icges(5,1) * t475 - Icges(5,4) * t474 - Icges(5,5) * t500;
t421 = Icges(5,1) * t473 - Icges(5,4) * t472 - Icges(5,5) * t498;
t420 = Icges(5,4) * t475 - Icges(5,2) * t474 - Icges(5,6) * t500;
t419 = Icges(5,4) * t473 - Icges(5,2) * t472 - Icges(5,6) * t498;
t418 = Icges(5,5) * t475 - Icges(5,6) * t474 - Icges(5,3) * t500;
t417 = Icges(5,5) * t473 - Icges(5,6) * t472 - Icges(5,3) * t498;
t416 = -pkin(5) * t580 + pkin(11) * t503 + t504 * t587;
t415 = rSges(7,1) * t467 + rSges(7,2) * t466 + rSges(7,3) * t503;
t414 = Icges(7,1) * t467 + Icges(7,4) * t466 + Icges(7,5) * t503;
t413 = Icges(7,4) * t467 + Icges(7,2) * t466 + Icges(7,6) * t503;
t412 = Icges(7,5) * t467 + Icges(7,6) * t466 + Icges(7,3) * t503;
t411 = qJD(6) * t472 + t434;
t410 = qJD(6) * t474 + t433;
t408 = t455 * t531 + (t544 * t567 - t574) * t540 + t577;
t407 = (-t454 - t496) * t531 + t567 * t578 + t565;
t406 = rSges(6,1) * t443 + rSges(6,2) * t442 + rSges(6,3) * t474;
t405 = rSges(6,1) * t441 + rSges(6,2) * t440 + rSges(6,3) * t472;
t404 = Icges(6,1) * t443 + Icges(6,4) * t442 + Icges(6,5) * t474;
t403 = Icges(6,1) * t441 + Icges(6,4) * t440 + Icges(6,5) * t472;
t402 = Icges(6,4) * t443 + Icges(6,2) * t442 + Icges(6,6) * t474;
t401 = Icges(6,4) * t441 + Icges(6,2) * t440 + Icges(6,6) * t472;
t400 = Icges(6,5) * t443 + Icges(6,6) * t442 + Icges(6,3) * t474;
t399 = Icges(6,5) * t441 + Icges(6,6) * t440 + Icges(6,3) * t472;
t398 = rSges(7,1) * t438 + rSges(7,2) * t437 + rSges(7,3) * t474;
t397 = rSges(7,1) * t436 + rSges(7,2) * t435 + rSges(7,3) * t472;
t396 = Icges(7,1) * t438 + Icges(7,4) * t437 + Icges(7,5) * t474;
t395 = Icges(7,1) * t436 + Icges(7,4) * t435 + Icges(7,5) * t472;
t394 = Icges(7,4) * t438 + Icges(7,2) * t437 + Icges(7,6) * t474;
t393 = Icges(7,4) * t436 + Icges(7,2) * t435 + Icges(7,6) * t472;
t392 = Icges(7,5) * t438 + Icges(7,6) * t437 + Icges(7,3) * t474;
t391 = Icges(7,5) * t436 + Icges(7,6) * t435 + Icges(7,3) * t472;
t390 = -pkin(5) * t581 + pkin(11) * t474 + t475 * t587;
t389 = -pkin(5) * t582 + pkin(11) * t472 + t473 * t587;
t388 = (t454 * t544 + t455 * t546) * t575 + t573;
t387 = t424 * t505 - t461 * t476 + t553;
t386 = -t423 * t505 + t461 * t477 + t551;
t385 = t423 * t476 - t424 * t477 + t558;
t384 = t406 * t465 - t428 * t433 + t550;
t383 = -t405 * t465 + t428 * t434 + t549;
t382 = t405 * t433 - t406 * t434 + t555;
t381 = t390 * t465 + t398 * t446 - t410 * t415 - t416 * t433 + t550;
t380 = -t389 * t465 - t397 * t446 + t411 * t415 + t416 * t434 + t549;
t379 = t389 * t433 - t390 * t434 + t397 * t410 - t398 * t411 + t555;
t1 = t476 * ((-t418 * t500 - t420 * t474 + t422 * t475) * t476 + (-t417 * t500 - t419 * t474 + t421 * t475) * t477 + (-t458 * t500 - t459 * t474 + t460 * t475) * t505) / 0.2e1 + t477 * ((-t418 * t498 - t420 * t472 + t422 * t473) * t476 + (-t417 * t498 - t419 * t472 + t421 * t473) * t477 + (-t458 * t498 - t459 * t472 + t460 * t473) * t505) / 0.2e1 + t505 * ((-t418 * t512 - t420 * t503 + t422 * t504) * t476 + (-t417 * t512 - t419 * t503 + t421 * t504) * t477 + (-t458 * t512 - t459 * t503 + t460 * t504) * t505) / 0.2e1 + t433 * ((t400 * t474 + t402 * t442 + t404 * t443) * t433 + (t399 * t474 + t401 * t442 + t403 * t443) * t434 + (t425 * t474 + t426 * t442 + t427 * t443) * t465) / 0.2e1 + t434 * ((t400 * t472 + t402 * t440 + t404 * t441) * t433 + (t399 * t472 + t401 * t440 + t403 * t441) * t434 + (t425 * t472 + t426 * t440 + t427 * t441) * t465) / 0.2e1 + t465 * ((t400 * t503 + t402 * t468 + t404 * t469) * t433 + (t399 * t503 + t401 * t468 + t403 * t469) * t434 + (t425 * t503 + t426 * t468 + t427 * t469) * t465) / 0.2e1 + t410 * ((t474 * t392 + t437 * t394 + t438 * t396) * t410 + (t391 * t474 + t393 * t437 + t395 * t438) * t411 + (t412 * t474 + t413 * t437 + t414 * t438) * t446) / 0.2e1 + t411 * ((t392 * t472 + t394 * t435 + t396 * t436) * t410 + (t391 * t472 + t393 * t435 + t395 * t436) * t411 + (t412 * t472 + t413 * t435 + t414 * t436) * t446) / 0.2e1 + t446 * ((t392 * t503 + t394 * t466 + t396 * t467) * t410 + (t391 * t503 + t393 * t466 + t395 * t467) * t411 + (t412 * t503 + t413 * t466 + t414 * t467) * t446) / 0.2e1 + m(3) * (t439 ^ 2 + t444 ^ 2 + t445 ^ 2) / 0.2e1 + m(4) * (t388 ^ 2 + t407 ^ 2 + t408 ^ 2) / 0.2e1 + m(5) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(6) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 + m(7) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + ((t585 * t486 + (t488 * t590 + t490 * t543) * t540) * t530 - (t585 * t485 + (t487 * t590 + t489 * t543) * t540) * t571 + ((t448 * t585 + t512 * t450 + t513 * t452) * t544 - (t447 * t585 + t512 * t449 + t513 * t451) * t546) * t575 + ((t509 * t590 + t510 * t543) * t540 + t512 * t482 + t513 * t483 + t594 * t585) * t531) * t531 / 0.2e1 + (Icges(2,3) + m(2) * (t525 ^ 2 + t526 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t449 * t500 - t451 * t501 - t487 * t518 - t489 * t519) * t546 + (t500 * t450 + t501 * t452 + t518 * t488 + t519 * t490 - t595) * t544) * t575 + (t482 * t500 + t483 * t501 + t509 * t518 + t510 * t519 + t579 * t594) * t531) * t530 / 0.2e1 - (((-t498 * t449 - t499 * t451 - t516 * t487 - t517 * t489 + t595) * t546 + (t450 * t498 + t452 * t499 + t488 * t516 + t490 * t517) * t544) * t575 + (t482 * t498 + t483 * t499 + t509 * t516 + t510 * t517 - t578 * t594) * t531) * t571 / 0.2e1;
T  = t1;
