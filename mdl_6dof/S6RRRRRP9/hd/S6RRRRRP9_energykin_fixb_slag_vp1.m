% Calculate kinetic energy for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:16
% EndTime: 2019-03-10 02:00:19
% DurationCPUTime: 3.08s
% Computational Cost: add. (3157->337), mult. (6859->517), div. (0->0), fcn. (8417->12), ass. (0->157)
t582 = Icges(6,1) + Icges(7,1);
t581 = Icges(6,4) + Icges(7,4);
t580 = Icges(6,5) + Icges(7,5);
t579 = Icges(6,2) + Icges(7,2);
t578 = Icges(6,6) + Icges(7,6);
t577 = Icges(6,3) + Icges(7,3);
t576 = rSges(7,3) + qJ(6);
t522 = sin(qJ(2));
t523 = sin(qJ(1));
t525 = cos(qJ(2));
t526 = cos(qJ(1));
t557 = cos(pkin(6));
t538 = t526 * t557;
t495 = t522 * t523 - t525 * t538;
t496 = t522 * t538 + t523 * t525;
t519 = sin(pkin(6));
t552 = t519 * t526;
t459 = Icges(3,5) * t496 - Icges(3,6) * t495 - Icges(3,3) * t552;
t539 = t523 * t557;
t497 = t526 * t522 + t525 * t539;
t498 = -t522 * t539 + t526 * t525;
t554 = t519 * t523;
t460 = Icges(3,5) * t498 - Icges(3,6) * t497 + Icges(3,3) * t554;
t575 = t519 * (t459 * t526 - t460 * t523);
t521 = sin(qJ(3));
t561 = cos(qJ(3));
t481 = t496 * t561 - t521 * t552;
t518 = qJ(4) + qJ(5);
t514 = sin(t518);
t515 = cos(t518);
t446 = -t481 * t514 + t495 * t515;
t447 = t481 * t515 + t495 * t514;
t542 = t519 * t561;
t480 = t496 * t521 + t526 * t542;
t574 = t446 * t578 + t447 * t580 + t480 * t577;
t483 = t498 * t561 + t521 * t554;
t448 = -t483 * t514 + t497 * t515;
t449 = t483 * t515 + t497 * t514;
t482 = t498 * t521 - t523 * t542;
t573 = t448 * t578 + t449 * t580 + t482 * t577;
t572 = t446 * t579 + t447 * t581 + t480 * t578;
t571 = t448 * t579 + t449 * t581 + t482 * t578;
t570 = t581 * t446 + t447 * t582 + t580 * t480;
t569 = t581 * t448 + t449 * t582 + t580 * t482;
t494 = t521 * t557 + t522 * t542;
t553 = t519 * t525;
t474 = -t494 * t514 - t515 * t553;
t475 = t494 * t515 - t514 * t553;
t493 = t519 * t521 * t522 - t557 * t561;
t568 = t474 * t578 + t475 * t580 + t493 * t577;
t567 = t474 * t579 + t475 * t581 + t493 * t578;
t566 = t581 * t474 + t475 * t582 + t580 * t493;
t524 = cos(qJ(4));
t559 = t524 * pkin(4);
t520 = sin(qJ(4));
t556 = t495 * t520;
t555 = t497 * t520;
t540 = pkin(5) * t514;
t547 = pkin(5) * t515;
t551 = rSges(7,1) * t447 + rSges(7,2) * t446 + t576 * t480 + t481 * t547 + t495 * t540;
t550 = rSges(7,1) * t449 + rSges(7,2) * t448 + t576 * t482 + t483 * t547 + t497 * t540;
t549 = rSges(7,1) * t475 + rSges(7,2) * t474 + t576 * t493 + t494 * t547 - t540 * t553;
t471 = pkin(2) * t496 + pkin(9) * t495;
t472 = pkin(2) * t498 + pkin(9) * t497;
t544 = qJD(2) * t519;
t509 = t523 * t544;
t541 = t526 * t544;
t548 = t471 * t509 + t472 * t541;
t484 = qJD(3) * t497 + t509;
t545 = qJD(1) * (pkin(1) * t523 - pkin(8) * t552);
t510 = qJD(2) * t557 + qJD(1);
t543 = t520 * t553;
t444 = qJD(4) * t482 + t484;
t485 = qJD(3) * t495 - t541;
t442 = pkin(3) * t481 + pkin(10) * t480;
t443 = pkin(3) * t483 + pkin(10) * t482;
t536 = t484 * t442 - t443 * t485 + t548;
t445 = qJD(4) * t480 + t485;
t500 = -qJD(3) * t553 + t510;
t499 = (pkin(2) * t522 - pkin(9) * t525) * t519;
t501 = qJD(1) * (pkin(1) * t526 + pkin(8) * t554);
t535 = t510 * t472 - t499 * t509 + t501;
t473 = qJD(4) * t493 + t500;
t385 = pkin(4) * t556 + pkin(11) * t480 + t481 * t559;
t386 = pkin(4) * t555 + pkin(11) * t482 + t483 * t559;
t534 = t444 * t385 - t386 * t445 + t536;
t533 = -t471 * t510 - t499 * t541 - t545;
t470 = pkin(3) * t494 + pkin(10) * t493;
t532 = t500 * t443 - t470 * t484 + t535;
t531 = -t442 * t500 + t485 * t470 + t533;
t427 = -pkin(4) * t543 + pkin(11) * t493 + t494 * t559;
t530 = t473 * t386 - t427 * t444 + t532;
t529 = -t385 * t473 + t445 * t427 + t531;
t506 = rSges(2,1) * t526 - rSges(2,2) * t523;
t505 = rSges(2,1) * t523 + rSges(2,2) * t526;
t489 = t557 * rSges(3,3) + (rSges(3,1) * t522 + rSges(3,2) * t525) * t519;
t488 = Icges(3,5) * t557 + (Icges(3,1) * t522 + Icges(3,4) * t525) * t519;
t487 = Icges(3,6) * t557 + (Icges(3,4) * t522 + Icges(3,2) * t525) * t519;
t486 = Icges(3,3) * t557 + (Icges(3,5) * t522 + Icges(3,6) * t525) * t519;
t479 = t494 * t524 - t543;
t478 = -t494 * t520 - t524 * t553;
t467 = rSges(3,1) * t498 - rSges(3,2) * t497 + rSges(3,3) * t554;
t466 = rSges(3,1) * t496 - rSges(3,2) * t495 - rSges(3,3) * t552;
t464 = Icges(3,1) * t498 - Icges(3,4) * t497 + Icges(3,5) * t554;
t463 = Icges(3,1) * t496 - Icges(3,4) * t495 - Icges(3,5) * t552;
t462 = Icges(3,4) * t498 - Icges(3,2) * t497 + Icges(3,6) * t554;
t461 = Icges(3,4) * t496 - Icges(3,2) * t495 - Icges(3,6) * t552;
t458 = rSges(4,1) * t494 - rSges(4,2) * t493 - rSges(4,3) * t553;
t457 = Icges(4,1) * t494 - Icges(4,4) * t493 - Icges(4,5) * t553;
t456 = Icges(4,4) * t494 - Icges(4,2) * t493 - Icges(4,6) * t553;
t455 = Icges(4,5) * t494 - Icges(4,6) * t493 - Icges(4,3) * t553;
t454 = qJD(5) * t493 + t473;
t453 = t483 * t524 + t555;
t452 = -t483 * t520 + t497 * t524;
t451 = t481 * t524 + t556;
t450 = -t481 * t520 + t495 * t524;
t439 = rSges(4,1) * t483 - rSges(4,2) * t482 + rSges(4,3) * t497;
t438 = rSges(4,1) * t481 - rSges(4,2) * t480 + rSges(4,3) * t495;
t437 = Icges(4,1) * t483 - Icges(4,4) * t482 + Icges(4,5) * t497;
t436 = Icges(4,1) * t481 - Icges(4,4) * t480 + Icges(4,5) * t495;
t435 = Icges(4,4) * t483 - Icges(4,2) * t482 + Icges(4,6) * t497;
t434 = Icges(4,4) * t481 - Icges(4,2) * t480 + Icges(4,6) * t495;
t433 = Icges(4,5) * t483 - Icges(4,6) * t482 + Icges(4,3) * t497;
t432 = Icges(4,5) * t481 - Icges(4,6) * t480 + Icges(4,3) * t495;
t431 = rSges(5,1) * t479 + rSges(5,2) * t478 + rSges(5,3) * t493;
t430 = Icges(5,1) * t479 + Icges(5,4) * t478 + Icges(5,5) * t493;
t429 = Icges(5,4) * t479 + Icges(5,2) * t478 + Icges(5,6) * t493;
t428 = Icges(5,5) * t479 + Icges(5,6) * t478 + Icges(5,3) * t493;
t426 = rSges(6,1) * t475 + rSges(6,2) * t474 + rSges(6,3) * t493;
t424 = qJD(5) * t480 + t445;
t423 = qJD(5) * t482 + t444;
t415 = t467 * t510 - t489 * t509 + t501;
t414 = -t466 * t510 - t489 * t541 - t545;
t413 = (t466 * t523 + t467 * t526) * t544;
t411 = rSges(5,1) * t453 + rSges(5,2) * t452 + rSges(5,3) * t482;
t410 = rSges(5,1) * t451 + rSges(5,2) * t450 + rSges(5,3) * t480;
t409 = Icges(5,1) * t453 + Icges(5,4) * t452 + Icges(5,5) * t482;
t408 = Icges(5,1) * t451 + Icges(5,4) * t450 + Icges(5,5) * t480;
t407 = Icges(5,4) * t453 + Icges(5,2) * t452 + Icges(5,6) * t482;
t406 = Icges(5,4) * t451 + Icges(5,2) * t450 + Icges(5,6) * t480;
t405 = Icges(5,5) * t453 + Icges(5,6) * t452 + Icges(5,3) * t482;
t404 = Icges(5,5) * t451 + Icges(5,6) * t450 + Icges(5,3) * t480;
t402 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t482;
t400 = rSges(6,1) * t447 + rSges(6,2) * t446 + rSges(6,3) * t480;
t380 = t439 * t500 - t458 * t484 + t535;
t379 = -t438 * t500 + t458 * t485 + t533;
t378 = t438 * t484 - t439 * t485 + t548;
t377 = t411 * t473 - t431 * t444 + t532;
t376 = -t410 * t473 + t431 * t445 + t531;
t375 = t410 * t444 - t411 * t445 + t536;
t374 = t402 * t454 - t423 * t426 + t530;
t373 = -t400 * t454 + t424 * t426 + t529;
t372 = t400 * t423 - t402 * t424 + t534;
t371 = qJD(6) * t480 - t423 * t549 + t454 * t550 + t530;
t370 = qJD(6) * t482 + t424 * t549 - t454 * t551 + t529;
t369 = qJD(6) * t493 + t423 * t551 - t424 * t550 + t534;
t1 = m(3) * (t413 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + m(4) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(5) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(6) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(7) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + t510 * ((t557 * t460 + (t462 * t525 + t464 * t522) * t519) * t509 - (t557 * t459 + (t461 * t525 + t463 * t522) * t519) * t541 + (t557 * t486 + (t487 * t525 + t488 * t522) * t519) * t510) / 0.2e1 + t484 * ((t497 * t433 - t482 * t435 + t483 * t437) * t484 + (t432 * t497 - t434 * t482 + t436 * t483) * t485 + (t455 * t497 - t456 * t482 + t457 * t483) * t500) / 0.2e1 + t485 * ((t433 * t495 - t435 * t480 + t437 * t481) * t484 + (t495 * t432 - t480 * t434 + t481 * t436) * t485 + (t455 * t495 - t456 * t480 + t457 * t481) * t500) / 0.2e1 + t500 * ((-t433 * t553 - t435 * t493 + t437 * t494) * t484 + (-t432 * t553 - t434 * t493 + t436 * t494) * t485 + (-t455 * t553 - t493 * t456 + t494 * t457) * t500) / 0.2e1 + t444 * ((t482 * t405 + t452 * t407 + t453 * t409) * t444 + (t404 * t482 + t406 * t452 + t408 * t453) * t445 + (t428 * t482 + t429 * t452 + t430 * t453) * t473) / 0.2e1 + t445 * ((t405 * t480 + t407 * t450 + t409 * t451) * t444 + (t480 * t404 + t450 * t406 + t451 * t408) * t445 + (t428 * t480 + t429 * t450 + t430 * t451) * t473) / 0.2e1 + t473 * ((t405 * t493 + t407 * t478 + t409 * t479) * t444 + (t404 * t493 + t406 * t478 + t408 * t479) * t445 + (t493 * t428 + t478 * t429 + t479 * t430) * t473) / 0.2e1 - ((-t486 * t552 - t487 * t495 + t488 * t496) * t510 + ((-t462 * t495 + t464 * t496) * t523 + (t495 * t461 - t496 * t463 + t575) * t526) * t544) * t541 / 0.2e1 + ((t486 * t554 - t487 * t497 + t488 * t498) * t510 + (-(-t461 * t497 + t463 * t498) * t526 + (-t497 * t462 + t498 * t464 - t575) * t523) * t544) * t509 / 0.2e1 + ((t448 * t567 + t449 * t566 + t482 * t568) * t454 + (t572 * t448 + t570 * t449 + t574 * t482) * t424 + (t571 * t448 + t569 * t449 + t573 * t482) * t423) * t423 / 0.2e1 + ((t446 * t567 + t447 * t566 + t480 * t568) * t454 + (t572 * t446 + t570 * t447 + t574 * t480) * t424 + (t571 * t446 + t569 * t447 + t573 * t480) * t423) * t424 / 0.2e1 + ((t474 * t567 + t475 * t566 + t493 * t568) * t454 + (t572 * t474 + t570 * t475 + t574 * t493) * t424 + (t571 * t474 + t569 * t475 + t573 * t493) * t423) * t454 / 0.2e1 + (m(2) * (t505 ^ 2 + t506 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
