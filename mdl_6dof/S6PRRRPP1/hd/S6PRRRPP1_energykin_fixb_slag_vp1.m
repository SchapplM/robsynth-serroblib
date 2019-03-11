% Calculate kinetic energy for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:49
% EndTime: 2019-03-08 22:42:51
% DurationCPUTime: 2.58s
% Computational Cost: add. (2941->318), mult. (6508->470), div. (0->0), fcn. (8021->12), ass. (0->147)
t581 = Icges(6,1) + Icges(7,1);
t580 = -Icges(6,4) + Icges(7,5);
t579 = Icges(7,4) + Icges(6,5);
t578 = Icges(6,2) + Icges(7,3);
t577 = -Icges(6,6) + Icges(7,6);
t576 = Icges(7,2) + Icges(5,3) + Icges(6,3);
t575 = rSges(7,1) + pkin(5);
t574 = rSges(7,3) + qJ(6);
t521 = sin(pkin(10));
t523 = cos(pkin(10));
t530 = cos(qJ(2));
t524 = cos(pkin(6));
t528 = sin(qJ(2));
t551 = t524 * t528;
t507 = t521 * t530 + t523 * t551;
t522 = sin(pkin(6));
t527 = sin(qJ(3));
t553 = t522 * t527;
t560 = cos(qJ(3));
t490 = t507 * t560 - t523 * t553;
t550 = t524 * t530;
t506 = t521 * t528 - t523 * t550;
t544 = qJ(4) + pkin(11);
t520 = sin(t544);
t539 = cos(t544);
t458 = t490 * t520 - t506 * t539;
t459 = t490 * t539 + t506 * t520;
t542 = t522 * t560;
t489 = t507 * t527 + t523 * t542;
t573 = t458 * t578 + t459 * t580 + t489 * t577;
t509 = -t521 * t551 + t523 * t530;
t492 = t509 * t560 + t521 * t553;
t508 = t521 * t550 + t523 * t528;
t460 = t492 * t520 - t508 * t539;
t461 = t492 * t539 + t508 * t520;
t491 = t509 * t527 - t521 * t542;
t572 = t460 * t578 + t461 * t580 + t491 * t577;
t571 = t580 * t458 + t459 * t581 + t579 * t489;
t570 = t580 * t460 + t461 * t581 + t579 * t491;
t511 = t524 * t527 + t528 * t542;
t552 = t522 * t530;
t484 = t511 * t520 + t539 * t552;
t485 = t511 * t539 - t520 * t552;
t510 = -t524 * t560 + t528 * t553;
t569 = t484 * t578 + t485 * t580 + t510 * t577;
t568 = t580 * t484 + t485 * t581 + t579 * t510;
t526 = sin(qJ(4));
t529 = cos(qJ(4));
t462 = -t490 * t526 + t506 * t529;
t557 = t506 * t526;
t463 = t490 * t529 + t557;
t567 = Icges(5,5) * t463 + Icges(5,6) * t462 + t458 * t577 + t459 * t579 + t489 * t576;
t464 = -t492 * t526 + t508 * t529;
t556 = t508 * t526;
t465 = t492 * t529 + t556;
t566 = Icges(5,5) * t465 + Icges(5,6) * t464 + t460 * t577 + t461 * t579 + t491 * t576;
t493 = -t511 * t526 - t529 * t552;
t543 = t526 * t552;
t494 = t511 * t529 - t543;
t565 = Icges(5,5) * t494 + Icges(5,6) * t493 + t484 * t577 + t485 * t579 + t510 * t576;
t564 = qJD(2) ^ 2;
t559 = pkin(4) * t529;
t555 = t521 * t522;
t554 = t522 * t523;
t549 = rSges(7,2) * t489 + t574 * t458 + t575 * t459;
t548 = rSges(7,2) * t491 + t574 * t460 + t575 * t461;
t547 = rSges(7,2) * t510 + t574 * t484 + t575 * t485;
t546 = qJD(2) * t522;
t516 = t521 * t546;
t495 = qJD(3) * t508 + t516;
t519 = qJD(2) * t524;
t541 = t523 * t546;
t481 = pkin(2) * t507 + pkin(8) * t506;
t482 = pkin(2) * t509 + pkin(8) * t508;
t540 = t481 * t516 + t482 * t541 + qJD(1);
t496 = qJD(3) * t506 - t541;
t513 = -qJD(3) * t552 + t519;
t512 = (pkin(2) * t528 - pkin(8) * t530) * t522;
t538 = t482 * t519 - t512 * t516;
t454 = pkin(3) * t490 + pkin(9) * t489;
t455 = pkin(3) * t492 + pkin(9) * t491;
t537 = t495 * t454 - t455 * t496 + t540;
t536 = (-t481 * t524 - t512 * t554) * qJD(2);
t397 = pkin(4) * t557 + qJ(5) * t489 + t490 * t559;
t456 = qJD(4) * t491 + t495;
t535 = qJD(5) * t510 + t456 * t397 + t537;
t483 = t511 * pkin(3) + t510 * pkin(9);
t534 = t513 * t455 - t483 * t495 + t538;
t398 = pkin(4) * t556 + qJ(5) * t491 + t492 * t559;
t488 = qJD(4) * t510 + t513;
t533 = qJD(5) * t489 + t488 * t398 + t534;
t532 = -t454 * t513 + t496 * t483 + t536;
t436 = -pkin(4) * t543 + qJ(5) * t510 + t511 * t559;
t457 = qJD(4) * t489 + t496;
t531 = qJD(5) * t491 + t457 * t436 + t532;
t502 = t524 * rSges(3,3) + (rSges(3,1) * t528 + rSges(3,2) * t530) * t522;
t501 = Icges(3,5) * t524 + (Icges(3,1) * t528 + Icges(3,4) * t530) * t522;
t500 = Icges(3,6) * t524 + (Icges(3,4) * t528 + Icges(3,2) * t530) * t522;
t499 = Icges(3,3) * t524 + (Icges(3,5) * t528 + Icges(3,6) * t530) * t522;
t479 = t511 * rSges(4,1) - t510 * rSges(4,2) - rSges(4,3) * t552;
t478 = Icges(4,1) * t511 - Icges(4,4) * t510 - Icges(4,5) * t552;
t477 = Icges(4,4) * t511 - Icges(4,2) * t510 - Icges(4,6) * t552;
t476 = Icges(4,5) * t511 - Icges(4,6) * t510 - Icges(4,3) * t552;
t473 = rSges(3,1) * t509 - rSges(3,2) * t508 + rSges(3,3) * t555;
t472 = rSges(3,1) * t507 - rSges(3,2) * t506 - rSges(3,3) * t554;
t471 = Icges(3,1) * t509 - Icges(3,4) * t508 + Icges(3,5) * t555;
t470 = Icges(3,1) * t507 - Icges(3,4) * t506 - Icges(3,5) * t554;
t469 = Icges(3,4) * t509 - Icges(3,2) * t508 + Icges(3,6) * t555;
t468 = Icges(3,4) * t507 - Icges(3,2) * t506 - Icges(3,6) * t554;
t467 = Icges(3,5) * t509 - Icges(3,6) * t508 + Icges(3,3) * t555;
t466 = Icges(3,5) * t507 - Icges(3,6) * t506 - Icges(3,3) * t554;
t450 = (-t472 * t524 - t502 * t554) * qJD(2);
t449 = (t473 * t524 - t502 * t555) * qJD(2);
t448 = rSges(5,1) * t494 + rSges(5,2) * t493 + rSges(5,3) * t510;
t447 = Icges(5,1) * t494 + Icges(5,4) * t493 + Icges(5,5) * t510;
t446 = Icges(5,4) * t494 + Icges(5,2) * t493 + Icges(5,6) * t510;
t444 = rSges(4,1) * t492 - rSges(4,2) * t491 + rSges(4,3) * t508;
t443 = rSges(4,1) * t490 - rSges(4,2) * t489 + rSges(4,3) * t506;
t442 = Icges(4,1) * t492 - Icges(4,4) * t491 + Icges(4,5) * t508;
t441 = Icges(4,1) * t490 - Icges(4,4) * t489 + Icges(4,5) * t506;
t440 = Icges(4,4) * t492 - Icges(4,2) * t491 + Icges(4,6) * t508;
t439 = Icges(4,4) * t490 - Icges(4,2) * t489 + Icges(4,6) * t506;
t438 = Icges(4,5) * t492 - Icges(4,6) * t491 + Icges(4,3) * t508;
t437 = Icges(4,5) * t490 - Icges(4,6) * t489 + Icges(4,3) * t506;
t435 = rSges(6,1) * t485 - rSges(6,2) * t484 + rSges(6,3) * t510;
t426 = qJD(1) + (t472 * t521 + t473 * t523) * t546;
t423 = rSges(5,1) * t465 + rSges(5,2) * t464 + rSges(5,3) * t491;
t422 = rSges(5,1) * t463 + rSges(5,2) * t462 + rSges(5,3) * t489;
t421 = Icges(5,1) * t465 + Icges(5,4) * t464 + Icges(5,5) * t491;
t420 = Icges(5,1) * t463 + Icges(5,4) * t462 + Icges(5,5) * t489;
t419 = Icges(5,4) * t465 + Icges(5,2) * t464 + Icges(5,6) * t491;
t418 = Icges(5,4) * t463 + Icges(5,2) * t462 + Icges(5,6) * t489;
t414 = rSges(6,1) * t461 - rSges(6,2) * t460 + rSges(6,3) * t491;
t412 = rSges(6,1) * t459 - rSges(6,2) * t458 + rSges(6,3) * t489;
t394 = -t443 * t513 + t479 * t496 + t536;
t393 = t444 * t513 - t479 * t495 + t538;
t392 = t443 * t495 - t444 * t496 + t540;
t391 = -t422 * t488 + t448 * t457 + t532;
t390 = t423 * t488 - t448 * t456 + t534;
t389 = t422 * t456 - t423 * t457 + t537;
t388 = t435 * t457 + (-t397 - t412) * t488 + t531;
t387 = t414 * t488 + (-t435 - t436) * t456 + t533;
t386 = t412 * t456 + (-t398 - t414) * t457 + t535;
t385 = qJD(6) * t460 + t547 * t457 + (-t397 - t549) * t488 + t531;
t384 = qJD(6) * t458 + t548 * t488 + (-t436 - t547) * t456 + t533;
t383 = qJD(6) * t484 + t549 * t456 + (-t398 - t548) * t457 + t535;
t1 = m(7) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(5) * (t389 ^ 2 + t390 ^ 2 + t391 ^ 2) / 0.2e1 + m(4) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + m(3) * (t426 ^ 2 + t449 ^ 2 + t450 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(6) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + t495 * ((t438 * t508 - t440 * t491 + t442 * t492) * t495 + (t437 * t508 - t439 * t491 + t441 * t492) * t496 + (t476 * t508 - t477 * t491 + t478 * t492) * t513) / 0.2e1 + t496 * ((t438 * t506 - t440 * t489 + t442 * t490) * t495 + (t437 * t506 - t439 * t489 + t441 * t490) * t496 + (t476 * t506 - t477 * t489 + t478 * t490) * t513) / 0.2e1 + t513 * ((-t438 * t552 - t510 * t440 + t511 * t442) * t495 + (-t437 * t552 - t510 * t439 + t511 * t441) * t496 + (-t476 * t552 - t510 * t477 + t511 * t478) * t513) / 0.2e1 - t564 * ((-t467 * t554 - t469 * t506 + t471 * t507) * t555 - (-t466 * t554 - t468 * t506 + t470 * t507) * t554 + (-t499 * t554 - t500 * t506 + t501 * t507) * t524) * t554 / 0.2e1 + (t524 * (t524 ^ 2 * t499 + (((t469 * t530 + t471 * t528) * t521 - (t468 * t530 + t470 * t528) * t523) * t522 + (-t466 * t523 + t467 * t521 + t500 * t530 + t501 * t528) * t524) * t522) + ((t467 * t555 - t469 * t508 + t471 * t509) * t555 - (t466 * t555 - t468 * t508 + t470 * t509) * t554 + (t499 * t555 - t500 * t508 + t501 * t509) * t524) * t555) * t564 / 0.2e1 + ((t446 * t464 + t447 * t465 + t460 * t569 + t461 * t568 + t491 * t565) * t488 + (t418 * t464 + t420 * t465 + t460 * t573 + t571 * t461 + t567 * t491) * t457 + (t419 * t464 + t421 * t465 + t460 * t572 + t461 * t570 + t491 * t566) * t456) * t456 / 0.2e1 + ((t446 * t462 + t447 * t463 + t458 * t569 + t459 * t568 + t489 * t565) * t488 + (t418 * t462 + t420 * t463 + t458 * t573 + t571 * t459 + t567 * t489) * t457 + (t419 * t462 + t421 * t463 + t458 * t572 + t459 * t570 + t489 * t566) * t456) * t457 / 0.2e1 + ((t446 * t493 + t447 * t494 + t484 * t569 + t485 * t568 + t510 * t565) * t488 + (t418 * t493 + t420 * t494 + t484 * t573 + t571 * t485 + t567 * t510) * t457 + (t419 * t493 + t421 * t494 + t484 * t572 + t485 * t570 + t510 * t566) * t456) * t488 / 0.2e1;
T  = t1;
