% Calculate kinetic energy for
% S6RRPRPR12
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:14
% EndTime: 2019-03-09 11:17:17
% DurationCPUTime: 3.15s
% Computational Cost: add. (2313->335), mult. (4517->493), div. (0->0), fcn. (5282->12), ass. (0->155)
t585 = Icges(3,1) + Icges(4,2);
t584 = Icges(3,4) + Icges(4,6);
t583 = Icges(3,5) - Icges(4,4);
t582 = Icges(3,2) + Icges(4,3);
t581 = Icges(3,6) - Icges(4,5);
t580 = Icges(3,3) + Icges(4,1);
t579 = Icges(5,3) + Icges(6,3);
t515 = sin(pkin(6));
t516 = cos(pkin(6));
t524 = cos(qJ(2));
t525 = cos(qJ(1));
t549 = t524 * t525;
t520 = sin(qJ(2));
t521 = sin(qJ(1));
t552 = t520 * t521;
t495 = -t516 * t549 + t552;
t550 = t521 * t524;
t551 = t520 * t525;
t496 = t516 * t551 + t550;
t497 = t516 * t550 + t551;
t498 = -t516 * t552 + t549;
t554 = t515 * t525;
t556 = t515 * t521;
t567 = (t495 * t581 - t496 * t583 + t554 * t580) * t525 + (-t497 * t581 + t498 * t583 + t556 * t580) * t521;
t578 = t515 * t567;
t544 = qJ(4) + pkin(11);
t514 = sin(t544);
t541 = cos(t544);
t461 = -t497 * t541 + t514 * t556;
t535 = t515 * t541;
t462 = t497 * t514 + t521 * t535;
t519 = sin(qJ(4));
t523 = cos(qJ(4));
t468 = t497 * t523 - t519 * t556;
t558 = t497 * t519;
t469 = t523 * t556 + t558;
t577 = Icges(5,5) * t469 + Icges(6,5) * t462 + Icges(5,6) * t468 - Icges(6,6) * t461 + t498 * t579;
t463 = t495 * t541 + t514 * t554;
t464 = t495 * t514 - t525 * t535;
t470 = t495 * t523 + t519 * t554;
t559 = t495 * t519;
t471 = -t523 * t554 + t559;
t576 = Icges(5,5) * t471 + Icges(6,5) * t464 + Icges(5,6) * t470 + Icges(6,6) * t463 + t496 * t579;
t483 = t516 * t514 + t524 * t535;
t555 = t515 * t524;
t484 = -t514 * t555 + t516 * t541;
t493 = -t516 * t519 - t523 * t555;
t553 = t519 * t524;
t494 = -t515 * t553 + t516 * t523;
t557 = t515 * t520;
t575 = Icges(5,5) * t494 + Icges(6,5) * t484 + Icges(5,6) * t493 - Icges(6,6) * t483 + t557 * t579;
t574 = t497 * t582 - t498 * t584 - t556 * t581;
t573 = t495 * t582 - t496 * t584 + t554 * t581;
t572 = -t584 * t497 + t498 * t585 + t583 * t556;
t571 = t584 * t495 - t496 * t585 + t583 * t554;
t570 = t581 * t516 + (t520 * t584 + t524 * t582) * t515;
t569 = t583 * t516 + (t520 * t585 + t584 * t524) * t515;
t568 = t580 * t516 + (t520 * t583 + t524 * t581) * t515;
t561 = pkin(4) * t523;
t457 = pkin(2) * t496 + qJ(3) * t495;
t458 = pkin(2) * t498 + qJ(3) * t497;
t546 = qJD(2) * t515;
t510 = t521 * t546;
t542 = t525 * t546;
t548 = t457 * t510 + t458 * t542;
t472 = qJD(4) * t498 + t510;
t547 = qJD(1) * (pkin(1) * t521 - pkin(8) * t554);
t545 = qJD(3) * t524;
t511 = qJD(2) * t516 + qJD(1);
t501 = qJD(1) * (pkin(1) * t525 + pkin(8) * t556);
t543 = qJD(3) * t495 + t511 * t458 + t501;
t500 = qJD(4) * t557 + t511;
t540 = qJD(3) * t497 - t547;
t499 = (pkin(2) * t520 - qJ(3) * t524) * t515;
t537 = (-rSges(4,1) * t516 - (-rSges(4,2) * t520 - rSges(4,3) * t524) * t515 - t499) * t546;
t536 = (-pkin(3) * t516 - pkin(9) * t557 - t499) * t546;
t473 = qJD(4) * t496 - t542;
t474 = pkin(3) * t556 + pkin(9) * t498;
t475 = -pkin(3) * t554 + pkin(9) * t496;
t532 = t474 * t542 + t475 * t510 - t515 * t545 + t548;
t417 = pkin(4) * t559 + qJ(5) * t496 - t554 * t561;
t531 = qJD(5) * t557 + t472 * t417 + t532;
t530 = t511 * t474 + t521 * t536 + t543;
t416 = pkin(4) * t558 + qJ(5) * t498 + t556 * t561;
t529 = qJD(5) * t496 + t500 * t416 + t530;
t528 = (-t457 - t475) * t511 + t525 * t536 + t540;
t455 = t561 * t516 + (-pkin(4) * t553 + qJ(5) * t520) * t515;
t527 = qJD(5) * t498 + t473 * t455 + t528;
t522 = cos(qJ(6));
t518 = sin(qJ(6));
t505 = rSges(2,1) * t525 - rSges(2,2) * t521;
t504 = rSges(2,1) * t521 + rSges(2,2) * t525;
t485 = t516 * rSges(3,3) + (rSges(3,1) * t520 + rSges(3,2) * t524) * t515;
t460 = t484 * t522 + t518 * t557;
t459 = -t484 * t518 + t522 * t557;
t456 = qJD(6) * t483 + t500;
t454 = rSges(3,1) * t498 - rSges(3,2) * t497 + rSges(3,3) * t556;
t453 = rSges(3,1) * t496 - rSges(3,2) * t495 - rSges(3,3) * t554;
t452 = -rSges(4,1) * t554 - rSges(4,2) * t496 + rSges(4,3) * t495;
t451 = rSges(4,1) * t556 - rSges(4,2) * t498 + rSges(4,3) * t497;
t435 = rSges(5,1) * t494 + rSges(5,2) * t493 + rSges(5,3) * t557;
t434 = Icges(5,1) * t494 + Icges(5,4) * t493 + Icges(5,5) * t557;
t433 = Icges(5,4) * t494 + Icges(5,2) * t493 + Icges(5,6) * t557;
t431 = pkin(5) * t484 + pkin(10) * t483;
t430 = rSges(6,1) * t484 - rSges(6,2) * t483 + rSges(6,3) * t557;
t429 = Icges(6,1) * t484 - Icges(6,4) * t483 + Icges(6,5) * t557;
t428 = Icges(6,4) * t484 - Icges(6,2) * t483 + Icges(6,6) * t557;
t426 = t464 * t522 + t496 * t518;
t425 = -t464 * t518 + t496 * t522;
t424 = t462 * t522 + t498 * t518;
t423 = -t462 * t518 + t498 * t522;
t422 = -qJD(6) * t463 + t473;
t421 = qJD(6) * t461 + t472;
t419 = pkin(5) * t464 - pkin(10) * t463;
t418 = pkin(5) * t462 + pkin(10) * t461;
t415 = rSges(5,1) * t471 + rSges(5,2) * t470 + rSges(5,3) * t496;
t414 = rSges(5,1) * t469 + rSges(5,2) * t468 + rSges(5,3) * t498;
t413 = Icges(5,1) * t471 + Icges(5,4) * t470 + Icges(5,5) * t496;
t412 = Icges(5,1) * t469 + Icges(5,4) * t468 + Icges(5,5) * t498;
t411 = Icges(5,4) * t471 + Icges(5,2) * t470 + Icges(5,6) * t496;
t410 = Icges(5,4) * t469 + Icges(5,2) * t468 + Icges(5,6) * t498;
t407 = rSges(6,1) * t464 + rSges(6,2) * t463 + rSges(6,3) * t496;
t406 = rSges(6,1) * t462 - rSges(6,2) * t461 + rSges(6,3) * t498;
t405 = Icges(6,1) * t464 + Icges(6,4) * t463 + Icges(6,5) * t496;
t404 = Icges(6,1) * t462 - Icges(6,4) * t461 + Icges(6,5) * t498;
t403 = Icges(6,4) * t464 + Icges(6,2) * t463 + Icges(6,6) * t496;
t402 = Icges(6,4) * t462 - Icges(6,2) * t461 + Icges(6,6) * t498;
t398 = t454 * t511 - t485 * t510 + t501;
t397 = -t453 * t511 - t485 * t542 - t547;
t396 = rSges(7,1) * t460 + rSges(7,2) * t459 + rSges(7,3) * t483;
t395 = Icges(7,1) * t460 + Icges(7,4) * t459 + Icges(7,5) * t483;
t394 = Icges(7,4) * t460 + Icges(7,2) * t459 + Icges(7,6) * t483;
t393 = Icges(7,5) * t460 + Icges(7,6) * t459 + Icges(7,3) * t483;
t392 = (t453 * t521 + t454 * t525) * t546;
t390 = rSges(7,1) * t426 + rSges(7,2) * t425 - rSges(7,3) * t463;
t389 = rSges(7,1) * t424 + rSges(7,2) * t423 + rSges(7,3) * t461;
t388 = Icges(7,1) * t426 + Icges(7,4) * t425 - Icges(7,5) * t463;
t387 = Icges(7,1) * t424 + Icges(7,4) * t423 + Icges(7,5) * t461;
t386 = Icges(7,4) * t426 + Icges(7,2) * t425 - Icges(7,6) * t463;
t385 = Icges(7,4) * t424 + Icges(7,2) * t423 + Icges(7,6) * t461;
t384 = Icges(7,5) * t426 + Icges(7,6) * t425 - Icges(7,3) * t463;
t383 = Icges(7,5) * t424 + Icges(7,6) * t423 + Icges(7,3) * t461;
t382 = t451 * t511 + t521 * t537 + t543;
t381 = (-t452 - t457) * t511 + t525 * t537 + t540;
t380 = (-t545 + (t451 * t525 + t452 * t521) * qJD(2)) * t515 + t548;
t379 = t414 * t500 - t435 * t472 + t530;
t378 = -t415 * t500 + t435 * t473 + t528;
t377 = -t414 * t473 + t415 * t472 + t532;
t376 = t406 * t500 + (-t430 - t455) * t472 + t529;
t375 = t430 * t473 + (-t407 - t417) * t500 + t527;
t374 = t407 * t472 + (-t406 - t416) * t473 + t531;
t373 = t389 * t456 - t396 * t421 + t418 * t500 + (-t431 - t455) * t472 + t529;
t372 = -t390 * t456 + t396 * t422 + t431 * t473 + (-t417 - t419) * t500 + t527;
t371 = -t389 * t422 + t390 * t421 + t419 * t472 + (-t416 - t418) * t473 + t531;
t1 = m(4) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(6) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(5) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + t421 * ((t461 * t383 + t423 * t385 + t424 * t387) * t421 + (t384 * t461 + t386 * t423 + t388 * t424) * t422 + (t393 * t461 + t394 * t423 + t395 * t424) * t456) / 0.2e1 + t456 * ((t383 * t483 + t385 * t459 + t387 * t460) * t421 + (t384 * t483 + t386 * t459 + t388 * t460) * t422 + (t483 * t393 + t459 * t394 + t460 * t395) * t456) / 0.2e1 + t422 * ((-t383 * t463 + t385 * t425 + t387 * t426) * t421 + (-t463 * t384 + t425 * t386 + t426 * t388) * t422 + (-t393 * t463 + t394 * t425 + t395 * t426) * t456) / 0.2e1 + m(3) * (t392 ^ 2 + t397 ^ 2 + t398 ^ 2) / 0.2e1 + m(7) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + ((-t428 * t461 + t429 * t462 + t433 * t468 + t434 * t469 + t498 * t575) * t500 + (-t403 * t461 + t405 * t462 + t411 * t468 + t413 * t469 + t498 * t576) * t473 + (-t461 * t402 + t462 * t404 + t468 * t410 + t469 * t412 + t577 * t498) * t472) * t472 / 0.2e1 + ((t428 * t463 + t429 * t464 + t433 * t470 + t434 * t471 + t496 * t575) * t500 + (t463 * t403 + t464 * t405 + t470 * t411 + t471 * t413 + t576 * t496) * t473 + (t402 * t463 + t404 * t464 + t410 * t470 + t412 * t471 + t496 * t577) * t472) * t473 / 0.2e1 + ((-t483 * t428 + t484 * t429 + t493 * t433 + t494 * t434 + t575 * t557) * t500 + (-t403 * t483 + t405 * t484 + t411 * t493 + t413 * t494 + t557 * t576) * t473 + (-t402 * t483 + t404 * t484 + t410 * t493 + t412 * t494 + t557 * t577) * t472) * t500 / 0.2e1 + ((t567 * t516 + ((t520 * t571 + t524 * t573) * t525 + (t520 * t572 - t524 * t574) * t521) * t515) * t546 + (t568 * t516 + (t520 * t569 + t524 * t570) * t515) * t511) * t511 / 0.2e1 + (Icges(2,3) + m(2) * (t504 ^ 2 + t505 ^ 2)) * qJD(1) ^ 2 / 0.2e1 - (((-t495 * t573 + t496 * t571 - t578) * t525 + (t495 * t574 + t496 * t572) * t521) * t546 + (-t495 * t570 + t496 * t569 - t554 * t568) * t511) * t542 / 0.2e1 + (((-t497 * t573 + t498 * t571) * t525 + (t497 * t574 + t498 * t572 + t578) * t521) * t546 + (-t497 * t570 + t498 * t569 + t556 * t568) * t511) * t510 / 0.2e1;
T  = t1;
