% Calculate kinetic energy for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:40
% EndTime: 2019-03-09 11:11:44
% DurationCPUTime: 3.96s
% Computational Cost: add. (1405->307), mult. (2173->469), div. (0->0), fcn. (2130->10), ass. (0->156)
t590 = Icges(3,4) + Icges(4,6);
t589 = Icges(3,1) + Icges(4,2);
t588 = -Icges(3,2) - Icges(4,3);
t503 = cos(qJ(2));
t587 = t590 * t503;
t500 = sin(qJ(2));
t586 = t590 * t500;
t585 = -Icges(4,4) + Icges(3,5);
t584 = Icges(4,5) - Icges(3,6);
t583 = t588 * t500 + t587;
t582 = -t589 * t503 + t586;
t581 = Icges(4,1) + Icges(3,3);
t580 = Icges(6,3) + Icges(5,3);
t501 = sin(qJ(1));
t504 = cos(qJ(1));
t579 = t583 * t501 + t584 * t504;
t578 = -t584 * t501 + t583 * t504;
t577 = t582 * t501 + t585 * t504;
t576 = t585 * t501 - t582 * t504;
t575 = t588 * t503 - t586;
t574 = t589 * t500 + t587;
t573 = t584 * t500 + t585 * t503;
t497 = qJ(4) + pkin(10);
t489 = sin(t497);
t490 = cos(t497);
t548 = t500 * t504;
t422 = -t489 * t501 + t490 * t548;
t423 = t489 * t548 + t490 * t501;
t499 = sin(qJ(4));
t502 = cos(qJ(4));
t545 = t502 * t504;
t452 = -t499 * t501 + t500 * t545;
t547 = t501 * t502;
t453 = t499 * t548 + t547;
t544 = t503 * t504;
t572 = Icges(5,5) * t453 + Icges(6,5) * t423 + Icges(5,6) * t452 + Icges(6,6) * t422 + t544 * t580;
t549 = t500 * t501;
t424 = t489 * t504 + t490 * t549;
t425 = t489 * t549 - t490 * t504;
t454 = t499 * t504 + t500 * t547;
t455 = t499 * t549 - t545;
t546 = t501 * t503;
t571 = Icges(5,5) * t455 + Icges(6,5) * t425 + Icges(5,6) * t454 + Icges(6,6) * t424 + t546 * t580;
t570 = (-Icges(5,5) * t499 - Icges(6,5) * t489 - Icges(5,6) * t502 - Icges(6,6) * t490) * t503 + t580 * t500;
t569 = t501 * t581 + t573 * t504;
t568 = t573 * t501 - t504 * t581;
t567 = t585 * t500 - t584 * t503;
t566 = t500 * t575 + t503 * t574;
t565 = -t500 * t578 + t503 * t576;
t564 = t500 * t579 + t503 * t577;
t557 = pkin(4) * t499;
t555 = t502 * pkin(4);
t527 = pkin(2) * t503 + qJ(3) * t500;
t456 = t527 * t501;
t478 = pkin(1) * t501 - pkin(7) * t504;
t543 = -t456 - t478;
t542 = pkin(5) * t490;
t494 = qJD(2) * t501;
t538 = qJD(4) * t503;
t459 = t504 * t538 + t494;
t540 = qJD(2) * t504;
t539 = qJD(3) * t500;
t537 = qJD(5) * t503;
t536 = qJD(6) * t503;
t485 = qJD(4) * t500 + qJD(1);
t457 = t527 * t504;
t466 = qJD(1) * (pkin(1) * t504 + pkin(7) * t501);
t535 = qJD(1) * t457 + t501 * t539 + t466;
t532 = pkin(5) * t489;
t473 = pkin(2) * t500 - qJ(3) * t503;
t531 = qJD(2) * (rSges(4,2) * t500 + rSges(4,3) * t503 - t473);
t460 = t501 * t538 - t540;
t530 = -qJD(3) * t503 + t456 * t494 + t457 * t540;
t529 = rSges(3,1) * t503 - rSges(3,2) * t500;
t528 = -rSges(4,2) * t503 + rSges(4,3) * t500;
t526 = qJD(2) * (-pkin(8) * t500 - t473);
t462 = pkin(3) * t501 + pkin(8) * t544;
t463 = -pkin(3) * t504 + pkin(8) * t546;
t513 = t462 * t540 + t463 * t494 + t530;
t512 = qJ(5) * t503 + t500 * t557;
t405 = t501 * t512 - t504 * t555;
t511 = qJD(5) * t500 + t459 * t405 + t513;
t510 = pkin(9) * t503 + t500 * t532;
t509 = qJD(1) * t462 + t501 * t526 + t535;
t404 = t501 * t555 + t504 * t512;
t508 = t485 * t404 + t501 * t537 + t509;
t484 = t504 * t539;
t507 = t484 + (-t463 + t543) * qJD(1) + t504 * t526;
t448 = qJ(5) * t500 - t503 * t557;
t506 = t460 * t448 + t504 * t537 + t507;
t491 = qJ(6) + t497;
t487 = cos(t491);
t486 = sin(t491);
t477 = rSges(2,1) * t504 - rSges(2,2) * t501;
t476 = rSges(2,1) * t501 + rSges(2,2) * t504;
t475 = rSges(3,1) * t500 + rSges(3,2) * t503;
t464 = qJD(6) * t500 + t485;
t447 = -rSges(4,1) * t504 + t501 * t528;
t446 = rSges(4,1) * t501 + t504 * t528;
t445 = rSges(3,3) * t501 + t504 * t529;
t444 = rSges(5,3) * t500 + (-rSges(5,1) * t499 - rSges(5,2) * t502) * t503;
t443 = -rSges(3,3) * t504 + t501 * t529;
t432 = Icges(5,5) * t500 + (-Icges(5,1) * t499 - Icges(5,4) * t502) * t503;
t429 = Icges(5,6) * t500 + (-Icges(5,4) * t499 - Icges(5,2) * t502) * t503;
t421 = t501 * t536 + t460;
t420 = t504 * t536 + t459;
t419 = t486 * t549 - t487 * t504;
t418 = t486 * t504 + t487 * t549;
t417 = t486 * t548 + t487 * t501;
t416 = -t486 * t501 + t487 * t548;
t415 = rSges(6,3) * t500 + (-rSges(6,1) * t489 - rSges(6,2) * t490) * t503;
t414 = Icges(6,5) * t500 + (-Icges(6,1) * t489 - Icges(6,4) * t490) * t503;
t413 = Icges(6,6) * t500 + (-Icges(6,4) * t489 - Icges(6,2) * t490) * t503;
t411 = rSges(7,3) * t500 + (-rSges(7,1) * t486 - rSges(7,2) * t487) * t503;
t410 = Icges(7,5) * t500 + (-Icges(7,1) * t486 - Icges(7,4) * t487) * t503;
t409 = Icges(7,6) * t500 + (-Icges(7,4) * t486 - Icges(7,2) * t487) * t503;
t408 = Icges(7,3) * t500 + (-Icges(7,5) * t486 - Icges(7,6) * t487) * t503;
t406 = pkin(9) * t500 - t503 * t532;
t403 = rSges(5,1) * t455 + rSges(5,2) * t454 + rSges(5,3) * t546;
t402 = rSges(5,1) * t453 + rSges(5,2) * t452 + rSges(5,3) * t544;
t401 = Icges(5,1) * t455 + Icges(5,4) * t454 + Icges(5,5) * t546;
t400 = Icges(5,1) * t453 + Icges(5,4) * t452 + Icges(5,5) * t544;
t399 = Icges(5,4) * t455 + Icges(5,2) * t454 + Icges(5,6) * t546;
t398 = Icges(5,4) * t453 + Icges(5,2) * t452 + Icges(5,6) * t544;
t394 = qJD(1) * t445 - t475 * t494 + t466;
t393 = -t475 * t540 + (-t443 - t478) * qJD(1);
t392 = (t443 * t501 + t445 * t504) * qJD(2);
t391 = rSges(6,1) * t425 + rSges(6,2) * t424 + rSges(6,3) * t546;
t390 = rSges(6,1) * t423 + rSges(6,2) * t422 + rSges(6,3) * t544;
t389 = Icges(6,1) * t425 + Icges(6,4) * t424 + Icges(6,5) * t546;
t388 = Icges(6,1) * t423 + Icges(6,4) * t422 + Icges(6,5) * t544;
t387 = Icges(6,4) * t425 + Icges(6,2) * t424 + Icges(6,6) * t546;
t386 = Icges(6,4) * t423 + Icges(6,2) * t422 + Icges(6,6) * t544;
t382 = rSges(7,1) * t419 + rSges(7,2) * t418 + rSges(7,3) * t546;
t381 = rSges(7,1) * t417 + rSges(7,2) * t416 + rSges(7,3) * t544;
t380 = Icges(7,1) * t419 + Icges(7,4) * t418 + Icges(7,5) * t546;
t379 = Icges(7,1) * t417 + Icges(7,4) * t416 + Icges(7,5) * t544;
t378 = Icges(7,4) * t419 + Icges(7,2) * t418 + Icges(7,6) * t546;
t377 = Icges(7,4) * t417 + Icges(7,2) * t416 + Icges(7,6) * t544;
t376 = Icges(7,5) * t419 + Icges(7,6) * t418 + Icges(7,3) * t546;
t375 = Icges(7,5) * t417 + Icges(7,6) * t416 + Icges(7,3) * t544;
t374 = t501 * t510 - t504 * t542;
t373 = t501 * t542 + t504 * t510;
t372 = qJD(1) * t446 + t501 * t531 + t535;
t371 = t484 + t504 * t531 + (-t447 + t543) * qJD(1);
t370 = (t446 * t504 + t447 * t501) * qJD(2) + t530;
t369 = t402 * t485 - t444 * t459 + t509;
t368 = -t403 * t485 + t444 * t460 + t507;
t367 = -t402 * t460 + t403 * t459 + t513;
t366 = t390 * t485 + (-t415 - t448) * t459 + t508;
t365 = t415 * t460 + (-t391 - t405) * t485 + t506;
t364 = t391 * t459 + (-t390 - t404) * t460 + t511;
t363 = t373 * t485 + t381 * t464 - t411 * t420 + (-t406 - t448) * t459 + t508;
t362 = -t382 * t464 + t406 * t460 + t411 * t421 + (-t374 - t405) * t485 + t506;
t361 = t374 * t459 - t381 * t421 + t382 * t420 + (-t373 - t404) * t460 + t511;
t1 = m(7) * (t361 ^ 2 + t362 ^ 2 + t363 ^ 2) / 0.2e1 + m(6) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + m(5) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(4) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(3) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + t421 * ((t375 * t546 + t377 * t418 + t379 * t419) * t420 + (t376 * t546 + t418 * t378 + t419 * t380) * t421 + (t408 * t546 + t409 * t418 + t410 * t419) * t464) / 0.2e1 + t464 * ((t375 * t420 + t376 * t421 + t408 * t464) * t500 + ((-t377 * t487 - t379 * t486) * t420 + (-t378 * t487 - t380 * t486) * t421 + (-t409 * t487 - t410 * t486) * t464) * t503) / 0.2e1 + t420 * ((t375 * t544 + t416 * t377 + t417 * t379) * t420 + (t376 * t544 + t378 * t416 + t380 * t417) * t421 + (t408 * t544 + t409 * t416 + t410 * t417) * t464) / 0.2e1 + ((t413 * t422 + t414 * t423 + t429 * t452 + t432 * t453 + t544 * t570) * t485 + (t387 * t422 + t389 * t423 + t399 * t452 + t401 * t453 + t544 * t571) * t460 + (t422 * t386 + t423 * t388 + t452 * t398 + t453 * t400 + t544 * t572) * t459) * t459 / 0.2e1 + ((t413 * t424 + t414 * t425 + t429 * t454 + t432 * t455 + t546 * t570) * t485 + (t424 * t387 + t425 * t389 + t454 * t399 + t455 * t401 + t546 * t571) * t460 + (t386 * t424 + t388 * t425 + t398 * t454 + t400 * t455 + t546 * t572) * t459) * t460 / 0.2e1 + (((-t413 * t490 - t414 * t489 - t429 * t502 - t432 * t499) * t485 + (-t387 * t490 - t389 * t489 - t399 * t502 - t401 * t499) * t460 + (-t386 * t490 - t388 * t489 - t398 * t502 - t400 * t499) * t459) * t503 + (t459 * t572 + t571 * t460 + t570 * t485) * t500) * t485 / 0.2e1 + (Icges(2,3) + m(2) * (t476 ^ 2 + t477 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t500 * t577 - t503 * t579) * t504 + (t500 * t576 + t503 * t578) * t501) * qJD(2) + (t574 * t500 - t575 * t503) * qJD(1)) * qJD(1) / 0.2e1 + ((t569 * t501 ^ 2 + (t564 * t504 + (t565 - t568) * t501) * t504) * qJD(2) + (t501 * t567 + t504 * t566) * qJD(1)) * t494 / 0.2e1 - ((t568 * t504 ^ 2 + (t565 * t501 + (t564 - t569) * t504) * t501) * qJD(2) + (t501 * t566 - t504 * t567) * qJD(1)) * t540 / 0.2e1;
T  = t1;
