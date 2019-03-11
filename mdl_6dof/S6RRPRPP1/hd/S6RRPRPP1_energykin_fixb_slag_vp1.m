% Calculate kinetic energy for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:06
% EndTime: 2019-03-09 09:45:08
% DurationCPUTime: 2.31s
% Computational Cost: add. (1918->273), mult. (2141->398), div. (0->0), fcn. (2108->10), ass. (0->147)
t557 = Icges(6,1) + Icges(7,1);
t556 = -Icges(6,4) + Icges(7,5);
t555 = Icges(7,4) + Icges(6,5);
t554 = Icges(6,2) + Icges(7,3);
t553 = -Icges(7,6) + Icges(6,6);
t552 = Icges(3,3) + Icges(4,3);
t462 = qJ(2) + pkin(9);
t457 = sin(t462);
t459 = cos(t462);
t466 = sin(qJ(2));
t469 = cos(qJ(2));
t551 = Icges(3,5) * t469 + Icges(4,5) * t459 - Icges(3,6) * t466 - Icges(4,6) * t457;
t550 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t549 = rSges(7,1) + pkin(5);
t548 = rSges(7,3) + qJ(6);
t461 = qJ(4) + pkin(10);
t458 = cos(t461);
t470 = cos(qJ(1));
t456 = sin(t461);
t467 = sin(qJ(1));
t516 = t456 * t467;
t416 = t458 * t470 + t459 * t516;
t510 = t467 * t458;
t417 = -t456 * t470 + t459 * t510;
t515 = t457 * t467;
t547 = t554 * t416 + t556 * t417 - t553 * t515;
t513 = t459 * t470;
t418 = t456 * t513 - t510;
t419 = t458 * t513 + t516;
t514 = t457 * t470;
t546 = t554 * t418 + t556 * t419 - t553 * t514;
t545 = t556 * t416 + t557 * t417 + t555 * t515;
t544 = t556 * t418 + t557 * t419 + t555 * t514;
t543 = t553 * t459 + (t554 * t456 + t556 * t458) * t457;
t542 = -t555 * t459 + (t556 * t456 + t557 * t458) * t457;
t541 = t551 * t467 - t552 * t470;
t540 = t552 * t467 + t551 * t470;
t539 = Icges(3,5) * t466 + Icges(4,5) * t457 + Icges(3,6) * t469 + Icges(4,6) * t459;
t518 = Icges(4,4) * t457;
t438 = Icges(4,2) * t459 + t518;
t517 = Icges(4,4) * t459;
t439 = Icges(4,1) * t457 + t517;
t520 = Icges(3,4) * t466;
t444 = Icges(3,2) * t469 + t520;
t519 = Icges(3,4) * t469;
t445 = Icges(3,1) * t466 + t519;
t538 = -t438 * t457 + t439 * t459 - t444 * t466 + t445 * t469;
t486 = -Icges(4,2) * t457 + t517;
t408 = Icges(4,6) * t467 + t470 * t486;
t488 = Icges(4,1) * t459 - t518;
t410 = Icges(4,5) * t467 + t470 * t488;
t487 = -Icges(3,2) * t466 + t519;
t424 = Icges(3,6) * t467 + t470 * t487;
t489 = Icges(3,1) * t469 - t520;
t426 = Icges(3,5) * t467 + t470 * t489;
t537 = -t408 * t457 + t410 * t459 - t424 * t466 + t426 * t469;
t407 = -Icges(4,6) * t470 + t467 * t486;
t409 = -Icges(4,5) * t470 + t467 * t488;
t423 = -Icges(3,6) * t470 + t467 * t487;
t425 = -Icges(3,5) * t470 + t467 * t489;
t536 = t407 * t457 - t409 * t459 + t423 * t466 - t425 * t469;
t468 = cos(qJ(4));
t508 = t468 * t470;
t465 = sin(qJ(4));
t512 = t465 * t467;
t431 = -t459 * t512 - t508;
t509 = t467 * t468;
t511 = t465 * t470;
t432 = t459 * t509 - t511;
t535 = Icges(5,5) * t432 + Icges(5,6) * t431 - t553 * t416 + t555 * t417 - t550 * t515;
t433 = -t459 * t511 + t509;
t434 = t459 * t508 + t512;
t534 = Icges(5,5) * t434 + Icges(5,6) * t433 - t553 * t418 + t555 * t419 - t550 * t514;
t533 = t550 * t459 + (Icges(5,5) * t468 - Icges(5,6) * t465 - t553 * t456 + t555 * t458) * t457;
t526 = pkin(2) * t466;
t524 = pkin(2) * t469;
t523 = pkin(4) * t468;
t507 = rSges(7,2) * t515 + t548 * t416 + t417 * t549;
t506 = rSges(7,2) * t514 + t548 * t418 + t419 * t549;
t505 = -rSges(7,2) * t459 + (t548 * t456 + t458 * t549) * t457;
t403 = -qJ(3) * t470 + t467 * t524;
t404 = qJ(3) * t467 + t470 * t524;
t501 = qJD(2) * t470;
t502 = qJD(2) * t467;
t504 = t403 * t502 + t404 * t501;
t451 = pkin(1) * t467 - pkin(7) * t470;
t503 = -t403 - t451;
t500 = qJD(4) * t457;
t499 = qJD(5) * t457;
t495 = pkin(3) * t459 + pkin(8) * t457;
t429 = t495 * t467;
t430 = t495 * t470;
t496 = t429 * t502 + t430 * t501 + t504;
t442 = qJD(1) * (pkin(1) * t470 + pkin(7) * t467);
t494 = qJD(1) * t404 - qJD(3) * t470 + t442;
t493 = rSges(3,1) * t469 - rSges(3,2) * t466;
t492 = rSges(4,1) * t459 - rSges(4,2) * t457;
t491 = qJD(2) * (-rSges(4,1) * t457 - rSges(4,2) * t459 - t526);
t490 = qJD(2) * (-pkin(3) * t457 + pkin(8) * t459 - t526);
t477 = qJ(5) * t457 + t459 * t523;
t371 = -pkin(4) * t511 + t467 * t477;
t435 = t470 * t500 + t502;
t476 = -qJD(5) * t459 + t435 * t371 + t496;
t475 = qJD(1) * t430 + t467 * t490 + t494;
t460 = qJD(3) * t467;
t474 = t460 + (-t429 + t503) * qJD(1) + t470 * t490;
t372 = pkin(4) * t512 + t470 * t477;
t452 = -qJD(4) * t459 + qJD(1);
t473 = t452 * t372 + t467 * t499 + t475;
t387 = -qJ(5) * t459 + t457 * t523;
t436 = t467 * t500 - t501;
t472 = t436 * t387 + t470 * t499 + t474;
t450 = rSges(2,1) * t470 - rSges(2,2) * t467;
t449 = rSges(2,1) * t467 + rSges(2,2) * t470;
t448 = rSges(3,1) * t466 + rSges(3,2) * t469;
t428 = t467 * rSges(3,3) + t470 * t493;
t427 = -t470 * rSges(3,3) + t467 * t493;
t413 = t467 * rSges(4,3) + t470 * t492;
t412 = -t470 * rSges(4,3) + t467 * t492;
t402 = -t459 * rSges(5,3) + (rSges(5,1) * t468 - rSges(5,2) * t465) * t457;
t400 = -Icges(5,5) * t459 + (Icges(5,1) * t468 - Icges(5,4) * t465) * t457;
t399 = -Icges(5,6) * t459 + (Icges(5,4) * t468 - Icges(5,2) * t465) * t457;
t395 = -t459 * rSges(6,3) + (rSges(6,1) * t458 - rSges(6,2) * t456) * t457;
t386 = qJD(1) * t428 - t448 * t502 + t442;
t385 = -t448 * t501 + (-t427 - t451) * qJD(1);
t382 = (t427 * t467 + t428 * t470) * qJD(2);
t381 = rSges(5,1) * t434 + rSges(5,2) * t433 + rSges(5,3) * t514;
t380 = rSges(5,1) * t432 + rSges(5,2) * t431 + rSges(5,3) * t515;
t379 = Icges(5,1) * t434 + Icges(5,4) * t433 + Icges(5,5) * t514;
t378 = Icges(5,1) * t432 + Icges(5,4) * t431 + Icges(5,5) * t515;
t377 = Icges(5,4) * t434 + Icges(5,2) * t433 + Icges(5,6) * t514;
t376 = Icges(5,4) * t432 + Icges(5,2) * t431 + Icges(5,6) * t515;
t370 = rSges(6,1) * t419 - rSges(6,2) * t418 + rSges(6,3) * t514;
t368 = rSges(6,1) * t417 - rSges(6,2) * t416 + rSges(6,3) * t515;
t352 = qJD(1) * t413 + t467 * t491 + t494;
t351 = t460 + t470 * t491 + (-t412 + t503) * qJD(1);
t350 = (t412 * t467 + t413 * t470) * qJD(2) + t504;
t349 = t381 * t452 - t402 * t435 + t475;
t348 = -t380 * t452 + t402 * t436 + t474;
t347 = t380 * t435 - t381 * t436 + t496;
t346 = t370 * t452 + (-t387 - t395) * t435 + t473;
t345 = t395 * t436 + (-t368 - t371) * t452 + t472;
t344 = t368 * t435 + (-t370 - t372) * t436 + t476;
t343 = qJD(6) * t416 + t506 * t452 + (-t387 - t505) * t435 + t473;
t342 = qJD(6) * t418 + t505 * t436 + (-t371 - t507) * t452 + t472;
t341 = qJD(6) * t456 * t457 + t507 * t435 + (-t372 - t506) * t436 + t476;
t1 = m(3) * (t382 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(4) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(5) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(6) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(7) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t449 ^ 2 + t450 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t407 * t459 - t409 * t457 - t423 * t469 - t425 * t466) * t470 + (t408 * t459 + t410 * t457 + t424 * t469 + t426 * t466) * t467) * qJD(2) + (t459 * t438 + t457 * t439 + t469 * t444 + t466 * t445) * qJD(1)) * qJD(1) / 0.2e1 + ((t540 * t467 ^ 2 + (t536 * t470 + (t537 - t541) * t467) * t470) * qJD(2) + (t539 * t467 + t538 * t470) * qJD(1)) * t502 / 0.2e1 - ((t541 * t470 ^ 2 + (t537 * t467 + (t536 - t540) * t470) * t467) * qJD(2) + (t538 * t467 - t539 * t470) * qJD(1)) * t501 / 0.2e1 + ((t399 * t433 + t400 * t434 + t543 * t418 + t542 * t419 + t533 * t514) * t452 + (t376 * t433 + t378 * t434 + t547 * t418 + t545 * t419 + t535 * t514) * t436 + (t433 * t377 + t434 * t379 + t546 * t418 + t544 * t419 + t534 * t514) * t435) * t435 / 0.2e1 + ((t399 * t431 + t400 * t432 + t543 * t416 + t542 * t417 + t533 * t515) * t452 + (t431 * t376 + t432 * t378 + t547 * t416 + t545 * t417 + t535 * t515) * t436 + (t377 * t431 + t379 * t432 + t546 * t416 + t544 * t417 + t534 * t515) * t435) * t436 / 0.2e1 + ((-t534 * t435 - t535 * t436 - t533 * t452) * t459 + ((-t399 * t465 + t400 * t468 + t543 * t456 + t542 * t458) * t452 + (-t376 * t465 + t378 * t468 + t547 * t456 + t545 * t458) * t436 + (-t377 * t465 + t379 * t468 + t546 * t456 + t544 * t458) * t435) * t457) * t452 / 0.2e1;
T  = t1;
