% Calculate kinetic energy for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:02
% EndTime: 2019-03-09 08:06:04
% DurationCPUTime: 2.61s
% Computational Cost: add. (1675->288), mult. (2392->422), div. (0->0), fcn. (2513->10), ass. (0->149)
t539 = Icges(5,1) + Icges(6,1);
t538 = Icges(5,4) - Icges(6,5);
t537 = Icges(6,4) + Icges(5,5);
t536 = Icges(5,2) + Icges(6,3);
t535 = -Icges(6,6) + Icges(5,6);
t534 = Icges(3,3) + Icges(4,3);
t533 = -Icges(5,3) - Icges(6,2);
t451 = qJ(2) + pkin(9);
t448 = sin(t451);
t449 = cos(t451);
t456 = sin(qJ(2));
t459 = cos(qJ(2));
t532 = Icges(3,5) * t459 + Icges(4,5) * t449 - Icges(3,6) * t456 - Icges(4,6) * t448;
t453 = cos(pkin(10));
t460 = cos(qJ(1));
t501 = t453 * t460;
t452 = sin(pkin(10));
t457 = sin(qJ(1));
t503 = t452 * t457;
t422 = t449 * t503 + t501;
t500 = t457 * t453;
t502 = t452 * t460;
t423 = t449 * t500 - t502;
t505 = t448 * t457;
t531 = -t536 * t422 + t538 * t423 + t535 * t505;
t424 = t449 * t502 - t500;
t425 = t449 * t501 + t503;
t504 = t448 * t460;
t530 = t536 * t424 - t538 * t425 - t535 * t504;
t529 = t535 * t422 - t537 * t423 + t533 * t505;
t528 = -t535 * t424 + t537 * t425 - t533 * t504;
t527 = t538 * t422 - t539 * t423 - t537 * t505;
t526 = -t538 * t424 + t539 * t425 + t537 * t504;
t525 = t535 * t449 + (t536 * t452 - t538 * t453) * t448;
t524 = t533 * t449 + (-t535 * t452 + t537 * t453) * t448;
t523 = -t537 * t449 + (-t538 * t452 + t539 * t453) * t448;
t522 = t532 * t457 - t534 * t460;
t521 = t534 * t457 + t532 * t460;
t520 = Icges(3,5) * t456 + Icges(4,5) * t448 + Icges(3,6) * t459 + Icges(4,6) * t449;
t507 = Icges(4,4) * t448;
t430 = Icges(4,2) * t449 + t507;
t506 = Icges(4,4) * t449;
t431 = Icges(4,1) * t448 + t506;
t509 = Icges(3,4) * t456;
t437 = Icges(3,2) * t459 + t509;
t508 = Icges(3,4) * t459;
t438 = Icges(3,1) * t456 + t508;
t519 = -t430 * t448 + t431 * t449 - t437 * t456 + t438 * t459;
t477 = -Icges(4,2) * t448 + t506;
t398 = -Icges(4,6) * t460 + t457 * t477;
t479 = Icges(4,1) * t449 - t507;
t400 = -Icges(4,5) * t460 + t457 * t479;
t478 = -Icges(3,2) * t456 + t508;
t414 = -Icges(3,6) * t460 + t457 * t478;
t480 = Icges(3,1) * t459 - t509;
t416 = -Icges(3,5) * t460 + t457 * t480;
t518 = t398 * t448 - t400 * t449 + t414 * t456 - t416 * t459;
t399 = Icges(4,6) * t457 + t460 * t477;
t401 = Icges(4,5) * t457 + t460 * t479;
t415 = Icges(3,6) * t457 + t460 * t478;
t417 = Icges(3,5) * t457 + t460 * t480;
t517 = -t399 * t448 + t401 * t449 - t415 * t456 + t417 * t459;
t513 = pkin(2) * t456;
t511 = pkin(2) * t459;
t394 = -qJ(3) * t460 + t457 * t511;
t395 = qJ(3) * t457 + t460 * t511;
t495 = qJD(2) * t460;
t496 = qJD(2) * t457;
t499 = t394 * t496 + t395 * t495;
t444 = pkin(1) * t457 - pkin(7) * t460;
t498 = -t394 - t444;
t450 = qJD(3) * t457;
t494 = qJD(4) * t448;
t497 = t460 * t494 + t450;
t493 = qJD(6) * t448;
t482 = pkin(3) * t449 + qJ(4) * t448;
t418 = t482 * t457;
t492 = -t418 + t498;
t491 = qJD(5) * t424 + t497;
t488 = -pkin(3) * t448 + qJ(4) * t449 - t513;
t379 = pkin(4) * t423 + qJ(5) * t422;
t487 = -t379 + t492;
t486 = -(pkin(4) * t453 + qJ(5) * t452) * t448 + t488;
t434 = qJD(1) * (pkin(1) * t460 + pkin(7) * t457);
t485 = qJD(1) * t395 - qJD(3) * t460 + t434;
t484 = rSges(3,1) * t459 - rSges(3,2) * t456;
t483 = rSges(4,1) * t449 - rSges(4,2) * t448;
t481 = qJD(2) * (-rSges(4,1) * t448 - rSges(4,2) * t449 - t513);
t468 = qJD(2) * (rSges(5,3) * t449 - (rSges(5,1) * t453 - rSges(5,2) * t452) * t448 + t488);
t419 = t482 * t460;
t467 = qJD(1) * t419 + t457 * t494 + t485;
t466 = -qJD(4) * t449 + t418 * t496 + t419 * t495 + t499;
t465 = qJD(2) * (rSges(6,2) * t449 - (rSges(6,1) * t453 + rSges(6,3) * t452) * t448 + t486);
t464 = qJD(2) * (-pkin(5) * t448 * t453 - pkin(8) * t449 + t486);
t380 = pkin(4) * t425 + qJ(5) * t424;
t463 = qJD(1) * t380 + qJD(5) * t422 + t467;
t462 = qJD(5) * t448 * t452 + t379 * t496 + t380 * t495 + t466;
t458 = cos(qJ(6));
t455 = sin(qJ(6));
t445 = qJD(6) * t449 + qJD(1);
t443 = rSges(2,1) * t460 - rSges(2,2) * t457;
t442 = rSges(2,1) * t457 + rSges(2,2) * t460;
t441 = rSges(3,1) * t456 + rSges(3,2) * t459;
t428 = -t457 * t493 - t495;
t427 = -t460 * t493 + t496;
t421 = rSges(3,3) * t457 + t460 * t484;
t420 = -rSges(3,3) * t460 + t457 * t484;
t407 = (t452 * t455 + t453 * t458) * t448;
t406 = (t452 * t458 - t453 * t455) * t448;
t405 = rSges(4,3) * t457 + t460 * t483;
t404 = -rSges(4,3) * t460 + t457 * t483;
t382 = pkin(5) * t425 - pkin(8) * t504;
t381 = pkin(5) * t423 - pkin(8) * t505;
t377 = t424 * t455 + t425 * t458;
t376 = t424 * t458 - t425 * t455;
t375 = t422 * t455 + t423 * t458;
t374 = t422 * t458 - t423 * t455;
t371 = qJD(1) * t421 - t441 * t496 + t434;
t370 = -t441 * t495 + (-t420 - t444) * qJD(1);
t369 = (t420 * t457 + t421 * t460) * qJD(2);
t368 = rSges(5,1) * t425 - rSges(5,2) * t424 + rSges(5,3) * t504;
t367 = rSges(6,1) * t425 + rSges(6,2) * t504 + rSges(6,3) * t424;
t366 = rSges(5,1) * t423 - rSges(5,2) * t422 + rSges(5,3) * t505;
t365 = rSges(6,1) * t423 + rSges(6,2) * t505 + rSges(6,3) * t422;
t352 = rSges(7,1) * t407 + rSges(7,2) * t406 + rSges(7,3) * t449;
t351 = Icges(7,1) * t407 + Icges(7,4) * t406 + Icges(7,5) * t449;
t350 = Icges(7,4) * t407 + Icges(7,2) * t406 + Icges(7,6) * t449;
t349 = Icges(7,5) * t407 + Icges(7,6) * t406 + Icges(7,3) * t449;
t348 = qJD(1) * t405 + t457 * t481 + t485;
t347 = t450 + t460 * t481 + (-t404 + t498) * qJD(1);
t346 = rSges(7,1) * t377 + rSges(7,2) * t376 - rSges(7,3) * t504;
t345 = rSges(7,1) * t375 + rSges(7,2) * t374 - rSges(7,3) * t505;
t344 = Icges(7,1) * t377 + Icges(7,4) * t376 - Icges(7,5) * t504;
t343 = Icges(7,1) * t375 + Icges(7,4) * t374 - Icges(7,5) * t505;
t342 = Icges(7,4) * t377 + Icges(7,2) * t376 - Icges(7,6) * t504;
t341 = Icges(7,4) * t375 + Icges(7,2) * t374 - Icges(7,6) * t505;
t340 = Icges(7,5) * t377 + Icges(7,6) * t376 - Icges(7,3) * t504;
t339 = Icges(7,5) * t375 + Icges(7,6) * t374 - Icges(7,3) * t505;
t338 = (t404 * t457 + t405 * t460) * qJD(2) + t499;
t337 = qJD(1) * t368 + t457 * t468 + t467;
t336 = t460 * t468 + (-t366 + t492) * qJD(1) + t497;
t335 = (t366 * t457 + t368 * t460) * qJD(2) + t466;
t334 = qJD(1) * t367 + t457 * t465 + t463;
t333 = t460 * t465 + (-t365 + t487) * qJD(1) + t491;
t332 = (t365 * t457 + t367 * t460) * qJD(2) + t462;
t331 = qJD(1) * t382 + t346 * t445 - t352 * t427 + t457 * t464 + t463;
t330 = -t345 * t445 + t352 * t428 + t460 * t464 + (-t381 + t487) * qJD(1) + t491;
t329 = t345 * t427 - t346 * t428 + (t381 * t457 + t382 * t460) * qJD(2) + t462;
t1 = t427 * ((-t340 * t504 + t376 * t342 + t377 * t344) * t427 + (-t339 * t504 + t341 * t376 + t343 * t377) * t428 + (-t349 * t504 + t350 * t376 + t351 * t377) * t445) / 0.2e1 + t428 * ((-t340 * t505 + t342 * t374 + t344 * t375) * t427 + (-t339 * t505 + t374 * t341 + t375 * t343) * t428 + (-t349 * t505 + t350 * t374 + t351 * t375) * t445) / 0.2e1 + t445 * ((t340 * t449 + t342 * t406 + t344 * t407) * t427 + (t339 * t449 + t341 * t406 + t343 * t407) * t428 + (t449 * t349 + t406 * t350 + t407 * t351) * t445) / 0.2e1 + m(5) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(6) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(7) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(3) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(4) * (t338 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + (m(2) * (t442 ^ 2 + t443 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t415 * t459 + t417 * t456) * t457 - (t414 * t459 + t416 * t456) * t460 + ((-t398 - t529) * t460 + (t399 - t528) * t457) * t449 + ((t531 * t452 + t527 * t453 - t400) * t460 + (t530 * t452 + t526 * t453 + t401) * t457) * t448) * qJD(2) + (t459 * t437 + t456 * t438 + (t430 - t524) * t449 + (t525 * t452 + t523 * t453 + t431) * t448) * qJD(1)) * qJD(1) / 0.2e1 + (((t531 * t424 + t527 * t425 + t518 * t460 + t529 * t504) * t460 + (t530 * t424 + t526 * t425 + t528 * t504 + (t517 - t522) * t460 + t521 * t457) * t457) * qJD(2) + (t525 * t424 + t523 * t425 + t520 * t457 + t519 * t460 + t524 * t504) * qJD(1)) * t496 / 0.2e1 - (((t530 * t422 + t526 * t423 + t517 * t457 + t528 * t505) * t457 + (t531 * t422 + t527 * t423 + t529 * t505 + (t518 - t521) * t457 + t522 * t460) * t460) * qJD(2) + (t525 * t422 + t523 * t423 + t519 * t457 - t520 * t460 + t524 * t505) * qJD(1)) * t495 / 0.2e1;
T  = t1;
