% Calculate kinetic energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:14
% EndTime: 2019-03-09 08:17:18
% DurationCPUTime: 3.50s
% Computational Cost: add. (972->269), mult. (2378->393), div. (0->0), fcn. (2499->8), ass. (0->139)
t553 = Icges(3,4) + Icges(4,6);
t552 = Icges(3,1) + Icges(4,2);
t551 = Icges(3,2) + Icges(4,3);
t462 = cos(qJ(2));
t550 = t553 * t462;
t459 = sin(qJ(2));
t549 = t553 * t459;
t548 = -Icges(4,4) + Icges(3,5);
t547 = Icges(4,5) - Icges(3,6);
t546 = t551 * t459 - t550;
t545 = -t552 * t462 + t549;
t544 = Icges(4,1) + Icges(3,3);
t543 = Icges(5,1) + Icges(6,1);
t542 = Icges(5,4) - Icges(6,5);
t541 = Icges(6,4) + Icges(5,5);
t540 = Icges(5,2) + Icges(6,3);
t539 = Icges(6,6) - Icges(5,6);
t538 = Icges(5,3) + Icges(6,2);
t460 = sin(qJ(1));
t463 = cos(qJ(1));
t537 = -t546 * t460 + t547 * t463;
t536 = t547 * t460 + t546 * t463;
t535 = t545 * t460 + t548 * t463;
t534 = t548 * t460 - t545 * t463;
t533 = -t551 * t462 - t549;
t532 = t552 * t459 + t550;
t531 = t547 * t459 + t548 * t462;
t456 = sin(pkin(9));
t457 = cos(pkin(9));
t506 = t459 * t463;
t422 = t456 * t460 - t457 * t506;
t423 = t456 * t506 + t457 * t460;
t504 = t462 * t463;
t530 = t540 * t422 - t542 * t423 + t539 * t504;
t507 = t459 * t460;
t424 = t456 * t463 + t457 * t507;
t425 = t456 * t507 - t457 * t463;
t505 = t460 * t462;
t529 = -t540 * t424 - t542 * t425 + t539 * t505;
t528 = t539 * t422 + t541 * t423 + t538 * t504;
t527 = t539 * t424 - t541 * t425 - t538 * t505;
t526 = -t542 * t422 + t543 * t423 + t541 * t504;
t525 = t542 * t424 + t543 * t425 + t541 * t505;
t524 = (t542 * t456 + t540 * t457) * t462 + t539 * t459;
t523 = (-t541 * t456 + t539 * t457) * t462 + t538 * t459;
t522 = (-t543 * t456 - t542 * t457) * t462 + t541 * t459;
t521 = t544 * t460 + t531 * t463;
t520 = t531 * t460 - t544 * t463;
t519 = t548 * t459 - t547 * t462;
t518 = t533 * t459 + t532 * t462;
t517 = t537 * t459 + t535 * t462;
t516 = t536 * t459 + t534 * t462;
t484 = pkin(2) * t462 + qJ(3) * t459;
t427 = t484 * t460;
t447 = pkin(1) * t460 - pkin(7) * t463;
t503 = -t427 - t447;
t499 = qJD(3) * t459;
t453 = t463 * t499;
t498 = qJD(4) * t462;
t502 = t463 * t498 + t453;
t501 = qJD(2) * t460;
t500 = qJD(2) * t463;
t497 = qJD(6) * t462;
t496 = qJD(5) * t422 + t502;
t428 = t484 * t463;
t435 = qJD(1) * (pkin(1) * t463 + pkin(7) * t460);
t495 = qJD(1) * t428 + t460 * t499 + t435;
t434 = -pkin(3) * t463 + qJ(4) * t505;
t494 = -t434 + t503;
t442 = pkin(2) * t459 - qJ(3) * t462;
t491 = -qJ(4) * t459 - t442;
t490 = qJD(2) * (rSges(4,2) * t459 + rSges(4,3) * t462 - t442);
t386 = pkin(4) * t425 - qJ(5) * t424;
t489 = -t386 + t494;
t488 = -(-pkin(4) * t456 + qJ(5) * t457) * t462 + t491;
t487 = -qJD(3) * t462 + t427 * t501 + t428 * t500;
t486 = rSges(3,1) * t462 - rSges(3,2) * t459;
t485 = -rSges(4,2) * t462 + rSges(4,3) * t459;
t433 = pkin(3) * t460 + qJ(4) * t504;
t483 = qJD(1) * t433 + t460 * t498 + t495;
t470 = qJD(2) * (-rSges(5,3) * t459 - (-rSges(5,1) * t456 - rSges(5,2) * t457) * t462 + t491);
t385 = pkin(4) * t423 + qJ(5) * t422;
t469 = qJD(1) * t385 - qJD(5) * t424 + t483;
t468 = qJD(2) * (-rSges(6,2) * t459 - (-rSges(6,1) * t456 + rSges(6,3) * t457) * t462 + t488);
t467 = qJD(2) * (pkin(5) * t456 * t462 + pkin(8) * t459 + t488);
t466 = qJD(4) * t459 + t433 * t500 + t434 * t501 + t487;
t465 = qJD(5) * t462 * t457 + t385 * t500 + t386 * t501 + t466;
t461 = cos(qJ(6));
t458 = sin(qJ(6));
t454 = -qJD(6) * t459 + qJD(1);
t446 = rSges(2,1) * t463 - rSges(2,2) * t460;
t445 = rSges(2,1) * t460 + rSges(2,2) * t463;
t444 = rSges(3,1) * t459 + rSges(3,2) * t462;
t431 = -t460 * t497 - t500;
t430 = -t463 * t497 + t501;
t416 = (-t456 * t461 + t457 * t458) * t462;
t415 = (t456 * t458 + t457 * t461) * t462;
t414 = -rSges(4,1) * t463 + t460 * t485;
t413 = rSges(4,1) * t460 + t463 * t485;
t412 = rSges(3,3) * t460 + t463 * t486;
t411 = -rSges(3,3) * t463 + t460 * t486;
t388 = pkin(5) * t425 - pkin(8) * t505;
t387 = pkin(5) * t423 - pkin(8) * t504;
t383 = -t424 * t458 + t425 * t461;
t382 = -t424 * t461 - t425 * t458;
t381 = t422 * t458 + t423 * t461;
t380 = t422 * t461 - t423 * t458;
t377 = rSges(5,1) * t425 + rSges(5,2) * t424 + rSges(5,3) * t505;
t376 = rSges(6,1) * t425 + rSges(6,2) * t505 - rSges(6,3) * t424;
t375 = rSges(5,1) * t423 - rSges(5,2) * t422 + rSges(5,3) * t504;
t374 = rSges(6,1) * t423 + rSges(6,2) * t504 + rSges(6,3) * t422;
t361 = rSges(7,1) * t416 + rSges(7,2) * t415 - rSges(7,3) * t459;
t360 = Icges(7,1) * t416 + Icges(7,4) * t415 - Icges(7,5) * t459;
t359 = Icges(7,4) * t416 + Icges(7,2) * t415 - Icges(7,6) * t459;
t358 = Icges(7,5) * t416 + Icges(7,6) * t415 - Icges(7,3) * t459;
t357 = qJD(1) * t412 - t444 * t501 + t435;
t356 = -t444 * t500 + (-t411 - t447) * qJD(1);
t355 = (t411 * t460 + t412 * t463) * qJD(2);
t354 = qJD(1) * t413 + t460 * t490 + t495;
t353 = t453 + t463 * t490 + (-t414 + t503) * qJD(1);
t352 = rSges(7,1) * t383 + rSges(7,2) * t382 - rSges(7,3) * t505;
t351 = rSges(7,1) * t381 + rSges(7,2) * t380 - rSges(7,3) * t504;
t350 = Icges(7,1) * t383 + Icges(7,4) * t382 - Icges(7,5) * t505;
t349 = Icges(7,1) * t381 + Icges(7,4) * t380 - Icges(7,5) * t504;
t348 = Icges(7,4) * t383 + Icges(7,2) * t382 - Icges(7,6) * t505;
t347 = Icges(7,4) * t381 + Icges(7,2) * t380 - Icges(7,6) * t504;
t346 = Icges(7,5) * t383 + Icges(7,6) * t382 - Icges(7,3) * t505;
t345 = Icges(7,5) * t381 + Icges(7,6) * t380 - Icges(7,3) * t504;
t344 = (t413 * t463 + t414 * t460) * qJD(2) + t487;
t343 = qJD(1) * t375 + t460 * t470 + t483;
t342 = t463 * t470 + (-t377 + t494) * qJD(1) + t502;
t341 = (t375 * t463 + t377 * t460) * qJD(2) + t466;
t340 = qJD(1) * t374 + t460 * t468 + t469;
t339 = t463 * t468 + (-t376 + t489) * qJD(1) + t496;
t338 = (t374 * t463 + t376 * t460) * qJD(2) + t465;
t337 = qJD(1) * t387 + t351 * t454 - t361 * t430 + t460 * t467 + t469;
t336 = -t352 * t454 + t361 * t431 + t463 * t467 + (-t388 + t489) * qJD(1) + t496;
t335 = -t351 * t431 + t352 * t430 + (t387 * t463 + t388 * t460) * qJD(2) + t465;
t1 = t431 * ((-t345 * t505 + t347 * t382 + t349 * t383) * t430 + (-t346 * t505 + t348 * t382 + t350 * t383) * t431 + (-t358 * t505 + t359 * t382 + t360 * t383) * t454) / 0.2e1 + t454 * ((-t345 * t459 + t347 * t415 + t349 * t416) * t430 + (-t346 * t459 + t348 * t415 + t350 * t416) * t431 + (-t459 * t358 + t415 * t359 + t416 * t360) * t454) / 0.2e1 + t430 * ((-t345 * t504 + t380 * t347 + t381 * t349) * t430 + (-t346 * t504 + t348 * t380 + t350 * t381) * t431 + (-t358 * t504 + t359 * t380 + t360 * t381) * t454) / 0.2e1 + m(4) * (t344 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(5) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(6) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(3) * (t355 ^ 2 + t356 ^ 2 + t357 ^ 2) / 0.2e1 + m(7) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + (m(2) * (t445 ^ 2 + t446 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((((t456 * t525 - t457 * t529 - t537) * t463 + (-t526 * t456 + t457 * t530 - t536) * t460) * t462 + ((t527 + t535) * t463 + (t528 + t534) * t460) * t459) * qJD(2) + ((-t456 * t522 + t457 * t524 - t533) * t462 + (t523 + t532) * t459) * qJD(1)) * qJD(1) / 0.2e1 + (((-t529 * t422 - t525 * t423 + t517 * t463 + t527 * t504) * t463 + (t530 * t422 + t526 * t423 + t528 * t504 + (t516 - t520) * t463 + t521 * t460) * t460) * qJD(2) + (t422 * t524 + t423 * t522 + t460 * t519 + t463 * t518 + t504 * t523) * qJD(1)) * t501 / 0.2e1 - (((-t530 * t424 + t526 * t425 + t516 * t460 + t528 * t505) * t460 + (t529 * t424 - t525 * t425 + t527 * t505 + (t517 - t521) * t460 + t520 * t463) * t463) * qJD(2) + (-t424 * t524 + t425 * t522 + t460 * t518 - t463 * t519 + t505 * t523) * qJD(1)) * t500 / 0.2e1;
T  = t1;
