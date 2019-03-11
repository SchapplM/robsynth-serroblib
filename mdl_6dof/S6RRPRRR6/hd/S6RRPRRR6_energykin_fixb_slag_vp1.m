% Calculate kinetic energy for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:27
% EndTime: 2019-03-09 13:50:31
% DurationCPUTime: 3.74s
% Computational Cost: add. (1547->294), mult. (2509->457), div. (0->0), fcn. (2666->10), ass. (0->155)
t560 = Icges(3,4) - Icges(4,5);
t559 = Icges(3,1) + Icges(4,1);
t558 = Icges(3,2) + Icges(4,3);
t483 = cos(qJ(2));
t557 = t560 * t483;
t479 = sin(qJ(2));
t556 = t560 * t479;
t555 = Icges(4,4) + Icges(3,5);
t554 = Icges(3,6) - Icges(4,6);
t553 = t558 * t479 - t557;
t552 = t559 * t483 - t556;
t551 = Icges(4,2) + Icges(3,3);
t480 = sin(qJ(1));
t484 = cos(qJ(1));
t550 = t553 * t480 + t554 * t484;
t549 = -t554 * t480 + t553 * t484;
t548 = -t552 * t480 + t555 * t484;
t547 = t555 * t480 + t552 * t484;
t546 = -t558 * t483 - t556;
t545 = t559 * t479 + t557;
t544 = -t554 * t479 + t555 * t483;
t543 = t544 * t480 - t551 * t484;
t542 = t551 * t480 + t544 * t484;
t541 = t555 * t479 + t554 * t483;
t540 = t546 * t479 + t545 * t483;
t539 = t549 * t479 + t547 * t483;
t538 = -t550 * t479 + t548 * t483;
t523 = qJ(4) + qJ(5);
t476 = sin(t523);
t517 = cos(t523);
t512 = t479 * t517;
t441 = -t483 * t476 + t512;
t482 = cos(qJ(4));
t533 = pkin(4) * t482;
t478 = sin(qJ(4));
t527 = t478 * t479;
t526 = t478 * t483;
t524 = t483 * t484;
t508 = pkin(2) * t483 + qJ(3) * t479;
t446 = t508 * t480;
t469 = pkin(1) * t480 - pkin(7) * t484;
t522 = -t446 - t469;
t475 = qJD(2) * t480;
t521 = qJD(2) * t484;
t520 = qJD(3) * t479;
t447 = t508 * t484;
t453 = qJD(1) * (pkin(1) * t484 + pkin(7) * t480);
t519 = qJD(1) * t447 + t480 * t520 + t453;
t451 = pkin(3) * t480 * t483 + pkin(8) * t484;
t518 = -t451 + t522;
t464 = pkin(2) * t479 - qJ(3) * t483;
t514 = qJD(2) * (-rSges(4,1) * t479 + rSges(4,3) * t483 - t464);
t457 = qJD(4) * t484 - t521;
t492 = pkin(4) * t527 + t483 * t533;
t393 = pkin(9) * t484 + t480 * t492;
t513 = -t393 + t518;
t445 = qJD(5) * t484 + t457;
t511 = -qJD(3) * t483 + t446 * t475 + t447 * t521;
t510 = rSges(3,1) * t483 - rSges(3,2) * t479;
t509 = rSges(4,1) * t483 + rSges(4,3) * t479;
t507 = qJD(2) * (-pkin(3) * t479 - t464);
t450 = t479 * t482 - t526;
t494 = t482 * t483 + t527;
t444 = t475 + (-qJD(4) - qJD(5)) * t480;
t452 = pkin(3) * t524 - pkin(8) * t480;
t493 = t451 * t475 + t452 * t521 + t511;
t471 = t484 * t520;
t491 = t484 * t507 + t471;
t440 = t479 * t476 + t483 * t517;
t415 = -pkin(4) * t526 + t479 * t533;
t490 = t457 * t415 + t491;
t394 = -pkin(9) * t480 + t484 * t492;
t456 = -qJD(4) * t480 + t475;
t489 = t456 * t393 - t394 * t457 + t493;
t488 = qJD(1) * t452 + t480 * t507 + t519;
t487 = qJD(1) * t394 - t415 * t456 + t488;
t481 = cos(qJ(6));
t477 = sin(qJ(6));
t468 = rSges(2,1) * t484 - rSges(2,2) * t480;
t467 = rSges(2,1) * t480 + rSges(2,2) * t484;
t466 = rSges(3,1) * t479 + rSges(3,2) * t483;
t438 = t494 * t484;
t437 = t450 * t484;
t436 = t494 * t480;
t435 = t450 * t480;
t434 = rSges(3,3) * t480 + t484 * t510;
t433 = rSges(4,2) * t480 + t484 * t509;
t432 = -rSges(3,3) * t484 + t480 * t510;
t431 = -rSges(4,2) * t484 + t480 * t509;
t416 = qJD(6) * t440 + qJD(1);
t414 = t440 * t484;
t413 = t476 * t524 - t484 * t512;
t412 = t440 * t480;
t411 = t441 * t480;
t410 = rSges(5,1) * t450 - rSges(5,2) * t494;
t409 = Icges(5,1) * t450 - Icges(5,4) * t494;
t408 = Icges(5,4) * t450 - Icges(5,2) * t494;
t407 = Icges(5,5) * t450 - Icges(5,6) * t494;
t405 = t414 * t481 - t477 * t480;
t404 = -t414 * t477 - t480 * t481;
t403 = t412 * t481 + t477 * t484;
t402 = -t412 * t477 + t481 * t484;
t401 = pkin(5) * t441 + pkin(10) * t440;
t400 = rSges(6,1) * t441 - rSges(6,2) * t440;
t399 = Icges(6,1) * t441 - Icges(6,4) * t440;
t398 = Icges(6,4) * t441 - Icges(6,2) * t440;
t397 = Icges(6,5) * t441 - Icges(6,6) * t440;
t396 = -qJD(6) * t411 + t445;
t395 = qJD(6) * t413 + t444;
t391 = rSges(5,1) * t438 + rSges(5,2) * t437 - rSges(5,3) * t480;
t390 = rSges(5,1) * t436 + rSges(5,2) * t435 + rSges(5,3) * t484;
t389 = Icges(5,1) * t438 + Icges(5,4) * t437 - Icges(5,5) * t480;
t388 = Icges(5,1) * t436 + Icges(5,4) * t435 + Icges(5,5) * t484;
t387 = Icges(5,4) * t438 + Icges(5,2) * t437 - Icges(5,6) * t480;
t386 = Icges(5,4) * t436 + Icges(5,2) * t435 + Icges(5,6) * t484;
t385 = Icges(5,5) * t438 + Icges(5,6) * t437 - Icges(5,3) * t480;
t384 = Icges(5,5) * t436 + Icges(5,6) * t435 + Icges(5,3) * t484;
t383 = pkin(5) * t414 + pkin(10) * t413;
t382 = pkin(5) * t412 - pkin(10) * t411;
t381 = qJD(1) * t434 - t466 * t475 + t453;
t380 = -t466 * t521 + (-t432 - t469) * qJD(1);
t379 = (t432 * t480 + t434 * t484) * qJD(2);
t377 = rSges(6,1) * t414 - rSges(6,2) * t413 - rSges(6,3) * t480;
t376 = rSges(6,1) * t412 + rSges(6,2) * t411 + rSges(6,3) * t484;
t375 = Icges(6,1) * t414 - Icges(6,4) * t413 - Icges(6,5) * t480;
t374 = Icges(6,1) * t412 + Icges(6,4) * t411 + Icges(6,5) * t484;
t373 = Icges(6,4) * t414 - Icges(6,2) * t413 - Icges(6,6) * t480;
t372 = Icges(6,4) * t412 + Icges(6,2) * t411 + Icges(6,6) * t484;
t371 = Icges(6,5) * t414 - Icges(6,6) * t413 - Icges(6,3) * t480;
t370 = Icges(6,5) * t412 + Icges(6,6) * t411 + Icges(6,3) * t484;
t369 = rSges(7,3) * t440 + (rSges(7,1) * t481 - rSges(7,2) * t477) * t441;
t368 = Icges(7,5) * t440 + (Icges(7,1) * t481 - Icges(7,4) * t477) * t441;
t367 = Icges(7,6) * t440 + (Icges(7,4) * t481 - Icges(7,2) * t477) * t441;
t366 = Icges(7,3) * t440 + (Icges(7,5) * t481 - Icges(7,6) * t477) * t441;
t365 = qJD(1) * t433 + t480 * t514 + t519;
t364 = t471 + t484 * t514 + (-t431 + t522) * qJD(1);
t363 = rSges(7,1) * t405 + rSges(7,2) * t404 + rSges(7,3) * t413;
t362 = rSges(7,1) * t403 + rSges(7,2) * t402 - rSges(7,3) * t411;
t361 = Icges(7,1) * t405 + Icges(7,4) * t404 + Icges(7,5) * t413;
t360 = Icges(7,1) * t403 + Icges(7,4) * t402 - Icges(7,5) * t411;
t359 = Icges(7,4) * t405 + Icges(7,2) * t404 + Icges(7,6) * t413;
t358 = Icges(7,4) * t403 + Icges(7,2) * t402 - Icges(7,6) * t411;
t357 = Icges(7,5) * t405 + Icges(7,6) * t404 + Icges(7,3) * t413;
t356 = Icges(7,5) * t403 + Icges(7,6) * t402 - Icges(7,3) * t411;
t355 = (t431 * t480 + t433 * t484) * qJD(2) + t511;
t354 = qJD(1) * t391 - t410 * t456 + t488;
t353 = t410 * t457 + (-t390 + t518) * qJD(1) + t491;
t352 = t390 * t456 - t391 * t457 + t493;
t351 = qJD(1) * t377 - t400 * t444 + t487;
t350 = t400 * t445 + (-t376 + t513) * qJD(1) + t490;
t349 = t376 * t444 - t377 * t445 + t489;
t348 = qJD(1) * t383 + t363 * t416 - t369 * t395 - t401 * t444 + t487;
t347 = -t362 * t416 + t369 * t396 + t401 * t445 + (-t382 + t513) * qJD(1) + t490;
t346 = t362 * t395 - t363 * t396 + t382 * t444 - t383 * t445 + t489;
t1 = m(3) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(7) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(6) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(5) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(4) * (t355 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + t416 * ((t356 * t396 + t357 * t395 + t366 * t416) * t440 + ((-t359 * t477 + t361 * t481) * t395 + (-t358 * t477 + t360 * t481) * t396 + (-t367 * t477 + t368 * t481) * t416) * t441) / 0.2e1 + t395 * ((t413 * t357 + t404 * t359 + t405 * t361) * t395 + (t356 * t413 + t358 * t404 + t360 * t405) * t396 + (t366 * t413 + t367 * t404 + t368 * t405) * t416) / 0.2e1 + t396 * ((-t357 * t411 + t359 * t402 + t361 * t403) * t395 + (-t411 * t356 + t402 * t358 + t403 * t360) * t396 + (-t366 * t411 + t367 * t402 + t368 * t403) * t416) / 0.2e1 + t444 * ((-t480 * t371 - t413 * t373 + t414 * t375) * t444 + (-t370 * t480 - t372 * t413 + t374 * t414) * t445 + (-t397 * t480 - t398 * t413 + t399 * t414) * qJD(1)) / 0.2e1 + t445 * ((t371 * t484 + t373 * t411 + t375 * t412) * t444 + (t484 * t370 + t411 * t372 + t412 * t374) * t445 + (t397 * t484 + t398 * t411 + t399 * t412) * qJD(1)) / 0.2e1 + t456 * ((-t480 * t385 + t437 * t387 + t438 * t389) * t456 + (-t384 * t480 + t386 * t437 + t388 * t438) * t457 + (-t407 * t480 + t408 * t437 + t409 * t438) * qJD(1)) / 0.2e1 + t457 * ((t385 * t484 + t387 * t435 + t389 * t436) * t456 + (t484 * t384 + t435 * t386 + t436 * t388) * t457 + (t407 * t484 + t408 * t435 + t409 * t436) * qJD(1)) / 0.2e1 + (Icges(2,3) + m(2) * (t467 ^ 2 + t468 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t542 * t480 ^ 2 + (t538 * t484 + (t539 - t543) * t480) * t484) * qJD(2) + (t541 * t480 + t540 * t484) * qJD(1)) * t475 / 0.2e1 - ((t543 * t484 ^ 2 + (t539 * t480 + (t538 - t542) * t484) * t480) * qJD(2) + (t540 * t480 - t541 * t484) * qJD(1)) * t521 / 0.2e1 + ((-t373 * t440 + t375 * t441) * t444 + (-t372 * t440 + t374 * t441) * t445 + (-t387 * t494 + t389 * t450) * t456 + (-t386 * t494 + t388 * t450) * t457 + ((t548 * t479 + t550 * t483) * t484 + (t547 * t479 - t549 * t483) * t480) * qJD(2) + (-t440 * t398 + t441 * t399 - t494 * t408 + t450 * t409 + t545 * t479 - t546 * t483) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
