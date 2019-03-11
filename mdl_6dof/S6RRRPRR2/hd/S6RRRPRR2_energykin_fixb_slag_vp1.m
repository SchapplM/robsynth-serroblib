% Calculate kinetic energy for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:15
% EndTime: 2019-03-09 18:06:18
% DurationCPUTime: 2.73s
% Computational Cost: add. (2118->310), mult. (1953->484), div. (0->0), fcn. (1854->12), ass. (0->170)
t548 = Icges(4,3) + Icges(5,3);
t472 = qJ(2) + qJ(3);
t461 = pkin(11) + t472;
t457 = sin(t461);
t458 = cos(t461);
t466 = sin(t472);
t468 = cos(t472);
t547 = Icges(4,5) * t468 + Icges(5,5) * t458 - Icges(4,6) * t466 - Icges(5,6) * t457;
t475 = sin(qJ(1));
t478 = cos(qJ(1));
t530 = Icges(5,4) * t458;
t499 = -Icges(5,2) * t457 + t530;
t394 = -Icges(5,6) * t478 + t499 * t475;
t395 = Icges(5,6) * t475 + t499 * t478;
t531 = Icges(5,4) * t457;
t502 = Icges(5,1) * t458 - t531;
t396 = -Icges(5,5) * t478 + t502 * t475;
t397 = Icges(5,5) * t475 + t502 * t478;
t532 = Icges(4,4) * t468;
t500 = -Icges(4,2) * t466 + t532;
t405 = -Icges(4,6) * t478 + t500 * t475;
t406 = Icges(4,6) * t475 + t500 * t478;
t533 = Icges(4,4) * t466;
t503 = Icges(4,1) * t468 - t533;
t407 = -Icges(4,5) * t478 + t503 * t475;
t408 = Icges(4,5) * t475 + t503 * t478;
t434 = Icges(5,2) * t458 + t531;
t435 = Icges(5,1) * t457 + t530;
t440 = Icges(4,2) * t468 + t533;
t441 = Icges(4,1) * t466 + t532;
t464 = qJD(2) * t475;
t448 = qJD(3) * t475 + t464;
t449 = (-qJD(2) - qJD(3)) * t478;
t546 = (-t394 * t457 + t396 * t458 - t405 * t466 + t407 * t468) * t449 + (-t395 * t457 + t397 * t458 - t406 * t466 + t408 * t468) * t448 + (-t434 * t457 + t435 * t458 - t440 * t466 + t441 * t468) * qJD(1);
t545 = (t547 * t475 - t548 * t478) * t449 + (t548 * t475 + t547 * t478) * t448 + (Icges(4,5) * t466 + Icges(5,5) * t457 + Icges(4,6) * t468 + Icges(5,6) * t458) * qJD(1);
t541 = pkin(3) * t466;
t477 = cos(qJ(2));
t539 = t477 * pkin(2);
t476 = cos(qJ(5));
t538 = pkin(5) * t476;
t474 = sin(qJ(2));
t535 = Icges(3,4) * t474;
t534 = Icges(3,4) * t477;
t529 = t457 * t475;
t528 = t457 * t478;
t471 = qJ(5) + qJ(6);
t465 = sin(t471);
t527 = t465 * t475;
t526 = t465 * t478;
t467 = cos(t471);
t525 = t467 * t475;
t524 = t467 * t478;
t473 = sin(qJ(5));
t523 = t473 * t475;
t522 = t473 * t478;
t521 = t475 * t476;
t520 = t476 * t478;
t401 = -pkin(8) * t478 + t539 * t475;
t402 = pkin(8) * t475 + t539 * t478;
t515 = qJD(2) * t478;
t519 = t401 * t464 + t402 * t515;
t456 = pkin(1) * t475 - pkin(7) * t478;
t518 = -t401 - t456;
t517 = pkin(3) * t468;
t514 = qJD(5) * t457;
t513 = qJD(6) * t457;
t512 = pkin(2) * qJD(2) * t474;
t376 = -qJ(4) * t478 + t517 * t475;
t511 = t448 * t376 + t519;
t510 = -t376 + t518;
t418 = t478 * t514 + t448;
t509 = t478 * t512;
t508 = pkin(4) * t458 + pkin(9) * t457;
t507 = rSges(3,1) * t477 - rSges(3,2) * t474;
t506 = rSges(4,1) * t468 - rSges(4,2) * t466;
t505 = rSges(5,1) * t458 - rSges(5,2) * t457;
t504 = Icges(3,1) * t477 - t535;
t501 = -Icges(3,2) * t474 + t534;
t498 = Icges(3,5) * t477 - Icges(3,6) * t474;
t422 = -Icges(3,6) * t478 + t501 * t475;
t424 = -Icges(3,5) * t478 + t504 * t475;
t495 = t422 * t474 - t424 * t477;
t423 = Icges(3,6) * t475 + t501 * t478;
t425 = Icges(3,5) * t475 + t504 * t478;
t494 = -t423 * t474 + t425 * t477;
t451 = Icges(3,2) * t477 + t535;
t452 = Icges(3,1) * t474 + t534;
t493 = -t451 * t474 + t452 * t477;
t419 = t475 * t514 + t449;
t444 = qJD(1) * (pkin(1) * t478 + pkin(7) * t475);
t492 = qJD(1) * t402 - t475 * t512 + t444;
t491 = qJD(4) * t475 + t449 * t541 - t509;
t490 = pkin(10) * t457 + t538 * t458;
t377 = qJ(4) * t475 + t517 * t478;
t416 = t508 * t475;
t417 = t508 * t478;
t487 = t448 * t416 + (-t377 - t417) * t449 + t511;
t486 = qJD(1) * t377 - qJD(4) * t478 + t492;
t438 = pkin(4) * t457 - pkin(9) * t458;
t485 = t449 * t438 + (-t416 + t510) * qJD(1) + t491;
t484 = qJD(1) * t417 + (-t438 - t541) * t448 + t486;
t455 = rSges(2,1) * t478 - rSges(2,2) * t475;
t454 = rSges(2,1) * t475 + rSges(2,2) * t478;
t453 = rSges(3,1) * t474 + rSges(3,2) * t477;
t450 = Icges(3,5) * t474 + Icges(3,6) * t477;
t447 = -qJD(5) * t458 + qJD(1);
t442 = rSges(4,1) * t466 + rSges(4,2) * t468;
t436 = rSges(5,1) * t457 + rSges(5,2) * t458;
t432 = qJD(1) + (-qJD(5) - qJD(6)) * t458;
t431 = t458 * t520 + t523;
t430 = -t458 * t522 + t521;
t429 = t458 * t521 - t522;
t428 = -t458 * t523 - t520;
t427 = rSges(3,3) * t475 + t507 * t478;
t426 = -rSges(3,3) * t478 + t507 * t475;
t421 = Icges(3,3) * t475 + t498 * t478;
t420 = -Icges(3,3) * t478 + t498 * t475;
t415 = t458 * t524 + t527;
t414 = -t458 * t526 + t525;
t413 = t458 * t525 - t526;
t412 = -t458 * t527 - t524;
t411 = rSges(4,3) * t475 + t506 * t478;
t410 = -rSges(4,3) * t478 + t506 * t475;
t399 = rSges(5,3) * t475 + t505 * t478;
t398 = -rSges(5,3) * t478 + t505 * t475;
t388 = -rSges(6,3) * t458 + (rSges(6,1) * t476 - rSges(6,2) * t473) * t457;
t387 = -Icges(6,5) * t458 + (Icges(6,1) * t476 - Icges(6,4) * t473) * t457;
t386 = -Icges(6,6) * t458 + (Icges(6,4) * t476 - Icges(6,2) * t473) * t457;
t385 = -Icges(6,3) * t458 + (Icges(6,5) * t476 - Icges(6,6) * t473) * t457;
t384 = t475 * t513 + t419;
t383 = t478 * t513 + t418;
t381 = -rSges(7,3) * t458 + (rSges(7,1) * t467 - rSges(7,2) * t465) * t457;
t380 = -Icges(7,5) * t458 + (Icges(7,1) * t467 - Icges(7,4) * t465) * t457;
t379 = -Icges(7,6) * t458 + (Icges(7,4) * t467 - Icges(7,2) * t465) * t457;
t378 = -Icges(7,3) * t458 + (Icges(7,5) * t467 - Icges(7,6) * t465) * t457;
t375 = -pkin(10) * t458 + t538 * t457;
t373 = qJD(1) * t427 - t453 * t464 + t444;
t372 = -t453 * t515 + (-t426 - t456) * qJD(1);
t371 = (t426 * t475 + t427 * t478) * qJD(2);
t369 = rSges(6,1) * t431 + rSges(6,2) * t430 + rSges(6,3) * t528;
t368 = rSges(6,1) * t429 + rSges(6,2) * t428 + rSges(6,3) * t529;
t367 = Icges(6,1) * t431 + Icges(6,4) * t430 + Icges(6,5) * t528;
t366 = Icges(6,1) * t429 + Icges(6,4) * t428 + Icges(6,5) * t529;
t365 = Icges(6,4) * t431 + Icges(6,2) * t430 + Icges(6,6) * t528;
t364 = Icges(6,4) * t429 + Icges(6,2) * t428 + Icges(6,6) * t529;
t363 = Icges(6,5) * t431 + Icges(6,6) * t430 + Icges(6,3) * t528;
t362 = Icges(6,5) * t429 + Icges(6,6) * t428 + Icges(6,3) * t529;
t361 = pkin(5) * t523 + t490 * t478;
t360 = -pkin(5) * t522 + t490 * t475;
t359 = rSges(7,1) * t415 + rSges(7,2) * t414 + rSges(7,3) * t528;
t358 = rSges(7,1) * t413 + rSges(7,2) * t412 + rSges(7,3) * t529;
t357 = Icges(7,1) * t415 + Icges(7,4) * t414 + Icges(7,5) * t528;
t356 = Icges(7,1) * t413 + Icges(7,4) * t412 + Icges(7,5) * t529;
t355 = Icges(7,4) * t415 + Icges(7,2) * t414 + Icges(7,6) * t528;
t354 = Icges(7,4) * t413 + Icges(7,2) * t412 + Icges(7,6) * t529;
t353 = Icges(7,5) * t415 + Icges(7,6) * t414 + Icges(7,3) * t528;
t352 = Icges(7,5) * t413 + Icges(7,6) * t412 + Icges(7,3) * t529;
t351 = qJD(1) * t411 - t442 * t448 + t492;
t350 = -t509 + t442 * t449 + (-t410 + t518) * qJD(1);
t349 = t410 * t448 - t411 * t449 + t519;
t348 = qJD(1) * t399 + (-t436 - t541) * t448 + t486;
t347 = t436 * t449 + (-t398 + t510) * qJD(1) + t491;
t346 = t398 * t448 + (-t377 - t399) * t449 + t511;
t345 = t369 * t447 - t388 * t418 + t484;
t344 = -t368 * t447 + t388 * t419 + t485;
t343 = t368 * t418 - t369 * t419 + t487;
t342 = t359 * t432 + t361 * t447 - t375 * t418 - t381 * t383 + t484;
t341 = -t358 * t432 - t360 * t447 + t375 * t419 + t381 * t384 + t485;
t340 = t358 * t383 - t359 * t384 + t360 * t418 - t361 * t419 + t487;
t1 = -((-t478 * t450 + t493 * t475) * qJD(1) + (t478 ^ 2 * t420 + (t494 * t475 + (-t421 + t495) * t478) * t475) * qJD(2)) * t515 / 0.2e1 + ((t475 * t450 + t493 * t478) * qJD(1) + (t475 ^ 2 * t421 + (t495 * t478 + (-t420 + t494) * t475) * t478) * qJD(2)) * t464 / 0.2e1 + t418 * ((t363 * t528 + t430 * t365 + t431 * t367) * t418 + (t362 * t528 + t364 * t430 + t366 * t431) * t419 + (t385 * t528 + t386 * t430 + t387 * t431) * t447) / 0.2e1 + t419 * ((t363 * t529 + t365 * t428 + t367 * t429) * t418 + (t362 * t529 + t428 * t364 + t429 * t366) * t419 + (t385 * t529 + t386 * t428 + t387 * t429) * t447) / 0.2e1 + t447 * ((-t362 * t419 - t363 * t418 - t385 * t447) * t458 + ((-t365 * t473 + t367 * t476) * t418 + (-t364 * t473 + t366 * t476) * t419 + (-t386 * t473 + t387 * t476) * t447) * t457) / 0.2e1 + t383 * ((t353 * t528 + t414 * t355 + t415 * t357) * t383 + (t352 * t528 + t354 * t414 + t356 * t415) * t384 + (t378 * t528 + t379 * t414 + t380 * t415) * t432) / 0.2e1 + t384 * ((t353 * t529 + t355 * t412 + t357 * t413) * t383 + (t352 * t529 + t412 * t354 + t413 * t356) * t384 + (t378 * t529 + t379 * t412 + t380 * t413) * t432) / 0.2e1 + t432 * ((-t352 * t384 - t353 * t383 - t378 * t432) * t458 + ((-t355 * t465 + t357 * t467) * t383 + (-t354 * t465 + t356 * t467) * t384 + (-t379 * t465 + t380 * t467) * t432) * t457) / 0.2e1 + m(3) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(4) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(5) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(6) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(7) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + (t545 * t475 + t546 * t478) * t448 / 0.2e1 + (t546 * t475 - t545 * t478) * t449 / 0.2e1 + (Icges(2,3) + m(2) * (t454 ^ 2 + t455 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t423 * t477 + t425 * t474) * t475 - (t422 * t477 + t424 * t474) * t478) * qJD(2) + (t394 * t458 + t396 * t457 + t405 * t468 + t407 * t466) * t449 + (t395 * t458 + t397 * t457 + t406 * t468 + t408 * t466) * t448 + (t458 * t434 + t457 * t435 + t468 * t440 + t466 * t441 + t477 * t451 + t474 * t452) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
