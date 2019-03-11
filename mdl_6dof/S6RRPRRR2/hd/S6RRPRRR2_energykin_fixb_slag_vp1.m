% Calculate kinetic energy for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:38
% EndTime: 2019-03-09 13:16:41
% DurationCPUTime: 2.53s
% Computational Cost: add. (2096->308), mult. (1931->487), div. (0->0), fcn. (1832->12), ass. (0->170)
t555 = Icges(3,3) + Icges(4,3);
t471 = qJ(2) + pkin(11);
t461 = sin(t471);
t462 = cos(t471);
t475 = sin(qJ(2));
t478 = cos(qJ(2));
t554 = Icges(3,5) * t478 + Icges(4,5) * t462 - Icges(3,6) * t475 - Icges(4,6) * t461;
t476 = sin(qJ(1));
t479 = cos(qJ(1));
t553 = t554 * t476 - t555 * t479;
t552 = t555 * t476 + t554 * t479;
t551 = Icges(3,5) * t475 + Icges(4,5) * t461 + Icges(3,6) * t478 + Icges(4,6) * t462;
t536 = Icges(4,4) * t461;
t440 = Icges(4,2) * t462 + t536;
t535 = Icges(4,4) * t462;
t441 = Icges(4,1) * t461 + t535;
t538 = Icges(3,4) * t475;
t451 = Icges(3,2) * t478 + t538;
t537 = Icges(3,4) * t478;
t452 = Icges(3,1) * t475 + t537;
t550 = -t440 * t461 + t441 * t462 - t451 * t475 + t452 * t478;
t501 = -Icges(4,2) * t461 + t535;
t407 = Icges(4,6) * t476 + t479 * t501;
t504 = Icges(4,1) * t462 - t536;
t409 = Icges(4,5) * t476 + t479 * t504;
t502 = -Icges(3,2) * t475 + t537;
t424 = Icges(3,6) * t476 + t479 * t502;
t505 = Icges(3,1) * t478 - t538;
t426 = Icges(3,5) * t476 + t479 * t505;
t549 = -t407 * t461 + t409 * t462 - t424 * t475 + t426 * t478;
t406 = -Icges(4,6) * t479 + t476 * t501;
t408 = -Icges(4,5) * t479 + t476 * t504;
t423 = -Icges(3,6) * t479 + t476 * t502;
t425 = -Icges(3,5) * t479 + t476 * t505;
t548 = t406 * t461 - t408 * t462 + t423 * t475 - t425 * t478;
t544 = pkin(2) * t475;
t542 = t478 * pkin(2);
t477 = cos(qJ(5));
t541 = pkin(5) * t477;
t463 = qJ(4) + t471;
t457 = sin(t463);
t534 = Icges(5,4) * t457;
t458 = cos(t463);
t533 = Icges(5,4) * t458;
t532 = t457 * t476;
t531 = t457 * t479;
t472 = qJ(5) + qJ(6);
t467 = sin(t472);
t530 = t467 * t476;
t529 = t467 * t479;
t468 = cos(t472);
t528 = t468 * t476;
t527 = t468 * t479;
t474 = sin(qJ(5));
t526 = t474 * t476;
t525 = t474 * t479;
t524 = t476 * t477;
t523 = t477 * t479;
t402 = -qJ(3) * t479 + t476 * t542;
t403 = qJ(3) * t476 + t479 * t542;
t466 = qJD(2) * t476;
t518 = qJD(2) * t479;
t522 = t402 * t466 + t403 * t518;
t456 = pkin(1) * t476 - pkin(7) * t479;
t521 = -t402 - t456;
t520 = pkin(3) * t462;
t448 = qJD(4) * t476 + t466;
t517 = qJD(5) * t457;
t516 = qJD(6) * t457;
t376 = -pkin(8) * t479 + t476 * t520;
t515 = -t376 + t521;
t419 = t479 * t517 + t448;
t449 = (-qJD(2) - qJD(4)) * t479;
t377 = pkin(8) * t476 + t479 * t520;
t512 = t376 * t466 + t377 * t518 + t522;
t511 = pkin(4) * t458 + pkin(9) * t457;
t444 = qJD(1) * (pkin(1) * t479 + pkin(7) * t476);
t510 = qJD(1) * t403 - qJD(3) * t479 + t444;
t509 = rSges(3,1) * t478 - rSges(3,2) * t475;
t508 = rSges(4,1) * t462 - rSges(4,2) * t461;
t507 = rSges(5,1) * t458 - rSges(5,2) * t457;
t506 = qJD(2) * (-rSges(4,1) * t461 - rSges(4,2) * t462 - t544);
t503 = Icges(5,1) * t458 - t534;
t500 = -Icges(5,2) * t457 + t533;
t497 = Icges(5,5) * t458 - Icges(5,6) * t457;
t420 = t476 * t517 + t449;
t490 = qJD(2) * (-pkin(3) * t461 - t544);
t417 = t511 * t476;
t418 = t511 * t479;
t489 = t448 * t417 - t418 * t449 + t512;
t488 = pkin(10) * t457 + t458 * t541;
t487 = qJD(1) * (Icges(5,5) * t457 + Icges(5,6) * t458) + (-Icges(5,3) * t479 + t476 * t497) * t449 + (Icges(5,3) * t476 + t479 * t497) * t448;
t465 = qJD(3) * t476;
t486 = t479 * t490 + t465;
t485 = qJD(1) * t377 + t476 * t490 + t510;
t438 = pkin(4) * t457 - pkin(9) * t458;
t484 = t449 * t438 + (-t417 + t515) * qJD(1) + t486;
t483 = qJD(1) * t418 - t438 * t448 + t485;
t395 = -Icges(5,6) * t479 + t476 * t500;
t396 = Icges(5,6) * t476 + t479 * t500;
t397 = -Icges(5,5) * t479 + t476 * t503;
t398 = Icges(5,5) * t476 + t479 * t503;
t435 = Icges(5,2) * t458 + t534;
t436 = Icges(5,1) * t457 + t533;
t482 = (-t396 * t457 + t398 * t458) * t448 + (-t395 * t457 + t397 * t458) * t449 + (-t435 * t457 + t436 * t458) * qJD(1);
t455 = rSges(2,1) * t479 - rSges(2,2) * t476;
t454 = rSges(2,1) * t476 + rSges(2,2) * t479;
t453 = rSges(3,1) * t475 + rSges(3,2) * t478;
t447 = -qJD(5) * t458 + qJD(1);
t437 = rSges(5,1) * t457 + rSges(5,2) * t458;
t433 = qJD(1) + (-qJD(5) - qJD(6)) * t458;
t432 = t458 * t523 + t526;
t431 = -t458 * t525 + t524;
t430 = t458 * t524 - t525;
t429 = -t458 * t526 - t523;
t428 = rSges(3,3) * t476 + t479 * t509;
t427 = -rSges(3,3) * t479 + t476 * t509;
t416 = t458 * t527 + t530;
t415 = -t458 * t529 + t528;
t414 = t458 * t528 - t529;
t413 = -t458 * t530 - t527;
t411 = rSges(4,3) * t476 + t479 * t508;
t410 = -rSges(4,3) * t479 + t476 * t508;
t401 = rSges(5,3) * t476 + t479 * t507;
t400 = -rSges(5,3) * t479 + t476 * t507;
t389 = -rSges(6,3) * t458 + (rSges(6,1) * t477 - rSges(6,2) * t474) * t457;
t388 = -Icges(6,5) * t458 + (Icges(6,1) * t477 - Icges(6,4) * t474) * t457;
t387 = -Icges(6,6) * t458 + (Icges(6,4) * t477 - Icges(6,2) * t474) * t457;
t386 = -Icges(6,3) * t458 + (Icges(6,5) * t477 - Icges(6,6) * t474) * t457;
t385 = t476 * t516 + t420;
t384 = t479 * t516 + t419;
t382 = -rSges(7,3) * t458 + (rSges(7,1) * t468 - rSges(7,2) * t467) * t457;
t381 = -Icges(7,5) * t458 + (Icges(7,1) * t468 - Icges(7,4) * t467) * t457;
t380 = -Icges(7,6) * t458 + (Icges(7,4) * t468 - Icges(7,2) * t467) * t457;
t379 = -Icges(7,3) * t458 + (Icges(7,5) * t468 - Icges(7,6) * t467) * t457;
t378 = -pkin(10) * t458 + t457 * t541;
t372 = qJD(1) * t428 - t453 * t466 + t444;
t371 = -t453 * t518 + (-t427 - t456) * qJD(1);
t370 = (t427 * t476 + t428 * t479) * qJD(2);
t369 = rSges(6,1) * t432 + rSges(6,2) * t431 + rSges(6,3) * t531;
t368 = rSges(6,1) * t430 + rSges(6,2) * t429 + rSges(6,3) * t532;
t367 = Icges(6,1) * t432 + Icges(6,4) * t431 + Icges(6,5) * t531;
t366 = Icges(6,1) * t430 + Icges(6,4) * t429 + Icges(6,5) * t532;
t365 = Icges(6,4) * t432 + Icges(6,2) * t431 + Icges(6,6) * t531;
t364 = Icges(6,4) * t430 + Icges(6,2) * t429 + Icges(6,6) * t532;
t363 = Icges(6,5) * t432 + Icges(6,6) * t431 + Icges(6,3) * t531;
t362 = Icges(6,5) * t430 + Icges(6,6) * t429 + Icges(6,3) * t532;
t361 = pkin(5) * t526 + t479 * t488;
t360 = -pkin(5) * t525 + t476 * t488;
t359 = rSges(7,1) * t416 + rSges(7,2) * t415 + rSges(7,3) * t531;
t358 = rSges(7,1) * t414 + rSges(7,2) * t413 + rSges(7,3) * t532;
t357 = Icges(7,1) * t416 + Icges(7,4) * t415 + Icges(7,5) * t531;
t356 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t532;
t355 = Icges(7,4) * t416 + Icges(7,2) * t415 + Icges(7,6) * t531;
t354 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t532;
t353 = Icges(7,5) * t416 + Icges(7,6) * t415 + Icges(7,3) * t531;
t352 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t532;
t351 = qJD(1) * t411 + t476 * t506 + t510;
t350 = t465 + t479 * t506 + (-t410 + t521) * qJD(1);
t349 = (t410 * t476 + t411 * t479) * qJD(2) + t522;
t348 = qJD(1) * t401 - t437 * t448 + t485;
t347 = t437 * t449 + (-t400 + t515) * qJD(1) + t486;
t346 = t400 * t448 - t401 * t449 + t512;
t345 = t369 * t447 - t389 * t419 + t483;
t344 = -t368 * t447 + t389 * t420 + t484;
t343 = t368 * t419 - t369 * t420 + t489;
t342 = t359 * t433 + t361 * t447 - t378 * t419 - t382 * t384 + t483;
t341 = -t358 * t433 - t360 * t447 + t378 * t420 + t382 * t385 + t484;
t340 = t358 * t384 - t359 * t385 + t360 * t419 - t361 * t420 + t489;
t1 = m(6) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(7) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + t448 * (t487 * t476 + t482 * t479) / 0.2e1 + t449 * (t482 * t476 - t487 * t479) / 0.2e1 + t419 * ((t363 * t531 + t431 * t365 + t432 * t367) * t419 + (t362 * t531 + t364 * t431 + t366 * t432) * t420 + (t386 * t531 + t387 * t431 + t388 * t432) * t447) / 0.2e1 + t420 * ((t363 * t532 + t365 * t429 + t367 * t430) * t419 + (t362 * t532 + t429 * t364 + t430 * t366) * t420 + (t386 * t532 + t387 * t429 + t388 * t430) * t447) / 0.2e1 + t447 * ((-t362 * t420 - t363 * t419 - t386 * t447) * t458 + ((-t365 * t474 + t367 * t477) * t419 + (-t364 * t474 + t366 * t477) * t420 + (-t387 * t474 + t388 * t477) * t447) * t457) / 0.2e1 + t384 * ((t353 * t531 + t415 * t355 + t416 * t357) * t384 + (t352 * t531 + t354 * t415 + t356 * t416) * t385 + (t379 * t531 + t380 * t415 + t381 * t416) * t433) / 0.2e1 + t385 * ((t353 * t532 + t355 * t413 + t357 * t414) * t384 + (t352 * t532 + t413 * t354 + t414 * t356) * t385 + (t379 * t532 + t380 * t413 + t381 * t414) * t433) / 0.2e1 + t433 * ((-t352 * t385 - t353 * t384 - t379 * t433) * t458 + ((-t355 * t467 + t357 * t468) * t384 + (-t354 * t467 + t356 * t468) * t385 + (-t380 * t467 + t381 * t468) * t433) * t457) / 0.2e1 + m(3) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(4) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(5) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t454 ^ 2 + t455 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t552 * t476 ^ 2 + (t548 * t479 + (t549 - t553) * t476) * t479) * qJD(2) + (t476 * t551 + t479 * t550) * qJD(1)) * t466 / 0.2e1 - ((t553 * t479 ^ 2 + (t549 * t476 + (t548 - t552) * t479) * t476) * qJD(2) + (t476 * t550 - t479 * t551) * qJD(1)) * t518 / 0.2e1 + ((t396 * t458 + t398 * t457) * t448 + (t395 * t458 + t397 * t457) * t449 + ((-t406 * t462 - t408 * t461 - t423 * t478 - t425 * t475) * t479 + (t407 * t462 + t461 * t409 + t424 * t478 + t426 * t475) * t476) * qJD(2) + (t458 * t435 + t457 * t436 + t462 * t440 + t461 * t441 + t478 * t451 + t475 * t452) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
