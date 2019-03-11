% Calculate kinetic energy for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:09
% EndTime: 2019-03-10 03:33:11
% DurationCPUTime: 2.31s
% Computational Cost: add. (2160->312), mult. (1995->506), div. (0->0), fcn. (1896->12), ass. (0->175)
t471 = qJ(2) + qJ(3);
t464 = sin(t471);
t541 = pkin(3) * t464;
t476 = cos(qJ(2));
t539 = t476 * pkin(2);
t475 = cos(qJ(5));
t538 = pkin(5) * t475;
t473 = sin(qJ(2));
t535 = Icges(3,4) * t473;
t534 = Icges(3,4) * t476;
t533 = Icges(4,4) * t464;
t466 = cos(t471);
t532 = Icges(4,4) * t466;
t467 = qJ(4) + t471;
t456 = sin(t467);
t531 = Icges(5,4) * t456;
t457 = cos(t467);
t530 = Icges(5,4) * t457;
t474 = sin(qJ(1));
t529 = t456 * t474;
t477 = cos(qJ(1));
t528 = t456 * t477;
t470 = qJ(5) + qJ(6);
t463 = sin(t470);
t527 = t463 * t474;
t526 = t463 * t477;
t465 = cos(t470);
t525 = t465 * t474;
t524 = t465 * t477;
t472 = sin(qJ(5));
t523 = t472 * t474;
t522 = t472 * t477;
t521 = t474 * t475;
t520 = t475 * t477;
t400 = -pkin(8) * t477 + t474 * t539;
t401 = pkin(8) * t474 + t477 * t539;
t462 = qJD(2) * t474;
t515 = qJD(2) * t477;
t519 = t400 * t462 + t401 * t515;
t455 = pkin(1) * t474 - pkin(7) * t477;
t518 = -t400 - t455;
t517 = pkin(3) * t466;
t446 = qJD(3) * t474 + t462;
t514 = qJD(5) * t456;
t513 = qJD(6) * t456;
t512 = -qJD(2) - qJD(3);
t511 = pkin(2) * qJD(2) * t473;
t372 = -pkin(9) * t477 + t474 * t517;
t510 = -t372 + t518;
t436 = qJD(4) * t474 + t446;
t509 = t477 * t511;
t396 = t477 * t514 + t436;
t508 = pkin(4) * t457 + pkin(10) * t456;
t507 = rSges(3,1) * t476 - rSges(3,2) * t473;
t506 = rSges(4,1) * t466 - rSges(4,2) * t464;
t505 = rSges(5,1) * t457 - rSges(5,2) * t456;
t437 = (-qJD(4) + t512) * t477;
t504 = Icges(3,1) * t476 - t535;
t503 = Icges(4,1) * t466 - t533;
t502 = Icges(5,1) * t457 - t531;
t501 = -Icges(3,2) * t473 + t534;
t500 = -Icges(4,2) * t464 + t532;
t499 = -Icges(5,2) * t456 + t530;
t498 = Icges(3,5) * t476 - Icges(3,6) * t473;
t497 = Icges(4,5) * t466 - Icges(4,6) * t464;
t496 = Icges(5,5) * t457 - Icges(5,6) * t456;
t419 = -Icges(3,6) * t477 + t474 * t501;
t421 = -Icges(3,5) * t477 + t474 * t504;
t495 = t419 * t473 - t421 * t476;
t420 = Icges(3,6) * t474 + t477 * t501;
t422 = Icges(3,5) * t474 + t477 * t504;
t494 = -t420 * t473 + t422 * t476;
t450 = Icges(3,2) * t476 + t535;
t451 = Icges(3,1) * t473 + t534;
t493 = -t450 * t473 + t451 * t476;
t447 = t512 * t477;
t492 = t447 * t541 - t509;
t373 = pkin(9) * t474 + t477 * t517;
t491 = t446 * t372 - t373 * t447 + t519;
t443 = qJD(1) * (pkin(1) * t477 + pkin(7) * t474);
t490 = qJD(1) * t401 - t474 * t511 + t443;
t397 = t474 * t514 + t437;
t489 = pkin(11) * t456 + t457 * t538;
t488 = (Icges(5,5) * t456 + Icges(5,6) * t457) * qJD(1) + (-Icges(5,3) * t477 + t474 * t496) * t437 + (Icges(5,3) * t474 + t477 * t496) * t436;
t487 = (Icges(4,5) * t464 + Icges(4,6) * t466) * qJD(1) + (-Icges(4,3) * t477 + t474 * t497) * t447 + (Icges(4,3) * t474 + t477 * t497) * t446;
t415 = t508 * t474;
t416 = t508 * t477;
t486 = t436 * t415 - t416 * t437 + t491;
t485 = qJD(1) * t373 - t446 * t541 + t490;
t435 = pkin(4) * t456 - pkin(10) * t457;
t484 = t437 * t435 + (-t415 + t510) * qJD(1) + t492;
t483 = qJD(1) * t416 - t435 * t436 + t485;
t392 = -Icges(5,6) * t477 + t474 * t499;
t393 = Icges(5,6) * t474 + t477 * t499;
t394 = -Icges(5,5) * t477 + t474 * t502;
t395 = Icges(5,5) * t474 + t477 * t502;
t432 = Icges(5,2) * t457 + t531;
t433 = Icges(5,1) * t456 + t530;
t482 = (-t393 * t456 + t395 * t457) * t436 + (-t392 * t456 + t394 * t457) * t437 + (-t432 * t456 + t433 * t457) * qJD(1);
t404 = -Icges(4,6) * t477 + t474 * t500;
t405 = Icges(4,6) * t474 + t477 * t500;
t406 = -Icges(4,5) * t477 + t474 * t503;
t407 = Icges(4,5) * t474 + t477 * t503;
t439 = Icges(4,2) * t466 + t533;
t440 = Icges(4,1) * t464 + t532;
t481 = (-t405 * t464 + t407 * t466) * t446 + (-t404 * t464 + t406 * t466) * t447 + (-t439 * t464 + t440 * t466) * qJD(1);
t454 = rSges(2,1) * t477 - rSges(2,2) * t474;
t453 = rSges(2,1) * t474 + rSges(2,2) * t477;
t452 = rSges(3,1) * t473 + rSges(3,2) * t476;
t449 = Icges(3,5) * t473 + Icges(3,6) * t476;
t448 = -qJD(5) * t457 + qJD(1);
t441 = rSges(4,1) * t464 + rSges(4,2) * t466;
t434 = rSges(5,1) * t456 + rSges(5,2) * t457;
t429 = qJD(1) + (-qJD(5) - qJD(6)) * t457;
t428 = t457 * t520 + t523;
t427 = -t457 * t522 + t521;
t426 = t457 * t521 - t522;
t425 = -t457 * t523 - t520;
t424 = rSges(3,3) * t474 + t477 * t507;
t423 = -rSges(3,3) * t477 + t474 * t507;
t418 = Icges(3,3) * t474 + t477 * t498;
t417 = -Icges(3,3) * t477 + t474 * t498;
t414 = t457 * t524 + t527;
t413 = -t457 * t526 + t525;
t412 = t457 * t525 - t526;
t411 = -t457 * t527 - t524;
t409 = rSges(4,3) * t474 + t477 * t506;
t408 = -rSges(4,3) * t477 + t474 * t506;
t399 = rSges(5,3) * t474 + t477 * t505;
t398 = -rSges(5,3) * t477 + t474 * t505;
t388 = -rSges(6,3) * t457 + (rSges(6,1) * t475 - rSges(6,2) * t472) * t456;
t385 = -Icges(6,5) * t457 + (Icges(6,1) * t475 - Icges(6,4) * t472) * t456;
t384 = -Icges(6,6) * t457 + (Icges(6,4) * t475 - Icges(6,2) * t472) * t456;
t383 = -Icges(6,3) * t457 + (Icges(6,5) * t475 - Icges(6,6) * t472) * t456;
t381 = -rSges(7,3) * t457 + (rSges(7,1) * t465 - rSges(7,2) * t463) * t456;
t380 = -Icges(7,5) * t457 + (Icges(7,1) * t465 - Icges(7,4) * t463) * t456;
t379 = -Icges(7,6) * t457 + (Icges(7,4) * t465 - Icges(7,2) * t463) * t456;
t378 = -Icges(7,3) * t457 + (Icges(7,5) * t465 - Icges(7,6) * t463) * t456;
t377 = t474 * t513 + t397;
t376 = t477 * t513 + t396;
t374 = -pkin(11) * t457 + t456 * t538;
t370 = qJD(1) * t424 - t452 * t462 + t443;
t369 = -t452 * t515 + (-t423 - t455) * qJD(1);
t367 = (t423 * t474 + t424 * t477) * qJD(2);
t366 = rSges(6,1) * t428 + rSges(6,2) * t427 + rSges(6,3) * t528;
t365 = rSges(6,1) * t426 + rSges(6,2) * t425 + rSges(6,3) * t529;
t364 = Icges(6,1) * t428 + Icges(6,4) * t427 + Icges(6,5) * t528;
t363 = Icges(6,1) * t426 + Icges(6,4) * t425 + Icges(6,5) * t529;
t362 = Icges(6,4) * t428 + Icges(6,2) * t427 + Icges(6,6) * t528;
t361 = Icges(6,4) * t426 + Icges(6,2) * t425 + Icges(6,6) * t529;
t360 = Icges(6,5) * t428 + Icges(6,6) * t427 + Icges(6,3) * t528;
t359 = Icges(6,5) * t426 + Icges(6,6) * t425 + Icges(6,3) * t529;
t358 = pkin(5) * t523 + t477 * t489;
t357 = -pkin(5) * t522 + t474 * t489;
t356 = rSges(7,1) * t414 + rSges(7,2) * t413 + rSges(7,3) * t528;
t355 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t529;
t354 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t528;
t353 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t529;
t352 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t528;
t351 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t529;
t350 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t528;
t349 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t529;
t348 = qJD(1) * t409 - t441 * t446 + t490;
t347 = -t509 + t441 * t447 + (-t408 + t518) * qJD(1);
t346 = t408 * t446 - t409 * t447 + t519;
t345 = qJD(1) * t399 - t434 * t436 + t485;
t344 = t434 * t437 + (-t398 + t510) * qJD(1) + t492;
t343 = t398 * t436 - t399 * t437 + t491;
t342 = t366 * t448 - t388 * t396 + t483;
t341 = -t365 * t448 + t388 * t397 + t484;
t340 = t365 * t396 - t366 * t397 + t486;
t339 = t356 * t429 + t358 * t448 - t374 * t396 - t376 * t381 + t483;
t338 = -t355 * t429 - t357 * t448 + t374 * t397 + t377 * t381 + t484;
t337 = t355 * t376 - t356 * t377 + t357 * t396 - t358 * t397 + t486;
t1 = t436 * (t474 * t488 + t477 * t482) / 0.2e1 + t437 * (t474 * t482 - t488 * t477) / 0.2e1 + t396 * ((t360 * t528 + t427 * t362 + t428 * t364) * t396 + (t359 * t528 + t361 * t427 + t363 * t428) * t397 + (t383 * t528 + t384 * t427 + t385 * t428) * t448) / 0.2e1 + t397 * ((t360 * t529 + t362 * t425 + t364 * t426) * t396 + (t359 * t529 + t425 * t361 + t426 * t363) * t397 + (t383 * t529 + t384 * t425 + t385 * t426) * t448) / 0.2e1 + t448 * ((-t359 * t397 - t360 * t396 - t383 * t448) * t457 + ((-t362 * t472 + t364 * t475) * t396 + (-t361 * t472 + t363 * t475) * t397 + (-t384 * t472 + t385 * t475) * t448) * t456) / 0.2e1 + t376 * ((t350 * t528 + t413 * t352 + t414 * t354) * t376 + (t349 * t528 + t351 * t413 + t353 * t414) * t377 + (t378 * t528 + t379 * t413 + t380 * t414) * t429) / 0.2e1 + t377 * ((t350 * t529 + t352 * t411 + t354 * t412) * t376 + (t349 * t529 + t411 * t351 + t412 * t353) * t377 + (t378 * t529 + t379 * t411 + t380 * t412) * t429) / 0.2e1 + t429 * ((-t349 * t377 - t350 * t376 - t378 * t429) * t457 + ((-t352 * t463 + t354 * t465) * t376 + (-t351 * t463 + t353 * t465) * t377 + (-t379 * t463 + t380 * t465) * t429) * t456) / 0.2e1 + m(3) * (t367 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(5) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(6) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(7) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + t446 * (t474 * t487 + t477 * t481) / 0.2e1 + t447 * (t474 * t481 - t477 * t487) / 0.2e1 - ((-t477 * t449 + t474 * t493) * qJD(1) + (t477 ^ 2 * t417 + (t494 * t474 + (-t418 + t495) * t477) * t474) * qJD(2)) * t515 / 0.2e1 + ((t474 * t449 + t477 * t493) * qJD(1) + (t474 ^ 2 * t418 + (t495 * t477 + (-t417 + t494) * t474) * t477) * qJD(2)) * t462 / 0.2e1 + (m(2) * (t453 ^ 2 + t454 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t405 * t466 + t407 * t464) * t446 + (t404 * t466 + t406 * t464) * t447 + (t393 * t457 + t395 * t456) * t436 + (t392 * t457 + t394 * t456) * t437 + ((t420 * t476 + t422 * t473) * t474 - (t419 * t476 + t421 * t473) * t477) * qJD(2) + (t457 * t432 + t456 * t433 + t466 * t439 + t464 * t440 + t476 * t450 + t473 * t451) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
