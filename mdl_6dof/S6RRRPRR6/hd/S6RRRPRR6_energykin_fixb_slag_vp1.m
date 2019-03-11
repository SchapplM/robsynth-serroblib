% Calculate kinetic energy for
% S6RRRPRR6
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:05
% EndTime: 2019-03-09 18:25:08
% DurationCPUTime: 3.52s
% Computational Cost: add. (2217->356), mult. (2682->557), div. (0->0), fcn. (2695->12), ass. (0->170)
t547 = -Icges(5,3) - Icges(4,3);
t489 = qJ(3) + pkin(11);
t481 = sin(t489);
t482 = cos(t489);
t496 = cos(qJ(1));
t493 = sin(qJ(1));
t495 = cos(qJ(2));
t531 = t493 * t495;
t425 = -t481 * t531 - t482 * t496;
t426 = -t481 * t496 + t482 * t531;
t491 = sin(qJ(3));
t494 = cos(qJ(3));
t446 = -t491 * t531 - t494 * t496;
t534 = t491 * t496;
t447 = t494 * t531 - t534;
t492 = sin(qJ(2));
t533 = t492 * t493;
t546 = Icges(4,5) * t447 + Icges(5,5) * t426 + Icges(4,6) * t446 + Icges(5,6) * t425 - t547 * t533;
t530 = t495 * t496;
t427 = -t481 * t530 + t482 * t493;
t428 = t481 * t493 + t482 * t530;
t448 = -t491 * t530 + t493 * t494;
t535 = t491 * t493;
t449 = t494 * t530 + t535;
t532 = t492 * t496;
t545 = Icges(4,5) * t449 + Icges(5,5) * t428 + Icges(4,6) * t448 + Icges(5,6) * t427 - t547 * t532;
t544 = t547 * t495 + (Icges(4,5) * t494 + Icges(5,5) * t482 - Icges(4,6) * t491 - Icges(5,6) * t481) * t492;
t539 = t494 * pkin(3);
t537 = Icges(3,4) * t492;
t536 = Icges(3,4) * t495;
t483 = qJ(5) + t489;
t478 = cos(t483);
t529 = pkin(5) * t478;
t477 = sin(t483);
t528 = pkin(5) * t477;
t516 = pkin(2) * t495 + pkin(8) * t492;
t450 = t516 * t493;
t451 = t516 * t496;
t484 = qJD(2) * t493;
t523 = qJD(2) * t496;
t527 = t450 * t484 + t451 * t523;
t526 = pkin(4) * t482;
t522 = qJD(3) * t492;
t452 = t496 * t522 + t484;
t521 = qJD(4) * t492;
t520 = qJD(5) * t492;
t519 = qJD(6) * t492;
t518 = -qJD(3) - qJD(5);
t422 = t496 * t520 + t452;
t517 = pkin(4) * t481;
t453 = t493 * t522 - t523;
t423 = t493 * t520 + t453;
t515 = rSges(3,1) * t495 - rSges(3,2) * t492;
t514 = Icges(3,1) * t495 - t537;
t513 = -Icges(3,2) * t492 + t536;
t512 = Icges(3,5) * t495 - Icges(3,6) * t492;
t433 = -Icges(3,6) * t496 + t493 * t513;
t436 = -Icges(3,5) * t496 + t493 * t514;
t511 = t433 * t492 - t436 * t495;
t434 = Icges(3,6) * t493 + t496 * t513;
t437 = Icges(3,5) * t493 + t496 * t514;
t510 = -t434 * t492 + t437 * t495;
t459 = Icges(3,2) * t495 + t537;
t460 = Icges(3,1) * t492 + t536;
t509 = -t459 * t492 + t460 * t495;
t457 = qJD(1) * (pkin(1) * t496 + pkin(7) * t493);
t464 = pkin(2) * t492 - pkin(8) * t495;
t508 = qJD(1) * t451 - t464 * t484 + t457;
t506 = qJ(4) * t492 + t495 * t539;
t393 = -pkin(3) * t534 + t493 * t506;
t507 = -qJD(4) * t495 + t452 * t393 + t527;
t394 = pkin(3) * t535 + t496 * t506;
t474 = -qJD(3) * t495 + qJD(1);
t505 = t474 * t394 + t493 * t521 + t508;
t465 = pkin(1) * t493 - pkin(7) * t496;
t504 = (-t450 - t465) * qJD(1) - t464 * t523;
t503 = pkin(10) * t492 + t495 * t529;
t502 = pkin(9) * t492 + t495 * t526;
t413 = -qJ(4) * t495 + t492 * t539;
t501 = t453 * t413 + t496 * t521 + t504;
t354 = t493 * t502 - t496 * t517;
t355 = t493 * t517 + t496 * t502;
t500 = t452 * t354 + (-t355 - t394) * t453 + t507;
t397 = -pkin(9) * t495 + t492 * t526;
t499 = t474 * t355 + (-t397 - t413) * t452 + t505;
t498 = t453 * t397 + (-t354 - t393) * t474 + t501;
t479 = qJ(6) + t483;
t467 = cos(t479);
t466 = sin(t479);
t463 = rSges(2,1) * t496 - rSges(2,2) * t493;
t462 = rSges(2,1) * t493 + rSges(2,2) * t496;
t461 = rSges(3,1) * t492 + rSges(3,2) * t495;
t458 = Icges(3,5) * t492 + Icges(3,6) * t495;
t455 = t495 * t518 + qJD(1);
t445 = qJD(1) + (-qJD(6) + t518) * t495;
t440 = rSges(3,3) * t493 + t496 * t515;
t439 = -rSges(3,3) * t496 + t493 * t515;
t438 = -rSges(4,3) * t495 + (rSges(4,1) * t494 - rSges(4,2) * t491) * t492;
t435 = -Icges(4,5) * t495 + (Icges(4,1) * t494 - Icges(4,4) * t491) * t492;
t432 = -Icges(4,6) * t495 + (Icges(4,4) * t494 - Icges(4,2) * t491) * t492;
t431 = Icges(3,3) * t493 + t496 * t512;
t430 = -Icges(3,3) * t496 + t493 * t512;
t421 = t477 * t493 + t478 * t530;
t420 = -t477 * t530 + t478 * t493;
t419 = -t477 * t496 + t478 * t531;
t418 = -t477 * t531 - t478 * t496;
t417 = -rSges(5,3) * t495 + (rSges(5,1) * t482 - rSges(5,2) * t481) * t492;
t416 = -Icges(5,5) * t495 + (Icges(5,1) * t482 - Icges(5,4) * t481) * t492;
t415 = -Icges(5,6) * t495 + (Icges(5,4) * t482 - Icges(5,2) * t481) * t492;
t412 = t466 * t493 + t467 * t530;
t411 = -t466 * t530 + t467 * t493;
t410 = -t466 * t496 + t467 * t531;
t409 = -t466 * t531 - t467 * t496;
t408 = -rSges(6,3) * t495 + (rSges(6,1) * t478 - rSges(6,2) * t477) * t492;
t407 = -Icges(6,5) * t495 + (Icges(6,1) * t478 - Icges(6,4) * t477) * t492;
t406 = -Icges(6,6) * t495 + (Icges(6,4) * t478 - Icges(6,2) * t477) * t492;
t405 = -Icges(6,3) * t495 + (Icges(6,5) * t478 - Icges(6,6) * t477) * t492;
t404 = t493 * t519 + t423;
t403 = t496 * t519 + t422;
t402 = -rSges(7,3) * t495 + (rSges(7,1) * t467 - rSges(7,2) * t466) * t492;
t401 = -Icges(7,5) * t495 + (Icges(7,1) * t467 - Icges(7,4) * t466) * t492;
t400 = -Icges(7,6) * t495 + (Icges(7,4) * t467 - Icges(7,2) * t466) * t492;
t399 = -Icges(7,3) * t495 + (Icges(7,5) * t467 - Icges(7,6) * t466) * t492;
t396 = rSges(4,1) * t449 + rSges(4,2) * t448 + rSges(4,3) * t532;
t395 = rSges(4,1) * t447 + rSges(4,2) * t446 + rSges(4,3) * t533;
t392 = Icges(4,1) * t449 + Icges(4,4) * t448 + Icges(4,5) * t532;
t391 = Icges(4,1) * t447 + Icges(4,4) * t446 + Icges(4,5) * t533;
t390 = Icges(4,4) * t449 + Icges(4,2) * t448 + Icges(4,6) * t532;
t389 = Icges(4,4) * t447 + Icges(4,2) * t446 + Icges(4,6) * t533;
t386 = qJD(1) * t440 - t461 * t484 + t457;
t385 = -t461 * t523 + (-t439 - t465) * qJD(1);
t384 = (t439 * t493 + t440 * t496) * qJD(2);
t382 = rSges(5,1) * t428 + rSges(5,2) * t427 + rSges(5,3) * t532;
t381 = rSges(5,1) * t426 + rSges(5,2) * t425 + rSges(5,3) * t533;
t380 = Icges(5,1) * t428 + Icges(5,4) * t427 + Icges(5,5) * t532;
t379 = Icges(5,1) * t426 + Icges(5,4) * t425 + Icges(5,5) * t533;
t378 = Icges(5,4) * t428 + Icges(5,2) * t427 + Icges(5,6) * t532;
t377 = Icges(5,4) * t426 + Icges(5,2) * t425 + Icges(5,6) * t533;
t373 = -pkin(10) * t495 + t492 * t529;
t371 = rSges(6,1) * t421 + rSges(6,2) * t420 + rSges(6,3) * t532;
t370 = rSges(6,1) * t419 + rSges(6,2) * t418 + rSges(6,3) * t533;
t369 = Icges(6,1) * t421 + Icges(6,4) * t420 + Icges(6,5) * t532;
t368 = Icges(6,1) * t419 + Icges(6,4) * t418 + Icges(6,5) * t533;
t367 = Icges(6,4) * t421 + Icges(6,2) * t420 + Icges(6,6) * t532;
t366 = Icges(6,4) * t419 + Icges(6,2) * t418 + Icges(6,6) * t533;
t365 = Icges(6,5) * t421 + Icges(6,6) * t420 + Icges(6,3) * t532;
t364 = Icges(6,5) * t419 + Icges(6,6) * t418 + Icges(6,3) * t533;
t363 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t532;
t362 = rSges(7,1) * t410 + rSges(7,2) * t409 + rSges(7,3) * t533;
t361 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t532;
t360 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t533;
t359 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t532;
t358 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t533;
t357 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t532;
t356 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t533;
t351 = t493 * t528 + t496 * t503;
t350 = t493 * t503 - t496 * t528;
t349 = t396 * t474 - t438 * t452 + t508;
t348 = -t395 * t474 + t438 * t453 + t504;
t347 = t395 * t452 - t396 * t453 + t527;
t346 = t382 * t474 + (-t413 - t417) * t452 + t505;
t345 = t417 * t453 + (-t381 - t393) * t474 + t501;
t344 = t381 * t452 + (-t382 - t394) * t453 + t507;
t343 = t371 * t455 - t408 * t422 + t499;
t342 = -t370 * t455 + t408 * t423 + t498;
t341 = t370 * t422 - t371 * t423 + t500;
t340 = t351 * t455 + t363 * t445 - t373 * t422 - t402 * t403 + t499;
t339 = -t350 * t455 - t362 * t445 + t373 * t423 + t402 * t404 + t498;
t338 = t350 * t422 - t351 * t423 + t362 * t403 - t363 * t404 + t500;
t1 = ((t493 * t458 + t496 * t509) * qJD(1) + (t493 ^ 2 * t431 + (t511 * t496 + (-t430 + t510) * t493) * t496) * qJD(2)) * t484 / 0.2e1 - ((-t496 * t458 + t493 * t509) * qJD(1) + (t496 ^ 2 * t430 + (t510 * t493 + (-t431 + t511) * t496) * t493) * qJD(2)) * t523 / 0.2e1 + qJD(1) * ((t459 * t495 + t460 * t492) * qJD(1) + ((t434 * t495 + t437 * t492) * t493 - (t433 * t495 + t436 * t492) * t496) * qJD(2)) / 0.2e1 + m(7) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(6) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(4) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(5) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + t403 * ((t357 * t532 + t359 * t411 + t361 * t412) * t403 + (t356 * t532 + t358 * t411 + t360 * t412) * t404 + (t399 * t532 + t400 * t411 + t401 * t412) * t445) / 0.2e1 + t445 * ((-t356 * t404 - t357 * t403 - t399 * t445) * t495 + ((-t359 * t466 + t361 * t467) * t403 + (-t358 * t466 + t360 * t467) * t404 + (-t400 * t466 + t401 * t467) * t445) * t492) / 0.2e1 + t404 * ((t357 * t533 + t359 * t409 + t361 * t410) * t403 + (t356 * t533 + t358 * t409 + t360 * t410) * t404 + (t399 * t533 + t400 * t409 + t401 * t410) * t445) / 0.2e1 + t422 * ((t365 * t532 + t367 * t420 + t369 * t421) * t422 + (t364 * t532 + t366 * t420 + t368 * t421) * t423 + (t405 * t532 + t406 * t420 + t407 * t421) * t455) / 0.2e1 + t423 * ((t365 * t533 + t367 * t418 + t369 * t419) * t422 + (t364 * t533 + t366 * t418 + t368 * t419) * t423 + (t405 * t533 + t406 * t418 + t407 * t419) * t455) / 0.2e1 + t455 * ((-t364 * t423 - t365 * t422 - t405 * t455) * t495 + ((-t367 * t477 + t369 * t478) * t422 + (-t366 * t477 + t368 * t478) * t423 + (-t406 * t477 + t407 * t478) * t455) * t492) / 0.2e1 + ((t415 * t427 + t416 * t428 + t432 * t448 + t435 * t449 + t532 * t544) * t474 + (t377 * t427 + t379 * t428 + t389 * t448 + t391 * t449 + t532 * t546) * t453 + (t378 * t427 + t380 * t428 + t390 * t448 + t392 * t449 + t545 * t532) * t452) * t452 / 0.2e1 + ((t415 * t425 + t416 * t426 + t432 * t446 + t435 * t447 + t533 * t544) * t474 + (t377 * t425 + t379 * t426 + t389 * t446 + t391 * t447 + t533 * t546) * t453 + (t378 * t425 + t380 * t426 + t390 * t446 + t392 * t447 + t533 * t545) * t452) * t453 / 0.2e1 + ((-t545 * t452 - t453 * t546 - t544 * t474) * t495 + ((-t415 * t481 + t416 * t482 - t432 * t491 + t435 * t494) * t474 + (-t377 * t481 + t379 * t482 - t389 * t491 + t391 * t494) * t453 + (-t378 * t481 + t380 * t482 - t390 * t491 + t392 * t494) * t452) * t492) * t474 / 0.2e1 + (Icges(2,3) + m(2) * (t462 ^ 2 + t463 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
