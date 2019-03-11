% Calculate kinetic energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:13
% EndTime: 2019-03-08 22:30:15
% DurationCPUTime: 2.63s
% Computational Cost: add. (2432->334), mult. (5743->509), div. (0->0), fcn. (6987->12), ass. (0->154)
t561 = Icges(4,1) + Icges(5,2);
t560 = Icges(5,1) + Icges(4,3);
t559 = -Icges(4,4) - Icges(5,6);
t558 = Icges(5,4) - Icges(4,5);
t557 = Icges(5,5) - Icges(4,6);
t556 = Icges(4,2) + Icges(5,3);
t508 = sin(pkin(11));
t510 = cos(pkin(11));
t515 = cos(qJ(2));
t511 = cos(pkin(6));
t513 = sin(qJ(2));
t532 = t511 * t513;
t491 = t508 * t515 + t510 * t532;
t509 = sin(pkin(6));
t542 = cos(qJ(3));
t528 = t509 * t542;
t541 = sin(qJ(3));
t473 = t491 * t541 + t510 * t528;
t527 = t509 * t541;
t474 = t491 * t542 - t510 * t527;
t531 = t511 * t515;
t490 = t508 * t513 - t510 * t531;
t555 = t556 * t473 + t559 * t474 + t557 * t490;
t493 = -t508 * t532 + t510 * t515;
t475 = t493 * t541 - t508 * t528;
t476 = t493 * t542 + t508 * t527;
t492 = t508 * t531 + t510 * t513;
t554 = t556 * t475 + t559 * t476 + t557 * t492;
t553 = t557 * t473 - t558 * t474 + t560 * t490;
t552 = t557 * t475 - t558 * t476 + t560 * t492;
t551 = t559 * t473 + t561 * t474 - t558 * t490;
t550 = t559 * t475 + t561 * t476 - t558 * t492;
t494 = -t511 * t542 + t513 * t527;
t495 = t511 * t541 + t513 * t528;
t533 = t509 * t515;
t549 = t556 * t494 + t559 * t495 - t557 * t533;
t548 = t559 * t494 + t561 * t495 + t558 * t533;
t547 = t557 * t494 - t558 * t495 - t560 * t533;
t546 = qJD(2) ^ 2;
t514 = cos(qJ(5));
t540 = pkin(5) * t514;
t512 = sin(qJ(5));
t538 = t473 * t512;
t537 = t475 * t512;
t536 = t494 * t512;
t535 = t508 * t509;
t534 = t509 * t510;
t530 = qJD(2) * t509;
t501 = t508 * t530;
t479 = qJD(3) * t492 + t501;
t504 = qJD(2) * t511;
t430 = qJD(5) * t476 + t479;
t526 = t510 * t530;
t463 = pkin(2) * t491 + pkin(8) * t490;
t464 = pkin(2) * t493 + pkin(8) * t492;
t525 = t463 * t501 + t464 * t526 + qJD(1);
t480 = qJD(3) * t490 - t526;
t497 = -qJD(3) * t533 + t504;
t428 = pkin(3) * t474 + qJ(4) * t473;
t524 = qJD(4) * t494 + t479 * t428 + t525;
t496 = (pkin(2) * t513 - pkin(8) * t515) * t509;
t523 = t464 * t504 - t496 * t501;
t431 = qJD(5) * t474 + t480;
t472 = qJD(5) * t495 + t497;
t429 = pkin(3) * t476 + qJ(4) * t475;
t522 = qJD(4) * t473 + t497 * t429 + t523;
t521 = (-t463 * t511 - t496 * t534) * qJD(2);
t441 = t490 * pkin(4) + t474 * pkin(9);
t442 = t492 * pkin(4) + t476 * pkin(9);
t520 = t479 * t441 + (-t429 - t442) * t480 + t524;
t465 = pkin(3) * t495 + qJ(4) * t494;
t519 = qJD(4) * t475 + t480 * t465 + t521;
t481 = -pkin(4) * t533 + t495 * pkin(9);
t518 = t497 * t442 + (-t465 - t481) * t479 + t522;
t517 = t480 * t481 + (-t428 - t441) * t497 + t519;
t507 = qJ(5) + qJ(6);
t506 = cos(t507);
t505 = sin(t507);
t485 = rSges(3,3) * t511 + (rSges(3,1) * t513 + rSges(3,2) * t515) * t509;
t484 = Icges(3,5) * t511 + (Icges(3,1) * t513 + Icges(3,4) * t515) * t509;
t483 = Icges(3,6) * t511 + (Icges(3,4) * t513 + Icges(3,2) * t515) * t509;
t482 = Icges(3,3) * t511 + (Icges(3,5) * t513 + Icges(3,6) * t515) * t509;
t478 = -t514 * t533 + t536;
t477 = t494 * t514 + t512 * t533;
t467 = t494 * t505 - t506 * t533;
t466 = t494 * t506 + t505 * t533;
t461 = rSges(4,1) * t495 - rSges(4,2) * t494 - rSges(4,3) * t533;
t460 = -rSges(5,1) * t533 - rSges(5,2) * t495 + rSges(5,3) * t494;
t451 = rSges(3,1) * t493 - rSges(3,2) * t492 + rSges(3,3) * t535;
t450 = rSges(3,1) * t491 - rSges(3,2) * t490 - rSges(3,3) * t534;
t449 = Icges(3,1) * t493 - Icges(3,4) * t492 + Icges(3,5) * t535;
t448 = Icges(3,1) * t491 - Icges(3,4) * t490 - Icges(3,5) * t534;
t447 = Icges(3,4) * t493 - Icges(3,2) * t492 + Icges(3,6) * t535;
t446 = Icges(3,4) * t491 - Icges(3,2) * t490 - Icges(3,6) * t534;
t445 = Icges(3,5) * t493 - Icges(3,6) * t492 + Icges(3,3) * t535;
t444 = Icges(3,5) * t491 - Icges(3,6) * t490 - Icges(3,3) * t534;
t443 = qJD(6) * t495 + t472;
t439 = t492 * t514 + t537;
t438 = t475 * t514 - t492 * t512;
t437 = t490 * t514 + t538;
t436 = t473 * t514 - t490 * t512;
t435 = t475 * t505 + t492 * t506;
t434 = t475 * t506 - t492 * t505;
t433 = t473 * t505 + t490 * t506;
t432 = t473 * t506 - t490 * t505;
t425 = pkin(5) * t536 + pkin(10) * t495 - t533 * t540;
t422 = (-t450 * t511 - t485 * t534) * qJD(2);
t421 = (t451 * t511 - t485 * t535) * qJD(2);
t420 = rSges(6,1) * t478 + rSges(6,2) * t477 + rSges(6,3) * t495;
t419 = Icges(6,1) * t478 + Icges(6,4) * t477 + Icges(6,5) * t495;
t418 = Icges(6,4) * t478 + Icges(6,2) * t477 + Icges(6,6) * t495;
t417 = Icges(6,5) * t478 + Icges(6,6) * t477 + Icges(6,3) * t495;
t416 = rSges(4,1) * t476 - rSges(4,2) * t475 + rSges(4,3) * t492;
t415 = rSges(4,1) * t474 - rSges(4,2) * t473 + rSges(4,3) * t490;
t414 = rSges(5,1) * t492 - rSges(5,2) * t476 + rSges(5,3) * t475;
t413 = rSges(5,1) * t490 - rSges(5,2) * t474 + rSges(5,3) * t473;
t400 = rSges(7,1) * t467 + rSges(7,2) * t466 + rSges(7,3) * t495;
t399 = Icges(7,1) * t467 + Icges(7,4) * t466 + Icges(7,5) * t495;
t398 = Icges(7,4) * t467 + Icges(7,2) * t466 + Icges(7,6) * t495;
t397 = Icges(7,5) * t467 + Icges(7,6) * t466 + Icges(7,3) * t495;
t396 = qJD(6) * t474 + t431;
t395 = qJD(6) * t476 + t430;
t393 = qJD(1) + (t450 * t508 + t451 * t510) * t530;
t392 = pkin(5) * t537 + pkin(10) * t476 + t492 * t540;
t391 = pkin(5) * t538 + pkin(10) * t474 + t490 * t540;
t390 = rSges(6,1) * t439 + rSges(6,2) * t438 + rSges(6,3) * t476;
t389 = rSges(6,1) * t437 + rSges(6,2) * t436 + rSges(6,3) * t474;
t388 = Icges(6,1) * t439 + Icges(6,4) * t438 + Icges(6,5) * t476;
t387 = Icges(6,1) * t437 + Icges(6,4) * t436 + Icges(6,5) * t474;
t386 = Icges(6,4) * t439 + Icges(6,2) * t438 + Icges(6,6) * t476;
t385 = Icges(6,4) * t437 + Icges(6,2) * t436 + Icges(6,6) * t474;
t384 = Icges(6,5) * t439 + Icges(6,6) * t438 + Icges(6,3) * t476;
t383 = Icges(6,5) * t437 + Icges(6,6) * t436 + Icges(6,3) * t474;
t382 = rSges(7,1) * t435 + rSges(7,2) * t434 + rSges(7,3) * t476;
t381 = rSges(7,1) * t433 + rSges(7,2) * t432 + rSges(7,3) * t474;
t380 = Icges(7,1) * t435 + Icges(7,4) * t434 + Icges(7,5) * t476;
t379 = Icges(7,1) * t433 + Icges(7,4) * t432 + Icges(7,5) * t474;
t378 = Icges(7,4) * t435 + Icges(7,2) * t434 + Icges(7,6) * t476;
t377 = Icges(7,4) * t433 + Icges(7,2) * t432 + Icges(7,6) * t474;
t376 = Icges(7,5) * t435 + Icges(7,6) * t434 + Icges(7,3) * t476;
t375 = Icges(7,5) * t433 + Icges(7,6) * t432 + Icges(7,3) * t474;
t374 = -t415 * t497 + t461 * t480 + t521;
t373 = t416 * t497 - t461 * t479 + t523;
t372 = t415 * t479 - t416 * t480 + t525;
t371 = t460 * t480 + (-t413 - t428) * t497 + t519;
t370 = t414 * t497 + (-t460 - t465) * t479 + t522;
t369 = t413 * t479 + (-t414 - t429) * t480 + t524;
t368 = -t389 * t472 + t420 * t431 + t517;
t367 = t390 * t472 - t420 * t430 + t518;
t366 = t389 * t430 - t390 * t431 + t520;
t365 = -t381 * t443 - t391 * t472 + t396 * t400 + t425 * t431 + t517;
t364 = t382 * t443 + t392 * t472 - t395 * t400 - t425 * t430 + t518;
t363 = t381 * t395 - t382 * t396 + t391 * t430 - t392 * t431 + t520;
t1 = -t546 * ((-t445 * t534 - t447 * t490 + t449 * t491) * t535 - (-t444 * t534 - t446 * t490 + t448 * t491) * t534 + (-t482 * t534 - t483 * t490 + t484 * t491) * t511) * t534 / 0.2e1 + m(7) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + m(6) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(5) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(4) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(3) * (t393 ^ 2 + t421 ^ 2 + t422 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t395 * ((t476 * t376 + t434 * t378 + t435 * t380) * t395 + (t375 * t476 + t377 * t434 + t379 * t435) * t396 + (t397 * t476 + t398 * t434 + t399 * t435) * t443) / 0.2e1 + t443 * ((t376 * t495 + t378 * t466 + t380 * t467) * t395 + (t375 * t495 + t377 * t466 + t379 * t467) * t396 + (t495 * t397 + t466 * t398 + t467 * t399) * t443) / 0.2e1 + t396 * ((t376 * t474 + t378 * t432 + t380 * t433) * t395 + (t474 * t375 + t432 * t377 + t433 * t379) * t396 + (t397 * t474 + t398 * t432 + t399 * t433) * t443) / 0.2e1 + t430 * ((t476 * t384 + t438 * t386 + t439 * t388) * t430 + (t383 * t476 + t385 * t438 + t387 * t439) * t431 + (t417 * t476 + t418 * t438 + t419 * t439) * t472) / 0.2e1 + t431 * ((t384 * t474 + t386 * t436 + t388 * t437) * t430 + (t474 * t383 + t436 * t385 + t437 * t387) * t431 + (t417 * t474 + t418 * t436 + t419 * t437) * t472) / 0.2e1 + t472 * ((t384 * t495 + t386 * t477 + t388 * t478) * t430 + (t383 * t495 + t385 * t477 + t387 * t478) * t431 + (t495 * t417 + t477 * t418 + t478 * t419) * t472) / 0.2e1 + ((t475 * t549 + t476 * t548 + t492 * t547) * t497 + (t475 * t555 + t551 * t476 + t553 * t492) * t480 + (t554 * t475 + t550 * t476 + t552 * t492) * t479) * t479 / 0.2e1 + ((t473 * t549 + t474 * t548 + t490 * t547) * t497 + (t555 * t473 + t551 * t474 + t553 * t490) * t480 + (t473 * t554 + t474 * t550 + t490 * t552) * t479) * t480 / 0.2e1 + ((t549 * t494 + t548 * t495 - t547 * t533) * t497 + (t494 * t555 + t551 * t495 - t553 * t533) * t480 + (t494 * t554 + t495 * t550 - t552 * t533) * t479) * t497 / 0.2e1 + (((t445 * t535 - t447 * t492 + t449 * t493) * t535 - (t444 * t535 - t446 * t492 + t448 * t493) * t534 + (t482 * t535 - t483 * t492 + t484 * t493) * t511) * t535 + t511 * (t511 ^ 2 * t482 + (((t447 * t515 + t449 * t513) * t508 - (t446 * t515 + t448 * t513) * t510) * t509 + (-t444 * t510 + t445 * t508 + t483 * t515 + t484 * t513) * t511) * t509)) * t546 / 0.2e1;
T  = t1;
