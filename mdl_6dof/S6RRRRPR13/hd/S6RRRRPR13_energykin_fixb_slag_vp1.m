% Calculate kinetic energy for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR13_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:04
% EndTime: 2019-03-09 23:51:07
% DurationCPUTime: 3.12s
% Computational Cost: add. (2975->338), mult. (7485->514), div. (0->0), fcn. (9356->12), ass. (0->152)
t560 = Icges(5,1) + Icges(6,1);
t559 = -Icges(5,4) + Icges(6,5);
t558 = Icges(6,4) + Icges(5,5);
t557 = Icges(5,2) + Icges(6,3);
t556 = Icges(6,2) + Icges(5,3);
t555 = -Icges(5,6) + Icges(6,6);
t513 = sin(qJ(2));
t514 = sin(qJ(1));
t516 = cos(qJ(2));
t517 = cos(qJ(1));
t538 = cos(pkin(6));
t528 = t517 * t538;
t493 = t513 * t514 - t516 * t528;
t494 = t513 * t528 + t514 * t516;
t509 = sin(pkin(6));
t535 = t509 * t517;
t456 = Icges(3,5) * t494 - Icges(3,6) * t493 - Icges(3,3) * t535;
t529 = t514 * t538;
t495 = t517 * t513 + t516 * t529;
t496 = -t513 * t529 + t517 * t516;
t537 = t509 * t514;
t457 = Icges(3,5) * t496 - Icges(3,6) * t495 + Icges(3,3) * t537;
t554 = t509 * (t456 * t517 - t457 * t514);
t512 = sin(qJ(3));
t540 = cos(qJ(3));
t477 = t494 * t540 - t512 * t535;
t511 = sin(qJ(4));
t539 = cos(qJ(4));
t446 = t477 * t511 - t493 * t539;
t447 = t477 * t539 + t493 * t511;
t531 = t509 * t540;
t476 = t494 * t512 + t517 * t531;
t553 = t557 * t446 + t559 * t447 + t555 * t476;
t479 = t496 * t540 + t512 * t537;
t448 = t479 * t511 - t495 * t539;
t449 = t479 * t539 + t495 * t511;
t478 = t496 * t512 - t514 * t531;
t552 = t557 * t448 + t559 * t449 + t555 * t478;
t551 = t555 * t446 + t558 * t447 + t556 * t476;
t550 = t555 * t448 + t558 * t449 + t556 * t478;
t549 = t559 * t446 + t560 * t447 + t558 * t476;
t548 = t559 * t448 + t560 * t449 + t558 * t478;
t492 = t512 * t538 + t513 * t531;
t536 = t509 * t516;
t474 = t492 * t511 + t536 * t539;
t475 = t492 * t539 - t511 * t536;
t491 = t509 * t512 * t513 - t538 * t540;
t547 = t557 * t474 + t559 * t475 + t555 * t491;
t546 = t555 * t474 + t558 * t475 + t556 * t491;
t545 = t559 * t474 + t560 * t475 + t558 * t491;
t468 = pkin(2) * t494 + pkin(9) * t493;
t469 = pkin(2) * t496 + pkin(9) * t495;
t532 = qJD(2) * t509;
t505 = t514 * t532;
t530 = t517 * t532;
t534 = t468 * t505 + t469 * t530;
t480 = qJD(3) * t495 + t505;
t533 = qJD(1) * (pkin(1) * t514 - pkin(8) * t535);
t506 = qJD(2) * t538 + qJD(1);
t442 = qJD(4) * t478 + t480;
t481 = qJD(3) * t493 - t530;
t440 = pkin(3) * t477 + pkin(10) * t476;
t441 = pkin(3) * t479 + pkin(10) * t478;
t526 = t480 * t440 - t441 * t481 + t534;
t443 = qJD(4) * t476 + t481;
t498 = -qJD(3) * t536 + t506;
t497 = (pkin(2) * t513 - pkin(9) * t516) * t509;
t499 = qJD(1) * (pkin(1) * t517 + pkin(8) * t537);
t525 = t506 * t469 - t497 * t505 + t499;
t470 = qJD(4) * t491 + t498;
t409 = pkin(4) * t447 + qJ(5) * t446;
t524 = qJD(5) * t474 + t442 * t409 + t526;
t523 = -t468 * t506 - t497 * t530 - t533;
t467 = pkin(3) * t492 + pkin(10) * t491;
t522 = t498 * t441 - t467 * t480 + t525;
t410 = pkin(4) * t449 + qJ(5) * t448;
t521 = qJD(5) * t446 + t470 * t410 + t522;
t520 = -t440 * t498 + t481 * t467 + t523;
t439 = pkin(4) * t475 + qJ(5) * t474;
t519 = qJD(5) * t448 + t443 * t439 + t520;
t515 = cos(qJ(6));
t510 = sin(qJ(6));
t502 = rSges(2,1) * t517 - rSges(2,2) * t514;
t501 = rSges(2,1) * t514 + rSges(2,2) * t517;
t485 = t538 * rSges(3,3) + (rSges(3,1) * t513 + rSges(3,2) * t516) * t509;
t484 = Icges(3,5) * t538 + (Icges(3,1) * t513 + Icges(3,4) * t516) * t509;
t483 = Icges(3,6) * t538 + (Icges(3,4) * t513 + Icges(3,2) * t516) * t509;
t482 = Icges(3,3) * t538 + (Icges(3,5) * t513 + Icges(3,6) * t516) * t509;
t464 = rSges(3,1) * t496 - rSges(3,2) * t495 + rSges(3,3) * t537;
t463 = rSges(3,1) * t494 - rSges(3,2) * t493 - rSges(3,3) * t535;
t461 = Icges(3,1) * t496 - Icges(3,4) * t495 + Icges(3,5) * t537;
t460 = Icges(3,1) * t494 - Icges(3,4) * t493 - Icges(3,5) * t535;
t459 = Icges(3,4) * t496 - Icges(3,2) * t495 + Icges(3,6) * t537;
t458 = Icges(3,4) * t494 - Icges(3,2) * t493 - Icges(3,6) * t535;
t455 = rSges(4,1) * t492 - rSges(4,2) * t491 - rSges(4,3) * t536;
t454 = Icges(4,1) * t492 - Icges(4,4) * t491 - Icges(4,5) * t536;
t453 = Icges(4,4) * t492 - Icges(4,2) * t491 - Icges(4,6) * t536;
t452 = Icges(4,5) * t492 - Icges(4,6) * t491 - Icges(4,3) * t536;
t451 = pkin(5) * t475 - pkin(11) * t491;
t450 = -qJD(6) * t491 + t470;
t437 = t474 * t510 + t475 * t515;
t436 = t474 * t515 - t475 * t510;
t434 = rSges(4,1) * t479 - rSges(4,2) * t478 + rSges(4,3) * t495;
t433 = rSges(4,1) * t477 - rSges(4,2) * t476 + rSges(4,3) * t493;
t432 = Icges(4,1) * t479 - Icges(4,4) * t478 + Icges(4,5) * t495;
t431 = Icges(4,1) * t477 - Icges(4,4) * t476 + Icges(4,5) * t493;
t430 = Icges(4,4) * t479 - Icges(4,2) * t478 + Icges(4,6) * t495;
t429 = Icges(4,4) * t477 - Icges(4,2) * t476 + Icges(4,6) * t493;
t428 = Icges(4,5) * t479 - Icges(4,6) * t478 + Icges(4,3) * t495;
t427 = Icges(4,5) * t477 - Icges(4,6) * t476 + Icges(4,3) * t493;
t426 = rSges(5,1) * t475 - rSges(5,2) * t474 + rSges(5,3) * t491;
t425 = rSges(6,1) * t475 + rSges(6,2) * t491 + rSges(6,3) * t474;
t418 = pkin(5) * t449 - pkin(11) * t478;
t417 = pkin(5) * t447 - pkin(11) * t476;
t416 = -qJD(6) * t476 + t443;
t415 = -qJD(6) * t478 + t442;
t413 = t464 * t506 - t485 * t505 + t499;
t412 = -t463 * t506 - t485 * t530 - t533;
t411 = (t463 * t514 + t464 * t517) * t532;
t408 = t448 * t510 + t449 * t515;
t407 = t448 * t515 - t449 * t510;
t406 = t446 * t510 + t447 * t515;
t405 = t446 * t515 - t447 * t510;
t403 = rSges(5,1) * t449 - rSges(5,2) * t448 + rSges(5,3) * t478;
t402 = rSges(6,1) * t449 + rSges(6,2) * t478 + rSges(6,3) * t448;
t401 = rSges(5,1) * t447 - rSges(5,2) * t446 + rSges(5,3) * t476;
t400 = rSges(6,1) * t447 + rSges(6,2) * t476 + rSges(6,3) * t446;
t386 = rSges(7,1) * t437 + rSges(7,2) * t436 - rSges(7,3) * t491;
t385 = Icges(7,1) * t437 + Icges(7,4) * t436 - Icges(7,5) * t491;
t384 = Icges(7,4) * t437 + Icges(7,2) * t436 - Icges(7,6) * t491;
t383 = Icges(7,5) * t437 + Icges(7,6) * t436 - Icges(7,3) * t491;
t381 = t434 * t498 - t455 * t480 + t525;
t380 = -t433 * t498 + t455 * t481 + t523;
t379 = rSges(7,1) * t408 + rSges(7,2) * t407 - rSges(7,3) * t478;
t378 = rSges(7,1) * t406 + rSges(7,2) * t405 - rSges(7,3) * t476;
t377 = Icges(7,1) * t408 + Icges(7,4) * t407 - Icges(7,5) * t478;
t376 = Icges(7,1) * t406 + Icges(7,4) * t405 - Icges(7,5) * t476;
t375 = Icges(7,4) * t408 + Icges(7,2) * t407 - Icges(7,6) * t478;
t374 = Icges(7,4) * t406 + Icges(7,2) * t405 - Icges(7,6) * t476;
t373 = Icges(7,5) * t408 + Icges(7,6) * t407 - Icges(7,3) * t478;
t372 = Icges(7,5) * t406 + Icges(7,6) * t405 - Icges(7,3) * t476;
t371 = t433 * t480 - t434 * t481 + t534;
t370 = t403 * t470 - t426 * t442 + t522;
t369 = -t401 * t470 + t426 * t443 + t520;
t368 = t401 * t442 - t403 * t443 + t526;
t367 = t402 * t470 + (-t425 - t439) * t442 + t521;
t366 = t425 * t443 + (-t400 - t409) * t470 + t519;
t365 = t400 * t442 + (-t402 - t410) * t443 + t524;
t364 = t379 * t450 - t386 * t415 + t418 * t470 + (-t439 - t451) * t442 + t521;
t363 = -t378 * t450 + t386 * t416 + t443 * t451 + (-t409 - t417) * t470 + t519;
t362 = t378 * t415 - t379 * t416 + t417 * t442 + (-t410 - t418) * t443 + t524;
t1 = ((t482 * t537 - t483 * t495 + t484 * t496) * t506 + (-(-t458 * t495 + t460 * t496) * t517 + (-t495 * t459 + t496 * t461 - t554) * t514) * t532) * t505 / 0.2e1 - ((-t482 * t535 - t483 * t493 + t484 * t494) * t506 + ((-t459 * t493 + t461 * t494) * t514 + (t493 * t458 - t494 * t460 + t554) * t517) * t532) * t530 / 0.2e1 + t450 * ((-t373 * t491 + t375 * t436 + t377 * t437) * t415 + (-t372 * t491 + t374 * t436 + t376 * t437) * t416 + (-t383 * t491 + t384 * t436 + t385 * t437) * t450) / 0.2e1 + t415 * ((-t478 * t373 + t407 * t375 + t408 * t377) * t415 + (-t372 * t478 + t374 * t407 + t376 * t408) * t416 + (-t383 * t478 + t384 * t407 + t385 * t408) * t450) / 0.2e1 + t506 * ((t538 * t457 + (t459 * t516 + t461 * t513) * t509) * t505 - (t538 * t456 + (t458 * t516 + t460 * t513) * t509) * t530 + (t538 * t482 + (t483 * t516 + t484 * t513) * t509) * t506) / 0.2e1 + m(6) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(5) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t371 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(7) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(3) * (t411 ^ 2 + t412 ^ 2 + t413 ^ 2) / 0.2e1 + t416 * ((-t373 * t476 + t375 * t405 + t377 * t406) * t415 + (-t476 * t372 + t405 * t374 + t406 * t376) * t416 + (-t383 * t476 + t384 * t405 + t385 * t406) * t450) / 0.2e1 + t480 * ((t428 * t495 - t430 * t478 + t432 * t479) * t480 + (t427 * t495 - t429 * t478 + t431 * t479) * t481 + (t452 * t495 - t453 * t478 + t454 * t479) * t498) / 0.2e1 + t481 * ((t428 * t493 - t430 * t476 + t432 * t477) * t480 + (t427 * t493 - t429 * t476 + t431 * t477) * t481 + (t452 * t493 - t453 * t476 + t454 * t477) * t498) / 0.2e1 + t498 * ((-t428 * t536 - t430 * t491 + t432 * t492) * t480 + (-t427 * t536 - t429 * t491 + t431 * t492) * t481 + (-t452 * t536 - t453 * t491 + t454 * t492) * t498) / 0.2e1 + ((t448 * t547 + t449 * t545 + t478 * t546) * t470 + (t448 * t553 + t449 * t549 + t478 * t551) * t443 + (t552 * t448 + t548 * t449 + t550 * t478) * t442) * t442 / 0.2e1 + ((t446 * t547 + t447 * t545 + t476 * t546) * t470 + (t553 * t446 + t549 * t447 + t551 * t476) * t443 + (t446 * t552 + t447 * t548 + t476 * t550) * t442) * t443 / 0.2e1 + ((t547 * t474 + t545 * t475 + t546 * t491) * t470 + (t474 * t553 + t475 * t549 + t491 * t551) * t443 + (t474 * t552 + t475 * t548 + t491 * t550) * t442) * t470 / 0.2e1 + (Icges(2,3) + m(2) * (t501 ^ 2 + t502 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
