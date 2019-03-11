% Calculate kinetic energy for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:09
% EndTime: 2019-03-08 21:18:11
% DurationCPUTime: 2.54s
% Computational Cost: add. (2336->328), mult. (5527->479), div. (0->0), fcn. (6711->12), ass. (0->151)
t561 = Icges(5,1) + Icges(4,3);
t560 = -Icges(4,4) - Icges(5,6);
t559 = Icges(5,4) - Icges(4,5);
t558 = Icges(5,5) - Icges(4,6);
t557 = Icges(4,2) + Icges(5,3);
t556 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t506 = sin(pkin(10));
t509 = cos(pkin(10));
t513 = cos(qJ(2));
t510 = cos(pkin(6));
t512 = sin(qJ(2));
t533 = t510 * t512;
t488 = t506 * t513 + t509 * t533;
t507 = sin(pkin(6));
t542 = cos(qJ(3));
t525 = t507 * t542;
t541 = sin(qJ(3));
t472 = t488 * t541 + t509 * t525;
t524 = t507 * t541;
t473 = t488 * t542 - t509 * t524;
t532 = t510 * t513;
t487 = t506 * t512 - t509 * t532;
t555 = t472 * t557 + t473 * t560 + t487 * t558;
t490 = -t506 * t533 + t509 * t513;
t474 = t490 * t541 - t506 * t525;
t475 = t490 * t542 + t506 * t524;
t489 = t506 * t532 + t509 * t512;
t554 = t474 * t557 + t475 * t560 + t489 * t558;
t553 = t472 * t558 - t473 * t559 + t487 * t561;
t552 = t474 * t558 - t475 * t559 + t489 * t561;
t491 = -t510 * t542 + t512 * t524;
t492 = t510 * t541 + t512 * t525;
t534 = t507 * t513;
t551 = t491 * t557 + t492 * t560 - t534 * t558;
t550 = t491 * t558 - t492 * t559 - t534 * t561;
t505 = sin(pkin(11));
t508 = cos(pkin(11));
t434 = t472 * t508 - t487 * t505;
t539 = t472 * t505;
t435 = t487 * t508 + t539;
t549 = Icges(6,5) * t435 + Icges(6,6) * t434 + t472 * t560 + t473 * t556 - t487 * t559;
t436 = t474 * t508 - t489 * t505;
t538 = t474 * t505;
t437 = t489 * t508 + t538;
t548 = Icges(6,5) * t437 + Icges(6,6) * t436 + t474 * t560 + t475 * t556 - t489 * t559;
t470 = t491 * t508 + t505 * t534;
t537 = t491 * t505;
t471 = -t508 * t534 + t537;
t547 = Icges(6,5) * t471 + Icges(6,6) * t470 + t491 * t560 + t492 * t556 + t534 * t559;
t546 = qJD(2) ^ 2;
t540 = pkin(5) * t508;
t536 = t506 * t507;
t535 = t507 * t509;
t426 = pkin(3) * t473 + qJ(4) * t472;
t439 = pkin(4) * t487 + qJ(5) * t473;
t530 = -t426 - t439;
t427 = pkin(3) * t475 + qJ(4) * t474;
t440 = pkin(4) * t489 + qJ(5) * t475;
t529 = -t427 - t440;
t462 = pkin(3) * t492 + qJ(4) * t491;
t478 = -pkin(4) * t534 + t492 * qJ(5);
t528 = -t462 - t478;
t527 = qJD(2) * t507;
t498 = t506 * t527;
t476 = qJD(3) * t489 + t498;
t501 = qJD(2) * t510;
t523 = t509 * t527;
t460 = pkin(2) * t488 + pkin(8) * t487;
t461 = pkin(2) * t490 + pkin(8) * t489;
t522 = t460 * t498 + t461 * t523 + qJD(1);
t477 = qJD(3) * t487 - t523;
t494 = -qJD(3) * t534 + t501;
t521 = qJD(4) * t491 + t426 * t476 + t522;
t493 = (pkin(2) * t512 - pkin(8) * t513) * t507;
t520 = t461 * t501 - t493 * t498;
t519 = qJD(5) * t492 + t439 * t476 + t521;
t518 = qJD(4) * t472 + t427 * t494 + t520;
t517 = (-t460 * t510 - t493 * t535) * qJD(2);
t516 = qJD(5) * t473 + t440 * t494 + t518;
t515 = qJD(4) * t474 + t462 * t477 + t517;
t514 = qJD(5) * t475 + t477 * t478 + t515;
t504 = pkin(11) + qJ(6);
t503 = cos(t504);
t502 = sin(t504);
t482 = t510 * rSges(3,3) + (rSges(3,1) * t512 + rSges(3,2) * t513) * t507;
t481 = Icges(3,5) * t510 + (Icges(3,1) * t512 + Icges(3,4) * t513) * t507;
t480 = Icges(3,6) * t510 + (Icges(3,4) * t512 + Icges(3,2) * t513) * t507;
t479 = Icges(3,3) * t510 + (Icges(3,5) * t512 + Icges(3,6) * t513) * t507;
t469 = qJD(6) * t492 + t494;
t464 = t491 * t502 - t503 * t534;
t463 = t491 * t503 + t502 * t534;
t458 = rSges(4,1) * t492 - rSges(4,2) * t491 - rSges(4,3) * t534;
t457 = -rSges(5,1) * t534 - rSges(5,2) * t492 + rSges(5,3) * t491;
t448 = rSges(3,1) * t490 - rSges(3,2) * t489 + rSges(3,3) * t536;
t447 = rSges(3,1) * t488 - rSges(3,2) * t487 - rSges(3,3) * t535;
t446 = Icges(3,1) * t490 - Icges(3,4) * t489 + Icges(3,5) * t536;
t445 = Icges(3,1) * t488 - Icges(3,4) * t487 - Icges(3,5) * t535;
t444 = Icges(3,4) * t490 - Icges(3,2) * t489 + Icges(3,6) * t536;
t443 = Icges(3,4) * t488 - Icges(3,2) * t487 - Icges(3,6) * t535;
t442 = Icges(3,5) * t490 - Icges(3,6) * t489 + Icges(3,3) * t536;
t441 = Icges(3,5) * t488 - Icges(3,6) * t487 - Icges(3,3) * t535;
t433 = t474 * t502 + t489 * t503;
t432 = t474 * t503 - t489 * t502;
t431 = t472 * t502 + t487 * t503;
t430 = t472 * t503 - t487 * t502;
t429 = qJD(6) * t473 + t477;
t428 = qJD(6) * t475 + t476;
t421 = pkin(5) * t537 + pkin(9) * t492 - t534 * t540;
t420 = (-t447 * t510 - t482 * t535) * qJD(2);
t419 = (t448 * t510 - t482 * t536) * qJD(2);
t418 = rSges(6,1) * t471 + rSges(6,2) * t470 + rSges(6,3) * t492;
t417 = rSges(4,1) * t475 - rSges(4,2) * t474 + rSges(4,3) * t489;
t416 = rSges(4,1) * t473 - rSges(4,2) * t472 + rSges(4,3) * t487;
t415 = rSges(5,1) * t489 - rSges(5,2) * t475 + rSges(5,3) * t474;
t414 = rSges(5,1) * t487 - rSges(5,2) * t473 + rSges(5,3) * t472;
t413 = Icges(6,1) * t471 + Icges(6,4) * t470 + Icges(6,5) * t492;
t412 = Icges(6,4) * t471 + Icges(6,2) * t470 + Icges(6,6) * t492;
t398 = rSges(7,1) * t464 + rSges(7,2) * t463 + rSges(7,3) * t492;
t397 = Icges(7,1) * t464 + Icges(7,4) * t463 + Icges(7,5) * t492;
t396 = Icges(7,4) * t464 + Icges(7,2) * t463 + Icges(7,6) * t492;
t395 = Icges(7,5) * t464 + Icges(7,6) * t463 + Icges(7,3) * t492;
t393 = qJD(1) + (t447 * t506 + t448 * t509) * t527;
t392 = pkin(5) * t538 + pkin(9) * t475 + t489 * t540;
t391 = pkin(5) * t539 + pkin(9) * t473 + t487 * t540;
t390 = rSges(6,1) * t437 + rSges(6,2) * t436 + rSges(6,3) * t475;
t389 = rSges(6,1) * t435 + rSges(6,2) * t434 + rSges(6,3) * t473;
t388 = Icges(6,1) * t437 + Icges(6,4) * t436 + Icges(6,5) * t475;
t387 = Icges(6,1) * t435 + Icges(6,4) * t434 + Icges(6,5) * t473;
t386 = Icges(6,4) * t437 + Icges(6,2) * t436 + Icges(6,6) * t475;
t385 = Icges(6,4) * t435 + Icges(6,2) * t434 + Icges(6,6) * t473;
t382 = rSges(7,1) * t433 + rSges(7,2) * t432 + rSges(7,3) * t475;
t381 = rSges(7,1) * t431 + rSges(7,2) * t430 + rSges(7,3) * t473;
t380 = Icges(7,1) * t433 + Icges(7,4) * t432 + Icges(7,5) * t475;
t379 = Icges(7,1) * t431 + Icges(7,4) * t430 + Icges(7,5) * t473;
t378 = Icges(7,4) * t433 + Icges(7,2) * t432 + Icges(7,6) * t475;
t377 = Icges(7,4) * t431 + Icges(7,2) * t430 + Icges(7,6) * t473;
t376 = Icges(7,5) * t433 + Icges(7,6) * t432 + Icges(7,3) * t475;
t375 = Icges(7,5) * t431 + Icges(7,6) * t430 + Icges(7,3) * t473;
t374 = -t416 * t494 + t458 * t477 + t517;
t373 = t417 * t494 - t458 * t476 + t520;
t372 = t416 * t476 - t417 * t477 + t522;
t371 = t457 * t477 + (-t414 - t426) * t494 + t515;
t370 = t415 * t494 + (-t457 - t462) * t476 + t518;
t369 = t414 * t476 + (-t415 - t427) * t477 + t521;
t368 = t418 * t477 + (-t389 + t530) * t494 + t514;
t367 = t390 * t494 + (-t418 + t528) * t476 + t516;
t366 = t389 * t476 + (-t390 + t529) * t477 + t519;
t365 = -t381 * t469 + t398 * t429 + t421 * t477 + (-t391 + t530) * t494 + t514;
t364 = t382 * t469 + t392 * t494 - t398 * t428 + (-t421 + t528) * t476 + t516;
t363 = t381 * t428 - t382 * t429 + t391 * t476 + (-t392 + t529) * t477 + t519;
t1 = t428 * ((t475 * t376 + t378 * t432 + t433 * t380) * t428 + (t375 * t475 + t377 * t432 + t379 * t433) * t429 + (t395 * t475 + t396 * t432 + t397 * t433) * t469) / 0.2e1 + t429 * ((t376 * t473 + t378 * t430 + t380 * t431) * t428 + (t375 * t473 + t430 * t377 + t431 * t379) * t429 + (t395 * t473 + t396 * t430 + t397 * t431) * t469) / 0.2e1 + t469 * ((t376 * t492 + t378 * t463 + t380 * t464) * t428 + (t375 * t492 + t377 * t463 + t379 * t464) * t429 + (t395 * t492 + t396 * t463 + t397 * t464) * t469) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t393 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 + m(7) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + m(5) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(6) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(4) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 - t546 * ((-t442 * t535 - t444 * t487 + t446 * t488) * t536 - (-t441 * t535 - t443 * t487 + t445 * t488) * t535 + (-t479 * t535 - t480 * t487 + t481 * t488) * t510) * t535 / 0.2e1 + (t510 * (t510 ^ 2 * t479 + (((t444 * t513 + t446 * t512) * t506 - (t443 * t513 + t445 * t512) * t509) * t507 + (-t441 * t509 + t442 * t506 + t480 * t513 + t481 * t512) * t510) * t507) + ((t442 * t536 - t444 * t489 + t446 * t490) * t536 - (t441 * t536 - t443 * t489 + t445 * t490) * t535 + (t479 * t536 - t480 * t489 + t481 * t490) * t510) * t536) * t546 / 0.2e1 + ((t412 * t436 + t413 * t437 + t474 * t551 + t475 * t547 + t489 * t550) * t494 + (t385 * t436 + t387 * t437 + t474 * t555 + t475 * t549 + t489 * t553) * t477 + (t436 * t386 + t437 * t388 + t554 * t474 + t548 * t475 + t552 * t489) * t476) * t476 / 0.2e1 + ((t412 * t434 + t413 * t435 + t472 * t551 + t473 * t547 + t487 * t550) * t494 + (t434 * t385 + t435 * t387 + t555 * t472 + t549 * t473 + t553 * t487) * t477 + (t386 * t434 + t388 * t435 + t472 * t554 + t473 * t548 + t487 * t552) * t476) * t477 / 0.2e1 + ((t470 * t412 + t471 * t413 + t551 * t491 + t547 * t492 - t550 * t534) * t494 + (t385 * t470 + t387 * t471 + t491 * t555 + t492 * t549 - t534 * t553) * t477 + (t386 * t470 + t388 * t471 + t491 * t554 + t492 * t548 - t534 * t552) * t476) * t494 / 0.2e1;
T  = t1;
