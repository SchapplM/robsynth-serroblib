% Calculate kinetic energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:26
% EndTime: 2019-03-09 16:21:29
% DurationCPUTime: 2.57s
% Computational Cost: add. (2405->336), mult. (5579->488), div. (0->0), fcn. (6745->12), ass. (0->155)
t568 = Icges(5,1) + Icges(4,3);
t567 = -Icges(4,4) - Icges(5,6);
t566 = Icges(5,4) - Icges(4,5);
t565 = Icges(5,5) - Icges(4,6);
t564 = Icges(4,2) + Icges(5,3);
t563 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t513 = sin(qJ(2));
t514 = sin(qJ(1));
t515 = cos(qJ(2));
t516 = cos(qJ(1));
t545 = cos(pkin(6));
t527 = t516 * t545;
t488 = t513 * t514 - t515 * t527;
t489 = t513 * t527 + t514 * t515;
t510 = sin(pkin(6));
t539 = t510 * t516;
t448 = Icges(3,5) * t489 - Icges(3,6) * t488 - Icges(3,3) * t539;
t528 = t514 * t545;
t490 = t516 * t513 + t515 * t528;
t491 = -t513 * t528 + t516 * t515;
t541 = t510 * t514;
t449 = Icges(3,5) * t491 - Icges(3,6) * t490 + Icges(3,3) * t541;
t562 = (t448 * t516 - t449 * t514) * t510;
t548 = cos(qJ(3));
t531 = t510 * t548;
t547 = sin(qJ(3));
t471 = t489 * t547 + t516 * t531;
t530 = t510 * t547;
t472 = t489 * t548 - t516 * t530;
t561 = t564 * t471 + t567 * t472 + t565 * t488;
t473 = t491 * t547 - t514 * t531;
t474 = t491 * t548 + t514 * t530;
t560 = t564 * t473 + t567 * t474 + t565 * t490;
t559 = t565 * t471 - t566 * t472 + t568 * t488;
t558 = t565 * t473 - t566 * t474 + t568 * t490;
t486 = t513 * t530 - t545 * t548;
t487 = t513 * t531 + t545 * t547;
t540 = t510 * t515;
t557 = t564 * t486 + t567 * t487 - t565 * t540;
t556 = t565 * t486 - t566 * t487 - t568 * t540;
t509 = sin(pkin(11));
t511 = cos(pkin(11));
t434 = t471 * t511 - t488 * t509;
t544 = t471 * t509;
t435 = t488 * t511 + t544;
t555 = Icges(6,5) * t435 + Icges(6,6) * t434 + t567 * t471 + t563 * t472 - t566 * t488;
t436 = t473 * t511 - t490 * t509;
t543 = t473 * t509;
t437 = t490 * t511 + t543;
t554 = Icges(6,5) * t437 + Icges(6,6) * t436 + t567 * t473 + t563 * t474 - t566 * t490;
t469 = t486 * t511 + t509 * t540;
t542 = t486 * t509;
t470 = -t511 * t540 + t542;
t553 = Icges(6,5) * t470 + Icges(6,6) * t469 + t567 * t486 + t563 * t487 + t566 * t540;
t546 = pkin(5) * t511;
t425 = pkin(3) * t472 + qJ(4) * t471;
t438 = pkin(4) * t488 + qJ(5) * t472;
t537 = -t425 - t438;
t426 = pkin(3) * t474 + qJ(4) * t473;
t439 = pkin(4) * t490 + qJ(5) * t474;
t536 = -t426 - t439;
t460 = pkin(2) * t489 + pkin(9) * t488;
t461 = pkin(2) * t491 + pkin(9) * t490;
t532 = qJD(2) * t510;
t501 = t514 * t532;
t529 = t516 * t532;
t535 = t460 * t501 + t461 * t529;
t459 = pkin(3) * t487 + qJ(4) * t486;
t477 = -pkin(4) * t540 + qJ(5) * t487;
t534 = -t459 - t477;
t475 = qJD(3) * t490 + t501;
t533 = qJD(1) * (pkin(1) * t514 - pkin(8) * t539);
t502 = qJD(2) * t545 + qJD(1);
t526 = qJD(4) * t486 + t475 * t425 + t535;
t476 = qJD(3) * t488 - t529;
t493 = -qJD(3) * t540 + t502;
t524 = qJD(5) * t487 + t475 * t438 + t526;
t492 = (pkin(2) * t513 - pkin(9) * t515) * t510;
t494 = qJD(1) * (pkin(1) * t516 + pkin(8) * t541);
t523 = t502 * t461 - t492 * t501 + t494;
t522 = qJD(4) * t471 + t493 * t426 + t523;
t521 = -t460 * t502 - t492 * t529 - t533;
t520 = qJD(5) * t472 + t493 * t439 + t522;
t519 = qJD(4) * t473 + t476 * t459 + t521;
t518 = qJD(5) * t474 + t476 * t477 + t519;
t508 = pkin(11) + qJ(6);
t507 = cos(t508);
t506 = sin(t508);
t497 = rSges(2,1) * t516 - rSges(2,2) * t514;
t496 = rSges(2,1) * t514 + rSges(2,2) * t516;
t481 = t545 * rSges(3,3) + (rSges(3,1) * t513 + rSges(3,2) * t515) * t510;
t480 = Icges(3,5) * t545 + (Icges(3,1) * t513 + Icges(3,4) * t515) * t510;
t479 = Icges(3,6) * t545 + (Icges(3,4) * t513 + Icges(3,2) * t515) * t510;
t478 = Icges(3,3) * t545 + (Icges(3,5) * t513 + Icges(3,6) * t515) * t510;
t464 = qJD(6) * t487 + t493;
t463 = t486 * t506 - t507 * t540;
t462 = t486 * t507 + t506 * t540;
t456 = rSges(3,1) * t491 - rSges(3,2) * t490 + rSges(3,3) * t541;
t455 = rSges(3,1) * t489 - rSges(3,2) * t488 - rSges(3,3) * t539;
t453 = Icges(3,1) * t491 - Icges(3,4) * t490 + Icges(3,5) * t541;
t452 = Icges(3,1) * t489 - Icges(3,4) * t488 - Icges(3,5) * t539;
t451 = Icges(3,4) * t491 - Icges(3,2) * t490 + Icges(3,6) * t541;
t450 = Icges(3,4) * t489 - Icges(3,2) * t488 - Icges(3,6) * t539;
t447 = rSges(4,1) * t487 - rSges(4,2) * t486 - rSges(4,3) * t540;
t446 = -rSges(5,1) * t540 - rSges(5,2) * t487 + rSges(5,3) * t486;
t432 = t473 * t506 + t490 * t507;
t431 = t473 * t507 - t490 * t506;
t430 = t471 * t506 + t488 * t507;
t429 = t471 * t507 - t488 * t506;
t428 = qJD(6) * t472 + t476;
t427 = qJD(6) * t474 + t475;
t420 = pkin(5) * t542 + pkin(10) * t487 - t540 * t546;
t419 = rSges(4,1) * t474 - rSges(4,2) * t473 + rSges(4,3) * t490;
t418 = rSges(4,1) * t472 - rSges(4,2) * t471 + rSges(4,3) * t488;
t417 = rSges(5,1) * t490 - rSges(5,2) * t474 + rSges(5,3) * t473;
t416 = rSges(5,1) * t488 - rSges(5,2) * t472 + rSges(5,3) * t471;
t403 = rSges(6,1) * t470 + rSges(6,2) * t469 + rSges(6,3) * t487;
t402 = Icges(6,1) * t470 + Icges(6,4) * t469 + Icges(6,5) * t487;
t401 = Icges(6,4) * t470 + Icges(6,2) * t469 + Icges(6,6) * t487;
t398 = rSges(7,1) * t463 + rSges(7,2) * t462 + rSges(7,3) * t487;
t397 = Icges(7,1) * t463 + Icges(7,4) * t462 + Icges(7,5) * t487;
t396 = Icges(7,4) * t463 + Icges(7,2) * t462 + Icges(7,6) * t487;
t395 = Icges(7,5) * t463 + Icges(7,6) * t462 + Icges(7,3) * t487;
t394 = t456 * t502 - t481 * t501 + t494;
t393 = -t455 * t502 - t481 * t529 - t533;
t392 = (t455 * t514 + t456 * t516) * t532;
t391 = pkin(5) * t543 + pkin(10) * t474 + t490 * t546;
t390 = pkin(5) * t544 + pkin(10) * t472 + t488 * t546;
t389 = rSges(6,1) * t437 + rSges(6,2) * t436 + rSges(6,3) * t474;
t388 = rSges(6,1) * t435 + rSges(6,2) * t434 + rSges(6,3) * t472;
t387 = Icges(6,1) * t437 + Icges(6,4) * t436 + Icges(6,5) * t474;
t386 = Icges(6,1) * t435 + Icges(6,4) * t434 + Icges(6,5) * t472;
t385 = Icges(6,4) * t437 + Icges(6,2) * t436 + Icges(6,6) * t474;
t384 = Icges(6,4) * t435 + Icges(6,2) * t434 + Icges(6,6) * t472;
t381 = rSges(7,1) * t432 + rSges(7,2) * t431 + rSges(7,3) * t474;
t380 = rSges(7,1) * t430 + rSges(7,2) * t429 + rSges(7,3) * t472;
t379 = Icges(7,1) * t432 + Icges(7,4) * t431 + Icges(7,5) * t474;
t378 = Icges(7,1) * t430 + Icges(7,4) * t429 + Icges(7,5) * t472;
t377 = Icges(7,4) * t432 + Icges(7,2) * t431 + Icges(7,6) * t474;
t376 = Icges(7,4) * t430 + Icges(7,2) * t429 + Icges(7,6) * t472;
t375 = Icges(7,5) * t432 + Icges(7,6) * t431 + Icges(7,3) * t474;
t374 = Icges(7,5) * t430 + Icges(7,6) * t429 + Icges(7,3) * t472;
t373 = t419 * t493 - t447 * t475 + t523;
t372 = -t418 * t493 + t447 * t476 + t521;
t371 = t418 * t475 - t419 * t476 + t535;
t370 = t417 * t493 + (-t446 - t459) * t475 + t522;
t369 = t446 * t476 + (-t416 - t425) * t493 + t519;
t368 = t416 * t475 + (-t417 - t426) * t476 + t526;
t367 = t389 * t493 + (-t403 + t534) * t475 + t520;
t366 = t403 * t476 + (-t388 + t537) * t493 + t518;
t365 = t388 * t475 + (-t389 + t536) * t476 + t524;
t364 = t381 * t464 + t391 * t493 - t398 * t427 + (-t420 + t534) * t475 + t520;
t363 = -t380 * t464 + t398 * t428 + t420 * t476 + (-t390 + t537) * t493 + t518;
t362 = t380 * t427 - t381 * t428 + t390 * t475 + (-t391 + t536) * t476 + t524;
t1 = t427 * ((t474 * t375 + t431 * t377 + t432 * t379) * t427 + (t374 * t474 + t376 * t431 + t378 * t432) * t428 + (t395 * t474 + t396 * t431 + t397 * t432) * t464) / 0.2e1 + t428 * ((t375 * t472 + t377 * t429 + t379 * t430) * t427 + (t472 * t374 + t429 * t376 + t430 * t378) * t428 + (t395 * t472 + t396 * t429 + t397 * t430) * t464) / 0.2e1 + t502 * ((t545 * t449 + (t451 * t515 + t453 * t513) * t510) * t501 - (t545 * t448 + (t450 * t515 + t452 * t513) * t510) * t529 + (t545 * t478 + (t479 * t515 + t480 * t513) * t510) * t502) / 0.2e1 + m(5) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(6) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(7) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(4) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(3) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + t464 * ((t375 * t487 + t377 * t462 + t379 * t463) * t427 + (t374 * t487 + t376 * t462 + t378 * t463) * t428 + (t487 * t395 + t462 * t396 + t463 * t397) * t464) / 0.2e1 - ((-t478 * t539 - t479 * t488 + t480 * t489) * t502 + ((-t451 * t488 + t453 * t489) * t514 + (t488 * t450 - t489 * t452 + t562) * t516) * t532) * t529 / 0.2e1 + ((t478 * t541 - t479 * t490 + t480 * t491) * t502 + (-(-t450 * t490 + t452 * t491) * t516 + (-t490 * t451 + t491 * t453 - t562) * t514) * t532) * t501 / 0.2e1 + (Icges(2,3) + m(2) * (t496 ^ 2 + t497 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t401 * t436 + t402 * t437 + t473 * t557 + t474 * t553 + t490 * t556) * t493 + (t384 * t436 + t386 * t437 + t473 * t561 + t474 * t555 + t490 * t559) * t476 + (t385 * t436 + t387 * t437 + t560 * t473 + t554 * t474 + t558 * t490) * t475) * t475 / 0.2e1 + ((t401 * t434 + t402 * t435 + t471 * t557 + t472 * t553 + t488 * t556) * t493 + (t384 * t434 + t386 * t435 + t561 * t471 + t555 * t472 + t559 * t488) * t476 + (t385 * t434 + t387 * t435 + t471 * t560 + t472 * t554 + t488 * t558) * t475) * t476 / 0.2e1 + ((t401 * t469 + t402 * t470 + t557 * t486 + t553 * t487 - t556 * t540) * t493 + (t384 * t469 + t386 * t470 + t486 * t561 + t487 * t555 - t540 * t559) * t476 + (t385 * t469 + t387 * t470 + t486 * t560 + t487 * t554 - t540 * t558) * t475) * t493 / 0.2e1;
T  = t1;
