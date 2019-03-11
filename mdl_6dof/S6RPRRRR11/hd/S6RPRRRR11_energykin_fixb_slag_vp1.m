% Calculate kinetic energy for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:16
% EndTime: 2019-03-09 07:38:18
% DurationCPUTime: 2.22s
% Computational Cost: add. (5374->348), mult. (14219->546), div. (0->0), fcn. (18501->16), ass. (0->163)
t552 = cos(pkin(6));
t590 = cos(pkin(13));
t594 = sin(qJ(1));
t576 = t594 * t590;
t550 = sin(pkin(13));
t557 = cos(qJ(1));
t583 = t557 * t550;
t569 = t552 * t576 + t583;
t551 = sin(pkin(6));
t589 = sin(pkin(7));
t579 = t551 * t589;
t591 = cos(pkin(7));
t598 = t569 * t591 - t579 * t594;
t596 = cos(qJ(3));
t595 = cos(qJ(4));
t556 = cos(qJ(5));
t593 = pkin(5) * t556;
t535 = t552 * t583 + t576;
t555 = sin(qJ(3));
t578 = t557 * t590;
t581 = t594 * t550;
t568 = -t552 * t578 + t581;
t562 = t568 * t591;
t577 = t596 * t589;
t584 = t551 * t557;
t514 = t535 * t555 + t562 * t596 + t577 * t584;
t553 = sin(qJ(5));
t588 = t514 * t553;
t536 = -t552 * t581 + t578;
t516 = t536 * t555 + t596 * t598;
t587 = t516 * t553;
t574 = t591 * t590;
t585 = t551 * t550;
t526 = -t551 * t574 * t596 - t552 * t577 + t555 * t585;
t586 = t526 * t553;
t580 = t551 * t591;
t528 = -t557 * t580 + t568 * t589;
t524 = qJD(3) * t528;
t500 = qJD(4) * t514 + t524;
t529 = t569 * t589 + t580 * t594;
t525 = qJD(3) * t529;
t501 = qJD(4) * t516 + t525;
t565 = t552 * t591 - t579 * t590;
t532 = qJD(3) * t565 + qJD(1);
t515 = t535 * t596 + (-t557 * t579 - t562) * t555;
t554 = sin(qJ(4));
t496 = t515 * t554 - t528 * t595;
t460 = qJD(5) * t496 + t500;
t517 = t536 * t596 - t555 * t598;
t498 = t517 * t554 - t529 * t595;
t461 = qJD(5) * t498 + t501;
t582 = t551 * t594;
t518 = qJD(4) * t526 + t532;
t527 = t552 * t589 * t555 + (t550 * t596 + t555 * t574) * t551;
t512 = t527 * t554 - t565 * t595;
t488 = qJD(5) * t512 + t518;
t575 = -qJD(2) * t584 + qJD(1) * (pkin(1) * t557 + qJ(2) * t582);
t538 = pkin(1) * t594 - qJ(2) * t584;
t544 = qJD(2) * t582;
t573 = t544 + (-t535 * pkin(2) - pkin(9) * t528 - t538) * qJD(1);
t571 = qJD(1) * (t536 * pkin(2) + pkin(9) * t529) + t575;
t486 = pkin(3) * t515 + pkin(10) * t514;
t487 = pkin(3) * t517 + pkin(10) * t516;
t546 = qJD(2) * t552;
t570 = t486 * t525 - t487 * t524 + t546;
t508 = pkin(3) * t527 + pkin(10) * t526;
t567 = -t486 * t532 + t508 * t524 + t573;
t497 = t515 * t595 + t528 * t554;
t458 = pkin(4) * t497 + pkin(11) * t496;
t499 = t517 * t595 + t529 * t554;
t459 = pkin(4) * t499 + pkin(11) * t498;
t566 = t458 * t501 - t459 * t500 + t570;
t564 = t487 * t532 - t508 * t525 + t571;
t513 = t527 * t595 + t554 * t565;
t485 = pkin(4) * t513 + pkin(11) * t512;
t561 = -t458 * t518 + t485 * t500 + t567;
t560 = t459 * t518 - t485 * t501 + t564;
t549 = qJ(5) + qJ(6);
t548 = cos(t549);
t547 = sin(t549);
t542 = rSges(2,1) * t557 - rSges(2,2) * t594;
t541 = rSges(2,1) * t594 + rSges(2,2) * t557;
t507 = qJD(1) * (t536 * rSges(3,1) - rSges(3,2) * t569 + rSges(3,3) * t582) + t575;
t506 = t544 + (-t535 * rSges(3,1) + rSges(3,2) * t568 + rSges(3,3) * t584 - t538) * qJD(1);
t505 = rSges(4,1) * t527 - rSges(4,2) * t526 + rSges(4,3) * t565;
t504 = Icges(4,1) * t527 - Icges(4,4) * t526 + Icges(4,5) * t565;
t503 = Icges(4,4) * t527 - Icges(4,2) * t526 + Icges(4,6) * t565;
t502 = Icges(4,5) * t527 - Icges(4,6) * t526 + Icges(4,3) * t565;
t493 = t513 * t556 + t586;
t492 = -t513 * t553 + t526 * t556;
t491 = t513 * t548 + t526 * t547;
t490 = -t513 * t547 + t526 * t548;
t482 = rSges(4,1) * t517 - rSges(4,2) * t516 + rSges(4,3) * t529;
t481 = rSges(4,1) * t515 - rSges(4,2) * t514 + rSges(4,3) * t528;
t480 = Icges(4,1) * t517 - Icges(4,4) * t516 + Icges(4,5) * t529;
t479 = Icges(4,1) * t515 - Icges(4,4) * t514 + Icges(4,5) * t528;
t478 = Icges(4,4) * t517 - Icges(4,2) * t516 + Icges(4,6) * t529;
t477 = Icges(4,4) * t515 - Icges(4,2) * t514 + Icges(4,6) * t528;
t476 = Icges(4,5) * t517 - Icges(4,6) * t516 + Icges(4,3) * t529;
t475 = Icges(4,5) * t515 - Icges(4,6) * t514 + Icges(4,3) * t528;
t474 = rSges(5,1) * t513 - rSges(5,2) * t512 + rSges(5,3) * t526;
t473 = Icges(5,1) * t513 - Icges(5,4) * t512 + Icges(5,5) * t526;
t472 = Icges(5,4) * t513 - Icges(5,2) * t512 + Icges(5,6) * t526;
t471 = Icges(5,5) * t513 - Icges(5,6) * t512 + Icges(5,3) * t526;
t470 = t499 * t556 + t587;
t469 = -t499 * t553 + t516 * t556;
t468 = t497 * t556 + t588;
t467 = -t497 * t553 + t514 * t556;
t466 = t499 * t548 + t516 * t547;
t465 = -t499 * t547 + t516 * t548;
t464 = t497 * t548 + t514 * t547;
t463 = -t497 * t547 + t514 * t548;
t462 = qJD(6) * t512 + t488;
t455 = rSges(5,1) * t499 - rSges(5,2) * t498 + rSges(5,3) * t516;
t454 = rSges(5,1) * t497 - rSges(5,2) * t496 + rSges(5,3) * t514;
t453 = Icges(5,1) * t499 - Icges(5,4) * t498 + Icges(5,5) * t516;
t452 = Icges(5,1) * t497 - Icges(5,4) * t496 + Icges(5,5) * t514;
t451 = Icges(5,4) * t499 - Icges(5,2) * t498 + Icges(5,6) * t516;
t450 = Icges(5,4) * t497 - Icges(5,2) * t496 + Icges(5,6) * t514;
t449 = Icges(5,5) * t499 - Icges(5,6) * t498 + Icges(5,3) * t516;
t448 = Icges(5,5) * t497 - Icges(5,6) * t496 + Icges(5,3) * t514;
t447 = rSges(6,1) * t493 + rSges(6,2) * t492 + rSges(6,3) * t512;
t446 = Icges(6,1) * t493 + Icges(6,4) * t492 + Icges(6,5) * t512;
t445 = Icges(6,4) * t493 + Icges(6,2) * t492 + Icges(6,6) * t512;
t444 = Icges(6,5) * t493 + Icges(6,6) * t492 + Icges(6,3) * t512;
t443 = rSges(7,1) * t491 + rSges(7,2) * t490 + rSges(7,3) * t512;
t442 = Icges(7,1) * t491 + Icges(7,4) * t490 + Icges(7,5) * t512;
t441 = Icges(7,4) * t491 + Icges(7,2) * t490 + Icges(7,6) * t512;
t440 = Icges(7,5) * t491 + Icges(7,6) * t490 + Icges(7,3) * t512;
t439 = pkin(5) * t586 + pkin(12) * t512 + t513 * t593;
t438 = qJD(6) * t498 + t461;
t437 = qJD(6) * t496 + t460;
t435 = t482 * t532 - t505 * t525 + t571;
t434 = -t481 * t532 + t505 * t524 + t573;
t433 = t546 + (t481 * t529 - t482 * t528) * qJD(3);
t432 = rSges(6,1) * t470 + rSges(6,2) * t469 + rSges(6,3) * t498;
t431 = rSges(6,1) * t468 + rSges(6,2) * t467 + rSges(6,3) * t496;
t430 = Icges(6,1) * t470 + Icges(6,4) * t469 + Icges(6,5) * t498;
t429 = Icges(6,1) * t468 + Icges(6,4) * t467 + Icges(6,5) * t496;
t428 = Icges(6,4) * t470 + Icges(6,2) * t469 + Icges(6,6) * t498;
t427 = Icges(6,4) * t468 + Icges(6,2) * t467 + Icges(6,6) * t496;
t426 = Icges(6,5) * t470 + Icges(6,6) * t469 + Icges(6,3) * t498;
t425 = Icges(6,5) * t468 + Icges(6,6) * t467 + Icges(6,3) * t496;
t424 = rSges(7,1) * t466 + rSges(7,2) * t465 + rSges(7,3) * t498;
t423 = rSges(7,1) * t464 + rSges(7,2) * t463 + rSges(7,3) * t496;
t422 = Icges(7,1) * t466 + Icges(7,4) * t465 + Icges(7,5) * t498;
t421 = Icges(7,1) * t464 + Icges(7,4) * t463 + Icges(7,5) * t496;
t420 = Icges(7,4) * t466 + Icges(7,2) * t465 + Icges(7,6) * t498;
t419 = Icges(7,4) * t464 + Icges(7,2) * t463 + Icges(7,6) * t496;
t418 = Icges(7,5) * t466 + Icges(7,6) * t465 + Icges(7,3) * t498;
t417 = Icges(7,5) * t464 + Icges(7,6) * t463 + Icges(7,3) * t496;
t416 = pkin(5) * t587 + pkin(12) * t498 + t499 * t593;
t415 = pkin(5) * t588 + pkin(12) * t496 + t497 * t593;
t414 = t455 * t518 - t474 * t501 + t564;
t413 = -t454 * t518 + t474 * t500 + t567;
t412 = t454 * t501 - t455 * t500 + t570;
t411 = t432 * t488 - t447 * t461 + t560;
t410 = -t431 * t488 + t447 * t460 + t561;
t409 = t431 * t461 - t432 * t460 + t566;
t408 = t416 * t488 + t424 * t462 - t438 * t443 - t439 * t461 + t560;
t407 = -t415 * t488 - t423 * t462 + t437 * t443 + t439 * t460 + t561;
t406 = t415 * t461 - t416 * t460 + t423 * t438 - t424 * t437 + t566;
t1 = ((t502 * t529 - t503 * t516 + t504 * t517) * t532 + ((t476 * t529 - t478 * t516 + t480 * t517) * t529 + (t475 * t529 - t477 * t516 + t479 * t517) * t528) * qJD(3)) * t525 / 0.2e1 + m(7) * (t406 ^ 2 + t407 ^ 2 + t408 ^ 2) / 0.2e1 + m(6) * (t409 ^ 2 + t410 ^ 2 + t411 ^ 2) / 0.2e1 + m(5) * (t412 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(4) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t552 ^ 2 + t506 ^ 2 + t507 ^ 2) / 0.2e1 + ((t502 * t528 - t503 * t514 + t504 * t515) * t532 + ((t476 * t528 - t478 * t514 + t480 * t515) * t529 + (t475 * t528 - t477 * t514 + t479 * t515) * t528) * qJD(3)) * t524 / 0.2e1 + t532 * ((t502 * t565 - t526 * t503 + t527 * t504) * t532 + ((t476 * t565 - t526 * t478 + t527 * t480) * t529 + (t475 * t565 - t526 * t477 + t527 * t479) * t528) * qJD(3)) / 0.2e1 + t437 * ((t418 * t496 + t420 * t463 + t422 * t464) * t438 + (t496 * t417 + t463 * t419 + t464 * t421) * t437 + (t440 * t496 + t441 * t463 + t442 * t464) * t462) / 0.2e1 + t462 * ((t418 * t512 + t420 * t490 + t422 * t491) * t438 + (t417 * t512 + t419 * t490 + t421 * t491) * t437 + (t512 * t440 + t490 * t441 + t491 * t442) * t462) / 0.2e1 + t438 * ((t498 * t418 + t465 * t420 + t466 * t422) * t438 + (t417 * t498 + t419 * t465 + t421 * t466) * t437 + (t440 * t498 + t441 * t465 + t442 * t466) * t462) / 0.2e1 + t461 * ((t498 * t426 + t469 * t428 + t470 * t430) * t461 + (t425 * t498 + t427 * t469 + t429 * t470) * t460 + (t444 * t498 + t445 * t469 + t446 * t470) * t488) / 0.2e1 + t460 * ((t426 * t496 + t428 * t467 + t430 * t468) * t461 + (t496 * t425 + t467 * t427 + t468 * t429) * t460 + (t444 * t496 + t445 * t467 + t446 * t468) * t488) / 0.2e1 + t488 * ((t426 * t512 + t428 * t492 + t430 * t493) * t461 + (t425 * t512 + t427 * t492 + t429 * t493) * t460 + (t512 * t444 + t492 * t445 + t493 * t446) * t488) / 0.2e1 + t500 * ((t449 * t514 - t451 * t496 + t453 * t497) * t501 + (t514 * t448 - t496 * t450 + t497 * t452) * t500 + (t471 * t514 - t472 * t496 + t473 * t497) * t518) / 0.2e1 + t518 * ((t449 * t526 - t451 * t512 + t453 * t513) * t501 + (t448 * t526 - t450 * t512 + t452 * t513) * t500 + (t526 * t471 - t512 * t472 + t513 * t473) * t518) / 0.2e1 + t501 * ((t516 * t449 - t498 * t451 + t499 * t453) * t501 + (t448 * t516 - t450 * t498 + t452 * t499) * t500 + (t471 * t516 - t472 * t498 + t473 * t499) * t518) / 0.2e1 + (m(2) * (t541 ^ 2 + t542 ^ 2) + Icges(2,3) + (Icges(3,5) * t552 + (Icges(3,1) * t550 + Icges(3,4) * t590) * t551) * t585 + t551 * t590 * (Icges(3,6) * t552 + (Icges(3,4) * t550 + Icges(3,2) * t590) * t551) + t552 * (Icges(3,3) * t552 + (Icges(3,5) * t550 + Icges(3,6) * t590) * t551)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
