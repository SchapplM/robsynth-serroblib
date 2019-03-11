% Calculate kinetic energy for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:38:57
% EndTime: 2019-03-09 11:39:00
% DurationCPUTime: 2.32s
% Computational Cost: add. (1911->263), mult. (1906->397), div. (0->0), fcn. (1807->10), ass. (0->151)
t562 = Icges(6,1) + Icges(7,1);
t561 = Icges(6,4) + Icges(7,4);
t560 = -Icges(7,5) - Icges(6,5);
t559 = Icges(6,2) + Icges(7,2);
t558 = -Icges(7,6) - Icges(6,6);
t557 = Icges(3,3) + Icges(4,3);
t556 = -Icges(7,3) - Icges(6,3);
t459 = qJ(2) + pkin(10);
t451 = sin(t459);
t452 = cos(t459);
t463 = sin(qJ(2));
t466 = cos(qJ(2));
t555 = Icges(3,5) * t466 + Icges(4,5) * t452 - Icges(3,6) * t463 - Icges(4,6) * t451;
t453 = qJ(4) + t459;
t448 = cos(t453);
t465 = cos(qJ(5));
t467 = cos(qJ(1));
t513 = t465 * t467;
t462 = sin(qJ(5));
t464 = sin(qJ(1));
t516 = t462 * t464;
t422 = -t448 * t516 - t513;
t514 = t464 * t465;
t515 = t462 * t467;
t423 = t448 * t514 - t515;
t447 = sin(t453);
t518 = t447 * t464;
t554 = -t558 * t422 - t560 * t423 - t556 * t518;
t424 = -t448 * t515 + t514;
t425 = t448 * t513 + t516;
t517 = t447 * t467;
t553 = -t558 * t424 - t560 * t425 - t556 * t517;
t552 = t559 * t422 + t561 * t423 - t558 * t518;
t551 = t559 * t424 + t561 * t425 - t558 * t517;
t550 = t561 * t422 + t562 * t423 - t560 * t518;
t549 = t561 * t424 + t562 * t425 - t560 * t517;
t548 = t556 * t448 + (t558 * t462 - t560 * t465) * t447;
t547 = t558 * t448 + (-t559 * t462 + t561 * t465) * t447;
t546 = t560 * t448 + (-t561 * t462 + t562 * t465) * t447;
t499 = pkin(4) * t448 + pkin(9) * t447;
t410 = t499 * t464;
t430 = t447 * pkin(4) - t448 * pkin(9);
t439 = (-qJD(2) - qJD(4)) * t467;
t507 = pkin(3) * t452;
t375 = -pkin(8) * t467 + t464 * t507;
t529 = t466 * pkin(2);
t399 = -qJ(3) * t467 + t464 * t529;
t446 = t464 * pkin(1) - t467 * pkin(7);
t508 = -t399 - t446;
t503 = -t375 + t508;
t545 = t439 * t430 + (-t410 + t503) * qJD(1);
t544 = t555 * t464 - t557 * t467;
t411 = t499 * t467;
t456 = qJD(2) * t464;
t438 = qJD(4) * t464 + t456;
t543 = qJD(1) * t411 - t438 * t430;
t542 = t557 * t464 + t555 * t467;
t541 = Icges(3,5) * t463 + Icges(4,5) * t451 + Icges(3,6) * t466 + Icges(4,6) * t452;
t523 = Icges(4,4) * t451;
t432 = Icges(4,2) * t452 + t523;
t522 = Icges(4,4) * t452;
t433 = Icges(4,1) * t451 + t522;
t525 = Icges(3,4) * t463;
t441 = Icges(3,2) * t466 + t525;
t524 = Icges(3,4) * t466;
t442 = Icges(3,1) * t463 + t524;
t540 = -t432 * t451 + t433 * t452 - t441 * t463 + t442 * t466;
t489 = -Icges(4,2) * t451 + t522;
t403 = -Icges(4,6) * t467 + t464 * t489;
t492 = Icges(4,1) * t452 - t523;
t405 = -Icges(4,5) * t467 + t464 * t492;
t490 = -Icges(3,2) * t463 + t524;
t416 = -Icges(3,6) * t467 + t464 * t490;
t493 = Icges(3,1) * t466 - t525;
t418 = -Icges(3,5) * t467 + t464 * t493;
t539 = t403 * t451 - t405 * t452 + t416 * t463 - t418 * t466;
t404 = Icges(4,6) * t464 + t467 * t489;
t406 = Icges(4,5) * t464 + t467 * t492;
t417 = Icges(3,6) * t464 + t467 * t490;
t419 = Icges(3,5) * t464 + t467 * t493;
t538 = -t404 * t451 + t406 * t452 - t417 * t463 + t419 * t466;
t531 = pkin(2) * t463;
t528 = pkin(5) * t465;
t521 = Icges(5,4) * t447;
t520 = Icges(5,4) * t448;
t474 = qJ(6) * t447 + t448 * t528;
t512 = rSges(7,1) * t423 + rSges(7,2) * t422 + rSges(7,3) * t518 - pkin(5) * t515 + t464 * t474;
t511 = rSges(7,1) * t425 + rSges(7,2) * t424 + rSges(7,3) * t517 + pkin(5) * t516 + t467 * t474;
t510 = (-qJ(6) - rSges(7,3)) * t448 + (rSges(7,1) * t465 - rSges(7,2) * t462 + t528) * t447;
t400 = qJ(3) * t464 + t467 * t529;
t505 = qJD(2) * t467;
t509 = t399 * t456 + t400 * t505;
t504 = qJD(5) * t447;
t376 = pkin(8) * t464 + t467 * t507;
t500 = t375 * t456 + t376 * t505 + t509;
t436 = qJD(1) * (t467 * pkin(1) + t464 * pkin(7));
t498 = qJD(1) * t400 - qJD(3) * t467 + t436;
t497 = rSges(3,1) * t466 - rSges(3,2) * t463;
t496 = rSges(4,1) * t452 - rSges(4,2) * t451;
t495 = rSges(5,1) * t448 - rSges(5,2) * t447;
t494 = qJD(2) * (-rSges(4,1) * t451 - rSges(4,2) * t452 - t531);
t491 = Icges(5,1) * t448 - t521;
t488 = -Icges(5,2) * t447 + t520;
t485 = Icges(5,5) * t448 - Icges(5,6) * t447;
t477 = qJD(1) * t376 + t498;
t476 = (-pkin(3) * t451 - t531) * qJD(2);
t475 = t438 * t410 - t439 * t411 + t500;
t473 = qJD(1) * (Icges(5,5) * t447 + Icges(5,6) * t448) + (-Icges(5,3) * t467 + t464 * t485) * t439 + (Icges(5,3) * t464 + t467 * t485) * t438;
t455 = qJD(3) * t464;
t472 = t467 * t476 + t455;
t471 = qJD(6) * t447 + t476;
t470 = t464 * t476 + t477;
t392 = -Icges(5,6) * t467 + t464 * t488;
t393 = Icges(5,6) * t464 + t467 * t488;
t394 = -Icges(5,5) * t467 + t464 * t491;
t395 = Icges(5,5) * t464 + t467 * t491;
t427 = Icges(5,2) * t448 + t521;
t428 = Icges(5,1) * t447 + t520;
t469 = (-t393 * t447 + t395 * t448) * t438 + (-t392 * t447 + t394 * t448) * t439 + (-t427 * t447 + t428 * t448) * qJD(1);
t445 = rSges(2,1) * t467 - rSges(2,2) * t464;
t444 = rSges(2,1) * t464 + rSges(2,2) * t467;
t443 = rSges(3,1) * t463 + rSges(3,2) * t466;
t437 = -qJD(5) * t448 + qJD(1);
t429 = rSges(5,1) * t447 + rSges(5,2) * t448;
t421 = t464 * rSges(3,3) + t467 * t497;
t420 = -t467 * rSges(3,3) + t464 * t497;
t413 = t464 * t504 + t439;
t412 = t467 * t504 + t438;
t408 = t464 * rSges(4,3) + t467 * t496;
t407 = -t467 * rSges(4,3) + t464 * t496;
t398 = t464 * rSges(5,3) + t467 * t495;
t397 = -t467 * rSges(5,3) + t464 * t495;
t386 = -t448 * rSges(6,3) + (rSges(6,1) * t465 - rSges(6,2) * t462) * t447;
t371 = qJD(1) * t421 - t443 * t456 + t436;
t370 = -t443 * t505 + (-t420 - t446) * qJD(1);
t369 = (t420 * t464 + t421 * t467) * qJD(2);
t368 = rSges(6,1) * t425 + rSges(6,2) * t424 + rSges(6,3) * t517;
t366 = rSges(6,1) * t423 + rSges(6,2) * t422 + rSges(6,3) * t518;
t350 = qJD(1) * t408 + t464 * t494 + t498;
t349 = t455 + t467 * t494 + (-t407 + t508) * qJD(1);
t348 = (t407 * t464 + t408 * t467) * qJD(2) + t509;
t347 = qJD(1) * t398 - t438 * t429 + t470;
t346 = t439 * t429 + (-t397 + t503) * qJD(1) + t472;
t345 = t397 * t438 - t398 * t439 + t500;
t344 = t437 * t368 - t412 * t386 + t470 + t543;
t343 = -t437 * t366 + t413 * t386 + t472 + t545;
t342 = t366 * t412 - t368 * t413 + t475;
t341 = -t412 * t510 + t437 * t511 + t464 * t471 + t477 + t543;
t340 = t413 * t510 - t437 * t512 + t467 * t471 + t455 + t545;
t339 = -qJD(6) * t448 + t412 * t512 - t413 * t511 + t475;
t1 = m(4) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(5) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(6) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(7) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(3) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + t438 * (t473 * t464 + t469 * t467) / 0.2e1 + t439 * (t469 * t464 - t473 * t467) / 0.2e1 + ((t547 * t424 + t546 * t425 + t548 * t517) * t437 + (t552 * t424 + t550 * t425 + t554 * t517) * t413 + (t551 * t424 + t549 * t425 + t553 * t517) * t412) * t412 / 0.2e1 + ((t547 * t422 + t546 * t423 + t548 * t518) * t437 + (t552 * t422 + t550 * t423 + t554 * t518) * t413 + (t551 * t422 + t549 * t423 + t553 * t518) * t412) * t413 / 0.2e1 + ((-t553 * t412 - t554 * t413 - t548 * t437) * t448 + ((-t547 * t462 + t546 * t465) * t437 + (-t552 * t462 + t550 * t465) * t413 + (-t551 * t462 + t549 * t465) * t412) * t447) * t437 / 0.2e1 + (Icges(2,3) + m(2) * (t444 ^ 2 + t445 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t542 * t464 ^ 2 + (t539 * t467 + (t538 - t544) * t464) * t467) * qJD(2) + (t541 * t464 + t540 * t467) * qJD(1)) * t456 / 0.2e1 - ((t544 * t467 ^ 2 + (t538 * t464 + (t539 - t542) * t467) * t464) * qJD(2) + (t540 * t464 - t541 * t467) * qJD(1)) * t505 / 0.2e1 + ((t393 * t448 + t395 * t447) * t438 + (t392 * t448 + t394 * t447) * t439 + ((-t403 * t452 - t405 * t451 - t416 * t466 - t418 * t463) * t467 + (t404 * t452 + t406 * t451 + t417 * t466 + t419 * t463) * t464) * qJD(2) + (t448 * t427 + t447 * t428 + t452 * t432 + t451 * t433 + t466 * t441 + t463 * t442) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
