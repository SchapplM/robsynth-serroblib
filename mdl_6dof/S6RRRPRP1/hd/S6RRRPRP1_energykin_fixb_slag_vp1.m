% Calculate kinetic energy for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:41
% EndTime: 2019-03-09 16:30:44
% DurationCPUTime: 2.47s
% Computational Cost: add. (1933->266), mult. (1928->394), div. (0->0), fcn. (1829->10), ass. (0->152)
t555 = Icges(6,1) + Icges(7,1);
t554 = Icges(6,4) + Icges(7,4);
t553 = -Icges(7,5) - Icges(6,5);
t552 = Icges(6,2) + Icges(7,2);
t551 = -Icges(7,6) - Icges(6,6);
t550 = Icges(4,3) + Icges(5,3);
t549 = -Icges(7,3) - Icges(6,3);
t459 = qJ(2) + qJ(3);
t451 = pkin(10) + t459;
t447 = sin(t451);
t448 = cos(t451);
t455 = sin(t459);
t456 = cos(t459);
t548 = Icges(4,5) * t456 + Icges(5,5) * t448 - Icges(4,6) * t455 - Icges(5,6) * t447;
t464 = cos(qJ(5));
t466 = cos(qJ(1));
t512 = t464 * t466;
t461 = sin(qJ(5));
t463 = sin(qJ(1));
t515 = t461 * t463;
t421 = -t448 * t515 - t512;
t513 = t463 * t464;
t514 = t461 * t466;
t422 = t448 * t513 - t514;
t517 = t447 * t463;
t547 = -t551 * t421 - t553 * t422 - t549 * t517;
t423 = -t448 * t514 + t513;
t424 = t448 * t512 + t515;
t516 = t447 * t466;
t546 = -t551 * t423 - t553 * t424 - t549 * t516;
t545 = t552 * t421 + t554 * t422 - t551 * t517;
t544 = t552 * t423 + t554 * t424 - t551 * t516;
t543 = t554 * t421 + t555 * t422 - t553 * t517;
t542 = t554 * t423 + t555 * t424 - t553 * t516;
t504 = pkin(3) * t456;
t376 = qJ(4) * t463 + t466 * t504;
t541 = qJD(1) * t376 - qJD(4) * t466;
t540 = t549 * t448 + (t551 * t461 - t553 * t464) * t447;
t539 = t551 * t448 + (-t552 * t461 + t554 * t464) * t447;
t538 = t553 * t448 + (-t554 * t461 + t555 * t464) * t447;
t518 = Icges(5,4) * t448;
t485 = -Icges(5,2) * t447 + t518;
t391 = -Icges(5,6) * t466 + t463 * t485;
t392 = Icges(5,6) * t463 + t466 * t485;
t519 = Icges(5,4) * t447;
t488 = Icges(5,1) * t448 - t519;
t393 = -Icges(5,5) * t466 + t463 * t488;
t394 = Icges(5,5) * t463 + t466 * t488;
t520 = Icges(4,4) * t456;
t486 = -Icges(4,2) * t455 + t520;
t402 = -Icges(4,6) * t466 + t463 * t486;
t403 = Icges(4,6) * t463 + t466 * t486;
t521 = Icges(4,4) * t455;
t489 = Icges(4,1) * t456 - t521;
t404 = -Icges(4,5) * t466 + t463 * t489;
t405 = Icges(4,5) * t463 + t466 * t489;
t426 = Icges(5,2) * t448 + t519;
t427 = Icges(5,1) * t447 + t518;
t432 = Icges(4,2) * t456 + t521;
t433 = Icges(4,1) * t455 + t520;
t454 = qJD(2) * t463;
t438 = qJD(3) * t463 + t454;
t439 = (-qJD(2) - qJD(3)) * t466;
t537 = (-t391 * t447 + t393 * t448 - t402 * t455 + t404 * t456) * t439 + (-t392 * t447 + t394 * t448 - t403 * t455 + t405 * t456) * t438 + (-t426 * t447 + t427 * t448 - t432 * t455 + t433 * t456) * qJD(1);
t536 = (t548 * t463 - t550 * t466) * t439 + (t550 * t463 + t548 * t466) * t438 + (Icges(4,5) * t455 + Icges(5,5) * t447 + Icges(4,6) * t456 + Icges(5,6) * t448) * qJD(1);
t529 = pkin(3) * t455;
t465 = cos(qJ(2));
t527 = t465 * pkin(2);
t526 = pkin(5) * t464;
t462 = sin(qJ(2));
t523 = Icges(3,4) * t462;
t522 = Icges(3,4) * t465;
t476 = qJ(6) * t447 + t448 * t526;
t511 = rSges(7,1) * t422 + rSges(7,2) * t421 + rSges(7,3) * t517 - pkin(5) * t514 + t463 * t476;
t510 = rSges(7,1) * t424 + rSges(7,2) * t423 + rSges(7,3) * t516 + pkin(5) * t515 + t466 * t476;
t509 = (-qJ(6) - rSges(7,3)) * t448 + (rSges(7,1) * t464 - rSges(7,2) * t461 + t526) * t447;
t398 = -pkin(8) * t466 + t463 * t527;
t399 = pkin(8) * t463 + t466 * t527;
t502 = qJD(2) * t466;
t508 = t398 * t454 + t399 * t502;
t436 = qJD(1) * (pkin(1) * t466 + pkin(7) * t463);
t507 = qJD(1) * t399 + t436;
t446 = pkin(1) * t463 - pkin(7) * t466;
t506 = -t398 - t446;
t505 = qJD(4) * t463 + t439 * t529;
t500 = qJD(5) * t447;
t499 = pkin(2) * qJD(2) * t462;
t375 = -qJ(4) * t466 + t463 * t504;
t498 = t438 * t375 + t508;
t497 = -t375 + t506;
t496 = t463 * t499;
t495 = t466 * t499;
t494 = pkin(4) * t448 + pkin(9) * t447;
t493 = rSges(3,1) * t465 - rSges(3,2) * t462;
t492 = rSges(4,1) * t456 - rSges(4,2) * t455;
t491 = rSges(5,1) * t448 - rSges(5,2) * t447;
t490 = Icges(3,1) * t465 - t523;
t487 = -Icges(3,2) * t462 + t522;
t484 = Icges(3,5) * t465 - Icges(3,6) * t462;
t415 = -Icges(3,6) * t466 + t463 * t487;
t417 = -Icges(3,5) * t466 + t463 * t490;
t481 = t415 * t462 - t417 * t465;
t416 = Icges(3,6) * t463 + t466 * t487;
t418 = Icges(3,5) * t463 + t466 * t490;
t480 = -t416 * t462 + t418 * t465;
t441 = Icges(3,2) * t465 + t523;
t442 = Icges(3,1) * t462 + t522;
t479 = -t441 * t462 + t442 * t465;
t478 = qJD(6) * t447 - t499;
t477 = -t496 + t507;
t409 = t494 * t463;
t410 = t494 * t466;
t473 = t438 * t409 + (-t376 - t410) * t439 + t498;
t430 = pkin(4) * t447 - pkin(9) * t448;
t472 = t439 * t430 + (-t409 + t497) * qJD(1) + t505;
t471 = qJD(1) * t410 + (-t430 - t529) * t438 + t507 + t541;
t445 = rSges(2,1) * t466 - rSges(2,2) * t463;
t444 = rSges(2,1) * t463 + rSges(2,2) * t466;
t443 = rSges(3,1) * t462 + rSges(3,2) * t465;
t440 = Icges(3,5) * t462 + Icges(3,6) * t465;
t437 = -qJD(5) * t448 + qJD(1);
t434 = rSges(4,1) * t455 + rSges(4,2) * t456;
t428 = rSges(5,1) * t447 + rSges(5,2) * t448;
t420 = rSges(3,3) * t463 + t466 * t493;
t419 = -rSges(3,3) * t466 + t463 * t493;
t414 = Icges(3,3) * t463 + t466 * t484;
t413 = -Icges(3,3) * t466 + t463 * t484;
t412 = t463 * t500 + t439;
t411 = t466 * t500 + t438;
t408 = rSges(4,3) * t463 + t466 * t492;
t407 = -rSges(4,3) * t466 + t463 * t492;
t396 = rSges(5,3) * t463 + t466 * t491;
t395 = -rSges(5,3) * t466 + t463 * t491;
t385 = -rSges(6,3) * t448 + (rSges(6,1) * t464 - rSges(6,2) * t461) * t447;
t372 = qJD(1) * t420 - t443 * t454 + t436;
t371 = -t443 * t502 + (-t419 - t446) * qJD(1);
t370 = (t419 * t463 + t420 * t466) * qJD(2);
t368 = rSges(6,1) * t424 + rSges(6,2) * t423 + rSges(6,3) * t516;
t366 = rSges(6,1) * t422 + rSges(6,2) * t421 + rSges(6,3) * t517;
t350 = qJD(1) * t408 - t434 * t438 + t477;
t349 = -t495 + t434 * t439 + (-t407 + t506) * qJD(1);
t348 = t407 * t438 - t408 * t439 + t508;
t347 = qJD(1) * t396 + (-t428 - t529) * t438 + t477 + t541;
t346 = -t495 + t428 * t439 + (-t395 + t497) * qJD(1) + t505;
t345 = t395 * t438 + (-t376 - t396) * t439 + t498;
t344 = t368 * t437 - t385 * t411 + t471 - t496;
t343 = -t366 * t437 + t385 * t412 + t472 - t495;
t342 = t366 * t411 - t368 * t412 + t473;
t341 = -t411 * t509 + t437 * t510 + t463 * t478 + t471;
t340 = t412 * t509 - t437 * t511 + t466 * t478 + t472;
t339 = -qJD(6) * t448 + t411 * t511 - t412 * t510 + t473;
t1 = m(7) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(5) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(6) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(4) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(3) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + ((t463 * t440 + t466 * t479) * qJD(1) + (t463 ^ 2 * t414 + (t481 * t466 + (-t413 + t480) * t463) * t466) * qJD(2)) * t454 / 0.2e1 - ((-t466 * t440 + t463 * t479) * qJD(1) + (t466 ^ 2 * t413 + (t480 * t463 + (-t414 + t481) * t466) * t463) * qJD(2)) * t502 / 0.2e1 + ((t423 * t539 + t424 * t538 + t516 * t540) * t437 + (t545 * t423 + t543 * t424 + t516 * t547) * t412 + (t423 * t544 + t424 * t542 + t516 * t546) * t411) * t411 / 0.2e1 + ((t421 * t539 + t422 * t538 + t517 * t540) * t437 + (t545 * t421 + t543 * t422 + t517 * t547) * t412 + (t421 * t544 + t422 * t542 + t546 * t517) * t411) * t412 / 0.2e1 + ((-t546 * t411 - t412 * t547 - t540 * t437) * t448 + ((-t461 * t539 + t464 * t538) * t437 + (-t461 * t545 + t464 * t543) * t412 + (-t461 * t544 + t464 * t542) * t411) * t447) * t437 / 0.2e1 + (t463 * t536 + t466 * t537) * t438 / 0.2e1 + (t463 * t537 - t466 * t536) * t439 / 0.2e1 + (m(2) * (t444 ^ 2 + t445 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t416 * t465 + t418 * t462) * t463 - (t415 * t465 + t417 * t462) * t466) * qJD(2) + (t391 * t448 + t393 * t447 + t402 * t456 + t404 * t455) * t439 + (t392 * t448 + t394 * t447 + t403 * t456 + t405 * t455) * t438 + (t448 * t426 + t447 * t427 + t456 * t432 + t455 * t433 + t465 * t441 + t462 * t442) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
