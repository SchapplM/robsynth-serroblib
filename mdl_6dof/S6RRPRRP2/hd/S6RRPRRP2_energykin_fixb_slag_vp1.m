% Calculate kinetic energy for
% S6RRPRRP2
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:03
% EndTime: 2019-03-09 11:43:05
% DurationCPUTime: 2.21s
% Computational Cost: add. (1869->255), mult. (1897->390), div. (0->0), fcn. (1808->10), ass. (0->149)
t555 = Icges(6,1) + Icges(7,1);
t554 = -Icges(6,4) + Icges(7,5);
t553 = Icges(7,4) + Icges(6,5);
t552 = Icges(6,2) + Icges(7,3);
t551 = -Icges(7,6) + Icges(6,6);
t550 = Icges(3,3) + Icges(4,3);
t549 = -Icges(6,3) - Icges(7,2);
t458 = qJ(2) + pkin(10);
t450 = sin(t458);
t451 = cos(t458);
t461 = sin(qJ(2));
t464 = cos(qJ(2));
t548 = Icges(3,5) * t464 + Icges(4,5) * t451 - Icges(3,6) * t461 - Icges(4,6) * t450;
t547 = rSges(7,1) + pkin(5);
t546 = rSges(7,3) + qJ(6);
t452 = qJ(4) + t458;
t447 = cos(t452);
t463 = cos(qJ(5));
t465 = cos(qJ(1));
t509 = t463 * t465;
t460 = sin(qJ(5));
t462 = sin(qJ(1));
t512 = t460 * t462;
t421 = t447 * t512 + t509;
t510 = t462 * t463;
t511 = t460 * t465;
t422 = t447 * t510 - t511;
t446 = sin(t452);
t514 = t446 * t462;
t545 = t421 * t552 + t422 * t554 - t514 * t551;
t423 = t447 * t511 - t510;
t424 = t447 * t509 + t512;
t513 = t446 * t465;
t544 = t423 * t552 + t424 * t554 - t513 * t551;
t543 = -t421 * t551 + t422 * t553 - t514 * t549;
t542 = -t423 * t551 + t424 * t553 - t513 * t549;
t541 = t554 * t421 + t422 * t555 + t553 * t514;
t540 = t554 * t423 + t424 * t555 + t553 * t513;
t539 = t551 * t447 + (t460 * t552 + t463 * t554) * t446;
t538 = t549 * t447 + (-t460 * t551 + t463 * t553) * t446;
t537 = -t553 * t447 + (t554 * t460 + t463 * t555) * t446;
t536 = t462 * t548 - t465 * t550;
t535 = t462 * t550 + t465 * t548;
t534 = Icges(3,5) * t461 + Icges(4,5) * t450 + Icges(3,6) * t464 + Icges(4,6) * t451;
t518 = Icges(4,4) * t450;
t431 = Icges(4,2) * t451 + t518;
t517 = Icges(4,4) * t451;
t432 = Icges(4,1) * t450 + t517;
t520 = Icges(3,4) * t461;
t440 = Icges(3,2) * t464 + t520;
t519 = Icges(3,4) * t464;
t441 = Icges(3,1) * t461 + t519;
t533 = -t431 * t450 + t432 * t451 - t440 * t461 + t441 * t464;
t485 = -Icges(4,2) * t450 + t517;
t402 = Icges(4,6) * t462 + t465 * t485;
t488 = Icges(4,1) * t451 - t518;
t404 = Icges(4,5) * t462 + t465 * t488;
t486 = -Icges(3,2) * t461 + t519;
t416 = Icges(3,6) * t462 + t465 * t486;
t489 = Icges(3,1) * t464 - t520;
t418 = Icges(3,5) * t462 + t465 * t489;
t532 = -t402 * t450 + t404 * t451 - t416 * t461 + t418 * t464;
t401 = -Icges(4,6) * t465 + t462 * t485;
t403 = -Icges(4,5) * t465 + t462 * t488;
t415 = -Icges(3,6) * t465 + t462 * t486;
t417 = -Icges(3,5) * t465 + t462 * t489;
t531 = t401 * t450 - t403 * t451 + t415 * t461 - t417 * t464;
t524 = pkin(2) * t461;
t522 = t464 * pkin(2);
t516 = Icges(5,4) * t446;
t515 = Icges(5,4) * t447;
t508 = rSges(7,2) * t514 + t546 * t421 + t547 * t422;
t507 = rSges(7,2) * t513 + t546 * t423 + t547 * t424;
t506 = -rSges(7,2) * t447 + (t546 * t460 + t547 * t463) * t446;
t397 = -qJ(3) * t465 + t462 * t522;
t398 = qJ(3) * t462 + t465 * t522;
t455 = qJD(2) * t462;
t501 = qJD(2) * t465;
t505 = t397 * t455 + t398 * t501;
t445 = pkin(1) * t462 - pkin(7) * t465;
t504 = -t397 - t445;
t503 = pkin(3) * t451;
t437 = qJD(4) * t462 + t455;
t500 = qJD(5) * t446;
t374 = -pkin(8) * t465 + t462 * t503;
t499 = -t374 + t504;
t438 = (-qJD(2) - qJD(4)) * t465;
t375 = pkin(8) * t462 + t465 * t503;
t496 = t374 * t455 + t375 * t501 + t505;
t495 = pkin(4) * t447 + pkin(9) * t446;
t435 = qJD(1) * (pkin(1) * t465 + pkin(7) * t462);
t494 = qJD(1) * t398 - qJD(3) * t465 + t435;
t493 = rSges(3,1) * t464 - rSges(3,2) * t461;
t492 = rSges(4,1) * t451 - rSges(4,2) * t450;
t491 = rSges(5,1) * t447 - rSges(5,2) * t446;
t490 = qJD(2) * (-rSges(4,1) * t450 - rSges(4,2) * t451 - t524);
t487 = Icges(5,1) * t447 - t516;
t484 = -Icges(5,2) * t446 + t515;
t481 = Icges(5,5) * t447 - Icges(5,6) * t446;
t474 = qJD(2) * (-pkin(3) * t450 - t524);
t409 = t495 * t462;
t410 = t495 * t465;
t473 = t437 * t409 - t410 * t438 + t496;
t472 = qJD(1) * (Icges(5,5) * t446 + Icges(5,6) * t447) + (-Icges(5,3) * t465 + t462 * t481) * t438 + (Icges(5,3) * t462 + t465 * t481) * t437;
t454 = qJD(3) * t462;
t471 = t465 * t474 + t454;
t470 = qJD(1) * t375 + t462 * t474 + t494;
t429 = pkin(4) * t446 - pkin(9) * t447;
t469 = t438 * t429 + (-t409 + t499) * qJD(1) + t471;
t468 = qJD(1) * t410 - t429 * t437 + t470;
t390 = -Icges(5,6) * t465 + t462 * t484;
t391 = Icges(5,6) * t462 + t465 * t484;
t392 = -Icges(5,5) * t465 + t462 * t487;
t393 = Icges(5,5) * t462 + t465 * t487;
t426 = Icges(5,2) * t447 + t516;
t427 = Icges(5,1) * t446 + t515;
t467 = (-t391 * t446 + t393 * t447) * t437 + (-t390 * t446 + t392 * t447) * t438 + (-t426 * t446 + t427 * t447) * qJD(1);
t444 = rSges(2,1) * t465 - rSges(2,2) * t462;
t443 = rSges(2,1) * t462 + rSges(2,2) * t465;
t442 = rSges(3,1) * t461 + rSges(3,2) * t464;
t436 = -qJD(5) * t447 + qJD(1);
t428 = rSges(5,1) * t446 + rSges(5,2) * t447;
t420 = rSges(3,3) * t462 + t465 * t493;
t419 = -rSges(3,3) * t465 + t462 * t493;
t412 = t462 * t500 + t438;
t411 = t465 * t500 + t437;
t406 = rSges(4,3) * t462 + t465 * t492;
t405 = -rSges(4,3) * t465 + t462 * t492;
t396 = rSges(5,3) * t462 + t465 * t491;
t395 = -rSges(5,3) * t465 + t462 * t491;
t384 = -rSges(6,3) * t447 + (rSges(6,1) * t463 - rSges(6,2) * t460) * t446;
t368 = qJD(1) * t420 - t442 * t455 + t435;
t367 = -t442 * t501 + (-t419 - t445) * qJD(1);
t366 = (t419 * t462 + t420 * t465) * qJD(2);
t365 = rSges(6,1) * t424 - rSges(6,2) * t423 + rSges(6,3) * t513;
t363 = rSges(6,1) * t422 - rSges(6,2) * t421 + rSges(6,3) * t514;
t349 = qJD(1) * t406 + t462 * t490 + t494;
t348 = t454 + t465 * t490 + (-t405 + t504) * qJD(1);
t347 = (t405 * t462 + t406 * t465) * qJD(2) + t505;
t346 = qJD(1) * t396 - t428 * t437 + t470;
t345 = t428 * t438 + (-t395 + t499) * qJD(1) + t471;
t344 = t395 * t437 - t396 * t438 + t496;
t343 = t365 * t436 - t384 * t411 + t468;
t342 = -t363 * t436 + t384 * t412 + t469;
t341 = t363 * t411 - t365 * t412 + t473;
t340 = qJD(6) * t421 - t411 * t506 + t436 * t507 + t468;
t339 = qJD(6) * t423 + t412 * t506 - t436 * t508 + t469;
t338 = qJD(6) * t446 * t460 + t411 * t508 - t412 * t507 + t473;
t1 = m(7) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(4) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(5) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + t437 * (t472 * t462 + t467 * t465) / 0.2e1 + t438 * (t467 * t462 - t472 * t465) / 0.2e1 + m(6) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + ((t539 * t423 + t537 * t424 + t538 * t513) * t436 + (t545 * t423 + t541 * t424 + t543 * t513) * t412 + (t544 * t423 + t540 * t424 + t542 * t513) * t411) * t411 / 0.2e1 + ((t539 * t421 + t537 * t422 + t538 * t514) * t436 + (t545 * t421 + t541 * t422 + t543 * t514) * t412 + (t544 * t421 + t540 * t422 + t542 * t514) * t411) * t412 / 0.2e1 + ((-t542 * t411 - t543 * t412 - t538 * t436) * t447 + ((t539 * t460 + t537 * t463) * t436 + (t545 * t460 + t541 * t463) * t412 + (t544 * t460 + t540 * t463) * t411) * t446) * t436 / 0.2e1 + (m(2) * (t443 ^ 2 + t444 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t535 * t462 ^ 2 + (t531 * t465 + (t532 - t536) * t462) * t465) * qJD(2) + (t534 * t462 + t533 * t465) * qJD(1)) * t455 / 0.2e1 - ((t536 * t465 ^ 2 + (t532 * t462 + (t531 - t535) * t465) * t462) * qJD(2) + (t533 * t462 - t534 * t465) * qJD(1)) * t501 / 0.2e1 + ((t391 * t447 + t393 * t446) * t437 + (t390 * t447 + t392 * t446) * t438 + ((-t401 * t451 - t403 * t450 - t415 * t464 - t417 * t461) * t465 + (t402 * t451 + t404 * t450 + t416 * t464 + t418 * t461) * t462) * qJD(2) + (t447 * t426 + t446 * t427 + t451 * t431 + t450 * t432 + t464 * t440 + t461 * t441) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
