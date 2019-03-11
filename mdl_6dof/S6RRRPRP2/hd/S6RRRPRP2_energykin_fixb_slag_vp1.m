% Calculate kinetic energy for
% S6RRRPRP2
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:05
% EndTime: 2019-03-09 16:35:07
% DurationCPUTime: 2.39s
% Computational Cost: add. (1891->257), mult. (1919->387), div. (0->0), fcn. (1830->10), ass. (0->149)
t548 = Icges(6,1) + Icges(7,1);
t547 = -Icges(6,4) + Icges(7,5);
t546 = Icges(7,4) + Icges(6,5);
t545 = Icges(6,2) + Icges(7,3);
t544 = -Icges(7,6) + Icges(6,6);
t543 = Icges(4,3) + Icges(5,3);
t542 = -Icges(6,3) - Icges(7,2);
t458 = qJ(2) + qJ(3);
t450 = pkin(10) + t458;
t446 = sin(t450);
t447 = cos(t450);
t454 = sin(t458);
t455 = cos(t458);
t541 = Icges(4,5) * t455 + Icges(5,5) * t447 - Icges(4,6) * t454 - Icges(5,6) * t446;
t540 = rSges(7,1) + pkin(5);
t539 = rSges(7,3) + qJ(6);
t462 = cos(qJ(5));
t464 = cos(qJ(1));
t506 = t462 * t464;
t459 = sin(qJ(5));
t461 = sin(qJ(1));
t509 = t459 * t461;
t420 = t447 * t509 + t506;
t507 = t461 * t462;
t508 = t459 * t464;
t421 = t447 * t507 - t508;
t511 = t446 * t461;
t538 = t545 * t420 + t547 * t421 - t544 * t511;
t422 = t447 * t508 - t507;
t423 = t447 * t506 + t509;
t510 = t446 * t464;
t537 = t545 * t422 + t547 * t423 - t544 * t510;
t536 = -t544 * t420 + t546 * t421 - t542 * t511;
t535 = -t544 * t422 + t546 * t423 - t542 * t510;
t534 = t547 * t420 + t548 * t421 + t546 * t511;
t533 = t547 * t422 + t548 * t423 + t546 * t510;
t532 = t544 * t447 + (t545 * t459 + t547 * t462) * t446;
t531 = t542 * t447 + (-t544 * t459 + t546 * t462) * t446;
t530 = -t546 * t447 + (t547 * t459 + t548 * t462) * t446;
t512 = Icges(5,4) * t447;
t483 = -Icges(5,2) * t446 + t512;
t389 = -Icges(5,6) * t464 + t461 * t483;
t390 = Icges(5,6) * t461 + t464 * t483;
t513 = Icges(5,4) * t446;
t486 = Icges(5,1) * t447 - t513;
t391 = -Icges(5,5) * t464 + t461 * t486;
t392 = Icges(5,5) * t461 + t464 * t486;
t514 = Icges(4,4) * t455;
t484 = -Icges(4,2) * t454 + t514;
t400 = -Icges(4,6) * t464 + t461 * t484;
t401 = Icges(4,6) * t461 + t464 * t484;
t515 = Icges(4,4) * t454;
t487 = Icges(4,1) * t455 - t515;
t402 = -Icges(4,5) * t464 + t461 * t487;
t403 = Icges(4,5) * t461 + t464 * t487;
t425 = Icges(5,2) * t447 + t513;
t426 = Icges(5,1) * t446 + t512;
t431 = Icges(4,2) * t455 + t515;
t432 = Icges(4,1) * t454 + t514;
t453 = qJD(2) * t461;
t437 = qJD(3) * t461 + t453;
t438 = (-qJD(2) - qJD(3)) * t464;
t529 = (-t389 * t446 + t391 * t447 - t400 * t454 + t402 * t455) * t438 + (-t390 * t446 + t392 * t447 - t401 * t454 + t403 * t455) * t437 + (-t425 * t446 + t426 * t447 - t431 * t454 + t432 * t455) * qJD(1);
t528 = (t541 * t461 - t543 * t464) * t438 + (t543 * t461 + t541 * t464) * t437 + (Icges(4,5) * t454 + Icges(5,5) * t446 + Icges(4,6) * t455 + Icges(5,6) * t447) * qJD(1);
t521 = pkin(3) * t454;
t463 = cos(qJ(2));
t519 = t463 * pkin(2);
t460 = sin(qJ(2));
t517 = Icges(3,4) * t460;
t516 = Icges(3,4) * t463;
t505 = rSges(7,2) * t511 + t539 * t420 + t540 * t421;
t504 = rSges(7,2) * t510 + t539 * t422 + t540 * t423;
t503 = -rSges(7,2) * t447 + (t539 * t459 + t540 * t462) * t446;
t396 = -pkin(8) * t464 + t461 * t519;
t397 = pkin(8) * t461 + t464 * t519;
t498 = qJD(2) * t464;
t502 = t396 * t453 + t397 * t498;
t445 = pkin(1) * t461 - pkin(7) * t464;
t501 = -t396 - t445;
t500 = pkin(3) * t455;
t497 = qJD(5) * t446;
t496 = pkin(2) * qJD(2) * t460;
t373 = -qJ(4) * t464 + t461 * t500;
t495 = t437 * t373 + t502;
t494 = -t373 + t501;
t493 = t464 * t496;
t492 = pkin(4) * t447 + pkin(9) * t446;
t491 = rSges(3,1) * t463 - rSges(3,2) * t460;
t490 = rSges(4,1) * t455 - rSges(4,2) * t454;
t489 = rSges(5,1) * t447 - rSges(5,2) * t446;
t488 = Icges(3,1) * t463 - t517;
t485 = -Icges(3,2) * t460 + t516;
t482 = Icges(3,5) * t463 - Icges(3,6) * t460;
t414 = -Icges(3,6) * t464 + t461 * t485;
t416 = -Icges(3,5) * t464 + t461 * t488;
t479 = t414 * t460 - t416 * t463;
t415 = Icges(3,6) * t461 + t464 * t485;
t417 = Icges(3,5) * t461 + t464 * t488;
t478 = -t415 * t460 + t417 * t463;
t440 = Icges(3,2) * t463 + t517;
t441 = Icges(3,1) * t460 + t516;
t477 = -t440 * t460 + t441 * t463;
t435 = qJD(1) * (pkin(1) * t464 + pkin(7) * t461);
t476 = qJD(1) * t397 - t461 * t496 + t435;
t475 = qJD(4) * t461 + t438 * t521 - t493;
t374 = qJ(4) * t461 + t464 * t500;
t408 = t492 * t461;
t409 = t492 * t464;
t472 = t437 * t408 + (-t374 - t409) * t438 + t495;
t471 = qJD(1) * t374 - qJD(4) * t464 + t476;
t429 = pkin(4) * t446 - pkin(9) * t447;
t470 = t438 * t429 + (-t408 + t494) * qJD(1) + t475;
t469 = qJD(1) * t409 + (-t429 - t521) * t437 + t471;
t444 = rSges(2,1) * t464 - rSges(2,2) * t461;
t443 = rSges(2,1) * t461 + rSges(2,2) * t464;
t442 = rSges(3,1) * t460 + rSges(3,2) * t463;
t439 = Icges(3,5) * t460 + Icges(3,6) * t463;
t436 = -qJD(5) * t447 + qJD(1);
t433 = rSges(4,1) * t454 + rSges(4,2) * t455;
t427 = rSges(5,1) * t446 + rSges(5,2) * t447;
t419 = rSges(3,3) * t461 + t464 * t491;
t418 = -rSges(3,3) * t464 + t461 * t491;
t413 = Icges(3,3) * t461 + t464 * t482;
t412 = -Icges(3,3) * t464 + t461 * t482;
t411 = t461 * t497 + t438;
t410 = t464 * t497 + t437;
t406 = rSges(4,3) * t461 + t464 * t490;
t405 = -rSges(4,3) * t464 + t461 * t490;
t394 = rSges(5,3) * t461 + t464 * t489;
t393 = -rSges(5,3) * t464 + t461 * t489;
t383 = -rSges(6,3) * t447 + (rSges(6,1) * t462 - rSges(6,2) * t459) * t446;
t369 = qJD(1) * t419 - t442 * t453 + t435;
t368 = -t442 * t498 + (-t418 - t445) * qJD(1);
t367 = (t418 * t461 + t419 * t464) * qJD(2);
t365 = rSges(6,1) * t423 - rSges(6,2) * t422 + rSges(6,3) * t510;
t363 = rSges(6,1) * t421 - rSges(6,2) * t420 + rSges(6,3) * t511;
t349 = qJD(1) * t406 - t433 * t437 + t476;
t348 = -t493 + t433 * t438 + (-t405 + t501) * qJD(1);
t347 = t405 * t437 - t406 * t438 + t502;
t346 = qJD(1) * t394 + (-t427 - t521) * t437 + t471;
t345 = t427 * t438 + (-t393 + t494) * qJD(1) + t475;
t344 = t393 * t437 + (-t374 - t394) * t438 + t495;
t343 = t365 * t436 - t383 * t410 + t469;
t342 = -t363 * t436 + t383 * t411 + t470;
t341 = t363 * t410 - t365 * t411 + t472;
t340 = qJD(6) * t420 - t410 * t503 + t436 * t504 + t469;
t339 = qJD(6) * t422 + t411 * t503 - t436 * t505 + t470;
t338 = qJD(6) * t446 * t459 + t410 * t505 - t411 * t504 + t472;
t1 = ((t461 * t439 + t464 * t477) * qJD(1) + (t461 ^ 2 * t413 + (t479 * t464 + (-t412 + t478) * t461) * t464) * qJD(2)) * t453 / 0.2e1 - ((-t464 * t439 + t461 * t477) * qJD(1) + (t464 ^ 2 * t412 + (t478 * t461 + (-t413 + t479) * t464) * t461) * qJD(2)) * t498 / 0.2e1 + m(7) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(5) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(6) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(4) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(3) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + ((t532 * t422 + t530 * t423 + t531 * t510) * t436 + (t538 * t422 + t534 * t423 + t536 * t510) * t411 + (t537 * t422 + t533 * t423 + t535 * t510) * t410) * t410 / 0.2e1 + ((t532 * t420 + t530 * t421 + t531 * t511) * t436 + (t538 * t420 + t534 * t421 + t536 * t511) * t411 + (t537 * t420 + t533 * t421 + t535 * t511) * t410) * t411 / 0.2e1 + ((-t535 * t410 - t536 * t411 - t531 * t436) * t447 + ((t532 * t459 + t530 * t462) * t436 + (t538 * t459 + t534 * t462) * t411 + (t537 * t459 + t533 * t462) * t410) * t446) * t436 / 0.2e1 + (t528 * t461 + t529 * t464) * t437 / 0.2e1 + (t529 * t461 - t528 * t464) * t438 / 0.2e1 + (Icges(2,3) + m(2) * (t443 ^ 2 + t444 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t415 * t463 + t417 * t460) * t461 - (t414 * t463 + t416 * t460) * t464) * qJD(2) + (t389 * t447 + t391 * t446 + t400 * t455 + t402 * t454) * t438 + (t390 * t447 + t392 * t446 + t401 * t455 + t403 * t454) * t437 + (t447 * t425 + t446 * t426 + t455 * t431 + t454 * t432 + t463 * t440 + t460 * t441) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
