% Calculate kinetic energy for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:50
% EndTime: 2019-03-09 09:22:53
% DurationCPUTime: 2.71s
% Computational Cost: add. (1453->321), mult. (3121->484), div. (0->0), fcn. (3464->10), ass. (0->157)
t535 = Icges(4,1) + Icges(5,1);
t534 = Icges(4,4) - Icges(5,5);
t533 = Icges(5,4) + Icges(4,5);
t532 = Icges(4,2) + Icges(5,3);
t531 = -Icges(5,6) + Icges(4,6);
t530 = -Icges(4,3) - Icges(5,2);
t469 = sin(pkin(10));
t470 = cos(pkin(10));
t476 = cos(qJ(1));
t473 = sin(qJ(1));
t475 = cos(qJ(2));
t509 = t475 * t473;
t438 = t469 * t509 + t470 * t476;
t439 = -t469 * t476 + t470 * t509;
t472 = sin(qJ(2));
t511 = t472 * t473;
t529 = -t438 * t532 + t439 * t534 + t511 * t531;
t508 = t475 * t476;
t440 = t469 * t508 - t470 * t473;
t441 = t469 * t473 + t470 * t508;
t510 = t472 * t476;
t528 = t440 * t532 - t441 * t534 - t510 * t531;
t527 = t438 * t531 - t439 * t533 + t511 * t530;
t526 = -t440 * t531 + t441 * t533 - t510 * t530;
t525 = t438 * t534 - t439 * t535 - t511 * t533;
t524 = -t440 * t534 + t441 * t535 + t510 * t533;
t523 = t531 * t475 + (t469 * t532 - t470 * t534) * t472;
t522 = t530 * t475 + (-t469 * t531 + t470 * t533) * t472;
t521 = -t533 * t475 + (-t469 * t534 + t470 * t535) * t472;
t474 = cos(qJ(5));
t518 = pkin(5) * t474;
t516 = Icges(3,4) * t472;
t515 = Icges(3,4) * t475;
t471 = sin(qJ(5));
t514 = t438 * t471;
t513 = t440 * t471;
t512 = t469 * t471;
t503 = qJD(3) * t472;
t460 = t476 * t503;
t507 = qJD(4) * t440 + t460;
t453 = pkin(2) * t472 - qJ(3) * t475;
t506 = -(pkin(3) * t470 + qJ(4) * t469) * t472 - t453;
t490 = pkin(2) * t475 + qJ(3) * t472;
t443 = t490 * t473;
t457 = pkin(1) * t473 - pkin(7) * t476;
t505 = -t443 - t457;
t465 = qJD(2) * t473;
t504 = qJD(2) * t476;
t502 = qJD(5) * t472;
t461 = qJD(5) * t475 + qJD(1);
t404 = pkin(3) * t439 + qJ(4) * t438;
t501 = -t404 + t505;
t444 = t490 * t476;
t449 = qJD(1) * (pkin(1) * t476 + pkin(7) * t473);
t500 = qJD(1) * t444 + t473 * t503 + t449;
t499 = t472 * pkin(9);
t496 = qJD(2) * (rSges(4,3) * t475 - (rSges(4,1) * t470 - rSges(4,2) * t469) * t472 - t453);
t495 = t472 * (-qJD(5) - qJD(6));
t494 = qJD(2) * (rSges(5,2) * t475 - (rSges(5,1) * t470 + rSges(5,3) * t469) * t472 + t506);
t493 = qJD(2) * (-pkin(4) * t470 * t472 - pkin(8) * t475 + t506);
t492 = -qJD(3) * t475 + t443 * t465 + t444 * t504;
t491 = rSges(3,1) * t475 - rSges(3,2) * t472;
t405 = pkin(3) * t441 + qJ(4) * t440;
t489 = qJD(1) * t405 + qJD(4) * t438 + t500;
t488 = Icges(3,1) * t475 - t516;
t487 = -Icges(3,2) * t472 + t515;
t486 = Icges(3,5) * t475 - Icges(3,6) * t472;
t425 = -Icges(3,6) * t476 + t473 * t487;
t427 = -Icges(3,5) * t476 + t473 * t488;
t485 = t425 * t472 - t427 * t475;
t426 = Icges(3,6) * t473 + t476 * t487;
t428 = Icges(3,5) * t473 + t476 * t488;
t484 = -t426 * t472 + t428 * t475;
t451 = Icges(3,2) * t475 + t516;
t452 = Icges(3,1) * t472 + t515;
t483 = -t451 * t472 + t452 * t475;
t482 = qJD(4) * t469 * t472 + t404 * t465 + t405 * t504 + t492;
t409 = pkin(4) * t439 - pkin(8) * t511;
t410 = pkin(4) * t441 - pkin(8) * t510;
t481 = t409 * t465 + t410 * t504 + t482;
t480 = qJD(1) * t410 + t473 * t493 + t489;
t479 = (-t409 + t501) * qJD(1) + t476 * t493 + t507;
t468 = qJ(5) + qJ(6);
t467 = cos(t468);
t466 = sin(t468);
t456 = rSges(2,1) * t476 - rSges(2,2) * t473;
t455 = rSges(2,1) * t473 + rSges(2,2) * t476;
t454 = rSges(3,1) * t472 + rSges(3,2) * t475;
t450 = Icges(3,5) * t472 + Icges(3,6) * t475;
t448 = qJD(6) * t475 + t461;
t446 = -t473 * t502 - t504;
t445 = -t476 * t502 + t465;
t434 = (t470 * t474 + t512) * t472;
t433 = (t469 * t474 - t470 * t471) * t472;
t432 = rSges(3,3) * t473 + t476 * t491;
t431 = -rSges(3,3) * t476 + t473 * t491;
t424 = Icges(3,3) * t473 + t476 * t486;
t423 = -Icges(3,3) * t476 + t473 * t486;
t420 = t473 * t495 - t504;
t419 = t476 * t495 + t465;
t412 = (t466 * t469 + t467 * t470) * t472;
t411 = (-t466 * t470 + t467 * t469) * t472;
t402 = t441 * t474 + t513;
t401 = t440 * t474 - t441 * t471;
t400 = t439 * t474 + t514;
t399 = t438 * t474 - t439 * t471;
t398 = pkin(9) * t475 + (pkin(5) * t512 + t470 * t518) * t472;
t395 = t440 * t466 + t441 * t467;
t394 = t440 * t467 - t441 * t466;
t393 = t438 * t466 + t439 * t467;
t392 = t438 * t467 - t439 * t466;
t391 = rSges(4,1) * t441 - rSges(4,2) * t440 + rSges(4,3) * t510;
t390 = rSges(5,1) * t441 + rSges(5,2) * t510 + rSges(5,3) * t440;
t389 = rSges(4,1) * t439 - rSges(4,2) * t438 + rSges(4,3) * t511;
t388 = rSges(5,1) * t439 + rSges(5,2) * t511 + rSges(5,3) * t438;
t375 = rSges(6,1) * t434 + rSges(6,2) * t433 + rSges(6,3) * t475;
t374 = Icges(6,1) * t434 + Icges(6,4) * t433 + Icges(6,5) * t475;
t373 = Icges(6,4) * t434 + Icges(6,2) * t433 + Icges(6,6) * t475;
t372 = Icges(6,5) * t434 + Icges(6,6) * t433 + Icges(6,3) * t475;
t371 = qJD(1) * t432 - t454 * t465 + t449;
t370 = -t454 * t504 + (-t431 - t457) * qJD(1);
t369 = (t431 * t473 + t432 * t476) * qJD(2);
t368 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t475;
t367 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t475;
t366 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t475;
t365 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t475;
t364 = pkin(5) * t513 + t441 * t518 - t476 * t499;
t363 = pkin(5) * t514 + t439 * t518 - t473 * t499;
t362 = rSges(6,1) * t402 + rSges(6,2) * t401 - rSges(6,3) * t510;
t361 = rSges(6,1) * t400 + rSges(6,2) * t399 - rSges(6,3) * t511;
t360 = Icges(6,1) * t402 + Icges(6,4) * t401 - Icges(6,5) * t510;
t359 = Icges(6,1) * t400 + Icges(6,4) * t399 - Icges(6,5) * t511;
t358 = Icges(6,4) * t402 + Icges(6,2) * t401 - Icges(6,6) * t510;
t357 = Icges(6,4) * t400 + Icges(6,2) * t399 - Icges(6,6) * t511;
t356 = Icges(6,5) * t402 + Icges(6,6) * t401 - Icges(6,3) * t510;
t355 = Icges(6,5) * t400 + Icges(6,6) * t399 - Icges(6,3) * t511;
t354 = rSges(7,1) * t395 + rSges(7,2) * t394 - rSges(7,3) * t510;
t353 = rSges(7,1) * t393 + rSges(7,2) * t392 - rSges(7,3) * t511;
t352 = Icges(7,1) * t395 + Icges(7,4) * t394 - Icges(7,5) * t510;
t351 = Icges(7,1) * t393 + Icges(7,4) * t392 - Icges(7,5) * t511;
t350 = Icges(7,4) * t395 + Icges(7,2) * t394 - Icges(7,6) * t510;
t349 = Icges(7,4) * t393 + Icges(7,2) * t392 - Icges(7,6) * t511;
t348 = Icges(7,5) * t395 + Icges(7,6) * t394 - Icges(7,3) * t510;
t347 = Icges(7,5) * t393 + Icges(7,6) * t392 - Icges(7,3) * t511;
t346 = qJD(1) * t391 + t473 * t496 + t500;
t345 = t460 + t476 * t496 + (-t389 + t505) * qJD(1);
t344 = (t389 * t473 + t391 * t476) * qJD(2) + t492;
t343 = qJD(1) * t390 + t473 * t494 + t489;
t342 = t476 * t494 + (-t388 + t501) * qJD(1) + t507;
t341 = (t388 * t473 + t390 * t476) * qJD(2) + t482;
t340 = t362 * t461 - t375 * t445 + t480;
t339 = -t361 * t461 + t375 * t446 + t479;
t338 = t361 * t445 - t362 * t446 + t481;
t337 = t354 * t448 + t364 * t461 - t368 * t419 - t398 * t445 + t480;
t336 = -t353 * t448 - t363 * t461 + t368 * t420 + t398 * t446 + t479;
t335 = t353 * t419 - t354 * t420 + t363 * t445 - t364 * t446 + t481;
t1 = m(7) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(6) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(4) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(5) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(3) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + t419 * ((-t348 * t510 + t394 * t350 + t395 * t352) * t419 + (-t347 * t510 + t349 * t394 + t351 * t395) * t420 + (-t365 * t510 + t366 * t394 + t367 * t395) * t448) / 0.2e1 + t448 * ((t348 * t475 + t350 * t411 + t352 * t412) * t419 + (t347 * t475 + t349 * t411 + t351 * t412) * t420 + (t475 * t365 + t411 * t366 + t412 * t367) * t448) / 0.2e1 + t420 * ((-t348 * t511 + t350 * t392 + t352 * t393) * t419 + (-t347 * t511 + t392 * t349 + t393 * t351) * t420 + (-t365 * t511 + t366 * t392 + t367 * t393) * t448) / 0.2e1 + t445 * ((-t356 * t510 + t401 * t358 + t402 * t360) * t445 + (-t355 * t510 + t357 * t401 + t359 * t402) * t446 + (-t372 * t510 + t373 * t401 + t374 * t402) * t461) / 0.2e1 + t446 * ((-t356 * t511 + t358 * t399 + t360 * t400) * t445 + (-t355 * t511 + t399 * t357 + t400 * t359) * t446 + (-t372 * t511 + t373 * t399 + t374 * t400) * t461) / 0.2e1 + t461 * ((t356 * t475 + t358 * t433 + t360 * t434) * t445 + (t355 * t475 + t357 * t433 + t359 * t434) * t446 + (t475 * t372 + t433 * t373 + t434 * t374) * t461) / 0.2e1 + (Icges(2,3) + m(2) * (t455 ^ 2 + t456 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((((-t425 - t527) * t476 + (t426 - t526) * t473) * t475 + ((t469 * t529 + t470 * t525 - t427) * t476 + (t469 * t528 + t470 * t524 + t428) * t473) * t472) * qJD(2) + ((t451 - t522) * t475 + (t469 * t523 + t470 * t521 + t452) * t472) * qJD(1)) * qJD(1) / 0.2e1 + (((t440 * t529 + t525 * t441 + t485 * t476 + t527 * t510) * t476 + ((-t423 + t484) * t476 + t424 * t473 + t526 * t510 + t524 * t441 + t528 * t440) * t473) * qJD(2) + (t440 * t523 + t441 * t521 + t473 * t450 + t476 * t483 + t510 * t522) * qJD(1)) * t465 / 0.2e1 - (((t423 * t476 + t438 * t529 + t525 * t439 + t527 * t511) * t476 + (t484 * t473 + (-t424 + t485) * t476 + t526 * t511 + t524 * t439 + t528 * t438) * t473) * qJD(2) + (t438 * t523 + t439 * t521 - t476 * t450 + t473 * t483 + t511 * t522) * qJD(1)) * t504 / 0.2e1;
T  = t1;
