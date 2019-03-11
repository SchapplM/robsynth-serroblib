% Calculate kinetic energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:11
% EndTime: 2019-03-09 09:07:14
% DurationCPUTime: 2.83s
% Computational Cost: add. (1602->288), mult. (4014->440), div. (0->0), fcn. (4635->10), ass. (0->136)
t540 = Icges(3,1) + Icges(4,1) + Icges(5,1);
t539 = Icges(3,4) - Icges(5,4) - Icges(4,5);
t538 = -Icges(5,5) + Icges(4,4) + Icges(3,5);
t537 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t536 = Icges(4,6) - Icges(5,6) - Icges(3,6);
t535 = Icges(3,3) + Icges(4,2) + Icges(5,3);
t481 = sin(pkin(6));
t482 = cos(pkin(6));
t488 = cos(qJ(2));
t489 = cos(qJ(1));
t515 = t488 * t489;
t485 = sin(qJ(2));
t486 = sin(qJ(1));
t518 = t485 * t486;
t462 = -t482 * t515 + t518;
t516 = t486 * t488;
t517 = t485 * t489;
t463 = t482 * t517 + t516;
t464 = t482 * t516 + t517;
t465 = -t482 * t518 + t515;
t519 = t481 * t489;
t521 = t481 * t486;
t526 = (-t536 * t462 - t538 * t463 + t535 * t519) * t489 + (t536 * t464 + t538 * t465 + t535 * t521) * t486;
t534 = t526 * t481;
t533 = t537 * t462 - t539 * t463 - t536 * t519;
t532 = t537 * t464 - t539 * t465 + t536 * t521;
t531 = t539 * t462 - t540 * t463 + t538 * t519;
t530 = -t539 * t464 + t540 * t465 + t538 * t521;
t529 = t536 * t482 + (-t539 * t485 - t537 * t488) * t481;
t528 = t535 * t482 + (t538 * t485 - t536 * t488) * t481;
t527 = t538 * t482 + (t540 * t485 + t539 * t488) * t481;
t523 = cos(qJ(5));
t522 = t481 * t485;
t520 = t481 * t488;
t426 = pkin(2) * t463 + qJ(3) * t462;
t427 = pkin(2) * t465 + qJ(3) * t464;
t510 = qJD(2) * t481;
t478 = t486 * t510;
t505 = t489 * t510;
t514 = t426 * t478 + t427 * t505;
t442 = pkin(3) * t463 + qJ(4) * t519;
t513 = -t426 - t442;
t466 = (pkin(2) * t485 - qJ(3) * t488) * t481;
t512 = -pkin(3) * t522 + qJ(4) * t482 - t466;
t440 = -qJD(5) * t464 + t478;
t511 = qJD(1) * (pkin(1) * t486 - pkin(8) * t519);
t509 = qJD(3) * t488;
t508 = qJD(4) * t486;
t479 = qJD(2) * t482 + qJD(1);
t469 = qJD(1) * (pkin(1) * t489 + pkin(8) * t521);
t507 = qJD(3) * t462 + t479 * t427 + t469;
t506 = t481 * t523;
t468 = qJD(5) * t520 + t479;
t504 = qJD(3) * t464 - t511;
t501 = qJD(2) * (rSges(5,3) * t482 - (rSges(5,1) * t485 - rSges(5,2) * t488) * t481 + t512);
t500 = qJD(2) * (-(pkin(4) * t485 + pkin(9) * t488) * t481 + t512);
t499 = (-rSges(4,2) * t482 - (rSges(4,1) * t485 - rSges(4,3) * t488) * t481 - t466) * t510;
t441 = -qJD(5) * t462 - t505;
t443 = pkin(3) * t465 - qJ(4) * t521;
t498 = qJD(4) * t519 + t479 * t443 + t507;
t494 = -qJD(4) * t482 + t442 * t478 + t443 * t505 + t514;
t428 = pkin(4) * t463 - pkin(9) * t462;
t429 = pkin(4) * t465 - pkin(9) * t464;
t493 = t428 * t478 + t429 * t505 - t481 * t509 + t494;
t492 = t479 * t429 + t500 * t521 + t498;
t491 = (t489 * t500 - t508) * t481 + (-t428 + t513) * t479 + t504;
t487 = cos(qJ(6));
t484 = sin(qJ(5));
t483 = sin(qJ(6));
t473 = rSges(2,1) * t489 - rSges(2,2) * t486;
t472 = rSges(2,1) * t486 + rSges(2,2) * t489;
t461 = -t482 * t484 + t485 * t506;
t460 = t482 * t523 + t484 * t522;
t455 = rSges(3,3) * t482 + (rSges(3,1) * t485 + rSges(3,2) * t488) * t481;
t439 = t465 * t523 - t484 * t521;
t438 = t465 * t484 + t486 * t506;
t437 = t463 * t523 + t484 * t519;
t436 = t463 * t484 - t489 * t506;
t435 = t461 * t487 + t483 * t520;
t434 = -t461 * t483 + t487 * t520;
t433 = qJD(6) * t460 + t468;
t425 = pkin(5) * t461 + pkin(10) * t460;
t422 = rSges(3,1) * t465 - rSges(3,2) * t464 + rSges(3,3) * t521;
t421 = rSges(4,1) * t465 + rSges(4,2) * t521 + rSges(4,3) * t464;
t420 = rSges(5,1) * t465 + rSges(5,2) * t464 - rSges(5,3) * t521;
t419 = rSges(3,1) * t463 - rSges(3,2) * t462 - rSges(3,3) * t519;
t418 = rSges(4,1) * t463 - rSges(4,2) * t519 + rSges(4,3) * t462;
t417 = rSges(5,1) * t463 + rSges(5,2) * t462 + rSges(5,3) * t519;
t394 = rSges(6,1) * t461 - rSges(6,2) * t460 + rSges(6,3) * t520;
t393 = Icges(6,1) * t461 - Icges(6,4) * t460 + Icges(6,5) * t520;
t392 = Icges(6,4) * t461 - Icges(6,2) * t460 + Icges(6,6) * t520;
t391 = Icges(6,5) * t461 - Icges(6,6) * t460 + Icges(6,3) * t520;
t390 = t439 * t487 - t464 * t483;
t389 = -t439 * t483 - t464 * t487;
t388 = t437 * t487 - t462 * t483;
t387 = -t437 * t483 - t462 * t487;
t386 = qJD(6) * t436 + t441;
t385 = qJD(6) * t438 + t440;
t384 = pkin(5) * t439 + pkin(10) * t438;
t383 = pkin(5) * t437 + pkin(10) * t436;
t382 = rSges(6,1) * t439 - rSges(6,2) * t438 - rSges(6,3) * t464;
t381 = rSges(6,1) * t437 - rSges(6,2) * t436 - rSges(6,3) * t462;
t380 = Icges(6,1) * t439 - Icges(6,4) * t438 - Icges(6,5) * t464;
t379 = Icges(6,1) * t437 - Icges(6,4) * t436 - Icges(6,5) * t462;
t378 = Icges(6,4) * t439 - Icges(6,2) * t438 - Icges(6,6) * t464;
t377 = Icges(6,4) * t437 - Icges(6,2) * t436 - Icges(6,6) * t462;
t376 = Icges(6,5) * t439 - Icges(6,6) * t438 - Icges(6,3) * t464;
t375 = Icges(6,5) * t437 - Icges(6,6) * t436 - Icges(6,3) * t462;
t374 = rSges(7,1) * t435 + rSges(7,2) * t434 + rSges(7,3) * t460;
t373 = Icges(7,1) * t435 + Icges(7,4) * t434 + Icges(7,5) * t460;
t372 = Icges(7,4) * t435 + Icges(7,2) * t434 + Icges(7,6) * t460;
t371 = Icges(7,5) * t435 + Icges(7,6) * t434 + Icges(7,3) * t460;
t370 = t422 * t479 - t455 * t478 + t469;
t369 = -t419 * t479 - t455 * t505 - t511;
t368 = (t419 * t486 + t422 * t489) * t510;
t367 = rSges(7,1) * t390 + rSges(7,2) * t389 + rSges(7,3) * t438;
t366 = rSges(7,1) * t388 + rSges(7,2) * t387 + rSges(7,3) * t436;
t365 = Icges(7,1) * t390 + Icges(7,4) * t389 + Icges(7,5) * t438;
t364 = Icges(7,1) * t388 + Icges(7,4) * t387 + Icges(7,5) * t436;
t363 = Icges(7,4) * t390 + Icges(7,2) * t389 + Icges(7,6) * t438;
t362 = Icges(7,4) * t388 + Icges(7,2) * t387 + Icges(7,6) * t436;
t361 = Icges(7,5) * t390 + Icges(7,6) * t389 + Icges(7,3) * t438;
t360 = Icges(7,5) * t388 + Icges(7,6) * t387 + Icges(7,3) * t436;
t359 = t421 * t479 + t486 * t499 + t507;
t358 = (-t418 - t426) * t479 + t489 * t499 + t504;
t357 = (-t509 + (t418 * t486 + t421 * t489) * qJD(2)) * t481 + t514;
t356 = t420 * t479 + t501 * t521 + t498;
t355 = (-t417 + t513) * t479 + (t489 * t501 - t508) * t481 + t504;
t354 = (-t509 + (t417 * t486 + t420 * t489) * qJD(2)) * t481 + t494;
t353 = t382 * t468 - t394 * t440 + t492;
t352 = -t381 * t468 + t394 * t441 + t491;
t351 = t381 * t440 - t382 * t441 + t493;
t350 = t367 * t433 - t374 * t385 + t384 * t468 - t425 * t440 + t492;
t349 = -t366 * t433 + t374 * t386 - t383 * t468 + t425 * t441 + t491;
t348 = t366 * t385 - t367 * t386 + t383 * t440 - t384 * t441 + t493;
t1 = m(5) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(6) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(7) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + t468 * ((t376 * t520 - t378 * t460 + t380 * t461) * t440 + (t375 * t520 - t377 * t460 + t379 * t461) * t441 + (t391 * t520 - t392 * t460 + t393 * t461) * t468) / 0.2e1 + t440 * ((-t376 * t464 - t378 * t438 + t380 * t439) * t440 + (-t375 * t464 - t377 * t438 + t379 * t439) * t441 + (-t391 * t464 - t392 * t438 + t393 * t439) * t468) / 0.2e1 + t441 * ((-t376 * t462 - t378 * t436 + t380 * t437) * t440 + (-t375 * t462 - t377 * t436 + t379 * t437) * t441 + (-t391 * t462 - t392 * t436 + t393 * t437) * t468) / 0.2e1 + t385 * ((t438 * t361 + t389 * t363 + t390 * t365) * t385 + (t360 * t438 + t362 * t389 + t364 * t390) * t386 + (t371 * t438 + t372 * t389 + t373 * t390) * t433) / 0.2e1 + t386 * ((t361 * t436 + t363 * t387 + t365 * t388) * t385 + (t436 * t360 + t387 * t362 + t388 * t364) * t386 + (t371 * t436 + t372 * t387 + t373 * t388) * t433) / 0.2e1 + t433 * ((t361 * t460 + t363 * t434 + t365 * t435) * t385 + (t360 * t460 + t362 * t434 + t364 * t435) * t386 + (t371 * t460 + t372 * t434 + t373 * t435) * t433) / 0.2e1 + m(3) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t472 ^ 2 + t473 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t526 * t482 + ((t531 * t485 + t533 * t488) * t489 + (t530 * t485 - t532 * t488) * t486) * t481) * t510 + (t528 * t482 + (t527 * t485 - t529 * t488) * t481) * t479) * t479 / 0.2e1 + (((-t533 * t464 + t531 * t465) * t489 + (t532 * t464 + t530 * t465 + t534) * t486) * t510 + (t529 * t464 + t527 * t465 + t528 * t521) * t479) * t478 / 0.2e1 - (((-t533 * t462 + t531 * t463 - t534) * t489 + (t532 * t462 + t530 * t463) * t486) * t510 + (t529 * t462 + t527 * t463 - t528 * t519) * t479) * t505 / 0.2e1;
T  = t1;
