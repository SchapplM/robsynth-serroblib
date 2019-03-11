% Calculate kinetic energy for
% S6RRPPRR7
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:08
% EndTime: 2019-03-09 09:17:10
% DurationCPUTime: 2.76s
% Computational Cost: add. (1602->288), mult. (4014->440), div. (0->0), fcn. (4635->10), ass. (0->135)
t532 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t531 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t530 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t529 = Icges(4,6) - Icges(3,6) - Icges(5,5);
t528 = Icges(5,6) + Icges(4,4) + Icges(3,5);
t527 = -Icges(5,3) - Icges(3,3) - Icges(4,2);
t474 = sin(pkin(6));
t475 = cos(pkin(6));
t481 = cos(qJ(2));
t482 = cos(qJ(1));
t509 = t481 * t482;
t478 = sin(qJ(2));
t479 = sin(qJ(1));
t513 = t478 * t479;
t453 = -t475 * t509 + t513;
t510 = t479 * t481;
t512 = t478 * t482;
t454 = t475 * t512 + t510;
t455 = t475 * t510 + t512;
t456 = -t475 * t513 + t509;
t508 = t482 * t474;
t511 = t479 * t474;
t518 = (-t529 * t453 - t528 * t454 - t527 * t508) * t482 + (t529 * t455 + t528 * t456 - t527 * t511) * t479;
t526 = t518 * t474;
t525 = t531 * t453 - t530 * t454 - t529 * t508;
t524 = t531 * t455 - t530 * t456 + t529 * t511;
t523 = t530 * t453 - t532 * t454 + t528 * t508;
t522 = -t530 * t455 + t532 * t456 + t528 * t511;
t521 = t527 * t475 + (-t528 * t478 + t529 * t481) * t474;
t520 = t529 * t475 + (-t530 * t478 - t531 * t481) * t474;
t519 = t528 * t475 + (t532 * t478 + t530 * t481) * t474;
t515 = cos(qJ(5));
t514 = t474 * t478;
t417 = pkin(2) * t454 + qJ(3) * t453;
t418 = pkin(2) * t456 + qJ(3) * t455;
t503 = qJD(2) * t474;
t470 = t479 * t503;
t498 = t482 * t503;
t507 = t417 * t470 + t418 * t498;
t433 = pkin(3) * t454 + qJ(4) * t508;
t506 = -t417 - t433;
t457 = (pkin(2) * t478 - qJ(3) * t481) * t474;
t505 = -pkin(3) * t514 + qJ(4) * t475 - t457;
t431 = qJD(5) * t456 + t470;
t504 = qJD(1) * (pkin(1) * t479 - pkin(8) * t508);
t502 = qJD(3) * t481;
t501 = qJD(4) * t479;
t471 = qJD(2) * t475 + qJD(1);
t460 = qJD(1) * (pkin(1) * t482 + pkin(8) * t511);
t500 = qJD(3) * t453 + t471 * t418 + t460;
t499 = t474 * t515;
t459 = qJD(5) * t514 + t471;
t497 = qJD(3) * t455 - t504;
t494 = qJD(2) * (rSges(5,3) * t475 - (-rSges(5,1) * t481 - rSges(5,2) * t478) * t474 + t505);
t493 = qJD(2) * (-(-pkin(4) * t481 + pkin(9) * t478) * t474 + t505);
t492 = (-rSges(4,2) * t475 - (rSges(4,1) * t478 - rSges(4,3) * t481) * t474 - t457) * t503;
t432 = qJD(5) * t454 - t498;
t434 = pkin(3) * t456 - qJ(4) * t511;
t491 = qJD(4) * t508 + t471 * t434 + t500;
t487 = -qJD(4) * t475 + t433 * t470 + t434 * t498 + t507;
t419 = pkin(4) * t453 + pkin(9) * t454;
t420 = pkin(4) * t455 + pkin(9) * t456;
t486 = t419 * t470 + t420 * t498 - t474 * t502 + t487;
t485 = t471 * t420 + t493 * t511 + t491;
t484 = (t482 * t493 - t501) * t474 + (-t419 + t506) * t471 + t497;
t480 = cos(qJ(6));
t477 = sin(qJ(5));
t476 = sin(qJ(6));
t464 = rSges(2,1) * t482 - rSges(2,2) * t479;
t463 = rSges(2,1) * t479 + rSges(2,2) * t482;
t452 = -t475 * t477 - t481 * t499;
t451 = t474 * t477 * t481 - t475 * t515;
t445 = rSges(3,3) * t475 + (rSges(3,1) * t478 + rSges(3,2) * t481) * t474;
t430 = t455 * t515 - t477 * t511;
t429 = t455 * t477 + t479 * t499;
t428 = t453 * t515 + t477 * t508;
t427 = t453 * t477 - t482 * t499;
t426 = t452 * t480 + t476 * t514;
t425 = -t452 * t476 + t480 * t514;
t424 = -qJD(6) * t451 + t459;
t416 = pkin(5) * t452 - pkin(10) * t451;
t413 = rSges(3,1) * t456 - rSges(3,2) * t455 + rSges(3,3) * t511;
t412 = rSges(4,1) * t456 + rSges(4,2) * t511 + rSges(4,3) * t455;
t411 = rSges(5,1) * t455 - rSges(5,2) * t456 - rSges(5,3) * t511;
t410 = rSges(3,1) * t454 - rSges(3,2) * t453 - rSges(3,3) * t508;
t409 = rSges(4,1) * t454 - rSges(4,2) * t508 + rSges(4,3) * t453;
t408 = rSges(5,1) * t453 - rSges(5,2) * t454 + rSges(5,3) * t508;
t385 = rSges(6,1) * t452 + rSges(6,2) * t451 + rSges(6,3) * t514;
t384 = Icges(6,1) * t452 + Icges(6,4) * t451 + Icges(6,5) * t514;
t383 = Icges(6,4) * t452 + Icges(6,2) * t451 + Icges(6,6) * t514;
t382 = Icges(6,5) * t452 + Icges(6,6) * t451 + Icges(6,3) * t514;
t381 = t430 * t480 + t456 * t476;
t380 = -t430 * t476 + t456 * t480;
t379 = t428 * t480 + t454 * t476;
t378 = -t428 * t476 + t454 * t480;
t377 = qJD(6) * t427 + t432;
t376 = qJD(6) * t429 + t431;
t375 = pkin(5) * t430 + pkin(10) * t429;
t374 = pkin(5) * t428 + pkin(10) * t427;
t373 = rSges(6,1) * t430 - rSges(6,2) * t429 + rSges(6,3) * t456;
t372 = rSges(6,1) * t428 - rSges(6,2) * t427 + rSges(6,3) * t454;
t371 = Icges(6,1) * t430 - Icges(6,4) * t429 + Icges(6,5) * t456;
t370 = Icges(6,1) * t428 - Icges(6,4) * t427 + Icges(6,5) * t454;
t369 = Icges(6,4) * t430 - Icges(6,2) * t429 + Icges(6,6) * t456;
t368 = Icges(6,4) * t428 - Icges(6,2) * t427 + Icges(6,6) * t454;
t367 = Icges(6,5) * t430 - Icges(6,6) * t429 + Icges(6,3) * t456;
t366 = Icges(6,5) * t428 - Icges(6,6) * t427 + Icges(6,3) * t454;
t365 = rSges(7,1) * t426 + rSges(7,2) * t425 - rSges(7,3) * t451;
t364 = Icges(7,1) * t426 + Icges(7,4) * t425 - Icges(7,5) * t451;
t363 = Icges(7,4) * t426 + Icges(7,2) * t425 - Icges(7,6) * t451;
t362 = Icges(7,5) * t426 + Icges(7,6) * t425 - Icges(7,3) * t451;
t361 = t413 * t471 - t445 * t470 + t460;
t360 = -t410 * t471 - t445 * t498 - t504;
t359 = (t410 * t479 + t413 * t482) * t503;
t358 = rSges(7,1) * t381 + rSges(7,2) * t380 + rSges(7,3) * t429;
t357 = rSges(7,1) * t379 + rSges(7,2) * t378 + rSges(7,3) * t427;
t356 = Icges(7,1) * t381 + Icges(7,4) * t380 + Icges(7,5) * t429;
t355 = Icges(7,1) * t379 + Icges(7,4) * t378 + Icges(7,5) * t427;
t354 = Icges(7,4) * t381 + Icges(7,2) * t380 + Icges(7,6) * t429;
t353 = Icges(7,4) * t379 + Icges(7,2) * t378 + Icges(7,6) * t427;
t352 = Icges(7,5) * t381 + Icges(7,6) * t380 + Icges(7,3) * t429;
t351 = Icges(7,5) * t379 + Icges(7,6) * t378 + Icges(7,3) * t427;
t350 = t412 * t471 + t479 * t492 + t500;
t349 = (-t409 - t417) * t471 + t482 * t492 + t497;
t348 = (-t502 + (t409 * t479 + t412 * t482) * qJD(2)) * t474 + t507;
t347 = t411 * t471 + t494 * t511 + t491;
t346 = (-t408 + t506) * t471 + (t482 * t494 - t501) * t474 + t497;
t345 = (-t502 + (t408 * t479 + t411 * t482) * qJD(2)) * t474 + t487;
t344 = t373 * t459 - t385 * t431 + t485;
t343 = -t372 * t459 + t385 * t432 + t484;
t342 = t372 * t431 - t373 * t432 + t486;
t341 = t358 * t424 - t365 * t376 + t375 * t459 - t416 * t431 + t485;
t340 = -t357 * t424 + t365 * t377 - t374 * t459 + t416 * t432 + t484;
t339 = t357 * t376 - t358 * t377 + t374 * t431 - t375 * t432 + t486;
t1 = t431 * ((t456 * t367 - t429 * t369 + t430 * t371) * t431 + (t366 * t456 - t368 * t429 + t370 * t430) * t432 + (t382 * t456 - t383 * t429 + t384 * t430) * t459) / 0.2e1 + t432 * ((t367 * t454 - t369 * t427 + t371 * t428) * t431 + (t454 * t366 - t427 * t368 + t428 * t370) * t432 + (t382 * t454 - t383 * t427 + t384 * t428) * t459) / 0.2e1 + t459 * ((t367 * t514 + t369 * t451 + t371 * t452) * t431 + (t366 * t514 + t368 * t451 + t370 * t452) * t432 + (t382 * t514 + t383 * t451 + t384 * t452) * t459) / 0.2e1 + t376 * ((t429 * t352 + t380 * t354 + t381 * t356) * t376 + (t351 * t429 + t353 * t380 + t355 * t381) * t377 + (t362 * t429 + t363 * t380 + t364 * t381) * t424) / 0.2e1 + t377 * ((t352 * t427 + t354 * t378 + t356 * t379) * t376 + (t427 * t351 + t378 * t353 + t379 * t355) * t377 + (t362 * t427 + t363 * t378 + t364 * t379) * t424) / 0.2e1 + t424 * ((-t352 * t451 + t354 * t425 + t356 * t426) * t376 + (-t351 * t451 + t353 * t425 + t355 * t426) * t377 + (-t451 * t362 + t425 * t363 + t426 * t364) * t424) / 0.2e1 + m(3) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(4) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(5) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(6) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(7) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + (m(2) * (t463 ^ 2 + t464 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t518 * t475 + ((t523 * t478 + t525 * t481) * t482 + (t522 * t478 - t524 * t481) * t479) * t474) * t503 + (-t521 * t475 + (t519 * t478 - t520 * t481) * t474) * t471) * t471 / 0.2e1 + (((-t525 * t455 + t523 * t456) * t482 + (t524 * t455 + t522 * t456 + t526) * t479) * t503 + (t520 * t455 + t519 * t456 - t521 * t511) * t471) * t470 / 0.2e1 - (((-t525 * t453 + t523 * t454 - t526) * t482 + (t524 * t453 + t522 * t454) * t479) * t503 + (t520 * t453 + t519 * t454 + t521 * t508) * t471) * t498 / 0.2e1;
T  = t1;
