% Calculate kinetic energy for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:49
% EndTime: 2019-03-10 03:27:51
% DurationCPUTime: 1.96s
% Computational Cost: add. (2017->287), mult. (1761->466), div. (0->0), fcn. (1606->12), ass. (0->171)
t453 = qJ(2) + qJ(3);
t447 = sin(t453);
t527 = pkin(3) * t447;
t449 = qJ(4) + t453;
t439 = sin(t449);
t526 = pkin(4) * t439;
t458 = cos(qJ(2));
t524 = t458 * pkin(2);
t455 = sin(qJ(2));
t522 = Icges(3,4) * t455;
t521 = Icges(3,4) * t458;
t520 = Icges(4,4) * t447;
t448 = cos(t453);
t519 = Icges(4,4) * t448;
t518 = Icges(5,4) * t439;
t440 = cos(t449);
t517 = Icges(5,4) * t440;
t442 = qJ(5) + t449;
t436 = sin(t442);
t516 = Icges(6,4) * t436;
t437 = cos(t442);
t515 = Icges(6,4) * t437;
t456 = sin(qJ(1));
t514 = t436 * t456;
t459 = cos(qJ(1));
t513 = t436 * t459;
t454 = sin(qJ(6));
t512 = t454 * t456;
t511 = t454 * t459;
t457 = cos(qJ(6));
t510 = t456 * t457;
t509 = t457 * t459;
t380 = -pkin(8) * t459 + t456 * t524;
t381 = pkin(8) * t456 + t459 * t524;
t446 = qJD(2) * t456;
t502 = qJD(2) * t459;
t508 = t380 * t446 + t381 * t502;
t435 = pkin(1) * t456 - pkin(7) * t459;
t507 = -t380 - t435;
t506 = pkin(4) * t440;
t505 = pkin(3) * t448;
t427 = qJD(3) * t456 + t446;
t501 = qJD(6) * t436;
t500 = -qJD(2) - qJD(3);
t499 = pkin(2) * qJD(2) * t455;
t353 = -pkin(9) * t459 + t456 * t505;
t498 = -t353 + t507;
t418 = qJD(4) * t456 + t427;
t497 = -qJD(4) + t500;
t496 = t459 * t499;
t346 = -pkin(10) * t459 + t456 * t506;
t495 = -t346 + t498;
t393 = qJD(5) * t456 + t418;
t494 = pkin(5) * t437 + pkin(11) * t436;
t493 = rSges(3,1) * t458 - rSges(3,2) * t455;
t492 = rSges(4,1) * t448 - rSges(4,2) * t447;
t491 = rSges(5,1) * t440 - rSges(5,2) * t439;
t490 = rSges(6,1) * t437 - rSges(6,2) * t436;
t489 = Icges(3,1) * t458 - t522;
t488 = Icges(4,1) * t448 - t520;
t487 = Icges(5,1) * t440 - t518;
t486 = Icges(6,1) * t437 - t516;
t485 = -Icges(3,2) * t455 + t521;
t484 = -Icges(4,2) * t447 + t519;
t483 = -Icges(5,2) * t439 + t517;
t482 = -Icges(6,2) * t436 + t515;
t481 = Icges(3,5) * t458 - Icges(3,6) * t455;
t480 = Icges(4,5) * t448 - Icges(4,6) * t447;
t479 = Icges(5,5) * t440 - Icges(5,6) * t439;
t478 = Icges(6,5) * t437 - Icges(6,6) * t436;
t397 = -Icges(3,6) * t459 + t456 * t485;
t399 = -Icges(3,5) * t459 + t456 * t489;
t477 = t397 * t455 - t399 * t458;
t398 = Icges(3,6) * t456 + t459 * t485;
t400 = Icges(3,5) * t456 + t459 * t489;
t476 = -t398 * t455 + t400 * t458;
t430 = Icges(3,2) * t458 + t522;
t431 = Icges(3,1) * t455 + t521;
t475 = -t430 * t455 + t431 * t458;
t428 = t500 * t459;
t474 = t428 * t527 - t496;
t354 = pkin(9) * t456 + t459 * t505;
t473 = t427 * t353 - t354 * t428 + t508;
t394 = (-qJD(5) + t497) * t459;
t426 = qJD(1) * (pkin(1) * t459 + pkin(7) * t456);
t472 = qJD(1) * t381 - t456 * t499 + t426;
t419 = t497 * t459;
t471 = t419 * t526 + t474;
t470 = (Icges(6,5) * t436 + Icges(6,6) * t437) * qJD(1) + (-Icges(6,3) * t459 + t456 * t478) * t394 + (Icges(6,3) * t456 + t459 * t478) * t393;
t469 = (Icges(5,5) * t439 + Icges(5,6) * t440) * qJD(1) + (-Icges(5,3) * t459 + t456 * t479) * t419 + (Icges(5,3) * t456 + t459 * t479) * t418;
t468 = (Icges(4,5) * t447 + Icges(4,6) * t448) * qJD(1) + (-Icges(4,3) * t459 + t456 * t480) * t428 + (Icges(4,3) * t456 + t459 * t480) * t427;
t347 = pkin(10) * t456 + t459 * t506;
t467 = t418 * t346 - t347 * t419 + t473;
t466 = qJD(1) * t354 - t427 * t527 + t472;
t465 = qJD(1) * t347 - t418 * t526 + t466;
t365 = -Icges(6,6) * t459 + t456 * t482;
t366 = Icges(6,6) * t456 + t459 * t482;
t367 = -Icges(6,5) * t459 + t456 * t486;
t368 = Icges(6,5) * t456 + t459 * t486;
t409 = Icges(6,2) * t437 + t516;
t410 = Icges(6,1) * t436 + t515;
t464 = (-t366 * t436 + t368 * t437) * t393 + (-t365 * t436 + t367 * t437) * t394 + (-t409 * t436 + t410 * t437) * qJD(1);
t374 = -Icges(5,6) * t459 + t456 * t483;
t375 = Icges(5,6) * t456 + t459 * t483;
t376 = -Icges(5,5) * t459 + t456 * t487;
t377 = Icges(5,5) * t456 + t459 * t487;
t415 = Icges(5,2) * t440 + t518;
t416 = Icges(5,1) * t439 + t517;
t463 = (-t375 * t439 + t377 * t440) * t418 + (-t374 * t439 + t376 * t440) * t419 + (-t415 * t439 + t416 * t440) * qJD(1);
t384 = -Icges(4,6) * t459 + t456 * t484;
t385 = Icges(4,6) * t456 + t459 * t484;
t386 = -Icges(4,5) * t459 + t456 * t488;
t387 = Icges(4,5) * t456 + t459 * t488;
t421 = Icges(4,2) * t448 + t520;
t422 = Icges(4,1) * t447 + t519;
t462 = (-t385 * t447 + t387 * t448) * t427 + (-t384 * t447 + t386 * t448) * t428 + (-t421 * t447 + t422 * t448) * qJD(1);
t434 = rSges(2,1) * t459 - rSges(2,2) * t456;
t433 = rSges(2,1) * t456 + rSges(2,2) * t459;
t432 = rSges(3,1) * t455 + rSges(3,2) * t458;
t429 = Icges(3,5) * t455 + Icges(3,6) * t458;
t425 = -qJD(6) * t437 + qJD(1);
t423 = rSges(4,1) * t447 + rSges(4,2) * t448;
t417 = rSges(5,1) * t439 + rSges(5,2) * t440;
t412 = pkin(5) * t436 - pkin(11) * t437;
t411 = rSges(6,1) * t436 + rSges(6,2) * t437;
t407 = rSges(3,3) * t456 + t459 * t493;
t406 = -rSges(3,3) * t459 + t456 * t493;
t404 = t437 * t509 + t512;
t403 = -t437 * t511 + t510;
t402 = t437 * t510 - t511;
t401 = -t437 * t512 - t509;
t396 = Icges(3,3) * t456 + t459 * t481;
t395 = -Icges(3,3) * t459 + t456 * t481;
t391 = rSges(4,3) * t456 + t459 * t492;
t390 = -rSges(4,3) * t459 + t456 * t492;
t389 = t494 * t459;
t388 = t494 * t456;
t379 = rSges(5,3) * t456 + t459 * t491;
t378 = -rSges(5,3) * t459 + t456 * t491;
t370 = rSges(6,3) * t456 + t459 * t490;
t369 = -rSges(6,3) * t459 + t456 * t490;
t360 = t456 * t501 + t394;
t359 = t459 * t501 + t393;
t358 = -rSges(7,3) * t437 + (rSges(7,1) * t457 - rSges(7,2) * t454) * t436;
t357 = -Icges(7,5) * t437 + (Icges(7,1) * t457 - Icges(7,4) * t454) * t436;
t356 = -Icges(7,6) * t437 + (Icges(7,4) * t457 - Icges(7,2) * t454) * t436;
t355 = -Icges(7,3) * t437 + (Icges(7,5) * t457 - Icges(7,6) * t454) * t436;
t351 = qJD(1) * t407 - t432 * t446 + t426;
t350 = -t432 * t502 + (-t406 - t435) * qJD(1);
t348 = (t406 * t456 + t407 * t459) * qJD(2);
t344 = rSges(7,1) * t404 + rSges(7,2) * t403 + rSges(7,3) * t513;
t343 = rSges(7,1) * t402 + rSges(7,2) * t401 + rSges(7,3) * t514;
t342 = Icges(7,1) * t404 + Icges(7,4) * t403 + Icges(7,5) * t513;
t341 = Icges(7,1) * t402 + Icges(7,4) * t401 + Icges(7,5) * t514;
t340 = Icges(7,4) * t404 + Icges(7,2) * t403 + Icges(7,6) * t513;
t339 = Icges(7,4) * t402 + Icges(7,2) * t401 + Icges(7,6) * t514;
t338 = Icges(7,5) * t404 + Icges(7,6) * t403 + Icges(7,3) * t513;
t337 = Icges(7,5) * t402 + Icges(7,6) * t401 + Icges(7,3) * t514;
t335 = qJD(1) * t391 - t423 * t427 + t472;
t334 = -t496 + t423 * t428 + (-t390 + t507) * qJD(1);
t333 = t390 * t427 - t391 * t428 + t508;
t332 = qJD(1) * t379 - t417 * t418 + t466;
t331 = t417 * t419 + (-t378 + t498) * qJD(1) + t474;
t330 = t378 * t418 - t379 * t419 + t473;
t329 = qJD(1) * t370 - t393 * t411 + t465;
t328 = t394 * t411 + (-t369 + t495) * qJD(1) + t471;
t327 = qJD(1) * t389 + t344 * t425 - t358 * t359 - t393 * t412 + t465;
t326 = -t343 * t425 + t358 * t360 + t394 * t412 + (-t388 + t495) * qJD(1) + t471;
t325 = t369 * t393 - t370 * t394 + t467;
t324 = t343 * t359 - t344 * t360 + t388 * t393 - t389 * t394 + t467;
t1 = m(7) * (t324 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(6) * (t325 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(5) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(4) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(3) * (t348 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + t360 * ((t338 * t514 + t340 * t401 + t342 * t402) * t359 + (t337 * t514 + t401 * t339 + t402 * t341) * t360 + (t355 * t514 + t356 * t401 + t357 * t402) * t425) / 0.2e1 + t359 * ((t338 * t513 + t403 * t340 + t404 * t342) * t359 + (t337 * t513 + t339 * t403 + t341 * t404) * t360 + (t355 * t513 + t356 * t403 + t357 * t404) * t425) / 0.2e1 + t425 * ((-t337 * t360 - t338 * t359 - t355 * t425) * t437 + ((-t340 * t454 + t342 * t457) * t359 + (-t339 * t454 + t341 * t457) * t360 + (-t356 * t454 + t357 * t457) * t425) * t436) / 0.2e1 + t393 * (t470 * t456 + t464 * t459) / 0.2e1 + t394 * (t464 * t456 - t470 * t459) / 0.2e1 + t418 * (t469 * t456 + t463 * t459) / 0.2e1 + t419 * (t463 * t456 - t469 * t459) / 0.2e1 + t427 * (t468 * t456 + t462 * t459) / 0.2e1 + t428 * (t462 * t456 - t468 * t459) / 0.2e1 - ((-t459 * t429 + t456 * t475) * qJD(1) + (t459 ^ 2 * t395 + (t476 * t456 + (-t396 + t477) * t459) * t456) * qJD(2)) * t502 / 0.2e1 + ((t456 * t429 + t459 * t475) * qJD(1) + (t456 ^ 2 * t396 + (t477 * t459 + (-t395 + t476) * t456) * t459) * qJD(2)) * t446 / 0.2e1 + (m(2) * (t433 ^ 2 + t434 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t366 * t437 + t368 * t436) * t393 + (t365 * t437 + t367 * t436) * t394 + (t375 * t440 + t377 * t439) * t418 + (t374 * t440 + t376 * t439) * t419 + ((t398 * t458 + t400 * t455) * t456 - (t397 * t458 + t399 * t455) * t459) * qJD(2) + (t385 * t448 + t387 * t447) * t427 + (t384 * t448 + t386 * t447) * t428 + (t437 * t409 + t436 * t410 + t440 * t415 + t439 * t416 + t448 * t421 + t447 * t422 + t458 * t430 + t455 * t431) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
