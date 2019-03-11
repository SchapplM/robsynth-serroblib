% Calculate kinetic energy for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:52
% EndTime: 2019-03-09 18:10:55
% DurationCPUTime: 3.47s
% Computational Cost: add. (1713->281), mult. (2274->434), div. (0->0), fcn. (2353->10), ass. (0->152)
t539 = Icges(4,4) - Icges(5,5);
t538 = Icges(4,1) + Icges(5,1);
t537 = Icges(4,2) + Icges(5,3);
t459 = qJ(2) + qJ(3);
t458 = cos(t459);
t536 = t539 * t458;
t457 = sin(t459);
t535 = t539 * t457;
t534 = Icges(5,4) + Icges(4,5);
t533 = Icges(4,6) - Icges(5,6);
t532 = t537 * t457 - t536;
t531 = t538 * t458 - t535;
t530 = Icges(5,2) + Icges(4,3);
t463 = sin(qJ(1));
t466 = cos(qJ(1));
t529 = t532 * t463 + t533 * t466;
t528 = -t533 * t463 + t532 * t466;
t527 = t531 * t463 - t534 * t466;
t526 = t534 * t463 + t531 * t466;
t525 = -t537 * t458 - t535;
t524 = t538 * t457 + t536;
t523 = -t533 * t457 + t534 * t458;
t456 = qJD(2) * t463;
t442 = qJD(3) * t463 + t456;
t443 = (-qJD(2) - qJD(3)) * t466;
t522 = (t457 * t529 + t458 * t527) * t443 + (t457 * t528 + t458 * t526) * t442 + (t457 * t525 + t458 * t524) * qJD(1);
t521 = (t523 * t463 - t466 * t530) * t443 + (t463 * t530 + t523 * t466) * t442 + (t534 * t457 + t533 * t458) * qJD(1);
t517 = cos(qJ(5));
t516 = pkin(4) * t457;
t465 = cos(qJ(2));
t514 = pkin(2) * t465;
t462 = sin(qJ(2));
t512 = Icges(3,4) * t462;
t511 = Icges(3,4) * t465;
t506 = t458 * t463;
t505 = t458 * t466;
t388 = -pkin(8) * t466 + t463 * t514;
t389 = pkin(8) * t463 + t466 * t514;
t502 = qJD(2) * t466;
t504 = t388 * t456 + t389 * t502;
t452 = pkin(1) * t463 - pkin(7) * t466;
t503 = -t388 - t452;
t501 = qJD(4) * t457;
t500 = pkin(2) * qJD(2) * t462;
t492 = pkin(3) * t458 + qJ(4) * t457;
t420 = t492 * t463;
t499 = -t420 + t503;
t498 = t457 * t517;
t497 = t466 * t500;
t428 = pkin(4) * t506 + pkin(9) * t466;
t496 = -t428 + t499;
t426 = -qJD(5) * t463 + t442;
t495 = rSges(3,1) * t465 - rSges(3,2) * t462;
t494 = rSges(4,1) * t458 - rSges(4,2) * t457;
t493 = rSges(5,1) * t458 + rSges(5,3) * t457;
t491 = Icges(3,1) * t465 - t512;
t488 = -Icges(3,2) * t462 + t511;
t485 = Icges(3,5) * t465 - Icges(3,6) * t462;
t414 = -Icges(3,6) * t466 + t463 * t488;
t416 = -Icges(3,5) * t466 + t463 * t491;
t482 = t414 * t462 - t416 * t465;
t415 = Icges(3,6) * t463 + t466 * t488;
t417 = Icges(3,5) * t463 + t466 * t491;
t481 = -t415 * t462 + t417 * t465;
t445 = Icges(3,2) * t465 + t512;
t446 = Icges(3,1) * t462 + t511;
t480 = -t445 * t462 + t446 * t465;
t427 = qJD(5) * t466 + t443;
t479 = -qJD(4) * t458 + t442 * t420 + t504;
t461 = sin(qJ(5));
t422 = t457 * t461 + t458 * t517;
t439 = qJD(1) * (pkin(1) * t466 + pkin(7) * t463);
t478 = qJD(1) * t389 - t463 * t500 + t439;
t436 = pkin(3) * t457 - qJ(4) * t458;
t477 = t443 * t436 + t466 * t501 - t497;
t476 = t443 * t516 + t477;
t421 = t492 * t466;
t473 = qJD(1) * t421 + t463 * t501 + t478;
t429 = pkin(4) * t505 - pkin(9) * t463;
t472 = t442 * t428 + (-t421 - t429) * t443 + t479;
t471 = qJD(1) * t429 + (-t436 - t516) * t442 + t473;
t464 = cos(qJ(6));
t460 = sin(qJ(6));
t449 = rSges(2,1) * t466 - rSges(2,2) * t463;
t448 = rSges(2,1) * t463 + rSges(2,2) * t466;
t447 = rSges(3,1) * t462 + rSges(3,2) * t465;
t444 = Icges(3,5) * t462 + Icges(3,6) * t465;
t438 = rSges(4,1) * t457 + rSges(4,2) * t458;
t437 = rSges(5,1) * t457 - rSges(5,3) * t458;
t423 = -t458 * t461 + t498;
t419 = rSges(3,3) * t463 + t466 * t495;
t418 = -rSges(3,3) * t466 + t463 * t495;
t413 = Icges(3,3) * t463 + t466 * t485;
t412 = -Icges(3,3) * t466 + t463 * t485;
t410 = qJD(6) * t422 + qJD(1);
t409 = t422 * t466;
t408 = t461 * t505 - t466 * t498;
t407 = t422 * t463;
t406 = t461 * t506 - t463 * t498;
t405 = rSges(4,3) * t463 + t466 * t494;
t404 = rSges(5,2) * t463 + t466 * t493;
t403 = -rSges(4,3) * t466 + t463 * t494;
t402 = -rSges(5,2) * t466 + t463 * t493;
t381 = t409 * t464 - t460 * t463;
t380 = -t409 * t460 - t463 * t464;
t379 = t407 * t464 + t460 * t466;
t378 = -t407 * t460 + t464 * t466;
t377 = pkin(5) * t423 + pkin(10) * t422;
t376 = rSges(6,1) * t423 - rSges(6,2) * t422;
t375 = Icges(6,1) * t423 - Icges(6,4) * t422;
t374 = Icges(6,4) * t423 - Icges(6,2) * t422;
t373 = Icges(6,5) * t423 - Icges(6,6) * t422;
t372 = qJD(6) * t406 + t427;
t371 = qJD(6) * t408 + t426;
t370 = pkin(5) * t409 + pkin(10) * t408;
t369 = pkin(5) * t407 + pkin(10) * t406;
t368 = qJD(1) * t419 - t447 * t456 + t439;
t367 = -t447 * t502 + (-t418 - t452) * qJD(1);
t366 = (t418 * t463 + t419 * t466) * qJD(2);
t365 = rSges(6,1) * t409 - rSges(6,2) * t408 - rSges(6,3) * t463;
t364 = rSges(6,1) * t407 - rSges(6,2) * t406 + rSges(6,3) * t466;
t363 = Icges(6,1) * t409 - Icges(6,4) * t408 - Icges(6,5) * t463;
t362 = Icges(6,1) * t407 - Icges(6,4) * t406 + Icges(6,5) * t466;
t361 = Icges(6,4) * t409 - Icges(6,2) * t408 - Icges(6,6) * t463;
t360 = Icges(6,4) * t407 - Icges(6,2) * t406 + Icges(6,6) * t466;
t359 = Icges(6,5) * t409 - Icges(6,6) * t408 - Icges(6,3) * t463;
t358 = Icges(6,5) * t407 - Icges(6,6) * t406 + Icges(6,3) * t466;
t357 = rSges(7,3) * t422 + (rSges(7,1) * t464 - rSges(7,2) * t460) * t423;
t356 = Icges(7,5) * t422 + (Icges(7,1) * t464 - Icges(7,4) * t460) * t423;
t355 = Icges(7,6) * t422 + (Icges(7,4) * t464 - Icges(7,2) * t460) * t423;
t354 = Icges(7,3) * t422 + (Icges(7,5) * t464 - Icges(7,6) * t460) * t423;
t353 = rSges(7,1) * t381 + rSges(7,2) * t380 + rSges(7,3) * t408;
t352 = rSges(7,1) * t379 + rSges(7,2) * t378 + rSges(7,3) * t406;
t351 = Icges(7,1) * t381 + Icges(7,4) * t380 + Icges(7,5) * t408;
t350 = Icges(7,1) * t379 + Icges(7,4) * t378 + Icges(7,5) * t406;
t349 = Icges(7,4) * t381 + Icges(7,2) * t380 + Icges(7,6) * t408;
t348 = Icges(7,4) * t379 + Icges(7,2) * t378 + Icges(7,6) * t406;
t347 = Icges(7,5) * t381 + Icges(7,6) * t380 + Icges(7,3) * t408;
t346 = Icges(7,5) * t379 + Icges(7,6) * t378 + Icges(7,3) * t406;
t345 = qJD(1) * t405 - t438 * t442 + t478;
t344 = -t497 + t438 * t443 + (-t403 + t503) * qJD(1);
t343 = t403 * t442 - t405 * t443 + t504;
t342 = qJD(1) * t404 + (-t436 - t437) * t442 + t473;
t341 = t437 * t443 + (-t402 + t499) * qJD(1) + t477;
t340 = t402 * t442 + (-t404 - t421) * t443 + t479;
t339 = qJD(1) * t365 - t376 * t426 + t471;
t338 = t376 * t427 + (-t364 + t496) * qJD(1) + t476;
t337 = t364 * t426 - t365 * t427 + t472;
t336 = qJD(1) * t370 + t353 * t410 - t357 * t371 - t377 * t426 + t471;
t335 = -t352 * t410 + t357 * t372 + t377 * t427 + (-t369 + t496) * qJD(1) + t476;
t334 = t352 * t371 - t353 * t372 + t369 * t426 - t370 * t427 + t472;
t1 = m(7) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(6) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(5) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(4) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(3) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + t371 * ((t408 * t347 + t380 * t349 + t381 * t351) * t371 + (t346 * t408 + t348 * t380 + t350 * t381) * t372 + (t354 * t408 + t355 * t380 + t356 * t381) * t410) / 0.2e1 + t372 * ((t347 * t406 + t349 * t378 + t351 * t379) * t371 + (t406 * t346 + t378 * t348 + t379 * t350) * t372 + (t354 * t406 + t355 * t378 + t356 * t379) * t410) / 0.2e1 + t410 * ((t346 * t372 + t347 * t371 + t354 * t410) * t422 + ((-t349 * t460 + t351 * t464) * t371 + (-t348 * t460 + t350 * t464) * t372 + (-t355 * t460 + t356 * t464) * t410) * t423) / 0.2e1 + t426 * ((-t463 * t359 - t408 * t361 + t409 * t363) * t426 + (-t358 * t463 - t360 * t408 + t362 * t409) * t427 + (-t373 * t463 - t374 * t408 + t375 * t409) * qJD(1)) / 0.2e1 + t427 * ((t359 * t466 - t361 * t406 + t363 * t407) * t426 + (t466 * t358 - t406 * t360 + t407 * t362) * t427 + (t373 * t466 - t374 * t406 + t375 * t407) * qJD(1)) / 0.2e1 - ((-t466 * t444 + t463 * t480) * qJD(1) + (t466 ^ 2 * t412 + (t481 * t463 + (-t413 + t482) * t466) * t463) * qJD(2)) * t502 / 0.2e1 + ((t463 * t444 + t466 * t480) * qJD(1) + (t463 ^ 2 * t413 + (t482 * t466 + (-t412 + t481) * t463) * t466) * qJD(2)) * t456 / 0.2e1 + (t521 * t463 + t522 * t466) * t442 / 0.2e1 + (t522 * t463 - t521 * t466) * t443 / 0.2e1 + (m(2) * (t448 ^ 2 + t449 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t415 * t465 + t417 * t462) * t463 - (t414 * t465 + t416 * t462) * t466) * qJD(2) + (-t361 * t422 + t363 * t423) * t426 + (-t360 * t422 + t362 * t423) * t427 + (t457 * t527 - t458 * t529) * t443 + (t457 * t526 - t458 * t528) * t442 + (-t422 * t374 + t423 * t375 + t465 * t445 + t462 * t446 + t524 * t457 - t525 * t458) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
