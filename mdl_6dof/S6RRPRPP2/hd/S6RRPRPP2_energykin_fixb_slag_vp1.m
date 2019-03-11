% Calculate kinetic energy for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:41
% EndTime: 2019-03-09 09:49:43
% DurationCPUTime: 2.10s
% Computational Cost: add. (1547->234), mult. (2108->343), div. (0->0), fcn. (2089->8), ass. (0->131)
t523 = Icges(3,3) + Icges(4,3);
t433 = qJ(2) + pkin(9);
t430 = sin(t433);
t431 = cos(t433);
t436 = sin(qJ(2));
t439 = cos(qJ(2));
t522 = Icges(3,5) * t439 + Icges(4,5) * t431 - Icges(3,6) * t436 - Icges(4,6) * t430;
t521 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t520 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t519 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t518 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t517 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t516 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t515 = rSges(7,1) + pkin(5);
t514 = rSges(7,3) + qJ(6);
t437 = sin(qJ(1));
t438 = cos(qJ(4));
t478 = t437 * t438;
t435 = sin(qJ(4));
t440 = cos(qJ(1));
t479 = t435 * t440;
t408 = t431 * t479 - t478;
t477 = t438 * t440;
t480 = t435 * t437;
t409 = t431 * t477 + t480;
t361 = pkin(4) * t409 + qJ(5) * t408;
t406 = t431 * t480 + t477;
t427 = -qJD(4) * t431 + qJD(1);
t513 = qJD(5) * t406 + t427 * t361;
t401 = (pkin(4) * t438 + qJ(5) * t435) * t430;
t469 = qJD(4) * t430;
t470 = qJD(2) * t440;
t412 = t437 * t469 - t470;
t512 = qJD(5) * t408 + t412 * t401;
t511 = t522 * t437 - t440 * t523;
t510 = t437 * t523 + t522 * t440;
t509 = Icges(3,5) * t436 + Icges(4,5) * t430 + Icges(3,6) * t439 + Icges(4,6) * t431;
t484 = Icges(4,4) * t430;
t414 = Icges(4,2) * t431 + t484;
t483 = Icges(4,4) * t431;
t415 = Icges(4,1) * t430 + t483;
t486 = Icges(3,4) * t436;
t420 = Icges(3,2) * t439 + t486;
t485 = Icges(3,4) * t439;
t421 = Icges(3,1) * t436 + t485;
t508 = -t414 * t430 + t415 * t431 - t420 * t436 + t421 * t439;
t456 = -Icges(4,2) * t430 + t483;
t384 = -Icges(4,6) * t440 + t437 * t456;
t458 = Icges(4,1) * t431 - t484;
t386 = -Icges(4,5) * t440 + t437 * t458;
t457 = -Icges(3,2) * t436 + t485;
t395 = -Icges(3,6) * t440 + t437 * t457;
t459 = Icges(3,1) * t439 - t486;
t397 = -Icges(3,5) * t440 + t437 * t459;
t507 = t384 * t430 - t386 * t431 + t395 * t436 - t397 * t439;
t385 = Icges(4,6) * t437 + t440 * t456;
t387 = Icges(4,5) * t437 + t440 * t458;
t396 = Icges(3,6) * t437 + t440 * t457;
t398 = Icges(3,5) * t437 + t440 * t459;
t506 = -t385 * t430 + t387 * t431 - t396 * t436 + t398 * t439;
t407 = t431 * t478 - t479;
t482 = t430 * t437;
t505 = t406 * t517 + t407 * t519 - t482 * t516;
t481 = t430 * t440;
t504 = t408 * t517 + t409 * t519 - t481 * t516;
t503 = t406 * t518 + t407 * t520 - t482 * t517;
t502 = t408 * t518 + t409 * t520 - t481 * t517;
t501 = t406 * t520 + t407 * t521 - t482 * t519;
t500 = t408 * t520 + t409 * t521 - t481 * t519;
t499 = t516 * t431 + (t435 * t517 + t438 * t519) * t430;
t498 = t517 * t431 + (t435 * t518 + t438 * t520) * t430;
t497 = t519 * t431 + (t435 * t520 + t438 * t521) * t430;
t490 = pkin(2) * t436;
t488 = pkin(2) * t439;
t476 = rSges(7,2) * t406 + t515 * t407 - t514 * t482;
t475 = rSges(7,2) * t408 + t515 * t409 - t514 * t481;
t380 = -qJ(3) * t440 + t437 * t488;
t381 = qJ(3) * t437 + t440 * t488;
t471 = qJD(2) * t437;
t474 = t380 * t471 + t381 * t470;
t473 = t514 * t431 + (rSges(7,2) * t435 + t515 * t438) * t430;
t426 = pkin(1) * t437 - pkin(7) * t440;
t472 = -t380 - t426;
t465 = pkin(3) * t431 + pkin(8) * t430;
t404 = t465 * t437;
t405 = t465 * t440;
t466 = t404 * t471 + t405 * t470 + t474;
t418 = qJD(1) * (pkin(1) * t440 + pkin(7) * t437);
t464 = qJD(1) * t381 - qJD(3) * t440 + t418;
t463 = rSges(3,1) * t439 - rSges(3,2) * t436;
t462 = rSges(4,1) * t431 - rSges(4,2) * t430;
t461 = qJD(2) * (-rSges(4,1) * t430 - rSges(4,2) * t431 - t490);
t460 = (-pkin(3) * t430 + pkin(8) * t431 - t490) * qJD(2);
t447 = qJD(1) * t405 + t464;
t360 = pkin(4) * t407 + qJ(5) * t406;
t411 = t440 * t469 + t471;
t446 = qJD(5) * t430 * t435 + t411 * t360 + t466;
t432 = qJD(3) * t437;
t445 = t432 + (-t404 + t472) * qJD(1);
t444 = -qJD(6) * t430 + t460;
t443 = t437 * t460 + t447;
t442 = t440 * t460 + t445;
t425 = rSges(2,1) * t440 - rSges(2,2) * t437;
t424 = rSges(2,1) * t437 + rSges(2,2) * t440;
t423 = rSges(3,1) * t436 + rSges(3,2) * t439;
t403 = rSges(3,3) * t437 + t440 * t463;
t402 = -rSges(3,3) * t440 + t437 * t463;
t389 = rSges(4,3) * t437 + t440 * t462;
t388 = -rSges(4,3) * t440 + t437 * t462;
t379 = -rSges(5,3) * t431 + (rSges(5,1) * t438 - rSges(5,2) * t435) * t430;
t378 = -rSges(6,2) * t431 + (rSges(6,1) * t438 + rSges(6,3) * t435) * t430;
t359 = qJD(1) * t403 - t423 * t471 + t418;
t358 = -t423 * t470 + (-t402 - t426) * qJD(1);
t357 = (t402 * t437 + t403 * t440) * qJD(2);
t356 = rSges(5,1) * t409 - rSges(5,2) * t408 + rSges(5,3) * t481;
t355 = rSges(6,1) * t409 + rSges(6,2) * t481 + rSges(6,3) * t408;
t353 = rSges(5,1) * t407 - rSges(5,2) * t406 + rSges(5,3) * t482;
t352 = rSges(6,1) * t407 + rSges(6,2) * t482 + rSges(6,3) * t406;
t330 = qJD(1) * t389 + t437 * t461 + t464;
t329 = t432 + t440 * t461 + (-t388 + t472) * qJD(1);
t328 = (t388 * t437 + t389 * t440) * qJD(2) + t474;
t327 = t356 * t427 - t379 * t411 + t443;
t326 = -t353 * t427 + t379 * t412 + t442;
t325 = t353 * t411 - t356 * t412 + t466;
t324 = t355 * t427 + (-t378 - t401) * t411 + t443 + t513;
t323 = t378 * t412 + (-t352 - t360) * t427 + t442 + t512;
t322 = t352 * t411 + (-t355 - t361) * t412 + t446;
t321 = t475 * t427 + t444 * t437 + (-t401 - t473) * t411 + t447 + t513;
t320 = t473 * t412 + t444 * t440 + (-t360 - t476) * t427 + t445 + t512;
t319 = qJD(6) * t431 + t476 * t411 + (-t361 - t475) * t412 + t446;
t1 = m(6) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(7) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(5) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(4) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(3) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t424 ^ 2 + t425 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t384 * t431 - t386 * t430 - t395 * t439 - t397 * t436) * t440 + (t385 * t431 + t387 * t430 + t396 * t439 + t398 * t436) * t437) * qJD(2) + (t414 * t431 + t415 * t430 + t420 * t439 + t421 * t436) * qJD(1)) * qJD(1) / 0.2e1 + ((t510 * t437 ^ 2 + (t507 * t440 + (t506 - t511) * t437) * t440) * qJD(2) + (t437 * t509 + t440 * t508) * qJD(1)) * t471 / 0.2e1 - ((t511 * t440 ^ 2 + (t506 * t437 + (t507 - t510) * t440) * t437) * qJD(2) + (t437 * t508 - t440 * t509) * qJD(1)) * t470 / 0.2e1 + ((t408 * t498 + t409 * t497 - t481 * t499) * t427 + (t408 * t503 + t409 * t501 - t481 * t505) * t412 + (t408 * t502 + t409 * t500 - t481 * t504) * t411) * t411 / 0.2e1 + ((t406 * t498 + t407 * t497 - t482 * t499) * t427 + (t406 * t503 + t407 * t501 - t482 * t505) * t412 + (t406 * t502 + t407 * t500 - t482 * t504) * t411) * t412 / 0.2e1 + ((t411 * t504 + t412 * t505 + t427 * t499) * t431 + ((t435 * t498 + t438 * t497) * t427 + (t435 * t503 + t438 * t501) * t412 + (t435 * t502 + t438 * t500) * t411) * t430) * t427 / 0.2e1;
T  = t1;
