% Calculate kinetic energy for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:23
% EndTime: 2019-03-09 05:55:25
% DurationCPUTime: 1.84s
% Computational Cost: add. (1826->218), mult. (1571->342), div. (0->0), fcn. (1514->10), ass. (0->125)
t475 = Icges(6,1) + Icges(7,1);
t474 = -Icges(6,4) + Icges(7,5);
t473 = Icges(7,4) + Icges(6,5);
t472 = Icges(6,2) + Icges(7,3);
t471 = -Icges(7,6) + Icges(6,6);
t470 = -Icges(6,3) - Icges(7,2);
t469 = rSges(7,1) + pkin(5);
t468 = rSges(7,3) + qJ(6);
t404 = qJ(1) + pkin(10);
t400 = sin(t404);
t401 = cos(t404);
t409 = cos(qJ(5));
t405 = qJ(3) + qJ(4);
t403 = cos(t405);
t406 = sin(qJ(5));
t444 = t403 * t406;
t374 = t400 * t444 + t401 * t409;
t443 = t403 * t409;
t375 = t400 * t443 - t401 * t406;
t402 = sin(t405);
t446 = t400 * t402;
t467 = t472 * t374 + t474 * t375 - t471 * t446;
t376 = -t400 * t409 + t401 * t444;
t377 = t400 * t406 + t401 * t443;
t445 = t401 * t402;
t466 = t472 * t376 + t474 * t377 - t471 * t445;
t465 = -t471 * t374 + t473 * t375 - t470 * t446;
t464 = -t471 * t376 + t473 * t377 - t470 * t445;
t463 = t474 * t374 + t475 * t375 + t473 * t446;
t462 = t474 * t376 + t475 * t377 + t473 * t445;
t461 = t471 * t403 + (t472 * t406 + t474 * t409) * t402;
t460 = t470 * t403 + (-t471 * t406 + t473 * t409) * t402;
t459 = -t473 * t403 + (t474 * t406 + t475 * t409) * t402;
t408 = sin(qJ(1));
t454 = pkin(1) * t408;
t410 = cos(qJ(3));
t452 = pkin(3) * t410;
t407 = sin(qJ(3));
t450 = Icges(4,4) * t407;
t449 = Icges(4,4) * t410;
t448 = Icges(5,4) * t402;
t447 = Icges(5,4) * t403;
t442 = rSges(7,2) * t446 + t468 * t374 + t469 * t375;
t441 = rSges(7,2) * t445 + t468 * t376 + t469 * t377;
t440 = -rSges(7,2) * t403 + (t468 * t406 + t469 * t409) * t402;
t411 = cos(qJ(1));
t399 = qJD(1) * t411 * pkin(1);
t439 = qJD(1) * (pkin(2) * t401 + pkin(7) * t400) + t399;
t397 = qJD(3) * t400;
t380 = qJD(4) * t400 + t397;
t438 = qJD(3) * t401;
t437 = qJD(5) * t402;
t436 = pkin(3) * qJD(3) * t407;
t340 = -pkin(8) * t401 + t400 * t452;
t341 = pkin(8) * t400 + t401 * t452;
t435 = t340 * t397 + t341 * t438 + qJD(2);
t434 = -pkin(2) * t400 + pkin(7) * t401 - t454;
t381 = (-qJD(3) - qJD(4)) * t401;
t433 = t401 * t436;
t432 = -t340 + t434;
t431 = pkin(4) * t403 + pkin(9) * t402;
t430 = rSges(4,1) * t410 - rSges(4,2) * t407;
t429 = rSges(5,1) * t403 - rSges(5,2) * t402;
t428 = Icges(4,1) * t410 - t450;
t427 = Icges(5,1) * t403 - t448;
t426 = -Icges(4,2) * t407 + t449;
t425 = -Icges(5,2) * t402 + t447;
t424 = Icges(4,5) * t410 - Icges(4,6) * t407;
t423 = Icges(5,5) * t403 - Icges(5,6) * t402;
t353 = -Icges(4,6) * t401 + t400 * t426;
t355 = -Icges(4,5) * t401 + t400 * t428;
t422 = t353 * t407 - t355 * t410;
t354 = Icges(4,6) * t400 + t401 * t426;
t356 = Icges(4,5) * t400 + t401 * t428;
t421 = -t354 * t407 + t356 * t410;
t389 = Icges(4,2) * t410 + t450;
t390 = Icges(4,1) * t407 + t449;
t420 = -t389 * t407 + t390 * t410;
t370 = t431 * t400;
t371 = t431 * t401;
t419 = t380 * t370 - t371 * t381 + t435;
t418 = qJD(1) * t341 - t400 * t436 + t439;
t417 = (Icges(5,5) * t402 + Icges(5,6) * t403) * qJD(1) + (-Icges(5,3) * t401 + t400 * t423) * t381 + (Icges(5,3) * t400 + t401 * t423) * t380;
t387 = pkin(4) * t402 - pkin(9) * t403;
t416 = qJD(1) * t371 - t380 * t387 + t418;
t415 = t381 * t387 + (-t370 + t432) * qJD(1) - t433;
t344 = -Icges(5,6) * t401 + t400 * t425;
t345 = Icges(5,6) * t400 + t401 * t425;
t346 = -Icges(5,5) * t401 + t400 * t427;
t347 = Icges(5,5) * t400 + t401 * t427;
t384 = Icges(5,2) * t403 + t448;
t385 = Icges(5,1) * t402 + t447;
t414 = (-t345 * t402 + t347 * t403) * t380 + (-t344 * t402 + t346 * t403) * t381 + (-t384 * t402 + t385 * t403) * qJD(1);
t394 = -qJD(5) * t403 + qJD(1);
t393 = rSges(2,1) * t411 - rSges(2,2) * t408;
t392 = rSges(2,1) * t408 + rSges(2,2) * t411;
t391 = rSges(4,1) * t407 + rSges(4,2) * t410;
t388 = Icges(4,5) * t407 + Icges(4,6) * t410;
t386 = rSges(5,1) * t402 + rSges(5,2) * t403;
t373 = t399 + qJD(1) * (rSges(3,1) * t401 - rSges(3,2) * t400);
t372 = (-rSges(3,1) * t400 - rSges(3,2) * t401 - t454) * qJD(1);
t368 = -rSges(6,3) * t403 + (rSges(6,1) * t409 - rSges(6,2) * t406) * t402;
t360 = rSges(4,3) * t400 + t401 * t430;
t359 = -rSges(4,3) * t401 + t400 * t430;
t358 = t400 * t437 + t381;
t357 = t401 * t437 + t380;
t352 = Icges(4,3) * t400 + t401 * t424;
t351 = -Icges(4,3) * t401 + t400 * t424;
t349 = rSges(5,3) * t400 + t401 * t429;
t348 = -rSges(5,3) * t401 + t400 * t429;
t333 = rSges(6,1) * t377 - rSges(6,2) * t376 + rSges(6,3) * t445;
t331 = rSges(6,1) * t375 - rSges(6,2) * t374 + rSges(6,3) * t446;
t317 = qJD(1) * t360 - t391 * t397 + t439;
t316 = -t391 * t438 + (-t359 + t434) * qJD(1);
t315 = qJD(2) + (t359 * t400 + t360 * t401) * qJD(3);
t314 = qJD(1) * t349 - t380 * t386 + t418;
t313 = -t433 + t381 * t386 + (-t348 + t432) * qJD(1);
t312 = t348 * t380 - t349 * t381 + t435;
t311 = t333 * t394 - t357 * t368 + t416;
t310 = -t331 * t394 + t358 * t368 + t415;
t309 = t331 * t357 - t333 * t358 + t419;
t308 = qJD(6) * t374 - t357 * t440 + t394 * t441 + t416;
t307 = qJD(6) * t376 + t358 * t440 - t394 * t442 + t415;
t306 = qJD(6) * t402 * t406 + t357 * t442 - t358 * t441 + t419;
t1 = m(7) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + m(6) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(4) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(5) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 - ((-t401 * t388 + t400 * t420) * qJD(1) + (t401 ^ 2 * t351 + (t421 * t400 + (-t352 + t422) * t401) * t400) * qJD(3)) * t438 / 0.2e1 + t380 * (t417 * t400 + t414 * t401) / 0.2e1 + t381 * (t414 * t400 - t417 * t401) / 0.2e1 + ((t400 * t388 + t401 * t420) * qJD(1) + (t400 ^ 2 * t352 + (t422 * t401 + (-t351 + t421) * t400) * t401) * qJD(3)) * t397 / 0.2e1 + ((t376 * t461 + t377 * t459 + t445 * t460) * t394 + (t376 * t467 + t463 * t377 + t465 * t445) * t358 + (t466 * t376 + t462 * t377 + t464 * t445) * t357) * t357 / 0.2e1 + ((t374 * t461 + t375 * t459 + t446 * t460) * t394 + (t467 * t374 + t463 * t375 + t465 * t446) * t358 + (t374 * t466 + t375 * t462 + t446 * t464) * t357) * t358 / 0.2e1 + ((-t357 * t464 - t358 * t465 - t394 * t460) * t403 + ((t406 * t461 + t409 * t459) * t394 + (t406 * t467 + t463 * t409) * t358 + (t406 * t466 + t409 * t462) * t357) * t402) * t394 / 0.2e1 + (((t354 * t410 + t356 * t407) * t400 - (t353 * t410 + t355 * t407) * t401) * qJD(3) + (t345 * t403 + t347 * t402) * t380 + (t344 * t403 + t346 * t402) * t381 + (t403 * t384 + t402 * t385 + t410 * t389 + t407 * t390) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t392 ^ 2 + t393 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
