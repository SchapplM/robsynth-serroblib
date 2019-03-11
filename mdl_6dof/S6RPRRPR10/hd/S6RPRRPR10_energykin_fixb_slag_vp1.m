% Calculate kinetic energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:21
% EndTime: 2019-03-09 05:34:23
% DurationCPUTime: 2.20s
% Computational Cost: add. (924->255), mult. (2206->391), div. (0->0), fcn. (2360->8), ass. (0->125)
t469 = Icges(5,1) + Icges(6,1);
t468 = Icges(5,4) - Icges(6,5);
t467 = Icges(6,4) + Icges(5,5);
t466 = Icges(5,2) + Icges(6,3);
t465 = Icges(6,6) - Icges(5,6);
t464 = Icges(5,3) + Icges(6,2);
t419 = sin(qJ(3));
t422 = cos(qJ(4));
t424 = cos(qJ(1));
t443 = t424 * t422;
t418 = sin(qJ(4));
t420 = sin(qJ(1));
t448 = t418 * t420;
t387 = t419 * t448 - t443;
t446 = t420 * t422;
t447 = t418 * t424;
t388 = t419 * t446 + t447;
t423 = cos(qJ(3));
t445 = t420 * t423;
t463 = t466 * t387 - t468 * t388 - t465 * t445;
t389 = t419 * t447 + t446;
t390 = -t419 * t443 + t448;
t444 = t423 * t424;
t462 = -t466 * t389 - t468 * t390 + t465 * t444;
t461 = t465 * t387 + t467 * t388 - t464 * t445;
t460 = -t465 * t389 + t467 * t390 + t464 * t444;
t459 = -t468 * t387 + t469 * t388 - t467 * t445;
t458 = t468 * t389 + t469 * t390 + t467 * t444;
t457 = (t466 * t418 - t468 * t422) * t423 + t465 * t419;
t456 = (t465 * t418 + t467 * t422) * t423 + t464 * t419;
t455 = (-t468 * t418 + t469 * t422) * t423 + t467 * t419;
t450 = Icges(4,4) * t419;
t449 = Icges(4,4) * t423;
t399 = qJD(1) * (pkin(1) * t424 + qJ(2) * t420);
t442 = qJD(1) * t424 * pkin(7) + t399;
t414 = qJD(3) * t420;
t441 = qJD(4) * t423;
t395 = t424 * t441 + t414;
t415 = qJD(3) * t424;
t410 = qJD(4) * t419 + qJD(1);
t403 = pkin(1) * t420 - qJ(2) * t424;
t440 = -pkin(7) * t420 - t403;
t439 = pkin(3) * t419 - pkin(8) * t423;
t392 = t439 * t420;
t393 = t439 * t424;
t438 = -t392 * t414 - t393 * t415;
t437 = rSges(4,1) * t419 + rSges(4,2) * t423;
t436 = Icges(4,1) * t419 + t449;
t435 = Icges(4,2) * t423 + t450;
t434 = Icges(4,5) * t419 + Icges(4,6) * t423;
t371 = Icges(4,6) * t424 + t420 * t435;
t375 = Icges(4,5) * t424 + t420 * t436;
t433 = -t371 * t423 - t375 * t419;
t372 = Icges(4,6) * t420 - t424 * t435;
t376 = Icges(4,5) * t420 - t424 * t436;
t432 = t372 * t423 + t376 * t419;
t401 = -Icges(4,2) * t419 + t449;
t402 = Icges(4,1) * t423 - t450;
t431 = t401 * t423 + t402 * t419;
t357 = pkin(4) * t390 - qJ(5) * t389;
t396 = -t420 * t441 + t415;
t430 = qJD(5) * t423 * t418 + t396 * t357 + t438;
t407 = pkin(3) * t423 + pkin(8) * t419;
t416 = qJD(2) * t420;
t429 = t407 * t414 + t416 + (t393 + t440) * qJD(1);
t428 = qJD(1) * t392 + (-qJD(3) * t407 - qJD(2)) * t424 + t442;
t391 = (pkin(4) * t422 + qJ(5) * t418) * t423;
t427 = qJD(5) * t387 + t395 * t391 + t429;
t356 = pkin(4) * t388 + qJ(5) * t387;
t426 = -qJD(5) * t389 + t410 * t356 + t428;
t421 = cos(qJ(6));
t417 = sin(qJ(6));
t406 = rSges(2,1) * t424 - rSges(2,2) * t420;
t405 = rSges(4,1) * t423 - rSges(4,2) * t419;
t404 = rSges(2,1) * t420 + rSges(2,2) * t424;
t400 = Icges(4,5) * t423 - Icges(4,6) * t419;
t398 = -qJD(6) * t419 + t410;
t397 = pkin(5) * t422 * t423 - pkin(9) * t419;
t383 = (t417 * t418 + t421 * t422) * t423;
t382 = (-t417 * t422 + t418 * t421) * t423;
t380 = rSges(4,3) * t420 - t424 * t437;
t379 = rSges(5,3) * t419 + (rSges(5,1) * t422 - rSges(5,2) * t418) * t423;
t378 = rSges(6,2) * t419 + (rSges(6,1) * t422 + rSges(6,3) * t418) * t423;
t377 = rSges(4,3) * t424 + t420 * t437;
t368 = Icges(4,3) * t420 - t424 * t434;
t367 = Icges(4,3) * t424 + t420 * t434;
t364 = t415 + (-qJD(4) + qJD(6)) * t445;
t363 = -qJD(6) * t444 + t395;
t362 = pkin(5) * t390 - pkin(9) * t444;
t361 = pkin(5) * t388 + pkin(9) * t445;
t359 = t399 - qJD(2) * t424 + qJD(1) * (-rSges(3,2) * t424 + rSges(3,3) * t420);
t358 = t416 + (rSges(3,2) * t420 + rSges(3,3) * t424 - t403) * qJD(1);
t355 = -t389 * t417 + t390 * t421;
t354 = -t389 * t421 - t390 * t417;
t353 = t387 * t417 + t388 * t421;
t352 = t387 * t421 - t388 * t417;
t351 = rSges(5,1) * t390 + rSges(5,2) * t389 + rSges(5,3) * t444;
t350 = rSges(6,1) * t390 + rSges(6,2) * t444 - rSges(6,3) * t389;
t349 = rSges(5,1) * t388 - rSges(5,2) * t387 - rSges(5,3) * t445;
t348 = rSges(6,1) * t388 - rSges(6,2) * t445 + rSges(6,3) * t387;
t334 = rSges(7,1) * t383 + rSges(7,2) * t382 - rSges(7,3) * t419;
t333 = Icges(7,1) * t383 + Icges(7,4) * t382 - Icges(7,5) * t419;
t332 = Icges(7,4) * t383 + Icges(7,2) * t382 - Icges(7,6) * t419;
t331 = Icges(7,5) * t383 + Icges(7,6) * t382 - Icges(7,3) * t419;
t329 = (-t377 * t420 + t380 * t424) * qJD(3);
t328 = qJD(1) * t377 + (-qJD(3) * t405 - qJD(2)) * t424 + t442;
t327 = t405 * t414 + t416 + (-t380 + t440) * qJD(1);
t326 = rSges(7,1) * t355 + rSges(7,2) * t354 - rSges(7,3) * t444;
t325 = rSges(7,1) * t353 + rSges(7,2) * t352 + rSges(7,3) * t445;
t324 = Icges(7,1) * t355 + Icges(7,4) * t354 - Icges(7,5) * t444;
t323 = Icges(7,1) * t353 + Icges(7,4) * t352 + Icges(7,5) * t445;
t322 = Icges(7,4) * t355 + Icges(7,2) * t354 - Icges(7,6) * t444;
t321 = Icges(7,4) * t353 + Icges(7,2) * t352 + Icges(7,6) * t445;
t320 = Icges(7,5) * t355 + Icges(7,6) * t354 - Icges(7,3) * t444;
t319 = Icges(7,5) * t353 + Icges(7,6) * t352 + Icges(7,3) * t445;
t318 = t349 * t410 - t379 * t396 + t428;
t317 = -t351 * t410 + t379 * t395 + t429;
t316 = -t349 * t395 + t351 * t396 + t438;
t315 = t348 * t410 + (-t378 - t391) * t396 + t426;
t314 = t378 * t395 + (-t350 - t357) * t410 + t427;
t313 = t350 * t396 + (-t348 - t356) * t395 + t430;
t312 = t325 * t398 - t334 * t364 + t361 * t410 + (-t391 - t397) * t396 + t426;
t311 = -t326 * t398 + t334 * t363 + t395 * t397 + (-t357 - t362) * t410 + t427;
t310 = -t325 * t363 + t326 * t364 + t362 * t396 + (-t356 - t361) * t395 + t430;
t1 = ((t424 * t400 + t420 * t431) * qJD(1) + (t424 ^ 2 * t367 + (t432 * t420 + (t368 - t433) * t424) * t420) * qJD(3)) * t415 / 0.2e1 + t364 * ((t319 * t445 + t352 * t321 + t353 * t323) * t364 + (t320 * t445 + t322 * t352 + t324 * t353) * t363 + (t331 * t445 + t332 * t352 + t333 * t353) * t398) / 0.2e1 + m(3) * (t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(6) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(7) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + m(5) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(4) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + ((t420 * t400 - t424 * t431) * qJD(1) + (t420 ^ 2 * t368 + (t433 * t424 + (t367 - t432) * t420) * t424) * qJD(3)) * t414 / 0.2e1 + qJD(1) * ((-t419 * t401 + t423 * t402) * qJD(1) + ((-t371 * t419 + t375 * t423) * t424 + (-t372 * t419 + t376 * t423) * t420) * qJD(3)) / 0.2e1 + t363 * ((-t319 * t444 + t321 * t354 + t323 * t355) * t364 + (-t320 * t444 + t354 * t322 + t355 * t324) * t363 + (-t331 * t444 + t332 * t354 + t333 * t355) * t398) / 0.2e1 + t398 * ((-t319 * t419 + t321 * t382 + t323 * t383) * t364 + (-t320 * t419 + t322 * t382 + t324 * t383) * t363 + (-t419 * t331 + t382 * t332 + t383 * t333) * t398) / 0.2e1 + ((-t457 * t389 + t455 * t390 + t456 * t444) * t410 + (-t463 * t389 + t459 * t390 + t461 * t444) * t396 + (-t462 * t389 + t458 * t390 + t460 * t444) * t395) * t395 / 0.2e1 + ((t457 * t387 + t455 * t388 - t456 * t445) * t410 + (t463 * t387 + t459 * t388 - t461 * t445) * t396 + (t462 * t387 + t458 * t388 - t460 * t445) * t395) * t396 / 0.2e1 + (((t457 * t418 + t455 * t422) * t410 + (t463 * t418 + t459 * t422) * t396 + (t462 * t418 + t458 * t422) * t395) * t423 + (t460 * t395 + t461 * t396 + t456 * t410) * t419) * t410 / 0.2e1 + (m(2) * (t404 ^ 2 + t406 ^ 2) + Icges(3,1) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
