% Calculate kinetic energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:07
% EndTime: 2019-12-31 22:08:09
% DurationCPUTime: 2.06s
% Computational Cost: add. (1826->255), mult. (4551->396), div. (0->0), fcn. (5547->10), ass. (0->122)
t472 = Icges(5,1) + Icges(6,1);
t471 = Icges(5,4) + Icges(6,4);
t470 = Icges(5,5) + Icges(6,5);
t469 = Icges(5,2) + Icges(6,2);
t468 = Icges(5,6) + Icges(6,6);
t467 = Icges(5,3) + Icges(6,3);
t466 = rSges(6,3) + qJ(5);
t420 = sin(qJ(2));
t421 = sin(qJ(1));
t423 = cos(qJ(2));
t424 = cos(qJ(1));
t448 = cos(pkin(5));
t432 = t424 * t448;
t399 = t420 * t421 - t423 * t432;
t400 = t420 * t432 + t421 * t423;
t416 = sin(pkin(5));
t443 = t416 * t424;
t368 = Icges(3,5) * t400 - Icges(3,6) * t399 - Icges(3,3) * t443;
t433 = t421 * t448;
t401 = t424 * t420 + t423 * t433;
t402 = -t420 * t433 + t424 * t423;
t445 = t416 * t421;
t369 = Icges(3,5) * t402 - Icges(3,6) * t401 + Icges(3,3) * t445;
t465 = (t368 * t424 - t369 * t421) * t416;
t419 = sin(qJ(3));
t451 = cos(qJ(3));
t386 = t400 * t451 - t419 * t443;
t418 = sin(qJ(4));
t422 = cos(qJ(4));
t360 = -t386 * t418 + t399 * t422;
t447 = t399 * t418;
t361 = t386 * t422 + t447;
t435 = t416 * t451;
t385 = t400 * t419 + t424 * t435;
t464 = t468 * t360 + t470 * t361 + t467 * t385;
t388 = t402 * t451 + t419 * t445;
t362 = -t388 * t418 + t401 * t422;
t446 = t401 * t418;
t363 = t388 * t422 + t446;
t387 = t402 * t419 - t421 * t435;
t463 = t468 * t362 + t470 * t363 + t467 * t387;
t462 = t469 * t360 + t471 * t361 + t468 * t385;
t461 = t469 * t362 + t471 * t363 + t468 * t387;
t460 = t471 * t360 + t472 * t361 + t470 * t385;
t459 = t471 * t362 + t472 * t363 + t470 * t387;
t398 = t419 * t448 + t420 * t435;
t444 = t416 * t423;
t383 = -t398 * t418 - t422 * t444;
t436 = t418 * t444;
t384 = t398 * t422 - t436;
t397 = t416 * t419 * t420 - t448 * t451;
t458 = t468 * t383 + t470 * t384 + t467 * t397;
t457 = t469 * t383 + t471 * t384 + t468 * t397;
t456 = t471 * t383 + t472 * t384 + t470 * t397;
t450 = pkin(4) * t422;
t442 = rSges(6,1) * t361 + rSges(6,2) * t360 + pkin(4) * t447 + t466 * t385 + t386 * t450;
t441 = rSges(6,1) * t363 + rSges(6,2) * t362 + pkin(4) * t446 + t466 * t387 + t388 * t450;
t440 = rSges(6,1) * t384 + rSges(6,2) * t383 - pkin(4) * t436 + t466 * t397 + t398 * t450;
t380 = pkin(2) * t400 + pkin(8) * t399;
t381 = pkin(2) * t402 + pkin(8) * t401;
t437 = qJD(2) * t416;
t411 = t421 * t437;
t434 = t424 * t437;
t439 = t380 * t411 + t381 * t434;
t389 = qJD(3) * t401 + t411;
t438 = qJD(1) * (pkin(1) * t421 - pkin(7) * t443);
t412 = qJD(2) * t448 + qJD(1);
t390 = qJD(3) * t399 - t434;
t356 = pkin(3) * t386 + pkin(9) * t385;
t357 = pkin(3) * t388 + pkin(9) * t387;
t430 = t389 * t356 - t357 * t390 + t439;
t404 = -qJD(3) * t444 + t412;
t403 = (pkin(2) * t420 - pkin(8) * t423) * t416;
t405 = qJD(1) * (pkin(1) * t424 + pkin(7) * t445);
t429 = t412 * t381 - t403 * t411 + t405;
t428 = -t380 * t412 - t403 * t434 - t438;
t379 = pkin(3) * t398 + pkin(9) * t397;
t427 = t404 * t357 - t379 * t389 + t429;
t426 = -t356 * t404 + t390 * t379 + t428;
t408 = rSges(2,1) * t424 - rSges(2,2) * t421;
t407 = rSges(2,1) * t421 + rSges(2,2) * t424;
t394 = t448 * rSges(3,3) + (rSges(3,1) * t420 + rSges(3,2) * t423) * t416;
t393 = Icges(3,5) * t448 + (Icges(3,1) * t420 + Icges(3,4) * t423) * t416;
t392 = Icges(3,6) * t448 + (Icges(3,4) * t420 + Icges(3,2) * t423) * t416;
t391 = Icges(3,3) * t448 + (Icges(3,5) * t420 + Icges(3,6) * t423) * t416;
t382 = qJD(4) * t397 + t404;
t376 = rSges(3,1) * t402 - rSges(3,2) * t401 + rSges(3,3) * t445;
t375 = rSges(3,1) * t400 - rSges(3,2) * t399 - rSges(3,3) * t443;
t373 = Icges(3,1) * t402 - Icges(3,4) * t401 + Icges(3,5) * t445;
t372 = Icges(3,1) * t400 - Icges(3,4) * t399 - Icges(3,5) * t443;
t371 = Icges(3,4) * t402 - Icges(3,2) * t401 + Icges(3,6) * t445;
t370 = Icges(3,4) * t400 - Icges(3,2) * t399 - Icges(3,6) * t443;
t367 = rSges(4,1) * t398 - rSges(4,2) * t397 - rSges(4,3) * t444;
t366 = Icges(4,1) * t398 - Icges(4,4) * t397 - Icges(4,5) * t444;
t365 = Icges(4,4) * t398 - Icges(4,2) * t397 - Icges(4,6) * t444;
t364 = Icges(4,5) * t398 - Icges(4,6) * t397 - Icges(4,3) * t444;
t359 = qJD(4) * t385 + t390;
t358 = qJD(4) * t387 + t389;
t353 = rSges(4,1) * t388 - rSges(4,2) * t387 + rSges(4,3) * t401;
t352 = rSges(4,1) * t386 - rSges(4,2) * t385 + rSges(4,3) * t399;
t351 = Icges(4,1) * t388 - Icges(4,4) * t387 + Icges(4,5) * t401;
t350 = Icges(4,1) * t386 - Icges(4,4) * t385 + Icges(4,5) * t399;
t349 = Icges(4,4) * t388 - Icges(4,2) * t387 + Icges(4,6) * t401;
t348 = Icges(4,4) * t386 - Icges(4,2) * t385 + Icges(4,6) * t399;
t347 = Icges(4,5) * t388 - Icges(4,6) * t387 + Icges(4,3) * t401;
t346 = Icges(4,5) * t386 - Icges(4,6) * t385 + Icges(4,3) * t399;
t345 = rSges(5,1) * t384 + rSges(5,2) * t383 + rSges(5,3) * t397;
t335 = t376 * t412 - t394 * t411 + t405;
t334 = -t375 * t412 - t394 * t434 - t438;
t333 = (t375 * t421 + t376 * t424) * t437;
t332 = rSges(5,1) * t363 + rSges(5,2) * t362 + rSges(5,3) * t387;
t330 = rSges(5,1) * t361 + rSges(5,2) * t360 + rSges(5,3) * t385;
t314 = t353 * t404 - t367 * t389 + t429;
t313 = -t352 * t404 + t367 * t390 + t428;
t312 = t352 * t389 - t353 * t390 + t439;
t311 = t332 * t382 - t345 * t358 + t427;
t310 = -t330 * t382 + t345 * t359 + t426;
t309 = t330 * t358 - t332 * t359 + t430;
t308 = qJD(5) * t385 - t358 * t440 + t382 * t441 + t427;
t307 = qJD(5) * t387 + t359 * t440 - t382 * t442 + t426;
t306 = qJD(5) * t397 + t358 * t442 - t359 * t441 + t430;
t1 = m(3) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + ((t391 * t445 - t392 * t401 + t393 * t402) * t412 + (-(-t370 * t401 + t372 * t402) * t424 + (-t401 * t371 + t402 * t373 - t465) * t421) * t437) * t411 / 0.2e1 - ((-t391 * t443 - t392 * t399 + t393 * t400) * t412 + ((-t371 * t399 + t373 * t400) * t421 + (t399 * t370 - t400 * t372 + t465) * t424) * t437) * t434 / 0.2e1 + t412 * ((t448 * t369 + (t371 * t423 + t373 * t420) * t416) * t411 - (t448 * t368 + (t370 * t423 + t372 * t420) * t416) * t434 + (t448 * t391 + (t392 * t423 + t393 * t420) * t416) * t412) / 0.2e1 + m(4) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + t389 * ((t401 * t347 - t387 * t349 + t388 * t351) * t389 + (t346 * t401 - t348 * t387 + t350 * t388) * t390 + (t364 * t401 - t365 * t387 + t366 * t388) * t404) / 0.2e1 + t390 * ((t347 * t399 - t349 * t385 + t351 * t386) * t389 + (t399 * t346 - t385 * t348 + t386 * t350) * t390 + (t364 * t399 - t365 * t385 + t366 * t386) * t404) / 0.2e1 + t404 * ((-t347 * t444 - t349 * t397 + t351 * t398) * t389 + (-t346 * t444 - t348 * t397 + t350 * t398) * t390 + (-t364 * t444 - t397 * t365 + t398 * t366) * t404) / 0.2e1 + m(5) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(6) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + ((t457 * t362 + t456 * t363 + t458 * t387) * t382 + (t462 * t362 + t460 * t363 + t464 * t387) * t359 + (t461 * t362 + t459 * t363 + t463 * t387) * t358) * t358 / 0.2e1 + ((t457 * t360 + t456 * t361 + t458 * t385) * t382 + (t462 * t360 + t460 * t361 + t464 * t385) * t359 + (t461 * t360 + t459 * t361 + t463 * t385) * t358) * t359 / 0.2e1 + ((t457 * t383 + t456 * t384 + t458 * t397) * t382 + (t462 * t383 + t460 * t384 + t464 * t397) * t359 + (t461 * t383 + t459 * t384 + t463 * t397) * t358) * t382 / 0.2e1 + (m(2) * (t407 ^ 2 + t408 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
