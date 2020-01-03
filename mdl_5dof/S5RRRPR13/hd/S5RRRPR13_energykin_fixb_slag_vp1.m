% Calculate kinetic energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR13_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR13_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:16
% EndTime: 2019-12-31 21:43:18
% DurationCPUTime: 2.04s
% Computational Cost: add. (1528->260), mult. (3829->398), div. (0->0), fcn. (4578->10), ass. (0->121)
t462 = Icges(4,1) + Icges(5,2);
t461 = Icges(5,1) + Icges(4,3);
t460 = -Icges(4,4) - Icges(5,6);
t459 = Icges(5,4) - Icges(4,5);
t458 = Icges(5,5) - Icges(4,6);
t457 = Icges(4,2) + Icges(5,3);
t417 = sin(qJ(2));
t418 = sin(qJ(1));
t420 = cos(qJ(2));
t421 = cos(qJ(1));
t440 = cos(pkin(5));
t429 = t421 * t440;
t398 = t417 * t418 - t420 * t429;
t399 = t417 * t429 + t418 * t420;
t415 = sin(pkin(5));
t437 = t415 * t421;
t363 = Icges(3,5) * t399 - Icges(3,6) * t398 - Icges(3,3) * t437;
t430 = t418 * t440;
t400 = t421 * t417 + t420 * t430;
t401 = -t417 * t430 + t421 * t420;
t439 = t415 * t418;
t364 = Icges(3,5) * t401 - Icges(3,6) * t400 + Icges(3,3) * t439;
t456 = (t363 * t421 - t364 * t418) * t415;
t442 = cos(qJ(3));
t433 = t415 * t442;
t441 = sin(qJ(3));
t382 = t399 * t441 + t421 * t433;
t432 = t415 * t441;
t383 = t399 * t442 - t421 * t432;
t455 = t457 * t382 + t460 * t383 + t458 * t398;
t384 = t401 * t441 - t418 * t433;
t385 = t401 * t442 + t418 * t432;
t454 = t457 * t384 + t460 * t385 + t458 * t400;
t453 = t458 * t382 - t459 * t383 + t461 * t398;
t452 = t458 * t384 - t459 * t385 + t461 * t400;
t451 = t460 * t382 + t462 * t383 - t459 * t398;
t450 = t460 * t384 + t462 * t385 - t459 * t400;
t396 = t417 * t432 - t440 * t442;
t397 = t417 * t433 + t440 * t441;
t438 = t415 * t420;
t449 = t457 * t396 + t460 * t397 - t458 * t438;
t448 = t460 * t396 + t462 * t397 + t459 * t438;
t447 = t458 * t396 - t459 * t397 - t461 * t438;
t375 = pkin(2) * t399 + pkin(8) * t398;
t376 = pkin(2) * t401 + pkin(8) * t400;
t434 = qJD(2) * t415;
t411 = t418 * t434;
t431 = t421 * t434;
t436 = t375 * t411 + t376 * t431;
t386 = qJD(3) * t400 + t411;
t435 = qJD(1) * (pkin(1) * t418 - pkin(7) * t437);
t412 = qJD(2) * t440 + qJD(1);
t345 = pkin(3) * t383 + qJ(4) * t382;
t428 = qJD(4) * t396 + t386 * t345 + t436;
t387 = qJD(3) * t398 - t431;
t403 = -qJD(3) * t438 + t412;
t402 = (pkin(2) * t417 - pkin(8) * t420) * t415;
t404 = qJD(1) * (pkin(1) * t421 + pkin(7) * t439);
t426 = t412 * t376 - t402 * t411 + t404;
t346 = pkin(3) * t385 + qJ(4) * t384;
t425 = qJD(4) * t382 + t403 * t346 + t426;
t424 = -t375 * t412 - t402 * t431 - t435;
t374 = pkin(3) * t397 + qJ(4) * t396;
t423 = qJD(4) * t384 + t387 * t374 + t424;
t419 = cos(qJ(5));
t416 = sin(qJ(5));
t407 = rSges(2,1) * t421 - rSges(2,2) * t418;
t406 = rSges(2,1) * t418 + rSges(2,2) * t421;
t392 = t440 * rSges(3,3) + (rSges(3,1) * t417 + rSges(3,2) * t420) * t415;
t391 = Icges(3,5) * t440 + (Icges(3,1) * t417 + Icges(3,4) * t420) * t415;
t390 = Icges(3,6) * t440 + (Icges(3,4) * t417 + Icges(3,2) * t420) * t415;
t389 = Icges(3,3) * t440 + (Icges(3,5) * t417 + Icges(3,6) * t420) * t415;
t388 = -pkin(4) * t438 + pkin(9) * t397;
t381 = t396 * t416 - t419 * t438;
t380 = t396 * t419 + t416 * t438;
t377 = qJD(5) * t397 + t403;
t371 = rSges(3,1) * t401 - rSges(3,2) * t400 + rSges(3,3) * t439;
t370 = rSges(3,1) * t399 - rSges(3,2) * t398 - rSges(3,3) * t437;
t368 = Icges(3,1) * t401 - Icges(3,4) * t400 + Icges(3,5) * t439;
t367 = Icges(3,1) * t399 - Icges(3,4) * t398 - Icges(3,5) * t437;
t366 = Icges(3,4) * t401 - Icges(3,2) * t400 + Icges(3,6) * t439;
t365 = Icges(3,4) * t399 - Icges(3,2) * t398 - Icges(3,6) * t437;
t362 = rSges(4,1) * t397 - rSges(4,2) * t396 - rSges(4,3) * t438;
t361 = -rSges(5,1) * t438 - rSges(5,2) * t397 + rSges(5,3) * t396;
t354 = pkin(4) * t400 + pkin(9) * t385;
t353 = pkin(4) * t398 + pkin(9) * t383;
t352 = t384 * t416 + t400 * t419;
t351 = t384 * t419 - t400 * t416;
t350 = t382 * t416 + t398 * t419;
t349 = t382 * t419 - t398 * t416;
t348 = qJD(5) * t383 + t387;
t347 = qJD(5) * t385 + t386;
t342 = rSges(4,1) * t385 - rSges(4,2) * t384 + rSges(4,3) * t400;
t341 = rSges(4,1) * t383 - rSges(4,2) * t382 + rSges(4,3) * t398;
t340 = rSges(5,1) * t400 - rSges(5,2) * t385 + rSges(5,3) * t384;
t339 = rSges(5,1) * t398 - rSges(5,2) * t383 + rSges(5,3) * t382;
t326 = rSges(6,1) * t381 + rSges(6,2) * t380 + rSges(6,3) * t397;
t325 = Icges(6,1) * t381 + Icges(6,4) * t380 + Icges(6,5) * t397;
t324 = Icges(6,4) * t381 + Icges(6,2) * t380 + Icges(6,6) * t397;
t323 = Icges(6,5) * t381 + Icges(6,6) * t380 + Icges(6,3) * t397;
t321 = t371 * t412 - t392 * t411 + t404;
t320 = -t370 * t412 - t392 * t431 - t435;
t319 = (t370 * t418 + t371 * t421) * t434;
t318 = rSges(6,1) * t352 + rSges(6,2) * t351 + rSges(6,3) * t385;
t317 = rSges(6,1) * t350 + rSges(6,2) * t349 + rSges(6,3) * t383;
t316 = Icges(6,1) * t352 + Icges(6,4) * t351 + Icges(6,5) * t385;
t315 = Icges(6,1) * t350 + Icges(6,4) * t349 + Icges(6,5) * t383;
t314 = Icges(6,4) * t352 + Icges(6,2) * t351 + Icges(6,6) * t385;
t313 = Icges(6,4) * t350 + Icges(6,2) * t349 + Icges(6,6) * t383;
t312 = Icges(6,5) * t352 + Icges(6,6) * t351 + Icges(6,3) * t385;
t311 = Icges(6,5) * t350 + Icges(6,6) * t349 + Icges(6,3) * t383;
t310 = t342 * t403 - t362 * t386 + t426;
t309 = -t341 * t403 + t362 * t387 + t424;
t308 = t341 * t386 - t342 * t387 + t436;
t307 = t340 * t403 + (-t361 - t374) * t386 + t425;
t306 = t361 * t387 + (-t339 - t345) * t403 + t423;
t305 = t339 * t386 + (-t340 - t346) * t387 + t428;
t304 = t318 * t377 - t326 * t347 + t354 * t403 + (-t374 - t388) * t386 + t425;
t303 = -t317 * t377 + t326 * t348 + t387 * t388 + (-t345 - t353) * t403 + t423;
t302 = t317 * t347 - t318 * t348 + t353 * t386 + (-t346 - t354) * t387 + t428;
t1 = m(3) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + ((t389 * t439 - t390 * t400 + t391 * t401) * t412 + (-(-t365 * t400 + t401 * t367) * t421 + (-t400 * t366 + t401 * t368 - t456) * t418) * t434) * t411 / 0.2e1 - ((-t389 * t437 - t390 * t398 + t391 * t399) * t412 + ((-t366 * t398 + t368 * t399) * t418 + (t398 * t365 - t399 * t367 + t456) * t421) * t434) * t431 / 0.2e1 + t412 * ((t440 * t364 + (t366 * t420 + t368 * t417) * t415) * t411 - (t440 * t363 + (t365 * t420 + t367 * t417) * t415) * t431 + (t440 * t389 + (t390 * t420 + t391 * t417) * t415) * t412) / 0.2e1 + m(4) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(5) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + m(6) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + t347 * ((t385 * t312 + t351 * t314 + t352 * t316) * t347 + (t311 * t385 + t313 * t351 + t315 * t352) * t348 + (t323 * t385 + t324 * t351 + t325 * t352) * t377) / 0.2e1 + t348 * ((t312 * t383 + t314 * t349 + t316 * t350) * t347 + (t383 * t311 + t349 * t313 + t350 * t315) * t348 + (t323 * t383 + t324 * t349 + t325 * t350) * t377) / 0.2e1 + t377 * ((t312 * t397 + t314 * t380 + t316 * t381) * t347 + (t311 * t397 + t313 * t380 + t315 * t381) * t348 + (t397 * t323 + t380 * t324 + t381 * t325) * t377) / 0.2e1 + ((t449 * t384 + t448 * t385 + t447 * t400) * t403 + (t455 * t384 + t451 * t385 + t453 * t400) * t387 + (t454 * t384 + t450 * t385 + t452 * t400) * t386) * t386 / 0.2e1 + ((t449 * t382 + t448 * t383 + t447 * t398) * t403 + (t455 * t382 + t451 * t383 + t453 * t398) * t387 + (t454 * t382 + t450 * t383 + t452 * t398) * t386) * t387 / 0.2e1 + ((t449 * t396 + t448 * t397 - t447 * t438) * t403 + (t455 * t396 + t451 * t397 - t453 * t438) * t387 + (t454 * t396 + t450 * t397 - t452 * t438) * t386) * t403 / 0.2e1 + (m(2) * (t406 ^ 2 + t407 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
