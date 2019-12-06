% Calculate kinetic energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:35
% EndTime: 2019-12-05 15:18:36
% DurationCPUTime: 1.16s
% Computational Cost: add. (3043->239), mult. (8424->386), div. (0->0), fcn. (10950->14), ass. (0->116)
t479 = cos(qJ(3));
t478 = cos(qJ(4));
t477 = cos(pkin(5));
t476 = cos(pkin(6));
t475 = cos(pkin(11));
t474 = sin(pkin(6));
t473 = sin(pkin(11));
t448 = sin(pkin(10));
t450 = cos(pkin(10));
t465 = t477 * t473;
t439 = t448 * t475 + t450 * t465;
t453 = sin(qJ(3));
t467 = t477 * t475;
t459 = t448 * t473 - t450 * t467;
t457 = t459 * t476;
t449 = sin(pkin(5));
t468 = t449 * t474;
t463 = t479 * t468;
t418 = t439 * t453 + t450 * t463 + t479 * t457;
t469 = t449 * t476;
t432 = -t450 * t469 + t459 * t474;
t428 = qJD(3) * t432;
t409 = qJD(4) * t418 + t428;
t440 = -t448 * t465 + t450 * t475;
t458 = t448 * t467 + t450 * t473;
t456 = t458 * t476;
t420 = t440 * t453 - t448 * t463 + t479 * t456;
t433 = t448 * t469 + t458 * t474;
t429 = qJD(3) * t433;
t410 = qJD(4) * t420 + t429;
t464 = t476 * t475;
t466 = t477 * t474;
t430 = -t479 * t466 + (t453 * t473 - t464 * t479) * t449;
t438 = -t475 * t468 + t477 * t476;
t437 = qJD(3) * t438;
t424 = qJD(4) * t430 + t437;
t472 = qJD(2) * t449;
t445 = qJD(2) * t477 + qJD(1);
t470 = t450 * t472;
t419 = t439 * t479 + (-t450 * t468 - t457) * t453;
t398 = pkin(3) * t419 + pkin(8) * t418;
t431 = t453 * t466 + (t453 * t464 + t473 * t479) * t449;
t415 = pkin(3) * t431 + pkin(8) * t430;
t444 = t448 * t472;
t462 = -t398 * t437 + t415 * t428 + t444;
t421 = t440 * t479 + (t448 * t468 - t456) * t453;
t399 = pkin(3) * t421 + pkin(8) * t420;
t461 = t398 * t429 - t399 * t428 + t445;
t460 = t399 * t437 - t415 * t429 - t470;
t454 = cos(qJ(5));
t452 = sin(qJ(4));
t451 = sin(qJ(5));
t423 = t431 * t478 + t438 * t452;
t422 = t431 * t452 - t438 * t478;
t414 = rSges(4,1) * t431 - rSges(4,2) * t430 + rSges(4,3) * t438;
t413 = Icges(4,1) * t431 - Icges(4,4) * t430 + Icges(4,5) * t438;
t412 = Icges(4,4) * t431 - Icges(4,2) * t430 + Icges(4,6) * t438;
t411 = Icges(4,5) * t431 - Icges(4,6) * t430 + Icges(4,3) * t438;
t408 = t423 * t454 + t430 * t451;
t407 = -t423 * t451 + t430 * t454;
t406 = t421 * t478 + t433 * t452;
t405 = t421 * t452 - t433 * t478;
t404 = t419 * t478 + t432 * t452;
t403 = t419 * t452 - t432 * t478;
t401 = qJD(5) * t422 + t424;
t400 = pkin(4) * t423 + pkin(9) * t422;
t396 = rSges(5,1) * t423 - rSges(5,2) * t422 + rSges(5,3) * t430;
t395 = Icges(5,1) * t423 - Icges(5,4) * t422 + Icges(5,5) * t430;
t394 = Icges(5,4) * t423 - Icges(5,2) * t422 + Icges(5,6) * t430;
t393 = Icges(5,5) * t423 - Icges(5,6) * t422 + Icges(5,3) * t430;
t391 = rSges(4,1) * t421 - rSges(4,2) * t420 + rSges(4,3) * t433;
t390 = rSges(4,1) * t419 - rSges(4,2) * t418 + rSges(4,3) * t432;
t389 = Icges(4,1) * t421 - Icges(4,4) * t420 + Icges(4,5) * t433;
t388 = Icges(4,1) * t419 - Icges(4,4) * t418 + Icges(4,5) * t432;
t387 = Icges(4,4) * t421 - Icges(4,2) * t420 + Icges(4,6) * t433;
t386 = Icges(4,4) * t419 - Icges(4,2) * t418 + Icges(4,6) * t432;
t385 = Icges(4,5) * t421 - Icges(4,6) * t420 + Icges(4,3) * t433;
t384 = Icges(4,5) * t419 - Icges(4,6) * t418 + Icges(4,3) * t432;
t383 = t406 * t454 + t420 * t451;
t382 = -t406 * t451 + t420 * t454;
t381 = t404 * t454 + t418 * t451;
t380 = -t404 * t451 + t418 * t454;
t379 = qJD(5) * t405 + t410;
t378 = qJD(5) * t403 + t409;
t377 = pkin(4) * t406 + pkin(9) * t405;
t376 = pkin(4) * t404 + pkin(9) * t403;
t375 = rSges(6,1) * t408 + rSges(6,2) * t407 + rSges(6,3) * t422;
t374 = Icges(6,1) * t408 + Icges(6,4) * t407 + Icges(6,5) * t422;
t373 = Icges(6,4) * t408 + Icges(6,2) * t407 + Icges(6,6) * t422;
t372 = Icges(6,5) * t408 + Icges(6,6) * t407 + Icges(6,3) * t422;
t371 = rSges(5,1) * t406 - rSges(5,2) * t405 + rSges(5,3) * t420;
t370 = rSges(5,1) * t404 - rSges(5,2) * t403 + rSges(5,3) * t418;
t369 = Icges(5,1) * t406 - Icges(5,4) * t405 + Icges(5,5) * t420;
t368 = Icges(5,1) * t404 - Icges(5,4) * t403 + Icges(5,5) * t418;
t367 = Icges(5,4) * t406 - Icges(5,2) * t405 + Icges(5,6) * t420;
t366 = Icges(5,4) * t404 - Icges(5,2) * t403 + Icges(5,6) * t418;
t365 = Icges(5,5) * t406 - Icges(5,6) * t405 + Icges(5,3) * t420;
t364 = Icges(5,5) * t404 - Icges(5,6) * t403 + Icges(5,3) * t418;
t363 = -t470 + (t391 * t438 - t414 * t433) * qJD(3);
t362 = t444 + (-t390 * t438 + t414 * t432) * qJD(3);
t361 = (t390 * t433 - t391 * t432) * qJD(3) + t445;
t360 = rSges(6,1) * t383 + rSges(6,2) * t382 + rSges(6,3) * t405;
t359 = rSges(6,1) * t381 + rSges(6,2) * t380 + rSges(6,3) * t403;
t358 = Icges(6,1) * t383 + Icges(6,4) * t382 + Icges(6,5) * t405;
t357 = Icges(6,1) * t381 + Icges(6,4) * t380 + Icges(6,5) * t403;
t356 = Icges(6,4) * t383 + Icges(6,2) * t382 + Icges(6,6) * t405;
t355 = Icges(6,4) * t381 + Icges(6,2) * t380 + Icges(6,6) * t403;
t354 = Icges(6,5) * t383 + Icges(6,6) * t382 + Icges(6,3) * t405;
t353 = Icges(6,5) * t381 + Icges(6,6) * t380 + Icges(6,3) * t403;
t352 = t371 * t424 - t396 * t410 + t460;
t351 = -t370 * t424 + t396 * t409 + t462;
t350 = t370 * t410 - t371 * t409 + t461;
t349 = t360 * t401 - t375 * t379 + t377 * t424 - t400 * t410 + t460;
t348 = -t359 * t401 + t375 * t378 - t376 * t424 + t400 * t409 + t462;
t347 = t359 * t379 - t360 * t378 + t376 * t410 - t377 * t409 + t461;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t445 ^ 2 + (t448 ^ 2 + t450 ^ 2) * qJD(2) ^ 2 * t449 ^ 2) / 0.2e1 + m(4) * (t361 ^ 2 + t362 ^ 2 + t363 ^ 2) / 0.2e1 + m(5) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + t410 * ((t420 * t365 - t405 * t367 + t406 * t369) * t410 + (t364 * t420 - t366 * t405 + t368 * t406) * t409 + (t393 * t420 - t394 * t405 + t395 * t406) * t424) / 0.2e1 + t409 * ((t365 * t418 - t367 * t403 + t369 * t404) * t410 + (t418 * t364 - t403 * t366 + t404 * t368) * t409 + (t393 * t418 - t394 * t403 + t395 * t404) * t424) / 0.2e1 + t424 * ((t365 * t430 - t367 * t422 + t369 * t423) * t410 + (t364 * t430 - t366 * t422 + t368 * t423) * t409 + (t430 * t393 - t422 * t394 + t423 * t395) * t424) / 0.2e1 + m(6) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + t379 * ((t405 * t354 + t382 * t356 + t383 * t358) * t379 + (t353 * t405 + t355 * t382 + t357 * t383) * t378 + (t372 * t405 + t373 * t382 + t374 * t383) * t401) / 0.2e1 + t378 * ((t354 * t403 + t356 * t380 + t358 * t381) * t379 + (t403 * t353 + t380 * t355 + t381 * t357) * t378 + (t372 * t403 + t373 * t380 + t374 * t381) * t401) / 0.2e1 + t401 * ((t354 * t422 + t356 * t407 + t358 * t408) * t379 + (t353 * t422 + t355 * t407 + t357 * t408) * t378 + (t422 * t372 + t407 * t373 + t408 * t374) * t401) / 0.2e1 + (t433 * ((t385 * t433 - t387 * t420 + t389 * t421) * t433 + (t384 * t433 - t386 * t420 + t388 * t421) * t432 + (t411 * t433 - t412 * t420 + t413 * t421) * t438) + t432 * ((t385 * t432 - t387 * t418 + t389 * t419) * t433 + (t384 * t432 - t386 * t418 + t388 * t419) * t432 + (t411 * t432 - t412 * t418 + t413 * t419) * t438) + t438 * ((t385 * t438 - t387 * t430 + t389 * t431) * t433 + (t384 * t438 - t386 * t430 + t388 * t431) * t432 + (t411 * t438 - t412 * t430 + t413 * t431) * t438)) * qJD(3) ^ 2 / 0.2e1;
T = t1;
