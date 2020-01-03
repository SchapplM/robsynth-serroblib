% Calculate kinetic energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:23
% EndTime: 2019-12-31 20:09:25
% DurationCPUTime: 2.23s
% Computational Cost: add. (636->188), mult. (1489->287), div. (0->0), fcn. (1442->6), ass. (0->110)
t474 = Icges(3,4) + Icges(4,6);
t473 = Icges(3,1) + Icges(4,2);
t472 = -Icges(3,2) - Icges(4,3);
t380 = cos(qJ(2));
t471 = t474 * t380;
t377 = sin(qJ(2));
t470 = t474 * t377;
t469 = -Icges(4,4) + Icges(3,5);
t468 = Icges(4,5) - Icges(3,6);
t467 = t472 * t377 + t471;
t466 = -t473 * t380 + t470;
t465 = Icges(4,1) + Icges(3,3);
t464 = Icges(5,1) + Icges(6,1);
t463 = Icges(5,4) + Icges(6,4);
t462 = Icges(6,5) + Icges(5,5);
t461 = Icges(5,2) + Icges(6,2);
t460 = Icges(6,6) + Icges(5,6);
t459 = Icges(6,3) + Icges(5,3);
t378 = sin(qJ(1));
t381 = cos(qJ(1));
t458 = t467 * t378 + t468 * t381;
t457 = -t468 * t378 + t467 * t381;
t456 = t466 * t378 + t469 * t381;
t455 = t469 * t378 - t466 * t381;
t454 = t472 * t380 - t470;
t453 = t473 * t377 + t471;
t452 = t468 * t377 + t469 * t380;
t379 = cos(qJ(4));
t417 = t381 * t379;
t376 = sin(qJ(4));
t422 = t378 * t376;
t347 = t377 * t417 - t422;
t418 = t381 * t376;
t421 = t378 * t379;
t348 = t377 * t418 + t421;
t419 = t380 * t381;
t451 = t460 * t347 + t462 * t348 + t459 * t419;
t349 = t377 * t421 + t418;
t350 = t377 * t422 - t417;
t420 = t378 * t380;
t450 = t460 * t349 + t462 * t350 + t459 * t420;
t449 = t461 * t347 + t463 * t348 + t460 * t419;
t448 = t461 * t349 + t463 * t350 + t460 * t420;
t447 = t463 * t347 + t464 * t348 + t462 * t419;
t446 = t463 * t349 + t464 * t350 + t462 * t420;
t445 = (-t462 * t376 - t460 * t379) * t380 + t459 * t377;
t444 = t452 * t378 - t465 * t381;
t443 = t465 * t378 + t452 * t381;
t442 = (-t463 * t376 - t461 * t379) * t380 + t460 * t377;
t441 = (-t464 * t376 - t463 * t379) * t380 + t462 * t377;
t440 = t469 * t377 - t468 * t380;
t439 = t454 * t377 + t453 * t380;
t438 = -t457 * t377 + t455 * t380;
t437 = t458 * t377 + t456 * t380;
t430 = pkin(4) * t376;
t428 = t379 * pkin(4);
t384 = qJ(5) * t380 + t377 * t430;
t416 = rSges(6,1) * t348 + rSges(6,2) * t347 + rSges(6,3) * t419 + t378 * t428 + t381 * t384;
t415 = rSges(6,1) * t350 + rSges(6,2) * t349 + rSges(6,3) * t420 + t378 * t384 - t381 * t428;
t414 = (-rSges(6,1) * t376 - rSges(6,2) * t379 - t430) * t380 + (rSges(6,3) + qJ(5)) * t377;
t400 = pkin(2) * t380 + qJ(3) * t377;
t351 = t400 * t378;
t370 = t378 * pkin(1) - t381 * pkin(6);
t413 = -t351 - t370;
t412 = qJD(2) * t378;
t411 = qJD(2) * t381;
t410 = qJD(3) * t377;
t409 = qJD(4) * t380;
t352 = t400 * t381;
t358 = qJD(1) * (t381 * pkin(1) + t378 * pkin(6));
t408 = qJD(1) * t352 + t378 * t410 + t358;
t365 = t377 * pkin(2) - t380 * qJ(3);
t405 = qJD(2) * (t377 * rSges(4,2) + t380 * rSges(4,3) - t365);
t356 = t378 * pkin(3) + pkin(7) * t419;
t404 = qJD(1) * t356 + t408;
t403 = -qJD(3) * t380 + t351 * t412 + t352 * t411;
t402 = rSges(3,1) * t380 - rSges(3,2) * t377;
t401 = -rSges(4,2) * t380 + rSges(4,3) * t377;
t399 = (-pkin(7) * t377 - t365) * qJD(2);
t357 = -t381 * pkin(3) + pkin(7) * t420;
t372 = t381 * t410;
t386 = t372 + (-t357 + t413) * qJD(1);
t385 = t356 * t411 + t357 * t412 + t403;
t383 = qJD(5) * t380 + t399;
t373 = qJD(4) * t377 + qJD(1);
t369 = t381 * rSges(2,1) - t378 * rSges(2,2);
t368 = t378 * rSges(2,1) + t381 * rSges(2,2);
t367 = t377 * rSges(3,1) + t380 * rSges(3,2);
t355 = t378 * t409 - t411;
t354 = t381 * t409 + t412;
t342 = -t381 * rSges(4,1) + t378 * t401;
t341 = t378 * rSges(4,1) + t381 * t401;
t340 = t378 * rSges(3,3) + t381 * t402;
t339 = t377 * rSges(5,3) + (-rSges(5,1) * t376 - rSges(5,2) * t379) * t380;
t337 = -t381 * rSges(3,3) + t378 * t402;
t314 = rSges(5,1) * t350 + rSges(5,2) * t349 + rSges(5,3) * t420;
t312 = rSges(5,1) * t348 + rSges(5,2) * t347 + rSges(5,3) * t419;
t298 = qJD(1) * t340 - t367 * t412 + t358;
t297 = -t367 * t411 + (-t337 - t370) * qJD(1);
t296 = (t337 * t378 + t340 * t381) * qJD(2);
t295 = qJD(1) * t341 + t378 * t405 + t408;
t294 = t372 + t381 * t405 + (-t342 + t413) * qJD(1);
t293 = (t341 * t381 + t342 * t378) * qJD(2) + t403;
t292 = t373 * t312 - t354 * t339 + t378 * t399 + t404;
t291 = -t373 * t314 + t355 * t339 + t381 * t399 + t386;
t290 = -t312 * t355 + t314 * t354 + t385;
t289 = -t354 * t414 + t373 * t416 + t378 * t383 + t404;
t288 = t355 * t414 - t373 * t415 + t381 * t383 + t386;
t287 = qJD(5) * t377 + t354 * t415 - t355 * t416 + t385;
t1 = m(3) * (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + m(4) * (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) / 0.2e1 + m(5) * (t290 ^ 2 + t291 ^ 2 + t292 ^ 2) / 0.2e1 + m(6) * (t287 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + ((t442 * t347 + t441 * t348 + t445 * t419) * t373 + (t448 * t347 + t446 * t348 + t450 * t419) * t355 + (t449 * t347 + t447 * t348 + t451 * t419) * t354) * t354 / 0.2e1 + ((t442 * t349 + t441 * t350 + t445 * t420) * t373 + (t448 * t349 + t446 * t350 + t450 * t420) * t355 + (t449 * t349 + t447 * t350 + t451 * t420) * t354) * t355 / 0.2e1 + (((-t441 * t376 - t442 * t379) * t373 + (-t446 * t376 - t448 * t379) * t355 + (-t447 * t376 - t449 * t379) * t354) * t380 + (t451 * t354 + t450 * t355 + t445 * t373) * t377) * t373 / 0.2e1 + (m(2) * (t368 ^ 2 + t369 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t456 * t377 - t458 * t380) * t381 + (t455 * t377 + t457 * t380) * t378) * qJD(2) + (t453 * t377 - t454 * t380) * qJD(1)) * qJD(1) / 0.2e1 + ((t443 * t378 ^ 2 + (t437 * t381 + (t438 - t444) * t378) * t381) * qJD(2) + (t440 * t378 + t439 * t381) * qJD(1)) * t412 / 0.2e1 - ((t444 * t381 ^ 2 + (t438 * t378 + (t437 - t443) * t381) * t378) * qJD(2) + (t439 * t378 - t440 * t381) * qJD(1)) * t411 / 0.2e1;
T = t1;
