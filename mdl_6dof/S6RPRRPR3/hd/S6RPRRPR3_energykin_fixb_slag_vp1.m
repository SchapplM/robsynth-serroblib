% Calculate kinetic energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:24
% EndTime: 2019-03-09 05:04:26
% DurationCPUTime: 2.13s
% Computational Cost: add. (1755->253), mult. (2195->390), div. (0->0), fcn. (2348->10), ass. (0->126)
t460 = Icges(5,1) + Icges(6,1);
t459 = -Icges(5,4) + Icges(6,5);
t458 = Icges(6,4) + Icges(5,5);
t457 = Icges(5,2) + Icges(6,3);
t456 = -Icges(6,6) + Icges(5,6);
t455 = -Icges(5,3) - Icges(6,2);
t406 = qJ(1) + pkin(10);
t404 = sin(t406);
t405 = cos(t406);
t412 = cos(qJ(4));
t408 = sin(qJ(4));
t413 = cos(qJ(3));
t436 = t408 * t413;
t367 = t404 * t436 + t405 * t412;
t435 = t412 * t413;
t368 = t404 * t435 - t405 * t408;
t409 = sin(qJ(3));
t438 = t404 * t409;
t454 = t457 * t367 + t459 * t368 - t456 * t438;
t369 = -t404 * t412 + t405 * t436;
t370 = t404 * t408 + t405 * t435;
t437 = t405 * t409;
t453 = t457 * t369 + t459 * t370 - t456 * t437;
t452 = -t456 * t367 + t458 * t368 - t455 * t438;
t451 = -t456 * t369 + t458 * t370 - t455 * t437;
t450 = t459 * t367 + t460 * t368 + t458 * t438;
t449 = t459 * t369 + t460 * t370 + t458 * t437;
t448 = t456 * t413 + (t457 * t408 + t459 * t412) * t409;
t447 = t455 * t413 + (-t456 * t408 + t458 * t412) * t409;
t446 = -t458 * t413 + (t459 * t408 + t460 * t412) * t409;
t410 = sin(qJ(1));
t441 = pkin(1) * t410;
t440 = Icges(4,4) * t409;
t439 = Icges(4,4) * t413;
t414 = cos(qJ(1));
t403 = qJD(1) * t414 * pkin(1);
t434 = qJD(1) * (pkin(2) * t405 + pkin(7) * t404) + t403;
t401 = qJD(3) * t404;
t432 = qJD(4) * t409;
t383 = t405 * t432 + t401;
t433 = qJD(3) * t405;
t431 = qJD(6) * t409;
t428 = pkin(3) * t413 + pkin(8) * t409;
t379 = t428 * t404;
t380 = t428 * t405;
t430 = t379 * t401 + t380 * t433 + qJD(2);
t429 = -pkin(2) * t404 + pkin(7) * t405 - t441;
t384 = t404 * t432 - t433;
t427 = rSges(4,1) * t413 - rSges(4,2) * t409;
t345 = pkin(4) * t368 + qJ(5) * t367;
t426 = qJD(5) * t409 * t408 + t383 * t345 + t430;
t425 = Icges(4,1) * t413 - t440;
t424 = -Icges(4,2) * t409 + t439;
t423 = Icges(4,5) * t413 - Icges(4,6) * t409;
t354 = -Icges(4,6) * t405 + t404 * t424;
t356 = -Icges(4,5) * t405 + t404 * t425;
t422 = t354 * t409 - t356 * t413;
t355 = Icges(4,6) * t404 + t405 * t424;
t357 = Icges(4,5) * t404 + t405 * t425;
t421 = -t355 * t409 + t357 * t413;
t391 = Icges(4,2) * t413 + t440;
t392 = Icges(4,1) * t409 + t439;
t420 = -t391 * t409 + t392 * t413;
t398 = pkin(3) * t409 - pkin(8) * t413;
t419 = qJD(1) * t380 - t398 * t401 + t434;
t346 = pkin(4) * t370 + qJ(5) * t369;
t402 = -qJD(4) * t413 + qJD(1);
t418 = qJD(5) * t367 + t402 * t346 + t419;
t417 = (-t379 + t429) * qJD(1) - t398 * t433;
t386 = (pkin(4) * t412 + qJ(5) * t408) * t409;
t416 = qJD(5) * t369 + t384 * t386 + t417;
t411 = cos(qJ(6));
t407 = sin(qJ(6));
t397 = rSges(2,1) * t414 - rSges(2,2) * t410;
t396 = rSges(2,1) * t410 + rSges(2,2) * t414;
t395 = rSges(4,1) * t409 + rSges(4,2) * t413;
t390 = Icges(4,5) * t409 + Icges(4,6) * t413;
t389 = qJD(1) + (-qJD(4) + qJD(6)) * t413;
t388 = pkin(5) * t409 * t412 + pkin(9) * t413;
t382 = (t407 * t408 + t411 * t412) * t409;
t381 = (-t407 * t412 + t408 * t411) * t409;
t378 = -rSges(5,3) * t413 + (rSges(5,1) * t412 - rSges(5,2) * t408) * t409;
t377 = -rSges(6,2) * t413 + (rSges(6,1) * t412 + rSges(6,3) * t408) * t409;
t363 = t403 + qJD(1) * (rSges(3,1) * t405 - rSges(3,2) * t404);
t362 = (-rSges(3,1) * t404 - rSges(3,2) * t405 - t441) * qJD(1);
t359 = rSges(4,3) * t404 + t405 * t427;
t358 = -rSges(4,3) * t405 + t404 * t427;
t353 = Icges(4,3) * t404 + t405 * t423;
t352 = -Icges(4,3) * t405 + t404 * t423;
t351 = -t404 * t431 + t384;
t350 = -t405 * t431 + t383;
t349 = pkin(5) * t370 - pkin(9) * t437;
t348 = pkin(5) * t368 - pkin(9) * t438;
t344 = rSges(7,1) * t382 + rSges(7,2) * t381 + rSges(7,3) * t413;
t343 = Icges(7,1) * t382 + Icges(7,4) * t381 + Icges(7,5) * t413;
t342 = Icges(7,4) * t382 + Icges(7,2) * t381 + Icges(7,6) * t413;
t341 = Icges(7,5) * t382 + Icges(7,6) * t381 + Icges(7,3) * t413;
t340 = t369 * t407 + t370 * t411;
t339 = t369 * t411 - t370 * t407;
t338 = t367 * t407 + t368 * t411;
t337 = t367 * t411 - t368 * t407;
t335 = rSges(5,1) * t370 - rSges(5,2) * t369 + rSges(5,3) * t437;
t334 = rSges(6,1) * t370 + rSges(6,2) * t437 + rSges(6,3) * t369;
t333 = rSges(5,1) * t368 - rSges(5,2) * t367 + rSges(5,3) * t438;
t332 = rSges(6,1) * t368 + rSges(6,2) * t438 + rSges(6,3) * t367;
t318 = qJD(1) * t359 - t395 * t401 + t434;
t317 = -t395 * t433 + (-t358 + t429) * qJD(1);
t316 = qJD(2) + (t358 * t404 + t359 * t405) * qJD(3);
t315 = rSges(7,1) * t340 + rSges(7,2) * t339 - rSges(7,3) * t437;
t314 = rSges(7,1) * t338 + rSges(7,2) * t337 - rSges(7,3) * t438;
t313 = Icges(7,1) * t340 + Icges(7,4) * t339 - Icges(7,5) * t437;
t312 = Icges(7,1) * t338 + Icges(7,4) * t337 - Icges(7,5) * t438;
t311 = Icges(7,4) * t340 + Icges(7,2) * t339 - Icges(7,6) * t437;
t310 = Icges(7,4) * t338 + Icges(7,2) * t337 - Icges(7,6) * t438;
t309 = Icges(7,5) * t340 + Icges(7,6) * t339 - Icges(7,3) * t437;
t308 = Icges(7,5) * t338 + Icges(7,6) * t337 - Icges(7,3) * t438;
t307 = t335 * t402 - t378 * t383 + t419;
t306 = -t333 * t402 + t378 * t384 + t417;
t305 = t333 * t383 - t335 * t384 + t430;
t304 = t334 * t402 + (-t377 - t386) * t383 + t418;
t303 = t377 * t384 + (-t332 - t345) * t402 + t416;
t302 = t332 * t383 + (-t334 - t346) * t384 + t426;
t301 = t315 * t389 - t344 * t350 + t349 * t402 + (-t386 - t388) * t383 + t418;
t300 = -t314 * t389 + t344 * t351 + t384 * t388 + (-t345 - t348) * t402 + t416;
t299 = t314 * t350 - t315 * t351 + t348 * t383 + (-t346 - t349) * t384 + t426;
t1 = ((t404 * t390 + t405 * t420) * qJD(1) + (t404 ^ 2 * t353 + (t422 * t405 + (-t352 + t421) * t404) * t405) * qJD(3)) * t401 / 0.2e1 + m(7) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + m(6) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 - ((-t405 * t390 + t404 * t420) * qJD(1) + (t405 ^ 2 * t352 + (t421 * t404 + (-t353 + t422) * t405) * t404) * qJD(3)) * t433 / 0.2e1 + qJD(1) * ((t413 * t391 + t409 * t392) * qJD(1) + ((t355 * t413 + t357 * t409) * t404 - (t354 * t413 + t356 * t409) * t405) * qJD(3)) / 0.2e1 + t350 * ((-t309 * t437 + t339 * t311 + t340 * t313) * t350 + (-t308 * t437 + t310 * t339 + t312 * t340) * t351 + (t339 * t342 + t340 * t343 - t341 * t437) * t389) / 0.2e1 + t351 * ((-t309 * t438 + t311 * t337 + t313 * t338) * t350 + (-t308 * t438 + t337 * t310 + t338 * t312) * t351 + (t337 * t342 + t338 * t343 - t341 * t438) * t389) / 0.2e1 + t389 * ((t309 * t413 + t311 * t381 + t313 * t382) * t350 + (t308 * t413 + t310 * t381 + t312 * t382) * t351 + (t413 * t341 + t381 * t342 + t382 * t343) * t389) / 0.2e1 + m(5) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + m(4) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t362 ^ 2 + t363 ^ 2) / 0.2e1 + ((t448 * t369 + t446 * t370 + t447 * t437) * t402 + (t454 * t369 + t450 * t370 + t452 * t437) * t384 + (t453 * t369 + t449 * t370 + t451 * t437) * t383) * t383 / 0.2e1 + ((t448 * t367 + t446 * t368 + t447 * t438) * t402 + (t454 * t367 + t450 * t368 + t452 * t438) * t384 + (t453 * t367 + t449 * t368 + t451 * t438) * t383) * t384 / 0.2e1 + ((-t451 * t383 - t452 * t384 - t447 * t402) * t413 + ((t448 * t408 + t446 * t412) * t402 + (t454 * t408 + t450 * t412) * t384 + (t453 * t408 + t449 * t412) * t383) * t409) * t402 / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t396 ^ 2 + t397 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
