% Calculate kinetic energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:06
% EndTime: 2019-03-09 02:33:07
% DurationCPUTime: 1.10s
% Computational Cost: add. (1088->214), mult. (1091->343), div. (0->0), fcn. (998->10), ass. (0->118)
t403 = sin(qJ(1));
t405 = cos(qJ(1));
t451 = qJD(1) * t405 * qJ(3) + qJD(3) * t403;
t398 = pkin(10) + qJ(4);
t389 = qJ(5) + t398;
t385 = cos(t389);
t384 = sin(t389);
t444 = Icges(6,4) * t384;
t418 = Icges(6,2) * t385 + t444;
t342 = Icges(6,6) * t405 + t418 * t403;
t343 = Icges(6,6) * t403 - t418 * t405;
t443 = Icges(6,4) * t385;
t420 = Icges(6,1) * t384 + t443;
t344 = Icges(6,5) * t405 + t420 * t403;
t345 = Icges(6,5) * t403 - t420 * t405;
t367 = -Icges(6,2) * t384 + t443;
t368 = Icges(6,1) * t385 - t444;
t392 = qJD(4) * t403;
t379 = qJD(5) * t403 + t392;
t393 = qJD(4) * t405;
t380 = qJD(5) * t405 + t393;
t450 = (t342 * t385 + t344 * t384) * t380 + (t343 * t385 + t345 * t384) * t379 + (t367 * t385 + t368 * t384) * qJD(1);
t399 = sin(pkin(10));
t448 = pkin(3) * t399;
t387 = sin(t398);
t446 = Icges(5,4) * t387;
t388 = cos(t398);
t445 = Icges(5,4) * t388;
t442 = t385 * t403;
t441 = t385 * t405;
t402 = sin(qJ(6));
t440 = t402 * t403;
t439 = t402 * t405;
t404 = cos(qJ(6));
t438 = t403 * t404;
t437 = t404 * t405;
t396 = qJD(2) * t403;
t435 = qJD(3) * t405 + t396;
t433 = qJD(6) * t385;
t432 = pkin(4) * qJD(4) * t388;
t431 = t403 * t432 + t435;
t430 = pkin(4) * t387;
t376 = qJD(1) * (pkin(1) * t405 + qJ(2) * t403);
t429 = -qJD(2) * t405 + t376;
t428 = qJD(1) * (pkin(7) * t405 + t403 * t448) + t376 + t451;
t427 = pkin(5) * t384 - pkin(9) * t385;
t381 = pkin(1) * t403 - qJ(2) * t405;
t426 = t405 * t448 - t381 + (-pkin(7) - qJ(3)) * t403;
t332 = pkin(8) * t403 - t430 * t405;
t333 = pkin(8) * t405 + t430 * t403;
t425 = t332 * t393 - t333 * t392;
t400 = cos(pkin(10));
t424 = rSges(4,1) * t399 + rSges(4,2) * t400;
t423 = rSges(5,1) * t387 + rSges(5,2) * t388;
t422 = rSges(6,1) * t384 + rSges(6,2) * t385;
t421 = Icges(5,1) * t387 + t445;
t419 = Icges(5,2) * t388 + t446;
t417 = Icges(5,5) * t387 + Icges(5,6) * t388;
t416 = Icges(6,5) * t384 + Icges(6,6) * t385;
t350 = Icges(5,6) * t405 + t419 * t403;
t352 = Icges(5,5) * t405 + t421 * t403;
t413 = -t350 * t388 - t352 * t387;
t351 = Icges(5,6) * t403 - t419 * t405;
t353 = Icges(5,5) * t403 - t421 * t405;
t412 = t351 * t388 + t353 * t387;
t372 = -Icges(5,2) * t387 + t445;
t373 = Icges(5,1) * t388 - t446;
t410 = t372 * t388 + t373 * t387;
t409 = -t332 + t426;
t408 = qJD(1) * (Icges(6,5) * t385 - Icges(6,6) * t384) + (Icges(6,3) * t405 + t416 * t403) * t380 + (Icges(6,3) * t403 - t416 * t405) * t379;
t407 = qJD(1) * t333 + (-qJD(2) - t432) * t405 + t428;
t383 = rSges(2,1) * t405 - rSges(2,2) * t403;
t382 = rSges(2,1) * t403 + rSges(2,2) * t405;
t378 = qJD(6) * t384 + qJD(1);
t374 = rSges(5,1) * t388 - rSges(5,2) * t387;
t371 = Icges(5,5) * t388 - Icges(5,6) * t387;
t370 = pkin(5) * t385 + pkin(9) * t384;
t369 = rSges(6,1) * t385 - rSges(6,2) * t384;
t364 = -t384 * t437 + t440;
t363 = t384 * t439 + t438;
t362 = t384 * t438 + t439;
t361 = -t384 * t440 + t437;
t359 = -t403 * t433 + t380;
t358 = t405 * t433 + t379;
t357 = t427 * t405;
t356 = t427 * t403;
t355 = rSges(5,3) * t403 - t423 * t405;
t354 = rSges(5,3) * t405 + t423 * t403;
t349 = Icges(5,3) * t403 - t417 * t405;
t348 = Icges(5,3) * t405 + t417 * t403;
t347 = rSges(6,3) * t403 - t422 * t405;
t346 = rSges(6,3) * t405 + t422 * t403;
t339 = qJD(1) * (-rSges(3,2) * t405 + rSges(3,3) * t403) + t429;
t338 = t396 + (rSges(3,2) * t403 + rSges(3,3) * t405 - t381) * qJD(1);
t337 = rSges(7,3) * t384 + (rSges(7,1) * t404 - rSges(7,2) * t402) * t385;
t336 = Icges(7,5) * t384 + (Icges(7,1) * t404 - Icges(7,4) * t402) * t385;
t335 = Icges(7,6) * t384 + (Icges(7,4) * t404 - Icges(7,2) * t402) * t385;
t334 = Icges(7,3) * t384 + (Icges(7,5) * t404 - Icges(7,6) * t402) * t385;
t329 = qJD(1) * (rSges(4,3) * t405 + t424 * t403) + t429 + t451;
t328 = (-t381 + t424 * t405 + (-rSges(4,3) - qJ(3)) * t403) * qJD(1) + t435;
t327 = rSges(7,1) * t364 + rSges(7,2) * t363 + rSges(7,3) * t441;
t326 = rSges(7,1) * t362 + rSges(7,2) * t361 - rSges(7,3) * t442;
t325 = Icges(7,1) * t364 + Icges(7,4) * t363 + Icges(7,5) * t441;
t324 = Icges(7,1) * t362 + Icges(7,4) * t361 - Icges(7,5) * t442;
t323 = Icges(7,4) * t364 + Icges(7,2) * t363 + Icges(7,6) * t441;
t322 = Icges(7,4) * t362 + Icges(7,2) * t361 - Icges(7,6) * t442;
t321 = Icges(7,5) * t364 + Icges(7,6) * t363 + Icges(7,3) * t441;
t320 = Icges(7,5) * t362 + Icges(7,6) * t361 - Icges(7,3) * t442;
t319 = (-t354 * t403 + t355 * t405) * qJD(4);
t318 = qJD(1) * t354 + (-qJD(4) * t374 - qJD(2)) * t405 + t428;
t317 = t374 * t392 + (-t355 + t426) * qJD(1) + t435;
t316 = qJD(1) * t346 - t369 * t380 + t407;
t315 = t369 * t379 + (-t347 + t409) * qJD(1) + t431;
t314 = -t346 * t379 + t347 * t380 + t425;
t313 = qJD(1) * t356 + t326 * t378 - t337 * t359 - t370 * t380 + t407;
t312 = -t327 * t378 + t337 * t358 + t370 * t379 + (t357 + t409) * qJD(1) + t431;
t311 = -t326 * t358 + t327 * t359 - t356 * t379 - t357 * t380 + t425;
t1 = m(3) * (t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(4) * (t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(5) * (t317 ^ 2 + t318 ^ 2 + t319 ^ 2) / 0.2e1 + t359 * ((-t320 * t442 + t361 * t322 + t362 * t324) * t359 + (-t321 * t442 + t323 * t361 + t325 * t362) * t358 + (-t334 * t442 + t335 * t361 + t336 * t362) * t378) / 0.2e1 + t358 * ((t320 * t441 + t322 * t363 + t324 * t364) * t359 + (t321 * t441 + t363 * t323 + t364 * t325) * t358 + (t334 * t441 + t335 * t363 + t336 * t364) * t378) / 0.2e1 + t378 * ((t320 * t359 + t321 * t358 + t334 * t378) * t384 + ((-t322 * t402 + t324 * t404) * t359 + (-t323 * t402 + t325 * t404) * t358 + (-t335 * t402 + t336 * t404) * t378) * t385) / 0.2e1 + t380 * (t450 * t403 + t408 * t405) / 0.2e1 + t379 * (t408 * t403 - t450 * t405) / 0.2e1 + ((t403 * t371 - t410 * t405) * qJD(1) + (t403 ^ 2 * t349 + (t413 * t405 + (t348 - t412) * t403) * t405) * qJD(4)) * t392 / 0.2e1 + m(7) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + m(6) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + ((t405 * t371 + t410 * t403) * qJD(1) + (t405 ^ 2 * t348 + (t412 * t403 + (t349 - t413) * t405) * t403) * qJD(4)) * t393 / 0.2e1 + ((-t342 * t384 + t344 * t385) * t380 + (-t343 * t384 + t345 * t385) * t379 + ((-t350 * t387 + t352 * t388) * t405 + (-t351 * t387 + t353 * t388) * t403) * qJD(4) + (-t384 * t367 + t385 * t368 - t387 * t372 + t388 * t373) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t382 ^ 2 + t383 ^ 2) + Icges(3,1) + Icges(4,1) * t400 ^ 2 + (-0.2e1 * Icges(4,4) * t400 + Icges(4,2) * t399) * t399 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
