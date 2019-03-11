% Calculate kinetic energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:41
% EndTime: 2019-03-09 02:02:42
% DurationCPUTime: 1.29s
% Computational Cost: add. (1048->180), mult. (1267->277), div. (0->0), fcn. (1238->8), ass. (0->98)
t431 = Icges(6,1) + Icges(7,1);
t430 = Icges(6,4) - Icges(7,5);
t429 = Icges(7,4) + Icges(6,5);
t428 = Icges(6,2) + Icges(7,3);
t427 = Icges(7,6) - Icges(6,6);
t426 = Icges(6,3) + Icges(7,2);
t425 = rSges(7,1) + pkin(5);
t424 = rSges(7,3) + qJ(6);
t374 = qJ(1) + pkin(9);
t372 = sin(t374);
t373 = cos(t374);
t378 = cos(qJ(5));
t375 = sin(qJ(5));
t376 = sin(qJ(4));
t405 = t375 * t376;
t340 = t372 * t405 - t373 * t378;
t404 = t376 * t378;
t341 = t372 * t404 + t373 * t375;
t379 = cos(qJ(4));
t407 = t372 * t379;
t423 = t428 * t340 - t430 * t341 - t427 * t407;
t342 = t372 * t378 + t373 * t405;
t343 = t372 * t375 - t373 * t404;
t406 = t373 * t379;
t422 = -t428 * t342 - t430 * t343 + t427 * t406;
t421 = t427 * t340 + t429 * t341 - t426 * t407;
t420 = -t427 * t342 + t429 * t343 + t426 * t406;
t419 = -t430 * t340 + t431 * t341 - t429 * t407;
t418 = t430 * t342 + t431 * t343 + t429 * t406;
t417 = (t428 * t375 - t430 * t378) * t379 + t427 * t376;
t416 = (t427 * t375 + t429 * t378) * t379 + t426 * t376;
t415 = (-t430 * t375 + t431 * t378) * t379 + t429 * t376;
t377 = sin(qJ(1));
t410 = pkin(1) * t377;
t409 = Icges(5,4) * t376;
t408 = Icges(5,4) * t379;
t403 = -rSges(7,2) * t407 + t424 * t340 + t425 * t341;
t402 = rSges(7,2) * t406 - t424 * t342 + t425 * t343;
t401 = rSges(7,2) * t376 + (t424 * t375 + t425 * t378) * t379;
t380 = cos(qJ(1));
t371 = qJD(1) * t380 * pkin(1);
t400 = qJD(1) * (pkin(2) * t373 + qJ(3) * t372) + t371;
t399 = qJD(4) * t372;
t398 = qJD(4) * t373;
t397 = qJD(5) * t379;
t396 = qJD(1) * t373 * pkin(7) + t400;
t395 = -pkin(2) * t372 + qJ(3) * t373 - t410;
t394 = pkin(4) * t376 - pkin(8) * t379;
t393 = rSges(5,1) * t376 + rSges(5,2) * t379;
t392 = Icges(5,1) * t376 + t408;
t391 = Icges(5,2) * t379 + t409;
t390 = Icges(5,5) * t376 + Icges(5,6) * t379;
t330 = Icges(5,6) * t373 + t372 * t391;
t332 = Icges(5,5) * t373 + t372 * t392;
t389 = -t330 * t379 - t332 * t376;
t331 = Icges(5,6) * t372 - t373 * t391;
t333 = Icges(5,5) * t372 - t373 * t392;
t388 = t331 * t379 + t333 * t376;
t361 = -Icges(5,2) * t376 + t408;
t362 = Icges(5,1) * t379 - t409;
t387 = t361 * t379 + t362 * t376;
t386 = -pkin(7) * t372 + t395;
t352 = t394 * t372;
t353 = t394 * t373;
t385 = -t352 * t399 - t353 * t398 + qJD(2);
t366 = pkin(4) * t379 + pkin(8) * t376;
t384 = qJD(1) * t352 + (-qJD(4) * t366 - qJD(3)) * t373 + t396;
t369 = qJD(3) * t372;
t383 = t366 * t399 + t369 + (t353 + t386) * qJD(1);
t381 = qJD(2) ^ 2;
t370 = qJD(5) * t376 + qJD(1);
t365 = rSges(2,1) * t380 - rSges(2,2) * t377;
t364 = rSges(5,1) * t379 - rSges(5,2) * t376;
t363 = rSges(2,1) * t377 + rSges(2,2) * t380;
t360 = Icges(5,5) * t379 - Icges(5,6) * t376;
t355 = -t372 * t397 + t398;
t354 = t373 * t397 + t399;
t351 = rSges(6,3) * t376 + (rSges(6,1) * t378 - rSges(6,2) * t375) * t379;
t338 = t371 + qJD(1) * (rSges(3,1) * t373 - rSges(3,2) * t372);
t337 = (-rSges(3,1) * t372 - rSges(3,2) * t373 - t410) * qJD(1);
t335 = rSges(5,3) * t372 - t373 * t393;
t334 = rSges(5,3) * t373 + t372 * t393;
t329 = Icges(5,3) * t372 - t373 * t390;
t328 = Icges(5,3) * t373 + t372 * t390;
t325 = -qJD(3) * t373 + qJD(1) * (-rSges(4,2) * t373 + rSges(4,3) * t372) + t400;
t324 = t369 + (rSges(4,2) * t372 + rSges(4,3) * t373 + t395) * qJD(1);
t323 = rSges(6,1) * t343 + rSges(6,2) * t342 + rSges(6,3) * t406;
t321 = rSges(6,1) * t341 - rSges(6,2) * t340 - rSges(6,3) * t407;
t307 = qJD(2) + (-t334 * t372 + t335 * t373) * qJD(4);
t306 = qJD(1) * t334 + (-qJD(4) * t364 - qJD(3)) * t373 + t396;
t305 = t364 * t399 + t369 + (-t335 + t386) * qJD(1);
t304 = t321 * t370 - t351 * t355 + t384;
t303 = -t323 * t370 + t351 * t354 + t383;
t302 = -t321 * t354 + t323 * t355 + t385;
t301 = -qJD(6) * t342 - t355 * t401 + t370 * t403 + t384;
t300 = qJD(6) * t340 + t354 * t401 - t370 * t402 + t383;
t299 = qJD(6) * t375 * t379 - t354 * t403 + t355 * t402 + t385;
t1 = ((t372 * t360 - t387 * t373) * qJD(1) + (t372 ^ 2 * t329 + (t389 * t373 + (t328 - t388) * t372) * t373) * qJD(4)) * t399 / 0.2e1 + qJD(1) * ((-t376 * t361 + t379 * t362) * qJD(1) + ((-t330 * t376 + t332 * t379) * t373 + (-t331 * t376 + t333 * t379) * t372) * qJD(4)) / 0.2e1 + m(7) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + m(6) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(5) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + m(3) * (t337 ^ 2 + t338 ^ 2 + t381) / 0.2e1 + m(4) * (t324 ^ 2 + t325 ^ 2 + t381) / 0.2e1 + ((t373 * t360 + t387 * t372) * qJD(1) + (t373 ^ 2 * t328 + (t388 * t372 + (t329 - t389) * t373) * t372) * qJD(4)) * t398 / 0.2e1 + ((-t342 * t417 + t343 * t415 + t406 * t416) * t370 + (-t342 * t423 + t419 * t343 + t421 * t406) * t355 + (-t422 * t342 + t418 * t343 + t420 * t406) * t354) * t354 / 0.2e1 + ((t340 * t417 + t341 * t415 - t407 * t416) * t370 + (t423 * t340 + t419 * t341 - t421 * t407) * t355 + (t340 * t422 + t341 * t418 - t407 * t420) * t354) * t355 / 0.2e1 + (((t375 * t417 + t378 * t415) * t370 + (t375 * t423 + t419 * t378) * t355 + (t375 * t422 + t378 * t418) * t354) * t379 + (t354 * t420 + t355 * t421 + t370 * t416) * t376) * t370 / 0.2e1 + (Icges(4,1) + Icges(2,3) + Icges(3,3) + m(2) * (t363 ^ 2 + t365 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
