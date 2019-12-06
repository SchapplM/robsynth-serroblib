% Calculate kinetic energy for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:02
% EndTime: 2019-12-05 15:35:03
% DurationCPUTime: 1.22s
% Computational Cost: add. (888->128), mult. (1261->204), div. (0->0), fcn. (1264->8), ass. (0->80)
t421 = Icges(5,1) + Icges(6,1);
t420 = Icges(5,4) - Icges(6,5);
t419 = Icges(6,4) + Icges(5,5);
t418 = Icges(5,2) + Icges(6,3);
t417 = -Icges(5,6) + Icges(6,6);
t416 = Icges(3,3) + Icges(4,3);
t415 = -Icges(5,3) - Icges(6,2);
t343 = qJ(2) + pkin(8);
t341 = sin(t343);
t342 = cos(t343);
t348 = sin(qJ(2));
t350 = cos(qJ(2));
t414 = Icges(3,5) * t350 + Icges(4,5) * t342 - Icges(3,6) * t348 - Icges(4,6) * t341;
t345 = cos(pkin(7));
t394 = t345 ^ 2;
t344 = sin(pkin(7));
t395 = t344 ^ 2;
t397 = t394 + t395;
t396 = qJD(2) * t397;
t413 = rSges(6,1) + pkin(4);
t412 = rSges(6,3) + qJ(5);
t411 = t345 * t344;
t349 = cos(qJ(4));
t382 = t345 * t349;
t347 = sin(qJ(4));
t385 = t344 * t347;
t328 = t342 * t385 + t382;
t383 = t345 * t347;
t384 = t344 * t349;
t329 = t342 * t384 - t383;
t387 = t341 * t344;
t410 = t418 * t328 - t420 * t329 + t417 * t387;
t330 = t342 * t383 - t384;
t331 = t342 * t382 + t385;
t386 = t341 * t345;
t409 = t418 * t330 - t420 * t331 + t417 * t386;
t408 = t417 * t328 + t419 * t329 - t415 * t387;
t407 = t417 * t330 + t419 * t331 - t415 * t386;
t406 = -t420 * t328 + t421 * t329 + t419 * t387;
t405 = -t420 * t330 + t421 * t331 + t419 * t386;
t404 = t417 * t342 + (-t418 * t347 + t420 * t349) * t341;
t403 = t415 * t342 + (t417 * t347 + t419 * t349) * t341;
t402 = t419 * t342 + (t420 * t347 - t421 * t349) * t341;
t401 = t414 * t344 - t416 * t345;
t400 = t416 * t344 + t414 * t345;
t393 = qJD(2) ^ 2;
t390 = pkin(2) * t348;
t389 = t350 * pkin(2);
t381 = rSges(6,2) * t387 + t412 * t328 + t413 * t329;
t380 = -rSges(6,2) * t386 - t412 * t330 - t413 * t331;
t379 = -t342 * rSges(6,2) + (t412 * t347 + t413 * t349) * t341;
t378 = qJD(2) * t344;
t377 = qJD(2) * t345;
t376 = qJD(3) * t345;
t375 = qJD(4) * t341;
t374 = qJD(4) * t342;
t373 = qJD(1) + (-qJ(3) * t345 + t389 * t344) * t378 + (qJ(3) * t344 + t389 * t345) * t377;
t369 = qJD(2) * (-t341 * rSges(4,1) - rSges(4,2) * t342 - t390);
t368 = qJD(2) * (-t341 * pkin(3) + pkin(6) * t342 - t390);
t367 = t373 + (pkin(3) * t342 + pkin(6) * t341) * t396;
t340 = qJD(3) * t344;
t353 = t345 * t368 + t340;
t352 = t344 * t368 - t376;
t337 = t348 * rSges(3,1) + t350 * rSges(3,2);
t333 = t344 * t375 - t377;
t332 = t345 * t375 + t378;
t312 = -t342 * rSges(5,3) + (rSges(5,1) * t349 - rSges(5,2) * t347) * t341;
t302 = t345 * t369 + t340;
t301 = t344 * t369 - t376;
t298 = t331 * rSges(5,1) - t330 * rSges(5,2) + rSges(5,3) * t386;
t296 = t329 * rSges(5,1) - t328 * rSges(5,2) + rSges(5,3) * t387;
t282 = qJD(1) + (rSges(3,1) * t350 - rSges(3,2) * t348) * t396;
t281 = t373 + (rSges(4,1) * t342 - rSges(4,2) * t341) * t396;
t280 = t296 * t374 + t333 * t312 + t353;
t279 = -t298 * t374 - t332 * t312 + t352;
t278 = qJD(5) * t330 + t379 * t333 + t381 * t374 + t353;
t277 = qJD(5) * t328 - t379 * t332 + t380 * t374 + t352;
t276 = t296 * t332 - t298 * t333 + t367;
t275 = qJD(5) * t341 * t347 + t381 * t332 + t380 * t333 + t367;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t337 ^ 2 * t393 * t397 + t282 ^ 2) / 0.2e1 + m(4) * (t281 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + m(5) * (t276 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + m(6) * (t275 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + ((t330 * t404 + t331 * t402 - t386 * t403) * t374 + (t330 * t410 + t406 * t331 + t408 * t386) * t333 + (t409 * t330 + t405 * t331 + t407 * t386) * t332) * t332 / 0.2e1 + ((t328 * t404 + t329 * t402 - t403 * t387) * t374 + (t410 * t328 + t406 * t329 + t408 * t387) * t333 + (t328 * t409 + t329 * t405 + t387 * t407) * t332) * t333 / 0.2e1 - ((-t332 * t407 - t333 * t408 + t374 * t403) * t342 + ((t347 * t404 + t349 * t402) * t374 + (t347 * t410 + t406 * t349) * t333 + (t347 * t409 + t349 * t405) * t332) * t341) * t374 / 0.2e1 + (t400 * t395 - t401 * t411) * t344 * t393 / 0.2e1 - (t401 * t394 - t400 * t411) * t345 * t393 / 0.2e1;
T = t1;
