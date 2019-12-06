% Calculate kinetic energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:19
% EndTime: 2019-12-05 15:40:20
% DurationCPUTime: 1.23s
% Computational Cost: add. (504->134), mult. (1262->211), div. (0->0), fcn. (1265->6), ass. (0->74)
t416 = Icges(4,1) + Icges(3,3);
t415 = Icges(5,1) + Icges(6,1);
t414 = Icges(5,4) - Icges(6,5);
t413 = Icges(6,4) + Icges(5,5);
t412 = Icges(5,2) + Icges(6,3);
t411 = Icges(6,6) - Icges(5,6);
t410 = Icges(5,3) + Icges(6,2);
t349 = sin(qJ(2));
t351 = cos(qJ(2));
t409 = (-Icges(4,4) + Icges(3,5)) * t351 + (Icges(4,5) - Icges(3,6)) * t349;
t347 = cos(pkin(7));
t388 = t347 ^ 2;
t346 = sin(pkin(7));
t389 = t346 ^ 2;
t391 = t388 + t389;
t390 = qJD(2) * t391;
t408 = rSges(6,1) + pkin(4);
t407 = rSges(6,3) + qJ(5);
t406 = t347 * t346;
t348 = sin(qJ(4));
t350 = cos(qJ(4));
t381 = t349 * t350;
t332 = t346 * t348 - t347 * t381;
t382 = t348 * t349;
t333 = t346 * t350 + t347 * t382;
t383 = t347 * t351;
t405 = t332 * t412 - t333 * t414 + t383 * t411;
t334 = t346 * t381 + t347 * t348;
t335 = t346 * t382 - t347 * t350;
t384 = t346 * t351;
t404 = -t334 * t412 - t335 * t414 + t384 * t411;
t403 = t332 * t411 + t333 * t413 + t383 * t410;
t402 = -t334 * t411 + t335 * t413 + t384 * t410;
t401 = -t332 * t414 + t333 * t415 + t383 * t413;
t400 = t334 * t414 + t335 * t415 + t384 * t413;
t399 = t346 * t409 - t347 * t416;
t398 = t346 * t416 + t347 * t409;
t397 = (t348 * t414 + t350 * t412) * t351 + t411 * t349;
t396 = (-t348 * t413 + t350 * t411) * t351 + t410 * t349;
t395 = (-t348 * t415 - t350 * t414) * t351 + t413 * t349;
t376 = qJD(2) * t347;
t377 = qJD(2) * t346;
t394 = (-pkin(3) * t347 + pkin(6) * t384) * t377 + (pkin(3) * t346 + pkin(6) * t383) * t376;
t387 = qJD(2) ^ 2;
t380 = rSges(6,2) * t383 + t332 * t407 + t333 * t408;
t379 = rSges(6,2) * t384 - t334 * t407 + t335 * t408;
t378 = rSges(6,2) * t349 + (-t348 * t408 + t350 * t407) * t351;
t375 = qJD(3) * t349;
t374 = qJD(4) * t349;
t373 = qJD(4) * t351;
t372 = qJD(1) + (pkin(2) * t351 + qJ(3) * t349) * t390;
t340 = pkin(2) * t349 - qJ(3) * t351;
t368 = qJD(2) * (rSges(4,2) * t349 + rSges(4,3) * t351 - t340);
t367 = qJD(2) * (-pkin(6) * t349 - t340);
t356 = -qJD(3) * t351 + t372;
t345 = t347 * t375;
t344 = t346 * t375;
t342 = rSges(3,1) * t349 + rSges(3,2) * t351;
t338 = t346 * t373 - t376;
t337 = t347 * t373 + t377;
t329 = rSges(5,3) * t349 + (-rSges(5,1) * t348 - rSges(5,2) * t350) * t351;
t305 = t347 * t368 + t345;
t304 = t346 * t368 + t344;
t303 = rSges(5,1) * t335 + rSges(5,2) * t334 + rSges(5,3) * t384;
t301 = rSges(5,1) * t333 - rSges(5,2) * t332 + rSges(5,3) * t383;
t287 = qJD(1) + (rSges(3,1) * t351 - rSges(3,2) * t349) * t390;
t286 = t356 + (-rSges(4,2) * t351 + rSges(4,3) * t349) * t390;
t285 = -t303 * t374 + t329 * t338 + t347 * t367 + t345;
t284 = t301 * t374 - t329 * t337 + t346 * t367 + t344;
t283 = -t301 * t338 + t303 * t337 + t356 + t394;
t282 = -t340 * t376 + qJD(5) * t332 + t345 + t378 * t338 + (-pkin(6) * t376 - qJD(4) * t379) * t349;
t281 = -t340 * t377 - qJD(5) * t334 + t344 - t378 * t337 + (-pkin(6) * t377 + qJD(4) * t380) * t349;
t280 = (qJD(5) * t350 - qJD(3)) * t351 - t380 * t338 + t379 * t337 + t372 + t394;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t342 ^ 2 * t387 * t391 + t287 ^ 2) / 0.2e1 + m(4) * (t286 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(5) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(6) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + ((t332 * t397 + t333 * t395 + t383 * t396) * t374 + (t332 * t404 + t333 * t400 + t383 * t402) * t338 + (t405 * t332 + t401 * t333 + t403 * t383) * t337) * t337 / 0.2e1 + ((-t334 * t397 + t335 * t395 + t384 * t396) * t374 + (-t404 * t334 + t400 * t335 + t402 * t384) * t338 + (-t334 * t405 + t335 * t401 + t384 * t403) * t337) * t338 / 0.2e1 + (((-t348 * t395 + t350 * t397) * t374 + (-t348 * t400 + t350 * t404) * t338 + (-t348 * t401 + t350 * t405) * t337) * t351 + (t337 * t403 + t338 * t402 + t374 * t396) * t349) * t374 / 0.2e1 + (t398 * t389 - t399 * t406) * t346 * t387 / 0.2e1 - (t399 * t388 - t398 * t406) * t347 * t387 / 0.2e1;
T = t1;
