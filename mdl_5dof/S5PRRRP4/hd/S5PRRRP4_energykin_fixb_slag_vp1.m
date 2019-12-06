% Calculate kinetic energy for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:35
% EndTime: 2019-12-05 16:45:36
% DurationCPUTime: 1.42s
% Computational Cost: add. (928->146), mult. (1301->243), div. (0->0), fcn. (1304->8), ass. (0->85)
t409 = Icges(5,1) + Icges(6,1);
t408 = Icges(5,4) - Icges(6,5);
t407 = Icges(6,4) + Icges(5,5);
t406 = Icges(5,2) + Icges(6,3);
t405 = -Icges(5,6) + Icges(6,6);
t404 = -Icges(5,3) - Icges(6,2);
t403 = rSges(6,1) + pkin(4);
t402 = rSges(6,3) + qJ(5);
t344 = sin(pkin(8));
t345 = cos(pkin(8));
t401 = t344 * t345;
t343 = qJ(2) + qJ(3);
t342 = cos(t343);
t348 = cos(qJ(4));
t378 = t345 * t348;
t346 = sin(qJ(4));
t381 = t344 * t346;
t327 = t342 * t381 + t378;
t379 = t345 * t346;
t380 = t344 * t348;
t328 = t342 * t380 - t379;
t341 = sin(t343);
t383 = t341 * t344;
t400 = t406 * t327 - t408 * t328 + t405 * t383;
t329 = t342 * t379 - t380;
t330 = t342 * t378 + t381;
t382 = t341 * t345;
t399 = t406 * t329 - t408 * t330 + t405 * t382;
t398 = t405 * t327 + t407 * t328 - t404 * t383;
t397 = t405 * t329 + t407 * t330 - t404 * t382;
t396 = -t408 * t327 + t409 * t328 + t407 * t383;
t395 = -t408 * t329 + t409 * t330 + t407 * t382;
t394 = t405 * t342 + (-t406 * t346 + t408 * t348) * t341;
t393 = t404 * t342 + (t405 * t346 + t407 * t348) * t341;
t392 = t407 * t342 + (t408 * t346 - t409 * t348) * t341;
t389 = t345 ^ 2;
t390 = t344 ^ 2;
t391 = t389 + t390;
t388 = qJD(2) ^ 2;
t349 = cos(qJ(2));
t385 = pkin(2) * t349;
t377 = rSges(6,2) * t383 + t402 * t327 + t403 * t328;
t376 = -rSges(6,2) * t382 - t402 * t329 - t403 * t330;
t375 = -rSges(6,2) * t342 + (t402 * t346 + t403 * t348) * t341;
t340 = qJD(2) * t344;
t333 = qJD(3) * t344 + t340;
t374 = qJD(4) * t341;
t373 = qJD(4) * t342;
t347 = sin(qJ(2));
t372 = pkin(2) * qJD(2) * t347;
t371 = qJD(1) + (-pkin(6) * t345 + t385 * t344) * t340 + qJD(2) * t345 * (pkin(6) * t344 + t385 * t345);
t334 = (-qJD(2) - qJD(3)) * t345;
t369 = t344 * t372;
t368 = t345 * t372;
t366 = rSges(4,1) * t342 - rSges(4,2) * t341;
t364 = Icges(4,1) * t342 - Icges(4,4) * t341;
t362 = Icges(4,4) * t342 - Icges(4,2) * t341;
t361 = Icges(3,5) * t349 - Icges(3,6) * t347;
t360 = Icges(4,5) * t342 - Icges(4,6) * t341;
t359 = (-Icges(4,3) * t345 + t360 * t344) * t334 + (Icges(4,3) * t344 + t360 * t345) * t333;
t332 = pkin(3) * t341 - pkin(7) * t342;
t356 = t334 * t332 - t368;
t354 = t371 + (t333 * t344 - t334 * t345) * (pkin(3) * t342 + pkin(7) * t341);
t353 = -t332 * t333 - t369;
t352 = (-(Icges(4,6) * t344 + t362 * t345) * t341 + (Icges(4,5) * t344 + t364 * t345) * t342) * t333 + (-(-Icges(4,6) * t345 + t362 * t344) * t341 + (-Icges(4,5) * t345 + t364 * t344) * t342) * t334;
t336 = rSges(3,1) * t347 + rSges(3,2) * t349;
t331 = rSges(4,1) * t341 + rSges(4,2) * t342;
t324 = t344 * t374 + t334;
t323 = t345 * t374 + t333;
t318 = Icges(3,3) * t344 + t361 * t345;
t317 = -Icges(3,3) * t345 + t361 * t344;
t309 = -rSges(5,3) * t342 + (rSges(5,1) * t348 - rSges(5,2) * t346) * t341;
t298 = t331 * t334 - t368;
t297 = -t331 * t333 - t369;
t294 = rSges(5,1) * t330 - rSges(5,2) * t329 + rSges(5,3) * t382;
t292 = rSges(5,1) * t328 - rSges(5,2) * t327 + rSges(5,3) * t383;
t278 = qJD(1) + t391 * qJD(2) * (rSges(3,1) * t349 - rSges(3,2) * t347);
t277 = t333 * (-rSges(4,3) * t345 + t366 * t344) - t334 * (rSges(4,3) * t344 + t366 * t345) + t371;
t276 = t292 * t373 + t309 * t324 + t356;
t275 = -t294 * t373 - t309 * t323 + t353;
t274 = qJD(5) * t329 + t375 * t324 + t377 * t373 + t356;
t273 = qJD(5) * t327 - t375 * t323 + t376 * t373 + t353;
t272 = t292 * t323 - t294 * t324 + t354;
t271 = qJD(5) * t341 * t346 + t377 * t323 + t376 * t324 + t354;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t336 ^ 2 * t388 * t391 + t278 ^ 2) / 0.2e1 + t388 * t344 * (-t317 * t401 + t390 * t318) / 0.2e1 - t388 * t345 * (t389 * t317 - t318 * t401) / 0.2e1 + m(4) * (t277 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + t333 * (t359 * t344 + t352 * t345) / 0.2e1 + t334 * (t352 * t344 - t359 * t345) / 0.2e1 + m(5) * (t272 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + m(6) * (t271 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + ((t329 * t394 + t330 * t392 - t382 * t393) * t373 + (t329 * t400 + t330 * t396 + t382 * t398) * t324 + (t329 * t399 + t330 * t395 + t382 * t397) * t323) * t323 / 0.2e1 + ((t327 * t394 + t328 * t392 - t383 * t393) * t373 + (t327 * t400 + t328 * t396 + t383 * t398) * t324 + (t327 * t399 + t328 * t395 + t383 * t397) * t323) * t324 / 0.2e1 - ((-t397 * t323 - t398 * t324 + t393 * t373) * t342 + ((t346 * t394 + t348 * t392) * t373 + (t346 * t400 + t348 * t396) * t324 + (t346 * t399 + t348 * t395) * t323) * t341) * t373 / 0.2e1;
T = t1;
