% Calculate kinetic energy for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR13_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR13_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:55
% EndTime: 2019-12-31 18:31:57
% DurationCPUTime: 1.78s
% Computational Cost: add. (783->176), mult. (1008->282), div. (0->0), fcn. (934->8), ass. (0->101)
t418 = Icges(4,4) + Icges(5,6);
t417 = Icges(4,1) + Icges(5,2);
t416 = -Icges(4,2) - Icges(5,3);
t341 = pkin(8) + qJ(3);
t339 = cos(t341);
t415 = t418 * t339;
t338 = sin(t341);
t414 = t418 * t338;
t413 = -Icges(5,4) + Icges(4,5);
t412 = Icges(5,5) - Icges(4,6);
t411 = t416 * t338 + t415;
t410 = -t417 * t339 + t414;
t409 = Icges(5,1) + Icges(4,3);
t346 = sin(qJ(1));
t348 = cos(qJ(1));
t408 = t411 * t346 + t412 * t348;
t407 = -t412 * t346 + t411 * t348;
t406 = t410 * t346 + t413 * t348;
t405 = t413 * t346 - t410 * t348;
t404 = t416 * t339 - t414;
t403 = t417 * t338 + t415;
t402 = t412 * t338 + t413 * t339;
t401 = t402 * t346 - t409 * t348;
t400 = t409 * t346 + t402 * t348;
t399 = t413 * t338 - t412 * t339;
t398 = t404 * t338 + t403 * t339;
t397 = -t407 * t338 + t405 * t339;
t396 = t408 * t338 + t406 * t339;
t343 = cos(pkin(8));
t391 = pkin(2) * t343;
t386 = t339 * t346;
t385 = t339 * t348;
t345 = sin(qJ(5));
t384 = t345 * t346;
t383 = t345 * t348;
t347 = cos(qJ(5));
t382 = t346 * t347;
t381 = t347 * t348;
t333 = pkin(1) * t346 - qJ(2) * t348;
t379 = pkin(6) * t348 - t346 * t391 - t333;
t340 = qJD(2) * t346;
t375 = qJD(4) * t338;
t378 = t348 * t375 + t340;
t377 = qJD(3) * t346;
t376 = qJD(3) * t348;
t374 = qJD(5) * t339;
t364 = pkin(3) * t339 + qJ(4) * t338;
t311 = t364 * t346;
t373 = -t311 + t379;
t325 = pkin(3) * t338 - qJ(4) * t339;
t370 = qJD(3) * (rSges(5,2) * t338 + rSges(5,3) * t339 - t325);
t330 = qJD(1) * (pkin(1) * t348 + qJ(2) * t346);
t369 = -qJD(2) * t348 + qJD(1) * (pkin(6) * t346 + t348 * t391) + t330;
t312 = t364 * t348;
t368 = -qJD(4) * t339 + t311 * t377 + t312 * t376;
t342 = sin(pkin(8));
t367 = rSges(3,1) * t343 - rSges(3,2) * t342;
t366 = rSges(4,1) * t339 - rSges(4,2) * t338;
t365 = -rSges(5,2) * t339 + rSges(5,3) * t338;
t363 = qJD(3) * (-pkin(7) * t338 - t325);
t350 = qJD(1) * t312 + t346 * t375 + t369;
t336 = qJD(5) * t338 + qJD(1);
t335 = rSges(2,1) * t348 - rSges(2,2) * t346;
t334 = rSges(2,1) * t346 + rSges(2,2) * t348;
t329 = -pkin(4) * t348 + pkin(7) * t386;
t328 = pkin(4) * t346 + pkin(7) * t385;
t327 = rSges(4,1) * t338 + rSges(4,2) * t339;
t318 = t346 * t374 - t376;
t317 = t348 * t374 + t377;
t316 = t338 * t384 - t381;
t315 = t338 * t382 + t383;
t314 = t338 * t383 + t382;
t313 = t338 * t381 - t384;
t309 = -rSges(5,1) * t348 + t346 * t365;
t308 = rSges(5,1) * t346 + t348 * t365;
t307 = rSges(4,3) * t346 + t348 * t366;
t306 = -rSges(4,3) * t348 + t346 * t366;
t290 = rSges(6,3) * t338 + (-rSges(6,1) * t345 - rSges(6,2) * t347) * t339;
t289 = Icges(6,5) * t338 + (-Icges(6,1) * t345 - Icges(6,4) * t347) * t339;
t288 = Icges(6,6) * t338 + (-Icges(6,4) * t345 - Icges(6,2) * t347) * t339;
t287 = Icges(6,3) * t338 + (-Icges(6,5) * t345 - Icges(6,6) * t347) * t339;
t285 = qJD(1) * t346 * rSges(3,3) + t330 + (qJD(1) * t367 - qJD(2)) * t348;
t284 = t340 + (t348 * rSges(3,3) - t346 * t367 - t333) * qJD(1);
t283 = rSges(6,1) * t316 + rSges(6,2) * t315 + rSges(6,3) * t386;
t282 = rSges(6,1) * t314 + rSges(6,2) * t313 + rSges(6,3) * t385;
t281 = Icges(6,1) * t316 + Icges(6,4) * t315 + Icges(6,5) * t386;
t280 = Icges(6,1) * t314 + Icges(6,4) * t313 + Icges(6,5) * t385;
t279 = Icges(6,4) * t316 + Icges(6,2) * t315 + Icges(6,6) * t386;
t278 = Icges(6,4) * t314 + Icges(6,2) * t313 + Icges(6,6) * t385;
t277 = Icges(6,5) * t316 + Icges(6,6) * t315 + Icges(6,3) * t386;
t276 = Icges(6,5) * t314 + Icges(6,6) * t313 + Icges(6,3) * t385;
t275 = (t306 * t346 + t307 * t348) * qJD(3);
t274 = qJD(1) * t307 - t327 * t377 + t369;
t273 = -t327 * t376 + t340 + (-t306 + t379) * qJD(1);
t272 = (t308 * t348 + t309 * t346) * qJD(3) + t368;
t271 = qJD(1) * t308 + t346 * t370 + t350;
t270 = t348 * t370 + (-t309 + t373) * qJD(1) + t378;
t269 = qJD(1) * t328 + t282 * t336 - t290 * t317 + t346 * t363 + t350;
t268 = -t283 * t336 + t290 * t318 + t348 * t363 + (-t329 + t373) * qJD(1) + t378;
t267 = -t282 * t318 + t283 * t317 + (t328 * t348 + t329 * t346) * qJD(3) + t368;
t1 = m(3) * (t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(4) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(5) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + m(6) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + t317 * ((t276 * t385 + t313 * t278 + t314 * t280) * t317 + (t277 * t385 + t279 * t313 + t281 * t314) * t318 + (t287 * t385 + t288 * t313 + t289 * t314) * t336) / 0.2e1 + t318 * ((t276 * t386 + t278 * t315 + t280 * t316) * t317 + (t277 * t386 + t315 * t279 + t316 * t281) * t318 + (t287 * t386 + t288 * t315 + t289 * t316) * t336) / 0.2e1 + t336 * ((t276 * t317 + t277 * t318 + t287 * t336) * t338 + ((-t278 * t347 - t280 * t345) * t317 + (-t279 * t347 - t281 * t345) * t318 + (-t288 * t347 - t289 * t345) * t336) * t339) / 0.2e1 + (((t406 * t338 - t408 * t339) * t348 + (t405 * t338 + t407 * t339) * t346) * qJD(3) + (t403 * t338 - t404 * t339) * qJD(1)) * qJD(1) / 0.2e1 + ((t400 * t346 ^ 2 + (t396 * t348 + (t397 - t401) * t346) * t348) * qJD(3) + (t399 * t346 + t398 * t348) * qJD(1)) * t377 / 0.2e1 - ((t401 * t348 ^ 2 + (t397 * t346 + (t396 - t400) * t348) * t346) * qJD(3) + (t398 * t346 - t399 * t348) * qJD(1)) * t376 / 0.2e1 + (m(2) * (t334 ^ 2 + t335 ^ 2) + Icges(2,3) + Icges(3,2) * t343 ^ 2 + (Icges(3,1) * t342 + 0.2e1 * Icges(3,4) * t343) * t342) * qJD(1) ^ 2 / 0.2e1;
T = t1;
