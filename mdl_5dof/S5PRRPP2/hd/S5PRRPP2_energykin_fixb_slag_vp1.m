% Calculate kinetic energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:41
% EndTime: 2019-12-05 16:08:42
% DurationCPUTime: 1.98s
% Computational Cost: add. (1003->178), mult. (1670->271), div. (0->0), fcn. (1716->8), ass. (0->97)
t426 = Icges(5,1) + Icges(6,1);
t425 = Icges(5,4) - Icges(6,5);
t424 = Icges(6,4) + Icges(5,5);
t423 = Icges(5,2) + Icges(6,3);
t422 = -Icges(5,6) + Icges(6,6);
t361 = cos(pkin(7));
t405 = t361 ^ 2;
t360 = sin(pkin(7));
t406 = t360 ^ 2;
t410 = t405 + t406;
t421 = qJD(2) * t410;
t420 = -Icges(5,3) - Icges(6,2) - Icges(4,3);
t419 = rSges(6,1) + pkin(4);
t418 = rSges(6,3) + qJ(5);
t417 = t361 * t360;
t359 = qJ(3) + pkin(8);
t357 = sin(t359);
t358 = cos(t359);
t366 = cos(qJ(2));
t397 = t360 * t366;
t333 = t357 * t397 + t358 * t361;
t334 = -t357 * t361 + t358 * t397;
t364 = sin(qJ(2));
t398 = t360 * t364;
t416 = t423 * t333 - t425 * t334 + t422 * t398;
t394 = t361 * t366;
t335 = t357 * t394 - t360 * t358;
t336 = t357 * t360 + t358 * t394;
t395 = t361 * t364;
t415 = t423 * t335 - t425 * t336 + t422 * t395;
t414 = -t425 * t333 + t426 * t334 + t424 * t398;
t413 = -t425 * t335 + t426 * t336 + t424 * t395;
t412 = t422 * t366 + (-t423 * t357 + t425 * t358) * t364;
t411 = t424 * t366 + (t425 * t357 - t426 * t358) * t364;
t365 = cos(qJ(3));
t363 = sin(qJ(3));
t393 = t363 * t366;
t344 = -t360 * t393 - t361 * t365;
t392 = t365 * t366;
t396 = t361 * t363;
t345 = t360 * t392 - t396;
t409 = Icges(4,5) * t345 + Icges(4,6) * t344 + t422 * t333 + t424 * t334 - t420 * t398;
t346 = t360 * t365 - t361 * t393;
t399 = t360 * t363;
t347 = t361 * t392 + t399;
t408 = Icges(4,5) * t347 + Icges(4,6) * t346 + t422 * t335 + t424 * t336 - t420 * t395;
t407 = t420 * t366 + (Icges(4,5) * t365 - Icges(4,6) * t363 + t422 * t357 + t424 * t358) * t364;
t404 = qJD(2) ^ 2;
t401 = pkin(3) * t365;
t391 = rSges(6,2) * t398 + t418 * t333 + t419 * t334;
t368 = qJ(4) * t364 + t366 * t401;
t314 = pkin(3) * t399 + t361 * t368;
t390 = -rSges(5,1) * t336 + rSges(5,2) * t335 - rSges(5,3) * t395 - t314;
t389 = -rSges(6,2) * t366 + (t418 * t357 + t419 * t358) * t364;
t388 = qJD(2) * t360;
t387 = qJD(2) * t361;
t386 = qJD(3) * t364;
t385 = qJD(3) * t366;
t384 = qJD(4) * t364;
t383 = -rSges(6,2) * t395 - t418 * t335 - t419 * t336 - t314;
t352 = pkin(2) * t364 - pkin(6) * t366;
t382 = t352 * t388;
t381 = t352 * t387;
t380 = qJD(1) + (pkin(2) * t366 + pkin(6) * t364) * t421;
t378 = t360 * t384 - t382;
t375 = Icges(3,5) * t366 - Icges(3,6) * t364;
t313 = -pkin(3) * t396 + t360 * t368;
t318 = -qJ(4) * t366 + t364 * t401;
t349 = t360 * t386 - t387;
t370 = t313 * t385 + t349 * t318 + t361 * t384 - t381;
t348 = t361 * t386 + t388;
t369 = -qJD(4) * t366 + t348 * t313 + t380;
t351 = rSges(3,1) * t364 + rSges(3,2) * t366;
t343 = -rSges(4,3) * t366 + (rSges(4,1) * t365 - rSges(4,2) * t363) * t364;
t341 = -Icges(4,5) * t366 + (Icges(4,1) * t365 - Icges(4,4) * t363) * t364;
t340 = -Icges(4,6) * t366 + (Icges(4,4) * t365 - Icges(4,2) * t363) * t364;
t328 = Icges(3,3) * t360 + t361 * t375;
t327 = -Icges(3,3) * t361 + t360 * t375;
t326 = -rSges(5,3) * t366 + (rSges(5,1) * t358 - rSges(5,2) * t357) * t364;
t316 = rSges(4,1) * t347 + rSges(4,2) * t346 + rSges(4,3) * t395;
t315 = rSges(4,1) * t345 + rSges(4,2) * t344 + rSges(4,3) * t398;
t312 = Icges(4,1) * t347 + Icges(4,4) * t346 + Icges(4,5) * t395;
t311 = Icges(4,1) * t345 + Icges(4,4) * t344 + Icges(4,5) * t398;
t310 = Icges(4,4) * t347 + Icges(4,2) * t346 + Icges(4,6) * t395;
t309 = Icges(4,4) * t345 + Icges(4,2) * t344 + Icges(4,6) * t398;
t303 = qJD(1) + (rSges(3,1) * t366 - rSges(3,2) * t364) * t421;
t300 = rSges(5,1) * t334 - rSges(5,2) * t333 + rSges(5,3) * t398;
t285 = t315 * t385 + t343 * t349 - t381;
t284 = -t316 * t385 - t343 * t348 - t382;
t283 = t315 * t348 - t316 * t349 + t380;
t282 = t300 * t385 + t326 * t349 + t370;
t281 = (-t318 - t326) * t348 + t390 * t385 + t378;
t280 = t300 * t348 + t349 * t390 + t369;
t279 = qJD(5) * t335 + t349 * t389 + t385 * t391 + t370;
t278 = qJD(5) * t333 + (-t318 - t389) * t348 + t383 * t385 + t378;
t277 = qJD(5) * t357 * t364 + t348 * t391 + t349 * t383 + t369;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t351 ^ 2 * t404 * t410 + t303 ^ 2) / 0.2e1 + t404 * t360 * (-t327 * t417 + t406 * t328) / 0.2e1 - t404 * t361 * (t405 * t327 - t328 * t417) / 0.2e1 + m(4) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(5) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(6) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + ((t335 * t412 + t336 * t411 - t340 * t346 - t341 * t347 - t395 * t407) * t385 + (t309 * t346 + t311 * t347 + t335 * t416 + t336 * t414 + t395 * t409) * t349 + (t346 * t310 + t347 * t312 + t415 * t335 + t413 * t336 + t408 * t395) * t348) * t348 / 0.2e1 + ((t333 * t412 + t334 * t411 - t340 * t344 - t345 * t341 - t398 * t407) * t385 + (t344 * t309 + t345 * t311 + t416 * t333 + t414 * t334 + t409 * t398) * t349 + (t310 * t344 + t312 * t345 + t333 * t415 + t334 * t413 + t398 * t408) * t348) * t349 / 0.2e1 - ((-t348 * t408 - t349 * t409 + t385 * t407) * t366 + ((t340 * t363 - t341 * t365 + t357 * t412 + t358 * t411) * t385 + (-t309 * t363 + t311 * t365 + t357 * t416 + t358 * t414) * t349 + (-t310 * t363 + t312 * t365 + t357 * t415 + t358 * t413) * t348) * t364) * t385 / 0.2e1;
T = t1;
