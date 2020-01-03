% Calculate kinetic energy for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:53
% EndTime: 2019-12-31 18:50:55
% DurationCPUTime: 1.34s
% Computational Cost: add. (996->177), mult. (1277->280), div. (0->0), fcn. (1259->8), ass. (0->98)
t405 = Icges(5,1) + Icges(6,1);
t404 = Icges(5,4) + Icges(6,4);
t403 = -Icges(6,5) - Icges(5,5);
t402 = Icges(5,2) + Icges(6,2);
t401 = -Icges(6,6) - Icges(5,6);
t400 = -Icges(6,3) - Icges(5,3);
t343 = pkin(8) + qJ(3);
t341 = cos(t343);
t350 = cos(qJ(4));
t351 = cos(qJ(1));
t376 = t350 * t351;
t348 = sin(qJ(4));
t349 = sin(qJ(1));
t379 = t348 * t349;
t322 = -t341 * t379 - t376;
t377 = t349 * t350;
t378 = t348 * t351;
t323 = t341 * t377 - t378;
t340 = sin(t343);
t381 = t340 * t349;
t399 = -t401 * t322 - t403 * t323 - t400 * t381;
t324 = -t341 * t378 + t377;
t325 = t341 * t376 + t379;
t380 = t340 * t351;
t398 = -t401 * t324 - t403 * t325 - t400 * t380;
t397 = t402 * t322 + t404 * t323 - t401 * t381;
t396 = t402 * t324 + t404 * t325 - t401 * t380;
t395 = t404 * t322 + t405 * t323 - t403 * t381;
t394 = t404 * t324 + t405 * t325 - t403 * t380;
t393 = t400 * t341 + (t401 * t348 - t403 * t350) * t340;
t392 = t401 * t341 + (-t402 * t348 + t404 * t350) * t340;
t391 = t403 * t341 + (-t404 * t348 + t405 * t350) * t340;
t345 = cos(pkin(8));
t386 = pkin(2) * t345;
t385 = pkin(4) * t350;
t383 = Icges(4,4) * t340;
t382 = Icges(4,4) * t341;
t353 = qJ(5) * t340 + t385 * t341;
t374 = rSges(6,1) * t323 + rSges(6,2) * t322 + rSges(6,3) * t381 - pkin(4) * t378 + t353 * t349;
t373 = rSges(6,1) * t325 + rSges(6,2) * t324 + rSges(6,3) * t380 + pkin(4) * t379 + t353 * t351;
t372 = (-qJ(5) - rSges(6,3)) * t341 + (rSges(6,1) * t350 - rSges(6,2) * t348 + t385) * t340;
t334 = pkin(1) * t349 - qJ(2) * t351;
t371 = pkin(6) * t351 - t386 * t349 - t334;
t366 = pkin(3) * t341 + pkin(7) * t340;
t320 = t366 * t349;
t321 = t366 * t351;
t368 = qJD(3) * t351;
t369 = qJD(3) * t349;
t370 = t320 * t369 + t321 * t368;
t367 = qJD(4) * t340;
t333 = qJD(1) * (pkin(1) * t351 + qJ(2) * t349);
t365 = -qJD(2) * t351 + qJD(1) * (pkin(6) * t349 + t386 * t351) + t333;
t344 = sin(pkin(8));
t364 = rSges(3,1) * t345 - rSges(3,2) * t344;
t363 = rSges(4,1) * t341 - rSges(4,2) * t340;
t362 = Icges(4,1) * t341 - t383;
t361 = -Icges(4,2) * t340 + t382;
t360 = Icges(4,5) * t341 - Icges(4,6) * t340;
t311 = -Icges(4,6) * t351 + t361 * t349;
t313 = -Icges(4,5) * t351 + t362 * t349;
t359 = t311 * t340 - t313 * t341;
t312 = Icges(4,6) * t349 + t361 * t351;
t314 = Icges(4,5) * t349 + t362 * t351;
t358 = -t312 * t340 + t314 * t341;
t329 = Icges(4,2) * t341 + t383;
t330 = Icges(4,1) * t340 + t382;
t357 = -t329 * t340 + t330 * t341;
t332 = pkin(3) * t340 - pkin(7) * t341;
t356 = -qJD(3) * t332 + qJD(5) * t340;
t355 = qJD(1) * t321 + t365;
t342 = qJD(2) * t349;
t354 = t342 + (-t320 + t371) * qJD(1);
t337 = -qJD(4) * t341 + qJD(1);
t336 = rSges(2,1) * t351 - rSges(2,2) * t349;
t335 = rSges(2,1) * t349 + rSges(2,2) * t351;
t331 = rSges(4,1) * t340 + rSges(4,2) * t341;
t328 = Icges(4,5) * t340 + Icges(4,6) * t341;
t327 = t349 * t367 - t368;
t326 = t351 * t367 + t369;
t316 = rSges(4,3) * t349 + t363 * t351;
t315 = -rSges(4,3) * t351 + t363 * t349;
t310 = Icges(4,3) * t349 + t360 * t351;
t309 = -Icges(4,3) * t351 + t360 * t349;
t307 = -rSges(5,3) * t341 + (rSges(5,1) * t350 - rSges(5,2) * t348) * t340;
t297 = qJD(1) * t349 * rSges(3,3) + t333 + (qJD(1) * t364 - qJD(2)) * t351;
t296 = t342 + (t351 * rSges(3,3) - t364 * t349 - t334) * qJD(1);
t295 = rSges(5,1) * t325 + rSges(5,2) * t324 + rSges(5,3) * t380;
t293 = rSges(5,1) * t323 + rSges(5,2) * t322 + rSges(5,3) * t381;
t277 = (t315 * t349 + t316 * t351) * qJD(3);
t276 = qJD(1) * t316 - t331 * t369 + t365;
t275 = -t331 * t368 + t342 + (-t315 + t371) * qJD(1);
t274 = t293 * t326 - t295 * t327 + t370;
t273 = t295 * t337 - t307 * t326 - t332 * t369 + t355;
t272 = -t293 * t337 + t307 * t327 - t332 * t368 + t354;
t271 = -t372 * t326 + t373 * t337 + t356 * t349 + t355;
t270 = t372 * t327 - t374 * t337 + t356 * t351 + t354;
t269 = -qJD(5) * t341 + t374 * t326 - t373 * t327 + t370;
t1 = m(3) * (t296 ^ 2 + t297 ^ 2) / 0.2e1 + m(4) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + ((t349 * t328 + t357 * t351) * qJD(1) + (t349 ^ 2 * t310 + (t359 * t351 + (-t309 + t358) * t349) * t351) * qJD(3)) * t369 / 0.2e1 - ((-t351 * t328 + t357 * t349) * qJD(1) + (t351 ^ 2 * t309 + (t358 * t349 + (-t310 + t359) * t351) * t349) * qJD(3)) * t368 / 0.2e1 + qJD(1) * ((t341 * t329 + t340 * t330) * qJD(1) + ((t312 * t341 + t314 * t340) * t349 - (t311 * t341 + t313 * t340) * t351) * qJD(3)) / 0.2e1 + m(5) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + m(6) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + ((t324 * t392 + t325 * t391 + t380 * t393) * t337 + (t397 * t324 + t395 * t325 + t380 * t399) * t327 + (t396 * t324 + t394 * t325 + t398 * t380) * t326) * t326 / 0.2e1 + ((t322 * t392 + t323 * t391 + t381 * t393) * t337 + (t397 * t322 + t395 * t323 + t399 * t381) * t327 + (t322 * t396 + t323 * t394 + t381 * t398) * t326) * t327 / 0.2e1 + ((-t398 * t326 - t327 * t399 - t393 * t337) * t341 + ((-t348 * t392 + t350 * t391) * t337 + (-t348 * t397 + t350 * t395) * t327 + (-t348 * t396 + t350 * t394) * t326) * t340) * t337 / 0.2e1 + (m(2) * (t335 ^ 2 + t336 ^ 2) + Icges(2,3) + Icges(3,2) * t345 ^ 2 + (Icges(3,1) * t344 + 0.2e1 * Icges(3,4) * t345) * t344) * qJD(1) ^ 2 / 0.2e1;
T = t1;
