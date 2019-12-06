% Calculate kinetic energy for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:32
% EndTime: 2019-12-05 17:09:33
% DurationCPUTime: 1.42s
% Computational Cost: add. (1128->198), mult. (1335->342), div. (0->0), fcn. (1328->10), ass. (0->104)
t348 = sin(pkin(9));
t349 = cos(pkin(9));
t397 = t348 * t349;
t394 = t349 ^ 2;
t395 = t348 ^ 2;
t396 = t394 + t395;
t393 = qJD(2) ^ 2;
t353 = cos(qJ(2));
t392 = pkin(2) * t353;
t352 = cos(qJ(4));
t391 = pkin(4) * t352;
t347 = qJ(2) + qJ(3);
t343 = sin(t347);
t388 = t343 * t348;
t387 = t343 * t349;
t345 = cos(t347);
t386 = t345 * t348;
t385 = t345 * t349;
t350 = sin(qJ(4));
t384 = t348 * t350;
t383 = t348 * t352;
t382 = t349 * t350;
t381 = t349 * t352;
t341 = qJD(2) * t348;
t332 = qJD(3) * t348 + t341;
t380 = qJD(4) * t343;
t379 = qJD(4) * t345;
t378 = qJD(5) * t343;
t351 = sin(qJ(2));
t377 = pkin(2) * qJD(2) * t351;
t322 = t349 * t380 + t332;
t376 = qJD(1) + (-pkin(6) * t349 + t348 * t392) * t341 + qJD(2) * t349 * (pkin(6) * t348 + t349 * t392);
t333 = (-qJD(2) - qJD(3)) * t349;
t375 = t348 * t377;
t374 = t349 * t377;
t372 = rSges(4,1) * t345 - rSges(4,2) * t343;
t370 = Icges(4,1) * t345 - Icges(4,4) * t343;
t368 = Icges(4,4) * t345 - Icges(4,2) * t343;
t367 = Icges(3,5) * t353 - Icges(3,6) * t351;
t366 = Icges(4,5) * t345 - Icges(4,6) * t343;
t365 = (-Icges(4,3) * t349 + t348 * t366) * t333 + (Icges(4,3) * t348 + t349 * t366) * t332;
t323 = t348 * t380 + t333;
t331 = pkin(3) * t343 - pkin(7) * t345;
t362 = t333 * t331 - t374;
t360 = t376 + (t332 * t348 - t333 * t349) * (pkin(3) * t345 + pkin(7) * t343);
t359 = -t331 * t332 - t375;
t358 = pkin(8) * t343 + t345 * t391;
t357 = (-(Icges(4,6) * t348 + t349 * t368) * t343 + (Icges(4,5) * t348 + t349 * t370) * t345) * t332 + (-(-Icges(4,6) * t349 + t348 * t368) * t343 + (-Icges(4,5) * t349 + t348 * t370) * t345) * t333;
t346 = qJ(4) + qJ(5);
t344 = cos(t346);
t342 = sin(t346);
t337 = rSges(3,1) * t351 + rSges(3,2) * t353;
t330 = rSges(4,1) * t343 + rSges(4,2) * t345;
t329 = (-qJD(4) - qJD(5)) * t345;
t328 = t345 * t381 + t384;
t327 = -t345 * t382 + t383;
t326 = t345 * t383 - t382;
t325 = -t345 * t384 - t381;
t317 = Icges(3,3) * t348 + t349 * t367;
t316 = -Icges(3,3) * t349 + t348 * t367;
t315 = t342 * t348 + t344 * t385;
t314 = -t342 * t385 + t344 * t348;
t313 = -t342 * t349 + t344 * t386;
t312 = -t342 * t386 - t344 * t349;
t304 = -rSges(5,3) * t345 + (rSges(5,1) * t352 - rSges(5,2) * t350) * t343;
t303 = -Icges(5,5) * t345 + (Icges(5,1) * t352 - Icges(5,4) * t350) * t343;
t302 = -Icges(5,6) * t345 + (Icges(5,4) * t352 - Icges(5,2) * t350) * t343;
t301 = -Icges(5,3) * t345 + (Icges(5,5) * t352 - Icges(5,6) * t350) * t343;
t300 = -rSges(6,3) * t345 + (rSges(6,1) * t344 - rSges(6,2) * t342) * t343;
t299 = -Icges(6,5) * t345 + (Icges(6,1) * t344 - Icges(6,4) * t342) * t343;
t298 = -Icges(6,6) * t345 + (Icges(6,4) * t344 - Icges(6,2) * t342) * t343;
t297 = -Icges(6,3) * t345 + (Icges(6,5) * t344 - Icges(6,6) * t342) * t343;
t296 = t348 * t378 + t323;
t295 = t349 * t378 + t322;
t291 = -pkin(8) * t345 + t343 * t391;
t290 = t330 * t333 - t374;
t289 = -t330 * t332 - t375;
t288 = rSges(5,1) * t328 + rSges(5,2) * t327 + rSges(5,3) * t387;
t287 = rSges(5,1) * t326 + rSges(5,2) * t325 + rSges(5,3) * t388;
t286 = Icges(5,1) * t328 + Icges(5,4) * t327 + Icges(5,5) * t387;
t285 = Icges(5,1) * t326 + Icges(5,4) * t325 + Icges(5,5) * t388;
t284 = Icges(5,4) * t328 + Icges(5,2) * t327 + Icges(5,6) * t387;
t283 = Icges(5,4) * t326 + Icges(5,2) * t325 + Icges(5,6) * t388;
t282 = Icges(5,5) * t328 + Icges(5,6) * t327 + Icges(5,3) * t387;
t281 = Icges(5,5) * t326 + Icges(5,6) * t325 + Icges(5,3) * t388;
t280 = pkin(4) * t384 + t349 * t358;
t279 = -pkin(4) * t382 + t348 * t358;
t278 = qJD(1) + t396 * qJD(2) * (rSges(3,1) * t353 - rSges(3,2) * t351);
t277 = rSges(6,1) * t315 + rSges(6,2) * t314 + rSges(6,3) * t387;
t276 = rSges(6,1) * t313 + rSges(6,2) * t312 + rSges(6,3) * t388;
t275 = Icges(6,1) * t315 + Icges(6,4) * t314 + Icges(6,5) * t387;
t274 = Icges(6,1) * t313 + Icges(6,4) * t312 + Icges(6,5) * t388;
t273 = Icges(6,4) * t315 + Icges(6,2) * t314 + Icges(6,6) * t387;
t272 = Icges(6,4) * t313 + Icges(6,2) * t312 + Icges(6,6) * t388;
t271 = Icges(6,5) * t315 + Icges(6,6) * t314 + Icges(6,3) * t387;
t270 = Icges(6,5) * t313 + Icges(6,6) * t312 + Icges(6,3) * t388;
t269 = t332 * (-rSges(4,3) * t349 + t348 * t372) - t333 * (rSges(4,3) * t348 + t349 * t372) + t376;
t268 = t287 * t379 + t304 * t323 + t362;
t267 = -t288 * t379 - t304 * t322 + t359;
t266 = t287 * t322 - t288 * t323 + t360;
t265 = -t276 * t329 + t279 * t379 + t291 * t323 + t296 * t300 + t362;
t264 = t277 * t329 - t280 * t379 - t291 * t322 - t295 * t300 + t359;
t263 = t276 * t295 - t277 * t296 + t279 * t322 - t280 * t323 + t360;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t396 * t393 * t337 ^ 2 + t278 ^ 2) / 0.2e1 + t393 * t348 * (-t316 * t397 + t395 * t317) / 0.2e1 - t393 * t349 * (t394 * t316 - t317 * t397) / 0.2e1 + m(4) * (t269 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + t332 * (t348 * t365 + t349 * t357) / 0.2e1 + t333 * (t348 * t357 - t349 * t365) / 0.2e1 + m(5) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + t322 * ((t282 * t387 + t284 * t327 + t286 * t328) * t322 + (t281 * t387 + t283 * t327 + t285 * t328) * t323 - (t301 * t387 + t302 * t327 + t303 * t328) * t379) / 0.2e1 + t323 * ((t282 * t388 + t284 * t325 + t286 * t326) * t322 + (t281 * t388 + t283 * t325 + t285 * t326) * t323 - (t301 * t388 + t302 * t325 + t303 * t326) * t379) / 0.2e1 - ((-t281 * t323 - t282 * t322 + t301 * t379) * t345 + ((-t284 * t350 + t286 * t352) * t322 + (-t283 * t350 + t285 * t352) * t323 - (-t302 * t350 + t303 * t352) * t379) * t343) * t379 / 0.2e1 + m(6) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t295 * ((t271 * t387 + t314 * t273 + t315 * t275) * t295 + (t270 * t387 + t272 * t314 + t274 * t315) * t296 + (t297 * t387 + t298 * t314 + t299 * t315) * t329) / 0.2e1 + t296 * ((t271 * t388 + t273 * t312 + t275 * t313) * t295 + (t270 * t388 + t312 * t272 + t313 * t274) * t296 + (t297 * t388 + t298 * t312 + t299 * t313) * t329) / 0.2e1 + t329 * ((-t270 * t296 - t271 * t295 - t297 * t329) * t345 + ((-t273 * t342 + t275 * t344) * t295 + (-t272 * t342 + t274 * t344) * t296 + (-t298 * t342 + t299 * t344) * t329) * t343) / 0.2e1;
T = t1;
