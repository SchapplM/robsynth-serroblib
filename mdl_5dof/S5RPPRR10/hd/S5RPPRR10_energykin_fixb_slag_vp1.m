% Calculate kinetic energy for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:58
% EndTime: 2019-12-31 18:03:58
% DurationCPUTime: 0.74s
% Computational Cost: add. (597->165), mult. (1036->264), div. (0->0), fcn. (1085->8), ass. (0->86)
t328 = cos(pkin(8));
t359 = t328 ^ 2;
t357 = pkin(3) * t328;
t331 = cos(qJ(4));
t355 = t331 * pkin(4);
t327 = sin(pkin(8));
t329 = sin(qJ(4));
t353 = t327 * t329;
t352 = t328 * t329;
t330 = sin(qJ(1));
t332 = cos(qJ(1));
t317 = t330 * pkin(1) - t332 * qJ(2);
t339 = pkin(2) * t328 + qJ(3) * t327;
t351 = -t339 * t330 - t317;
t323 = qJD(2) * t330;
t349 = qJD(3) * t327;
t350 = t332 * t349 + t323;
t348 = qJD(3) * t328;
t347 = qJD(4) * t330;
t346 = qJD(4) * t332;
t345 = t330 * qJD(1);
t344 = qJD(4) + qJD(5);
t314 = qJD(1) * (t332 * pkin(1) + t330 * qJ(2));
t343 = qJD(1) * t339 * t332 + t330 * t349 + t314;
t342 = -t332 * pkin(6) - t330 * t357 + t351;
t341 = rSges(3,1) * t328 - rSges(3,2) * t327;
t340 = rSges(4,1) * t328 + rSges(4,3) * t327;
t326 = qJ(4) + qJ(5);
t324 = sin(t326);
t325 = cos(t326);
t308 = -t328 * t324 + t327 * t325;
t338 = t327 * t324 + t328 * t325;
t312 = t327 * t331 - t352;
t337 = t328 * t331 + t353;
t336 = -qJD(2) * t332 + qJD(1) * (-t330 * pkin(6) + t332 * t357) + t343;
t335 = pkin(4) * t353 + t355 * t328;
t319 = t332 * rSges(2,1) - t330 * rSges(2,2);
t318 = t330 * rSges(2,1) + t332 * rSges(2,2);
t316 = t344 * t332;
t315 = t344 * t330;
t305 = t337 * t332;
t304 = t312 * t332;
t303 = t337 * t330;
t302 = t312 * t330;
t301 = -pkin(4) * t352 + t355 * t327;
t300 = t338 * t332;
t299 = t308 * t332;
t298 = t338 * t330;
t297 = t308 * t330;
t296 = t312 * rSges(5,1) - rSges(5,2) * t337;
t295 = Icges(5,1) * t312 - Icges(5,4) * t337;
t294 = Icges(5,4) * t312 - Icges(5,2) * t337;
t293 = Icges(5,5) * t312 - Icges(5,6) * t337;
t292 = t308 * rSges(6,1) - rSges(6,2) * t338;
t291 = Icges(6,1) * t308 - Icges(6,4) * t338;
t290 = Icges(6,4) * t308 - Icges(6,2) * t338;
t289 = Icges(6,5) * t308 - Icges(6,6) * t338;
t288 = rSges(3,3) * t345 + t314 + (qJD(1) * t341 - qJD(2)) * t332;
t287 = t323 + (t332 * rSges(3,3) - t341 * t330 - t317) * qJD(1);
t286 = -pkin(7) * t330 + t335 * t332;
t285 = pkin(7) * t332 + t335 * t330;
t284 = t305 * rSges(5,1) + t304 * rSges(5,2) - t330 * rSges(5,3);
t283 = t303 * rSges(5,1) + t302 * rSges(5,2) + t332 * rSges(5,3);
t282 = Icges(5,1) * t305 + Icges(5,4) * t304 - Icges(5,5) * t330;
t281 = Icges(5,1) * t303 + Icges(5,4) * t302 + Icges(5,5) * t332;
t280 = Icges(5,4) * t305 + Icges(5,2) * t304 - Icges(5,6) * t330;
t279 = Icges(5,4) * t303 + Icges(5,2) * t302 + Icges(5,6) * t332;
t278 = Icges(5,5) * t305 + Icges(5,6) * t304 - Icges(5,3) * t330;
t277 = Icges(5,5) * t303 + Icges(5,6) * t302 + Icges(5,3) * t332;
t276 = t300 * rSges(6,1) + t299 * rSges(6,2) - t330 * rSges(6,3);
t275 = t298 * rSges(6,1) + t297 * rSges(6,2) + t332 * rSges(6,3);
t274 = Icges(6,1) * t300 + Icges(6,4) * t299 - Icges(6,5) * t330;
t273 = Icges(6,1) * t298 + Icges(6,4) * t297 + Icges(6,5) * t332;
t272 = Icges(6,4) * t300 + Icges(6,2) * t299 - Icges(6,6) * t330;
t271 = Icges(6,4) * t298 + Icges(6,2) * t297 + Icges(6,6) * t332;
t270 = Icges(6,5) * t300 + Icges(6,6) * t299 - Icges(6,3) * t330;
t269 = Icges(6,5) * t298 + Icges(6,6) * t297 + Icges(6,3) * t332;
t268 = rSges(4,2) * t345 + (qJD(1) * t340 - qJD(2)) * t332 + t343;
t267 = (t332 * rSges(4,2) - t340 * t330 + t351) * qJD(1) + t350;
t266 = -t348 + (-t283 * t330 - t284 * t332) * qJD(4);
t265 = qJD(1) * t284 + t296 * t347 + t336;
t264 = t296 * t346 + (-t283 + t342) * qJD(1) + t350;
t263 = -t348 - t315 * t275 - t316 * t276 + (-t285 * t330 - t286 * t332) * qJD(4);
t262 = t301 * t347 + t315 * t292 + (t276 + t286) * qJD(1) + t336;
t261 = t301 * t346 + t316 * t292 + (-t275 - t285 + t342) * qJD(1) + t350;
t1 = m(3) * (t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t359 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(5) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 - ((-t330 * t293 + t304 * t294 + t305 * t295) * qJD(1) + (-(-t330 * t278 + t304 * t280 + t305 * t282) * t330 + (-t330 * t277 + t304 * t279 + t305 * t281) * t332) * qJD(4)) * t347 / 0.2e1 + ((t332 * t293 + t302 * t294 + t303 * t295) * qJD(1) + (-(t332 * t278 + t302 * t280 + t303 * t282) * t330 + (t332 * t277 + t302 * t279 + t303 * t281) * t332) * qJD(4)) * t346 / 0.2e1 + m(6) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 - t315 * (-(-t330 * t270 + t299 * t272 + t300 * t274) * t315 + (-t330 * t269 + t299 * t271 + t300 * t273) * t316 + (-t330 * t289 + t299 * t290 + t300 * t291) * qJD(1)) / 0.2e1 + t316 * (-(t332 * t270 + t297 * t272 + t298 * t274) * t315 + (t332 * t269 + t297 * t271 + t298 * t273) * t316 + (t332 * t289 + t297 * t290 + t298 * t291) * qJD(1)) / 0.2e1 + ((-(-t280 * t337 + t312 * t282) * t330 + (-t279 * t337 + t312 * t281) * t332) * qJD(4) - (-t272 * t338 + t308 * t274) * t315 + (-t271 * t338 + t308 * t273) * t316 + (-t290 * t338 + t308 * t291 - t294 * t337 + t312 * t295) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t318 ^ 2 + t319 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t359 + ((Icges(3,1) + Icges(4,1)) * t327 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t328) * t327) * qJD(1) ^ 2 / 0.2e1;
T = t1;
