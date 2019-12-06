% Calculate kinetic energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:01
% EndTime: 2019-12-05 19:00:02
% DurationCPUTime: 0.85s
% Computational Cost: add. (1062->169), mult. (788->278), div. (0->0), fcn. (658->10), ass. (0->106)
t316 = qJ(1) + qJ(2);
t308 = sin(t316);
t310 = cos(t316);
t315 = qJ(3) + qJ(4);
t311 = qJ(5) + t315;
t304 = sin(t311);
t305 = cos(t311);
t356 = Icges(6,4) * t305;
t338 = -Icges(6,2) * t304 + t356;
t254 = Icges(6,6) * t310 - t308 * t338;
t255 = Icges(6,6) * t308 + t310 * t338;
t357 = Icges(6,4) * t304;
t341 = Icges(6,1) * t305 - t357;
t256 = Icges(6,5) * t310 - t308 * t341;
t257 = Icges(6,5) * t308 + t310 * t341;
t302 = qJD(3) * t308;
t285 = qJD(4) * t308 + t302;
t278 = qJD(5) * t308 + t285;
t303 = qJD(3) * t310;
t286 = qJD(4) * t310 + t303;
t279 = qJD(5) * t310 + t286;
t282 = Icges(6,2) * t305 + t357;
t283 = Icges(6,1) * t304 + t356;
t313 = qJD(1) + qJD(2);
t369 = (t282 * t304 - t283 * t305) * t313 + (t254 * t304 - t256 * t305) * t279 + (t255 * t304 - t257 * t305) * t278;
t307 = sin(t315);
t309 = cos(t315);
t358 = Icges(5,4) * t309;
t339 = -Icges(5,2) * t307 + t358;
t262 = Icges(5,6) * t310 - t308 * t339;
t263 = Icges(5,6) * t308 + t310 * t339;
t359 = Icges(5,4) * t307;
t342 = Icges(5,1) * t309 - t359;
t264 = Icges(5,5) * t310 - t308 * t342;
t265 = Icges(5,5) * t308 + t310 * t342;
t288 = Icges(5,2) * t309 + t359;
t289 = Icges(5,1) * t307 + t358;
t368 = (t288 * t307 - t289 * t309) * t313 + (t262 * t307 - t264 * t309) * t286 + (t263 * t307 - t265 * t309) * t285;
t365 = pkin(4) * t307;
t319 = cos(qJ(3));
t364 = t319 * pkin(3);
t362 = pkin(1) * qJD(1);
t317 = sin(qJ(3));
t361 = Icges(4,4) * t317;
t360 = Icges(4,4) * t319;
t251 = pkin(8) * t308 + t310 * t364;
t291 = t310 * pkin(2) + t308 * pkin(7);
t355 = -t251 - t291;
t354 = pkin(4) * t309;
t352 = pkin(3) * qJD(3) * t317;
t318 = sin(qJ(1));
t351 = t318 * t362;
t320 = cos(qJ(1));
t350 = t320 * t362;
t349 = t313 * (-t308 * pkin(2) + t310 * pkin(7)) - t351;
t348 = t308 * t352 - t350;
t250 = pkin(8) * t310 - t308 * t364;
t347 = -t250 * t302 + t251 * t303;
t346 = rSges(4,1) * t319 - rSges(4,2) * t317;
t345 = rSges(5,1) * t309 - rSges(5,2) * t307;
t344 = rSges(6,1) * t305 - rSges(6,2) * t304;
t343 = Icges(4,1) * t319 - t361;
t340 = -Icges(4,2) * t317 + t360;
t337 = Icges(4,5) * t319 - Icges(4,6) * t317;
t336 = Icges(5,5) * t309 - Icges(5,6) * t307;
t335 = Icges(6,5) * t305 - Icges(6,6) * t304;
t270 = Icges(4,6) * t310 - t308 * t340;
t272 = Icges(4,5) * t310 - t308 * t343;
t330 = -t270 * t317 + t272 * t319;
t271 = Icges(4,6) * t308 + t310 * t340;
t273 = Icges(4,5) * t308 + t310 * t343;
t329 = t271 * t317 - t273 * t319;
t295 = Icges(4,2) * t319 + t361;
t296 = Icges(4,1) * t317 + t360;
t326 = t295 * t317 - t296 * t319;
t325 = (Icges(6,3) * t310 - t308 * t335) * t279 + (Icges(6,3) * t308 + t310 * t335) * t278 + (Icges(6,5) * t304 + Icges(6,6) * t305) * t313;
t324 = (Icges(5,3) * t310 - t308 * t336) * t286 + (Icges(5,3) * t308 + t310 * t336) * t285 + (Icges(5,5) * t307 + Icges(5,6) * t309) * t313;
t323 = t313 * t250 - t310 * t352 + t349;
t299 = t320 * rSges(2,1) - t318 * rSges(2,2);
t298 = -t318 * rSges(2,1) - t320 * rSges(2,2);
t297 = t317 * rSges(4,1) + t319 * rSges(4,2);
t294 = Icges(4,5) * t317 + Icges(4,6) * t319;
t290 = t307 * rSges(5,1) + t309 * rSges(5,2);
t284 = t304 * rSges(6,1) + t305 * rSges(6,2);
t277 = -t350 - t313 * (t310 * rSges(3,1) - t308 * rSges(3,2));
t276 = -t351 + t313 * (-t308 * rSges(3,1) - t310 * rSges(3,2));
t275 = t308 * rSges(4,3) + t310 * t346;
t274 = t310 * rSges(4,3) - t308 * t346;
t269 = Icges(4,3) * t308 + t310 * t337;
t268 = Icges(4,3) * t310 - t308 * t337;
t267 = t308 * rSges(5,3) + t310 * t345;
t266 = t310 * rSges(5,3) - t308 * t345;
t259 = t308 * rSges(6,3) + t310 * t344;
t258 = t310 * rSges(6,3) - t308 * t344;
t247 = pkin(9) * t308 + t310 * t354;
t246 = pkin(9) * t310 - t308 * t354;
t245 = (-t274 * t308 + t275 * t310) * qJD(3);
t244 = -t350 + t297 * t302 + (-t275 - t291) * t313;
t243 = t313 * t274 - t297 * t303 + t349;
t242 = t285 * t290 + (-t267 + t355) * t313 + t348;
t241 = t313 * t266 - t286 * t290 + t323;
t240 = -t285 * t266 + t286 * t267 + t347;
t239 = t285 * t365 + t278 * t284 + (-t247 - t259 + t355) * t313 + t348;
t238 = -t286 * t365 - t279 * t284 + (t246 + t258) * t313 + t323;
t237 = -t285 * t246 + t286 * t247 - t278 * t258 + t279 * t259 + t347;
t1 = m(3) * (t276 ^ 2 + t277 ^ 2) / 0.2e1 + t313 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + ((t310 * t294 + t308 * t326) * t313 + (t310 ^ 2 * t268 + (t329 * t308 + (t269 - t330) * t310) * t308) * qJD(3)) * t303 / 0.2e1 + ((t308 * t294 - t310 * t326) * t313 + (t308 ^ 2 * t269 + (t330 * t310 + (t268 - t329) * t308) * t310) * qJD(3)) * t302 / 0.2e1 + m(5) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + t286 * (t368 * t308 + t324 * t310) / 0.2e1 + t285 * (t324 * t308 - t368 * t310) / 0.2e1 + m(6) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + t279 * (t369 * t308 + t325 * t310) / 0.2e1 + t278 * (t325 * t308 - t369 * t310) / 0.2e1 + (m(2) * (t298 ^ 2 + t299 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t319 * t270 + t317 * t272) * t310 + (t319 * t271 + t317 * t273) * t308) * qJD(3) + (t309 * t262 + t307 * t264) * t286 + (t309 * t263 + t307 * t265) * t285 + (t305 * t254 + t304 * t256) * t279 + (t305 * t255 + t304 * t257) * t278 + (t305 * t282 + t304 * t283 + t309 * t288 + t307 * t289 + t319 * t295 + t317 * t296) * t313) * t313 / 0.2e1;
T = t1;
