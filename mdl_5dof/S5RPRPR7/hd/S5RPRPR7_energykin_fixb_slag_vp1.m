% Calculate kinetic energy for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:52
% EndTime: 2019-12-31 18:18:54
% DurationCPUTime: 1.33s
% Computational Cost: add. (1084->188), mult. (967->297), div. (0->0), fcn. (892->10), ass. (0->106)
t377 = Icges(4,3) + Icges(5,3);
t315 = qJ(3) + pkin(9);
t311 = sin(t315);
t313 = cos(t315);
t319 = sin(qJ(3));
t322 = cos(qJ(3));
t376 = Icges(4,5) * t322 + Icges(5,5) * t313 - Icges(4,6) * t319 - Icges(5,6) * t311;
t316 = qJ(1) + pkin(8);
t312 = sin(t316);
t314 = cos(t316);
t375 = t376 * t312 - t377 * t314;
t374 = t377 * t312 + t376 * t314;
t373 = Icges(4,5) * t319 + Icges(5,5) * t311 + Icges(4,6) * t322 + Icges(5,6) * t313;
t359 = Icges(5,4) * t311;
t296 = Icges(5,2) * t313 + t359;
t358 = Icges(5,4) * t313;
t297 = Icges(5,1) * t311 + t358;
t361 = Icges(4,4) * t319;
t302 = Icges(4,2) * t322 + t361;
t360 = Icges(4,4) * t322;
t303 = Icges(4,1) * t319 + t360;
t372 = -t296 * t311 + t297 * t313 - t302 * t319 + t303 * t322;
t334 = -Icges(5,2) * t311 + t358;
t267 = Icges(5,6) * t312 + t314 * t334;
t336 = Icges(5,1) * t313 - t359;
t269 = Icges(5,5) * t312 + t314 * t336;
t335 = -Icges(4,2) * t319 + t360;
t277 = Icges(4,6) * t312 + t314 * t335;
t337 = Icges(4,1) * t322 - t361;
t280 = Icges(4,5) * t312 + t314 * t337;
t371 = -t267 * t311 + t269 * t313 - t277 * t319 + t280 * t322;
t266 = -Icges(5,6) * t314 + t312 * t334;
t268 = -Icges(5,5) * t314 + t312 * t336;
t276 = -Icges(4,6) * t314 + t312 * t335;
t279 = -Icges(4,5) * t314 + t312 * t337;
t370 = t266 * t311 - t268 * t313 + t276 * t319 - t279 * t322;
t320 = sin(qJ(1));
t366 = pkin(1) * t320;
t365 = pkin(3) * t319;
t363 = pkin(3) * t322;
t357 = t311 * t312;
t356 = t311 * t314;
t318 = sin(qJ(5));
t355 = t312 * t318;
t321 = cos(qJ(5));
t354 = t312 * t321;
t353 = t314 * t318;
t352 = t314 * t321;
t323 = cos(qJ(1));
t310 = qJD(1) * t323 * pkin(1);
t351 = qJD(1) * (pkin(2) * t314 + pkin(6) * t312) + t310;
t350 = qJD(3) * t312;
t349 = qJD(3) * t314;
t348 = qJD(5) * t311;
t262 = -qJ(4) * t314 + t312 * t363;
t263 = qJ(4) * t312 + t314 * t363;
t347 = t262 * t350 + t263 * t349 + qJD(2);
t344 = -pkin(2) * t312 + pkin(6) * t314 - t366;
t343 = -t262 + t344;
t342 = pkin(4) * t313 + pkin(7) * t311;
t341 = rSges(4,1) * t322 - rSges(4,2) * t319;
t340 = rSges(5,1) * t313 - rSges(5,2) * t311;
t339 = qJD(3) * (-rSges(5,1) * t311 - rSges(5,2) * t313 - t365);
t338 = qJD(3) * (-pkin(4) * t311 + pkin(7) * t313 - t365);
t325 = qJD(1) * t263 - qJD(4) * t314 + t351;
t308 = qJD(4) * t312;
t307 = -qJD(5) * t313 + qJD(1);
t306 = rSges(2,1) * t323 - rSges(2,2) * t320;
t305 = rSges(2,1) * t320 + rSges(2,2) * t323;
t304 = rSges(4,1) * t319 + rSges(4,2) * t322;
t293 = t312 * t348 - t349;
t292 = t314 * t348 + t350;
t291 = t313 * t352 + t355;
t290 = -t313 * t353 + t354;
t289 = t313 * t354 - t353;
t288 = -t313 * t355 - t352;
t287 = t310 + qJD(1) * (rSges(3,1) * t314 - rSges(3,2) * t312);
t286 = (-rSges(3,1) * t312 - rSges(3,2) * t314 - t366) * qJD(1);
t285 = t342 * t314;
t284 = t342 * t312;
t283 = rSges(4,3) * t312 + t314 * t341;
t282 = -rSges(4,3) * t314 + t312 * t341;
t281 = -rSges(6,3) * t313 + (rSges(6,1) * t321 - rSges(6,2) * t318) * t311;
t278 = -Icges(6,5) * t313 + (Icges(6,1) * t321 - Icges(6,4) * t318) * t311;
t275 = -Icges(6,6) * t313 + (Icges(6,4) * t321 - Icges(6,2) * t318) * t311;
t272 = -Icges(6,3) * t313 + (Icges(6,5) * t321 - Icges(6,6) * t318) * t311;
t271 = rSges(5,3) * t312 + t314 * t340;
t270 = -rSges(5,3) * t314 + t312 * t340;
t258 = rSges(6,1) * t291 + rSges(6,2) * t290 + rSges(6,3) * t356;
t257 = rSges(6,1) * t289 + rSges(6,2) * t288 + rSges(6,3) * t357;
t256 = Icges(6,1) * t291 + Icges(6,4) * t290 + Icges(6,5) * t356;
t255 = Icges(6,1) * t289 + Icges(6,4) * t288 + Icges(6,5) * t357;
t254 = Icges(6,4) * t291 + Icges(6,2) * t290 + Icges(6,6) * t356;
t253 = Icges(6,4) * t289 + Icges(6,2) * t288 + Icges(6,6) * t357;
t252 = Icges(6,5) * t291 + Icges(6,6) * t290 + Icges(6,3) * t356;
t251 = Icges(6,5) * t289 + Icges(6,6) * t288 + Icges(6,3) * t357;
t250 = qJD(1) * t283 - t304 * t350 + t351;
t249 = -t304 * t349 + (-t282 + t344) * qJD(1);
t248 = qJD(2) + (t282 * t312 + t283 * t314) * qJD(3);
t247 = qJD(1) * t271 + t312 * t339 + t325;
t246 = t308 + t314 * t339 + (-t270 + t343) * qJD(1);
t245 = (t270 * t312 + t271 * t314) * qJD(3) + t347;
t244 = qJD(1) * t285 + t258 * t307 - t281 * t292 + t312 * t338 + t325;
t243 = -t257 * t307 + t281 * t293 + t308 + t314 * t338 + (-t284 + t343) * qJD(1);
t242 = t257 * t292 - t258 * t293 + (t284 * t312 + t285 * t314) * qJD(3) + t347;
t1 = m(3) * (qJD(2) ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + m(4) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(5) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + m(6) * (t242 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + t292 * ((t252 * t356 + t290 * t254 + t291 * t256) * t292 + (t251 * t356 + t253 * t290 + t255 * t291) * t293 + (t272 * t356 + t275 * t290 + t278 * t291) * t307) / 0.2e1 + t293 * ((t252 * t357 + t254 * t288 + t256 * t289) * t292 + (t251 * t357 + t288 * t253 + t289 * t255) * t293 + (t272 * t357 + t275 * t288 + t278 * t289) * t307) / 0.2e1 + t307 * ((-t251 * t293 - t252 * t292 - t272 * t307) * t313 + ((-t254 * t318 + t256 * t321) * t292 + (-t253 * t318 + t255 * t321) * t293 + (-t275 * t318 + t278 * t321) * t307) * t311) / 0.2e1 + (((-t266 * t313 - t268 * t311 - t276 * t322 - t279 * t319) * t314 + (t267 * t313 + t269 * t311 + t277 * t322 + t280 * t319) * t312) * qJD(3) + (t313 * t296 + t311 * t297 + t322 * t302 + t319 * t303) * qJD(1)) * qJD(1) / 0.2e1 + ((t374 * t312 ^ 2 + (t370 * t314 + (t371 - t375) * t312) * t314) * qJD(3) + (t373 * t312 + t372 * t314) * qJD(1)) * t350 / 0.2e1 - ((t375 * t314 ^ 2 + (t371 * t312 + (t370 - t374) * t314) * t312) * qJD(3) + (t372 * t312 - t373 * t314) * qJD(1)) * t349 / 0.2e1 + (m(2) * (t305 ^ 2 + t306 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
