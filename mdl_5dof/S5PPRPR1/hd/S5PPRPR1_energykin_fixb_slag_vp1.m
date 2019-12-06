% Calculate kinetic energy for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:51
% EndTime: 2019-12-05 15:00:52
% DurationCPUTime: 0.92s
% Computational Cost: add. (779->140), mult. (864->234), div. (0->0), fcn. (874->8), ass. (0->79)
t323 = cos(pkin(7));
t359 = t323 ^ 2;
t321 = sin(pkin(7));
t360 = t321 ^ 2;
t362 = t359 + t360;
t361 = qJD(3) * t362;
t322 = cos(pkin(9));
t357 = pkin(4) * t322;
t319 = pkin(8) + qJ(3);
t315 = sin(t319);
t356 = t315 * t321;
t355 = t315 * t323;
t317 = cos(t319);
t354 = t317 * t321;
t353 = t317 * t323;
t320 = sin(pkin(9));
t352 = t320 * t321;
t351 = t320 * t323;
t350 = t321 * t322;
t349 = t322 * t323;
t313 = qJD(2) * t321;
t343 = qJD(4) * t315;
t347 = t323 * t343 + t313;
t346 = qJD(2) * t323;
t345 = qJD(3) * t321;
t344 = qJD(3) * t323;
t342 = qJD(5) * t315;
t341 = qJD(5) * t317;
t308 = pkin(3) * t315 - qJ(4) * t317;
t338 = qJD(3) * (pkin(6) * t317 - t357 * t315 - t308);
t337 = qJD(3) * (rSges(5,3) * t317 - (rSges(5,1) * t322 - rSges(5,2) * t320) * t315 - t308);
t336 = t321 * t343 - t346;
t333 = Icges(4,5) * t317 - Icges(4,6) * t315;
t330 = -qJD(4) * t317 + qJD(1) + (pkin(3) * t317 + qJ(4) * t315) * t361;
t326 = qJD(1) ^ 2;
t318 = pkin(9) + qJ(5);
t316 = cos(t318);
t314 = sin(t318);
t309 = rSges(4,1) * t315 + rSges(4,2) * t317;
t307 = t321 * t342 - t344;
t306 = t323 * t342 + t345;
t305 = t317 * t349 + t352;
t304 = -t317 * t351 + t350;
t303 = t317 * t350 - t351;
t302 = -t317 * t352 - t349;
t301 = t314 * t321 + t316 * t353;
t300 = -t314 * t353 + t316 * t321;
t299 = -t314 * t323 + t316 * t354;
t298 = -t314 * t354 - t316 * t323;
t297 = -t309 * t345 - t346;
t296 = -t309 * t344 + t313;
t289 = Icges(4,3) * t321 + t333 * t323;
t288 = -Icges(4,3) * t323 + t333 * t321;
t286 = -rSges(6,3) * t317 + (rSges(6,1) * t316 - rSges(6,2) * t314) * t315;
t285 = -Icges(6,5) * t317 + (Icges(6,1) * t316 - Icges(6,4) * t314) * t315;
t284 = -Icges(6,6) * t317 + (Icges(6,4) * t316 - Icges(6,2) * t314) * t315;
t283 = -Icges(6,3) * t317 + (Icges(6,5) * t316 - Icges(6,6) * t314) * t315;
t281 = Icges(5,1) * t305 + Icges(5,4) * t304 + Icges(5,5) * t355;
t280 = Icges(5,1) * t303 + Icges(5,4) * t302 + Icges(5,5) * t356;
t279 = Icges(5,4) * t305 + Icges(5,2) * t304 + Icges(5,6) * t355;
t278 = Icges(5,4) * t303 + Icges(5,2) * t302 + Icges(5,6) * t356;
t277 = Icges(5,5) * t305 + Icges(5,6) * t304 + Icges(5,3) * t355;
t276 = Icges(5,5) * t303 + Icges(5,6) * t302 + Icges(5,3) * t356;
t275 = rSges(6,1) * t301 + rSges(6,2) * t300 + rSges(6,3) * t355;
t274 = rSges(6,1) * t299 + rSges(6,2) * t298 + rSges(6,3) * t356;
t273 = Icges(6,1) * t301 + Icges(6,4) * t300 + Icges(6,5) * t355;
t272 = Icges(6,1) * t299 + Icges(6,4) * t298 + Icges(6,5) * t356;
t271 = Icges(6,4) * t301 + Icges(6,2) * t300 + Icges(6,6) * t355;
t270 = Icges(6,4) * t299 + Icges(6,2) * t298 + Icges(6,6) * t356;
t269 = Icges(6,5) * t301 + Icges(6,6) * t300 + Icges(6,3) * t355;
t268 = Icges(6,5) * t299 + Icges(6,6) * t298 + Icges(6,3) * t356;
t267 = qJD(1) + (rSges(4,1) * t317 - rSges(4,2) * t315) * t361;
t266 = t321 * t337 + t336;
t265 = t323 * t337 + t347;
t264 = (t321 * (rSges(5,1) * t303 + rSges(5,2) * t302 + rSges(5,3) * t356) + t323 * (rSges(5,1) * t305 + rSges(5,2) * t304 + rSges(5,3) * t355)) * qJD(3) + t330;
t263 = -t275 * t341 - t286 * t306 + t321 * t338 + t336;
t262 = t274 * t341 + t286 * t307 + t323 * t338 + t347;
t261 = t306 * t274 - t307 * t275 + t330 + (pkin(6) * t315 + t357 * t317) * t361;
t1 = m(2) * t326 / 0.2e1 + m(3) * (qJD(2) ^ 2 * t362 + t326) / 0.2e1 + m(4) * (t267 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + m(5) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(6) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t306 * ((t269 * t355 + t271 * t300 + t273 * t301) * t306 + (t268 * t355 + t270 * t300 + t301 * t272) * t307 - (t283 * t355 + t284 * t300 + t285 * t301) * t341) / 0.2e1 + t307 * ((t269 * t356 + t271 * t298 + t273 * t299) * t306 + (t268 * t356 + t270 * t298 + t272 * t299) * t307 - (t283 * t356 + t284 * t298 + t285 * t299) * t341) / 0.2e1 - ((-t268 * t307 - t269 * t306 + t283 * t341) * t317 + ((-t271 * t314 + t273 * t316) * t306 + (-t270 * t314 + t272 * t316) * t307 - (-t284 * t314 + t285 * t316) * t341) * t315) * t341 / 0.2e1 + ((t360 * t289 + (t277 * t355 + t279 * t304 + t281 * t305) * t321 + (-t276 * t355 - t278 * t304 - t280 * t305 - t321 * t288) * t323) * t321 / 0.2e1 - (t359 * t288 - (t276 * t356 + t278 * t302 + t280 * t303) * t323 + (t277 * t356 + t279 * t302 + t281 * t303 - t323 * t289) * t321) * t323 / 0.2e1) * qJD(3) ^ 2;
T = t1;
