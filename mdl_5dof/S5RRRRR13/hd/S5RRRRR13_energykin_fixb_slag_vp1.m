% Calculate kinetic energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR13_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:19
% EndTime: 2024-09-27 17:30:20
% DurationCPUTime: 0.19s
% Computational Cost: add. (1647->172), mult. (1071->284), div. (0->0), fcn. (994->16), ass. (0->97)
t327 = sin(pkin(5));
t326 = qJ(1) + qJ(2);
t323 = qJ(3) + t326;
t311 = sin(t323);
t312 = cos(t323);
t329 = sin(qJ(4));
t328 = cos(pkin(5));
t331 = cos(qJ(4));
t346 = t328 * t331;
t287 = -t311 * t329 + t312 * t346;
t347 = t328 * t329;
t288 = t311 * t331 + t312 * t347;
t289 = -t311 * t346 - t312 * t329;
t290 = -t311 * t347 + t312 * t331;
t348 = t312 * t327;
t349 = t311 * t327;
t337 = (Icges(5,5) * t288 + Icges(5,6) * t287 - Icges(5,3) * t348) * t312 - (Icges(5,5) * t290 + Icges(5,6) * t289 + Icges(5,3) * t349) * t311;
t354 = t337 * t327;
t324 = qJD(1) + qJD(2);
t352 = pkin(2) * t324;
t351 = t331 * pkin(4);
t350 = pkin(1) * qJD(1);
t325 = qJ(4) + qJ(5);
t332 = cos(qJ(1));
t314 = t332 * t350;
t322 = cos(t326);
t345 = t322 * t352 + t314;
t344 = qJD(4) * t327;
t343 = qJD(4) + qJD(5);
t330 = sin(qJ(1));
t342 = t330 * t350;
t316 = qJD(3) + t324;
t341 = t316 * (t312 * pkin(3) + pkin(9) * t349) + t345;
t340 = t311 * t344;
t339 = t312 * t344;
t338 = pkin(4) * t347 - pkin(10) * t327;
t303 = qJD(4) * t328 + t316;
t320 = sin(t326);
t336 = -t320 * t352 - t342;
t335 = -(t311 * pkin(3) - pkin(9) * t348) * t316 + t336;
t321 = cos(t325);
t319 = sin(t325);
t318 = pkin(5) - t325;
t317 = pkin(5) + t325;
t310 = cos(t317);
t309 = sin(t318);
t308 = cos(t318) / 0.2e1;
t307 = sin(t317) / 0.2e1;
t306 = t332 * rSges(2,1) - t330 * rSges(2,2);
t305 = t330 * rSges(2,1) + t332 * rSges(2,2);
t301 = t308 - t310 / 0.2e1;
t300 = t308 + t310 / 0.2e1;
t299 = t307 - t309 / 0.2e1;
t298 = t307 + t309 / 0.2e1;
t297 = qJD(5) * t328 + t303;
t295 = t327 * t329 * pkin(4) + pkin(10) * t328;
t294 = t328 * rSges(5,3) + (rSges(5,1) * t329 + rSges(5,2) * t331) * t327;
t293 = Icges(5,5) * t328 + (Icges(5,1) * t329 + Icges(5,4) * t331) * t327;
t292 = Icges(5,6) * t328 + (Icges(5,4) * t329 + Icges(5,2) * t331) * t327;
t291 = Icges(5,3) * t328 + (Icges(5,5) * t329 + Icges(5,6) * t331) * t327;
t286 = t343 * t348;
t285 = t343 * t349;
t284 = t314 + t324 * (rSges(3,1) * t322 - rSges(3,2) * t320);
t283 = -t342 - t324 * (rSges(3,1) * t320 + rSges(3,2) * t322);
t281 = -t299 * t311 + t312 * t321;
t280 = -t300 * t311 - t312 * t319;
t279 = t299 * t312 + t311 * t321;
t278 = t300 * t312 - t311 * t319;
t277 = t316 * (rSges(4,1) * t312 - rSges(4,2) * t311) + t345;
t276 = -t316 * (rSges(4,1) * t311 + rSges(4,2) * t312) + t336;
t275 = rSges(6,1) * t301 + rSges(6,2) * t298 + rSges(6,3) * t328;
t274 = Icges(6,1) * t301 + Icges(6,4) * t298 + Icges(6,5) * t328;
t273 = Icges(6,4) * t301 + Icges(6,2) * t298 + Icges(6,6) * t328;
t272 = Icges(6,5) * t301 + Icges(6,6) * t298 + Icges(6,3) * t328;
t271 = -t338 * t311 + t351 * t312;
t270 = t351 * t311 + t338 * t312;
t269 = rSges(5,1) * t290 + rSges(5,2) * t289 + rSges(5,3) * t349;
t268 = rSges(5,1) * t288 + rSges(5,2) * t287 - rSges(5,3) * t348;
t267 = Icges(5,1) * t290 + Icges(5,4) * t289 + Icges(5,5) * t349;
t266 = Icges(5,1) * t288 + Icges(5,4) * t287 - Icges(5,5) * t348;
t265 = Icges(5,4) * t290 + Icges(5,2) * t289 + Icges(5,6) * t349;
t264 = Icges(5,4) * t288 + Icges(5,2) * t287 - Icges(5,6) * t348;
t261 = rSges(6,1) * t281 + rSges(6,2) * t280 + rSges(6,3) * t349;
t260 = rSges(6,1) * t279 + rSges(6,2) * t278 - rSges(6,3) * t348;
t259 = Icges(6,1) * t281 + Icges(6,4) * t280 + Icges(6,5) * t349;
t258 = Icges(6,1) * t279 + Icges(6,4) * t278 - Icges(6,5) * t348;
t257 = Icges(6,4) * t281 + Icges(6,2) * t280 + Icges(6,6) * t349;
t256 = Icges(6,4) * t279 + Icges(6,2) * t278 - Icges(6,6) * t348;
t255 = Icges(6,5) * t281 + Icges(6,6) * t280 + Icges(6,3) * t349;
t254 = Icges(6,5) * t279 + Icges(6,6) * t278 - Icges(6,3) * t348;
t253 = (t268 * t311 + t269 * t312) * t344;
t252 = t269 * t303 - t294 * t340 + t341;
t251 = -t268 * t303 - t294 * t339 + t335;
t250 = t261 * t297 + t271 * t303 - t275 * t285 - t295 * t340 + t341;
t249 = -t260 * t297 - t270 * t303 - t275 * t286 - t295 * t339 + t335;
t248 = t285 * t260 + t286 * t261 + (t270 * t311 + t271 * t312) * t344;
t1 = m(3) * (t283 ^ 2 + t284 ^ 2) / 0.2e1 + t324 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t276 ^ 2 + t277 ^ 2) / 0.2e1 + t316 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + ((t289 * t292 + t290 * t293 + t291 * t349) * t303 + (-(t264 * t289 + t266 * t290) * t312 + (t289 * t265 + t290 * t267 - t354) * t311) * t344) * t340 / 0.2e1 - ((t287 * t292 + t288 * t293 - t291 * t348) * t303 + ((t265 * t287 + t267 * t288) * t311 + (-t287 * t264 - t288 * t266 + t354) * t312) * t344) * t339 / 0.2e1 + t303 * ((t328 * t291 + (t292 * t331 + t293 * t329) * t327) * t303 + (((t265 * t331 + t267 * t329) * t311 - (t264 * t331 + t266 * t329) * t312) * t327 - t337 * t328) * t344) / 0.2e1 + m(6) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + t285 * ((t255 * t349 + t280 * t257 + t281 * t259) * t285 - (t254 * t349 + t256 * t280 + t258 * t281) * t286 + (t272 * t349 + t273 * t280 + t274 * t281) * t297) / 0.2e1 - t286 * ((-t255 * t348 + t257 * t278 + t259 * t279) * t285 - (-t254 * t348 + t278 * t256 + t279 * t258) * t286 + (-t272 * t348 + t273 * t278 + t274 * t279) * t297) / 0.2e1 + t297 * ((t255 * t328 + t257 * t298 + t259 * t301) * t285 - (t254 * t328 + t256 * t298 + t258 * t301) * t286 + (t328 * t272 + t298 * t273 + t301 * t274) * t297) / 0.2e1 + (m(2) * (t305 ^ 2 + t306 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
