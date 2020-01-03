% Calculate kinetic energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:21
% EndTime: 2019-12-31 17:29:22
% DurationCPUTime: 1.16s
% Computational Cost: add. (1101->224), mult. (2769->366), div. (0->0), fcn. (3324->10), ass. (0->102)
t324 = sin(qJ(2));
t325 = sin(qJ(1));
t327 = cos(qJ(2));
t328 = cos(qJ(1));
t343 = cos(pkin(4));
t333 = t328 * t343;
t305 = t324 * t325 - t327 * t333;
t306 = t324 * t333 + t325 * t327;
t321 = sin(pkin(4));
t340 = t321 * t328;
t274 = Icges(3,5) * t306 - Icges(3,6) * t305 - Icges(3,3) * t340;
t334 = t325 * t343;
t307 = t328 * t324 + t327 * t334;
t308 = -t324 * t334 + t328 * t327;
t342 = t321 * t325;
t275 = Icges(3,5) * t308 - Icges(3,6) * t307 + Icges(3,3) * t342;
t346 = (t274 * t328 - t275 * t325) * t321;
t344 = cos(qJ(3));
t341 = t321 * t327;
t286 = pkin(2) * t306 + pkin(7) * t305;
t287 = pkin(2) * t308 + pkin(7) * t307;
t337 = qJD(2) * t321;
t317 = t325 * t337;
t335 = t328 * t337;
t339 = t286 * t317 + t287 * t335;
t295 = qJD(3) * t307 + t317;
t338 = qJD(1) * (pkin(1) * t325 - pkin(6) * t340);
t318 = qJD(2) * t343 + qJD(1);
t336 = t321 * t344;
t296 = qJD(3) * t305 - t335;
t310 = -qJD(3) * t341 + t318;
t309 = (pkin(2) * t324 - pkin(7) * t327) * t321;
t311 = qJD(1) * (pkin(1) * t328 + pkin(6) * t342);
t331 = t318 * t287 - t309 * t317 + t311;
t330 = -t286 * t318 - t309 * t335 - t338;
t326 = cos(qJ(4));
t323 = sin(qJ(3));
t322 = sin(qJ(4));
t314 = rSges(2,1) * t328 - rSges(2,2) * t325;
t313 = rSges(2,1) * t325 + rSges(2,2) * t328;
t304 = t323 * t343 + t324 * t336;
t303 = t321 * t323 * t324 - t343 * t344;
t300 = t343 * rSges(3,3) + (rSges(3,1) * t324 + rSges(3,2) * t327) * t321;
t299 = Icges(3,5) * t343 + (Icges(3,1) * t324 + Icges(3,4) * t327) * t321;
t298 = Icges(3,6) * t343 + (Icges(3,4) * t324 + Icges(3,2) * t327) * t321;
t297 = Icges(3,3) * t343 + (Icges(3,5) * t324 + Icges(3,6) * t327) * t321;
t294 = t308 * t344 + t323 * t342;
t293 = t308 * t323 - t325 * t336;
t292 = t306 * t344 - t323 * t340;
t291 = t306 * t323 + t328 * t336;
t290 = t304 * t326 - t322 * t341;
t289 = -t304 * t322 - t326 * t341;
t288 = qJD(4) * t303 + t310;
t285 = pkin(3) * t304 + pkin(8) * t303;
t282 = rSges(3,1) * t308 - rSges(3,2) * t307 + rSges(3,3) * t342;
t281 = rSges(3,1) * t306 - rSges(3,2) * t305 - rSges(3,3) * t340;
t279 = Icges(3,1) * t308 - Icges(3,4) * t307 + Icges(3,5) * t342;
t278 = Icges(3,1) * t306 - Icges(3,4) * t305 - Icges(3,5) * t340;
t277 = Icges(3,4) * t308 - Icges(3,2) * t307 + Icges(3,6) * t342;
t276 = Icges(3,4) * t306 - Icges(3,2) * t305 - Icges(3,6) * t340;
t273 = rSges(4,1) * t304 - rSges(4,2) * t303 - rSges(4,3) * t341;
t272 = Icges(4,1) * t304 - Icges(4,4) * t303 - Icges(4,5) * t341;
t271 = Icges(4,4) * t304 - Icges(4,2) * t303 - Icges(4,6) * t341;
t270 = Icges(4,5) * t304 - Icges(4,6) * t303 - Icges(4,3) * t341;
t269 = t294 * t326 + t307 * t322;
t268 = -t294 * t322 + t307 * t326;
t267 = t292 * t326 + t305 * t322;
t266 = -t292 * t322 + t305 * t326;
t265 = qJD(4) * t291 + t296;
t264 = qJD(4) * t293 + t295;
t263 = pkin(3) * t294 + pkin(8) * t293;
t262 = pkin(3) * t292 + pkin(8) * t291;
t261 = rSges(4,1) * t294 - rSges(4,2) * t293 + rSges(4,3) * t307;
t260 = rSges(4,1) * t292 - rSges(4,2) * t291 + rSges(4,3) * t305;
t259 = Icges(4,1) * t294 - Icges(4,4) * t293 + Icges(4,5) * t307;
t258 = Icges(4,1) * t292 - Icges(4,4) * t291 + Icges(4,5) * t305;
t257 = Icges(4,4) * t294 - Icges(4,2) * t293 + Icges(4,6) * t307;
t256 = Icges(4,4) * t292 - Icges(4,2) * t291 + Icges(4,6) * t305;
t255 = Icges(4,5) * t294 - Icges(4,6) * t293 + Icges(4,3) * t307;
t254 = Icges(4,5) * t292 - Icges(4,6) * t291 + Icges(4,3) * t305;
t253 = rSges(5,1) * t290 + rSges(5,2) * t289 + rSges(5,3) * t303;
t252 = Icges(5,1) * t290 + Icges(5,4) * t289 + Icges(5,5) * t303;
t251 = Icges(5,4) * t290 + Icges(5,2) * t289 + Icges(5,6) * t303;
t250 = Icges(5,5) * t290 + Icges(5,6) * t289 + Icges(5,3) * t303;
t249 = t282 * t318 - t300 * t317 + t311;
t248 = -t281 * t318 - t300 * t335 - t338;
t247 = (t281 * t325 + t282 * t328) * t337;
t246 = rSges(5,1) * t269 + rSges(5,2) * t268 + rSges(5,3) * t293;
t245 = rSges(5,1) * t267 + rSges(5,2) * t266 + rSges(5,3) * t291;
t244 = Icges(5,1) * t269 + Icges(5,4) * t268 + Icges(5,5) * t293;
t243 = Icges(5,1) * t267 + Icges(5,4) * t266 + Icges(5,5) * t291;
t242 = Icges(5,4) * t269 + Icges(5,2) * t268 + Icges(5,6) * t293;
t241 = Icges(5,4) * t267 + Icges(5,2) * t266 + Icges(5,6) * t291;
t240 = Icges(5,5) * t269 + Icges(5,6) * t268 + Icges(5,3) * t293;
t239 = Icges(5,5) * t267 + Icges(5,6) * t266 + Icges(5,3) * t291;
t238 = t261 * t310 - t273 * t295 + t331;
t237 = -t260 * t310 + t273 * t296 + t330;
t236 = t260 * t295 - t261 * t296 + t339;
t235 = t246 * t288 - t253 * t264 + t263 * t310 - t285 * t295 + t331;
t234 = -t245 * t288 + t253 * t265 - t262 * t310 + t285 * t296 + t330;
t233 = t245 * t264 - t246 * t265 + t262 * t295 - t263 * t296 + t339;
t1 = m(3) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + ((t297 * t342 - t298 * t307 + t299 * t308) * t318 + (-(-t276 * t307 + t278 * t308) * t328 + (-t307 * t277 + t308 * t279 - t346) * t325) * t337) * t317 / 0.2e1 - ((-t297 * t340 - t298 * t305 + t299 * t306) * t318 + ((-t277 * t305 + t279 * t306) * t325 + (t305 * t276 - t306 * t278 + t346) * t328) * t337) * t335 / 0.2e1 + t318 * ((t343 * t275 + (t277 * t327 + t279 * t324) * t321) * t317 - (t343 * t274 + (t276 * t327 + t278 * t324) * t321) * t335 + (t343 * t297 + (t298 * t327 + t299 * t324) * t321) * t318) / 0.2e1 + m(4) * (t236 ^ 2 + t237 ^ 2 + t238 ^ 2) / 0.2e1 + t295 * ((t255 * t307 - t257 * t293 + t259 * t294) * t295 + (t254 * t307 - t256 * t293 + t258 * t294) * t296 + (t270 * t307 - t271 * t293 + t272 * t294) * t310) / 0.2e1 + t296 * ((t255 * t305 - t257 * t291 + t259 * t292) * t295 + (t254 * t305 - t256 * t291 + t258 * t292) * t296 + (t270 * t305 - t271 * t291 + t272 * t292) * t310) / 0.2e1 + t310 * ((-t255 * t341 - t257 * t303 + t259 * t304) * t295 + (-t254 * t341 - t256 * t303 + t258 * t304) * t296 + (-t270 * t341 - t271 * t303 + t272 * t304) * t310) / 0.2e1 + m(5) * (t233 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + t264 * ((t240 * t293 + t242 * t268 + t244 * t269) * t264 + (t239 * t293 + t241 * t268 + t243 * t269) * t265 + (t250 * t293 + t251 * t268 + t252 * t269) * t288) / 0.2e1 + t265 * ((t240 * t291 + t242 * t266 + t244 * t267) * t264 + (t239 * t291 + t241 * t266 + t243 * t267) * t265 + (t250 * t291 + t251 * t266 + t252 * t267) * t288) / 0.2e1 + t288 * ((t240 * t303 + t242 * t289 + t244 * t290) * t264 + (t239 * t303 + t241 * t289 + t243 * t290) * t265 + (t250 * t303 + t251 * t289 + t252 * t290) * t288) / 0.2e1 + (m(2) * (t313 ^ 2 + t314 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
