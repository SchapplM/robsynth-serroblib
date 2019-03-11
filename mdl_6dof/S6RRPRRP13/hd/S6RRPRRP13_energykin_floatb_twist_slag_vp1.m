% Calculate kinetic energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP13_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP13_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:29
% EndTime: 2019-03-09 12:55:33
% DurationCPUTime: 3.95s
% Computational Cost: add. (2361->345), mult. (5348->482), div. (0->0), fcn. (6303->10), ass. (0->160)
t371 = Icges(3,1) + Icges(4,2);
t370 = Icges(6,1) + Icges(7,1);
t369 = Icges(3,4) + Icges(4,6);
t368 = Icges(6,4) + Icges(7,4);
t367 = Icges(3,5) - Icges(4,4);
t366 = Icges(6,5) + Icges(7,5);
t365 = Icges(3,2) + Icges(4,3);
t364 = Icges(6,2) + Icges(7,2);
t363 = Icges(3,6) - Icges(4,5);
t362 = Icges(6,6) + Icges(7,6);
t361 = Icges(3,3) + Icges(4,1);
t360 = Icges(6,3) + Icges(7,3);
t359 = rSges(7,3) + qJ(6);
t290 = cos(pkin(6));
t295 = sin(qJ(1));
t297 = cos(qJ(2));
t321 = t295 * t297;
t294 = sin(qJ(2));
t298 = cos(qJ(1));
t322 = t294 * t298;
t258 = t290 * t321 + t322;
t293 = sin(qJ(4));
t289 = sin(pkin(6));
t334 = cos(qJ(4));
t311 = t289 * t334;
t229 = t258 * t293 + t295 * t311;
t320 = t297 * t298;
t323 = t294 * t295;
t259 = -t290 * t323 + t320;
t292 = sin(qJ(5));
t296 = cos(qJ(5));
t192 = -t229 * t292 + t259 * t296;
t328 = t259 * t292;
t193 = t229 * t296 + t328;
t326 = t289 * t295;
t228 = -t258 * t334 + t293 * t326;
t358 = t362 * t192 + t366 * t193 + t360 * t228;
t256 = -t290 * t320 + t323;
t231 = t256 * t293 - t298 * t311;
t257 = t290 * t322 + t321;
t194 = -t231 * t292 + t257 * t296;
t329 = t257 * t292;
t195 = t231 * t296 + t329;
t324 = t289 * t298;
t230 = t256 * t334 + t293 * t324;
t357 = t362 * t194 + t366 * t195 - t360 * t230;
t356 = t364 * t192 + t368 * t193 + t362 * t228;
t355 = t364 * t194 + t368 * t195 - t362 * t230;
t354 = t368 * t192 + t370 * t193 + t366 * t228;
t353 = t368 * t194 + t370 * t195 - t366 * t230;
t325 = t289 * t297;
t255 = t290 * t334 - t293 * t325;
t327 = t289 * t294;
t224 = -t255 * t292 + t296 * t327;
t312 = t292 * t327;
t225 = t255 * t296 + t312;
t254 = t290 * t293 + t297 * t311;
t352 = t362 * t224 + t366 * t225 + t360 * t254;
t351 = t364 * t224 + t368 * t225 + t362 * t254;
t350 = t368 * t224 + t370 * t225 + t366 * t254;
t349 = t365 * t258 - t369 * t259 - t363 * t326;
t348 = t365 * t256 - t369 * t257 + t363 * t324;
t347 = -t369 * t258 + t371 * t259 + t367 * t326;
t346 = -t369 * t256 + t371 * t257 - t367 * t324;
t345 = -t363 * t258 + t367 * t259 + t361 * t326;
t344 = -t363 * t256 + t367 * t257 - t361 * t324;
t343 = t361 * t290 + (t367 * t294 + t363 * t297) * t289;
t342 = t363 * t290 + (t369 * t294 + t365 * t297) * t289;
t341 = t367 * t290 + (t371 * t294 + t369 * t297) * t289;
t333 = pkin(8) * t290;
t332 = pkin(5) * t296;
t330 = Icges(2,4) * t295;
t319 = rSges(7,1) * t193 + rSges(7,2) * t192 + pkin(5) * t328 + t228 * t359 + t229 * t332;
t318 = rSges(7,1) * t195 + rSges(7,2) * t194 + pkin(5) * t329 - t230 * t359 + t231 * t332;
t317 = rSges(7,1) * t225 + rSges(7,2) * t224 + pkin(5) * t312 + t254 * t359 + t255 * t332;
t316 = qJD(2) * t289;
t315 = V_base(5) * pkin(7) + V_base(1);
t269 = t295 * t316 + V_base(4);
t286 = V_base(6) + qJD(1);
t227 = qJD(4) * t259 + t269;
t270 = qJD(2) * t290 + t286;
t252 = qJD(4) * t327 + t270;
t268 = -t298 * t316 + V_base(5);
t263 = t295 * pkin(1) - pkin(8) * t324;
t310 = -t263 * t286 + V_base(5) * t333 + t315;
t264 = pkin(1) * t298 + pkin(8) * t326;
t309 = V_base(4) * t263 - t264 * V_base(5) + V_base(3);
t226 = qJD(4) * t257 + t268;
t260 = (pkin(2) * t294 - qJ(3) * t297) * t289;
t308 = qJD(3) * t258 + t268 * t260 + t310;
t307 = t286 * t264 + V_base(2) + (-pkin(7) - t333) * V_base(4);
t223 = pkin(2) * t259 + qJ(3) * t258;
t306 = qJD(3) * t256 + t270 * t223 + t307;
t222 = pkin(2) * t257 + qJ(3) * t256;
t305 = -qJD(3) * t325 + t269 * t222 + t309;
t237 = -pkin(3) * t324 + t257 * pkin(9);
t262 = pkin(3) * t290 + pkin(9) * t327;
t304 = t268 * t262 + (-t222 - t237) * t270 + t308;
t236 = pkin(3) * t326 + pkin(9) * t259;
t303 = t270 * t236 + (-t260 - t262) * t269 + t306;
t302 = t269 * t237 + (-t223 - t236) * t268 + t305;
t189 = pkin(4) * t231 - pkin(10) * t230;
t219 = pkin(4) * t255 + pkin(10) * t254;
t301 = -t189 * t252 + t226 * t219 + t304;
t188 = pkin(4) * t229 + pkin(10) * t228;
t300 = t252 * t188 - t219 * t227 + t303;
t299 = -t188 * t226 + t227 * t189 + t302;
t287 = Icges(2,4) * t298;
t278 = rSges(2,1) * t298 - t295 * rSges(2,2);
t277 = t295 * rSges(2,1) + rSges(2,2) * t298;
t276 = Icges(2,1) * t298 - t330;
t275 = Icges(2,1) * t295 + t287;
t274 = -Icges(2,2) * t295 + t287;
t273 = Icges(2,2) * t298 + t330;
t267 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t266 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t265 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t246 = rSges(4,1) * t290 + (-rSges(4,2) * t294 - rSges(4,3) * t297) * t289;
t245 = rSges(3,3) * t290 + (rSges(3,1) * t294 + rSges(3,2) * t297) * t289;
t235 = V_base(5) * rSges(2,3) - t277 * t286 + t315;
t234 = t278 * t286 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t232 = t277 * V_base(4) - t278 * V_base(5) + V_base(3);
t220 = qJD(5) * t254 + t252;
t217 = rSges(3,1) * t259 - rSges(3,2) * t258 + rSges(3,3) * t326;
t216 = t257 * rSges(3,1) - t256 * rSges(3,2) - rSges(3,3) * t324;
t215 = -rSges(4,1) * t324 - t257 * rSges(4,2) + t256 * rSges(4,3);
t214 = rSges(4,1) * t326 - rSges(4,2) * t259 + rSges(4,3) * t258;
t201 = rSges(5,1) * t255 - rSges(5,2) * t254 + rSges(5,3) * t327;
t200 = Icges(5,1) * t255 - Icges(5,4) * t254 + Icges(5,5) * t327;
t199 = Icges(5,4) * t255 - Icges(5,2) * t254 + Icges(5,6) * t327;
t198 = Icges(5,5) * t255 - Icges(5,6) * t254 + Icges(5,3) * t327;
t191 = qJD(5) * t228 + t227;
t190 = -qJD(5) * t230 + t226;
t186 = rSges(5,1) * t231 + rSges(5,2) * t230 + rSges(5,3) * t257;
t185 = rSges(5,1) * t229 - rSges(5,2) * t228 + rSges(5,3) * t259;
t183 = Icges(5,1) * t231 + Icges(5,4) * t230 + Icges(5,5) * t257;
t182 = Icges(5,1) * t229 - Icges(5,4) * t228 + Icges(5,5) * t259;
t181 = Icges(5,4) * t231 + Icges(5,2) * t230 + Icges(5,6) * t257;
t180 = Icges(5,4) * t229 - Icges(5,2) * t228 + Icges(5,6) * t259;
t179 = Icges(5,5) * t231 + Icges(5,6) * t230 + Icges(5,3) * t257;
t178 = Icges(5,5) * t229 - Icges(5,6) * t228 + Icges(5,3) * t259;
t177 = rSges(6,1) * t225 + rSges(6,2) * t224 + rSges(6,3) * t254;
t167 = -t216 * t270 + t245 * t268 + t310;
t166 = t217 * t270 - t245 * t269 + t307;
t165 = rSges(6,1) * t195 + rSges(6,2) * t194 - rSges(6,3) * t230;
t163 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t228;
t149 = t216 * t269 - t217 * t268 + t309;
t146 = t246 * t268 + (-t215 - t222) * t270 + t308;
t145 = t214 * t270 + (-t246 - t260) * t269 + t306;
t144 = t215 * t269 + (-t214 - t223) * t268 + t305;
t143 = -t186 * t252 + t201 * t226 + t304;
t142 = t185 * t252 - t201 * t227 + t303;
t141 = -t185 * t226 + t186 * t227 + t302;
t140 = -t165 * t220 + t177 * t190 + t301;
t139 = t163 * t220 - t177 * t191 + t300;
t138 = -t163 * t190 + t165 * t191 + t299;
t137 = qJD(6) * t228 + t190 * t317 - t220 * t318 + t301;
t136 = -qJD(6) * t230 - t191 * t317 + t220 * t319 + t300;
t135 = qJD(6) * t254 - t190 * t319 + t191 * t318 + t299;
t1 = m(2) * (t232 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + m(1) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t227 * ((t259 * t178 - t228 * t180 + t229 * t182) * t227 + (t179 * t259 - t181 * t228 + t183 * t229) * t226 + (t198 * t259 - t199 * t228 + t200 * t229) * t252) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t149 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t252 * ((t178 * t327 - t180 * t254 + t182 * t255) * t227 + (t179 * t327 - t181 * t254 + t183 * t255) * t226 + (t198 * t327 - t254 * t199 + t255 * t200) * t252) / 0.2e1 + t226 * ((t178 * t257 + t180 * t230 + t182 * t231) * t227 + (t257 * t179 + t230 * t181 + t231 * t183) * t226 + (t198 * t257 + t199 * t230 + t200 * t231) * t252) / 0.2e1 + ((t194 * t351 + t195 * t350 - t230 * t352) * t220 + (t356 * t194 + t354 * t195 - t230 * t358) * t191 + (t355 * t194 + t353 * t195 - t357 * t230) * t190) * t190 / 0.2e1 + ((t192 * t351 + t193 * t350 + t228 * t352) * t220 + (t356 * t192 + t354 * t193 + t358 * t228) * t191 + (t192 * t355 + t193 * t353 + t228 * t357) * t190) * t191 / 0.2e1 + ((t351 * t224 + t350 * t225 + t352 * t254) * t220 + (t356 * t224 + t354 * t225 + t254 * t358) * t191 + (t224 * t355 + t225 * t353 + t254 * t357) * t190) * t220 / 0.2e1 + ((-t256 * t342 + t257 * t341 - t324 * t343) * t270 + (t256 * t349 + t257 * t347 - t324 * t345) * t269 + (t348 * t256 + t346 * t257 - t344 * t324) * t268) * t268 / 0.2e1 + ((-t258 * t342 + t259 * t341 + t326 * t343) * t270 + (t349 * t258 + t347 * t259 + t345 * t326) * t269 + (t258 * t348 + t259 * t346 + t326 * t344) * t268) * t269 / 0.2e1 + ((t268 * t344 + t269 * t345 + t343 * t270) * t290 + ((t294 * t341 + t297 * t342) * t270 + (t294 * t347 - t297 * t349) * t269 + (t294 * t346 - t297 * t348) * t268) * t289) * t270 / 0.2e1 + ((-t295 * t273 + t275 * t298 + Icges(1,4)) * V_base(5) + (-t295 * t274 + t276 * t298 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t273 * t298 + t295 * t275 + Icges(1,2)) * V_base(5) + (t274 * t298 + t295 * t276 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t295 + Icges(2,6) * t298) * V_base(5) + (Icges(2,5) * t298 - Icges(2,6) * t295) * V_base(4) + Icges(2,3) * t286 / 0.2e1) * t286;
T  = t1;
