% Calculate kinetic energy for
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:28
% EndTime: 2019-03-09 17:41:32
% DurationCPUTime: 3.85s
% Computational Cost: add. (2567->345), mult. (5854->481), div. (0->0), fcn. (6982->10), ass. (0->160)
t370 = Icges(4,1) + Icges(5,2);
t369 = Icges(5,1) + Icges(4,3);
t368 = Icges(6,1) + Icges(7,1);
t367 = -Icges(4,4) - Icges(5,6);
t366 = Icges(5,4) - Icges(4,5);
t365 = Icges(6,4) + Icges(7,4);
t364 = Icges(5,5) - Icges(4,6);
t363 = Icges(6,5) + Icges(7,5);
t362 = Icges(4,2) + Icges(5,3);
t361 = Icges(6,2) + Icges(7,2);
t360 = Icges(6,6) + Icges(7,6);
t359 = Icges(6,3) + Icges(7,3);
t358 = rSges(7,3) + qJ(6);
t289 = cos(pkin(6));
t293 = sin(qJ(1));
t295 = cos(qJ(2));
t319 = t293 * t295;
t292 = sin(qJ(2));
t296 = cos(qJ(1));
t320 = t292 * t296;
t256 = t289 * t320 + t319;
t288 = sin(pkin(6));
t333 = cos(qJ(3));
t310 = t288 * t333;
t332 = sin(qJ(3));
t235 = t256 * t332 + t296 * t310;
t318 = t295 * t296;
t321 = t292 * t293;
t255 = -t289 * t318 + t321;
t291 = sin(qJ(5));
t294 = cos(qJ(5));
t201 = t235 * t294 - t255 * t291;
t327 = t235 * t291;
t202 = t255 * t294 + t327;
t309 = t288 * t332;
t236 = t256 * t333 - t296 * t309;
t357 = t360 * t201 + t363 * t202 + t359 * t236;
t258 = -t289 * t321 + t318;
t237 = t258 * t332 - t293 * t310;
t257 = t289 * t319 + t320;
t203 = t237 * t294 - t257 * t291;
t326 = t237 * t291;
t204 = t257 * t294 + t326;
t238 = t258 * t333 + t293 * t309;
t356 = t360 * t203 + t363 * t204 + t359 * t238;
t355 = t361 * t201 + t365 * t202 + t360 * t236;
t354 = t361 * t203 + t365 * t204 + t360 * t238;
t353 = t365 * t201 + t368 * t202 + t363 * t236;
t352 = t365 * t203 + t368 * t204 + t363 * t238;
t253 = -t289 * t333 + t292 * t309;
t323 = t288 * t295;
t231 = t253 * t294 + t291 * t323;
t325 = t253 * t291;
t232 = -t294 * t323 + t325;
t254 = t289 * t332 + t292 * t310;
t351 = t360 * t231 + t363 * t232 + t359 * t254;
t350 = t361 * t231 + t365 * t232 + t360 * t254;
t349 = t365 * t231 + t368 * t232 + t363 * t254;
t348 = t362 * t235 + t367 * t236 + t364 * t255;
t347 = t362 * t237 + t367 * t238 + t364 * t257;
t346 = t364 * t235 - t366 * t236 + t369 * t255;
t345 = t364 * t237 - t366 * t238 + t369 * t257;
t344 = t367 * t235 + t370 * t236 - t366 * t255;
t343 = t367 * t237 + t370 * t238 - t366 * t257;
t342 = t362 * t253 + t367 * t254 - t364 * t323;
t341 = t367 * t253 + t370 * t254 + t366 * t323;
t340 = t364 * t253 - t366 * t254 - t369 * t323;
t331 = pkin(8) * t289;
t330 = pkin(5) * t294;
t328 = Icges(2,4) * t293;
t324 = t288 * t293;
t322 = t288 * t296;
t317 = rSges(7,1) * t202 + rSges(7,2) * t201 + pkin(5) * t327 + t236 * t358 + t255 * t330;
t316 = rSges(7,1) * t204 + rSges(7,2) * t203 + pkin(5) * t326 + t238 * t358 + t257 * t330;
t315 = rSges(7,1) * t232 + rSges(7,2) * t231 + pkin(5) * t325 + t254 * t358 - t323 * t330;
t314 = qJD(2) * t288;
t313 = V_base(5) * pkin(7) + V_base(1);
t267 = t293 * t314 + V_base(4);
t285 = V_base(6) + qJD(1);
t234 = qJD(3) * t257 + t267;
t268 = qJD(2) * t289 + t285;
t266 = -t296 * t314 + V_base(5);
t261 = t293 * pkin(1) - pkin(8) * t322;
t308 = -t261 * t285 + V_base(5) * t331 + t313;
t262 = pkin(1) * t296 + pkin(8) * t324;
t307 = V_base(4) * t261 - t262 * V_base(5) + V_base(3);
t233 = qJD(3) * t255 + t266;
t251 = -qJD(3) * t323 + t268;
t306 = t285 * t262 + V_base(2) + (-pkin(7) - t331) * V_base(4);
t227 = pkin(2) * t256 + pkin(9) * t255;
t260 = (pkin(2) * t292 - pkin(9) * t295) * t288;
t305 = -t227 * t268 + t266 * t260 + t308;
t228 = pkin(2) * t258 + pkin(9) * t257;
t304 = t267 * t227 - t228 * t266 + t307;
t225 = pkin(3) * t254 + qJ(4) * t253;
t303 = qJD(4) * t237 + t233 * t225 + t305;
t196 = pkin(3) * t236 + qJ(4) * t235;
t302 = qJD(4) * t253 + t234 * t196 + t304;
t301 = t268 * t228 - t260 * t267 + t306;
t197 = pkin(3) * t238 + qJ(4) * t237;
t300 = qJD(4) * t235 + t251 * t197 + t301;
t205 = pkin(4) * t255 + pkin(10) * t236;
t243 = -pkin(4) * t323 + pkin(10) * t254;
t299 = t233 * t243 + (-t196 - t205) * t251 + t303;
t206 = pkin(4) * t257 + pkin(10) * t238;
t298 = t234 * t205 + (-t197 - t206) * t233 + t302;
t297 = t251 * t206 + (-t225 - t243) * t234 + t300;
t286 = Icges(2,4) * t296;
t276 = rSges(2,1) * t296 - t293 * rSges(2,2);
t275 = t293 * rSges(2,1) + rSges(2,2) * t296;
t274 = Icges(2,1) * t296 - t328;
t273 = Icges(2,1) * t293 + t286;
t272 = -Icges(2,2) * t293 + t286;
t271 = Icges(2,2) * t296 + t328;
t265 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t264 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t263 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t247 = rSges(3,3) * t289 + (rSges(3,1) * t292 + rSges(3,2) * t295) * t288;
t246 = Icges(3,5) * t289 + (Icges(3,1) * t292 + Icges(3,4) * t295) * t288;
t245 = Icges(3,6) * t289 + (Icges(3,4) * t292 + Icges(3,2) * t295) * t288;
t244 = Icges(3,3) * t289 + (Icges(3,5) * t292 + Icges(3,6) * t295) * t288;
t242 = V_base(5) * rSges(2,3) - t275 * t285 + t313;
t241 = t276 * t285 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t239 = t275 * V_base(4) - t276 * V_base(5) + V_base(3);
t226 = qJD(5) * t254 + t251;
t224 = rSges(3,1) * t258 - rSges(3,2) * t257 + rSges(3,3) * t324;
t223 = t256 * rSges(3,1) - t255 * rSges(3,2) - rSges(3,3) * t322;
t222 = Icges(3,1) * t258 - Icges(3,4) * t257 + Icges(3,5) * t324;
t221 = Icges(3,1) * t256 - Icges(3,4) * t255 - Icges(3,5) * t322;
t220 = Icges(3,4) * t258 - Icges(3,2) * t257 + Icges(3,6) * t324;
t219 = Icges(3,4) * t256 - Icges(3,2) * t255 - Icges(3,6) * t322;
t218 = Icges(3,5) * t258 - Icges(3,6) * t257 + Icges(3,3) * t324;
t217 = Icges(3,5) * t256 - Icges(3,6) * t255 - Icges(3,3) * t322;
t216 = rSges(4,1) * t254 - rSges(4,2) * t253 - rSges(4,3) * t323;
t215 = -rSges(5,1) * t323 - rSges(5,2) * t254 + rSges(5,3) * t253;
t199 = qJD(5) * t238 + t234;
t198 = qJD(5) * t236 + t233;
t191 = rSges(4,1) * t238 - rSges(4,2) * t237 + rSges(4,3) * t257;
t190 = rSges(4,1) * t236 - rSges(4,2) * t235 + rSges(4,3) * t255;
t189 = rSges(5,1) * t257 - rSges(5,2) * t238 + rSges(5,3) * t237;
t188 = rSges(5,1) * t255 - rSges(5,2) * t236 + rSges(5,3) * t235;
t174 = rSges(6,1) * t232 + rSges(6,2) * t231 + rSges(6,3) * t254;
t165 = -t223 * t268 + t247 * t266 + t308;
t164 = t224 * t268 - t247 * t267 + t306;
t161 = rSges(6,1) * t204 + rSges(6,2) * t203 + rSges(6,3) * t238;
t159 = rSges(6,1) * t202 + rSges(6,2) * t201 + rSges(6,3) * t236;
t145 = t223 * t267 - t224 * t266 + t307;
t144 = -t190 * t251 + t216 * t233 + t305;
t143 = t191 * t251 - t216 * t234 + t301;
t142 = t190 * t234 - t191 * t233 + t304;
t141 = t215 * t233 + (-t188 - t196) * t251 + t303;
t140 = t189 * t251 + (-t215 - t225) * t234 + t300;
t139 = t188 * t234 + (-t189 - t197) * t233 + t302;
t138 = -t159 * t226 + t174 * t198 + t299;
t137 = t161 * t226 - t174 * t199 + t297;
t136 = t159 * t199 - t161 * t198 + t298;
t135 = qJD(6) * t238 + t198 * t315 - t226 * t317 + t299;
t134 = qJD(6) * t236 - t199 * t315 + t226 * t316 + t297;
t133 = qJD(6) * t254 - t198 * t316 + t199 * t317 + t298;
t1 = m(3) * (t145 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t268 * ((t217 * t266 + t218 * t267 + t244 * t268) * t289 + ((t220 * t295 + t222 * t292) * t267 + (t219 * t295 + t221 * t292) * t266 + (t245 * t295 + t246 * t292) * t268) * t288) / 0.2e1 + m(2) * (t239 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + m(1) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t266 * ((-t218 * t322 - t255 * t220 + t256 * t222) * t267 + (-t217 * t322 - t255 * t219 + t256 * t221) * t266 + (-t244 * t322 - t255 * t245 + t256 * t246) * t268) / 0.2e1 + t267 * ((t218 * t324 - t257 * t220 + t258 * t222) * t267 + (t217 * t324 - t219 * t257 + t221 * t258) * t266 + (t244 * t324 - t245 * t257 + t246 * t258) * t268) / 0.2e1 + ((t201 * t350 + t202 * t349 + t236 * t351) * t226 + (t201 * t354 + t202 * t352 + t236 * t356) * t199 + (t355 * t201 + t353 * t202 + t357 * t236) * t198) * t198 / 0.2e1 + ((t203 * t350 + t204 * t349 + t238 * t351) * t226 + (t354 * t203 + t352 * t204 + t356 * t238) * t199 + (t355 * t203 + t353 * t204 + t238 * t357) * t198) * t199 / 0.2e1 + ((t231 * t350 + t232 * t349 + t254 * t351) * t226 + (t231 * t354 + t232 * t352 + t254 * t356) * t199 + (t355 * t231 + t353 * t232 + t254 * t357) * t198) * t226 / 0.2e1 + ((t235 * t342 + t236 * t341 + t255 * t340) * t251 + (t235 * t347 + t236 * t343 + t255 * t345) * t234 + (t235 * t348 + t236 * t344 + t255 * t346) * t233) * t233 / 0.2e1 + ((t237 * t342 + t238 * t341 + t257 * t340) * t251 + (t347 * t237 + t238 * t343 + t345 * t257) * t234 + (t237 * t348 + t238 * t344 + t257 * t346) * t233) * t234 / 0.2e1 + ((t253 * t342 + t254 * t341 - t323 * t340) * t251 + (t253 * t347 + t254 * t343 - t323 * t345) * t234 + (t253 * t348 + t254 * t344 - t323 * t346) * t233) * t251 / 0.2e1 + ((-t293 * t271 + t273 * t296 + Icges(1,4)) * V_base(5) + (-t293 * t272 + t296 * t274 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t271 + t293 * t273 + Icges(1,2)) * V_base(5) + (t272 * t296 + t293 * t274 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t293 + Icges(2,6) * t296) * V_base(5) + (Icges(2,5) * t296 - Icges(2,6) * t293) * V_base(4) + Icges(2,3) * t285 / 0.2e1) * t285;
T  = t1;
