% Calculate kinetic energy for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:47
% EndTime: 2019-03-09 17:21:50
% DurationCPUTime: 3.53s
% Computational Cost: add. (1623->314), mult. (3393->444), div. (0->0), fcn. (3696->8), ass. (0->150)
t356 = Icges(4,1) + Icges(5,1);
t355 = Icges(6,1) + Icges(7,1);
t354 = -Icges(4,4) + Icges(5,5);
t353 = Icges(5,4) + Icges(4,5);
t352 = -Icges(6,4) + Icges(7,5);
t351 = Icges(7,4) + Icges(6,5);
t350 = Icges(4,2) + Icges(5,3);
t349 = Icges(6,2) + Icges(7,3);
t348 = Icges(7,2) + Icges(6,3);
t347 = -Icges(5,6) + Icges(4,6);
t346 = Icges(6,6) - Icges(7,6);
t345 = -Icges(4,3) - Icges(5,2);
t344 = rSges(7,1) + pkin(5);
t343 = rSges(7,3) + qJ(6);
t276 = sin(qJ(3));
t279 = cos(qJ(3));
t281 = cos(qJ(1));
t278 = sin(qJ(1));
t280 = cos(qJ(2));
t309 = t278 * t280;
t234 = t276 * t309 + t279 * t281;
t235 = -t276 * t281 + t279 * t309;
t275 = sin(qJ(5));
t317 = cos(qJ(5));
t191 = -t234 * t317 + t235 * t275;
t192 = t234 * t275 + t235 * t317;
t277 = sin(qJ(2));
t312 = t277 * t278;
t342 = t349 * t191 + t352 * t192 + t346 * t312;
t308 = t280 * t281;
t236 = t276 * t308 - t278 * t279;
t237 = t278 * t276 + t279 * t308;
t193 = -t236 * t317 + t237 * t275;
t194 = t236 * t275 + t237 * t317;
t310 = t277 * t281;
t341 = t349 * t193 + t352 * t194 + t346 * t310;
t340 = -t346 * t191 + t351 * t192 - t348 * t312;
t339 = -t346 * t193 + t351 * t194 - t348 * t310;
t338 = t352 * t191 + t355 * t192 - t351 * t312;
t337 = t352 * t193 + t355 * t194 - t351 * t310;
t311 = t277 * t279;
t313 = t276 * t277;
t227 = t275 * t311 - t313 * t317;
t228 = (t275 * t276 + t279 * t317) * t277;
t336 = t349 * t227 + t352 * t228 - t346 * t280;
t335 = -t346 * t227 + t351 * t228 + t348 * t280;
t334 = t352 * t227 + t355 * t228 + t351 * t280;
t333 = t350 * t234 + t354 * t235 - t347 * t312;
t332 = t350 * t236 + t354 * t237 - t347 * t310;
t331 = -t347 * t234 + t353 * t235 - t345 * t312;
t330 = -t347 * t236 + t353 * t237 - t345 * t310;
t329 = t354 * t234 + t356 * t235 + t353 * t312;
t328 = t354 * t236 + t356 * t237 + t353 * t310;
t327 = t347 * t280 + (t350 * t276 + t354 * t279) * t277;
t326 = t345 * t280 + (-t347 * t276 + t353 * t279) * t277;
t325 = -t353 * t280 + (t354 * t276 + t356 * t279) * t277;
t316 = Icges(2,4) * t278;
t315 = Icges(3,4) * t277;
t314 = Icges(3,4) * t280;
t307 = -rSges(7,2) * t312 + t343 * t191 + t344 * t192;
t306 = -rSges(7,2) * t310 + t343 * t193 + t344 * t194;
t305 = rSges(7,2) * t280 + t343 * t227 + t344 * t228;
t304 = qJD(3) * t277;
t303 = qJD(5) * t277;
t302 = V_base(5) * pkin(6) + V_base(1);
t266 = qJD(2) * t278 + V_base(4);
t271 = V_base(6) + qJD(1);
t233 = t281 * t304 + t266;
t299 = pkin(2) * t280 + pkin(8) * t277;
t265 = -qJD(2) * t281 + V_base(5);
t298 = rSges(3,1) * t280 - rSges(3,2) * t277;
t297 = Icges(3,1) * t280 - t315;
t296 = -Icges(3,2) * t277 + t314;
t295 = Icges(3,5) * t280 - Icges(3,6) * t277;
t232 = t278 * t304 + t265;
t263 = pkin(1) * t281 + t278 * pkin(7);
t294 = -V_base(4) * pkin(6) + t271 * t263 + V_base(2);
t262 = t278 * pkin(1) - pkin(7) * t281;
t293 = V_base(4) * t262 - t263 * V_base(5) + V_base(3);
t240 = t299 * t278;
t261 = pkin(2) * t277 - pkin(8) * t280;
t292 = t265 * t261 + (-t240 - t262) * t271 + t302;
t291 = (-Icges(3,3) * t281 + t278 * t295) * t265 + (Icges(3,3) * t278 + t281 * t295) * t266 + (Icges(3,5) * t277 + Icges(3,6) * t280) * t271;
t241 = t299 * t281;
t290 = t271 * t241 - t261 * t266 + t294;
t238 = (pkin(3) * t279 + qJ(4) * t276) * t277;
t289 = qJD(4) * t236 + t232 * t238 + t292;
t288 = t266 * t240 - t241 * t265 + t293;
t196 = pkin(3) * t237 + qJ(4) * t236;
t257 = -qJD(3) * t280 + t271;
t287 = qJD(4) * t234 + t257 * t196 + t290;
t195 = pkin(3) * t235 + qJ(4) * t234;
t286 = qJD(4) * t313 + t233 * t195 + t288;
t203 = pkin(4) * t235 - pkin(9) * t312;
t244 = pkin(4) * t311 + pkin(9) * t280;
t285 = t232 * t244 + (-t195 - t203) * t257 + t289;
t204 = t237 * pkin(4) - pkin(9) * t310;
t284 = t257 * t204 + (-t238 - t244) * t233 + t287;
t283 = t233 * t203 + (-t196 - t204) * t232 + t286;
t215 = -Icges(3,6) * t281 + t278 * t296;
t216 = Icges(3,6) * t278 + t281 * t296;
t219 = -Icges(3,5) * t281 + t278 * t297;
t220 = Icges(3,5) * t278 + t281 * t297;
t251 = Icges(3,2) * t280 + t315;
t254 = Icges(3,1) * t277 + t314;
t282 = (-t216 * t277 + t220 * t280) * t266 + (-t215 * t277 + t219 * t280) * t265 + (-t251 * t277 + t254 * t280) * t271;
t273 = Icges(2,4) * t281;
t260 = rSges(2,1) * t281 - t278 * rSges(2,2);
t259 = t278 * rSges(2,1) + rSges(2,2) * t281;
t258 = rSges(3,1) * t277 + rSges(3,2) * t280;
t256 = Icges(2,1) * t281 - t316;
t255 = Icges(2,1) * t278 + t273;
t253 = -Icges(2,2) * t278 + t273;
t252 = Icges(2,2) * t281 + t316;
t247 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t246 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t245 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t242 = (-qJD(3) + qJD(5)) * t280 + t271;
t224 = t278 * rSges(3,3) + t281 * t298;
t223 = -rSges(3,3) * t281 + t278 * t298;
t222 = -rSges(4,3) * t280 + (rSges(4,1) * t279 - rSges(4,2) * t276) * t277;
t221 = -rSges(5,2) * t280 + (rSges(5,1) * t279 + rSges(5,3) * t276) * t277;
t207 = -t281 * t303 + t233;
t206 = -t278 * t303 + t232;
t202 = V_base(5) * rSges(2,3) - t259 * t271 + t302;
t201 = t260 * t271 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t199 = t259 * V_base(4) - t260 * V_base(5) + V_base(3);
t188 = t237 * rSges(4,1) - t236 * rSges(4,2) + rSges(4,3) * t310;
t187 = t237 * rSges(5,1) + rSges(5,2) * t310 + t236 * rSges(5,3);
t186 = rSges(4,1) * t235 - rSges(4,2) * t234 + rSges(4,3) * t312;
t185 = rSges(5,1) * t235 + rSges(5,2) * t312 + rSges(5,3) * t234;
t172 = rSges(6,1) * t228 - rSges(6,2) * t227 + rSges(6,3) * t280;
t162 = t258 * t265 + (-t223 - t262) * t271 + t302;
t161 = t224 * t271 - t258 * t266 + t294;
t158 = t223 * t266 - t224 * t265 + t293;
t157 = t194 * rSges(6,1) - t193 * rSges(6,2) - rSges(6,3) * t310;
t155 = rSges(6,1) * t192 - rSges(6,2) * t191 - rSges(6,3) * t312;
t141 = -t186 * t257 + t222 * t232 + t292;
t140 = t188 * t257 - t222 * t233 + t290;
t139 = t186 * t233 - t188 * t232 + t288;
t138 = t221 * t232 + (-t185 - t195) * t257 + t289;
t137 = t187 * t257 + (-t221 - t238) * t233 + t287;
t136 = t185 * t233 + (-t187 - t196) * t232 + t286;
t135 = -t155 * t242 + t172 * t206 + t285;
t134 = t157 * t242 - t172 * t207 + t284;
t133 = t155 * t207 - t157 * t206 + t283;
t132 = qJD(6) * t193 + t206 * t305 - t242 * t307 + t285;
t131 = qJD(6) * t191 - t207 * t305 + t242 * t306 + t284;
t130 = qJD(6) * t227 - t206 * t306 + t207 * t307 + t283;
t1 = t266 * (t278 * t291 + t281 * t282) / 0.2e1 + t265 * (t278 * t282 - t281 * t291) / 0.2e1 + m(6) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(5) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(4) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(7) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(3) * (t158 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(2) * (t199 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(1) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + ((t191 * t336 + t192 * t334 - t312 * t335) * t242 + (t191 * t341 + t192 * t337 - t312 * t339) * t207 + (t342 * t191 + t338 * t192 - t340 * t312) * t206) * t206 / 0.2e1 + ((t193 * t336 + t194 * t334 - t310 * t335) * t242 + (t341 * t193 + t337 * t194 - t339 * t310) * t207 + (t193 * t342 + t338 * t194 - t340 * t310) * t206) * t207 / 0.2e1 + ((t234 * t327 + t235 * t325 + t312 * t326) * t257 + (t234 * t332 + t235 * t328 + t312 * t330) * t233 + (t333 * t234 + t329 * t235 + t331 * t312) * t232) * t232 / 0.2e1 + ((t236 * t327 + t237 * t325 + t310 * t326) * t257 + (t332 * t236 + t328 * t237 + t330 * t310) * t233 + (t236 * t333 + t237 * t329 + t310 * t331) * t232) * t233 / 0.2e1 + ((t336 * t227 + t334 * t228 + t335 * t280) * t242 + (t227 * t341 + t228 * t337 + t280 * t339) * t207 + (t227 * t342 + t338 * t228 + t340 * t280) * t206) * t242 / 0.2e1 + ((-t232 * t331 - t233 * t330 - t257 * t326) * t280 + ((t276 * t327 + t279 * t325) * t257 + (t276 * t332 + t279 * t328) * t233 + (t276 * t333 + t279 * t329) * t232) * t277) * t257 / 0.2e1 + ((t216 * t280 + t220 * t277) * t266 + (t215 * t280 + t219 * t277) * t265 + (t280 * t251 + t277 * t254 + Icges(2,3)) * t271) * t271 / 0.2e1 + ((-t278 * t252 + t255 * t281 + Icges(1,4)) * V_base(5) + (-t278 * t253 + t281 * t256 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t281 * t252 + t278 * t255 + Icges(1,2)) * V_base(5) + (t253 * t281 + t278 * t256 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t271 * (Icges(2,5) * t281 - Icges(2,6) * t278) + V_base(5) * t271 * (Icges(2,5) * t278 + Icges(2,6) * t281) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
