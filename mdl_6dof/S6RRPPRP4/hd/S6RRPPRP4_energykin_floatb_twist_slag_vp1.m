% Calculate kinetic energy for
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:36:56
% EndTime: 2019-03-09 08:36:59
% DurationCPUTime: 3.13s
% Computational Cost: add. (1527->316), mult. (3233->427), div. (0->0), fcn. (3536->8), ass. (0->151)
t355 = Icges(4,1) + Icges(5,1);
t354 = Icges(6,1) + Icges(7,1);
t353 = -Icges(4,4) + Icges(5,5);
t352 = Icges(5,4) + Icges(4,5);
t351 = -Icges(6,4) + Icges(7,5);
t350 = Icges(7,4) + Icges(6,5);
t349 = Icges(4,2) + Icges(5,3);
t348 = Icges(6,2) + Icges(7,3);
t347 = Icges(7,2) + Icges(6,3);
t346 = -Icges(5,6) + Icges(4,6);
t345 = Icges(6,6) - Icges(7,6);
t344 = -Icges(4,3) - Icges(5,2);
t343 = rSges(7,1) + pkin(5);
t342 = rSges(7,3) + qJ(6);
t271 = sin(pkin(9));
t272 = cos(pkin(9));
t277 = cos(qJ(1));
t275 = sin(qJ(1));
t276 = cos(qJ(2));
t309 = t275 * t276;
t229 = t271 * t309 + t272 * t277;
t230 = -t271 * t277 + t272 * t309;
t273 = sin(qJ(5));
t317 = cos(qJ(5));
t189 = -t229 * t317 + t230 * t273;
t190 = t229 * t273 + t230 * t317;
t274 = sin(qJ(2));
t311 = t274 * t275;
t341 = t189 * t348 + t190 * t351 + t311 * t345;
t308 = t276 * t277;
t231 = t271 * t308 - t275 * t272;
t232 = t275 * t271 + t272 * t308;
t191 = -t231 * t317 + t232 * t273;
t192 = t231 * t273 + t232 * t317;
t310 = t274 * t277;
t340 = t191 * t348 + t192 * t351 + t310 * t345;
t339 = -t189 * t345 + t190 * t350 - t311 * t347;
t338 = -t191 * t345 + t192 * t350 - t310 * t347;
t337 = t189 * t351 + t190 * t354 - t311 * t350;
t336 = t191 * t351 + t192 * t354 - t310 * t350;
t312 = t272 * t274;
t313 = t271 * t274;
t224 = t273 * t312 - t313 * t317;
t225 = (t271 * t273 + t272 * t317) * t274;
t335 = t224 * t348 + t225 * t351 - t276 * t345;
t334 = -t224 * t345 + t225 * t350 + t276 * t347;
t333 = t224 * t351 + t225 * t354 + t276 * t350;
t332 = t229 * t349 + t230 * t353 - t311 * t346;
t331 = t231 * t349 + t232 * t353 - t310 * t346;
t330 = -t229 * t346 + t230 * t352 - t311 * t344;
t329 = -t231 * t346 + t232 * t352 - t310 * t344;
t328 = t353 * t229 + t230 * t355 + t352 * t311;
t327 = t353 * t231 + t232 * t355 + t352 * t310;
t326 = t346 * t276 + (t271 * t349 + t272 * t353) * t274;
t325 = t344 * t276 + (-t271 * t346 + t272 * t352) * t274;
t324 = -t352 * t276 + (t353 * t271 + t272 * t355) * t274;
t316 = Icges(2,4) * t275;
t315 = Icges(3,4) * t274;
t314 = Icges(3,4) * t276;
t307 = -rSges(7,2) * t311 + t342 * t189 + t343 * t190;
t306 = -rSges(7,2) * t310 + t342 * t191 + t343 * t192;
t305 = rSges(7,2) * t276 + t342 * t224 + t343 * t225;
t194 = pkin(3) * t232 + qJ(4) * t231;
t293 = pkin(2) * t276 + qJ(3) * t274;
t237 = t293 * t277;
t304 = -t194 - t237;
t235 = (pkin(3) * t272 + qJ(4) * t271) * t274;
t254 = pkin(2) * t274 - qJ(3) * t276;
t303 = -t235 - t254;
t236 = t293 * t275;
t258 = t275 * pkin(1) - pkin(7) * t277;
t302 = -t236 - t258;
t301 = qJD(3) * t274;
t300 = qJD(5) * t274;
t299 = V_base(5) * pkin(6) + V_base(1);
t193 = pkin(3) * t230 + qJ(4) * t229;
t296 = -t193 + t302;
t262 = qJD(2) * t275 + V_base(4);
t267 = V_base(6) + qJD(1);
t261 = -qJD(2) * t277 + V_base(5);
t295 = t261 * t254 + t277 * t301 + t299;
t294 = rSges(3,1) * t276 - rSges(3,2) * t274;
t292 = Icges(3,1) * t276 - t315;
t291 = -Icges(3,2) * t274 + t314;
t290 = Icges(3,5) * t276 - Icges(3,6) * t274;
t259 = pkin(1) * t277 + t275 * pkin(7);
t289 = -V_base(4) * pkin(6) + t267 * t259 + V_base(2);
t288 = V_base(4) * t258 - t259 * V_base(5) + V_base(3);
t287 = qJD(4) * t231 + t261 * t235 + t295;
t286 = (-Icges(3,3) * t277 + t275 * t290) * t261 + (Icges(3,3) * t275 + t277 * t290) * t262 + (Icges(3,5) * t274 + Icges(3,6) * t276) * t267;
t285 = t267 * t237 + t275 * t301 + t289;
t284 = -qJD(3) * t276 + t262 * t236 + t288;
t283 = qJD(4) * t229 + t267 * t194 + t285;
t282 = qJD(4) * t313 + t262 * t193 + t284;
t200 = pkin(4) * t230 - pkin(8) * t311;
t240 = pkin(4) * t312 + pkin(8) * t276;
t281 = t261 * t240 + (-t200 + t296) * t267 + t287;
t201 = t232 * pkin(4) - pkin(8) * t310;
t280 = t267 * t201 + (-t240 + t303) * t262 + t283;
t279 = t262 * t200 + (-t201 + t304) * t261 + t282;
t216 = -Icges(3,6) * t277 + t275 * t291;
t217 = Icges(3,6) * t275 + t277 * t291;
t218 = -Icges(3,5) * t277 + t275 * t292;
t219 = Icges(3,5) * t275 + t277 * t292;
t247 = Icges(3,2) * t276 + t315;
t250 = Icges(3,1) * t274 + t314;
t278 = (-t217 * t274 + t219 * t276) * t262 + (-t216 * t274 + t218 * t276) * t261 + (-t247 * t274 + t250 * t276) * t267;
t269 = Icges(2,4) * t277;
t257 = rSges(2,1) * t277 - t275 * rSges(2,2);
t256 = t275 * rSges(2,1) + rSges(2,2) * t277;
t255 = rSges(3,1) * t274 + rSges(3,2) * t276;
t253 = qJD(5) * t276 + t267;
t252 = Icges(2,1) * t277 - t316;
t251 = Icges(2,1) * t275 + t269;
t249 = -Icges(2,2) * t275 + t269;
t248 = Icges(2,2) * t277 + t316;
t243 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t242 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t241 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t234 = -t277 * t300 + t262;
t233 = -t275 * t300 + t261;
t223 = t275 * rSges(3,3) + t277 * t294;
t222 = -rSges(3,3) * t277 + t275 * t294;
t213 = -rSges(4,3) * t276 + (rSges(4,1) * t272 - rSges(4,2) * t271) * t274;
t212 = -rSges(5,2) * t276 + (rSges(5,1) * t272 + rSges(5,3) * t271) * t274;
t199 = V_base(5) * rSges(2,3) - t256 * t267 + t299;
t198 = t257 * t267 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t197 = t256 * V_base(4) - t257 * V_base(5) + V_base(3);
t186 = t232 * rSges(4,1) - t231 * rSges(4,2) + rSges(4,3) * t310;
t185 = t232 * rSges(5,1) + rSges(5,2) * t310 + t231 * rSges(5,3);
t184 = rSges(4,1) * t230 - rSges(4,2) * t229 + rSges(4,3) * t311;
t183 = rSges(5,1) * t230 + rSges(5,2) * t311 + rSges(5,3) * t229;
t169 = rSges(6,1) * t225 - rSges(6,2) * t224 + rSges(6,3) * t276;
t161 = t255 * t261 + (-t222 - t258) * t267 + t299;
t160 = t223 * t267 - t255 * t262 + t289;
t157 = t222 * t262 - t223 * t261 + t288;
t156 = t192 * rSges(6,1) - t191 * rSges(6,2) - rSges(6,3) * t310;
t154 = rSges(6,1) * t190 - rSges(6,2) * t189 - rSges(6,3) * t311;
t140 = t213 * t261 + (-t184 + t302) * t267 + t295;
t139 = t186 * t267 + (-t213 - t254) * t262 + t285;
t138 = t184 * t262 + (-t186 - t237) * t261 + t284;
t137 = t212 * t261 + (-t183 + t296) * t267 + t287;
t136 = t185 * t267 + (-t212 + t303) * t262 + t283;
t135 = t183 * t262 + (-t185 + t304) * t261 + t282;
t134 = -t154 * t253 + t169 * t233 + t281;
t133 = t156 * t253 - t169 * t234 + t280;
t132 = t154 * t234 - t156 * t233 + t279;
t131 = qJD(6) * t191 + t233 * t305 - t253 * t307 + t281;
t130 = qJD(6) * t189 - t234 * t305 + t253 * t306 + t280;
t129 = qJD(6) * t224 - t233 * t306 + t234 * t307 + t279;
t1 = m(5) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(4) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(6) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t157 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(1) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + m(2) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + ((t189 * t335 + t190 * t333 - t311 * t334) * t253 + (t189 * t340 + t190 * t336 - t311 * t338) * t234 + (t341 * t189 + t337 * t190 - t339 * t311) * t233) * t233 / 0.2e1 + ((t191 * t335 + t192 * t333 - t310 * t334) * t253 + (t340 * t191 + t336 * t192 - t338 * t310) * t234 + (t191 * t341 + t337 * t192 - t339 * t310) * t233) * t234 / 0.2e1 + ((t335 * t224 + t333 * t225 + t334 * t276) * t253 + (t224 * t340 + t225 * t336 + t276 * t338) * t234 + (t224 * t341 + t337 * t225 + t339 * t276) * t233) * t253 / 0.2e1 + ((-t275 * t248 + t251 * t277 + Icges(1,4)) * V_base(5) + (-t275 * t249 + t252 * t277 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t248 * t277 + t275 * t251 + Icges(1,2)) * V_base(5) + (t249 * t277 + t275 * t252 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t275 * t278 - t286 * t277 + (t229 * t326 + t230 * t324 + t311 * t325) * t267 + (t229 * t331 + t327 * t230 + t329 * t311) * t262 + (t332 * t229 + t328 * t230 + t330 * t311) * t261) * t261 / 0.2e1 + (t275 * t286 + t277 * t278 + (t231 * t326 + t232 * t324 + t310 * t325) * t267 + (t331 * t231 + t327 * t232 + t329 * t310) * t262 + (t231 * t332 + t232 * t328 + t310 * t330) * t261) * t262 / 0.2e1 + (((t217 - t329) * t262 + (t216 - t330) * t261) * t276 + ((t271 * t331 + t272 * t327 + t219) * t262 + (t271 * t332 + t272 * t328 + t218) * t261) * t274 + (Icges(2,3) + (t247 - t325) * t276 + (t271 * t326 + t272 * t324 + t250) * t274) * t267) * t267 / 0.2e1 + t267 * V_base(4) * (Icges(2,5) * t277 - Icges(2,6) * t275) + t267 * V_base(5) * (Icges(2,5) * t275 + Icges(2,6) * t277) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
