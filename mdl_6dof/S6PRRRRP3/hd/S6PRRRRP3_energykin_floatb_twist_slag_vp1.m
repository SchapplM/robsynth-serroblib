% Calculate kinetic energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:15
% EndTime: 2019-03-09 00:06:20
% DurationCPUTime: 4.44s
% Computational Cost: add. (3371->391), mult. (7008->568), div. (0->0), fcn. (8493->12), ass. (0->175)
t370 = Icges(6,1) + Icges(7,1);
t369 = Icges(6,4) + Icges(7,4);
t368 = Icges(6,5) + Icges(7,5);
t367 = Icges(6,2) + Icges(7,2);
t366 = Icges(6,6) + Icges(7,6);
t365 = Icges(6,3) + Icges(7,3);
t364 = rSges(7,3) + qJ(6);
t300 = sin(pkin(11));
t302 = cos(pkin(11));
t308 = cos(qJ(2));
t303 = cos(pkin(6));
t306 = sin(qJ(2));
t337 = t303 * t306;
t259 = t300 * t308 + t302 * t337;
t301 = sin(pkin(6));
t305 = sin(qJ(3));
t339 = t301 * t305;
t349 = cos(qJ(3));
t243 = t259 * t349 - t302 * t339;
t336 = t303 * t308;
t258 = t300 * t306 - t302 * t336;
t299 = qJ(4) + qJ(5);
t294 = sin(t299);
t295 = cos(t299);
t208 = -t243 * t294 + t258 * t295;
t209 = t243 * t295 + t258 * t294;
t323 = t301 * t349;
t242 = t259 * t305 + t302 * t323;
t361 = t366 * t208 + t368 * t209 + t365 * t242;
t261 = -t300 * t337 + t302 * t308;
t245 = t261 * t349 + t300 * t339;
t260 = t300 * t336 + t302 * t306;
t210 = -t245 * t294 + t260 * t295;
t211 = t245 * t295 + t260 * t294;
t244 = t261 * t305 - t300 * t323;
t360 = t366 * t210 + t368 * t211 + t365 * t244;
t359 = t367 * t208 + t369 * t209 + t366 * t242;
t358 = t367 * t210 + t369 * t211 + t366 * t244;
t357 = t369 * t208 + t370 * t209 + t368 * t242;
t356 = t369 * t210 + t370 * t211 + t368 * t244;
t266 = t303 * t305 + t306 * t323;
t338 = t301 * t308;
t235 = -t266 * t294 - t295 * t338;
t236 = t266 * t295 - t294 * t338;
t265 = -t303 * t349 + t306 * t339;
t355 = t366 * t235 + t368 * t236 + t365 * t265;
t354 = t367 * t235 + t369 * t236 + t366 * t265;
t353 = t369 * t235 + t370 * t236 + t368 * t265;
t347 = pkin(7) * t303;
t307 = cos(qJ(4));
t346 = t307 * pkin(4);
t344 = Icges(2,4) * t300;
t304 = sin(qJ(4));
t343 = t258 * t304;
t342 = t260 * t304;
t341 = t300 * t301;
t340 = t301 * t302;
t322 = pkin(5) * t294;
t332 = pkin(5) * t295;
t335 = rSges(7,1) * t209 + rSges(7,2) * t208 + t242 * t364 + t243 * t332 + t258 * t322;
t334 = rSges(7,1) * t211 + rSges(7,2) * t210 + t244 * t364 + t245 * t332 + t260 * t322;
t333 = rSges(7,1) * t236 + rSges(7,2) * t235 + t265 * t364 + t266 * t332 - t322 * t338;
t330 = qJD(2) * t301;
t329 = V_base(5) * qJ(1) + V_base(1);
t325 = qJD(1) + V_base(3);
t324 = t304 * t338;
t275 = t300 * t330 + V_base(4);
t287 = qJD(2) * t303 + V_base(6);
t241 = qJD(3) * t260 + t275;
t207 = qJD(4) * t244 + t241;
t274 = -t302 * t330 + V_base(5);
t240 = qJD(3) * t258 + t274;
t262 = -qJD(3) * t338 + t287;
t268 = pkin(1) * t300 - pkin(7) * t340;
t321 = -t268 * V_base(6) + V_base(5) * t347 + t329;
t269 = pkin(1) * t302 + pkin(7) * t341;
t320 = V_base(4) * t268 - t269 * V_base(5) + t325;
t206 = qJD(4) * t242 + t240;
t234 = qJD(4) * t265 + t262;
t319 = V_base(6) * t269 + V_base(2) + (-qJ(1) - t347) * V_base(4);
t231 = pkin(2) * t259 + pkin(8) * t258;
t267 = (pkin(2) * t306 - pkin(8) * t308) * t301;
t318 = -t231 * t287 + t274 * t267 + t321;
t232 = pkin(2) * t261 + pkin(8) * t260;
t317 = t275 * t231 - t232 * t274 + t320;
t316 = t287 * t232 - t267 * t275 + t319;
t204 = t243 * pkin(3) + t242 * pkin(9);
t233 = t266 * pkin(3) + t265 * pkin(9);
t315 = -t204 * t262 + t240 * t233 + t318;
t205 = t245 * pkin(3) + t244 * pkin(9);
t314 = t241 * t204 - t205 * t240 + t317;
t313 = t262 * t205 - t233 * t241 + t316;
t147 = pkin(4) * t343 + pkin(10) * t242 + t243 * t346;
t189 = -pkin(4) * t324 + pkin(10) * t265 + t266 * t346;
t312 = -t147 * t234 + t206 * t189 + t315;
t148 = pkin(4) * t342 + pkin(10) * t244 + t245 * t346;
t311 = t207 * t147 - t148 * t206 + t314;
t310 = t234 * t148 - t189 * t207 + t313;
t293 = Icges(2,4) * t302;
t284 = rSges(2,1) * t302 - rSges(2,2) * t300;
t283 = rSges(2,1) * t300 + rSges(2,2) * t302;
t282 = Icges(2,1) * t302 - t344;
t281 = Icges(2,1) * t300 + t293;
t280 = -Icges(2,2) * t300 + t293;
t279 = Icges(2,2) * t302 + t344;
t273 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t272 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t271 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t254 = rSges(3,3) * t303 + (rSges(3,1) * t306 + rSges(3,2) * t308) * t301;
t253 = Icges(3,5) * t303 + (Icges(3,1) * t306 + Icges(3,4) * t308) * t301;
t252 = Icges(3,6) * t303 + (Icges(3,4) * t306 + Icges(3,2) * t308) * t301;
t251 = Icges(3,3) * t303 + (Icges(3,5) * t306 + Icges(3,6) * t308) * t301;
t250 = V_base(5) * rSges(2,3) - t283 * V_base(6) + t329;
t249 = t284 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t247 = t266 * t307 - t324;
t246 = -t266 * t304 - t307 * t338;
t239 = t283 * V_base(4) - t284 * V_base(5) + t325;
t230 = rSges(4,1) * t266 - rSges(4,2) * t265 - rSges(4,3) * t338;
t229 = Icges(4,1) * t266 - Icges(4,4) * t265 - Icges(4,5) * t338;
t228 = Icges(4,4) * t266 - Icges(4,2) * t265 - Icges(4,6) * t338;
t227 = Icges(4,5) * t266 - Icges(4,6) * t265 - Icges(4,3) * t338;
t226 = rSges(3,1) * t261 - rSges(3,2) * t260 + rSges(3,3) * t341;
t225 = rSges(3,1) * t259 - rSges(3,2) * t258 - rSges(3,3) * t340;
t224 = Icges(3,1) * t261 - Icges(3,4) * t260 + Icges(3,5) * t341;
t223 = Icges(3,1) * t259 - Icges(3,4) * t258 - Icges(3,5) * t340;
t222 = Icges(3,4) * t261 - Icges(3,2) * t260 + Icges(3,6) * t341;
t221 = Icges(3,4) * t259 - Icges(3,2) * t258 - Icges(3,6) * t340;
t220 = Icges(3,5) * t261 - Icges(3,6) * t260 + Icges(3,3) * t341;
t219 = Icges(3,5) * t259 - Icges(3,6) * t258 - Icges(3,3) * t340;
t216 = qJD(5) * t265 + t234;
t215 = t245 * t307 + t342;
t214 = -t245 * t304 + t260 * t307;
t213 = t243 * t307 + t343;
t212 = -t243 * t304 + t258 * t307;
t202 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t265;
t200 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t265;
t199 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t265;
t198 = Icges(5,5) * t247 + Icges(5,6) * t246 + Icges(5,3) * t265;
t197 = rSges(4,1) * t245 - rSges(4,2) * t244 + rSges(4,3) * t260;
t196 = rSges(4,1) * t243 - rSges(4,2) * t242 + rSges(4,3) * t258;
t195 = Icges(4,1) * t245 - Icges(4,4) * t244 + Icges(4,5) * t260;
t194 = Icges(4,1) * t243 - Icges(4,4) * t242 + Icges(4,5) * t258;
t193 = Icges(4,4) * t245 - Icges(4,2) * t244 + Icges(4,6) * t260;
t192 = Icges(4,4) * t243 - Icges(4,2) * t242 + Icges(4,6) * t258;
t191 = Icges(4,5) * t245 - Icges(4,6) * t244 + Icges(4,3) * t260;
t190 = Icges(4,5) * t243 - Icges(4,6) * t242 + Icges(4,3) * t258;
t188 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t265;
t180 = qJD(5) * t244 + t207;
t179 = qJD(5) * t242 + t206;
t177 = -t225 * t287 + t254 * t274 + t321;
t176 = t226 * t287 - t254 * t275 + t319;
t174 = rSges(5,1) * t215 + rSges(5,2) * t214 + rSges(5,3) * t244;
t173 = rSges(5,1) * t213 + rSges(5,2) * t212 + rSges(5,3) * t242;
t172 = Icges(5,1) * t215 + Icges(5,4) * t214 + Icges(5,5) * t244;
t171 = Icges(5,1) * t213 + Icges(5,4) * t212 + Icges(5,5) * t242;
t170 = Icges(5,4) * t215 + Icges(5,2) * t214 + Icges(5,6) * t244;
t169 = Icges(5,4) * t213 + Icges(5,2) * t212 + Icges(5,6) * t242;
t168 = Icges(5,5) * t215 + Icges(5,6) * t214 + Icges(5,3) * t244;
t167 = Icges(5,5) * t213 + Icges(5,6) * t212 + Icges(5,3) * t242;
t165 = t225 * t275 - t226 * t274 + t320;
t164 = rSges(6,1) * t211 + rSges(6,2) * t210 + rSges(6,3) * t244;
t162 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t242;
t142 = -t196 * t262 + t230 * t240 + t318;
t141 = t197 * t262 - t230 * t241 + t316;
t140 = t196 * t241 - t197 * t240 + t317;
t139 = -t173 * t234 + t202 * t206 + t315;
t138 = t174 * t234 - t202 * t207 + t313;
t137 = t173 * t207 - t174 * t206 + t314;
t136 = -t162 * t216 + t179 * t188 + t312;
t135 = t164 * t216 - t180 * t188 + t310;
t134 = t162 * t180 - t164 * t179 + t311;
t133 = qJD(6) * t244 + t179 * t333 - t216 * t335 + t312;
t132 = qJD(6) * t242 - t180 * t333 + t216 * t334 + t310;
t131 = qJD(6) * t265 - t179 * t334 + t180 * t335 + t311;
t1 = t287 * ((t219 * t274 + t220 * t275 + t251 * t287) * t303 + ((t222 * t308 + t224 * t306) * t275 + (t221 * t308 + t223 * t306) * t274 + (t252 * t308 + t253 * t306) * t287) * t301) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t262 * ((-t191 * t338 - t193 * t265 + t195 * t266) * t241 + (-t190 * t338 - t192 * t265 + t194 * t266) * t240 + (-t227 * t338 - t228 * t265 + t229 * t266) * t262) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + t274 * ((-t220 * t340 - t222 * t258 + t224 * t259) * t275 + (-t219 * t340 - t221 * t258 + t223 * t259) * t274 + (-t251 * t340 - t252 * t258 + t253 * t259) * t287) / 0.2e1 + t275 * ((t220 * t341 - t222 * t260 + t224 * t261) * t275 + (t219 * t341 - t221 * t260 + t223 * t261) * t274 + (t251 * t341 - t252 * t260 + t253 * t261) * t287) / 0.2e1 + m(3) * (t165 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + t206 * ((t168 * t242 + t170 * t212 + t172 * t213) * t207 + (t167 * t242 + t169 * t212 + t171 * t213) * t206 + (t198 * t242 + t199 * t212 + t200 * t213) * t234) / 0.2e1 + t207 * ((t168 * t244 + t170 * t214 + t172 * t215) * t207 + (t167 * t244 + t169 * t214 + t171 * t215) * t206 + (t198 * t244 + t199 * t214 + t200 * t215) * t234) / 0.2e1 + m(2) * (t239 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + t240 * ((t191 * t258 - t193 * t242 + t195 * t243) * t241 + (t190 * t258 - t192 * t242 + t194 * t243) * t240 + (t227 * t258 - t228 * t242 + t229 * t243) * t262) / 0.2e1 + t241 * ((t191 * t260 - t193 * t244 + t195 * t245) * t241 + (t190 * t260 - t192 * t244 + t194 * t245) * t240 + (t227 * t260 - t228 * t244 + t229 * t245) * t262) / 0.2e1 + t234 * ((t168 * t265 + t170 * t246 + t172 * t247) * t207 + (t167 * t265 + t169 * t246 + t171 * t247) * t206 + (t198 * t265 + t199 * t246 + t200 * t247) * t234) / 0.2e1 + m(1) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + ((t208 * t354 + t209 * t353 + t242 * t355) * t216 + (t208 * t358 + t209 * t356 + t242 * t360) * t180 + (t359 * t208 + t357 * t209 + t361 * t242) * t179) * t179 / 0.2e1 + ((t210 * t354 + t211 * t353 + t244 * t355) * t216 + (t358 * t210 + t356 * t211 + t360 * t244) * t180 + (t210 * t359 + t211 * t357 + t244 * t361) * t179) * t180 / 0.2e1 + ((t354 * t235 + t353 * t236 + t355 * t265) * t216 + (t235 * t358 + t236 * t356 + t265 * t360) * t180 + (t235 * t359 + t236 * t357 + t265 * t361) * t179) * t216 / 0.2e1 + ((-t279 * t300 + t281 * t302 + Icges(1,4)) * V_base(5) + (-t280 * t300 + t282 * t302 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t279 * t302 + t281 * t300 + Icges(1,2)) * V_base(5) + (t280 * t302 + t282 * t300 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t300 + Icges(2,6) * t302 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t302 - Icges(2,6) * t300 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
