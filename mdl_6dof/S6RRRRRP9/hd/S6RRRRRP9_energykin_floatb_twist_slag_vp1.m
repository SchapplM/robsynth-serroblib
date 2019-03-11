% Calculate kinetic energy for
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:11
% EndTime: 2019-03-10 02:00:16
% DurationCPUTime: 4.54s
% Computational Cost: add. (3431->391), mult. (7008->570), div. (0->0), fcn. (8493->12), ass. (0->176)
t368 = Icges(6,1) + Icges(7,1);
t367 = Icges(6,4) + Icges(7,4);
t366 = Icges(6,5) + Icges(7,5);
t365 = Icges(6,2) + Icges(7,2);
t364 = Icges(6,6) + Icges(7,6);
t363 = Icges(6,3) + Icges(7,3);
t362 = rSges(7,3) + qJ(6);
t302 = cos(pkin(6));
t306 = sin(qJ(1));
t308 = cos(qJ(2));
t336 = t306 * t308;
t305 = sin(qJ(2));
t309 = cos(qJ(1));
t337 = t305 * t309;
t263 = t302 * t337 + t336;
t304 = sin(qJ(3));
t301 = sin(pkin(6));
t339 = t301 * t309;
t349 = cos(qJ(3));
t244 = t263 * t349 - t304 * t339;
t335 = t308 * t309;
t338 = t305 * t306;
t262 = -t302 * t335 + t338;
t300 = qJ(4) + qJ(5);
t294 = sin(t300);
t295 = cos(t300);
t208 = -t244 * t294 + t262 * t295;
t209 = t244 * t295 + t262 * t294;
t324 = t301 * t349;
t243 = t263 * t304 + t309 * t324;
t361 = t208 * t364 + t209 * t366 + t243 * t363;
t265 = -t302 * t338 + t335;
t341 = t301 * t306;
t246 = t265 * t349 + t304 * t341;
t264 = t302 * t336 + t337;
t210 = -t246 * t294 + t264 * t295;
t211 = t246 * t295 + t264 * t294;
t245 = t265 * t304 - t306 * t324;
t360 = t210 * t364 + t211 * t366 + t245 * t363;
t359 = t208 * t365 + t209 * t367 + t243 * t364;
t358 = t210 * t365 + t211 * t367 + t245 * t364;
t357 = t367 * t208 + t209 * t368 + t366 * t243;
t356 = t367 * t210 + t211 * t368 + t366 * t245;
t261 = t302 * t304 + t305 * t324;
t340 = t301 * t308;
t235 = -t261 * t294 - t295 * t340;
t236 = t261 * t295 - t294 * t340;
t260 = t301 * t304 * t305 - t302 * t349;
t355 = t235 * t364 + t236 * t366 + t260 * t363;
t354 = t235 * t365 + t236 * t367 + t260 * t364;
t353 = t367 * t235 + t236 * t368 + t366 * t260;
t347 = pkin(8) * t302;
t307 = cos(qJ(4));
t346 = t307 * pkin(4);
t344 = Icges(2,4) * t306;
t303 = sin(qJ(4));
t343 = t262 * t303;
t342 = t264 * t303;
t323 = pkin(5) * t294;
t331 = pkin(5) * t295;
t334 = rSges(7,1) * t209 + rSges(7,2) * t208 + t243 * t362 + t244 * t331 + t262 * t323;
t333 = rSges(7,1) * t211 + rSges(7,2) * t210 + t245 * t362 + t246 * t331 + t264 * t323;
t332 = rSges(7,1) * t236 + rSges(7,2) * t235 + t260 * t362 + t261 * t331 - t323 * t340;
t329 = qJD(2) * t301;
t328 = V_base(5) * pkin(7) + V_base(1);
t325 = t303 * t340;
t276 = t306 * t329 + V_base(4);
t293 = V_base(6) + qJD(1);
t242 = qJD(3) * t264 + t276;
t277 = qJD(2) * t302 + t293;
t207 = qJD(4) * t245 + t242;
t275 = -t309 * t329 + V_base(5);
t268 = pkin(1) * t306 - pkin(8) * t339;
t322 = -t268 * t293 + V_base(5) * t347 + t328;
t269 = pkin(1) * t309 + pkin(8) * t341;
t321 = V_base(4) * t268 - t269 * V_base(5) + V_base(3);
t241 = qJD(3) * t262 + t275;
t206 = qJD(4) * t243 + t241;
t258 = -qJD(3) * t340 + t277;
t320 = t293 * t269 + V_base(2) + (-pkin(7) - t347) * V_base(4);
t232 = qJD(4) * t260 + t258;
t233 = pkin(2) * t263 + pkin(9) * t262;
t267 = (pkin(2) * t305 - pkin(9) * t308) * t301;
t319 = -t233 * t277 + t275 * t267 + t322;
t234 = pkin(2) * t265 + pkin(9) * t264;
t318 = t276 * t233 - t234 * t275 + t321;
t317 = t277 * t234 - t267 * t276 + t320;
t204 = t244 * pkin(3) + t243 * pkin(10);
t231 = t261 * pkin(3) + t260 * pkin(10);
t316 = -t204 * t258 + t241 * t231 + t319;
t205 = t246 * pkin(3) + t245 * pkin(10);
t315 = t242 * t204 - t205 * t241 + t318;
t314 = t258 * t205 - t231 * t242 + t317;
t147 = pkin(4) * t343 + pkin(11) * t243 + t244 * t346;
t189 = -pkin(4) * t325 + pkin(11) * t260 + t261 * t346;
t313 = -t147 * t232 + t206 * t189 + t316;
t148 = pkin(4) * t342 + pkin(11) * t245 + t246 * t346;
t312 = t207 * t147 - t148 * t206 + t315;
t311 = t232 * t148 - t189 * t207 + t314;
t296 = Icges(2,4) * t309;
t285 = rSges(2,1) * t309 - rSges(2,2) * t306;
t284 = rSges(2,1) * t306 + rSges(2,2) * t309;
t283 = Icges(2,1) * t309 - t344;
t282 = Icges(2,1) * t306 + t296;
t281 = -Icges(2,2) * t306 + t296;
t280 = Icges(2,2) * t309 + t344;
t273 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t272 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t271 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t254 = rSges(3,3) * t302 + (rSges(3,1) * t305 + rSges(3,2) * t308) * t301;
t253 = Icges(3,5) * t302 + (Icges(3,1) * t305 + Icges(3,4) * t308) * t301;
t252 = Icges(3,6) * t302 + (Icges(3,4) * t305 + Icges(3,2) * t308) * t301;
t251 = Icges(3,3) * t302 + (Icges(3,5) * t305 + Icges(3,6) * t308) * t301;
t250 = V_base(5) * rSges(2,3) - t284 * t293 + t328;
t249 = t285 * t293 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t247 = t284 * V_base(4) - t285 * V_base(5) + V_base(3);
t240 = t261 * t307 - t325;
t239 = -t261 * t303 - t307 * t340;
t230 = rSges(3,1) * t265 - rSges(3,2) * t264 + rSges(3,3) * t341;
t229 = rSges(3,1) * t263 - rSges(3,2) * t262 - rSges(3,3) * t339;
t228 = Icges(3,1) * t265 - Icges(3,4) * t264 + Icges(3,5) * t341;
t227 = Icges(3,1) * t263 - Icges(3,4) * t262 - Icges(3,5) * t339;
t226 = Icges(3,4) * t265 - Icges(3,2) * t264 + Icges(3,6) * t341;
t225 = Icges(3,4) * t263 - Icges(3,2) * t262 - Icges(3,6) * t339;
t224 = Icges(3,5) * t265 - Icges(3,6) * t264 + Icges(3,3) * t341;
t223 = Icges(3,5) * t263 - Icges(3,6) * t262 - Icges(3,3) * t339;
t222 = rSges(4,1) * t261 - rSges(4,2) * t260 - rSges(4,3) * t340;
t221 = Icges(4,1) * t261 - Icges(4,4) * t260 - Icges(4,5) * t340;
t220 = Icges(4,4) * t261 - Icges(4,2) * t260 - Icges(4,6) * t340;
t219 = Icges(4,5) * t261 - Icges(4,6) * t260 - Icges(4,3) * t340;
t216 = t246 * t307 + t342;
t215 = -t246 * t303 + t264 * t307;
t214 = t244 * t307 + t343;
t213 = -t244 * t303 + t262 * t307;
t212 = qJD(5) * t260 + t232;
t202 = rSges(4,1) * t246 - rSges(4,2) * t245 + rSges(4,3) * t264;
t201 = rSges(4,1) * t244 - rSges(4,2) * t243 + rSges(4,3) * t262;
t199 = Icges(4,1) * t246 - Icges(4,4) * t245 + Icges(4,5) * t264;
t198 = Icges(4,1) * t244 - Icges(4,4) * t243 + Icges(4,5) * t262;
t197 = Icges(4,4) * t246 - Icges(4,2) * t245 + Icges(4,6) * t264;
t196 = Icges(4,4) * t244 - Icges(4,2) * t243 + Icges(4,6) * t262;
t195 = Icges(4,5) * t246 - Icges(4,6) * t245 + Icges(4,3) * t264;
t194 = Icges(4,5) * t244 - Icges(4,6) * t243 + Icges(4,3) * t262;
t193 = rSges(5,1) * t240 + rSges(5,2) * t239 + rSges(5,3) * t260;
t192 = Icges(5,1) * t240 + Icges(5,4) * t239 + Icges(5,5) * t260;
t191 = Icges(5,4) * t240 + Icges(5,2) * t239 + Icges(5,6) * t260;
t190 = Icges(5,5) * t240 + Icges(5,6) * t239 + Icges(5,3) * t260;
t188 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t260;
t180 = qJD(5) * t245 + t207;
t179 = qJD(5) * t243 + t206;
t177 = -t229 * t277 + t254 * t275 + t322;
t176 = t230 * t277 - t254 * t276 + t320;
t174 = rSges(5,1) * t216 + rSges(5,2) * t215 + rSges(5,3) * t245;
t173 = rSges(5,1) * t214 + rSges(5,2) * t213 + rSges(5,3) * t243;
t172 = Icges(5,1) * t216 + Icges(5,4) * t215 + Icges(5,5) * t245;
t171 = Icges(5,1) * t214 + Icges(5,4) * t213 + Icges(5,5) * t243;
t170 = Icges(5,4) * t216 + Icges(5,2) * t215 + Icges(5,6) * t245;
t169 = Icges(5,4) * t214 + Icges(5,2) * t213 + Icges(5,6) * t243;
t168 = Icges(5,5) * t216 + Icges(5,6) * t215 + Icges(5,3) * t245;
t167 = Icges(5,5) * t214 + Icges(5,6) * t213 + Icges(5,3) * t243;
t166 = t229 * t276 - t230 * t275 + t321;
t165 = rSges(6,1) * t211 + rSges(6,2) * t210 + rSges(6,3) * t245;
t163 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t243;
t142 = -t201 * t258 + t222 * t241 + t319;
t141 = t202 * t258 - t222 * t242 + t317;
t140 = t201 * t242 - t202 * t241 + t318;
t139 = -t173 * t232 + t193 * t206 + t316;
t138 = t174 * t232 - t193 * t207 + t314;
t137 = t173 * t207 - t174 * t206 + t315;
t136 = -t163 * t212 + t179 * t188 + t313;
t135 = t165 * t212 - t180 * t188 + t311;
t134 = t163 * t180 - t165 * t179 + t312;
t133 = qJD(6) * t245 + t179 * t332 - t212 * t334 + t313;
t132 = qJD(6) * t243 - t180 * t332 + t212 * t333 + t311;
t131 = qJD(6) * t260 - t179 * t333 + t180 * t334 + t312;
t1 = t276 * ((t224 * t341 - t264 * t226 + t265 * t228) * t276 + (t223 * t341 - t225 * t264 + t227 * t265) * t275 + (t251 * t341 - t252 * t264 + t253 * t265) * t277) / 0.2e1 + t275 * ((-t224 * t339 - t226 * t262 + t228 * t263) * t276 + (-t223 * t339 - t262 * t225 + t263 * t227) * t275 + (-t251 * t339 - t252 * t262 + t253 * t263) * t277) / 0.2e1 + t258 * ((-t195 * t340 - t197 * t260 + t199 * t261) * t242 + (-t194 * t340 - t196 * t260 + t198 * t261) * t241 + (-t219 * t340 - t220 * t260 + t261 * t221) * t258) / 0.2e1 + t277 * ((t223 * t275 + t224 * t276 + t251 * t277) * t302 + ((t226 * t308 + t228 * t305) * t276 + (t225 * t308 + t227 * t305) * t275 + (t252 * t308 + t253 * t305) * t277) * t301) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(3) * (t166 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + t206 * ((t168 * t243 + t170 * t213 + t172 * t214) * t207 + (t243 * t167 + t169 * t213 + t171 * t214) * t206 + (t190 * t243 + t191 * t213 + t192 * t214) * t232) / 0.2e1 + t207 * ((t245 * t168 + t215 * t170 + t216 * t172) * t207 + (t167 * t245 + t169 * t215 + t171 * t216) * t206 + (t190 * t245 + t191 * t215 + t192 * t216) * t232) / 0.2e1 + m(2) * (t247 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + t232 * ((t168 * t260 + t170 * t239 + t172 * t240) * t207 + (t167 * t260 + t169 * t239 + t171 * t240) * t206 + (t260 * t190 + t239 * t191 + t240 * t192) * t232) / 0.2e1 + t241 * ((t195 * t262 - t197 * t243 + t199 * t244) * t242 + (t262 * t194 - t243 * t196 + t244 * t198) * t241 + (t219 * t262 - t220 * t243 + t221 * t244) * t258) / 0.2e1 + t242 * ((t264 * t195 - t245 * t197 + t246 * t199) * t242 + (t194 * t264 - t196 * t245 + t198 * t246) * t241 + (t219 * t264 - t220 * t245 + t221 * t246) * t258) / 0.2e1 + m(1) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + ((t208 * t354 + t209 * t353 + t243 * t355) * t212 + (t208 * t358 + t209 * t356 + t243 * t360) * t180 + (t359 * t208 + t357 * t209 + t361 * t243) * t179) * t179 / 0.2e1 + ((t210 * t354 + t211 * t353 + t245 * t355) * t212 + (t358 * t210 + t356 * t211 + t360 * t245) * t180 + (t359 * t210 + t357 * t211 + t245 * t361) * t179) * t180 / 0.2e1 + ((t354 * t235 + t353 * t236 + t355 * t260) * t212 + (t235 * t358 + t236 * t356 + t260 * t360) * t180 + (t359 * t235 + t357 * t236 + t260 * t361) * t179) * t212 / 0.2e1 + ((-t280 * t306 + t282 * t309 + Icges(1,4)) * V_base(5) + (-t281 * t306 + t283 * t309 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t280 * t309 + t282 * t306 + Icges(1,2)) * V_base(5) + (t281 * t309 + t283 * t306 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t306 + Icges(2,6) * t309) * V_base(5) + (Icges(2,5) * t309 - Icges(2,6) * t306) * V_base(4) + Icges(2,3) * t293 / 0.2e1) * t293;
T  = t1;
