% Calculate kinetic energy for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:28
% EndTime: 2019-03-09 11:32:32
% DurationCPUTime: 4.01s
% Computational Cost: add. (2063->350), mult. (4626->485), div. (0->0), fcn. (5334->10), ass. (0->158)
t363 = Icges(3,1) + Icges(4,2);
t362 = Icges(5,1) + Icges(6,2);
t361 = Icges(6,1) + Icges(5,3);
t360 = Icges(3,4) + Icges(4,6);
t359 = Icges(5,4) + Icges(6,6);
t358 = Icges(6,4) - Icges(5,5);
t357 = Icges(3,5) - Icges(4,4);
t356 = Icges(6,5) - Icges(5,6);
t355 = Icges(3,2) + Icges(4,3);
t354 = Icges(5,2) + Icges(6,3);
t353 = Icges(3,6) - Icges(4,5);
t352 = Icges(3,3) + Icges(4,1);
t292 = cos(pkin(6));
t297 = cos(qJ(2));
t298 = cos(qJ(1));
t317 = t297 * t298;
t294 = sin(qJ(2));
t295 = sin(qJ(1));
t320 = t294 * t295;
t258 = -t292 * t317 + t320;
t291 = sin(pkin(6));
t326 = sin(qJ(4));
t311 = t291 * t326;
t327 = cos(qJ(4));
t229 = t258 * t327 + t298 * t311;
t312 = t291 * t327;
t230 = -t258 * t326 + t298 * t312;
t318 = t295 * t297;
t319 = t294 * t298;
t259 = t292 * t319 + t318;
t351 = -t354 * t229 + t359 * t230 + t356 * t259;
t260 = t292 * t318 + t319;
t227 = -t260 * t327 + t295 * t311;
t228 = t260 * t326 + t295 * t312;
t261 = -t292 * t320 + t317;
t350 = t354 * t227 - t359 * t228 + t356 * t261;
t349 = t356 * t227 - t358 * t228 + t361 * t261;
t348 = -t356 * t229 + t358 * t230 + t361 * t259;
t347 = -t359 * t229 + t362 * t230 + t358 * t259;
t346 = t359 * t227 - t362 * t228 + t358 * t261;
t256 = t292 * t326 + t297 * t312;
t257 = t292 * t327 - t297 * t311;
t323 = t291 * t294;
t345 = t354 * t256 - t359 * t257 + t356 * t323;
t344 = t359 * t256 - t362 * t257 + t358 * t323;
t343 = t356 * t256 - t358 * t257 + t361 * t323;
t322 = t291 * t295;
t342 = t355 * t260 - t360 * t261 - t353 * t322;
t321 = t291 * t298;
t341 = t355 * t258 - t360 * t259 + t353 * t321;
t340 = -t360 * t260 + t363 * t261 + t357 * t322;
t339 = -t360 * t258 + t363 * t259 - t357 * t321;
t338 = -t353 * t260 + t357 * t261 + t352 * t322;
t337 = -t353 * t258 + t357 * t259 - t352 * t321;
t336 = t352 * t292 + (t357 * t294 + t353 * t297) * t291;
t335 = t353 * t292 + (t360 * t294 + t355 * t297) * t291;
t334 = t357 * t292 + (t363 * t294 + t360 * t297) * t291;
t325 = pkin(8) * t292;
t324 = Icges(2,4) * t295;
t316 = qJD(2) * t291;
t315 = V_base(5) * pkin(7) + V_base(1);
t271 = t295 * t316 + V_base(4);
t288 = V_base(6) + qJD(1);
t226 = qJD(4) * t261 + t271;
t272 = qJD(2) * t292 + t288;
t254 = qJD(4) * t323 + t272;
t270 = -t298 * t316 + V_base(5);
t265 = t295 * pkin(1) - pkin(8) * t321;
t310 = -t265 * t288 + V_base(5) * t325 + t315;
t266 = pkin(1) * t298 + pkin(8) * t322;
t309 = V_base(4) * t265 - t266 * V_base(5) + V_base(3);
t225 = qJD(4) * t259 + t270;
t262 = (pkin(2) * t294 - qJ(3) * t297) * t291;
t308 = qJD(3) * t260 + t270 * t262 + t310;
t307 = t288 * t266 + V_base(2) + (-pkin(7) - t325) * V_base(4);
t220 = pkin(2) * t261 + qJ(3) * t260;
t306 = qJD(3) * t258 + t272 * t220 + t307;
t219 = pkin(2) * t259 + qJ(3) * t258;
t305 = -qJD(3) * t291 * t297 + t271 * t219 + t309;
t237 = -pkin(3) * t321 + t259 * pkin(9);
t264 = pkin(3) * t292 + pkin(9) * t323;
t304 = t270 * t264 + (-t219 - t237) * t272 + t308;
t216 = pkin(4) * t257 + qJ(5) * t256;
t303 = qJD(5) * t227 + t225 * t216 + t304;
t236 = pkin(3) * t322 + pkin(9) * t261;
t302 = t272 * t236 + (-t262 - t264) * t271 + t306;
t301 = t271 * t237 + (-t220 - t236) * t270 + t305;
t179 = pkin(4) * t228 + qJ(5) * t227;
t300 = -qJD(5) * t229 + t254 * t179 + t302;
t180 = -pkin(4) * t230 - qJ(5) * t229;
t299 = qJD(5) * t256 + t226 * t180 + t301;
t296 = cos(qJ(6));
t293 = sin(qJ(6));
t289 = Icges(2,4) * t298;
t280 = rSges(2,1) * t298 - t295 * rSges(2,2);
t279 = t295 * rSges(2,1) + rSges(2,2) * t298;
t278 = Icges(2,1) * t298 - t324;
t277 = Icges(2,1) * t295 + t289;
t276 = -Icges(2,2) * t295 + t289;
t275 = Icges(2,2) * t298 + t324;
t269 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t268 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t267 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t246 = rSges(4,1) * t292 + (-rSges(4,2) * t294 - rSges(4,3) * t297) * t291;
t245 = rSges(3,3) * t292 + (rSges(3,1) * t294 + rSges(3,2) * t297) * t291;
t235 = pkin(5) * t323 + pkin(10) * t257;
t234 = V_base(5) * rSges(2,3) - t279 * t288 + t315;
t233 = t280 * t288 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t231 = t279 * V_base(4) - t280 * V_base(5) + V_base(3);
t224 = t256 * t293 + t296 * t323;
t223 = t256 * t296 - t293 * t323;
t217 = qJD(6) * t257 + t254;
t214 = rSges(3,1) * t261 - rSges(3,2) * t260 + rSges(3,3) * t322;
t213 = t259 * rSges(3,1) - t258 * rSges(3,2) - rSges(3,3) * t321;
t212 = -rSges(4,1) * t321 - t259 * rSges(4,2) + t258 * rSges(4,3);
t211 = rSges(4,1) * t322 - rSges(4,2) * t261 + rSges(4,3) * t260;
t198 = rSges(5,1) * t257 - rSges(5,2) * t256 + rSges(5,3) * t323;
t197 = rSges(6,1) * t323 - rSges(6,2) * t257 + rSges(6,3) * t256;
t188 = pkin(5) * t261 + pkin(10) * t228;
t187 = pkin(5) * t259 - pkin(10) * t230;
t186 = t227 * t293 + t261 * t296;
t185 = t227 * t296 - t261 * t293;
t184 = -t229 * t293 + t259 * t296;
t183 = -t229 * t296 - t259 * t293;
t182 = qJD(6) * t228 + t226;
t181 = -qJD(6) * t230 + t225;
t177 = rSges(6,1) * t261 - rSges(6,2) * t228 + rSges(6,3) * t227;
t176 = rSges(6,1) * t259 + rSges(6,2) * t230 - rSges(6,3) * t229;
t175 = -rSges(5,1) * t230 + rSges(5,2) * t229 + rSges(5,3) * t259;
t174 = rSges(5,1) * t228 - rSges(5,2) * t227 + rSges(5,3) * t261;
t160 = rSges(7,1) * t224 + rSges(7,2) * t223 + rSges(7,3) * t257;
t159 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t257;
t158 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t257;
t157 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t257;
t155 = -t213 * t272 + t245 * t270 + t310;
t154 = t214 * t272 - t245 * t271 + t307;
t153 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t228;
t152 = rSges(7,1) * t184 + rSges(7,2) * t183 - rSges(7,3) * t230;
t151 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t228;
t150 = Icges(7,1) * t184 + Icges(7,4) * t183 - Icges(7,5) * t230;
t149 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t228;
t148 = Icges(7,4) * t184 + Icges(7,2) * t183 - Icges(7,6) * t230;
t147 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t228;
t146 = Icges(7,5) * t184 + Icges(7,6) * t183 - Icges(7,3) * t230;
t145 = t213 * t271 - t214 * t270 + t309;
t144 = t246 * t270 + (-t212 - t219) * t272 + t308;
t143 = t211 * t272 + (-t246 - t262) * t271 + t306;
t142 = t212 * t271 + (-t211 - t220) * t270 + t305;
t141 = -t175 * t254 + t198 * t225 + t304;
t140 = t174 * t254 - t198 * t226 + t302;
t139 = -t174 * t225 + t175 * t226 + t301;
t138 = t197 * t225 + (-t176 - t180) * t254 + t303;
t137 = t177 * t254 + (-t197 - t216) * t226 + t300;
t136 = t176 * t226 + (-t177 - t179) * t225 + t299;
t135 = t303 - t152 * t217 + t160 * t181 + t225 * t235 + (-t180 - t187) * t254;
t134 = t153 * t217 - t160 * t182 + t188 * t254 + (-t216 - t235) * t226 + t300;
t133 = t152 * t182 - t153 * t181 + t187 * t226 + (-t179 - t188) * t225 + t299;
t1 = m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t181 * ((-t147 * t230 + t149 * t183 + t151 * t184) * t182 + (-t230 * t146 + t183 * t148 + t184 * t150) * t181 + (-t157 * t230 + t158 * t183 + t159 * t184) * t217) / 0.2e1 + m(2) * (t231 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + t182 * ((t228 * t147 + t185 * t149 + t186 * t151) * t182 + (t146 * t228 + t148 * t185 + t150 * t186) * t181 + (t157 * t228 + t158 * t185 + t159 * t186) * t217) / 0.2e1 + m(3) * (t145 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(1) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + t217 * ((t147 * t257 + t149 * t223 + t151 * t224) * t182 + (t146 * t257 + t148 * t223 + t150 * t224) * t181 + (t157 * t257 + t158 * t223 + t159 * t224) * t217) / 0.2e1 + ((-t229 * t345 + t230 * t344 + t259 * t343) * t254 + (-t229 * t350 + t230 * t346 + t259 * t349) * t226 + (-t351 * t229 + t347 * t230 + t348 * t259) * t225) * t225 / 0.2e1 + ((t227 * t345 - t228 * t344 + t261 * t343) * t254 + (t350 * t227 - t346 * t228 + t349 * t261) * t226 + (t227 * t351 - t347 * t228 + t348 * t261) * t225) * t226 / 0.2e1 + ((t345 * t256 - t344 * t257 + t343 * t323) * t254 + (t256 * t350 - t257 * t346 + t323 * t349) * t226 + (t256 * t351 - t347 * t257 + t348 * t323) * t225) * t254 / 0.2e1 + ((-t258 * t335 + t259 * t334 - t321 * t336) * t272 + (t258 * t342 + t259 * t340 - t321 * t338) * t271 + (t341 * t258 + t339 * t259 - t337 * t321) * t270) * t270 / 0.2e1 + ((-t260 * t335 + t261 * t334 + t322 * t336) * t272 + (t342 * t260 + t340 * t261 + t338 * t322) * t271 + (t260 * t341 + t261 * t339 + t322 * t337) * t270) * t271 / 0.2e1 + ((t270 * t337 + t271 * t338 + t272 * t336) * t292 + ((t294 * t334 + t297 * t335) * t272 + (t294 * t340 - t297 * t342) * t271 + (t294 * t339 - t297 * t341) * t270) * t291) * t272 / 0.2e1 + ((-t295 * t275 + t277 * t298 + Icges(1,4)) * V_base(5) + (-t295 * t276 + t278 * t298 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t275 * t298 + t295 * t277 + Icges(1,2)) * V_base(5) + (t276 * t298 + t295 * t278 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t295 + Icges(2,6) * t298) * V_base(5) + (Icges(2,5) * t298 - Icges(2,6) * t295) * V_base(4) + Icges(2,3) * t288 / 0.2e1) * t288;
T  = t1;
