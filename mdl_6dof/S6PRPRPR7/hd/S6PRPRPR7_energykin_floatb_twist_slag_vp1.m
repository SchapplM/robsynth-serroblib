% Calculate kinetic energy for
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:06
% EndTime: 2019-03-08 19:51:10
% DurationCPUTime: 4.10s
% Computational Cost: add. (2003->350), mult. (4626->483), div. (0->0), fcn. (5334->10), ass. (0->156)
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
t289 = sin(pkin(10));
t291 = cos(pkin(10));
t294 = sin(qJ(2));
t292 = cos(pkin(6));
t296 = cos(qJ(2));
t317 = t292 * t296;
t253 = t289 * t294 - t291 * t317;
t290 = sin(pkin(6));
t324 = sin(qJ(4));
t309 = t290 * t324;
t325 = cos(qJ(4));
t227 = t253 * t325 + t291 * t309;
t310 = t290 * t325;
t228 = -t253 * t324 + t291 * t310;
t318 = t292 * t294;
t254 = t289 * t296 + t291 * t318;
t351 = -t227 * t354 + t228 * t359 + t254 * t356;
t255 = t289 * t317 + t291 * t294;
t225 = -t255 * t325 + t289 * t309;
t226 = t255 * t324 + t289 * t310;
t256 = -t289 * t318 + t291 * t296;
t350 = t225 * t354 - t226 * t359 + t256 * t356;
t349 = t225 * t356 - t226 * t358 + t256 * t361;
t348 = -t227 * t356 + t228 * t358 + t254 * t361;
t347 = -t227 * t359 + t228 * t362 + t254 * t358;
t346 = t225 * t359 - t226 * t362 + t256 * t358;
t321 = t289 * t290;
t345 = t255 * t355 - t256 * t360 - t321 * t353;
t320 = t290 * t291;
t344 = t253 * t355 - t254 * t360 + t320 * t353;
t343 = -t255 * t360 + t256 * t363 + t321 * t357;
t342 = -t253 * t360 + t254 * t363 - t320 * t357;
t341 = -t255 * t353 + t256 * t357 + t321 * t352;
t340 = -t253 * t353 + t254 * t357 - t320 * t352;
t260 = t292 * t324 + t296 * t310;
t261 = t292 * t325 - t296 * t309;
t319 = t290 * t294;
t339 = t260 * t354 - t261 * t359 + t319 * t356;
t338 = t260 * t359 - t261 * t362 + t319 * t358;
t337 = t260 * t356 - t261 * t358 + t319 * t361;
t336 = t352 * t292 + (t294 * t357 + t296 * t353) * t290;
t335 = t353 * t292 + (t294 * t360 + t296 * t355) * t290;
t334 = t357 * t292 + (t294 * t363 + t296 * t360) * t290;
t323 = pkin(7) * t292;
t322 = Icges(2,4) * t289;
t316 = qJD(2) * t290;
t315 = V_base(5) * qJ(1) + V_base(1);
t311 = qJD(1) + V_base(3);
t270 = t289 * t316 + V_base(4);
t282 = qJD(2) * t292 + V_base(6);
t224 = qJD(4) * t256 + t270;
t257 = qJD(4) * t319 + t282;
t269 = -t291 * t316 + V_base(5);
t223 = qJD(4) * t254 + t269;
t263 = pkin(1) * t289 - pkin(7) * t320;
t308 = -t263 * V_base(6) + t323 * V_base(5) + t315;
t264 = pkin(1) * t291 + pkin(7) * t321;
t307 = t263 * V_base(4) - t264 * V_base(5) + t311;
t306 = V_base(6) * t264 + V_base(2) + (-qJ(1) - t323) * V_base(4);
t262 = (pkin(2) * t294 - qJ(3) * t296) * t290;
t305 = qJD(3) * t255 + t262 * t269 + t308;
t216 = pkin(2) * t256 + qJ(3) * t255;
t304 = qJD(3) * t253 + t216 * t282 + t306;
t215 = pkin(2) * t254 + qJ(3) * t253;
t303 = -qJD(3) * t290 * t296 + t215 * t270 + t307;
t233 = -pkin(3) * t320 + pkin(8) * t254;
t265 = pkin(3) * t292 + pkin(8) * t319;
t302 = t269 * t265 + (-t215 - t233) * t282 + t305;
t232 = pkin(3) * t321 + pkin(8) * t256;
t301 = t282 * t232 + (-t262 - t265) * t270 + t304;
t217 = pkin(4) * t261 + qJ(5) * t260;
t300 = qJD(5) * t225 + t217 * t223 + t302;
t299 = t270 * t233 + (-t216 - t232) * t269 + t303;
t178 = pkin(4) * t226 + qJ(5) * t225;
t298 = -qJD(5) * t227 + t178 * t257 + t301;
t179 = -pkin(4) * t228 - qJ(5) * t227;
t297 = qJD(5) * t260 + t179 * t224 + t299;
t295 = cos(qJ(6));
t293 = sin(qJ(6));
t287 = Icges(2,4) * t291;
t278 = rSges(2,1) * t291 - rSges(2,2) * t289;
t277 = rSges(2,1) * t289 + rSges(2,2) * t291;
t276 = Icges(2,1) * t291 - t322;
t275 = Icges(2,1) * t289 + t287;
t274 = -Icges(2,2) * t289 + t287;
t273 = Icges(2,2) * t291 + t322;
t268 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t267 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t266 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t245 = t292 * rSges(4,1) + (-rSges(4,2) * t294 - rSges(4,3) * t296) * t290;
t244 = t292 * rSges(3,3) + (rSges(3,1) * t294 + rSges(3,2) * t296) * t290;
t236 = V_base(5) * rSges(2,3) - t277 * V_base(6) + t315;
t235 = t278 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t234 = pkin(5) * t319 + pkin(9) * t261;
t230 = t260 * t293 + t295 * t319;
t229 = t260 * t295 - t293 * t319;
t222 = t277 * V_base(4) - t278 * V_base(5) + t311;
t219 = qJD(6) * t261 + t257;
t213 = rSges(5,1) * t261 - rSges(5,2) * t260 + rSges(5,3) * t319;
t212 = rSges(6,1) * t319 - rSges(6,2) * t261 + rSges(6,3) * t260;
t205 = rSges(3,1) * t256 - rSges(3,2) * t255 + rSges(3,3) * t321;
t204 = rSges(3,1) * t254 - rSges(3,2) * t253 - rSges(3,3) * t320;
t203 = -rSges(4,1) * t320 - rSges(4,2) * t254 + rSges(4,3) * t253;
t202 = rSges(4,1) * t321 - rSges(4,2) * t256 + rSges(4,3) * t255;
t187 = pkin(5) * t256 + pkin(9) * t226;
t186 = pkin(5) * t254 - pkin(9) * t228;
t185 = t225 * t293 + t256 * t295;
t184 = t225 * t295 - t256 * t293;
t183 = -t227 * t293 + t254 * t295;
t182 = -t227 * t295 - t254 * t293;
t181 = qJD(6) * t226 + t224;
t180 = -qJD(6) * t228 + t223;
t176 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t261;
t174 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t261;
t173 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t261;
t172 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t261;
t171 = rSges(6,1) * t256 - rSges(6,2) * t226 + rSges(6,3) * t225;
t170 = rSges(6,1) * t254 + rSges(6,2) * t228 - rSges(6,3) * t227;
t169 = -rSges(5,1) * t228 + rSges(5,2) * t227 + rSges(5,3) * t254;
t168 = rSges(5,1) * t226 - rSges(5,2) * t225 + rSges(5,3) * t256;
t154 = -t204 * t282 + t244 * t269 + t308;
t153 = t205 * t282 - t244 * t270 + t306;
t152 = rSges(7,1) * t185 + rSges(7,2) * t184 + rSges(7,3) * t226;
t151 = rSges(7,1) * t183 + rSges(7,2) * t182 - rSges(7,3) * t228;
t150 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t226;
t149 = Icges(7,1) * t183 + Icges(7,4) * t182 - Icges(7,5) * t228;
t148 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t226;
t147 = Icges(7,4) * t183 + Icges(7,2) * t182 - Icges(7,6) * t228;
t146 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t226;
t145 = Icges(7,5) * t183 + Icges(7,6) * t182 - Icges(7,3) * t228;
t144 = t204 * t270 - t205 * t269 + t307;
t143 = t245 * t269 + (-t203 - t215) * t282 + t305;
t142 = t202 * t282 + (-t245 - t262) * t270 + t304;
t141 = t270 * t203 + (-t202 - t216) * t269 + t303;
t140 = -t169 * t257 + t213 * t223 + t302;
t139 = t168 * t257 - t213 * t224 + t301;
t138 = -t168 * t223 + t169 * t224 + t299;
t137 = t212 * t223 + (-t170 - t179) * t257 + t300;
t136 = t171 * t257 + (-t212 - t217) * t224 + t298;
t135 = t224 * t170 + (-t171 - t178) * t223 + t297;
t134 = t300 + (-t179 - t186) * t257 - t151 * t219 + t176 * t180 + t223 * t234;
t133 = t152 * t219 - t176 * t181 + t187 * t257 + (-t217 - t234) * t224 + t298;
t132 = t297 - t180 * t152 + t181 * t151 + t224 * t186 + (-t178 - t187) * t223;
t1 = t219 * ((t146 * t261 + t148 * t229 + t150 * t230) * t181 + (t145 * t261 + t147 * t229 + t149 * t230) * t180 + (t172 * t261 + t173 * t229 + t174 * t230) * t219) / 0.2e1 + m(4) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(1) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(5) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(6) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(2) * (t222 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + t180 * ((-t146 * t228 + t148 * t182 + t150 * t183) * t181 + (-t228 * t145 + t182 * t147 + t183 * t149) * t180 + (-t172 * t228 + t173 * t182 + t174 * t183) * t219) / 0.2e1 + t181 * ((t226 * t146 + t184 * t148 + t185 * t150) * t181 + (t145 * t226 + t147 * t184 + t149 * t185) * t180 + (t172 * t226 + t173 * t184 + t174 * t185) * t219) / 0.2e1 + m(3) * (t144 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + ((-t227 * t339 + t228 * t338 + t254 * t337) * t257 + (-t227 * t350 + t228 * t346 + t254 * t349) * t224 + (-t351 * t227 + t347 * t228 + t348 * t254) * t223) * t223 / 0.2e1 + ((t225 * t339 - t226 * t338 + t256 * t337) * t257 + (t350 * t225 - t346 * t226 + t349 * t256) * t224 + (t225 * t351 - t226 * t347 + t256 * t348) * t223) * t224 / 0.2e1 + ((t260 * t339 - t338 * t261 + t337 * t319) * t257 + (t260 * t350 - t261 * t346 + t319 * t349) * t224 + (t260 * t351 - t261 * t347 + t319 * t348) * t223) * t257 / 0.2e1 + ((-t253 * t335 + t254 * t334 - t320 * t336) * t282 + (t253 * t345 + t254 * t343 - t320 * t341) * t270 + (t344 * t253 + t342 * t254 - t340 * t320) * t269) * t269 / 0.2e1 + ((-t255 * t335 + t256 * t334 + t321 * t336) * t282 + (t345 * t255 + t343 * t256 + t341 * t321) * t270 + (t255 * t344 + t256 * t342 + t321 * t340) * t269) * t270 / 0.2e1 + ((t269 * t340 + t270 * t341 + t282 * t336) * t292 + ((t294 * t334 + t296 * t335) * t282 + (t294 * t343 - t296 * t345) * t270 + (t294 * t342 - t296 * t344) * t269) * t290) * t282 / 0.2e1 + ((-t273 * t289 + t275 * t291 + Icges(1,4)) * V_base(5) + (-t274 * t289 + t276 * t291 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t273 * t291 + t275 * t289 + Icges(1,2)) * V_base(5) + (t274 * t291 + t276 * t289 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t291 - Icges(2,6) * t289 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t289 + Icges(2,6) * t291 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
