% Calculate kinetic energy for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:03
% EndTime: 2019-03-08 21:08:08
% DurationCPUTime: 4.65s
% Computational Cost: add. (2215->344), mult. (5148->481), div. (0->0), fcn. (6035->10), ass. (0->153)
t350 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t349 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t348 = -Icges(4,4) - Icges(6,4) + Icges(5,5);
t347 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t346 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t345 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t291 = sin(pkin(10));
t293 = cos(pkin(10));
t297 = cos(qJ(2));
t295 = sin(qJ(2));
t328 = cos(pkin(6));
t311 = t295 * t328;
t256 = t291 * t297 + t293 * t311;
t292 = sin(pkin(6));
t330 = cos(qJ(3));
t314 = t292 * t330;
t329 = sin(qJ(3));
t236 = t256 * t329 + t293 * t314;
t313 = t292 * t329;
t237 = t256 * t330 - t293 * t313;
t310 = t297 * t328;
t255 = t291 * t295 - t293 * t310;
t342 = t236 * t346 - t237 * t347 - t255 * t345;
t258 = -t291 * t311 + t293 * t297;
t238 = t258 * t329 - t291 * t314;
t239 = t258 * t330 + t291 * t313;
t257 = t291 * t310 + t293 * t295;
t341 = t238 * t346 - t239 * t347 - t257 * t345;
t340 = t236 * t349 + t237 * t348 - t255 * t346;
t339 = t238 * t349 + t239 * t348 - t257 * t346;
t338 = t348 * t236 + t237 * t350 + t347 * t255;
t337 = t348 * t238 + t239 * t350 + t347 * t257;
t262 = t295 * t313 - t328 * t330;
t263 = t295 * t314 + t328 * t329;
t324 = t292 * t297;
t336 = t262 * t346 - t263 * t347 + t324 * t345;
t335 = t262 * t349 + t263 * t348 + t324 * t346;
t334 = t348 * t262 + t263 * t350 - t347 * t324;
t327 = Icges(2,4) * t291;
t326 = t291 * t292;
t325 = t292 * t293;
t191 = pkin(3) * t237 + qJ(4) * t236;
t202 = pkin(4) * t237 - qJ(5) * t255;
t323 = -t191 - t202;
t192 = pkin(3) * t239 + qJ(4) * t238;
t203 = pkin(4) * t239 - qJ(5) * t257;
t322 = -t192 - t203;
t228 = pkin(3) * t263 + qJ(4) * t262;
t243 = t263 * pkin(4) + qJ(5) * t324;
t321 = -t228 - t243;
t320 = qJD(2) * t292;
t319 = V_base(5) * qJ(1) + V_base(1);
t315 = qJD(1) + V_base(3);
t312 = t328 * pkin(7);
t271 = t291 * t320 + V_base(4);
t283 = qJD(2) * t328 + V_base(6);
t235 = qJD(3) * t257 + t271;
t270 = -t293 * t320 + V_base(5);
t234 = qJD(3) * t255 + t270;
t259 = -qJD(3) * t324 + t283;
t265 = pkin(1) * t291 - pkin(7) * t325;
t309 = -t265 * V_base(6) + V_base(5) * t312 + t319;
t266 = pkin(1) * t293 + pkin(7) * t326;
t308 = V_base(4) * t265 - t266 * V_base(5) + t315;
t307 = V_base(6) * t266 + V_base(2) + (-t312 - qJ(1)) * V_base(4);
t226 = pkin(2) * t256 + pkin(8) * t255;
t264 = (pkin(2) * t295 - pkin(8) * t297) * t292;
t306 = -t226 * t283 + t270 * t264 + t309;
t227 = pkin(2) * t258 + pkin(8) * t257;
t305 = t271 * t226 - t227 * t270 + t308;
t304 = qJD(4) * t238 + t234 * t228 + t306;
t303 = qJD(4) * t262 + t235 * t191 + t305;
t302 = t283 * t227 - t271 * t264 + t307;
t301 = -qJD(5) * t257 + t234 * t243 + t304;
t300 = qJD(5) * t324 + t235 * t202 + t303;
t299 = qJD(4) * t236 + t259 * t192 + t302;
t298 = -qJD(5) * t255 + t259 * t203 + t299;
t296 = cos(qJ(6));
t294 = sin(qJ(6));
t289 = Icges(2,4) * t293;
t279 = rSges(2,1) * t293 - rSges(2,2) * t291;
t278 = rSges(2,1) * t291 + rSges(2,2) * t293;
t277 = Icges(2,1) * t293 - t327;
t276 = Icges(2,1) * t291 + t289;
t275 = -Icges(2,2) * t291 + t289;
t274 = Icges(2,2) * t293 + t327;
t269 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t268 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t267 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t249 = t328 * rSges(3,3) + (rSges(3,1) * t295 + rSges(3,2) * t297) * t292;
t248 = Icges(3,5) * t328 + (Icges(3,1) * t295 + Icges(3,4) * t297) * t292;
t247 = Icges(3,6) * t328 + (Icges(3,4) * t295 + Icges(3,2) * t297) * t292;
t246 = Icges(3,3) * t328 + (Icges(3,5) * t295 + Icges(3,6) * t297) * t292;
t245 = V_base(5) * rSges(2,3) - t278 * V_base(6) + t319;
t244 = t279 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t241 = t262 * t296 + t294 * t324;
t240 = -t262 * t294 + t296 * t324;
t233 = t278 * V_base(4) - t279 * V_base(5) + t315;
t230 = qJD(6) * t263 + t259;
t229 = pkin(5) * t262 + pkin(9) * t263;
t225 = t263 * rSges(4,1) - t262 * rSges(4,2) - rSges(4,3) * t324;
t224 = t263 * rSges(5,1) - rSges(5,2) * t324 + t262 * rSges(5,3);
t223 = t262 * rSges(6,1) - t263 * rSges(6,2) + rSges(6,3) * t324;
t213 = rSges(3,1) * t258 - rSges(3,2) * t257 + rSges(3,3) * t326;
t212 = rSges(3,1) * t256 - rSges(3,2) * t255 - rSges(3,3) * t325;
t211 = Icges(3,1) * t258 - Icges(3,4) * t257 + Icges(3,5) * t326;
t210 = Icges(3,1) * t256 - Icges(3,4) * t255 - Icges(3,5) * t325;
t209 = Icges(3,4) * t258 - Icges(3,2) * t257 + Icges(3,6) * t326;
t208 = Icges(3,4) * t256 - Icges(3,2) * t255 - Icges(3,6) * t325;
t207 = Icges(3,5) * t258 - Icges(3,6) * t257 + Icges(3,3) * t326;
t206 = Icges(3,5) * t256 - Icges(3,6) * t255 - Icges(3,3) * t325;
t201 = t238 * t296 - t257 * t294;
t200 = -t238 * t294 - t257 * t296;
t199 = t236 * t296 - t255 * t294;
t198 = -t236 * t294 - t255 * t296;
t196 = qJD(6) * t239 + t235;
t195 = qJD(6) * t237 + t234;
t194 = pkin(5) * t238 + pkin(9) * t239;
t193 = pkin(5) * t236 + pkin(9) * t237;
t188 = rSges(7,1) * t241 + rSges(7,2) * t240 + rSges(7,3) * t263;
t186 = Icges(7,1) * t241 + Icges(7,4) * t240 + Icges(7,5) * t263;
t185 = Icges(7,4) * t241 + Icges(7,2) * t240 + Icges(7,6) * t263;
t184 = Icges(7,5) * t241 + Icges(7,6) * t240 + Icges(7,3) * t263;
t182 = rSges(4,1) * t239 - rSges(4,2) * t238 + rSges(4,3) * t257;
t181 = rSges(5,1) * t239 + rSges(5,2) * t257 + rSges(5,3) * t238;
t180 = rSges(6,1) * t238 - rSges(6,2) * t239 - rSges(6,3) * t257;
t179 = rSges(4,1) * t237 - rSges(4,2) * t236 + rSges(4,3) * t255;
t178 = rSges(5,1) * t237 + rSges(5,2) * t255 + rSges(5,3) * t236;
t177 = rSges(6,1) * t236 - rSges(6,2) * t237 - rSges(6,3) * t255;
t157 = -t212 * t283 + t249 * t270 + t309;
t156 = t283 * t213 - t271 * t249 + t307;
t155 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t239;
t154 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t237;
t153 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t239;
t152 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t237;
t151 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t239;
t150 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t237;
t149 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t239;
t148 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t237;
t147 = t212 * t271 - t213 * t270 + t308;
t146 = -t179 * t259 + t225 * t234 + t306;
t145 = t259 * t182 - t235 * t225 + t302;
t144 = t179 * t235 - t182 * t234 + t305;
t143 = t224 * t234 + (-t178 - t191) * t259 + t304;
t142 = t259 * t181 + (-t224 - t228) * t235 + t299;
t141 = t178 * t235 + (-t181 - t192) * t234 + t303;
t140 = t223 * t234 + (-t177 + t323) * t259 + t301;
t139 = t259 * t180 + (-t223 + t321) * t235 + t298;
t138 = t177 * t235 + (-t180 + t322) * t234 + t300;
t137 = t301 + (-t193 + t323) * t259 - t154 * t230 + t188 * t195 + t229 * t234;
t136 = t230 * t155 - t196 * t188 + t259 * t194 + (-t229 + t321) * t235 + t298;
t135 = t300 + (-t194 + t322) * t234 + t154 * t196 - t155 * t195 + t193 * t235;
t1 = t283 * (((t209 * t297 + t211 * t295) * t271 + (t208 * t297 + t210 * t295) * t270 + (t247 * t297 + t248 * t295) * t283) * t292 + (t206 * t270 + t207 * t271 + t246 * t283) * t328) / 0.2e1 + t270 * ((-t207 * t325 - t209 * t255 + t211 * t256) * t271 + (-t206 * t325 - t255 * t208 + t256 * t210) * t270 + (-t246 * t325 - t247 * t255 + t248 * t256) * t283) / 0.2e1 + t271 * ((t207 * t326 - t257 * t209 + t258 * t211) * t271 + (t206 * t326 - t208 * t257 + t210 * t258) * t270 + (t246 * t326 - t247 * t257 + t248 * t258) * t283) / 0.2e1 + m(3) * (t147 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t195 * ((t149 * t237 + t151 * t198 + t153 * t199) * t196 + (t237 * t148 + t198 * t150 + t199 * t152) * t195 + (t184 * t237 + t185 * t198 + t186 * t199) * t230) / 0.2e1 + t196 * ((t239 * t149 + t200 * t151 + t201 * t153) * t196 + (t148 * t239 + t150 * t200 + t152 * t201) * t195 + (t184 * t239 + t185 * t200 + t186 * t201) * t230) / 0.2e1 + m(2) * (t233 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + t230 * ((t149 * t263 + t151 * t240 + t153 * t241) * t196 + (t148 * t263 + t150 * t240 + t152 * t241) * t195 + (t263 * t184 + t240 * t185 + t241 * t186) * t230) / 0.2e1 + m(1) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + ((-t274 * t291 + t276 * t293 + Icges(1,4)) * V_base(5) + (-t291 * t275 + t293 * t277 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t293 * t274 + t291 * t276 + Icges(1,2)) * V_base(5) + (t275 * t293 + t277 * t291 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t236 * t335 + t237 * t334 - t255 * t336) * t259 + (t236 * t339 + t237 * t337 - t255 * t341) * t235 + (t340 * t236 + t338 * t237 - t342 * t255) * t234) * t234 / 0.2e1 + ((t238 * t335 + t239 * t334 - t257 * t336) * t259 + (t339 * t238 + t337 * t239 - t341 * t257) * t235 + (t238 * t340 + t239 * t338 - t257 * t342) * t234) * t235 / 0.2e1 + ((t335 * t262 + t334 * t263 + t336 * t324) * t259 + (t262 * t339 + t263 * t337 + t324 * t341) * t235 + (t262 * t340 + t263 * t338 + t324 * t342) * t234) * t259 / 0.2e1 + ((Icges(2,5) * t291 + Icges(2,6) * t293 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t293 - Icges(2,6) * t291 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
