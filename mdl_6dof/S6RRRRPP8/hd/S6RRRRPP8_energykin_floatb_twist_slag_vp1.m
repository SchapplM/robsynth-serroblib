% Calculate kinetic energy for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:13
% EndTime: 2019-03-09 21:30:17
% DurationCPUTime: 4.54s
% Computational Cost: add. (2845->334), mult. (6592->476), div. (0->0), fcn. (8001->10), ass. (0->148)
t352 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t351 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t350 = Icges(6,4) + Icges(5,5) - Icges(7,5);
t349 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t348 = Icges(6,2) + Icges(5,3) + Icges(7,3);
t347 = -Icges(5,6) + Icges(6,6) - Icges(7,6);
t346 = rSges(7,1) + pkin(5);
t345 = -rSges(7,3) - qJ(6);
t299 = sin(qJ(2));
t300 = sin(qJ(1));
t301 = cos(qJ(2));
t302 = cos(qJ(1));
t330 = cos(pkin(6));
t315 = t302 * t330;
t266 = t299 * t315 + t300 * t301;
t298 = sin(qJ(3));
t296 = sin(pkin(6));
t326 = t296 * t302;
t332 = cos(qJ(3));
t246 = t266 * t332 - t298 * t326;
t265 = t299 * t300 - t301 * t315;
t297 = sin(qJ(4));
t331 = cos(qJ(4));
t217 = t246 * t297 - t265 * t331;
t218 = t246 * t331 + t265 * t297;
t318 = t296 * t332;
t245 = t266 * t298 + t302 * t318;
t344 = t347 * t217 + t350 * t218 + t348 * t245;
t316 = t300 * t330;
t268 = -t299 * t316 + t302 * t301;
t328 = t296 * t300;
t248 = t268 * t332 + t298 * t328;
t267 = t302 * t299 + t301 * t316;
t219 = t248 * t297 - t267 * t331;
t220 = t248 * t331 + t267 * t297;
t247 = t268 * t298 - t300 * t318;
t343 = t347 * t219 + t350 * t220 + t348 * t247;
t342 = t349 * t217 + t351 * t218 + t347 * t245;
t341 = t349 * t219 + t351 * t220 + t347 * t247;
t340 = t351 * t217 + t352 * t218 + t350 * t245;
t339 = t351 * t219 + t352 * t220 + t350 * t247;
t264 = t298 * t330 + t299 * t318;
t327 = t296 * t301;
t241 = t264 * t297 + t327 * t331;
t242 = t264 * t331 - t297 * t327;
t263 = t296 * t298 * t299 - t330 * t332;
t338 = t347 * t241 + t350 * t242 + t348 * t263;
t337 = t349 * t241 + t351 * t242 + t347 * t263;
t336 = t351 * t241 + t352 * t242 + t350 * t263;
t329 = Icges(2,4) * t300;
t325 = rSges(7,2) * t217 + t218 * t346 + t345 * t245;
t324 = rSges(7,2) * t219 + t220 * t346 + t345 * t247;
t323 = rSges(7,2) * t241 + t242 * t346 + t345 * t263;
t322 = qJD(2) * t296;
t321 = V_base(5) * pkin(7) + V_base(1);
t317 = t330 * pkin(8);
t277 = t300 * t322 + V_base(4);
t293 = V_base(6) + qJD(1);
t244 = qJD(3) * t267 + t277;
t278 = qJD(2) * t330 + t293;
t276 = -t302 * t322 + V_base(5);
t271 = t300 * pkin(1) - pkin(8) * t326;
t314 = -t271 * t293 + V_base(5) * t317 + t321;
t272 = pkin(1) * t302 + pkin(8) * t328;
t313 = V_base(4) * t271 - t272 * V_base(5) + V_base(3);
t243 = qJD(3) * t265 + t276;
t261 = -qJD(3) * t327 + t278;
t238 = pkin(2) * t266 + pkin(9) * t265;
t270 = (pkin(2) * t299 - pkin(9) * t301) * t296;
t312 = -t238 * t278 + t276 * t270 + t314;
t239 = pkin(2) * t268 + pkin(9) * t267;
t311 = t277 * t238 - t239 * t276 + t313;
t310 = t293 * t272 + V_base(2) + (-t317 - pkin(7)) * V_base(4);
t211 = pkin(3) * t246 + pkin(10) * t245;
t236 = pkin(3) * t264 + pkin(10) * t263;
t309 = -t211 * t261 + t243 * t236 + t312;
t212 = pkin(3) * t248 + pkin(10) * t247;
t308 = t244 * t211 - t212 * t243 + t311;
t307 = t278 * t239 - t277 * t270 + t310;
t210 = pkin(4) * t242 + qJ(5) * t241;
t213 = qJD(4) * t245 + t243;
t306 = qJD(5) * t219 + t213 * t210 + t309;
t183 = pkin(4) * t218 + qJ(5) * t217;
t214 = qJD(4) * t247 + t244;
t305 = qJD(5) * t241 + t214 * t183 + t308;
t304 = t261 * t212 - t244 * t236 + t307;
t184 = pkin(4) * t220 + qJ(5) * t219;
t237 = qJD(4) * t263 + t261;
t303 = qJD(5) * t217 + t237 * t184 + t304;
t294 = Icges(2,4) * t302;
t286 = rSges(2,1) * t302 - t300 * rSges(2,2);
t285 = t300 * rSges(2,1) + rSges(2,2) * t302;
t284 = Icges(2,1) * t302 - t329;
t283 = Icges(2,1) * t300 + t294;
t282 = -Icges(2,2) * t300 + t294;
t281 = Icges(2,2) * t302 + t329;
t275 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t274 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t273 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t256 = t330 * rSges(3,3) + (rSges(3,1) * t299 + rSges(3,2) * t301) * t296;
t255 = Icges(3,5) * t330 + (Icges(3,1) * t299 + Icges(3,4) * t301) * t296;
t254 = Icges(3,6) * t330 + (Icges(3,4) * t299 + Icges(3,2) * t301) * t296;
t253 = Icges(3,3) * t330 + (Icges(3,5) * t299 + Icges(3,6) * t301) * t296;
t252 = V_base(5) * rSges(2,3) - t285 * t293 + t321;
t251 = t286 * t293 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t249 = t285 * V_base(4) - t286 * V_base(5) + V_base(3);
t235 = rSges(3,1) * t268 - rSges(3,2) * t267 + rSges(3,3) * t328;
t234 = t266 * rSges(3,1) - t265 * rSges(3,2) - rSges(3,3) * t326;
t233 = Icges(3,1) * t268 - Icges(3,4) * t267 + Icges(3,5) * t328;
t232 = Icges(3,1) * t266 - Icges(3,4) * t265 - Icges(3,5) * t326;
t231 = Icges(3,4) * t268 - Icges(3,2) * t267 + Icges(3,6) * t328;
t230 = Icges(3,4) * t266 - Icges(3,2) * t265 - Icges(3,6) * t326;
t229 = Icges(3,5) * t268 - Icges(3,6) * t267 + Icges(3,3) * t328;
t228 = Icges(3,5) * t266 - Icges(3,6) * t265 - Icges(3,3) * t326;
t227 = rSges(4,1) * t264 - rSges(4,2) * t263 - rSges(4,3) * t327;
t226 = Icges(4,1) * t264 - Icges(4,4) * t263 - Icges(4,5) * t327;
t225 = Icges(4,4) * t264 - Icges(4,2) * t263 - Icges(4,6) * t327;
t224 = Icges(4,5) * t264 - Icges(4,6) * t263 - Icges(4,3) * t327;
t208 = rSges(4,1) * t248 - rSges(4,2) * t247 + rSges(4,3) * t267;
t207 = rSges(4,1) * t246 - rSges(4,2) * t245 + rSges(4,3) * t265;
t205 = Icges(4,1) * t248 - Icges(4,4) * t247 + Icges(4,5) * t267;
t204 = Icges(4,1) * t246 - Icges(4,4) * t245 + Icges(4,5) * t265;
t203 = Icges(4,4) * t248 - Icges(4,2) * t247 + Icges(4,6) * t267;
t202 = Icges(4,4) * t246 - Icges(4,2) * t245 + Icges(4,6) * t265;
t201 = Icges(4,5) * t248 - Icges(4,6) * t247 + Icges(4,3) * t267;
t200 = Icges(4,5) * t246 - Icges(4,6) * t245 + Icges(4,3) * t265;
t199 = rSges(5,1) * t242 - rSges(5,2) * t241 + rSges(5,3) * t263;
t198 = rSges(6,1) * t242 + rSges(6,2) * t263 + rSges(6,3) * t241;
t182 = -t234 * t278 + t256 * t276 + t314;
t181 = t278 * t235 - t277 * t256 + t310;
t179 = rSges(5,1) * t220 - rSges(5,2) * t219 + rSges(5,3) * t247;
t178 = rSges(6,1) * t220 + rSges(6,2) * t247 + rSges(6,3) * t219;
t176 = rSges(5,1) * t218 - rSges(5,2) * t217 + rSges(5,3) * t245;
t175 = rSges(6,1) * t218 + rSges(6,2) * t245 + rSges(6,3) * t217;
t155 = t234 * t277 - t235 * t276 + t313;
t152 = -t207 * t261 + t227 * t243 + t312;
t151 = t261 * t208 - t244 * t227 + t307;
t150 = t207 * t244 - t208 * t243 + t311;
t149 = -t176 * t237 + t199 * t213 + t309;
t148 = t237 * t179 - t214 * t199 + t304;
t147 = t176 * t214 - t179 * t213 + t308;
t146 = t198 * t213 + (-t175 - t183) * t237 + t306;
t145 = t237 * t178 + (-t198 - t210) * t214 + t303;
t144 = t175 * t214 + (-t178 - t184) * t213 + t305;
t143 = -qJD(6) * t247 + t323 * t213 + (-t183 - t325) * t237 + t306;
t142 = -qJD(6) * t245 + t324 * t237 + (-t210 - t323) * t214 + t303;
t141 = -qJD(6) * t263 + t325 * t214 + (-t184 - t324) * t213 + t305;
t1 = m(1) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t244 * ((t201 * t267 - t203 * t247 + t205 * t248) * t244 + (t200 * t267 - t202 * t247 + t204 * t248) * t243 + (t224 * t267 - t225 * t247 + t226 * t248) * t261) / 0.2e1 + t243 * ((t201 * t265 - t203 * t245 + t205 * t246) * t244 + (t200 * t265 - t202 * t245 + t204 * t246) * t243 + (t224 * t265 - t225 * t245 + t226 * t246) * t261) / 0.2e1 + m(2) * (t249 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(3) * (t155 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t278 * (((t231 * t301 + t233 * t299) * t277 + (t230 * t301 + t232 * t299) * t276 + (t254 * t301 + t255 * t299) * t278) * t296 + (t228 * t276 + t229 * t277 + t253 * t278) * t330) / 0.2e1 + t276 * ((-t229 * t326 - t265 * t231 + t266 * t233) * t277 + (-t228 * t326 - t265 * t230 + t266 * t232) * t276 + (-t253 * t326 - t265 * t254 + t266 * t255) * t278) / 0.2e1 + t261 * ((-t201 * t327 - t203 * t263 + t205 * t264) * t244 + (-t200 * t327 - t202 * t263 + t204 * t264) * t243 + (-t224 * t327 - t225 * t263 + t226 * t264) * t261) / 0.2e1 + t277 * ((t229 * t328 - t231 * t267 + t233 * t268) * t277 + (t228 * t328 - t230 * t267 + t232 * t268) * t276 + (t253 * t328 - t254 * t267 + t255 * t268) * t278) / 0.2e1 + ((-t300 * t281 + t283 * t302 + Icges(1,4)) * V_base(5) + (-t300 * t282 + t284 * t302 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t281 * t302 + t300 * t283 + Icges(1,2)) * V_base(5) + (t282 * t302 + t300 * t284 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t217 * t337 + t218 * t336 + t245 * t338) * t237 + (t217 * t341 + t218 * t339 + t245 * t343) * t214 + (t342 * t217 + t340 * t218 + t344 * t245) * t213) * t213 / 0.2e1 + ((t219 * t337 + t220 * t336 + t247 * t338) * t237 + (t341 * t219 + t339 * t220 + t343 * t247) * t214 + (t342 * t219 + t220 * t340 + t344 * t247) * t213) * t214 / 0.2e1 + ((t337 * t241 + t336 * t242 + t338 * t263) * t237 + (t241 * t341 + t242 * t339 + t263 * t343) * t214 + (t241 * t342 + t242 * t340 + t263 * t344) * t213) * t237 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t300 + Icges(2,6) * t302) * V_base(5) + (Icges(2,5) * t302 - Icges(2,6) * t300) * V_base(4) + Icges(2,3) * t293 / 0.2e1) * t293;
T  = t1;
