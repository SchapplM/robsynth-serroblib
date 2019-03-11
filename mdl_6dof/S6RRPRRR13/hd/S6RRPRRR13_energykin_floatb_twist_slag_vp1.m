% Calculate kinetic energy for
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:51
% EndTime: 2019-03-09 14:44:56
% DurationCPUTime: 4.39s
% Computational Cost: add. (2569->396), mult. (5438->572), div. (0->0), fcn. (6418->12), ass. (0->177)
t360 = Icges(3,1) + Icges(4,2);
t359 = Icges(3,4) + Icges(4,6);
t358 = Icges(3,5) - Icges(4,4);
t357 = Icges(3,2) + Icges(4,3);
t356 = Icges(3,6) - Icges(4,5);
t355 = Icges(3,3) + Icges(4,1);
t301 = cos(pkin(6));
t305 = sin(qJ(1));
t307 = cos(qJ(2));
t329 = t305 * t307;
t304 = sin(qJ(2));
t308 = cos(qJ(1));
t330 = t304 * t308;
t266 = t301 * t329 + t330;
t328 = t307 * t308;
t331 = t304 * t305;
t267 = -t301 * t331 + t328;
t300 = sin(pkin(6));
t334 = t300 * t305;
t354 = t357 * t266 - t359 * t267 - t356 * t334;
t264 = -t301 * t328 + t331;
t265 = t301 * t330 + t329;
t332 = t300 * t308;
t353 = t357 * t264 - t359 * t265 + t356 * t332;
t352 = -t359 * t266 + t360 * t267 + t358 * t334;
t351 = -t359 * t264 + t360 * t265 - t358 * t332;
t350 = -t356 * t266 + t358 * t267 + t355 * t334;
t349 = -t356 * t264 + t358 * t265 - t355 * t332;
t348 = t355 * t301 + (t358 * t304 + t356 * t307) * t300;
t347 = t356 * t301 + (t359 * t304 + t357 * t307) * t300;
t346 = t358 * t301 + (t360 * t304 + t359 * t307) * t300;
t342 = cos(qJ(4));
t341 = pkin(8) * t301;
t306 = cos(qJ(5));
t340 = pkin(5) * t306;
t338 = Icges(2,4) * t305;
t302 = sin(qJ(5));
t337 = t265 * t302;
t336 = t267 * t302;
t335 = t300 * t304;
t333 = t300 * t307;
t327 = qJD(2) * t300;
t326 = V_base(5) * pkin(7) + V_base(1);
t323 = t302 * t335;
t322 = t300 * t342;
t277 = t305 * t327 + V_base(4);
t294 = V_base(6) + qJD(1);
t234 = qJD(4) * t267 + t277;
t278 = qJD(2) * t301 + t294;
t303 = sin(qJ(4));
t235 = -t266 * t342 + t303 * t334;
t189 = qJD(5) * t235 + t234;
t260 = qJD(4) * t335 + t278;
t276 = -t308 * t327 + V_base(5);
t271 = pkin(1) * t305 - pkin(8) * t332;
t321 = -t271 * t294 + V_base(5) * t341 + t326;
t262 = t301 * t303 + t307 * t322;
t223 = qJD(5) * t262 + t260;
t272 = pkin(1) * t308 + pkin(8) * t334;
t320 = V_base(4) * t271 - t272 * V_base(5) + V_base(3);
t233 = qJD(4) * t265 + t276;
t237 = t264 * t342 + t303 * t332;
t188 = -qJD(5) * t237 + t233;
t268 = (pkin(2) * t304 - qJ(3) * t307) * t300;
t319 = qJD(3) * t266 + t276 * t268 + t321;
t318 = t294 * t272 + V_base(2) + (-pkin(7) - t341) * V_base(4);
t226 = pkin(2) * t267 + qJ(3) * t266;
t317 = qJD(3) * t264 + t278 * t226 + t318;
t225 = pkin(2) * t265 + qJ(3) * t264;
t316 = -qJD(3) * t333 + t277 * t225 + t320;
t244 = -pkin(3) * t332 + pkin(9) * t265;
t270 = pkin(3) * t301 + pkin(9) * t335;
t315 = t276 * t270 + (-t225 - t244) * t278 + t319;
t243 = pkin(3) * t334 + pkin(9) * t267;
t314 = t278 * t243 + (-t268 - t270) * t277 + t317;
t313 = t277 * t244 + (-t226 - t243) * t276 + t316;
t238 = t264 * t303 - t308 * t322;
t187 = t238 * pkin(4) - t237 * pkin(10);
t263 = t301 * t342 - t303 * t333;
t222 = t263 * pkin(4) + t262 * pkin(10);
t312 = -t187 * t260 + t233 * t222 + t315;
t236 = t266 * t303 + t305 * t322;
t186 = t236 * pkin(4) + t235 * pkin(10);
t311 = t260 * t186 - t222 * t234 + t314;
t310 = -t186 * t233 + t234 * t187 + t313;
t299 = qJ(5) + qJ(6);
t297 = Icges(2,4) * t308;
t296 = cos(t299);
t295 = sin(t299);
t286 = rSges(2,1) * t308 - rSges(2,2) * t305;
t285 = rSges(2,1) * t305 + rSges(2,2) * t308;
t284 = Icges(2,1) * t308 - t338;
t283 = Icges(2,1) * t305 + t297;
t282 = -Icges(2,2) * t305 + t297;
t281 = Icges(2,2) * t308 + t338;
t275 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t274 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t273 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t253 = rSges(4,1) * t301 + (-rSges(4,2) * t304 - rSges(4,3) * t307) * t300;
t252 = rSges(3,3) * t301 + (rSges(3,1) * t304 + rSges(3,2) * t307) * t300;
t242 = V_base(5) * rSges(2,3) - t285 * t294 + t326;
t241 = t286 * t294 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t239 = t285 * V_base(4) - t286 * V_base(5) + V_base(3);
t232 = t263 * t306 + t323;
t231 = -t263 * t302 + t306 * t335;
t228 = t263 * t296 + t295 * t335;
t227 = -t263 * t295 + t296 * t335;
t220 = rSges(3,1) * t267 - rSges(3,2) * t266 + rSges(3,3) * t334;
t219 = rSges(3,1) * t265 - rSges(3,2) * t264 - rSges(3,3) * t332;
t218 = -rSges(4,1) * t332 - rSges(4,2) * t265 + rSges(4,3) * t264;
t217 = rSges(4,1) * t334 - rSges(4,2) * t267 + rSges(4,3) * t266;
t204 = rSges(5,1) * t263 - rSges(5,2) * t262 + rSges(5,3) * t335;
t203 = Icges(5,1) * t263 - Icges(5,4) * t262 + Icges(5,5) * t335;
t202 = Icges(5,4) * t263 - Icges(5,2) * t262 + Icges(5,6) * t335;
t201 = Icges(5,5) * t263 - Icges(5,6) * t262 + Icges(5,3) * t335;
t198 = t238 * t306 + t337;
t197 = -t238 * t302 + t265 * t306;
t196 = t236 * t306 + t336;
t195 = -t236 * t302 + t267 * t306;
t194 = qJD(6) * t262 + t223;
t193 = t238 * t296 + t265 * t295;
t192 = -t238 * t295 + t265 * t296;
t191 = t236 * t296 + t267 * t295;
t190 = -t236 * t295 + t267 * t296;
t184 = rSges(5,1) * t238 + rSges(5,2) * t237 + rSges(5,3) * t265;
t183 = rSges(5,1) * t236 - rSges(5,2) * t235 + rSges(5,3) * t267;
t181 = Icges(5,1) * t238 + Icges(5,4) * t237 + Icges(5,5) * t265;
t180 = Icges(5,1) * t236 - Icges(5,4) * t235 + Icges(5,5) * t267;
t179 = Icges(5,4) * t238 + Icges(5,2) * t237 + Icges(5,6) * t265;
t178 = Icges(5,4) * t236 - Icges(5,2) * t235 + Icges(5,6) * t267;
t177 = Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t265;
t176 = Icges(5,5) * t236 - Icges(5,6) * t235 + Icges(5,3) * t267;
t175 = rSges(6,1) * t232 + rSges(6,2) * t231 + rSges(6,3) * t262;
t174 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t262;
t173 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t262;
t172 = Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t262;
t171 = pkin(5) * t323 + pkin(11) * t262 + t263 * t340;
t170 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t262;
t169 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t262;
t168 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t262;
t167 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t262;
t166 = qJD(6) * t235 + t189;
t165 = -qJD(6) * t237 + t188;
t163 = -t219 * t278 + t252 * t276 + t321;
t162 = t220 * t278 - t252 * t277 + t318;
t161 = rSges(6,1) * t198 + rSges(6,2) * t197 - rSges(6,3) * t237;
t160 = rSges(6,1) * t196 + rSges(6,2) * t195 + rSges(6,3) * t235;
t159 = Icges(6,1) * t198 + Icges(6,4) * t197 - Icges(6,5) * t237;
t158 = Icges(6,1) * t196 + Icges(6,4) * t195 + Icges(6,5) * t235;
t157 = Icges(6,4) * t198 + Icges(6,2) * t197 - Icges(6,6) * t237;
t156 = Icges(6,4) * t196 + Icges(6,2) * t195 + Icges(6,6) * t235;
t155 = Icges(6,5) * t198 + Icges(6,6) * t197 - Icges(6,3) * t237;
t154 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t235;
t153 = t219 * t277 - t220 * t276 + t320;
t152 = rSges(7,1) * t193 + rSges(7,2) * t192 - rSges(7,3) * t237;
t151 = rSges(7,1) * t191 + rSges(7,2) * t190 + rSges(7,3) * t235;
t150 = Icges(7,1) * t193 + Icges(7,4) * t192 - Icges(7,5) * t237;
t149 = Icges(7,1) * t191 + Icges(7,4) * t190 + Icges(7,5) * t235;
t148 = Icges(7,4) * t193 + Icges(7,2) * t192 - Icges(7,6) * t237;
t147 = Icges(7,4) * t191 + Icges(7,2) * t190 + Icges(7,6) * t235;
t146 = Icges(7,5) * t193 + Icges(7,6) * t192 - Icges(7,3) * t237;
t145 = Icges(7,5) * t191 + Icges(7,6) * t190 + Icges(7,3) * t235;
t144 = pkin(5) * t337 - pkin(11) * t237 + t238 * t340;
t143 = pkin(5) * t336 + pkin(11) * t235 + t236 * t340;
t142 = t253 * t276 + (-t218 - t225) * t278 + t319;
t141 = t217 * t278 + (-t253 - t268) * t277 + t317;
t140 = t218 * t277 + (-t217 - t226) * t276 + t316;
t139 = -t184 * t260 + t204 * t233 + t315;
t138 = t183 * t260 - t204 * t234 + t314;
t137 = -t183 * t233 + t184 * t234 + t313;
t136 = -t161 * t223 + t175 * t188 + t312;
t135 = t160 * t223 - t175 * t189 + t311;
t134 = -t160 * t188 + t161 * t189 + t310;
t133 = -t144 * t223 - t152 * t194 + t165 * t170 + t171 * t188 + t312;
t132 = t143 * t223 + t151 * t194 - t166 * t170 - t171 * t189 + t311;
t131 = -t143 * t188 + t144 * t189 - t151 * t165 + t152 * t166 + t310;
t1 = m(1) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t234 * ((t267 * t176 - t235 * t178 + t236 * t180) * t234 + (t177 * t267 - t179 * t235 + t181 * t236) * t233 + (t201 * t267 - t202 * t235 + t203 * t236) * t260) / 0.2e1 + t233 * ((t176 * t265 + t178 * t237 + t180 * t238) * t234 + (t265 * t177 + t237 * t179 + t238 * t181) * t233 + (t201 * t265 + t202 * t237 + t203 * t238) * t260) / 0.2e1 + t223 * ((t154 * t262 + t156 * t231 + t158 * t232) * t189 + (t155 * t262 + t157 * t231 + t159 * t232) * t188 + (t262 * t172 + t231 * t173 + t232 * t174) * t223) / 0.2e1 + t194 * ((t145 * t262 + t147 * t227 + t149 * t228) * t166 + (t146 * t262 + t148 * t227 + t150 * t228) * t165 + (t262 * t167 + t227 * t168 + t228 * t169) * t194) / 0.2e1 + m(2) * (t239 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + t188 * ((-t154 * t237 + t156 * t197 + t158 * t198) * t189 + (-t237 * t155 + t197 * t157 + t198 * t159) * t188 + (-t172 * t237 + t173 * t197 + t174 * t198) * t223) / 0.2e1 + t165 * ((-t145 * t237 + t147 * t192 + t149 * t193) * t166 + (-t237 * t146 + t192 * t148 + t193 * t150) * t165 + (-t167 * t237 + t168 * t192 + t169 * t193) * t194) / 0.2e1 + t189 * ((t235 * t154 + t195 * t156 + t196 * t158) * t189 + (t155 * t235 + t157 * t195 + t159 * t196) * t188 + (t172 * t235 + t173 * t195 + t174 * t196) * t223) / 0.2e1 + t166 * ((t235 * t145 + t190 * t147 + t191 * t149) * t166 + (t146 * t235 + t148 * t190 + t150 * t191) * t165 + (t167 * t235 + t168 * t190 + t169 * t191) * t194) / 0.2e1 + m(3) * (t153 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + t260 * ((t176 * t335 - t178 * t262 + t180 * t263) * t234 + (t177 * t335 - t179 * t262 + t181 * t263) * t233 + (t201 * t335 - t202 * t262 + t203 * t263) * t260) / 0.2e1 + ((-t264 * t347 + t265 * t346 - t332 * t348) * t278 + (t264 * t354 + t352 * t265 - t350 * t332) * t277 + (t353 * t264 + t351 * t265 - t349 * t332) * t276) * t276 / 0.2e1 + ((-t266 * t347 + t267 * t346 + t334 * t348) * t278 + (t354 * t266 + t352 * t267 + t350 * t334) * t277 + (t266 * t353 + t267 * t351 + t334 * t349) * t276) * t277 / 0.2e1 + ((t276 * t349 + t277 * t350 + t278 * t348) * t301 + ((t304 * t346 + t307 * t347) * t278 + (t352 * t304 - t307 * t354) * t277 + (t304 * t351 - t307 * t353) * t276) * t300) * t278 / 0.2e1 + ((-t281 * t305 + t283 * t308 + Icges(1,4)) * V_base(5) + (-t282 * t305 + t284 * t308 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t281 * t308 + t283 * t305 + Icges(1,2)) * V_base(5) + (t282 * t308 + t284 * t305 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t305 + Icges(2,6) * t308) * V_base(5) + (Icges(2,5) * t308 - Icges(2,6) * t305) * V_base(4) + Icges(2,3) * t294 / 0.2e1) * t294;
T  = t1;
