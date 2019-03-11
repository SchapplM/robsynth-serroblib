% Calculate kinetic energy for
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:48
% EndTime: 2019-03-08 20:40:53
% DurationCPUTime: 4.40s
% Computational Cost: add. (2587->396), mult. (4786->574), div. (0->0), fcn. (5502->12), ass. (0->177)
t365 = Icges(3,1) + Icges(4,2);
t364 = Icges(3,4) + Icges(4,6);
t363 = Icges(3,5) - Icges(4,4);
t362 = Icges(3,2) + Icges(4,3);
t361 = Icges(3,6) - Icges(4,5);
t360 = Icges(3,3) + Icges(4,1);
t299 = sin(pkin(11));
t301 = cos(pkin(11));
t305 = sin(qJ(2));
t302 = cos(pkin(6));
t308 = cos(qJ(2));
t332 = t302 * t308;
t263 = t299 * t332 + t301 * t305;
t333 = t302 * t305;
t264 = -t299 * t333 + t301 * t308;
t300 = sin(pkin(6));
t339 = t299 * t300;
t358 = t362 * t263 - t364 * t264 - t361 * t339;
t261 = t299 * t305 - t301 * t332;
t262 = t299 * t308 + t301 * t333;
t338 = t300 * t301;
t357 = t362 * t261 - t364 * t262 + t361 * t338;
t356 = -t364 * t263 + t365 * t264 + t363 * t339;
t355 = -t364 * t261 + t365 * t262 - t363 * t338;
t354 = -t361 * t263 + t363 * t264 + t360 * t339;
t353 = -t361 * t261 + t363 * t262 - t360 * t338;
t352 = t360 * t302 + (t363 * t305 + t361 * t308) * t300;
t351 = t361 * t302 + (t364 * t305 + t362 * t308) * t300;
t350 = t363 * t302 + (t365 * t305 + t364 * t308) * t300;
t345 = pkin(7) * t302;
t307 = cos(qJ(4));
t344 = pkin(4) * t307;
t342 = Icges(2,4) * t299;
t304 = sin(qJ(4));
t341 = t261 * t304;
t340 = t263 * t304;
t337 = t300 * t304;
t336 = t300 * t305;
t335 = t300 * t307;
t334 = t300 * t308;
t331 = t304 * t308;
t330 = qJ(4) + qJ(5);
t329 = qJD(2) * t300;
t328 = V_base(5) * qJ(1) + V_base(1);
t324 = qJD(1) + V_base(3);
t278 = t299 * t329 + V_base(4);
t289 = qJD(2) * t302 + V_base(6);
t323 = cos(t330);
t232 = qJD(4) * t264 + t278;
t265 = qJD(4) * t336 + t289;
t322 = t300 * t323;
t195 = qJD(5) * t264 + t232;
t243 = qJD(5) * t336 + t265;
t277 = -t301 * t329 + V_base(5);
t231 = qJD(4) * t262 + t277;
t271 = pkin(1) * t299 - pkin(7) * t338;
t321 = -t271 * V_base(6) + V_base(5) * t345 + t328;
t272 = pkin(1) * t301 + pkin(7) * t339;
t320 = V_base(4) * t271 - t272 * V_base(5) + t324;
t194 = qJD(5) * t262 + t231;
t319 = V_base(6) * t272 + V_base(2) + (-qJ(1) - t345) * V_base(4);
t270 = (pkin(2) * t305 - qJ(3) * t308) * t300;
t318 = qJD(3) * t263 + t277 * t270 + t321;
t222 = pkin(2) * t264 + qJ(3) * t263;
t317 = qJD(3) * t261 + t289 * t222 + t319;
t221 = pkin(2) * t262 + qJ(3) * t261;
t316 = -qJD(3) * t334 + t278 * t221 + t320;
t239 = -pkin(3) * t338 + t262 * pkin(8);
t273 = t302 * pkin(3) + pkin(8) * t336;
t315 = t277 * t273 + (-t221 - t239) * t289 + t318;
t238 = pkin(3) * t339 + t264 * pkin(8);
t314 = t289 * t238 + (-t270 - t273) * t278 + t317;
t313 = t278 * t239 + (-t222 - t238) * t277 + t316;
t181 = pkin(4) * t341 + pkin(9) * t262 - t338 * t344;
t219 = t344 * t302 + (-pkin(4) * t331 + pkin(9) * t305) * t300;
t312 = -t181 * t265 + t231 * t219 + t315;
t180 = pkin(4) * t340 + pkin(9) * t264 + t339 * t344;
t311 = t265 * t180 - t219 * t232 + t314;
t310 = -t180 * t231 + t232 * t181 + t313;
t306 = cos(qJ(6));
t303 = sin(qJ(6));
t297 = sin(t330);
t296 = Icges(2,4) * t301;
t286 = rSges(2,1) * t301 - rSges(2,2) * t299;
t285 = rSges(2,1) * t299 + rSges(2,2) * t301;
t284 = Icges(2,1) * t301 - t342;
t283 = Icges(2,1) * t299 + t296;
t282 = -Icges(2,2) * t299 + t296;
t281 = Icges(2,2) * t301 + t342;
t276 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t275 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t274 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t269 = -t300 * t331 + t302 * t307;
t268 = -t302 * t304 - t307 * t334;
t254 = -t297 * t334 + t302 * t323;
t253 = t302 * t297 + t308 * t322;
t252 = rSges(4,1) * t302 + (-rSges(4,2) * t305 - rSges(4,3) * t308) * t300;
t251 = rSges(3,3) * t302 + (rSges(3,1) * t305 + rSges(3,2) * t308) * t300;
t241 = V_base(5) * rSges(2,3) - t285 * V_base(6) + t328;
t240 = t286 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t236 = -t301 * t335 + t341;
t235 = t261 * t307 + t301 * t337;
t234 = t299 * t335 + t340;
t233 = t263 * t307 - t299 * t337;
t230 = t285 * V_base(4) - t286 * V_base(5) + t324;
t229 = t254 * t306 + t303 * t336;
t228 = -t254 * t303 + t306 * t336;
t227 = t261 * t297 - t301 * t322;
t226 = t261 * t323 + t297 * t338;
t225 = t263 * t297 + t299 * t322;
t224 = -t263 * t323 + t297 * t339;
t218 = pkin(5) * t254 + pkin(10) * t253;
t217 = rSges(5,1) * t269 + rSges(5,2) * t268 + rSges(5,3) * t336;
t216 = Icges(5,1) * t269 + Icges(5,4) * t268 + Icges(5,5) * t336;
t215 = Icges(5,4) * t269 + Icges(5,2) * t268 + Icges(5,6) * t336;
t214 = Icges(5,5) * t269 + Icges(5,6) * t268 + Icges(5,3) * t336;
t213 = rSges(3,1) * t264 - rSges(3,2) * t263 + rSges(3,3) * t339;
t212 = rSges(3,1) * t262 - rSges(3,2) * t261 - rSges(3,3) * t338;
t211 = -rSges(4,1) * t338 - rSges(4,2) * t262 + rSges(4,3) * t261;
t210 = rSges(4,1) * t339 - rSges(4,2) * t264 + rSges(4,3) * t263;
t197 = qJD(6) * t253 + t243;
t192 = rSges(6,1) * t254 - rSges(6,2) * t253 + rSges(6,3) * t336;
t191 = Icges(6,1) * t254 - Icges(6,4) * t253 + Icges(6,5) * t336;
t190 = Icges(6,4) * t254 - Icges(6,2) * t253 + Icges(6,6) * t336;
t189 = Icges(6,5) * t254 - Icges(6,6) * t253 + Icges(6,3) * t336;
t188 = t227 * t306 + t262 * t303;
t187 = -t227 * t303 + t262 * t306;
t186 = t225 * t306 + t264 * t303;
t185 = -t225 * t303 + t264 * t306;
t184 = pkin(5) * t227 - pkin(10) * t226;
t183 = pkin(5) * t225 + pkin(10) * t224;
t179 = rSges(5,1) * t236 + rSges(5,2) * t235 + rSges(5,3) * t262;
t178 = rSges(5,1) * t234 + rSges(5,2) * t233 + rSges(5,3) * t264;
t177 = Icges(5,1) * t236 + Icges(5,4) * t235 + Icges(5,5) * t262;
t176 = Icges(5,1) * t234 + Icges(5,4) * t233 + Icges(5,5) * t264;
t175 = Icges(5,4) * t236 + Icges(5,2) * t235 + Icges(5,6) * t262;
t174 = Icges(5,4) * t234 + Icges(5,2) * t233 + Icges(5,6) * t264;
t173 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t262;
t172 = Icges(5,5) * t234 + Icges(5,6) * t233 + Icges(5,3) * t264;
t171 = qJD(6) * t224 + t195;
t170 = -qJD(6) * t226 + t194;
t169 = rSges(6,1) * t227 + rSges(6,2) * t226 + rSges(6,3) * t262;
t168 = rSges(6,1) * t225 - rSges(6,2) * t224 + rSges(6,3) * t264;
t167 = Icges(6,1) * t227 + Icges(6,4) * t226 + Icges(6,5) * t262;
t166 = Icges(6,1) * t225 - Icges(6,4) * t224 + Icges(6,5) * t264;
t165 = Icges(6,4) * t227 + Icges(6,2) * t226 + Icges(6,6) * t262;
t164 = Icges(6,4) * t225 - Icges(6,2) * t224 + Icges(6,6) * t264;
t163 = Icges(6,5) * t227 + Icges(6,6) * t226 + Icges(6,3) * t262;
t162 = Icges(6,5) * t225 - Icges(6,6) * t224 + Icges(6,3) * t264;
t161 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t253;
t160 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t253;
t159 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t253;
t158 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t253;
t155 = -t212 * t289 + t251 * t277 + t321;
t154 = t213 * t289 - t251 * t278 + t319;
t153 = t212 * t278 - t213 * t277 + t320;
t152 = rSges(7,1) * t188 + rSges(7,2) * t187 - rSges(7,3) * t226;
t151 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t224;
t150 = Icges(7,1) * t188 + Icges(7,4) * t187 - Icges(7,5) * t226;
t149 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t224;
t148 = Icges(7,4) * t188 + Icges(7,2) * t187 - Icges(7,6) * t226;
t147 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t224;
t146 = Icges(7,5) * t188 + Icges(7,6) * t187 - Icges(7,3) * t226;
t145 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t224;
t144 = t252 * t277 + (-t211 - t221) * t289 + t318;
t143 = t210 * t289 + (-t252 - t270) * t278 + t317;
t142 = t211 * t278 + (-t210 - t222) * t277 + t316;
t141 = -t179 * t265 + t217 * t231 + t315;
t140 = t178 * t265 - t217 * t232 + t314;
t139 = -t178 * t231 + t179 * t232 + t313;
t138 = -t169 * t243 + t192 * t194 + t312;
t137 = t168 * t243 - t192 * t195 + t311;
t136 = -t168 * t194 + t169 * t195 + t310;
t135 = -t152 * t197 + t161 * t170 - t184 * t243 + t194 * t218 + t312;
t134 = t151 * t197 - t161 * t171 + t183 * t243 - t195 * t218 + t311;
t133 = -t151 * t170 + t152 * t171 - t183 * t194 + t184 * t195 + t310;
t1 = t194 * ((t162 * t262 + t164 * t226 + t166 * t227) * t195 + (t262 * t163 + t226 * t165 + t227 * t167) * t194 + (t189 * t262 + t190 * t226 + t191 * t227) * t243) / 0.2e1 + t171 * ((t224 * t145 + t185 * t147 + t186 * t149) * t171 + (t146 * t224 + t148 * t185 + t150 * t186) * t170 + (t158 * t224 + t159 * t185 + t160 * t186) * t197) / 0.2e1 + t197 * ((t145 * t253 + t147 * t228 + t149 * t229) * t171 + (t146 * t253 + t148 * t228 + t150 * t229) * t170 + (t253 * t158 + t228 * t159 + t229 * t160) * t197) / 0.2e1 + t170 * ((-t145 * t226 + t147 * t187 + t149 * t188) * t171 + (-t226 * t146 + t187 * t148 + t188 * t150) * t170 + (-t158 * t226 + t159 * t187 + t160 * t188) * t197) / 0.2e1 + m(1) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + t195 * ((t264 * t162 - t224 * t164 + t225 * t166) * t195 + (t163 * t264 - t165 * t224 + t167 * t225) * t194 + (t189 * t264 - t190 * t224 + t191 * t225) * t243) / 0.2e1 + t231 * ((t172 * t262 + t174 * t235 + t176 * t236) * t232 + (t262 * t173 + t235 * t175 + t236 * t177) * t231 + (t214 * t262 + t215 * t235 + t216 * t236) * t265) / 0.2e1 + t232 * ((t264 * t172 + t233 * t174 + t234 * t176) * t232 + (t173 * t264 + t175 * t233 + t177 * t234) * t231 + (t214 * t264 + t215 * t233 + t216 * t234) * t265) / 0.2e1 + m(3) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(2) * (t230 ^ 2 + t240 ^ 2 + t241 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t265 * ((t172 * t336 + t174 * t268 + t176 * t269) * t232 + (t173 * t336 + t175 * t268 + t177 * t269) * t231 + (t214 * t336 + t268 * t215 + t269 * t216) * t265) / 0.2e1 + t243 * ((t162 * t336 - t164 * t253 + t166 * t254) * t195 + (t163 * t336 - t165 * t253 + t167 * t254) * t194 + (t189 * t336 - t253 * t190 + t254 * t191) * t243) / 0.2e1 + ((-t261 * t351 + t262 * t350 - t338 * t352) * t289 + (t261 * t358 + t262 * t356 - t338 * t354) * t278 + (t261 * t357 + t262 * t355 - t338 * t353) * t277) * t277 / 0.2e1 + ((-t263 * t351 + t264 * t350 + t339 * t352) * t289 + (t263 * t358 + t264 * t356 + t339 * t354) * t278 + (t263 * t357 + t264 * t355 + t339 * t353) * t277) * t278 / 0.2e1 + ((t277 * t353 + t278 * t354 + t289 * t352) * t302 + ((t305 * t350 + t308 * t351) * t289 + (t305 * t356 - t308 * t358) * t278 + (t305 * t355 - t357 * t308) * t277) * t300) * t289 / 0.2e1 + ((-t281 * t299 + t283 * t301 + Icges(1,4)) * V_base(5) + (-t282 * t299 + t284 * t301 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t281 * t301 + t283 * t299 + Icges(1,2)) * V_base(5) + (t282 * t301 + t284 * t299 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t299 + Icges(2,6) * t301 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t301 - Icges(2,6) * t299 + Icges(1,5)) * V_base(4) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
