% Calculate kinetic energy for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:14
% EndTime: 2019-03-09 14:08:18
% DurationCPUTime: 4.41s
% Computational Cost: add. (3637->443), mult. (5384->641), div. (0->0), fcn. (6256->14), ass. (0->194)
t320 = cos(pkin(6));
t324 = sin(qJ(1));
t326 = cos(qJ(2));
t354 = t324 * t326;
t323 = sin(qJ(2));
t327 = cos(qJ(1));
t355 = t323 * t327;
t278 = t320 * t355 + t354;
t317 = sin(pkin(12));
t319 = cos(pkin(12));
t318 = sin(pkin(6));
t357 = t318 * t327;
t248 = -t278 * t317 - t319 * t357;
t344 = t317 * t357;
t249 = t278 * t319 - t344;
t353 = t326 * t327;
t356 = t323 * t324;
t277 = -t320 * t353 + t356;
t193 = Icges(4,5) * t249 + Icges(4,6) * t248 + Icges(4,3) * t277;
t230 = Icges(3,4) * t278 - Icges(3,2) * t277 - Icges(3,6) * t357;
t371 = t193 - t230;
t280 = -t320 * t356 + t353;
t359 = t318 * t324;
t250 = -t280 * t317 + t319 * t359;
t345 = t317 * t359;
t251 = t280 * t319 + t345;
t279 = t320 * t354 + t355;
t194 = Icges(4,5) * t251 + Icges(4,6) * t250 + Icges(4,3) * t279;
t231 = Icges(3,4) * t280 - Icges(3,2) * t279 + Icges(3,6) * t359;
t370 = t194 - t231;
t360 = t318 * t323;
t275 = -t317 * t360 + t319 * t320;
t361 = t317 * t320;
t276 = t319 * t360 + t361;
t358 = t318 * t326;
t223 = Icges(4,5) * t276 + Icges(4,6) * t275 - Icges(4,3) * t358;
t262 = Icges(3,6) * t320 + (Icges(3,4) * t323 + Icges(3,2) * t326) * t318;
t369 = t223 - t262;
t364 = pkin(8) * t320;
t363 = t319 * pkin(3);
t362 = Icges(2,4) * t324;
t316 = pkin(12) + qJ(4);
t310 = cos(t316);
t351 = pkin(4) * t310;
t349 = qJD(2) * t318;
t348 = V_base(5) * pkin(7) + V_base(1);
t292 = t324 * t349 + V_base(4);
t311 = V_base(6) + qJD(1);
t343 = qJ(5) + t316;
t309 = sin(t316);
t342 = pkin(4) * t309;
t253 = qJD(4) * t279 + t292;
t293 = qJD(2) * t320 + t311;
t341 = cos(t343);
t221 = qJD(5) * t279 + t253;
t291 = -t327 * t349 + V_base(5);
t340 = t318 * t341;
t285 = t324 * pkin(1) - pkin(8) * t357;
t339 = -t285 * t311 + V_base(5) * t364 + t348;
t286 = pkin(1) * t327 + pkin(8) * t359;
t338 = V_base(4) * t285 - t286 * V_base(5) + V_base(3);
t252 = qJD(4) * t277 + t291;
t220 = qJD(5) * t277 + t252;
t281 = (pkin(2) * t323 - qJ(3) * t326) * t318;
t337 = qJD(3) * t279 + t291 * t281 + t339;
t336 = t311 * t286 + V_base(2) + (-pkin(7) - t364) * V_base(4);
t258 = (-qJD(4) - qJD(5)) * t358 + t293;
t243 = pkin(2) * t280 + qJ(3) * t279;
t335 = qJD(3) * t277 + t293 * t243 + t336;
t242 = t278 * pkin(2) + t277 * qJ(3);
t334 = -qJD(3) * t358 + t292 * t242 + t338;
t189 = -pkin(3) * t344 + pkin(9) * t277 + t278 * t363;
t227 = pkin(3) * t361 + (-pkin(9) * t326 + t323 * t363) * t318;
t333 = t291 * t227 + (-t189 - t242) * t293 + t337;
t190 = pkin(3) * t345 + pkin(9) * t279 + t280 * t363;
t332 = t293 * t190 + (-t227 - t281) * t292 + t335;
t331 = t292 * t189 + (-t190 - t243) * t291 + t334;
t162 = pkin(10) * t277 + t278 * t351 - t342 * t357;
t203 = t342 * t320 + (-pkin(10) * t326 + t323 * t351) * t318;
t273 = -qJD(4) * t358 + t293;
t330 = -t162 * t273 + t252 * t203 + t333;
t163 = pkin(10) * t279 + t280 * t351 + t342 * t359;
t329 = t273 * t163 - t203 * t253 + t332;
t328 = t253 * t162 - t163 * t252 + t331;
t325 = cos(qJ(6));
t322 = sin(qJ(6));
t313 = Icges(2,4) * t327;
t306 = sin(t343);
t302 = rSges(2,1) * t327 - t324 * rSges(2,2);
t301 = t324 * rSges(2,1) + rSges(2,2) * t327;
t299 = Icges(2,1) * t327 - t362;
t298 = Icges(2,1) * t324 + t313;
t297 = -Icges(2,2) * t324 + t313;
t296 = Icges(2,2) * t327 + t362;
t290 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t289 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t288 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t266 = rSges(3,3) * t320 + (rSges(3,1) * t323 + rSges(3,2) * t326) * t318;
t265 = t309 * t320 + t310 * t360;
t264 = -t309 * t360 + t310 * t320;
t263 = Icges(3,5) * t320 + (Icges(3,1) * t323 + Icges(3,4) * t326) * t318;
t261 = Icges(3,3) * t320 + (Icges(3,5) * t323 + Icges(3,6) * t326) * t318;
t260 = t320 * t306 + t323 * t340;
t259 = t306 * t360 - t320 * t341;
t257 = V_base(5) * rSges(2,3) - t301 * t311 + t348;
t256 = t302 * t311 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t254 = t301 * V_base(4) - t302 * V_base(5) + V_base(3);
t247 = t280 * t310 + t309 * t359;
t246 = -t280 * t309 + t310 * t359;
t245 = t278 * t310 - t309 * t357;
t244 = -t278 * t309 - t310 * t357;
t241 = t280 * t341 + t306 * t359;
t240 = t280 * t306 - t324 * t340;
t239 = t278 * t341 - t306 * t357;
t238 = t278 * t306 + t327 * t340;
t237 = t260 * t325 - t322 * t358;
t236 = -t260 * t322 - t325 * t358;
t235 = rSges(3,1) * t280 - rSges(3,2) * t279 + rSges(3,3) * t359;
t234 = t278 * rSges(3,1) - t277 * rSges(3,2) - rSges(3,3) * t357;
t233 = Icges(3,1) * t280 - Icges(3,4) * t279 + Icges(3,5) * t359;
t232 = Icges(3,1) * t278 - Icges(3,4) * t277 - Icges(3,5) * t357;
t229 = Icges(3,5) * t280 - Icges(3,6) * t279 + Icges(3,3) * t359;
t228 = Icges(3,5) * t278 - Icges(3,6) * t277 - Icges(3,3) * t357;
t226 = rSges(4,1) * t276 + rSges(4,2) * t275 - rSges(4,3) * t358;
t225 = Icges(4,1) * t276 + Icges(4,4) * t275 - Icges(4,5) * t358;
t224 = Icges(4,4) * t276 + Icges(4,2) * t275 - Icges(4,6) * t358;
t218 = pkin(5) * t260 + pkin(11) * t259;
t217 = qJD(6) * t259 + t258;
t216 = rSges(5,1) * t265 + rSges(5,2) * t264 - rSges(5,3) * t358;
t215 = Icges(5,1) * t265 + Icges(5,4) * t264 - Icges(5,5) * t358;
t214 = Icges(5,4) * t265 + Icges(5,2) * t264 - Icges(5,6) * t358;
t213 = Icges(5,5) * t265 + Icges(5,6) * t264 - Icges(5,3) * t358;
t211 = rSges(6,1) * t260 - rSges(6,2) * t259 - rSges(6,3) * t358;
t210 = Icges(6,1) * t260 - Icges(6,4) * t259 - Icges(6,5) * t358;
t209 = Icges(6,4) * t260 - Icges(6,2) * t259 - Icges(6,6) * t358;
t208 = Icges(6,5) * t260 - Icges(6,6) * t259 - Icges(6,3) * t358;
t207 = t241 * t325 + t279 * t322;
t206 = -t241 * t322 + t279 * t325;
t205 = t239 * t325 + t277 * t322;
t204 = -t239 * t322 + t277 * t325;
t202 = pkin(5) * t241 + pkin(11) * t240;
t201 = pkin(5) * t239 + pkin(11) * t238;
t200 = rSges(4,1) * t251 + rSges(4,2) * t250 + rSges(4,3) * t279;
t199 = rSges(4,1) * t249 + rSges(4,2) * t248 + rSges(4,3) * t277;
t198 = Icges(4,1) * t251 + Icges(4,4) * t250 + Icges(4,5) * t279;
t197 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t277;
t196 = Icges(4,4) * t251 + Icges(4,2) * t250 + Icges(4,6) * t279;
t195 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t277;
t192 = qJD(6) * t240 + t221;
t191 = qJD(6) * t238 + t220;
t188 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t279;
t187 = rSges(5,1) * t245 + rSges(5,2) * t244 + rSges(5,3) * t277;
t186 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t279;
t185 = Icges(5,1) * t245 + Icges(5,4) * t244 + Icges(5,5) * t277;
t184 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t279;
t183 = Icges(5,4) * t245 + Icges(5,2) * t244 + Icges(5,6) * t277;
t182 = Icges(5,5) * t247 + Icges(5,6) * t246 + Icges(5,3) * t279;
t181 = Icges(5,5) * t245 + Icges(5,6) * t244 + Icges(5,3) * t277;
t180 = rSges(6,1) * t241 - rSges(6,2) * t240 + rSges(6,3) * t279;
t179 = rSges(6,1) * t239 - rSges(6,2) * t238 + rSges(6,3) * t277;
t178 = Icges(6,1) * t241 - Icges(6,4) * t240 + Icges(6,5) * t279;
t177 = Icges(6,1) * t239 - Icges(6,4) * t238 + Icges(6,5) * t277;
t176 = Icges(6,4) * t241 - Icges(6,2) * t240 + Icges(6,6) * t279;
t175 = Icges(6,4) * t239 - Icges(6,2) * t238 + Icges(6,6) * t277;
t174 = Icges(6,5) * t241 - Icges(6,6) * t240 + Icges(6,3) * t279;
t173 = Icges(6,5) * t239 - Icges(6,6) * t238 + Icges(6,3) * t277;
t170 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t259;
t169 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t259;
t168 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t259;
t167 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t259;
t165 = -t234 * t293 + t266 * t291 + t339;
t164 = t235 * t293 - t266 * t292 + t336;
t160 = t234 * t292 - t235 * t291 + t338;
t158 = rSges(7,1) * t207 + rSges(7,2) * t206 + rSges(7,3) * t240;
t157 = rSges(7,1) * t205 + rSges(7,2) * t204 + rSges(7,3) * t238;
t156 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t240;
t155 = Icges(7,1) * t205 + Icges(7,4) * t204 + Icges(7,5) * t238;
t154 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t240;
t153 = Icges(7,4) * t205 + Icges(7,2) * t204 + Icges(7,6) * t238;
t152 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t240;
t151 = Icges(7,5) * t205 + Icges(7,6) * t204 + Icges(7,3) * t238;
t150 = t226 * t291 + (-t199 - t242) * t293 + t337;
t149 = t200 * t293 + (-t226 - t281) * t292 + t335;
t148 = t199 * t292 + (-t200 - t243) * t291 + t334;
t147 = -t187 * t273 + t216 * t252 + t333;
t146 = t188 * t273 - t216 * t253 + t332;
t145 = t187 * t253 - t188 * t252 + t331;
t144 = -t179 * t258 + t211 * t220 + t330;
t143 = t180 * t258 - t211 * t221 + t329;
t142 = t179 * t221 - t180 * t220 + t328;
t141 = -t157 * t217 + t170 * t191 - t201 * t258 + t218 * t220 + t330;
t140 = t158 * t217 - t170 * t192 + t202 * t258 - t218 * t221 + t329;
t139 = t157 * t192 - t158 * t191 + t201 * t221 - t202 * t220 + t328;
t1 = t191 * ((t152 * t238 + t154 * t204 + t156 * t205) * t192 + (t238 * t151 + t204 * t153 + t205 * t155) * t191 + (t167 * t238 + t168 * t204 + t169 * t205) * t217) / 0.2e1 + t192 * ((t240 * t152 + t206 * t154 + t207 * t156) * t192 + (t151 * t240 + t153 * t206 + t155 * t207) * t191 + (t167 * t240 + t168 * t206 + t169 * t207) * t217) / 0.2e1 + t273 * ((-t182 * t358 + t184 * t264 + t186 * t265) * t253 + (-t181 * t358 + t183 * t264 + t185 * t265) * t252 + (-t213 * t358 + t214 * t264 + t265 * t215) * t273) / 0.2e1 + t258 * ((-t174 * t358 - t176 * t259 + t178 * t260) * t221 + (-t173 * t358 - t175 * t259 + t177 * t260) * t220 + (-t208 * t358 - t209 * t259 + t210 * t260) * t258) / 0.2e1 + m(3) * (t160 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(2) * (t254 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t217 * ((t152 * t259 + t154 * t236 + t156 * t237) * t192 + (t151 * t259 + t153 * t236 + t155 * t237) * t191 + (t259 * t167 + t236 * t168 + t237 * t169) * t217) / 0.2e1 + t220 * ((t174 * t277 - t176 * t238 + t178 * t239) * t221 + (t277 * t173 - t238 * t175 + t239 * t177) * t220 + (t208 * t277 - t209 * t238 + t210 * t239) * t258) / 0.2e1 + t252 * ((t182 * t277 + t184 * t244 + t186 * t245) * t253 + (t181 * t277 + t183 * t244 + t185 * t245) * t252 + (t213 * t277 + t214 * t244 + t215 * t245) * t273) / 0.2e1 + t221 * ((t279 * t174 - t240 * t176 + t241 * t178) * t221 + (t173 * t279 - t175 * t240 + t177 * t241) * t220 + (t208 * t279 - t209 * t240 + t210 * t241) * t258) / 0.2e1 + t253 * ((t182 * t279 + t184 * t246 + t186 * t247) * t253 + (t181 * t279 + t183 * t246 + t185 * t247) * t252 + (t213 * t279 + t214 * t246 + t215 * t247) * t273) / 0.2e1 + m(1) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + ((t224 * t248 + t225 * t249 - t261 * t357 + t278 * t263 + t277 * t369) * t293 + (t196 * t248 + t198 * t249 - t229 * t357 + t278 * t233 + t277 * t370) * t292 + (t195 * t248 + t197 * t249 - t228 * t357 + t278 * t232 + t277 * t371) * t291) * t291 / 0.2e1 + ((t224 * t250 + t225 * t251 + t261 * t359 + t263 * t280 + t279 * t369) * t293 + (t196 * t250 + t198 * t251 + t229 * t359 + t233 * t280 + t279 * t370) * t292 + (t195 * t250 + t197 * t251 + t228 * t359 + t232 * t280 + t279 * t371) * t291) * t292 / 0.2e1 + ((t228 * t291 + t229 * t292 + t261 * t293) * t320 + ((t231 * t326 + t233 * t323) * t292 + (t230 * t326 + t232 * t323) * t291 + (t262 * t326 + t263 * t323) * t293) * t318 + (-t194 * t358 + t196 * t275 + t198 * t276) * t292 + (-t193 * t358 + t195 * t275 + t197 * t276) * t291 + (-t223 * t358 + t224 * t275 + t225 * t276) * t293) * t293 / 0.2e1 + ((-t324 * t296 + t298 * t327 + Icges(1,4)) * V_base(5) + (-t324 * t297 + t299 * t327 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t327 + t324 * t298 + Icges(1,2)) * V_base(5) + (t297 * t327 + t324 * t299 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t324 + Icges(2,6) * t327) * V_base(5) + (Icges(2,5) * t327 - Icges(2,6) * t324) * V_base(4) + Icges(2,3) * t311 / 0.2e1) * t311;
T  = t1;
