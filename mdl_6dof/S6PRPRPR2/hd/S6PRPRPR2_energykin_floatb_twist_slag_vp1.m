% Calculate kinetic energy for
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:26
% EndTime: 2019-03-08 19:29:31
% DurationCPUTime: 4.71s
% Computational Cost: add. (3825->434), mult. (8792->609), div. (0->0), fcn. (11044->14), ass. (0->190)
t378 = Icges(5,2) + Icges(6,3);
t318 = cos(pkin(6));
t321 = sin(qJ(2));
t357 = sin(pkin(11));
t358 = cos(pkin(11));
t363 = cos(qJ(2));
t333 = t321 * t358 + t357 * t363;
t271 = t333 * t318;
t280 = -t321 * t357 + t363 * t358;
t314 = sin(pkin(10));
t317 = cos(pkin(10));
t252 = t271 * t317 + t280 * t314;
t315 = sin(pkin(6));
t320 = sin(qJ(4));
t350 = t315 * t320;
t362 = cos(qJ(4));
t230 = t252 * t362 - t317 * t350;
t329 = t318 * t280;
t251 = -t314 * t333 + t317 * t329;
t313 = sin(pkin(12));
t316 = cos(pkin(12));
t197 = -t230 * t313 - t251 * t316;
t355 = t251 * t313;
t198 = t230 * t316 - t355;
t340 = t315 * t362;
t229 = t252 * t320 + t317 * t340;
t376 = -Icges(5,4) * t230 + Icges(6,5) * t198 + Icges(5,6) * t251 + Icges(6,6) * t197 + t378 * t229;
t254 = -t271 * t314 + t280 * t317;
t232 = t254 * t362 + t314 * t350;
t253 = -t314 * t329 - t317 * t333;
t199 = -t232 * t313 - t253 * t316;
t354 = t253 * t313;
t200 = t232 * t316 - t354;
t231 = t254 * t320 - t314 * t340;
t375 = -Icges(5,4) * t232 + Icges(6,5) * t200 + Icges(5,6) * t253 + Icges(6,6) * t199 + t378 * t231;
t270 = t333 * t315;
t259 = t270 * t362 + t318 * t320;
t269 = t280 * t315;
t221 = -t259 * t313 - t269 * t316;
t353 = t269 * t313;
t222 = t259 * t316 - t353;
t258 = t270 * t320 - t318 * t362;
t374 = -Icges(5,4) * t259 + Icges(6,5) * t222 + Icges(5,6) * t269 + Icges(6,6) * t221 + t378 * t258;
t351 = t315 * t317;
t203 = Icges(4,5) * t252 + Icges(4,6) * t251 - Icges(4,3) * t351;
t339 = t318 * t363;
t273 = -t314 * t321 + t317 * t339;
t349 = t318 * t321;
t274 = t314 * t363 + t317 * t349;
t238 = Icges(3,5) * t274 + Icges(3,6) * t273 - Icges(3,3) * t351;
t373 = t203 + t238;
t352 = t314 * t315;
t204 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t352;
t275 = -t314 * t339 - t317 * t321;
t276 = -t314 * t349 + t317 * t363;
t239 = Icges(3,5) * t276 + Icges(3,6) * t275 + Icges(3,3) * t352;
t372 = t204 + t239;
t234 = Icges(4,5) * t270 + Icges(4,6) * t269 + Icges(4,3) * t318;
t265 = Icges(3,3) * t318 + (Icges(3,5) * t321 + Icges(3,6) * t363) * t315;
t371 = t234 + t265;
t361 = pkin(7) * t318;
t360 = pkin(2) * t363;
t359 = pkin(5) * t316;
t356 = Icges(2,4) * t314;
t347 = qJD(2) * t315;
t346 = qJD(3) * t315;
t345 = V_base(5) * qJ(1) + V_base(1);
t341 = qJD(1) + V_base(3);
t288 = t314 * t347 + V_base(4);
t299 = qJD(2) * t318 + V_base(6);
t338 = pkin(2) * t349 - qJ(3) * t315;
t228 = -qJD(4) * t253 + t288;
t257 = -qJD(4) * t269 + t299;
t287 = -t317 * t347 + V_base(5);
t227 = -qJD(4) * t251 + t287;
t282 = pkin(1) * t314 - pkin(7) * t351;
t335 = -t282 * V_base(6) + V_base(5) * t361 + t345;
t283 = pkin(1) * t317 + pkin(7) * t352;
t334 = V_base(4) * t282 - t283 * V_base(5) + t341;
t332 = V_base(6) * t283 + V_base(2) + (-qJ(1) - t361) * V_base(4);
t281 = pkin(2) * t315 * t321 + t318 * qJ(3);
t331 = t287 * t281 + t314 * t346 + t335;
t246 = t314 * t360 + t317 * t338;
t330 = qJD(3) * t318 + t288 * t246 + t334;
t247 = -t314 * t338 + t317 * t360;
t328 = t299 * t247 - t317 * t346 + t332;
t211 = pkin(3) * t252 - pkin(8) * t251;
t250 = pkin(3) * t270 - pkin(8) * t269;
t327 = t287 * t250 + (-t211 - t246) * t299 + t331;
t212 = pkin(3) * t254 - pkin(8) * t253;
t326 = t288 * t211 + (-t212 - t247) * t287 + t330;
t217 = pkin(4) * t259 + qJ(5) * t258;
t325 = qJD(5) * t231 + t227 * t217 + t327;
t189 = pkin(4) * t230 + qJ(5) * t229;
t324 = qJD(5) * t258 + t228 * t189 + t326;
t323 = t299 * t212 + (-t250 - t281) * t288 + t328;
t190 = pkin(4) * t232 + qJ(5) * t231;
t322 = qJD(5) * t229 + t257 * t190 + t323;
t312 = pkin(12) + qJ(6);
t310 = Icges(2,4) * t317;
t309 = cos(t312);
t308 = sin(t312);
t296 = rSges(2,1) * t317 - rSges(2,2) * t314;
t295 = rSges(2,1) * t314 + rSges(2,2) * t317;
t294 = Icges(2,1) * t317 - t356;
t293 = Icges(2,1) * t314 + t310;
t292 = -Icges(2,2) * t314 + t310;
t291 = Icges(2,2) * t317 + t356;
t286 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t285 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t284 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t268 = t318 * rSges(3,3) + (rSges(3,1) * t321 + rSges(3,2) * t363) * t315;
t267 = Icges(3,5) * t318 + (Icges(3,1) * t321 + Icges(3,4) * t363) * t315;
t266 = Icges(3,6) * t318 + (Icges(3,4) * t321 + Icges(3,2) * t363) * t315;
t261 = V_base(5) * rSges(2,3) - t295 * V_base(6) + t345;
t260 = t296 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t256 = t295 * V_base(4) - t296 * V_base(5) + t341;
t245 = rSges(3,1) * t276 + rSges(3,2) * t275 + rSges(3,3) * t352;
t244 = rSges(3,1) * t274 + rSges(3,2) * t273 - rSges(3,3) * t351;
t243 = Icges(3,1) * t276 + Icges(3,4) * t275 + Icges(3,5) * t352;
t242 = Icges(3,1) * t274 + Icges(3,4) * t273 - Icges(3,5) * t351;
t241 = Icges(3,4) * t276 + Icges(3,2) * t275 + Icges(3,6) * t352;
t240 = Icges(3,4) * t274 + Icges(3,2) * t273 - Icges(3,6) * t351;
t237 = rSges(4,1) * t270 + rSges(4,2) * t269 + rSges(4,3) * t318;
t236 = Icges(4,1) * t270 + Icges(4,4) * t269 + Icges(4,5) * t318;
t235 = Icges(4,4) * t270 + Icges(4,2) * t269 + Icges(4,6) * t318;
t220 = t259 * t309 - t269 * t308;
t219 = -t259 * t308 - t269 * t309;
t218 = qJD(6) * t258 + t257;
t216 = rSges(5,1) * t259 - rSges(5,2) * t258 - rSges(5,3) * t269;
t215 = Icges(5,1) * t259 - Icges(5,4) * t258 - Icges(5,5) * t269;
t213 = Icges(5,5) * t259 - Icges(5,6) * t258 - Icges(5,3) * t269;
t210 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t352;
t209 = rSges(4,1) * t252 + rSges(4,2) * t251 - rSges(4,3) * t351;
t208 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t352;
t207 = Icges(4,1) * t252 + Icges(4,4) * t251 - Icges(4,5) * t351;
t206 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t352;
t205 = Icges(4,4) * t252 + Icges(4,2) * t251 - Icges(4,6) * t351;
t196 = t232 * t309 - t253 * t308;
t195 = -t232 * t308 - t253 * t309;
t194 = t230 * t309 - t251 * t308;
t193 = -t230 * t308 - t251 * t309;
t192 = qJD(6) * t231 + t228;
t191 = qJD(6) * t229 + t227;
t188 = -t244 * t299 + t268 * t287 + t335;
t187 = t245 * t299 - t268 * t288 + t332;
t184 = rSges(6,1) * t222 + rSges(6,2) * t221 + rSges(6,3) * t258;
t183 = t244 * t288 - t245 * t287 + t334;
t182 = Icges(6,1) * t222 + Icges(6,4) * t221 + Icges(6,5) * t258;
t181 = Icges(6,4) * t222 + Icges(6,2) * t221 + Icges(6,6) * t258;
t179 = rSges(5,1) * t232 - rSges(5,2) * t231 - rSges(5,3) * t253;
t178 = rSges(5,1) * t230 - rSges(5,2) * t229 - rSges(5,3) * t251;
t177 = Icges(5,1) * t232 - Icges(5,4) * t231 - Icges(5,5) * t253;
t176 = Icges(5,1) * t230 - Icges(5,4) * t229 - Icges(5,5) * t251;
t173 = Icges(5,5) * t232 - Icges(5,6) * t231 - Icges(5,3) * t253;
t172 = Icges(5,5) * t230 - Icges(5,6) * t229 - Icges(5,3) * t251;
t171 = -pkin(5) * t353 + pkin(9) * t258 + t259 * t359;
t170 = rSges(7,1) * t220 + rSges(7,2) * t219 + rSges(7,3) * t258;
t169 = Icges(7,1) * t220 + Icges(7,4) * t219 + Icges(7,5) * t258;
t168 = Icges(7,4) * t220 + Icges(7,2) * t219 + Icges(7,6) * t258;
t167 = Icges(7,5) * t220 + Icges(7,6) * t219 + Icges(7,3) * t258;
t165 = rSges(6,1) * t200 + rSges(6,2) * t199 + rSges(6,3) * t231;
t164 = rSges(6,1) * t198 + rSges(6,2) * t197 + rSges(6,3) * t229;
t163 = Icges(6,1) * t200 + Icges(6,4) * t199 + Icges(6,5) * t231;
t162 = Icges(6,1) * t198 + Icges(6,4) * t197 + Icges(6,5) * t229;
t161 = Icges(6,4) * t200 + Icges(6,2) * t199 + Icges(6,6) * t231;
t160 = Icges(6,4) * t198 + Icges(6,2) * t197 + Icges(6,6) * t229;
t157 = rSges(7,1) * t196 + rSges(7,2) * t195 + rSges(7,3) * t231;
t156 = rSges(7,1) * t194 + rSges(7,2) * t193 + rSges(7,3) * t229;
t155 = Icges(7,1) * t196 + Icges(7,4) * t195 + Icges(7,5) * t231;
t154 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t229;
t153 = Icges(7,4) * t196 + Icges(7,2) * t195 + Icges(7,6) * t231;
t152 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t229;
t151 = Icges(7,5) * t196 + Icges(7,6) * t195 + Icges(7,3) * t231;
t150 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t229;
t149 = -pkin(5) * t354 + pkin(9) * t231 + t232 * t359;
t148 = -pkin(5) * t355 + pkin(9) * t229 + t230 * t359;
t147 = t237 * t287 + (-t209 - t246) * t299 + t331;
t146 = t210 * t299 + (-t237 - t281) * t288 + t328;
t145 = t209 * t288 + (-t210 - t247) * t287 + t330;
t144 = -t178 * t257 + t216 * t227 + t327;
t143 = t179 * t257 - t216 * t228 + t323;
t142 = t178 * t228 - t179 * t227 + t326;
t141 = t184 * t227 + (-t164 - t189) * t257 + t325;
t140 = t165 * t257 + (-t184 - t217) * t228 + t322;
t139 = t164 * t228 + (-t165 - t190) * t227 + t324;
t138 = -t156 * t218 + t170 * t191 + t171 * t227 + (-t148 - t189) * t257 + t325;
t137 = t149 * t257 + t157 * t218 - t170 * t192 + (-t171 - t217) * t228 + t322;
t136 = t324 + (-t149 - t190) * t227 + t148 * t228 + t156 * t192 - t157 * t191;
t1 = m(3) * (t183 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(4) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(6) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(5) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(7) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + t191 * ((t151 * t229 + t153 * t193 + t155 * t194) * t192 + (t150 * t229 + t152 * t193 + t154 * t194) * t191 + (t167 * t229 + t168 * t193 + t169 * t194) * t218) / 0.2e1 + m(2) * (t256 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t218 * ((t151 * t258 + t153 * t219 + t155 * t220) * t192 + (t150 * t258 + t152 * t219 + t154 * t220) * t191 + (t167 * t258 + t168 * t219 + t169 * t220) * t218) / 0.2e1 + m(1) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + t192 * ((t151 * t231 + t153 * t195 + t155 * t196) * t192 + (t150 * t231 + t152 * t195 + t154 * t196) * t191 + (t167 * t231 + t168 * t195 + t169 * t196) * t218) / 0.2e1 + ((t181 * t197 + t182 * t198 - t213 * t251 + t215 * t230 + t229 * t374) * t257 + (t161 * t197 + t163 * t198 - t173 * t251 + t177 * t230 + t229 * t375) * t228 + (t160 * t197 + t162 * t198 - t172 * t251 + t176 * t230 + t376 * t229) * t227) * t227 / 0.2e1 + ((t181 * t199 + t182 * t200 - t213 * t253 + t215 * t232 + t231 * t374) * t257 + (t161 * t199 + t163 * t200 - t173 * t253 + t177 * t232 + t375 * t231) * t228 + (t160 * t199 + t162 * t200 - t172 * t253 + t176 * t232 + t231 * t376) * t227) * t228 / 0.2e1 + ((t181 * t221 + t182 * t222 - t213 * t269 + t215 * t259 + t258 * t374) * t257 + (t161 * t221 + t163 * t222 - t173 * t269 + t177 * t259 + t258 * t375) * t228 + (t160 * t221 + t162 * t222 - t172 * t269 + t176 * t259 + t258 * t376) * t227) * t257 / 0.2e1 + ((t235 * t251 + t236 * t252 + t266 * t273 + t267 * t274 - t351 * t371) * t299 + (t206 * t251 + t208 * t252 + t241 * t273 + t243 * t274 - t351 * t372) * t288 + (t205 * t251 + t207 * t252 + t240 * t273 + t242 * t274 - t351 * t373) * t287) * t287 / 0.2e1 + ((t235 * t253 + t236 * t254 + t266 * t275 + t267 * t276 + t352 * t371) * t299 + (t206 * t253 + t208 * t254 + t241 * t275 + t243 * t276 + t352 * t372) * t288 + (t205 * t253 + t207 * t254 + t240 * t275 + t242 * t276 + t373 * t352) * t287) * t288 / 0.2e1 + ((t238 * t287 + t239 * t288 + t265 * t299) * t318 + ((t241 * t363 + t243 * t321) * t288 + (t240 * t363 + t242 * t321) * t287 + (t266 * t363 + t267 * t321) * t299) * t315 + (t204 * t318 + t206 * t269 + t208 * t270) * t288 + (t203 * t318 + t205 * t269 + t207 * t270) * t287 + (t234 * t318 + t235 * t269 + t236 * t270) * t299) * t299 / 0.2e1 + ((-t291 * t314 + t293 * t317 + Icges(1,4)) * V_base(5) + (-t292 * t314 + t294 * t317 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t291 * t317 + t293 * t314 + Icges(1,2)) * V_base(5) + (t292 * t317 + t294 * t314 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t314 + Icges(2,6) * t317 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t317 - Icges(2,6) * t314 + Icges(1,5)) * V_base(4) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
