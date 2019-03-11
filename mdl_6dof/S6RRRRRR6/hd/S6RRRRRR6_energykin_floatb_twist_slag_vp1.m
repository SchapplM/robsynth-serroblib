% Calculate kinetic energy for
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:35
% EndTime: 2019-03-10 04:11:39
% DurationCPUTime: 4.03s
% Computational Cost: add. (3855->443), mult. (6268->661), div. (0->0), fcn. (7460->14), ass. (0->197)
t315 = cos(pkin(6));
t365 = pkin(8) * t315;
t321 = cos(qJ(3));
t364 = pkin(3) * t321;
t320 = cos(qJ(5));
t363 = pkin(5) * t320;
t319 = sin(qJ(1));
t360 = Icges(2,4) * t319;
t322 = cos(qJ(2));
t323 = cos(qJ(1));
t348 = t322 * t323;
t318 = sin(qJ(2));
t351 = t318 * t319;
t277 = -t315 * t348 + t351;
t316 = sin(qJ(5));
t359 = t277 * t316;
t349 = t319 * t322;
t350 = t318 * t323;
t279 = t315 * t349 + t350;
t358 = t279 * t316;
t314 = sin(pkin(6));
t357 = t314 * t318;
t356 = t314 * t319;
t355 = t314 * t321;
t354 = t314 * t322;
t353 = t314 * t323;
t317 = sin(qJ(3));
t352 = t315 * t317;
t347 = qJ(3) + qJ(4);
t346 = qJD(2) * t314;
t345 = V_base(5) * pkin(7) + V_base(1);
t342 = t316 * t354;
t341 = t317 * t356;
t340 = t317 * t353;
t289 = t319 * t346 + V_base(4);
t339 = cos(t347);
t307 = V_base(6) + qJD(1);
t252 = qJD(3) * t279 + t289;
t291 = qJD(2) * t315 + t307;
t338 = t314 * t339;
t223 = qJD(4) * t279 + t252;
t288 = -t323 * t346 + V_base(5);
t280 = -t315 * t351 + t348;
t309 = sin(t347);
t249 = t280 * t309 - t319 * t338;
t195 = qJD(5) * t249 + t223;
t283 = pkin(1) * t319 - pkin(8) * t353;
t337 = -t283 * t307 + t365 * V_base(5) + t345;
t284 = pkin(1) * t323 + pkin(8) * t356;
t336 = V_base(4) * t283 - t284 * V_base(5) + V_base(3);
t251 = qJD(3) * t277 + t288;
t222 = qJD(4) * t277 + t251;
t335 = t307 * t284 + V_base(2) + (-pkin(7) - t365) * V_base(4);
t278 = t315 * t350 + t349;
t247 = t278 * t309 + t323 * t338;
t194 = qJD(5) * t247 + t222;
t261 = (-qJD(3) - qJD(4)) * t354 + t291;
t241 = t278 * pkin(2) + t277 * pkin(9);
t282 = (pkin(2) * t318 - pkin(9) * t322) * t314;
t334 = -t241 * t291 + t288 * t282 + t337;
t242 = t280 * pkin(2) + t279 * pkin(9);
t333 = t289 * t241 - t242 * t288 + t336;
t267 = t309 * t357 - t315 * t339;
t220 = qJD(5) * t267 + t261;
t332 = t291 * t242 - t282 * t289 + t335;
t192 = -pkin(3) * t340 + pkin(10) * t277 + t278 * t364;
t236 = pkin(3) * t352 + (-pkin(10) * t322 + t318 * t364) * t314;
t273 = -qJD(3) * t354 + t291;
t331 = -t192 * t273 + t251 * t236 + t334;
t193 = pkin(3) * t341 + pkin(10) * t279 + t280 * t364;
t330 = t252 * t192 - t193 * t251 + t333;
t329 = t273 * t193 - t236 * t252 + t332;
t248 = t278 * t339 - t309 * t353;
t206 = pkin(4) * t248 + pkin(11) * t247;
t268 = t315 * t309 + t318 * t338;
t229 = pkin(4) * t268 + pkin(11) * t267;
t328 = -t206 * t261 + t222 * t229 + t331;
t250 = t280 * t339 + t309 * t356;
t207 = pkin(4) * t250 + pkin(11) * t249;
t327 = t223 * t206 - t207 * t222 + t330;
t326 = t261 * t207 - t223 * t229 + t329;
t313 = qJ(5) + qJ(6);
t311 = Icges(2,4) * t323;
t310 = cos(t313);
t308 = sin(t313);
t299 = rSges(2,1) * t323 - rSges(2,2) * t319;
t298 = rSges(2,1) * t319 + rSges(2,2) * t323;
t297 = Icges(2,1) * t323 - t360;
t296 = Icges(2,1) * t319 + t311;
t295 = -Icges(2,2) * t319 + t311;
t294 = Icges(2,2) * t323 + t360;
t287 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t286 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t285 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t276 = t318 * t355 + t352;
t275 = t315 * t321 - t317 * t357;
t266 = rSges(3,3) * t315 + (rSges(3,1) * t318 + rSges(3,2) * t322) * t314;
t265 = Icges(3,5) * t315 + (Icges(3,1) * t318 + Icges(3,4) * t322) * t314;
t264 = Icges(3,6) * t315 + (Icges(3,4) * t318 + Icges(3,2) * t322) * t314;
t263 = Icges(3,3) * t315 + (Icges(3,5) * t318 + Icges(3,6) * t322) * t314;
t260 = V_base(5) * rSges(2,3) - t298 * t307 + t345;
t259 = t299 * t307 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t257 = t298 * V_base(4) - t299 * V_base(5) + V_base(3);
t256 = t280 * t321 + t341;
t255 = -t280 * t317 + t319 * t355;
t254 = t278 * t321 - t340;
t253 = -t278 * t317 - t321 * t353;
t246 = t268 * t320 - t342;
t245 = -t268 * t316 - t320 * t354;
t240 = t268 * t310 - t308 * t354;
t239 = -t268 * t308 - t310 * t354;
t238 = rSges(3,1) * t280 - rSges(3,2) * t279 + rSges(3,3) * t356;
t237 = rSges(3,1) * t278 - rSges(3,2) * t277 - rSges(3,3) * t353;
t235 = Icges(3,1) * t280 - Icges(3,4) * t279 + Icges(3,5) * t356;
t234 = Icges(3,1) * t278 - Icges(3,4) * t277 - Icges(3,5) * t353;
t233 = Icges(3,4) * t280 - Icges(3,2) * t279 + Icges(3,6) * t356;
t232 = Icges(3,4) * t278 - Icges(3,2) * t277 - Icges(3,6) * t353;
t231 = Icges(3,5) * t280 - Icges(3,6) * t279 + Icges(3,3) * t356;
t230 = Icges(3,5) * t278 - Icges(3,6) * t277 - Icges(3,3) * t353;
t228 = rSges(4,1) * t276 + rSges(4,2) * t275 - rSges(4,3) * t354;
t227 = Icges(4,1) * t276 + Icges(4,4) * t275 - Icges(4,5) * t354;
t226 = Icges(4,4) * t276 + Icges(4,2) * t275 - Icges(4,6) * t354;
t225 = Icges(4,5) * t276 + Icges(4,6) * t275 - Icges(4,3) * t354;
t219 = rSges(5,1) * t268 - rSges(5,2) * t267 - rSges(5,3) * t354;
t218 = Icges(5,1) * t268 - Icges(5,4) * t267 - Icges(5,5) * t354;
t217 = Icges(5,4) * t268 - Icges(5,2) * t267 - Icges(5,6) * t354;
t216 = Icges(5,5) * t268 - Icges(5,6) * t267 - Icges(5,3) * t354;
t215 = t250 * t320 + t358;
t214 = -t250 * t316 + t279 * t320;
t213 = t248 * t320 + t359;
t212 = -t248 * t316 + t277 * t320;
t211 = t250 * t310 + t279 * t308;
t210 = -t250 * t308 + t279 * t310;
t209 = t248 * t310 + t277 * t308;
t208 = -t248 * t308 + t277 * t310;
t205 = qJD(6) * t267 + t220;
t203 = rSges(4,1) * t256 + rSges(4,2) * t255 + rSges(4,3) * t279;
t202 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t277;
t201 = Icges(4,1) * t256 + Icges(4,4) * t255 + Icges(4,5) * t279;
t200 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t277;
t199 = Icges(4,4) * t256 + Icges(4,2) * t255 + Icges(4,6) * t279;
t198 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t277;
t197 = Icges(4,5) * t256 + Icges(4,6) * t255 + Icges(4,3) * t279;
t196 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t277;
t191 = rSges(5,1) * t250 - rSges(5,2) * t249 + rSges(5,3) * t279;
t190 = rSges(5,1) * t248 - rSges(5,2) * t247 + rSges(5,3) * t277;
t189 = Icges(5,1) * t250 - Icges(5,4) * t249 + Icges(5,5) * t279;
t188 = Icges(5,1) * t248 - Icges(5,4) * t247 + Icges(5,5) * t277;
t187 = Icges(5,4) * t250 - Icges(5,2) * t249 + Icges(5,6) * t279;
t186 = Icges(5,4) * t248 - Icges(5,2) * t247 + Icges(5,6) * t277;
t185 = Icges(5,5) * t250 - Icges(5,6) * t249 + Icges(5,3) * t279;
t184 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t277;
t182 = rSges(6,1) * t246 + rSges(6,2) * t245 + rSges(6,3) * t267;
t181 = Icges(6,1) * t246 + Icges(6,4) * t245 + Icges(6,5) * t267;
t180 = Icges(6,4) * t246 + Icges(6,2) * t245 + Icges(6,6) * t267;
t179 = Icges(6,5) * t246 + Icges(6,6) * t245 + Icges(6,3) * t267;
t177 = rSges(7,1) * t240 + rSges(7,2) * t239 + rSges(7,3) * t267;
t176 = Icges(7,1) * t240 + Icges(7,4) * t239 + Icges(7,5) * t267;
t175 = Icges(7,4) * t240 + Icges(7,2) * t239 + Icges(7,6) * t267;
t174 = Icges(7,5) * t240 + Icges(7,6) * t239 + Icges(7,3) * t267;
t173 = -pkin(5) * t342 + pkin(12) * t267 + t268 * t363;
t170 = qJD(6) * t249 + t195;
t169 = qJD(6) * t247 + t194;
t167 = -t237 * t291 + t266 * t288 + t337;
t166 = t238 * t291 - t266 * t289 + t335;
t165 = t237 * t289 - t238 * t288 + t336;
t164 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t249;
t163 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t247;
t162 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t249;
t161 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t247;
t160 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t249;
t159 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t247;
t158 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t249;
t157 = Icges(6,5) * t213 + Icges(6,6) * t212 + Icges(6,3) * t247;
t156 = rSges(7,1) * t211 + rSges(7,2) * t210 + rSges(7,3) * t249;
t155 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t247;
t154 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t249;
t153 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t247;
t152 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t249;
t151 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t247;
t150 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t249;
t149 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t247;
t148 = pkin(5) * t358 + pkin(12) * t249 + t250 * t363;
t147 = pkin(5) * t359 + pkin(12) * t247 + t248 * t363;
t146 = -t202 * t273 + t228 * t251 + t334;
t145 = t203 * t273 - t228 * t252 + t332;
t144 = t202 * t252 - t203 * t251 + t333;
t143 = -t190 * t261 + t219 * t222 + t331;
t142 = t191 * t261 - t219 * t223 + t329;
t141 = t190 * t223 - t191 * t222 + t330;
t140 = -t163 * t220 + t182 * t194 + t328;
t139 = t164 * t220 - t182 * t195 + t326;
t138 = t163 * t195 - t164 * t194 + t327;
t137 = -t147 * t220 - t155 * t205 + t169 * t177 + t173 * t194 + t328;
t136 = t148 * t220 + t156 * t205 - t170 * t177 - t173 * t195 + t326;
t135 = t147 * t195 - t148 * t194 + t155 * t170 - t156 * t169 + t327;
t1 = t273 * ((-t197 * t354 + t199 * t275 + t201 * t276) * t252 + (-t196 * t354 + t198 * t275 + t200 * t276) * t251 + (-t225 * t354 + t226 * t275 + t227 * t276) * t273) / 0.2e1 + t261 * ((-t185 * t354 - t187 * t267 + t189 * t268) * t223 + (-t184 * t354 - t186 * t267 + t188 * t268) * t222 + (-t216 * t354 - t217 * t267 + t218 * t268) * t261) / 0.2e1 + t289 * ((t231 * t356 - t233 * t279 + t235 * t280) * t289 + (t230 * t356 - t232 * t279 + t234 * t280) * t288 + (t263 * t356 - t264 * t279 + t265 * t280) * t291) / 0.2e1 + t288 * ((-t231 * t353 - t233 * t277 + t235 * t278) * t289 + (-t230 * t353 - t232 * t277 + t234 * t278) * t288 + (-t263 * t353 - t264 * t277 + t265 * t278) * t291) / 0.2e1 + ((t294 * t323 + t296 * t319 + Icges(1,2)) * V_base(5) + (t295 * t323 + t297 * t319 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t291 * ((t230 * t288 + t231 * t289 + t263 * t291) * t315 + ((t233 * t322 + t235 * t318) * t289 + (t232 * t322 + t234 * t318) * t288 + (t264 * t322 + t265 * t318) * t291) * t314) / 0.2e1 + ((-t294 * t319 + t296 * t323 + Icges(1,4)) * V_base(5) + (-t295 * t319 + t297 * t323 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((Icges(2,5) * t319 + Icges(2,6) * t323) * V_base(5) + (Icges(2,5) * t323 - Icges(2,6) * t319) * V_base(4) + Icges(2,3) * t307 / 0.2e1) * t307 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + t169 * ((t150 * t247 + t152 * t208 + t154 * t209) * t170 + (t247 * t149 + t208 * t151 + t209 * t153) * t169 + (t174 * t247 + t175 * t208 + t176 * t209) * t205) / 0.2e1 + t194 * ((t158 * t247 + t160 * t212 + t162 * t213) * t195 + (t247 * t157 + t212 * t159 + t213 * t161) * t194 + (t179 * t247 + t180 * t212 + t181 * t213) * t220) / 0.2e1 + t170 * ((t249 * t150 + t210 * t152 + t211 * t154) * t170 + (t149 * t249 + t151 * t210 + t153 * t211) * t169 + (t174 * t249 + t175 * t210 + t176 * t211) * t205) / 0.2e1 + t195 * ((t249 * t158 + t214 * t160 + t215 * t162) * t195 + (t157 * t249 + t159 * t214 + t161 * t215) * t194 + (t179 * t249 + t180 * t214 + t181 * t215) * t220) / 0.2e1 + m(2) * (t257 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t205 * ((t150 * t267 + t152 * t239 + t154 * t240) * t170 + (t149 * t267 + t151 * t239 + t153 * t240) * t169 + (t267 * t174 + t239 * t175 + t240 * t176) * t205) / 0.2e1 + t220 * ((t158 * t267 + t160 * t245 + t162 * t246) * t195 + (t157 * t267 + t159 * t245 + t161 * t246) * t194 + (t179 * t267 + t180 * t245 + t181 * t246) * t220) / 0.2e1 + t222 * ((t185 * t277 - t187 * t247 + t189 * t248) * t223 + (t184 * t277 - t186 * t247 + t188 * t248) * t222 + (t216 * t277 - t217 * t247 + t218 * t248) * t261) / 0.2e1 + t251 * ((t197 * t277 + t199 * t253 + t201 * t254) * t252 + (t196 * t277 + t198 * t253 + t200 * t254) * t251 + (t225 * t277 + t226 * t253 + t227 * t254) * t273) / 0.2e1 + t223 * ((t185 * t279 - t187 * t249 + t189 * t250) * t223 + (t184 * t279 - t186 * t249 + t188 * t250) * t222 + (t216 * t279 - t217 * t249 + t218 * t250) * t261) / 0.2e1 + t252 * ((t197 * t279 + t199 * t255 + t201 * t256) * t252 + (t196 * t279 + t198 * t255 + t200 * t256) * t251 + (t225 * t279 + t226 * t255 + t227 * t256) * t273) / 0.2e1 + m(1) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1;
T  = t1;
