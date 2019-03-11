% Calculate kinetic energy for
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:35
% EndTime: 2019-03-08 20:26:39
% DurationCPUTime: 4.39s
% Computational Cost: add. (3981->439), mult. (9164->639), div. (0->0), fcn. (11536->14), ass. (0->195)
t315 = sin(pkin(11));
t317 = cos(pkin(11));
t321 = sin(qJ(2));
t358 = sin(pkin(12));
t359 = cos(pkin(12));
t365 = cos(qJ(2));
t282 = -t321 * t358 + t359 * t365;
t318 = cos(pkin(6));
t331 = t318 * t282;
t335 = t321 * t359 + t358 * t365;
t253 = -t315 * t335 + t317 * t331;
t273 = t335 * t318;
t254 = t273 * t317 + t282 * t315;
t316 = sin(pkin(6));
t352 = t316 * t317;
t204 = Icges(4,5) * t254 + Icges(4,6) * t253 - Icges(4,3) * t352;
t341 = t318 * t365;
t275 = -t315 * t321 + t317 * t341;
t350 = t318 * t321;
t276 = t315 * t365 + t317 * t350;
t240 = Icges(3,5) * t276 + Icges(3,6) * t275 - Icges(3,3) * t352;
t372 = t204 + t240;
t255 = -t315 * t331 - t317 * t335;
t256 = -t273 * t315 + t282 * t317;
t353 = t315 * t316;
t205 = Icges(4,5) * t256 + Icges(4,6) * t255 + Icges(4,3) * t353;
t277 = -t315 * t341 - t317 * t321;
t278 = -t315 * t350 + t317 * t365;
t241 = Icges(3,5) * t278 + Icges(3,6) * t277 + Icges(3,3) * t353;
t371 = t205 + t241;
t271 = t282 * t316;
t272 = t335 * t316;
t236 = Icges(4,5) * t272 + Icges(4,6) * t271 + Icges(4,3) * t318;
t267 = Icges(3,3) * t318 + (Icges(3,5) * t321 + Icges(3,6) * t365) * t316;
t370 = t236 + t267;
t364 = cos(qJ(4));
t363 = pkin(7) * t318;
t362 = pkin(2) * t365;
t322 = cos(qJ(5));
t361 = pkin(5) * t322;
t357 = Icges(2,4) * t315;
t319 = sin(qJ(5));
t356 = t253 * t319;
t355 = t255 * t319;
t354 = t271 * t319;
t320 = sin(qJ(4));
t351 = t316 * t320;
t349 = qJD(2) * t316;
t348 = qJD(3) * t316;
t347 = V_base(5) * qJ(1) + V_base(1);
t343 = qJD(1) + V_base(3);
t342 = t316 * t364;
t290 = t315 * t349 + V_base(4);
t301 = qJD(2) * t318 + V_base(6);
t340 = pkin(2) * t350 - qJ(3) * t316;
t230 = -qJD(4) * t255 + t290;
t259 = -qJD(4) * t271 + t301;
t233 = t256 * t320 - t315 * t342;
t193 = qJD(5) * t233 + t230;
t260 = t272 * t320 - t318 * t364;
t220 = qJD(5) * t260 + t259;
t289 = -t317 * t349 + V_base(5);
t229 = -qJD(4) * t253 + t289;
t284 = pkin(1) * t315 - pkin(7) * t352;
t337 = -t284 * V_base(6) + t363 * V_base(5) + t347;
t285 = pkin(1) * t317 + pkin(7) * t353;
t336 = t284 * V_base(4) - t285 * V_base(5) + t343;
t231 = t254 * t320 + t317 * t342;
t192 = qJD(5) * t231 + t229;
t334 = V_base(6) * t285 + V_base(2) + (-qJ(1) - t363) * V_base(4);
t283 = pkin(2) * t316 * t321 + qJ(3) * t318;
t333 = t283 * t289 + t315 * t348 + t337;
t248 = t315 * t362 + t317 * t340;
t332 = qJD(3) * t318 + t248 * t290 + t336;
t249 = -t315 * t340 + t317 * t362;
t330 = t249 * t301 - t317 * t348 + t334;
t213 = pkin(3) * t254 - pkin(8) * t253;
t252 = pkin(3) * t272 - pkin(8) * t271;
t329 = t289 * t252 + (-t213 - t248) * t301 + t333;
t214 = pkin(3) * t256 - pkin(8) * t255;
t328 = t290 * t213 + (-t214 - t249) * t289 + t332;
t232 = t254 * t364 - t317 * t351;
t190 = pkin(4) * t232 + pkin(9) * t231;
t261 = t272 * t364 + t318 * t320;
t219 = pkin(4) * t261 + pkin(9) * t260;
t327 = -t190 * t259 + t219 * t229 + t329;
t234 = t256 * t364 + t315 * t351;
t191 = pkin(4) * t234 + pkin(9) * t233;
t326 = t190 * t230 - t191 * t229 + t328;
t325 = t301 * t214 + (-t252 - t283) * t290 + t330;
t324 = t191 * t259 - t219 * t230 + t325;
t314 = qJ(5) + qJ(6);
t312 = cos(t314);
t311 = sin(t314);
t310 = Icges(2,4) * t317;
t298 = rSges(2,1) * t317 - rSges(2,2) * t315;
t297 = rSges(2,1) * t315 + rSges(2,2) * t317;
t296 = Icges(2,1) * t317 - t357;
t295 = Icges(2,1) * t315 + t310;
t294 = -Icges(2,2) * t315 + t310;
t293 = Icges(2,2) * t317 + t357;
t288 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t287 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t286 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t270 = t318 * rSges(3,3) + (rSges(3,1) * t321 + rSges(3,2) * t365) * t316;
t269 = Icges(3,5) * t318 + (Icges(3,1) * t321 + Icges(3,4) * t365) * t316;
t268 = Icges(3,6) * t318 + (Icges(3,4) * t321 + Icges(3,2) * t365) * t316;
t263 = V_base(5) * rSges(2,3) - t297 * V_base(6) + t347;
t262 = t298 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t258 = t297 * V_base(4) - t298 * V_base(5) + t343;
t247 = rSges(3,1) * t278 + rSges(3,2) * t277 + rSges(3,3) * t353;
t246 = rSges(3,1) * t276 + rSges(3,2) * t275 - rSges(3,3) * t352;
t245 = Icges(3,1) * t278 + Icges(3,4) * t277 + Icges(3,5) * t353;
t244 = Icges(3,1) * t276 + Icges(3,4) * t275 - Icges(3,5) * t352;
t243 = Icges(3,4) * t278 + Icges(3,2) * t277 + Icges(3,6) * t353;
t242 = Icges(3,4) * t276 + Icges(3,2) * t275 - Icges(3,6) * t352;
t239 = rSges(4,1) * t272 + rSges(4,2) * t271 + rSges(4,3) * t318;
t238 = Icges(4,1) * t272 + Icges(4,4) * t271 + Icges(4,5) * t318;
t237 = Icges(4,4) * t272 + Icges(4,2) * t271 + Icges(4,6) * t318;
t224 = t261 * t322 - t354;
t223 = -t261 * t319 - t271 * t322;
t222 = t261 * t312 - t271 * t311;
t221 = -t261 * t311 - t271 * t312;
t218 = rSges(5,1) * t261 - rSges(5,2) * t260 - rSges(5,3) * t271;
t217 = Icges(5,1) * t261 - Icges(5,4) * t260 - Icges(5,5) * t271;
t216 = Icges(5,4) * t261 - Icges(5,2) * t260 - Icges(5,6) * t271;
t215 = Icges(5,5) * t261 - Icges(5,6) * t260 - Icges(5,3) * t271;
t212 = rSges(4,1) * t256 + rSges(4,2) * t255 + rSges(4,3) * t353;
t211 = rSges(4,1) * t254 + rSges(4,2) * t253 - rSges(4,3) * t352;
t210 = qJD(6) * t260 + t220;
t209 = Icges(4,1) * t256 + Icges(4,4) * t255 + Icges(4,5) * t353;
t208 = Icges(4,1) * t254 + Icges(4,4) * t253 - Icges(4,5) * t352;
t207 = Icges(4,4) * t256 + Icges(4,2) * t255 + Icges(4,6) * t353;
t206 = Icges(4,4) * t254 + Icges(4,2) * t253 - Icges(4,6) * t352;
t201 = t234 * t322 - t355;
t200 = -t234 * t319 - t255 * t322;
t199 = t232 * t322 - t356;
t198 = -t232 * t319 - t253 * t322;
t197 = t234 * t312 - t255 * t311;
t196 = -t234 * t311 - t255 * t312;
t195 = t232 * t312 - t253 * t311;
t194 = -t232 * t311 - t253 * t312;
t189 = -t246 * t301 + t270 * t289 + t337;
t188 = t247 * t301 - t270 * t290 + t334;
t185 = rSges(6,1) * t224 + rSges(6,2) * t223 + rSges(6,3) * t260;
t184 = Icges(6,1) * t224 + Icges(6,4) * t223 + Icges(6,5) * t260;
t183 = Icges(6,4) * t224 + Icges(6,2) * t223 + Icges(6,6) * t260;
t182 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t260;
t181 = t246 * t290 - t247 * t289 + t336;
t180 = rSges(5,1) * t234 - rSges(5,2) * t233 - rSges(5,3) * t255;
t179 = rSges(5,1) * t232 - rSges(5,2) * t231 - rSges(5,3) * t253;
t178 = Icges(5,1) * t234 - Icges(5,4) * t233 - Icges(5,5) * t255;
t177 = Icges(5,1) * t232 - Icges(5,4) * t231 - Icges(5,5) * t253;
t176 = Icges(5,4) * t234 - Icges(5,2) * t233 - Icges(5,6) * t255;
t175 = Icges(5,4) * t232 - Icges(5,2) * t231 - Icges(5,6) * t253;
t174 = Icges(5,5) * t234 - Icges(5,6) * t233 - Icges(5,3) * t255;
t173 = Icges(5,5) * t232 - Icges(5,6) * t231 - Icges(5,3) * t253;
t172 = -pkin(5) * t354 + pkin(10) * t260 + t261 * t361;
t171 = rSges(7,1) * t222 + rSges(7,2) * t221 + rSges(7,3) * t260;
t170 = Icges(7,1) * t222 + Icges(7,4) * t221 + Icges(7,5) * t260;
t169 = Icges(7,4) * t222 + Icges(7,2) * t221 + Icges(7,6) * t260;
t168 = Icges(7,5) * t222 + Icges(7,6) * t221 + Icges(7,3) * t260;
t167 = qJD(6) * t233 + t193;
t166 = qJD(6) * t231 + t192;
t164 = rSges(6,1) * t201 + rSges(6,2) * t200 + rSges(6,3) * t233;
t163 = rSges(6,1) * t199 + rSges(6,2) * t198 + rSges(6,3) * t231;
t162 = Icges(6,1) * t201 + Icges(6,4) * t200 + Icges(6,5) * t233;
t161 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t231;
t160 = Icges(6,4) * t201 + Icges(6,2) * t200 + Icges(6,6) * t233;
t159 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t231;
t158 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t233;
t157 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t231;
t156 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t233;
t155 = rSges(7,1) * t195 + rSges(7,2) * t194 + rSges(7,3) * t231;
t154 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t233;
t153 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t231;
t152 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t233;
t151 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t231;
t150 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t233;
t149 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t231;
t148 = -pkin(5) * t355 + pkin(10) * t233 + t234 * t361;
t147 = -pkin(5) * t356 + pkin(10) * t231 + t232 * t361;
t146 = t239 * t289 + (-t211 - t248) * t301 + t333;
t145 = t212 * t301 + (-t239 - t283) * t290 + t330;
t144 = t211 * t290 + (-t212 - t249) * t289 + t332;
t143 = -t179 * t259 + t218 * t229 + t329;
t142 = t180 * t259 - t218 * t230 + t325;
t141 = t179 * t230 - t180 * t229 + t328;
t140 = -t163 * t220 + t185 * t192 + t327;
t139 = t164 * t220 - t185 * t193 + t324;
t138 = t163 * t193 - t164 * t192 + t326;
t137 = -t147 * t220 - t155 * t210 + t166 * t171 + t172 * t192 + t327;
t136 = t148 * t220 + t156 * t210 - t167 * t171 - t172 * t193 + t324;
t135 = t147 * t193 - t148 * t192 + t155 * t167 - t156 * t166 + t326;
t1 = m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t181 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + t166 * ((t150 * t231 + t152 * t194 + t154 * t195) * t167 + (t149 * t231 + t151 * t194 + t153 * t195) * t166 + (t168 * t231 + t169 * t194 + t170 * t195) * t210) / 0.2e1 + t192 * ((t158 * t231 + t160 * t198 + t162 * t199) * t193 + (t157 * t231 + t159 * t198 + t161 * t199) * t192 + (t182 * t231 + t183 * t198 + t184 * t199) * t220) / 0.2e1 + t167 * ((t150 * t233 + t152 * t196 + t154 * t197) * t167 + (t149 * t233 + t151 * t196 + t153 * t197) * t166 + (t168 * t233 + t169 * t196 + t170 * t197) * t210) / 0.2e1 + t193 * ((t158 * t233 + t160 * t200 + t162 * t201) * t193 + (t157 * t233 + t159 * t200 + t161 * t201) * t192 + (t182 * t233 + t183 * t200 + t184 * t201) * t220) / 0.2e1 + t229 * ((-t174 * t253 - t176 * t231 + t178 * t232) * t230 + (-t173 * t253 - t175 * t231 + t177 * t232) * t229 + (-t215 * t253 - t216 * t231 + t217 * t232) * t259) / 0.2e1 + t230 * ((-t174 * t255 - t176 * t233 + t178 * t234) * t230 + (-t173 * t255 - t175 * t233 + t177 * t234) * t229 + (-t215 * t255 - t216 * t233 + t217 * t234) * t259) / 0.2e1 + t210 * ((t150 * t260 + t152 * t221 + t154 * t222) * t167 + (t149 * t260 + t151 * t221 + t153 * t222) * t166 + (t168 * t260 + t169 * t221 + t170 * t222) * t210) / 0.2e1 + t220 * ((t158 * t260 + t160 * t223 + t162 * t224) * t193 + (t157 * t260 + t159 * t223 + t161 * t224) * t192 + (t182 * t260 + t183 * t223 + t184 * t224) * t220) / 0.2e1 + m(2) * (t258 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t259 * ((-t174 * t271 - t176 * t260 + t178 * t261) * t230 + (-t173 * t271 - t175 * t260 + t177 * t261) * t229 + (-t215 * t271 - t216 * t260 + t217 * t261) * t259) / 0.2e1 + m(1) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + ((t237 * t253 + t238 * t254 + t268 * t275 + t269 * t276 - t352 * t370) * t301 + (t207 * t253 + t209 * t254 + t243 * t275 + t245 * t276 - t352 * t371) * t290 + (t206 * t253 + t208 * t254 + t242 * t275 + t244 * t276 - t372 * t352) * t289) * t289 / 0.2e1 + ((t237 * t255 + t238 * t256 + t268 * t277 + t269 * t278 + t353 * t370) * t301 + (t207 * t255 + t209 * t256 + t243 * t277 + t245 * t278 + t371 * t353) * t290 + (t206 * t255 + t208 * t256 + t242 * t277 + t244 * t278 + t353 * t372) * t289) * t290 / 0.2e1 + ((t240 * t289 + t241 * t290 + t267 * t301) * t318 + ((t243 * t365 + t245 * t321) * t290 + (t242 * t365 + t244 * t321) * t289 + (t268 * t365 + t269 * t321) * t301) * t316 + (t205 * t318 + t207 * t271 + t209 * t272) * t290 + (t204 * t318 + t206 * t271 + t208 * t272) * t289 + (t236 * t318 + t237 * t271 + t238 * t272) * t301) * t301 / 0.2e1 + ((-t293 * t315 + t295 * t317 + Icges(1,4)) * V_base(5) + (-t294 * t315 + t296 * t317 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t293 * t317 + t295 * t315 + Icges(1,2)) * V_base(5) + (t294 * t317 + t296 * t315 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t317 - Icges(2,6) * t315 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t315 + Icges(2,6) * t317 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
