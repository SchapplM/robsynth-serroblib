% Calculate kinetic energy for
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:46
% EndTime: 2019-03-09 00:47:50
% DurationCPUTime: 4.02s
% Computational Cost: add. (3579->442), mult. (7098->658), div. (0->0), fcn. (8608->14), ass. (0->192)
t361 = cos(qJ(3));
t318 = cos(pkin(6));
t359 = pkin(7) * t318;
t322 = cos(qJ(4));
t358 = t322 * pkin(4);
t315 = sin(pkin(12));
t356 = Icges(2,4) * t315;
t317 = cos(pkin(12));
t321 = sin(qJ(2));
t323 = cos(qJ(2));
t348 = t318 * t323;
t270 = t315 * t321 - t317 * t348;
t319 = sin(qJ(4));
t355 = t270 * t319;
t272 = t315 * t348 + t317 * t321;
t354 = t272 * t319;
t316 = sin(pkin(6));
t353 = t315 * t316;
t352 = t316 * t317;
t320 = sin(qJ(3));
t351 = t316 * t320;
t350 = t316 * t323;
t349 = t318 * t321;
t314 = qJ(4) + qJ(5);
t309 = cos(t314);
t347 = pkin(5) * t309;
t345 = qJD(2) * t316;
t344 = V_base(5) * qJ(1) + V_base(1);
t340 = qJD(1) + V_base(3);
t339 = t319 * t350;
t338 = t316 * t361;
t287 = t315 * t345 + V_base(4);
t299 = qJD(2) * t318 + V_base(6);
t308 = sin(t314);
t337 = pkin(5) * t308;
t252 = qJD(3) * t272 + t287;
t273 = -t315 * t349 + t317 * t323;
t255 = t273 * t320 - t315 * t338;
t210 = qJD(4) * t255 + t252;
t286 = -t317 * t345 + V_base(5);
t182 = qJD(5) * t255 + t210;
t251 = qJD(3) * t270 + t286;
t274 = -qJD(3) * t350 + t299;
t280 = pkin(1) * t315 - pkin(7) * t352;
t336 = -t280 * V_base(6) + V_base(5) * t359 + t344;
t281 = pkin(1) * t317 + pkin(7) * t353;
t335 = V_base(4) * t280 - t281 * V_base(5) + t340;
t271 = t315 * t323 + t317 * t349;
t253 = t271 * t320 + t317 * t338;
t209 = qJD(4) * t253 + t251;
t277 = -t318 * t361 + t321 * t351;
t243 = qJD(4) * t277 + t274;
t181 = qJD(5) * t253 + t209;
t223 = qJD(5) * t277 + t243;
t334 = V_base(6) * t281 + V_base(2) + (-qJ(1) - t359) * V_base(4);
t238 = pkin(2) * t271 + pkin(8) * t270;
t279 = (pkin(2) * t321 - pkin(8) * t323) * t316;
t333 = -t238 * t299 + t286 * t279 + t336;
t239 = pkin(2) * t273 + pkin(8) * t272;
t332 = t287 * t238 - t239 * t286 + t335;
t331 = t299 * t239 - t279 * t287 + t334;
t254 = t271 * t361 - t317 * t351;
t207 = t254 * pkin(3) + t253 * pkin(9);
t278 = t318 * t320 + t321 * t338;
t240 = t278 * pkin(3) + t277 * pkin(9);
t330 = -t207 * t274 + t251 * t240 + t333;
t256 = t273 * t361 + t315 * t351;
t208 = t256 * pkin(3) + t255 * pkin(9);
t329 = t252 * t207 - t208 * t251 + t332;
t328 = t274 * t208 - t240 * t252 + t331;
t155 = pkin(4) * t355 + pkin(10) * t253 + t254 * t358;
t191 = -pkin(4) * t339 + pkin(10) * t277 + t278 * t358;
t327 = -t155 * t243 + t209 * t191 + t330;
t156 = pkin(4) * t354 + pkin(10) * t255 + t256 * t358;
t326 = t210 * t155 - t156 * t209 + t329;
t325 = t243 * t156 - t191 * t210 + t328;
t311 = qJ(6) + t314;
t307 = Icges(2,4) * t317;
t304 = cos(t311);
t303 = sin(t311);
t296 = rSges(2,1) * t317 - rSges(2,2) * t315;
t295 = rSges(2,1) * t315 + rSges(2,2) * t317;
t294 = Icges(2,1) * t317 - t356;
t293 = Icges(2,1) * t315 + t307;
t292 = -Icges(2,2) * t315 + t307;
t291 = Icges(2,2) * t317 + t356;
t285 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t284 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t283 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t265 = rSges(3,3) * t318 + (rSges(3,1) * t321 + rSges(3,2) * t323) * t316;
t264 = Icges(3,5) * t318 + (Icges(3,1) * t321 + Icges(3,4) * t323) * t316;
t263 = Icges(3,6) * t318 + (Icges(3,4) * t321 + Icges(3,2) * t323) * t316;
t262 = Icges(3,3) * t318 + (Icges(3,5) * t321 + Icges(3,6) * t323) * t316;
t261 = V_base(5) * rSges(2,3) - t295 * V_base(6) + t344;
t260 = t296 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t258 = t278 * t322 - t339;
t257 = -t278 * t319 - t322 * t350;
t250 = t295 * V_base(4) - t296 * V_base(5) + t340;
t245 = t278 * t309 - t308 * t350;
t244 = -t278 * t308 - t309 * t350;
t242 = t278 * t304 - t303 * t350;
t241 = -t278 * t303 - t304 * t350;
t237 = rSges(4,1) * t278 - rSges(4,2) * t277 - rSges(4,3) * t350;
t236 = Icges(4,1) * t278 - Icges(4,4) * t277 - Icges(4,5) * t350;
t235 = Icges(4,4) * t278 - Icges(4,2) * t277 - Icges(4,6) * t350;
t234 = Icges(4,5) * t278 - Icges(4,6) * t277 - Icges(4,3) * t350;
t233 = rSges(3,1) * t273 - rSges(3,2) * t272 + rSges(3,3) * t353;
t232 = rSges(3,1) * t271 - rSges(3,2) * t270 - rSges(3,3) * t352;
t231 = Icges(3,1) * t273 - Icges(3,4) * t272 + Icges(3,5) * t353;
t230 = Icges(3,1) * t271 - Icges(3,4) * t270 - Icges(3,5) * t352;
t229 = Icges(3,4) * t273 - Icges(3,2) * t272 + Icges(3,6) * t353;
t228 = Icges(3,4) * t271 - Icges(3,2) * t270 - Icges(3,6) * t352;
t227 = Icges(3,5) * t273 - Icges(3,6) * t272 + Icges(3,3) * t353;
t226 = Icges(3,5) * t271 - Icges(3,6) * t270 - Icges(3,3) * t352;
t222 = t256 * t322 + t354;
t221 = -t256 * t319 + t272 * t322;
t220 = t254 * t322 + t355;
t219 = -t254 * t319 + t270 * t322;
t218 = t256 * t309 + t272 * t308;
t217 = -t256 * t308 + t272 * t309;
t216 = t254 * t309 + t270 * t308;
t215 = -t254 * t308 + t270 * t309;
t214 = t256 * t304 + t272 * t303;
t213 = -t256 * t303 + t272 * t304;
t212 = t254 * t304 + t270 * t303;
t211 = -t254 * t303 + t270 * t304;
t205 = qJD(6) * t277 + t223;
t204 = rSges(5,1) * t258 + rSges(5,2) * t257 + rSges(5,3) * t277;
t202 = Icges(5,1) * t258 + Icges(5,4) * t257 + Icges(5,5) * t277;
t201 = Icges(5,4) * t258 + Icges(5,2) * t257 + Icges(5,6) * t277;
t200 = Icges(5,5) * t258 + Icges(5,6) * t257 + Icges(5,3) * t277;
t199 = rSges(4,1) * t256 - rSges(4,2) * t255 + rSges(4,3) * t272;
t198 = rSges(4,1) * t254 - rSges(4,2) * t253 + rSges(4,3) * t270;
t197 = Icges(4,1) * t256 - Icges(4,4) * t255 + Icges(4,5) * t272;
t196 = Icges(4,1) * t254 - Icges(4,4) * t253 + Icges(4,5) * t270;
t195 = Icges(4,4) * t256 - Icges(4,2) * t255 + Icges(4,6) * t272;
t194 = Icges(4,4) * t254 - Icges(4,2) * t253 + Icges(4,6) * t270;
t193 = Icges(4,5) * t256 - Icges(4,6) * t255 + Icges(4,3) * t272;
t192 = Icges(4,5) * t254 - Icges(4,6) * t253 + Icges(4,3) * t270;
t190 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t277;
t189 = Icges(6,1) * t245 + Icges(6,4) * t244 + Icges(6,5) * t277;
t188 = Icges(6,4) * t245 + Icges(6,2) * t244 + Icges(6,6) * t277;
t187 = Icges(6,5) * t245 + Icges(6,6) * t244 + Icges(6,3) * t277;
t186 = rSges(7,1) * t242 + rSges(7,2) * t241 + rSges(7,3) * t277;
t185 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t277;
t184 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t277;
t183 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t277;
t179 = -t232 * t299 + t265 * t286 + t336;
t178 = t233 * t299 - t265 * t287 + t334;
t177 = pkin(11) * t277 + t278 * t347 - t337 * t350;
t176 = qJD(6) * t255 + t182;
t175 = qJD(6) * t253 + t181;
t174 = rSges(5,1) * t222 + rSges(5,2) * t221 + rSges(5,3) * t255;
t173 = rSges(5,1) * t220 + rSges(5,2) * t219 + rSges(5,3) * t253;
t172 = Icges(5,1) * t222 + Icges(5,4) * t221 + Icges(5,5) * t255;
t171 = Icges(5,1) * t220 + Icges(5,4) * t219 + Icges(5,5) * t253;
t170 = Icges(5,4) * t222 + Icges(5,2) * t221 + Icges(5,6) * t255;
t169 = Icges(5,4) * t220 + Icges(5,2) * t219 + Icges(5,6) * t253;
t168 = Icges(5,5) * t222 + Icges(5,6) * t221 + Icges(5,3) * t255;
t167 = Icges(5,5) * t220 + Icges(5,6) * t219 + Icges(5,3) * t253;
t165 = t232 * t287 - t233 * t286 + t335;
t164 = rSges(6,1) * t218 + rSges(6,2) * t217 + rSges(6,3) * t255;
t163 = rSges(6,1) * t216 + rSges(6,2) * t215 + rSges(6,3) * t253;
t162 = Icges(6,1) * t218 + Icges(6,4) * t217 + Icges(6,5) * t255;
t161 = Icges(6,1) * t216 + Icges(6,4) * t215 + Icges(6,5) * t253;
t160 = Icges(6,4) * t218 + Icges(6,2) * t217 + Icges(6,6) * t255;
t159 = Icges(6,4) * t216 + Icges(6,2) * t215 + Icges(6,6) * t253;
t158 = Icges(6,5) * t218 + Icges(6,6) * t217 + Icges(6,3) * t255;
t157 = Icges(6,5) * t216 + Icges(6,6) * t215 + Icges(6,3) * t253;
t154 = rSges(7,1) * t214 + rSges(7,2) * t213 + rSges(7,3) * t255;
t153 = rSges(7,1) * t212 + rSges(7,2) * t211 + rSges(7,3) * t253;
t152 = Icges(7,1) * t214 + Icges(7,4) * t213 + Icges(7,5) * t255;
t151 = Icges(7,1) * t212 + Icges(7,4) * t211 + Icges(7,5) * t253;
t150 = Icges(7,4) * t214 + Icges(7,2) * t213 + Icges(7,6) * t255;
t149 = Icges(7,4) * t212 + Icges(7,2) * t211 + Icges(7,6) * t253;
t148 = Icges(7,5) * t214 + Icges(7,6) * t213 + Icges(7,3) * t255;
t147 = Icges(7,5) * t212 + Icges(7,6) * t211 + Icges(7,3) * t253;
t145 = pkin(11) * t255 + t256 * t347 + t272 * t337;
t144 = pkin(11) * t253 + t254 * t347 + t270 * t337;
t142 = -t198 * t274 + t237 * t251 + t333;
t141 = t199 * t274 - t237 * t252 + t331;
t140 = t198 * t252 - t199 * t251 + t332;
t139 = -t173 * t243 + t204 * t209 + t330;
t138 = t174 * t243 - t204 * t210 + t328;
t137 = t173 * t210 - t174 * t209 + t329;
t136 = -t163 * t223 + t181 * t190 + t327;
t135 = t164 * t223 - t182 * t190 + t325;
t134 = t163 * t182 - t164 * t181 + t326;
t133 = -t144 * t223 - t153 * t205 + t175 * t186 + t177 * t181 + t327;
t132 = t145 * t223 + t154 * t205 - t176 * t186 - t177 * t182 + t325;
t131 = t144 * t182 - t145 * t181 + t153 * t176 - t154 * t175 + t326;
t1 = ((Icges(2,5) * t315 + Icges(2,6) * t317 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t317 - Icges(2,6) * t315 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((t291 * t317 + t293 * t315 + Icges(1,2)) * V_base(5) + (t292 * t317 + t294 * t315 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t287 * ((t227 * t353 - t229 * t272 + t231 * t273) * t287 + (t226 * t353 - t228 * t272 + t230 * t273) * t286 + (t262 * t353 - t263 * t272 + t264 * t273) * t299) / 0.2e1 + t286 * ((-t227 * t352 - t229 * t270 + t231 * t271) * t287 + (-t226 * t352 - t228 * t270 + t230 * t271) * t286 + (-t262 * t352 - t263 * t270 + t264 * t271) * t299) / 0.2e1 + t274 * ((-t193 * t350 - t195 * t277 + t197 * t278) * t252 + (-t192 * t350 - t194 * t277 + t196 * t278) * t251 + (-t234 * t350 - t235 * t277 + t236 * t278) * t274) / 0.2e1 + ((-t291 * t315 + t293 * t317 + Icges(1,4)) * V_base(5) + (-t292 * t315 + t294 * t317 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t299 * ((t226 * t286 + t227 * t287 + t262 * t299) * t318 + ((t229 * t323 + t231 * t321) * t287 + (t228 * t323 + t230 * t321) * t286 + (t263 * t323 + t264 * t321) * t299) * t316) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(3) * (t165 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + t175 * ((t148 * t253 + t150 * t211 + t152 * t212) * t176 + (t147 * t253 + t149 * t211 + t151 * t212) * t175 + (t183 * t253 + t184 * t211 + t185 * t212) * t205) / 0.2e1 + t181 * ((t158 * t253 + t160 * t215 + t162 * t216) * t182 + (t157 * t253 + t159 * t215 + t161 * t216) * t181 + (t187 * t253 + t188 * t215 + t189 * t216) * t223) / 0.2e1 + t209 * ((t168 * t253 + t170 * t219 + t172 * t220) * t210 + (t167 * t253 + t169 * t219 + t171 * t220) * t209 + (t200 * t253 + t201 * t219 + t202 * t220) * t243) / 0.2e1 + t176 * ((t148 * t255 + t150 * t213 + t152 * t214) * t176 + (t147 * t255 + t149 * t213 + t151 * t214) * t175 + (t183 * t255 + t184 * t213 + t185 * t214) * t205) / 0.2e1 + t182 * ((t158 * t255 + t160 * t217 + t162 * t218) * t182 + (t157 * t255 + t159 * t217 + t161 * t218) * t181 + (t187 * t255 + t188 * t217 + t189 * t218) * t223) / 0.2e1 + t210 * ((t168 * t255 + t170 * t221 + t172 * t222) * t210 + (t167 * t255 + t169 * t221 + t171 * t222) * t209 + (t200 * t255 + t201 * t221 + t202 * t222) * t243) / 0.2e1 + m(2) * (t250 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t251 * ((t193 * t270 - t195 * t253 + t197 * t254) * t252 + (t192 * t270 - t194 * t253 + t196 * t254) * t251 + (t234 * t270 - t235 * t253 + t236 * t254) * t274) / 0.2e1 + t252 * ((t193 * t272 - t195 * t255 + t197 * t256) * t252 + (t192 * t272 - t194 * t255 + t196 * t256) * t251 + (t234 * t272 - t235 * t255 + t236 * t256) * t274) / 0.2e1 + t205 * ((t148 * t277 + t150 * t241 + t152 * t242) * t176 + (t147 * t277 + t149 * t241 + t151 * t242) * t175 + (t183 * t277 + t184 * t241 + t185 * t242) * t205) / 0.2e1 + t223 * ((t158 * t277 + t160 * t244 + t162 * t245) * t182 + (t157 * t277 + t159 * t244 + t161 * t245) * t181 + (t187 * t277 + t188 * t244 + t189 * t245) * t223) / 0.2e1 + t243 * ((t168 * t277 + t170 * t257 + t172 * t258) * t210 + (t167 * t277 + t169 * t257 + t171 * t258) * t209 + (t200 * t277 + t201 * t257 + t202 * t258) * t243) / 0.2e1 + m(1) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1;
T  = t1;
