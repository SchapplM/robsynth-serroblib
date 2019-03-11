% Calculate kinetic energy for
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:29
% EndTime: 2019-03-09 10:55:33
% DurationCPUTime: 4.56s
% Computational Cost: add. (3527->437), mult. (5732->609), div. (0->0), fcn. (6800->14), ass. (0->191)
t371 = Icges(5,2) + Icges(6,3);
t313 = cos(pkin(6));
t317 = sin(qJ(1));
t318 = cos(qJ(2));
t345 = t317 * t318;
t316 = sin(qJ(2));
t319 = cos(qJ(1));
t346 = t316 * t319;
t272 = t313 * t346 + t345;
t340 = pkin(11) + qJ(4);
t302 = sin(t340);
t333 = cos(t340);
t310 = sin(pkin(6));
t348 = t310 * t319;
t243 = t272 * t333 - t302 * t348;
t344 = t318 * t319;
t347 = t316 * t317;
t271 = -t313 * t344 + t347;
t308 = sin(pkin(12));
t311 = cos(pkin(12));
t208 = -t243 * t308 + t271 * t311;
t354 = t271 * t308;
t209 = t243 * t311 + t354;
t332 = t310 * t333;
t242 = t272 * t302 + t319 * t332;
t370 = -Icges(5,4) * t243 + Icges(6,5) * t209 - Icges(5,6) * t271 + Icges(6,6) * t208 + t371 * t242;
t274 = -t313 * t347 + t344;
t350 = t310 * t317;
t245 = t274 * t333 + t302 * t350;
t273 = t313 * t345 + t346;
t210 = -t245 * t308 + t273 * t311;
t353 = t273 * t308;
t211 = t245 * t311 + t353;
t244 = t274 * t302 - t317 * t332;
t369 = -Icges(5,4) * t245 + Icges(6,5) * t211 - Icges(5,6) * t273 + Icges(6,6) * t210 + t371 * t244;
t261 = t313 * t302 + t316 * t332;
t349 = t310 * t318;
t240 = -t261 * t308 - t311 * t349;
t336 = t308 * t349;
t241 = t261 * t311 - t336;
t351 = t310 * t316;
t260 = t302 * t351 - t313 * t333;
t368 = -Icges(5,4) * t261 + Icges(6,5) * t241 + Icges(5,6) * t349 + Icges(6,6) * t240 + t371 * t260;
t309 = sin(pkin(11));
t312 = cos(pkin(11));
t246 = -t272 * t309 - t312 * t348;
t334 = t309 * t348;
t247 = t272 * t312 - t334;
t191 = Icges(4,5) * t247 + Icges(4,6) * t246 + Icges(4,3) * t271;
t227 = Icges(3,4) * t272 - Icges(3,2) * t271 - Icges(3,6) * t348;
t367 = t191 - t227;
t248 = -t274 * t309 + t312 * t350;
t335 = t309 * t350;
t249 = t274 * t312 + t335;
t192 = Icges(4,5) * t249 + Icges(4,6) * t248 + Icges(4,3) * t273;
t228 = Icges(3,4) * t274 - Icges(3,2) * t273 + Icges(3,6) * t350;
t366 = t192 - t228;
t269 = -t309 * t351 + t312 * t313;
t352 = t309 * t313;
t270 = t312 * t351 + t352;
t219 = Icges(4,5) * t270 + Icges(4,6) * t269 - Icges(4,3) * t349;
t258 = Icges(3,6) * t313 + (Icges(3,4) * t316 + Icges(3,2) * t318) * t310;
t365 = t219 - t258;
t358 = pkin(8) * t313;
t357 = pkin(3) * t312;
t356 = pkin(5) * t311;
t355 = Icges(2,4) * t317;
t341 = qJD(2) * t310;
t339 = V_base(5) * pkin(7) + V_base(1);
t284 = t317 * t341 + V_base(4);
t304 = V_base(6) + qJD(1);
t251 = qJD(4) * t273 + t284;
t285 = qJD(2) * t313 + t304;
t283 = -t319 * t341 + V_base(5);
t277 = t317 * pkin(1) - pkin(8) * t348;
t331 = -t277 * t304 + V_base(5) * t358 + t339;
t278 = pkin(1) * t319 + pkin(8) * t350;
t330 = V_base(4) * t277 - t278 * V_base(5) + V_base(3);
t250 = qJD(4) * t271 + t283;
t267 = -qJD(4) * t349 + t285;
t275 = (pkin(2) * t316 - qJ(3) * t318) * t310;
t329 = qJD(3) * t273 + t283 * t275 + t331;
t328 = t304 * t278 + V_base(2) + (-pkin(7) - t358) * V_base(4);
t237 = pkin(2) * t274 + qJ(3) * t273;
t327 = qJD(3) * t271 + t285 * t237 + t328;
t236 = t272 * pkin(2) + t271 * qJ(3);
t326 = -qJD(3) * t349 + t284 * t236 + t330;
t189 = -pkin(3) * t334 + pkin(9) * t271 + t272 * t357;
t224 = pkin(3) * t352 + (-pkin(9) * t318 + t316 * t357) * t310;
t325 = t283 * t224 + (-t189 - t236) * t285 + t329;
t222 = pkin(4) * t261 + qJ(5) * t260;
t324 = qJD(5) * t244 + t250 * t222 + t325;
t190 = pkin(3) * t335 + pkin(9) * t273 + t274 * t357;
t323 = t285 * t190 + (-t224 - t275) * t284 + t327;
t322 = t284 * t189 + (-t190 - t237) * t283 + t326;
t201 = pkin(4) * t245 + qJ(5) * t244;
t321 = qJD(5) * t242 + t267 * t201 + t323;
t200 = pkin(4) * t243 + qJ(5) * t242;
t320 = qJD(5) * t260 + t251 * t200 + t322;
t307 = pkin(12) + qJ(6);
t305 = Icges(2,4) * t319;
t303 = cos(t307);
t301 = sin(t307);
t293 = rSges(2,1) * t319 - t317 * rSges(2,2);
t292 = t317 * rSges(2,1) + rSges(2,2) * t319;
t291 = Icges(2,1) * t319 - t355;
t290 = Icges(2,1) * t317 + t305;
t289 = -Icges(2,2) * t317 + t305;
t288 = Icges(2,2) * t319 + t355;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t262 = rSges(3,3) * t313 + (rSges(3,1) * t316 + rSges(3,2) * t318) * t310;
t259 = Icges(3,5) * t313 + (Icges(3,1) * t316 + Icges(3,4) * t318) * t310;
t257 = Icges(3,3) * t313 + (Icges(3,5) * t316 + Icges(3,6) * t318) * t310;
t255 = V_base(5) * rSges(2,3) - t292 * t304 + t339;
t254 = t293 * t304 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t252 = t292 * V_base(4) - t293 * V_base(5) + V_base(3);
t235 = t261 * t303 - t301 * t349;
t234 = -t261 * t301 - t303 * t349;
t233 = qJD(6) * t260 + t267;
t232 = rSges(3,1) * t274 - rSges(3,2) * t273 + rSges(3,3) * t350;
t231 = t272 * rSges(3,1) - t271 * rSges(3,2) - rSges(3,3) * t348;
t230 = Icges(3,1) * t274 - Icges(3,4) * t273 + Icges(3,5) * t350;
t229 = Icges(3,1) * t272 - Icges(3,4) * t271 - Icges(3,5) * t348;
t226 = Icges(3,5) * t274 - Icges(3,6) * t273 + Icges(3,3) * t350;
t225 = Icges(3,5) * t272 - Icges(3,6) * t271 - Icges(3,3) * t348;
t223 = rSges(4,1) * t270 + rSges(4,2) * t269 - rSges(4,3) * t349;
t221 = Icges(4,1) * t270 + Icges(4,4) * t269 - Icges(4,5) * t349;
t220 = Icges(4,4) * t270 + Icges(4,2) * t269 - Icges(4,6) * t349;
t216 = rSges(5,1) * t261 - rSges(5,2) * t260 - rSges(5,3) * t349;
t215 = Icges(5,1) * t261 - Icges(5,4) * t260 - Icges(5,5) * t349;
t213 = Icges(5,5) * t261 - Icges(5,6) * t260 - Icges(5,3) * t349;
t207 = t245 * t303 + t273 * t301;
t206 = -t245 * t301 + t273 * t303;
t205 = t243 * t303 + t271 * t301;
t204 = -t243 * t301 + t271 * t303;
t203 = qJD(6) * t244 + t251;
t202 = qJD(6) * t242 + t250;
t198 = rSges(4,1) * t249 + rSges(4,2) * t248 + rSges(4,3) * t273;
t197 = rSges(4,1) * t247 + rSges(4,2) * t246 + rSges(4,3) * t271;
t196 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t273;
t195 = Icges(4,1) * t247 + Icges(4,4) * t246 + Icges(4,5) * t271;
t194 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t273;
t193 = Icges(4,4) * t247 + Icges(4,2) * t246 + Icges(4,6) * t271;
t188 = rSges(5,1) * t245 - rSges(5,2) * t244 + rSges(5,3) * t273;
t187 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t271;
t186 = Icges(5,1) * t245 - Icges(5,4) * t244 + Icges(5,5) * t273;
t185 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t271;
t182 = Icges(5,5) * t245 - Icges(5,6) * t244 + Icges(5,3) * t273;
t181 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t271;
t179 = rSges(6,1) * t241 + rSges(6,2) * t240 + rSges(6,3) * t260;
t178 = Icges(6,1) * t241 + Icges(6,4) * t240 + Icges(6,5) * t260;
t177 = Icges(6,4) * t241 + Icges(6,2) * t240 + Icges(6,6) * t260;
t172 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t260;
t171 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t260;
t170 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t260;
t169 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t260;
t168 = -pkin(5) * t336 + pkin(10) * t260 + t261 * t356;
t167 = -t231 * t285 + t262 * t283 + t331;
t166 = t232 * t285 - t262 * t284 + t328;
t165 = t231 * t284 - t232 * t283 + t330;
t164 = rSges(6,1) * t211 + rSges(6,2) * t210 + rSges(6,3) * t244;
t163 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t242;
t162 = Icges(6,1) * t211 + Icges(6,4) * t210 + Icges(6,5) * t244;
t161 = Icges(6,1) * t209 + Icges(6,4) * t208 + Icges(6,5) * t242;
t160 = Icges(6,4) * t211 + Icges(6,2) * t210 + Icges(6,6) * t244;
t159 = Icges(6,4) * t209 + Icges(6,2) * t208 + Icges(6,6) * t242;
t156 = rSges(7,1) * t207 + rSges(7,2) * t206 + rSges(7,3) * t244;
t155 = rSges(7,1) * t205 + rSges(7,2) * t204 + rSges(7,3) * t242;
t154 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t244;
t153 = Icges(7,1) * t205 + Icges(7,4) * t204 + Icges(7,5) * t242;
t152 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t244;
t151 = Icges(7,4) * t205 + Icges(7,2) * t204 + Icges(7,6) * t242;
t150 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t244;
t149 = Icges(7,5) * t205 + Icges(7,6) * t204 + Icges(7,3) * t242;
t148 = pkin(5) * t353 + pkin(10) * t244 + t245 * t356;
t147 = pkin(5) * t354 + pkin(10) * t242 + t243 * t356;
t146 = t223 * t283 + (-t197 - t236) * t285 + t329;
t145 = t198 * t285 + (-t223 - t275) * t284 + t327;
t144 = t197 * t284 + (-t198 - t237) * t283 + t326;
t143 = -t187 * t267 + t216 * t250 + t325;
t142 = t188 * t267 - t216 * t251 + t323;
t141 = t187 * t251 - t188 * t250 + t322;
t140 = t179 * t250 + (-t163 - t200) * t267 + t324;
t139 = t164 * t267 + (-t179 - t222) * t251 + t321;
t138 = t163 * t251 + (-t164 - t201) * t250 + t320;
t137 = t324 - t155 * t233 + t168 * t250 + t172 * t202 + (-t147 - t200) * t267;
t136 = t148 * t267 + t156 * t233 - t172 * t203 + (-t168 - t222) * t251 + t321;
t135 = t147 * t251 + t155 * t203 - t156 * t202 + (-t148 - t201) * t250 + t320;
t1 = m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t233 * ((t150 * t260 + t152 * t234 + t154 * t235) * t203 + (t149 * t260 + t151 * t234 + t153 * t235) * t202 + (t169 * t260 + t170 * t234 + t171 * t235) * t233) / 0.2e1 + m(2) * (t252 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t202 * ((t150 * t242 + t152 * t204 + t154 * t205) * t203 + (t149 * t242 + t151 * t204 + t153 * t205) * t202 + (t169 * t242 + t170 * t204 + t171 * t205) * t233) / 0.2e1 + t203 * ((t150 * t244 + t152 * t206 + t154 * t207) * t203 + (t149 * t244 + t151 * t206 + t153 * t207) * t202 + (t169 * t244 + t170 * t206 + t171 * t207) * t233) / 0.2e1 + m(3) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + ((t177 * t208 + t178 * t209 + t213 * t271 + t215 * t243 + t242 * t368) * t267 + (t160 * t208 + t162 * t209 + t182 * t271 + t186 * t243 + t242 * t369) * t251 + (t159 * t208 + t161 * t209 + t181 * t271 + t185 * t243 + t242 * t370) * t250) * t250 / 0.2e1 + ((t177 * t210 + t178 * t211 + t213 * t273 + t215 * t245 + t244 * t368) * t267 + (t160 * t210 + t162 * t211 + t182 * t273 + t186 * t245 + t244 * t369) * t251 + (t159 * t210 + t161 * t211 + t181 * t273 + t185 * t245 + t244 * t370) * t250) * t251 / 0.2e1 + ((t177 * t240 + t178 * t241 - t213 * t349 + t215 * t261 + t260 * t368) * t267 + (t160 * t240 + t162 * t241 - t182 * t349 + t186 * t261 + t260 * t369) * t251 + (t159 * t240 + t161 * t241 - t181 * t349 + t185 * t261 + t260 * t370) * t250) * t267 / 0.2e1 + ((t220 * t246 + t221 * t247 - t257 * t348 + t272 * t259 + t271 * t365) * t285 + (t194 * t246 + t196 * t247 - t226 * t348 + t272 * t230 + t271 * t366) * t284 + (t193 * t246 + t195 * t247 - t225 * t348 + t272 * t229 + t271 * t367) * t283) * t283 / 0.2e1 + ((t220 * t248 + t221 * t249 + t257 * t350 + t259 * t274 + t273 * t365) * t285 + (t194 * t248 + t196 * t249 + t226 * t350 + t230 * t274 + t273 * t366) * t284 + (t193 * t248 + t195 * t249 + t225 * t350 + t229 * t274 + t273 * t367) * t283) * t284 / 0.2e1 + ((t225 * t283 + t226 * t284 + t257 * t285) * t313 + ((t228 * t318 + t230 * t316) * t284 + (t227 * t318 + t229 * t316) * t283 + (t258 * t318 + t259 * t316) * t285) * t310 + (-t192 * t349 + t194 * t269 + t196 * t270) * t284 + (-t191 * t349 + t193 * t269 + t195 * t270) * t283 + (-t219 * t349 + t220 * t269 + t221 * t270) * t285) * t285 / 0.2e1 + ((-t317 * t288 + t290 * t319 + Icges(1,4)) * V_base(5) + (-t317 * t289 + t291 * t319 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t288 * t319 + t317 * t290 + Icges(1,2)) * V_base(5) + (t289 * t319 + t317 * t291 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t317 + Icges(2,6) * t319) * V_base(5) + (Icges(2,5) * t319 - Icges(2,6) * t317) * V_base(4) + Icges(2,3) * t304 / 0.2e1) * t304;
T  = t1;
