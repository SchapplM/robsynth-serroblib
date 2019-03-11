% Calculate kinetic energy for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:21
% EndTime: 2019-03-08 19:34:27
% DurationCPUTime: 6.10s
% Computational Cost: add. (3245->395), mult. (7754->555), div. (0->0), fcn. (9624->12), ass. (0->178)
t382 = Icges(5,1) + Icges(6,2);
t381 = -Icges(6,1) - Icges(5,3);
t380 = -Icges(5,4) - Icges(6,6);
t379 = Icges(6,4) - Icges(5,5);
t378 = Icges(6,5) - Icges(5,6);
t377 = Icges(5,2) + Icges(6,3);
t314 = sin(qJ(2));
t350 = sin(pkin(11));
t352 = cos(pkin(6));
t330 = t352 * t350;
t351 = cos(pkin(11));
t331 = t352 * t351;
t356 = cos(qJ(2));
t271 = t314 * t331 + t330 * t356;
t280 = -t314 * t350 + t356 * t351;
t310 = sin(pkin(10));
t312 = cos(pkin(10));
t252 = t271 * t312 + t280 * t310;
t311 = sin(pkin(6));
t355 = cos(qJ(4));
t339 = t311 * t355;
t354 = sin(qJ(4));
t228 = t252 * t354 + t312 * t339;
t338 = t311 * t354;
t229 = t252 * t355 - t312 * t338;
t323 = -t314 * t330 + t331 * t356;
t327 = t314 * t351 + t350 * t356;
t251 = -t310 * t327 + t312 * t323;
t374 = t377 * t228 + t380 * t229 - t378 * t251;
t254 = -t271 * t310 + t280 * t312;
t230 = t254 * t354 - t310 * t339;
t231 = t254 * t355 + t310 * t338;
t253 = -t310 * t323 - t312 * t327;
t373 = t377 * t230 + t380 * t231 - t378 * t253;
t372 = t378 * t228 - t379 * t229 + t381 * t251;
t371 = t378 * t230 - t379 * t231 + t381 * t253;
t370 = t380 * t228 + t382 * t229 + t379 * t251;
t369 = t380 * t230 + t382 * t231 + t379 * t253;
t347 = t311 * t312;
t200 = Icges(4,5) * t252 + Icges(4,6) * t251 - Icges(4,3) * t347;
t333 = t352 * t356;
t273 = -t310 * t314 + t312 * t333;
t335 = t314 * t352;
t274 = t310 * t356 + t312 * t335;
t238 = Icges(3,5) * t274 + Icges(3,6) * t273 - Icges(3,3) * t347;
t368 = t200 + t238;
t348 = t310 * t311;
t201 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t348;
t275 = -t310 * t333 - t312 * t314;
t276 = -t310 * t335 + t312 * t356;
t239 = Icges(3,5) * t276 + Icges(3,6) * t275 + Icges(3,3) * t348;
t367 = t201 + t239;
t270 = t327 * t311;
t258 = t270 * t354 - t352 * t355;
t259 = t270 * t355 + t352 * t354;
t269 = t280 * t311;
t366 = t377 * t258 + t380 * t259 - t378 * t269;
t365 = t378 * t258 - t379 * t259 + t381 * t269;
t364 = t380 * t258 + t382 * t259 + t379 * t269;
t234 = Icges(4,5) * t270 + Icges(4,6) * t269 + Icges(4,3) * t352;
t265 = Icges(3,3) * t352 + (Icges(3,5) * t314 + Icges(3,6) * t356) * t311;
t363 = t234 + t265;
t353 = pkin(2) * t356;
t349 = Icges(2,4) * t310;
t346 = qJD(2) * t311;
t345 = qJD(3) * t311;
t344 = V_base(5) * qJ(1) + V_base(1);
t340 = qJD(1) + V_base(3);
t337 = t352 * pkin(7);
t288 = t310 * t346 + V_base(4);
t300 = qJD(2) * t352 + V_base(6);
t336 = pkin(2) * t335 - qJ(3) * t311;
t227 = -qJD(4) * t253 + t288;
t257 = -qJD(4) * t269 + t300;
t287 = -t312 * t346 + V_base(5);
t226 = -qJD(4) * t251 + t287;
t282 = pkin(1) * t310 - pkin(7) * t347;
t329 = -t282 * V_base(6) + V_base(5) * t337 + t344;
t283 = pkin(1) * t312 + pkin(7) * t348;
t328 = V_base(4) * t282 - t283 * V_base(5) + t340;
t281 = t311 * t314 * pkin(2) + qJ(3) * t352;
t326 = t287 * t281 + t310 * t345 + t329;
t246 = t310 * t353 + t312 * t336;
t325 = qJD(3) * t352 + t288 * t246 + t328;
t324 = V_base(6) * t283 + V_base(2) + (-t337 - qJ(1)) * V_base(4);
t208 = pkin(3) * t252 - pkin(8) * t251;
t250 = pkin(3) * t270 - pkin(8) * t269;
t322 = t287 * t250 + (-t208 - t246) * t300 + t326;
t209 = pkin(3) * t254 - pkin(8) * t253;
t247 = -t310 * t336 + t312 * t353;
t321 = t288 * t208 + (-t209 - t247) * t287 + t325;
t320 = t300 * t247 - t312 * t345 + t324;
t218 = pkin(4) * t259 + qJ(5) * t258;
t319 = qJD(5) * t230 + t226 * t218 + t322;
t188 = pkin(4) * t229 + qJ(5) * t228;
t318 = qJD(5) * t258 + t227 * t188 + t321;
t317 = t300 * t209 + (-t250 - t281) * t288 + t320;
t189 = pkin(4) * t231 + qJ(5) * t230;
t316 = qJD(5) * t228 + t257 * t189 + t317;
t315 = cos(qJ(6));
t313 = sin(qJ(6));
t308 = Icges(2,4) * t312;
t296 = rSges(2,1) * t312 - rSges(2,2) * t310;
t295 = rSges(2,1) * t310 + rSges(2,2) * t312;
t294 = Icges(2,1) * t312 - t349;
t293 = Icges(2,1) * t310 + t308;
t292 = -Icges(2,2) * t310 + t308;
t291 = Icges(2,2) * t312 + t349;
t286 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t285 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t284 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t268 = t352 * rSges(3,3) + (rSges(3,1) * t314 + rSges(3,2) * t356) * t311;
t267 = Icges(3,5) * t352 + (Icges(3,1) * t314 + Icges(3,4) * t356) * t311;
t266 = Icges(3,6) * t352 + (Icges(3,4) * t314 + Icges(3,2) * t356) * t311;
t261 = V_base(5) * rSges(2,3) - t295 * V_base(6) + t344;
t260 = t296 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t256 = t295 * V_base(4) - t296 * V_base(5) + t340;
t245 = rSges(3,1) * t276 + rSges(3,2) * t275 + rSges(3,3) * t348;
t244 = rSges(3,1) * t274 + rSges(3,2) * t273 - rSges(3,3) * t347;
t243 = Icges(3,1) * t276 + Icges(3,4) * t275 + Icges(3,5) * t348;
t242 = Icges(3,1) * t274 + Icges(3,4) * t273 - Icges(3,5) * t347;
t241 = Icges(3,4) * t276 + Icges(3,2) * t275 + Icges(3,6) * t348;
t240 = Icges(3,4) * t274 + Icges(3,2) * t273 - Icges(3,6) * t347;
t237 = t270 * rSges(4,1) + t269 * rSges(4,2) + rSges(4,3) * t352;
t236 = Icges(4,1) * t270 + Icges(4,4) * t269 + Icges(4,5) * t352;
t235 = Icges(4,4) * t270 + Icges(4,2) * t269 + Icges(4,6) * t352;
t232 = -pkin(5) * t269 + pkin(9) * t259;
t221 = t258 * t313 - t269 * t315;
t220 = t258 * t315 + t269 * t313;
t219 = qJD(6) * t259 + t257;
t217 = rSges(5,1) * t259 - rSges(5,2) * t258 - rSges(5,3) * t269;
t216 = -rSges(6,1) * t269 - rSges(6,2) * t259 + rSges(6,3) * t258;
t207 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t348;
t206 = rSges(4,1) * t252 + rSges(4,2) * t251 - rSges(4,3) * t347;
t205 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t348;
t204 = Icges(4,1) * t252 + Icges(4,4) * t251 - Icges(4,5) * t347;
t203 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t348;
t202 = Icges(4,4) * t252 + Icges(4,2) * t251 - Icges(4,6) * t347;
t197 = -pkin(5) * t253 + pkin(9) * t231;
t196 = -pkin(5) * t251 + pkin(9) * t229;
t195 = t230 * t313 - t253 * t315;
t194 = t230 * t315 + t253 * t313;
t193 = t228 * t313 - t251 * t315;
t192 = t228 * t315 + t251 * t313;
t191 = qJD(6) * t231 + t227;
t190 = qJD(6) * t229 + t226;
t187 = -t244 * t300 + t268 * t287 + t329;
t186 = t300 * t245 - t288 * t268 + t324;
t183 = rSges(7,1) * t221 + rSges(7,2) * t220 + rSges(7,3) * t259;
t182 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t259;
t181 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t259;
t180 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t259;
t179 = t244 * t288 - t245 * t287 + t328;
t178 = rSges(5,1) * t231 - rSges(5,2) * t230 - rSges(5,3) * t253;
t177 = rSges(5,1) * t229 - rSges(5,2) * t228 - rSges(5,3) * t251;
t176 = -rSges(6,1) * t253 - rSges(6,2) * t231 + rSges(6,3) * t230;
t175 = -rSges(6,1) * t251 - rSges(6,2) * t229 + rSges(6,3) * t228;
t161 = rSges(7,1) * t195 + rSges(7,2) * t194 + rSges(7,3) * t231;
t160 = rSges(7,1) * t193 + rSges(7,2) * t192 + rSges(7,3) * t229;
t159 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t231;
t158 = Icges(7,1) * t193 + Icges(7,4) * t192 + Icges(7,5) * t229;
t157 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t231;
t156 = Icges(7,4) * t193 + Icges(7,2) * t192 + Icges(7,6) * t229;
t155 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t231;
t154 = Icges(7,5) * t193 + Icges(7,6) * t192 + Icges(7,3) * t229;
t153 = t237 * t287 + (-t206 - t246) * t300 + t326;
t152 = t300 * t207 + (-t237 - t281) * t288 + t320;
t151 = t206 * t288 + (-t207 - t247) * t287 + t325;
t150 = -t177 * t257 + t217 * t226 + t322;
t149 = t257 * t178 - t227 * t217 + t317;
t148 = t177 * t227 - t178 * t226 + t321;
t147 = t216 * t226 + (-t175 - t188) * t257 + t319;
t146 = t257 * t176 + (-t216 - t218) * t227 + t316;
t145 = t175 * t227 + (-t176 - t189) * t226 + t318;
t144 = t319 - t160 * t219 + t183 * t190 + t226 * t232 + (-t188 - t196) * t257;
t143 = t219 * t161 - t191 * t183 + t257 * t197 + (-t218 - t232) * t227 + t316;
t142 = t318 + t160 * t191 - t161 * t190 + t196 * t227 + (-t189 - t197) * t226;
t1 = m(6) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(3) * (t179 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(5) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(4) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(7) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + t190 * ((t155 * t229 + t157 * t192 + t159 * t193) * t191 + (t229 * t154 + t156 * t192 + t193 * t158) * t190 + (t180 * t229 + t181 * t192 + t182 * t193) * t219) / 0.2e1 + t191 * ((t231 * t155 + t157 * t194 + t159 * t195) * t191 + (t154 * t231 + t156 * t194 + t158 * t195) * t190 + (t180 * t231 + t181 * t194 + t182 * t195) * t219) / 0.2e1 + t219 * ((t155 * t259 + t157 * t220 + t159 * t221) * t191 + (t154 * t259 + t156 * t220 + t158 * t221) * t190 + (t259 * t180 + t220 * t181 + t221 * t182) * t219) / 0.2e1 + m(2) * (t256 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + m(1) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + ((t228 * t366 + t229 * t364 - t251 * t365) * t257 + (t228 * t373 + t229 * t369 - t251 * t371) * t227 + (t374 * t228 + t370 * t229 - t372 * t251) * t226) * t226 / 0.2e1 + ((t230 * t366 + t231 * t364 - t253 * t365) * t257 + (t373 * t230 + t369 * t231 - t371 * t253) * t227 + (t230 * t374 + t231 * t370 - t253 * t372) * t226) * t227 / 0.2e1 + ((t366 * t258 + t364 * t259 - t365 * t269) * t257 + (t258 * t373 + t259 * t369 - t269 * t371) * t227 + (t258 * t374 + t259 * t370 - t269 * t372) * t226) * t257 / 0.2e1 + ((t235 * t251 + t236 * t252 + t266 * t273 + t267 * t274 - t347 * t363) * t300 + (t203 * t251 + t205 * t252 + t241 * t273 + t243 * t274 - t347 * t367) * t288 + (t251 * t202 + t252 * t204 + t273 * t240 + t274 * t242 - t368 * t347) * t287) * t287 / 0.2e1 + ((t235 * t253 + t236 * t254 + t266 * t275 + t267 * t276 + t348 * t363) * t300 + (t253 * t203 + t254 * t205 + t275 * t241 + t276 * t243 + t367 * t348) * t288 + (t202 * t253 + t204 * t254 + t240 * t275 + t242 * t276 + t348 * t368) * t287) * t288 / 0.2e1 + (((t241 * t356 + t243 * t314) * t288 + (t240 * t356 + t242 * t314) * t287 + (t266 * t356 + t267 * t314) * t300) * t311 + (t238 * t287 + t239 * t288 + t265 * t300) * t352 + (t201 * t352 + t269 * t203 + t270 * t205) * t288 + (t200 * t352 + t269 * t202 + t270 * t204) * t287 + (t234 * t352 + t269 * t235 + t270 * t236) * t300) * t300 / 0.2e1 + ((-t291 * t310 + t293 * t312 + Icges(1,4)) * V_base(5) + (-t310 * t292 + t312 * t294 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t312 * t291 + t310 * t293 + Icges(1,2)) * V_base(5) + (t292 * t312 + t294 * t310 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t310 + Icges(2,6) * t312 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t312 - Icges(2,6) * t310 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
