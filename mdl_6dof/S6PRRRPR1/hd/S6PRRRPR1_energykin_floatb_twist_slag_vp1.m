% Calculate kinetic energy for
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:05
% EndTime: 2019-03-08 23:00:10
% DurationCPUTime: 4.82s
% Computational Cost: add. (3677->438), mult. (5584->632), div. (0->0), fcn. (6496->14), ass. (0->189)
t375 = Icges(5,3) + Icges(6,3);
t316 = sin(pkin(11));
t318 = cos(pkin(11));
t325 = cos(qJ(2));
t319 = cos(pkin(6));
t322 = sin(qJ(2));
t354 = t319 * t322;
t274 = t316 * t325 + t318 * t354;
t315 = qJ(3) + qJ(4);
t342 = pkin(12) + t315;
t306 = sin(t342);
t317 = sin(pkin(6));
t340 = cos(t342);
t339 = t317 * t340;
t236 = t274 * t306 + t318 * t339;
t360 = t317 * t318;
t237 = t274 * t340 - t306 * t360;
t310 = sin(t315);
t311 = cos(t315);
t244 = -t274 * t310 - t311 * t360;
t245 = t274 * t311 - t310 * t360;
t353 = t319 * t325;
t273 = t316 * t322 - t318 * t353;
t373 = Icges(5,5) * t245 + Icges(6,5) * t237 + Icges(5,6) * t244 - Icges(6,6) * t236 + t273 * t375;
t276 = -t316 * t354 + t318 * t325;
t238 = t276 * t306 - t316 * t339;
t361 = t316 * t317;
t239 = t276 * t340 + t306 * t361;
t246 = -t276 * t310 + t311 * t361;
t247 = t276 * t311 + t310 * t361;
t275 = t316 * t353 + t318 * t322;
t372 = Icges(5,5) * t247 + Icges(6,5) * t239 + Icges(5,6) * t246 - Icges(6,6) * t238 + t375 * t275;
t358 = t317 * t322;
t259 = t306 * t358 - t319 * t340;
t260 = t306 * t319 + t322 * t339;
t265 = -t310 * t358 + t311 * t319;
t266 = t310 * t319 + t311 * t358;
t356 = t317 * t325;
t371 = Icges(5,5) * t266 + Icges(6,5) * t260 + Icges(5,6) * t265 - Icges(6,6) * t259 - t356 * t375;
t365 = pkin(7) * t319;
t324 = cos(qJ(3));
t364 = t324 * pkin(3);
t362 = Icges(2,4) * t316;
t321 = sin(qJ(3));
t359 = t317 * t321;
t357 = t317 * t324;
t355 = t319 * t321;
t352 = pkin(4) * t311;
t350 = qJD(2) * t317;
t349 = V_base(5) * qJ(1) + V_base(1);
t345 = qJD(1) + V_base(3);
t344 = t316 * t359;
t343 = t318 * t359;
t291 = t316 * t350 + V_base(4);
t303 = qJD(2) * t319 + V_base(6);
t341 = pkin(4) * t310;
t250 = qJD(3) * t275 + t291;
t221 = qJD(4) * t275 + t250;
t290 = -t318 * t350 + V_base(5);
t249 = qJD(3) * t273 + t290;
t283 = pkin(1) * t316 - pkin(7) * t360;
t338 = -t283 * V_base(6) + t365 * V_base(5) + t349;
t284 = pkin(1) * t318 + pkin(7) * t361;
t337 = t283 * V_base(4) - t284 * V_base(5) + t345;
t220 = qJD(4) * t273 + t249;
t258 = (-qJD(3) - qJD(4)) * t356 + t303;
t336 = V_base(6) * t284 + V_base(2) + (-qJ(1) - t365) * V_base(4);
t240 = pkin(2) * t274 + pkin(8) * t273;
t282 = (pkin(2) * t322 - pkin(8) * t325) * t317;
t335 = -t240 * t303 + t282 * t290 + t338;
t241 = pkin(2) * t276 + pkin(8) * t275;
t334 = t240 * t291 - t241 * t290 + t337;
t333 = t241 * t303 - t282 * t291 + t336;
t189 = -pkin(3) * t343 + pkin(9) * t273 + t274 * t364;
t235 = pkin(3) * t355 + (-pkin(9) * t325 + t322 * t364) * t317;
t277 = -qJD(3) * t356 + t303;
t332 = -t189 * t277 + t235 * t249 + t335;
t190 = pkin(3) * t344 + pkin(9) * t275 + t276 * t364;
t331 = t189 * t250 - t190 * t249 + t334;
t330 = t190 * t277 - t235 * t250 + t333;
t204 = t341 * t319 + (-qJ(5) * t325 + t322 * t352) * t317;
t329 = qJD(5) * t275 + t204 * t220 + t332;
t164 = qJ(5) * t275 + t276 * t352 + t341 * t361;
t328 = qJD(5) * t273 + t164 * t258 + t330;
t163 = qJ(5) * t273 + t274 * t352 - t341 * t360;
t327 = -qJD(5) * t356 + t163 * t221 + t331;
t323 = cos(qJ(6));
t320 = sin(qJ(6));
t309 = Icges(2,4) * t318;
t300 = rSges(2,1) * t318 - rSges(2,2) * t316;
t299 = rSges(2,1) * t316 + rSges(2,2) * t318;
t298 = Icges(2,1) * t318 - t362;
t297 = Icges(2,1) * t316 + t309;
t296 = -Icges(2,2) * t316 + t309;
t295 = Icges(2,2) * t318 + t362;
t289 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t288 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t287 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t281 = t322 * t357 + t355;
t280 = t319 * t324 - t321 * t358;
t264 = rSges(3,3) * t319 + (rSges(3,1) * t322 + rSges(3,2) * t325) * t317;
t263 = Icges(3,5) * t319 + (Icges(3,1) * t322 + Icges(3,4) * t325) * t317;
t262 = Icges(3,6) * t319 + (Icges(3,4) * t322 + Icges(3,2) * t325) * t317;
t261 = Icges(3,3) * t319 + (Icges(3,5) * t322 + Icges(3,6) * t325) * t317;
t257 = V_base(5) * rSges(2,3) - t299 * V_base(6) + t349;
t256 = t300 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t254 = t276 * t324 + t344;
t253 = -t276 * t321 + t316 * t357;
t252 = t274 * t324 - t343;
t251 = -t274 * t321 - t318 * t357;
t248 = t299 * V_base(4) - t300 * V_base(5) + t345;
t243 = t260 * t323 - t320 * t356;
t242 = -t260 * t320 - t323 * t356;
t234 = rSges(4,1) * t281 + rSges(4,2) * t280 - rSges(4,3) * t356;
t233 = Icges(4,1) * t281 + Icges(4,4) * t280 - Icges(4,5) * t356;
t232 = Icges(4,4) * t281 + Icges(4,2) * t280 - Icges(4,6) * t356;
t231 = Icges(4,5) * t281 + Icges(4,6) * t280 - Icges(4,3) * t356;
t230 = rSges(3,1) * t276 - rSges(3,2) * t275 + rSges(3,3) * t361;
t229 = rSges(3,1) * t274 - rSges(3,2) * t273 - rSges(3,3) * t360;
t228 = Icges(3,1) * t276 - Icges(3,4) * t275 + Icges(3,5) * t361;
t227 = Icges(3,1) * t274 - Icges(3,4) * t273 - Icges(3,5) * t360;
t226 = Icges(3,4) * t276 - Icges(3,2) * t275 + Icges(3,6) * t361;
t225 = Icges(3,4) * t274 - Icges(3,2) * t273 - Icges(3,6) * t360;
t224 = Icges(3,5) * t276 - Icges(3,6) * t275 + Icges(3,3) * t361;
t223 = Icges(3,5) * t274 - Icges(3,6) * t273 - Icges(3,3) * t360;
t218 = qJD(6) * t259 + t258;
t217 = pkin(5) * t260 + pkin(10) * t259;
t216 = rSges(5,1) * t266 + rSges(5,2) * t265 - rSges(5,3) * t356;
t215 = Icges(5,1) * t266 + Icges(5,4) * t265 - Icges(5,5) * t356;
t214 = Icges(5,4) * t266 + Icges(5,2) * t265 - Icges(5,6) * t356;
t212 = rSges(6,1) * t260 - rSges(6,2) * t259 - rSges(6,3) * t356;
t211 = Icges(6,1) * t260 - Icges(6,4) * t259 - Icges(6,5) * t356;
t210 = Icges(6,4) * t260 - Icges(6,2) * t259 - Icges(6,6) * t356;
t208 = t239 * t323 + t275 * t320;
t207 = -t239 * t320 + t275 * t323;
t206 = t237 * t323 + t273 * t320;
t205 = -t237 * t320 + t273 * t323;
t202 = pkin(5) * t239 + pkin(10) * t238;
t201 = pkin(5) * t237 + pkin(10) * t236;
t200 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t275;
t199 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t273;
t198 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t275;
t197 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t273;
t196 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t275;
t195 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t273;
t194 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t275;
t193 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t273;
t192 = qJD(6) * t238 + t221;
t191 = qJD(6) * t236 + t220;
t188 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t275;
t187 = rSges(5,1) * t245 + rSges(5,2) * t244 + rSges(5,3) * t273;
t186 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t275;
t185 = Icges(5,1) * t245 + Icges(5,4) * t244 + Icges(5,5) * t273;
t184 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t275;
t183 = Icges(5,4) * t245 + Icges(5,2) * t244 + Icges(5,6) * t273;
t180 = rSges(6,1) * t239 - rSges(6,2) * t238 + rSges(6,3) * t275;
t179 = rSges(6,1) * t237 - rSges(6,2) * t236 + rSges(6,3) * t273;
t178 = Icges(6,1) * t239 - Icges(6,4) * t238 + Icges(6,5) * t275;
t177 = Icges(6,1) * t237 - Icges(6,4) * t236 + Icges(6,5) * t273;
t176 = Icges(6,4) * t239 - Icges(6,2) * t238 + Icges(6,6) * t275;
t175 = Icges(6,4) * t237 - Icges(6,2) * t236 + Icges(6,6) * t273;
t172 = rSges(7,1) * t243 + rSges(7,2) * t242 + rSges(7,3) * t259;
t171 = Icges(7,1) * t243 + Icges(7,4) * t242 + Icges(7,5) * t259;
t170 = Icges(7,4) * t243 + Icges(7,2) * t242 + Icges(7,6) * t259;
t169 = Icges(7,5) * t243 + Icges(7,6) * t242 + Icges(7,3) * t259;
t167 = -t229 * t303 + t264 * t290 + t338;
t166 = t230 * t303 - t264 * t291 + t336;
t160 = t229 * t291 - t230 * t290 + t337;
t159 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t238;
t158 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t236;
t157 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t238;
t156 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t236;
t155 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t238;
t154 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t236;
t153 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t238;
t152 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t236;
t150 = -t199 * t277 + t234 * t249 + t335;
t149 = t200 * t277 - t234 * t250 + t333;
t148 = t199 * t250 - t200 * t249 + t334;
t147 = -t187 * t258 + t216 * t220 + t332;
t146 = t188 * t258 - t216 * t221 + t330;
t145 = t187 * t221 - t188 * t220 + t331;
t144 = t212 * t220 + (-t163 - t179) * t258 + t329;
t143 = t180 * t258 + (-t204 - t212) * t221 + t328;
t142 = t179 * t221 + (-t164 - t180) * t220 + t327;
t141 = (-t163 - t201) * t258 - t158 * t218 + t172 * t191 + t217 * t220 + t329;
t140 = t159 * t218 - t172 * t192 + t202 * t258 + (-t204 - t217) * t221 + t328;
t139 = t327 + t158 * t192 - t159 * t191 + t201 * t221 + (-t164 - t202) * t220;
t1 = t303 * ((t223 * t290 + t224 * t291 + t261 * t303) * t319 + ((t226 * t325 + t228 * t322) * t291 + (t225 * t325 + t227 * t322) * t290 + (t262 * t325 + t263 * t322) * t303) * t317) / 0.2e1 + m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(3) * (t160 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t277 * ((-t194 * t356 + t196 * t280 + t198 * t281) * t250 + (-t193 * t356 + t195 * t280 + t197 * t281) * t249 + (-t231 * t356 + t280 * t232 + t281 * t233) * t277) / 0.2e1 + t290 * ((-t224 * t360 - t226 * t273 + t228 * t274) * t291 + (-t223 * t360 - t273 * t225 + t274 * t227) * t290 + (-t261 * t360 - t262 * t273 + t263 * t274) * t303) / 0.2e1 + t291 * ((t224 * t361 - t275 * t226 + t276 * t228) * t291 + (t223 * t361 - t225 * t275 + t227 * t276) * t290 + (t261 * t361 - t262 * t275 + t263 * t276) * t303) / 0.2e1 + t191 * ((t153 * t236 + t155 * t205 + t157 * t206) * t192 + (t236 * t152 + t205 * t154 + t206 * t156) * t191 + (t169 * t236 + t170 * t205 + t171 * t206) * t218) / 0.2e1 + t192 * ((t238 * t153 + t207 * t155 + t208 * t157) * t192 + (t152 * t238 + t154 * t207 + t156 * t208) * t191 + (t169 * t238 + t170 * t207 + t171 * t208) * t218) / 0.2e1 + m(2) * (t248 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t218 * ((t153 * t259 + t155 * t242 + t157 * t243) * t192 + (t152 * t259 + t154 * t242 + t156 * t243) * t191 + (t259 * t169 + t242 * t170 + t243 * t171) * t218) / 0.2e1 + t249 * ((t194 * t273 + t196 * t251 + t198 * t252) * t250 + (t273 * t193 + t251 * t195 + t252 * t197) * t249 + (t231 * t273 + t232 * t251 + t233 * t252) * t277) / 0.2e1 + t250 * ((t275 * t194 + t253 * t196 + t254 * t198) * t250 + (t193 * t275 + t195 * t253 + t197 * t254) * t249 + (t231 * t275 + t232 * t253 + t233 * t254) * t277) / 0.2e1 + m(1) * (t287 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + ((-t210 * t236 + t211 * t237 + t214 * t244 + t215 * t245 + t273 * t371) * t258 + (-t176 * t236 + t178 * t237 + t184 * t244 + t186 * t245 + t273 * t372) * t221 + (-t236 * t175 + t237 * t177 + t244 * t183 + t245 * t185 + t373 * t273) * t220) * t220 / 0.2e1 + ((-t210 * t238 + t211 * t239 + t214 * t246 + t215 * t247 + t275 * t371) * t258 + (-t238 * t176 + t239 * t178 + t246 * t184 + t247 * t186 + t372 * t275) * t221 + (-t175 * t238 + t177 * t239 + t183 * t246 + t185 * t247 + t275 * t373) * t220) * t221 / 0.2e1 + ((-t259 * t210 + t260 * t211 + t265 * t214 + t266 * t215 - t356 * t371) * t258 + (-t176 * t259 + t178 * t260 + t184 * t265 + t186 * t266 - t356 * t372) * t221 + (-t175 * t259 + t177 * t260 + t183 * t265 + t185 * t266 - t356 * t373) * t220) * t258 / 0.2e1 + ((-t295 * t316 + t297 * t318 + Icges(1,4)) * V_base(5) + (-t316 * t296 + t318 * t298 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t318 * t295 + t316 * t297 + Icges(1,2)) * V_base(5) + (t296 * t318 + t298 * t316 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t316 + Icges(2,6) * t318 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t318 - Icges(2,6) * t316 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
