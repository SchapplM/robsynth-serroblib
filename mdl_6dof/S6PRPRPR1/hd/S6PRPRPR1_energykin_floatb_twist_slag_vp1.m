% Calculate kinetic energy for
% S6PRPRPR1
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:34
% EndTime: 2019-03-08 19:24:40
% DurationCPUTime: 6.69s
% Computational Cost: add. (3817->436), mult. (7918->614), div. (0->0), fcn. (9820->14), ass. (0->193)
t393 = -Icges(5,3) - Icges(6,3);
t330 = sin(qJ(2));
t372 = sin(pkin(11));
t374 = cos(pkin(6));
t347 = t374 * t372;
t373 = cos(pkin(11));
t348 = t374 * t373;
t378 = cos(qJ(2));
t284 = t330 * t348 + t347 * t378;
t293 = -t330 * t372 + t378 * t373;
t324 = sin(pkin(10));
t326 = cos(pkin(10));
t263 = t284 * t326 + t293 * t324;
t364 = qJ(4) + pkin(12);
t321 = sin(t364);
t325 = sin(pkin(6));
t351 = cos(t364);
t346 = t325 * t351;
t230 = t263 * t321 + t326 * t346;
t369 = t325 * t326;
t231 = t263 * t351 - t321 * t369;
t329 = sin(qJ(4));
t332 = cos(qJ(4));
t367 = t325 * t332;
t238 = -t263 * t329 - t326 * t367;
t368 = t325 * t329;
t357 = t326 * t368;
t239 = t263 * t332 - t357;
t292 = -t330 * t373 - t372 * t378;
t340 = -t330 * t347 + t348 * t378;
t262 = t324 * t292 + t326 * t340;
t392 = Icges(5,5) * t239 + Icges(6,5) * t231 + Icges(5,6) * t238 - Icges(6,6) * t230 + t393 * t262;
t265 = -t284 * t324 + t293 * t326;
t232 = t265 * t321 - t324 * t346;
t370 = t324 * t325;
t233 = t265 * t351 + t321 * t370;
t240 = -t265 * t329 + t324 * t367;
t358 = t324 * t368;
t241 = t265 * t332 + t358;
t264 = t292 * t326 - t324 * t340;
t391 = Icges(5,5) * t241 + Icges(6,5) * t233 + Icges(5,6) * t240 - Icges(6,6) * t232 + t393 * t264;
t283 = t292 * t325;
t266 = -t283 * t321 - t351 * t374;
t267 = -t283 * t351 + t321 * t374;
t270 = t283 * t329 + t332 * t374;
t352 = t374 * t329;
t271 = -t283 * t332 + t352;
t282 = t293 * t325;
t390 = Icges(5,5) * t271 + Icges(6,5) * t267 + Icges(5,6) * t270 - Icges(6,6) * t266 + t393 * t282;
t211 = Icges(4,5) * t263 + Icges(4,6) * t262 - Icges(4,3) * t369;
t350 = t374 * t378;
t286 = -t324 * t330 + t326 * t350;
t354 = t330 * t374;
t287 = t324 * t378 + t326 * t354;
t247 = Icges(3,5) * t287 + Icges(3,6) * t286 - Icges(3,3) * t369;
t389 = t211 + t247;
t212 = Icges(4,5) * t265 + Icges(4,6) * t264 + Icges(4,3) * t370;
t288 = -t324 * t350 - t326 * t330;
t289 = -t324 * t354 + t326 * t378;
t248 = Icges(3,5) * t289 + Icges(3,6) * t288 + Icges(3,3) * t370;
t388 = t212 + t248;
t243 = -Icges(4,5) * t283 + Icges(4,6) * t282 + Icges(4,3) * t374;
t278 = Icges(3,3) * t374 + (Icges(3,5) * t330 + Icges(3,6) * t378) * t325;
t387 = t243 + t278;
t377 = pkin(2) * t378;
t376 = pkin(4) * t332;
t371 = Icges(2,4) * t324;
t366 = qJD(2) * t325;
t365 = qJD(3) * t325;
t363 = V_base(5) * qJ(1) + V_base(1);
t359 = qJD(1) + V_base(3);
t356 = t374 * pkin(7);
t302 = t324 * t366 + V_base(4);
t313 = qJD(2) * t374 + V_base(6);
t355 = pkin(2) * t354 - qJ(3) * t325;
t237 = -qJD(4) * t264 + t302;
t269 = -qJD(4) * t282 + t313;
t301 = -t326 * t366 + V_base(5);
t236 = -qJD(4) * t262 + t301;
t295 = pkin(1) * t324 - pkin(7) * t369;
t345 = -t295 * V_base(6) + V_base(5) * t356 + t363;
t296 = pkin(1) * t326 + pkin(7) * t370;
t344 = V_base(4) * t295 - t296 * V_base(5) + t359;
t294 = t325 * t330 * pkin(2) + qJ(3) * t374;
t343 = t301 * t294 + t324 * t365 + t345;
t255 = t324 * t377 + t326 * t355;
t342 = qJD(3) * t374 + t302 * t255 + t344;
t341 = V_base(6) * t296 + V_base(2) + (-t356 - qJ(1)) * V_base(4);
t220 = pkin(3) * t263 - pkin(8) * t262;
t261 = -t283 * pkin(3) - t282 * pkin(8);
t339 = t301 * t261 + (-t220 - t255) * t313 + t343;
t221 = pkin(3) * t265 - pkin(8) * t264;
t256 = -t324 * t355 + t326 * t377;
t338 = t302 * t220 + (-t221 - t256) * t301 + t342;
t337 = t313 * t256 - t326 * t365 + t341;
t205 = pkin(4) * t352 - qJ(5) * t282 - t283 * t376;
t336 = -qJD(5) * t264 + t236 * t205 + t339;
t171 = -pkin(4) * t357 - qJ(5) * t262 + t263 * t376;
t335 = -qJD(5) * t282 + t237 * t171 + t338;
t334 = t313 * t221 + (-t261 - t294) * t302 + t337;
t172 = pkin(4) * t358 - qJ(5) * t264 + t265 * t376;
t333 = -qJD(5) * t262 + t269 * t172 + t334;
t331 = cos(qJ(6));
t328 = sin(qJ(6));
t322 = Icges(2,4) * t326;
t310 = rSges(2,1) * t326 - rSges(2,2) * t324;
t309 = rSges(2,1) * t324 + rSges(2,2) * t326;
t308 = Icges(2,1) * t326 - t371;
t307 = Icges(2,1) * t324 + t322;
t306 = -Icges(2,2) * t324 + t322;
t305 = Icges(2,2) * t326 + t371;
t300 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t299 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t298 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t281 = t374 * rSges(3,3) + (rSges(3,1) * t330 + rSges(3,2) * t378) * t325;
t280 = Icges(3,5) * t374 + (Icges(3,1) * t330 + Icges(3,4) * t378) * t325;
t279 = Icges(3,6) * t374 + (Icges(3,4) * t330 + Icges(3,2) * t378) * t325;
t273 = V_base(5) * rSges(2,3) - t309 * V_base(6) + t363;
t272 = t310 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t268 = t309 * V_base(4) - t310 * V_base(5) + t359;
t254 = rSges(3,1) * t289 + rSges(3,2) * t288 + rSges(3,3) * t370;
t253 = rSges(3,1) * t287 + rSges(3,2) * t286 - rSges(3,3) * t369;
t252 = Icges(3,1) * t289 + Icges(3,4) * t288 + Icges(3,5) * t370;
t251 = Icges(3,1) * t287 + Icges(3,4) * t286 - Icges(3,5) * t369;
t250 = Icges(3,4) * t289 + Icges(3,2) * t288 + Icges(3,6) * t370;
t249 = Icges(3,4) * t287 + Icges(3,2) * t286 - Icges(3,6) * t369;
t246 = -t283 * rSges(4,1) + t282 * rSges(4,2) + rSges(4,3) * t374;
t245 = -Icges(4,1) * t283 + Icges(4,4) * t282 + Icges(4,5) * t374;
t244 = -Icges(4,4) * t283 + Icges(4,2) * t282 + Icges(4,6) * t374;
t229 = t267 * t331 - t282 * t328;
t228 = -t267 * t328 - t282 * t331;
t227 = qJD(6) * t266 + t269;
t226 = pkin(5) * t267 + pkin(9) * t266;
t225 = rSges(5,1) * t271 + rSges(5,2) * t270 - rSges(5,3) * t282;
t224 = Icges(5,1) * t271 + Icges(5,4) * t270 - Icges(5,5) * t282;
t223 = Icges(5,4) * t271 + Icges(5,2) * t270 - Icges(5,6) * t282;
t219 = rSges(4,1) * t265 + rSges(4,2) * t264 + rSges(4,3) * t370;
t218 = rSges(4,1) * t263 + rSges(4,2) * t262 - rSges(4,3) * t369;
t217 = rSges(6,1) * t267 - rSges(6,2) * t266 - rSges(6,3) * t282;
t216 = Icges(4,1) * t265 + Icges(4,4) * t264 + Icges(4,5) * t370;
t215 = Icges(4,1) * t263 + Icges(4,4) * t262 - Icges(4,5) * t369;
t214 = Icges(4,4) * t265 + Icges(4,2) * t264 + Icges(4,6) * t370;
t213 = Icges(4,4) * t263 + Icges(4,2) * t262 - Icges(4,6) * t369;
t210 = Icges(6,1) * t267 - Icges(6,4) * t266 - Icges(6,5) * t282;
t209 = Icges(6,4) * t267 - Icges(6,2) * t266 - Icges(6,6) * t282;
t204 = t233 * t331 - t264 * t328;
t203 = -t233 * t328 - t264 * t331;
t202 = t231 * t331 - t262 * t328;
t201 = -t231 * t328 - t262 * t331;
t200 = qJD(6) * t232 + t237;
t199 = qJD(6) * t230 + t236;
t198 = -t253 * t313 + t281 * t301 + t345;
t197 = t313 * t254 - t302 * t281 + t341;
t196 = pkin(5) * t233 + pkin(9) * t232;
t195 = pkin(5) * t231 + pkin(9) * t230;
t194 = t253 * t302 - t254 * t301 + t344;
t192 = rSges(5,1) * t241 + rSges(5,2) * t240 - rSges(5,3) * t264;
t191 = rSges(5,1) * t239 + rSges(5,2) * t238 - rSges(5,3) * t262;
t190 = Icges(5,1) * t241 + Icges(5,4) * t240 - Icges(5,5) * t264;
t189 = Icges(5,1) * t239 + Icges(5,4) * t238 - Icges(5,5) * t262;
t188 = Icges(5,4) * t241 + Icges(5,2) * t240 - Icges(5,6) * t264;
t187 = Icges(5,4) * t239 + Icges(5,2) * t238 - Icges(5,6) * t262;
t184 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t266;
t183 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t266;
t182 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t266;
t181 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t266;
t180 = rSges(6,1) * t233 - rSges(6,2) * t232 - rSges(6,3) * t264;
t179 = rSges(6,1) * t231 - rSges(6,2) * t230 - rSges(6,3) * t262;
t178 = Icges(6,1) * t233 - Icges(6,4) * t232 - Icges(6,5) * t264;
t177 = Icges(6,1) * t231 - Icges(6,4) * t230 - Icges(6,5) * t262;
t176 = Icges(6,4) * t233 - Icges(6,2) * t232 - Icges(6,6) * t264;
t175 = Icges(6,4) * t231 - Icges(6,2) * t230 - Icges(6,6) * t262;
t168 = rSges(7,1) * t204 + rSges(7,2) * t203 + rSges(7,3) * t232;
t167 = rSges(7,1) * t202 + rSges(7,2) * t201 + rSges(7,3) * t230;
t166 = t246 * t301 + (-t218 - t255) * t313 + t343;
t165 = t313 * t219 + (-t246 - t294) * t302 + t337;
t164 = Icges(7,1) * t204 + Icges(7,4) * t203 + Icges(7,5) * t232;
t163 = Icges(7,1) * t202 + Icges(7,4) * t201 + Icges(7,5) * t230;
t162 = Icges(7,4) * t204 + Icges(7,2) * t203 + Icges(7,6) * t232;
t161 = Icges(7,4) * t202 + Icges(7,2) * t201 + Icges(7,6) * t230;
t160 = Icges(7,5) * t204 + Icges(7,6) * t203 + Icges(7,3) * t232;
t159 = Icges(7,5) * t202 + Icges(7,6) * t201 + Icges(7,3) * t230;
t158 = t218 * t302 + (-t219 - t256) * t301 + t342;
t157 = -t191 * t269 + t225 * t236 + t339;
t156 = t269 * t192 - t237 * t225 + t334;
t155 = t191 * t237 - t192 * t236 + t338;
t154 = t217 * t236 + (-t171 - t179) * t269 + t336;
t153 = t269 * t180 + (-t205 - t217) * t237 + t333;
t152 = t179 * t237 + (-t172 - t180) * t236 + t335;
t151 = (-t171 - t195) * t269 - t167 * t227 + t184 * t199 + t226 * t236 + t336;
t150 = t227 * t168 - t200 * t184 + t269 * t196 + (-t205 - t226) * t237 + t333;
t149 = t167 * t200 - t168 * t199 + t195 * t237 + (-t172 - t196) * t236 + t335;
t1 = m(7) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(1) * (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + m(2) * (t268 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + t227 * ((t160 * t266 + t162 * t228 + t164 * t229) * t200 + (t159 * t266 + t161 * t228 + t163 * t229) * t199 + (t266 * t181 + t228 * t182 + t229 * t183) * t227) / 0.2e1 + t199 * ((t160 * t230 + t162 * t201 + t164 * t202) * t200 + (t230 * t159 + t201 * t161 + t202 * t163) * t199 + (t181 * t230 + t182 * t201 + t183 * t202) * t227) / 0.2e1 + t200 * ((t232 * t160 + t203 * t162 + t204 * t164) * t200 + (t159 * t232 + t161 * t203 + t163 * t204) * t199 + (t181 * t232 + t182 * t203 + t183 * t204) * t227) / 0.2e1 + m(3) * (t194 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(4) * (t158 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(5) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(6) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + ((-t209 * t230 + t210 * t231 + t223 * t238 + t224 * t239 - t262 * t390) * t269 + (-t176 * t230 + t178 * t231 + t188 * t238 + t190 * t239 - t262 * t391) * t237 + (-t175 * t230 + t177 * t231 + t187 * t238 + t189 * t239 - t392 * t262) * t236) * t236 / 0.2e1 + ((-t209 * t232 + t210 * t233 + t223 * t240 + t224 * t241 - t264 * t390) * t269 + (-t176 * t232 + t178 * t233 + t188 * t240 + t190 * t241 - t391 * t264) * t237 + (-t175 * t232 + t177 * t233 + t187 * t240 + t189 * t241 - t264 * t392) * t236) * t237 / 0.2e1 + ((-t209 * t266 + t210 * t267 + t223 * t270 + t224 * t271 - t390 * t282) * t269 + (-t176 * t266 + t178 * t267 + t188 * t270 + t190 * t271 - t282 * t391) * t237 + (-t266 * t175 + t267 * t177 + t187 * t270 + t189 * t271 - t282 * t392) * t236) * t269 / 0.2e1 + ((t244 * t262 + t245 * t263 + t279 * t286 + t280 * t287 - t369 * t387) * t313 + (t214 * t262 + t216 * t263 + t250 * t286 + t252 * t287 - t369 * t388) * t302 + (t213 * t262 + t215 * t263 + t249 * t286 + t251 * t287 - t389 * t369) * t301) * t301 / 0.2e1 + ((t244 * t264 + t245 * t265 + t279 * t288 + t280 * t289 + t370 * t387) * t313 + (t214 * t264 + t216 * t265 + t250 * t288 + t252 * t289 + t388 * t370) * t302 + (t213 * t264 + t215 * t265 + t249 * t288 + t251 * t289 + t370 * t389) * t301) * t302 / 0.2e1 + ((t212 * t374 + t282 * t214 - t283 * t216) * t302 + (t211 * t374 + t282 * t213 - t283 * t215) * t301 + (t243 * t374 + t282 * t244 - t283 * t245) * t313 + ((t250 * t378 + t252 * t330) * t302 + (t249 * t378 + t251 * t330) * t301 + (t279 * t378 + t280 * t330) * t313) * t325 + (t247 * t301 + t248 * t302 + t278 * t313) * t374) * t313 / 0.2e1 + ((-t305 * t324 + t307 * t326 + Icges(1,4)) * V_base(5) + (-t306 * t324 + t308 * t326 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t305 * t326 + t307 * t324 + Icges(1,2)) * V_base(5) + (t306 * t326 + t308 * t324 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t326 - Icges(2,6) * t324 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t324 + Icges(2,6) * t326 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
