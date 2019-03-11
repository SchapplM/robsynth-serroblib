% Calculate kinetic energy for
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:02
% EndTime: 2019-03-09 10:38:08
% DurationCPUTime: 5.94s
% Computational Cost: add. (3305->395), mult. (7754->559), div. (0->0), fcn. (9624->12), ass. (0->178)
t384 = Icges(5,1) + Icges(6,2);
t383 = -Icges(6,1) - Icges(5,3);
t382 = -Icges(5,4) - Icges(6,6);
t381 = Icges(6,4) - Icges(5,5);
t380 = Icges(6,5) - Icges(5,6);
t379 = Icges(5,2) + Icges(6,3);
t318 = sin(qJ(2));
t354 = sin(pkin(11));
t356 = cos(pkin(6));
t336 = t356 * t354;
t355 = cos(pkin(11));
t337 = t356 * t355;
t360 = cos(qJ(2));
t276 = t318 * t337 + t336 * t360;
t285 = -t318 * t354 + t360 * t355;
t319 = sin(qJ(1));
t321 = cos(qJ(1));
t257 = t276 * t321 + t319 * t285;
t316 = sin(pkin(6));
t359 = cos(qJ(4));
t345 = t316 * t359;
t358 = sin(qJ(4));
t233 = t257 * t358 + t321 * t345;
t344 = t316 * t358;
t234 = t257 * t359 - t321 * t344;
t327 = -t318 * t336 + t337 * t360;
t331 = t318 * t355 + t354 * t360;
t256 = -t319 * t331 + t321 * t327;
t378 = t233 * t379 + t234 * t382 - t256 * t380;
t259 = -t319 * t276 + t285 * t321;
t235 = t259 * t358 - t319 * t345;
t236 = t259 * t359 + t319 * t344;
t258 = -t319 * t327 - t321 * t331;
t377 = t235 * t379 + t236 * t382 - t258 * t380;
t376 = t233 * t380 - t234 * t381 + t256 * t383;
t375 = t235 * t380 - t236 * t381 + t258 * t383;
t374 = t382 * t233 + t234 * t384 + t381 * t256;
t373 = t382 * t235 + t236 * t384 + t381 * t258;
t351 = t316 * t321;
t205 = Icges(4,5) * t257 + Icges(4,6) * t256 - Icges(4,3) * t351;
t339 = t356 * t360;
t279 = -t319 * t318 + t321 * t339;
t341 = t318 * t356;
t280 = t319 * t360 + t321 * t341;
t243 = Icges(3,5) * t280 + Icges(3,6) * t279 - Icges(3,3) * t351;
t372 = t205 + t243;
t352 = t316 * t319;
t206 = Icges(4,5) * t259 + Icges(4,6) * t258 + Icges(4,3) * t352;
t281 = -t321 * t318 - t319 * t339;
t282 = -t319 * t341 + t321 * t360;
t244 = Icges(3,5) * t282 + Icges(3,6) * t281 + Icges(3,3) * t352;
t371 = t206 + t244;
t275 = t331 * t316;
t262 = t275 * t358 - t356 * t359;
t263 = t275 * t359 + t356 * t358;
t274 = t285 * t316;
t370 = t262 * t379 + t263 * t382 - t274 * t380;
t369 = t262 * t380 - t263 * t381 + t274 * t383;
t368 = t382 * t262 + t263 * t384 + t381 * t274;
t239 = Icges(4,5) * t275 + Icges(4,6) * t274 + Icges(4,3) * t356;
t270 = Icges(3,3) * t356 + (Icges(3,5) * t318 + Icges(3,6) * t360) * t316;
t367 = t239 + t270;
t357 = pkin(2) * t360;
t353 = Icges(2,4) * t319;
t350 = qJD(2) * t316;
t349 = qJD(3) * t316;
t348 = V_base(5) * pkin(7) + V_base(1);
t343 = t356 * pkin(8);
t293 = t319 * t350 + V_base(4);
t313 = V_base(6) + qJD(1);
t342 = pkin(2) * t341 - qJ(3) * t316;
t232 = -qJD(4) * t258 + t293;
t294 = qJD(2) * t356 + t313;
t261 = -qJD(4) * t274 + t294;
t292 = -t321 * t350 + V_base(5);
t287 = t319 * pkin(1) - pkin(8) * t351;
t335 = -t287 * t313 + V_base(5) * t343 + t348;
t288 = pkin(1) * t321 + pkin(8) * t352;
t334 = V_base(4) * t287 - t288 * V_base(5) + V_base(3);
t231 = -qJD(4) * t256 + t292;
t286 = t316 * t318 * pkin(2) + qJ(3) * t356;
t333 = t292 * t286 + t319 * t349 + t335;
t254 = t319 * t357 + t321 * t342;
t332 = qJD(3) * t356 + t293 * t254 + t334;
t330 = t313 * t288 + V_base(2) + (-t343 - pkin(7)) * V_base(4);
t221 = pkin(3) * t257 - pkin(9) * t256;
t251 = pkin(3) * t275 - pkin(9) * t274;
t329 = t292 * t251 + (-t221 - t254) * t294 + t333;
t222 = pkin(3) * t259 - pkin(9) * t258;
t255 = -t319 * t342 + t321 * t357;
t328 = t293 * t221 + (-t222 - t255) * t292 + t332;
t326 = t294 * t255 - t321 * t349 + t330;
t223 = pkin(4) * t263 + qJ(5) * t262;
t325 = qJD(5) * t235 + t231 * t223 + t329;
t193 = pkin(4) * t234 + qJ(5) * t233;
t324 = qJD(5) * t262 + t232 * t193 + t328;
t323 = t294 * t222 + (-t251 - t286) * t293 + t326;
t194 = pkin(4) * t236 + qJ(5) * t235;
t322 = qJD(5) * t233 + t261 * t194 + t323;
t320 = cos(qJ(6));
t317 = sin(qJ(6));
t314 = Icges(2,4) * t321;
t302 = rSges(2,1) * t321 - t319 * rSges(2,2);
t301 = t319 * rSges(2,1) + rSges(2,2) * t321;
t300 = Icges(2,1) * t321 - t353;
t299 = Icges(2,1) * t319 + t314;
t298 = -Icges(2,2) * t319 + t314;
t297 = Icges(2,2) * t321 + t353;
t291 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t290 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t289 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t273 = t356 * rSges(3,3) + (rSges(3,1) * t318 + rSges(3,2) * t360) * t316;
t272 = Icges(3,5) * t356 + (Icges(3,1) * t318 + Icges(3,4) * t360) * t316;
t271 = Icges(3,6) * t356 + (Icges(3,4) * t318 + Icges(3,2) * t360) * t316;
t266 = V_base(5) * rSges(2,3) - t301 * t313 + t348;
t265 = t302 * t313 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t264 = t301 * V_base(4) - t302 * V_base(5) + V_base(3);
t250 = rSges(3,1) * t282 + rSges(3,2) * t281 + rSges(3,3) * t352;
t249 = t280 * rSges(3,1) + t279 * rSges(3,2) - rSges(3,3) * t351;
t248 = Icges(3,1) * t282 + Icges(3,4) * t281 + Icges(3,5) * t352;
t247 = Icges(3,1) * t280 + Icges(3,4) * t279 - Icges(3,5) * t351;
t246 = Icges(3,4) * t282 + Icges(3,2) * t281 + Icges(3,6) * t352;
t245 = Icges(3,4) * t280 + Icges(3,2) * t279 - Icges(3,6) * t351;
t242 = t275 * rSges(4,1) + t274 * rSges(4,2) + rSges(4,3) * t356;
t241 = Icges(4,1) * t275 + Icges(4,4) * t274 + Icges(4,5) * t356;
t240 = Icges(4,4) * t275 + Icges(4,2) * t274 + Icges(4,6) * t356;
t237 = -pkin(5) * t274 + pkin(10) * t263;
t226 = t262 * t317 - t274 * t320;
t225 = t262 * t320 + t274 * t317;
t224 = qJD(6) * t263 + t261;
t220 = rSges(5,1) * t263 - rSges(5,2) * t262 - rSges(5,3) * t274;
t219 = -rSges(6,1) * t274 - rSges(6,2) * t263 + rSges(6,3) * t262;
t212 = rSges(4,1) * t259 + rSges(4,2) * t258 + rSges(4,3) * t352;
t211 = t257 * rSges(4,1) + t256 * rSges(4,2) - rSges(4,3) * t351;
t210 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t352;
t209 = Icges(4,1) * t257 + Icges(4,4) * t256 - Icges(4,5) * t351;
t208 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t352;
t207 = Icges(4,4) * t257 + Icges(4,2) * t256 - Icges(4,6) * t351;
t202 = -pkin(5) * t258 + pkin(10) * t236;
t201 = -pkin(5) * t256 + pkin(10) * t234;
t200 = t235 * t317 - t258 * t320;
t199 = t235 * t320 + t258 * t317;
t198 = t233 * t317 - t256 * t320;
t197 = t233 * t320 + t256 * t317;
t196 = qJD(6) * t236 + t232;
t195 = qJD(6) * t234 + t231;
t192 = -t249 * t294 + t273 * t292 + t335;
t191 = t294 * t250 - t293 * t273 + t330;
t188 = t249 * t293 - t250 * t292 + t334;
t187 = rSges(7,1) * t226 + rSges(7,2) * t225 + rSges(7,3) * t263;
t186 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t263;
t185 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t263;
t184 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t263;
t183 = rSges(5,1) * t236 - rSges(5,2) * t235 - rSges(5,3) * t258;
t182 = rSges(5,1) * t234 - rSges(5,2) * t233 - rSges(5,3) * t256;
t181 = -rSges(6,1) * t258 - rSges(6,2) * t236 + rSges(6,3) * t235;
t180 = -rSges(6,1) * t256 - rSges(6,2) * t234 + rSges(6,3) * t233;
t166 = rSges(7,1) * t200 + rSges(7,2) * t199 + rSges(7,3) * t236;
t165 = rSges(7,1) * t198 + rSges(7,2) * t197 + rSges(7,3) * t234;
t164 = Icges(7,1) * t200 + Icges(7,4) * t199 + Icges(7,5) * t236;
t163 = Icges(7,1) * t198 + Icges(7,4) * t197 + Icges(7,5) * t234;
t162 = Icges(7,4) * t200 + Icges(7,2) * t199 + Icges(7,6) * t236;
t161 = Icges(7,4) * t198 + Icges(7,2) * t197 + Icges(7,6) * t234;
t160 = Icges(7,5) * t200 + Icges(7,6) * t199 + Icges(7,3) * t236;
t159 = Icges(7,5) * t198 + Icges(7,6) * t197 + Icges(7,3) * t234;
t158 = t242 * t292 + (-t211 - t254) * t294 + t333;
t157 = t294 * t212 + (-t242 - t286) * t293 + t326;
t156 = t211 * t293 + (-t212 - t255) * t292 + t332;
t155 = -t182 * t261 + t220 * t231 + t329;
t154 = t261 * t183 - t232 * t220 + t323;
t153 = t182 * t232 - t183 * t231 + t328;
t152 = t219 * t231 + (-t180 - t193) * t261 + t325;
t151 = t261 * t181 + (-t219 - t223) * t232 + t322;
t150 = t180 * t232 + (-t181 - t194) * t231 + t324;
t149 = t325 - t165 * t224 + t187 * t195 + t231 * t237 + (-t193 - t201) * t261;
t148 = t224 * t166 - t196 * t187 + t261 * t202 + (-t223 - t237) * t232 + t322;
t147 = t165 * t196 - t166 * t195 + t201 * t232 + (-t194 - t202) * t231 + t324;
t1 = m(4) * (t156 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(5) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(7) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t188 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + t195 * ((t160 * t234 + t162 * t197 + t164 * t198) * t196 + (t159 * t234 + t161 * t197 + t163 * t198) * t195 + (t184 * t234 + t185 * t197 + t186 * t198) * t224) / 0.2e1 + t196 * ((t160 * t236 + t162 * t199 + t164 * t200) * t196 + (t159 * t236 + t161 * t199 + t163 * t200) * t195 + (t184 * t236 + t185 * t199 + t186 * t200) * t224) / 0.2e1 + t224 * ((t160 * t263 + t162 * t225 + t164 * t226) * t196 + (t159 * t263 + t161 * t225 + t163 * t226) * t195 + (t184 * t263 + t185 * t225 + t186 * t226) * t224) / 0.2e1 + m(2) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(1) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + ((t233 * t370 + t234 * t368 - t256 * t369) * t261 + (t233 * t377 + t234 * t373 - t256 * t375) * t232 + (t378 * t233 + t374 * t234 - t376 * t256) * t231) * t231 / 0.2e1 + ((t235 * t370 + t236 * t368 - t258 * t369) * t261 + (t377 * t235 + t373 * t236 - t375 * t258) * t232 + (t235 * t378 + t374 * t236 - t376 * t258) * t231) * t232 / 0.2e1 + ((t370 * t262 + t368 * t263 - t369 * t274) * t261 + (t262 * t377 + t263 * t373 - t274 * t375) * t232 + (t262 * t378 + t374 * t263 - t376 * t274) * t231) * t261 / 0.2e1 + ((t256 * t240 + t257 * t241 + t279 * t271 + t280 * t272 - t351 * t367) * t294 + (t256 * t208 + t257 * t210 + t279 * t246 + t280 * t248 - t351 * t371) * t293 + (t256 * t207 + t257 * t209 + t279 * t245 + t280 * t247 - t372 * t351) * t292) * t292 / 0.2e1 + ((t240 * t258 + t241 * t259 + t271 * t281 + t272 * t282 + t352 * t367) * t294 + (t208 * t258 + t210 * t259 + t246 * t281 + t248 * t282 + t371 * t352) * t293 + (t207 * t258 + t209 * t259 + t245 * t281 + t247 * t282 + t352 * t372) * t292) * t293 / 0.2e1 + (((t246 * t360 + t248 * t318) * t293 + (t245 * t360 + t247 * t318) * t292 + (t271 * t360 + t272 * t318) * t294) * t316 + (t243 * t292 + t244 * t293 + t270 * t294) * t356 + (t206 * t356 + t274 * t208 + t275 * t210) * t293 + (t205 * t356 + t274 * t207 + t275 * t209) * t292 + (t239 * t356 + t274 * t240 + t275 * t241) * t294) * t294 / 0.2e1 + ((-t319 * t297 + t299 * t321 + Icges(1,4)) * V_base(5) + (-t319 * t298 + t300 * t321 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t297 * t321 + t319 * t299 + Icges(1,2)) * V_base(5) + (t298 * t321 + t319 * t300 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t319 + Icges(2,6) * t321) * V_base(5) + (Icges(2,5) * t321 - Icges(2,6) * t319) * V_base(4) + Icges(2,3) * t313 / 0.2e1) * t313;
T  = t1;
