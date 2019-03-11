% Calculate kinetic energy for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:05
% EndTime: 2019-03-09 09:01:09
% DurationCPUTime: 4.12s
% Computational Cost: add. (2981->387), mult. (6935->544), div. (0->0), fcn. (8508->12), ass. (0->175)
t374 = Icges(4,1) + Icges(5,2);
t373 = Icges(4,4) + Icges(5,6);
t372 = Icges(5,4) - Icges(4,5);
t371 = Icges(5,5) - Icges(4,6);
t370 = -Icges(4,2) - Icges(5,3);
t369 = Icges(5,1) + Icges(4,3) + Icges(3,3);
t313 = sin(qJ(1));
t316 = cos(qJ(1));
t312 = sin(qJ(2));
t315 = cos(qJ(2));
t351 = sin(pkin(11));
t352 = cos(pkin(11));
t280 = -t312 * t351 + t315 * t352;
t309 = cos(pkin(6));
t325 = t309 * t280;
t329 = t312 * t352 + t315 * t351;
t250 = -t313 * t325 - t316 * t329;
t324 = t309 * t329;
t251 = t280 * t316 - t313 * t324;
t308 = sin(pkin(6));
t349 = t308 * t313;
t368 = t373 * t250 + t374 * t251 - t372 * t349;
t248 = -t313 * t329 + t316 * t325;
t249 = t313 * t280 + t316 * t324;
t348 = t308 * t316;
t367 = t373 * t248 + t374 * t249 + t372 * t348;
t366 = t370 * t248 - t373 * t249 - t371 * t348;
t365 = t370 * t250 - t373 * t251 + t371 * t349;
t270 = t280 * t308;
t271 = t329 * t308;
t364 = t373 * t270 + t374 * t271 - t372 * t309;
t363 = t370 * t270 - t373 * t271 + t371 * t309;
t345 = t313 * t315;
t347 = t312 * t316;
t276 = -t309 * t345 - t347;
t344 = t315 * t316;
t346 = t313 * t312;
t277 = -t309 * t346 + t344;
t362 = Icges(3,5) * t277 + Icges(3,6) * t276 - t371 * t250 - t372 * t251 + t369 * t349;
t274 = t309 * t344 - t346;
t275 = t309 * t347 + t345;
t361 = Icges(3,5) * t275 + Icges(3,6) * t274 - t371 * t248 - t372 * t249 - t369 * t348;
t360 = (Icges(3,5) * t312 + Icges(3,6) * t315) * t308 - t372 * t271 - t371 * t270 + t369 * t309;
t356 = cos(qJ(5));
t355 = pkin(2) * t312;
t354 = pkin(8) * t309;
t353 = pkin(2) * t315;
t350 = Icges(2,4) * t313;
t205 = pkin(3) * t249 - qJ(4) * t248;
t334 = -qJ(3) * t308 + t309 * t355;
t246 = t313 * t353 + t316 * t334;
t343 = -t205 - t246;
t206 = pkin(3) * t251 - qJ(4) * t250;
t247 = -t313 * t334 + t316 * t353;
t342 = -t206 - t247;
t241 = pkin(3) * t271 - qJ(4) * t270;
t281 = qJ(3) * t309 + t308 * t355;
t341 = -t241 - t281;
t340 = qJD(2) * t308;
t339 = qJD(3) * t308;
t338 = V_base(5) * pkin(7) + V_base(1);
t335 = t308 * t356;
t288 = t313 * t340 + V_base(4);
t305 = V_base(6) + qJD(1);
t216 = qJD(5) * t251 + t288;
t289 = qJD(2) * t309 + t305;
t253 = qJD(5) * t271 + t289;
t287 = -t316 * t340 + V_base(5);
t282 = t313 * pkin(1) - pkin(8) * t348;
t331 = -t282 * t305 + V_base(5) * t354 + t338;
t283 = pkin(1) * t316 + pkin(8) * t349;
t330 = V_base(4) * t282 - t283 * V_base(5) + V_base(3);
t215 = qJD(5) * t249 + t287;
t328 = t287 * t281 + t313 * t339 + t331;
t327 = qJD(3) * t309 + t288 * t246 + t330;
t326 = t305 * t283 + V_base(2) + (-pkin(7) - t354) * V_base(4);
t323 = -qJD(4) * t250 + t287 * t241 + t328;
t322 = -qJD(4) * t270 + t288 * t205 + t327;
t321 = t289 * t247 - t316 * t339 + t326;
t320 = -qJD(4) * t248 + t289 * t206 + t321;
t223 = -pkin(4) * t348 + t249 * pkin(9);
t259 = pkin(4) * t309 + pkin(9) * t271;
t319 = t287 * t259 + (-t223 + t343) * t289 + t323;
t222 = pkin(4) * t349 + pkin(9) * t251;
t318 = t288 * t223 + (-t222 + t342) * t287 + t322;
t317 = t289 * t222 + (-t259 + t341) * t288 + t320;
t314 = cos(qJ(6));
t311 = sin(qJ(5));
t310 = sin(qJ(6));
t306 = Icges(2,4) * t316;
t297 = rSges(2,1) * t316 - t313 * rSges(2,2);
t296 = t313 * rSges(2,1) + rSges(2,2) * t316;
t295 = Icges(2,1) * t316 - t350;
t294 = Icges(2,1) * t313 + t306;
t293 = -Icges(2,2) * t313 + t306;
t292 = Icges(2,2) * t316 + t350;
t286 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t285 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t284 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t269 = rSges(3,3) * t309 + (rSges(3,1) * t312 + rSges(3,2) * t315) * t308;
t268 = Icges(3,5) * t309 + (Icges(3,1) * t312 + Icges(3,4) * t315) * t308;
t267 = Icges(3,6) * t309 + (Icges(3,4) * t312 + Icges(3,2) * t315) * t308;
t258 = V_base(5) * rSges(2,3) - t296 * t305 + t338;
t257 = t297 * t305 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t256 = t296 * V_base(4) - t297 * V_base(5) + V_base(3);
t255 = -t270 * t311 + t309 * t356;
t254 = t270 * t356 + t309 * t311;
t240 = rSges(3,1) * t277 + rSges(3,2) * t276 + rSges(3,3) * t349;
t239 = t275 * rSges(3,1) + t274 * rSges(3,2) - rSges(3,3) * t348;
t237 = Icges(3,1) * t277 + Icges(3,4) * t276 + Icges(3,5) * t349;
t236 = Icges(3,1) * t275 + Icges(3,4) * t274 - Icges(3,5) * t348;
t235 = Icges(3,4) * t277 + Icges(3,2) * t276 + Icges(3,6) * t349;
t234 = Icges(3,4) * t275 + Icges(3,2) * t274 - Icges(3,6) * t348;
t231 = rSges(4,1) * t271 + rSges(4,2) * t270 + rSges(4,3) * t309;
t230 = rSges(5,1) * t309 - rSges(5,2) * t271 - rSges(5,3) * t270;
t220 = -t248 * t311 - t316 * t335;
t219 = -t248 * t356 + t311 * t348;
t218 = -t250 * t311 + t313 * t335;
t217 = t250 * t356 + t311 * t349;
t212 = t255 * t314 + t271 * t310;
t211 = -t255 * t310 + t271 * t314;
t208 = qJD(6) * t254 + t253;
t207 = pkin(5) * t255 + pkin(10) * t254;
t204 = rSges(6,1) * t255 - rSges(6,2) * t254 + rSges(6,3) * t271;
t203 = Icges(6,1) * t255 - Icges(6,4) * t254 + Icges(6,5) * t271;
t202 = Icges(6,4) * t255 - Icges(6,2) * t254 + Icges(6,6) * t271;
t201 = Icges(6,5) * t255 - Icges(6,6) * t254 + Icges(6,3) * t271;
t200 = rSges(4,1) * t251 + rSges(4,2) * t250 + rSges(4,3) * t349;
t199 = t249 * rSges(4,1) + t248 * rSges(4,2) - rSges(4,3) * t348;
t198 = -rSges(5,1) * t348 - t249 * rSges(5,2) - t248 * rSges(5,3);
t197 = rSges(5,1) * t349 - rSges(5,2) * t251 - rSges(5,3) * t250;
t182 = t220 * t314 + t249 * t310;
t181 = -t220 * t310 + t249 * t314;
t180 = t218 * t314 + t251 * t310;
t179 = -t218 * t310 + t251 * t314;
t178 = qJD(6) * t217 + t216;
t177 = -qJD(6) * t219 + t215;
t176 = pkin(5) * t220 - pkin(10) * t219;
t175 = pkin(5) * t218 + pkin(10) * t217;
t174 = -t239 * t289 + t269 * t287 + t331;
t173 = t240 * t289 - t269 * t288 + t326;
t172 = t239 * t288 - t240 * t287 + t330;
t171 = rSges(7,1) * t212 + rSges(7,2) * t211 + rSges(7,3) * t254;
t170 = Icges(7,1) * t212 + Icges(7,4) * t211 + Icges(7,5) * t254;
t169 = Icges(7,4) * t212 + Icges(7,2) * t211 + Icges(7,6) * t254;
t168 = Icges(7,5) * t212 + Icges(7,6) * t211 + Icges(7,3) * t254;
t167 = rSges(6,1) * t220 + rSges(6,2) * t219 + rSges(6,3) * t249;
t166 = rSges(6,1) * t218 - rSges(6,2) * t217 + rSges(6,3) * t251;
t165 = Icges(6,1) * t220 + Icges(6,4) * t219 + Icges(6,5) * t249;
t164 = Icges(6,1) * t218 - Icges(6,4) * t217 + Icges(6,5) * t251;
t163 = Icges(6,4) * t220 + Icges(6,2) * t219 + Icges(6,6) * t249;
t162 = Icges(6,4) * t218 - Icges(6,2) * t217 + Icges(6,6) * t251;
t161 = Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t249;
t160 = Icges(6,5) * t218 - Icges(6,6) * t217 + Icges(6,3) * t251;
t159 = rSges(7,1) * t182 + rSges(7,2) * t181 - rSges(7,3) * t219;
t158 = rSges(7,1) * t180 + rSges(7,2) * t179 + rSges(7,3) * t217;
t157 = Icges(7,1) * t182 + Icges(7,4) * t181 - Icges(7,5) * t219;
t156 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t217;
t155 = Icges(7,4) * t182 + Icges(7,2) * t181 - Icges(7,6) * t219;
t154 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t217;
t153 = Icges(7,5) * t182 + Icges(7,6) * t181 - Icges(7,3) * t219;
t152 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t217;
t151 = t231 * t287 + (-t199 - t246) * t289 + t328;
t150 = t289 * t200 + (-t231 - t281) * t288 + t321;
t149 = t199 * t288 + (-t200 - t247) * t287 + t327;
t148 = t230 * t287 + (-t198 + t343) * t289 + t323;
t147 = t289 * t197 + (-t230 + t341) * t288 + t320;
t146 = t198 * t288 + (-t197 + t342) * t287 + t322;
t145 = -t167 * t253 + t204 * t215 + t319;
t144 = t253 * t166 - t216 * t204 + t317;
t143 = -t166 * t215 + t167 * t216 + t318;
t142 = -t159 * t208 + t171 * t177 - t176 * t253 + t207 * t215 + t319;
t141 = t208 * t158 - t178 * t171 + t253 * t175 - t216 * t207 + t317;
t140 = -t158 * t177 + t159 * t178 - t175 * t215 + t176 * t216 + t318;
t1 = m(7) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(5) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(4) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(3) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + t178 * ((t217 * t152 + t179 * t154 + t156 * t180) * t178 + (t153 * t217 + t155 * t179 + t157 * t180) * t177 + (t168 * t217 + t169 * t179 + t170 * t180) * t208) / 0.2e1 + t177 * ((-t152 * t219 + t154 * t181 + t156 * t182) * t178 + (-t219 * t153 + t181 * t155 + t182 * t157) * t177 + (-t168 * t219 + t169 * t181 + t170 * t182) * t208) / 0.2e1 + t215 * ((t160 * t249 + t162 * t219 + t164 * t220) * t216 + (t249 * t161 + t219 * t163 + t220 * t165) * t215 + (t201 * t249 + t202 * t219 + t203 * t220) * t253) / 0.2e1 + t216 * ((t160 * t251 - t217 * t162 + t164 * t218) * t216 + (t161 * t251 - t163 * t217 + t165 * t218) * t215 + (t201 * t251 - t202 * t217 + t203 * t218) * t253) / 0.2e1 + t208 * ((t152 * t254 + t154 * t211 + t156 * t212) * t178 + (t153 * t254 + t155 * t211 + t157 * t212) * t177 + (t168 * t254 + t169 * t211 + t212 * t170) * t208) / 0.2e1 + m(2) * (t256 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t253 * ((t160 * t271 - t162 * t254 + t164 * t255) * t216 + (t161 * t271 - t163 * t254 + t165 * t255) * t215 + (t271 * t201 - t202 * t254 + t255 * t203) * t253) / 0.2e1 + m(1) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + ((-t313 * t292 + t294 * t316 + Icges(1,4)) * V_base(5) + (-t313 * t293 + t316 * t295 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t316 * t292 + t313 * t294 + Icges(1,2)) * V_base(5) + (t293 * t316 + t313 * t295 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t248 * t363 + t249 * t364 + t274 * t267 + t275 * t268 - t348 * t360) * t289 + (t274 * t235 + t275 * t237 - t365 * t248 + t249 * t368 - t362 * t348) * t288 + (t274 * t234 + t275 * t236 - t366 * t248 + t367 * t249 - t361 * t348) * t287) * t287 / 0.2e1 + ((-t250 * t363 + t251 * t364 + t267 * t276 + t268 * t277 + t349 * t360) * t289 + (t276 * t235 + t277 * t237 - t365 * t250 + t368 * t251 + t362 * t349) * t288 + (t234 * t276 + t236 * t277 - t250 * t366 + t251 * t367 + t349 * t361) * t287) * t288 / 0.2e1 + (((t267 * t315 + t268 * t312) * t308 + t364 * t271 - t363 * t270 + t360 * t309) * t289 + ((t235 * t315 + t237 * t312) * t308 + t368 * t271 - t365 * t270 + t362 * t309) * t288 + ((t234 * t315 + t236 * t312) * t308 + t367 * t271 - t366 * t270 + t361 * t309) * t287) * t289 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t313 + Icges(2,6) * t316) * V_base(5) + (Icges(2,5) * t316 - Icges(2,6) * t313) * V_base(4) + Icges(2,3) * t305 / 0.2e1) * t305;
T  = t1;
