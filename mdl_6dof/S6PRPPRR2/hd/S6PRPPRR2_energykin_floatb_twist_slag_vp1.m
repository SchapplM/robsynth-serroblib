% Calculate kinetic energy for
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:36
% EndTime: 2019-03-08 19:17:40
% DurationCPUTime: 3.99s
% Computational Cost: add. (2921->387), mult. (6935->543), div. (0->0), fcn. (8508->12), ass. (0->173)
t370 = Icges(4,1) + Icges(5,2);
t369 = Icges(4,4) + Icges(5,6);
t368 = Icges(5,4) - Icges(4,5);
t367 = Icges(5,5) - Icges(4,6);
t366 = -Icges(4,2) - Icges(5,3);
t365 = Icges(5,1) + Icges(4,3) + Icges(3,3);
t302 = sin(pkin(10));
t304 = cos(pkin(10));
t308 = sin(qJ(2));
t310 = cos(qJ(2));
t346 = sin(pkin(11));
t347 = cos(pkin(11));
t275 = -t308 * t346 + t310 * t347;
t305 = cos(pkin(6));
t319 = t305 * t275;
t323 = t308 * t347 + t310 * t346;
t245 = -t302 * t319 - t304 * t323;
t318 = t305 * t323;
t246 = t275 * t304 - t302 * t318;
t303 = sin(pkin(6));
t344 = t302 * t303;
t364 = t369 * t245 + t370 * t246 - t368 * t344;
t243 = -t302 * t323 + t304 * t319;
t244 = t302 * t275 + t304 * t318;
t343 = t303 * t304;
t363 = t369 * t243 + t370 * t244 + t368 * t343;
t362 = t366 * t243 - t369 * t244 - t367 * t343;
t361 = t366 * t245 - t369 * t246 + t367 * t344;
t265 = t275 * t303;
t266 = t323 * t303;
t360 = t369 * t265 + t370 * t266 - t368 * t305;
t359 = t366 * t265 - t369 * t266 + t367 * t305;
t340 = t305 * t310;
t270 = -t302 * t340 - t304 * t308;
t341 = t305 * t308;
t271 = -t302 * t341 + t304 * t310;
t356 = Icges(3,5) * t271 + Icges(3,6) * t270 - t367 * t245 - t368 * t246 + t365 * t344;
t268 = -t302 * t308 + t304 * t340;
t269 = t302 * t310 + t304 * t341;
t355 = Icges(3,5) * t269 + Icges(3,6) * t268 - t367 * t243 - t368 * t244 - t365 * t343;
t354 = (Icges(3,5) * t308 + Icges(3,6) * t310) * t303 - t368 * t266 - t367 * t265 + t365 * t305;
t350 = cos(qJ(5));
t349 = pkin(7) * t305;
t348 = pkin(2) * t310;
t345 = Icges(2,4) * t302;
t307 = sin(qJ(5));
t342 = t303 * t307;
t196 = pkin(3) * t244 - qJ(4) * t243;
t328 = pkin(2) * t341 - qJ(3) * t303;
t236 = t302 * t348 + t304 * t328;
t339 = -t196 - t236;
t197 = pkin(3) * t246 - qJ(4) * t245;
t237 = -t302 * t328 + t304 * t348;
t338 = -t197 - t237;
t238 = pkin(3) * t266 - qJ(4) * t265;
t276 = pkin(2) * t303 * t308 + qJ(3) * t305;
t337 = -t238 - t276;
t336 = qJD(2) * t303;
t335 = qJD(3) * t303;
t334 = V_base(5) * qJ(1) + V_base(1);
t330 = qJD(1) + V_base(3);
t329 = t303 * t350;
t283 = t302 * t336 + V_base(4);
t293 = qJD(2) * t305 + V_base(6);
t211 = qJD(5) * t246 + t283;
t249 = qJD(5) * t266 + t293;
t282 = -t304 * t336 + V_base(5);
t210 = qJD(5) * t244 + t282;
t277 = pkin(1) * t302 - pkin(7) * t343;
t325 = -t277 * V_base(6) + V_base(5) * t349 + t334;
t278 = pkin(1) * t304 + pkin(7) * t344;
t324 = V_base(4) * t277 - t278 * V_base(5) + t330;
t322 = V_base(6) * t278 + V_base(2) + (-qJ(1) - t349) * V_base(4);
t321 = t282 * t276 + t302 * t335 + t325;
t320 = qJD(3) * t305 + t283 * t236 + t324;
t317 = -qJD(4) * t245 + t282 * t238 + t321;
t316 = -qJD(4) * t265 + t283 * t196 + t320;
t315 = t293 * t237 - t304 * t335 + t322;
t314 = -qJD(4) * t243 + t293 * t197 + t315;
t218 = -pkin(4) * t343 + pkin(8) * t244;
t252 = pkin(4) * t305 + pkin(8) * t266;
t313 = t282 * t252 + (-t218 + t339) * t293 + t317;
t217 = pkin(4) * t344 + pkin(8) * t246;
t312 = t283 * t218 + (-t217 + t338) * t282 + t316;
t311 = t293 * t217 + (-t252 + t337) * t283 + t314;
t309 = cos(qJ(6));
t306 = sin(qJ(6));
t300 = Icges(2,4) * t304;
t291 = rSges(2,1) * t304 - rSges(2,2) * t302;
t290 = rSges(2,1) * t302 + rSges(2,2) * t304;
t289 = Icges(2,1) * t304 - t345;
t288 = Icges(2,1) * t302 + t300;
t287 = -Icges(2,2) * t302 + t300;
t286 = Icges(2,2) * t304 + t345;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t264 = t305 * rSges(3,3) + (rSges(3,1) * t308 + rSges(3,2) * t310) * t303;
t263 = Icges(3,5) * t305 + (Icges(3,1) * t308 + Icges(3,4) * t310) * t303;
t262 = Icges(3,6) * t305 + (Icges(3,4) * t308 + Icges(3,2) * t310) * t303;
t254 = V_base(5) * rSges(2,3) - t290 * V_base(6) + t334;
t253 = t291 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t251 = -t265 * t307 + t305 * t350;
t250 = t265 * t350 + t305 * t307;
t248 = t290 * V_base(4) - t291 * V_base(5) + t330;
t234 = rSges(3,1) * t271 + rSges(3,2) * t270 + rSges(3,3) * t344;
t233 = rSges(3,1) * t269 + rSges(3,2) * t268 - rSges(3,3) * t343;
t232 = Icges(3,1) * t271 + Icges(3,4) * t270 + Icges(3,5) * t344;
t231 = Icges(3,1) * t269 + Icges(3,4) * t268 - Icges(3,5) * t343;
t230 = Icges(3,4) * t271 + Icges(3,2) * t270 + Icges(3,6) * t344;
t229 = Icges(3,4) * t269 + Icges(3,2) * t268 - Icges(3,6) * t343;
t226 = rSges(4,1) * t266 + rSges(4,2) * t265 + rSges(4,3) * t305;
t225 = rSges(5,1) * t305 - rSges(5,2) * t266 - rSges(5,3) * t265;
t215 = -t243 * t307 - t304 * t329;
t214 = -t243 * t350 + t304 * t342;
t213 = -t245 * t307 + t302 * t329;
t212 = t245 * t350 + t302 * t342;
t207 = t251 * t309 + t266 * t306;
t206 = -t251 * t306 + t266 * t309;
t204 = qJD(6) * t250 + t249;
t202 = pkin(5) * t251 + pkin(9) * t250;
t201 = rSges(6,1) * t251 - rSges(6,2) * t250 + rSges(6,3) * t266;
t200 = Icges(6,1) * t251 - Icges(6,4) * t250 + Icges(6,5) * t266;
t199 = Icges(6,4) * t251 - Icges(6,2) * t250 + Icges(6,6) * t266;
t198 = Icges(6,5) * t251 - Icges(6,6) * t250 + Icges(6,3) * t266;
t195 = rSges(4,1) * t246 + rSges(4,2) * t245 + rSges(4,3) * t344;
t194 = rSges(4,1) * t244 + rSges(4,2) * t243 - rSges(4,3) * t343;
t193 = -rSges(5,1) * t343 - rSges(5,2) * t244 - rSges(5,3) * t243;
t192 = rSges(5,1) * t344 - rSges(5,2) * t246 - rSges(5,3) * t245;
t177 = t215 * t309 + t244 * t306;
t176 = -t215 * t306 + t244 * t309;
t175 = t213 * t309 + t246 * t306;
t174 = -t213 * t306 + t246 * t309;
t173 = qJD(6) * t212 + t211;
t172 = -qJD(6) * t214 + t210;
t171 = pkin(5) * t215 - pkin(9) * t214;
t170 = pkin(5) * t213 + pkin(9) * t212;
t169 = -t233 * t293 + t264 * t282 + t325;
t168 = t234 * t293 - t264 * t283 + t322;
t167 = rSges(7,1) * t207 + rSges(7,2) * t206 + rSges(7,3) * t250;
t166 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t250;
t165 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t250;
t164 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t250;
t163 = t233 * t283 - t234 * t282 + t324;
t162 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t244;
t161 = rSges(6,1) * t213 - rSges(6,2) * t212 + rSges(6,3) * t246;
t160 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t244;
t159 = Icges(6,1) * t213 - Icges(6,4) * t212 + Icges(6,5) * t246;
t158 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t244;
t157 = Icges(6,4) * t213 - Icges(6,2) * t212 + Icges(6,6) * t246;
t156 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t244;
t155 = Icges(6,5) * t213 - Icges(6,6) * t212 + Icges(6,3) * t246;
t154 = rSges(7,1) * t177 + rSges(7,2) * t176 - rSges(7,3) * t214;
t153 = rSges(7,1) * t175 + rSges(7,2) * t174 + rSges(7,3) * t212;
t152 = Icges(7,1) * t177 + Icges(7,4) * t176 - Icges(7,5) * t214;
t151 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t212;
t150 = Icges(7,4) * t177 + Icges(7,2) * t176 - Icges(7,6) * t214;
t149 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t212;
t148 = Icges(7,5) * t177 + Icges(7,6) * t176 - Icges(7,3) * t214;
t147 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t212;
t146 = t226 * t282 + (-t194 - t236) * t293 + t321;
t145 = t195 * t293 + (-t226 - t276) * t283 + t315;
t144 = t194 * t283 + (-t195 - t237) * t282 + t320;
t143 = t225 * t282 + (-t193 + t339) * t293 + t317;
t142 = t192 * t293 + (-t225 + t337) * t283 + t314;
t141 = t193 * t283 + (-t192 + t338) * t282 + t316;
t140 = -t162 * t249 + t201 * t210 + t313;
t139 = t161 * t249 - t201 * t211 + t311;
t138 = -t161 * t210 + t162 * t211 + t312;
t137 = -t154 * t204 + t167 * t172 - t171 * t249 + t202 * t210 + t313;
t136 = t153 * t204 - t167 * t173 + t170 * t249 - t202 * t211 + t311;
t135 = -t153 * t172 + t154 * t173 - t170 * t210 + t171 * t211 + t312;
t1 = m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(3) * (t163 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + t210 * ((t155 * t244 + t157 * t214 + t159 * t215) * t211 + (t244 * t156 + t214 * t158 + t215 * t160) * t210 + (t198 * t244 + t199 * t214 + t200 * t215) * t249) / 0.2e1 + t211 * ((t246 * t155 - t212 * t157 + t213 * t159) * t211 + (t156 * t246 - t158 * t212 + t160 * t213) * t210 + (t198 * t246 - t199 * t212 + t200 * t213) * t249) / 0.2e1 + t204 * ((t147 * t250 + t149 * t206 + t151 * t207) * t173 + (t148 * t250 + t150 * t206 + t152 * t207) * t172 + (t250 * t164 + t206 * t165 + t207 * t166) * t204) / 0.2e1 + m(2) * (t248 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + t249 * ((t155 * t266 - t157 * t250 + t159 * t251) * t211 + (t156 * t266 - t158 * t250 + t160 * t251) * t210 + (t198 * t266 - t199 * t250 + t200 * t251) * t249) / 0.2e1 + m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t173 * ((t212 * t147 + t174 * t149 + t175 * t151) * t173 + (t148 * t212 + t150 * t174 + t152 * t175) * t172 + (t164 * t212 + t165 * t174 + t166 * t175) * t204) / 0.2e1 + t172 * ((-t147 * t214 + t149 * t176 + t151 * t177) * t173 + (-t214 * t148 + t176 * t150 + t177 * t152) * t172 + (-t164 * t214 + t165 * t176 + t166 * t177) * t204) / 0.2e1 + ((-t286 * t302 + t288 * t304 + Icges(1,4)) * V_base(5) + (-t287 * t302 + t289 * t304 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t286 * t304 + t288 * t302 + Icges(1,2)) * V_base(5) + (t287 * t304 + t289 * t302 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t243 * t359 + t244 * t360 + t262 * t268 + t263 * t269 - t343 * t354) * t293 + (t230 * t268 + t232 * t269 - t361 * t243 + t244 * t364 - t356 * t343) * t283 + (t229 * t268 + t231 * t269 - t243 * t362 + t244 * t363 - t343 * t355) * t282) * t282 / 0.2e1 + ((-t245 * t359 + t246 * t360 + t262 * t270 + t263 * t271 + t344 * t354) * t293 + (t230 * t270 + t232 * t271 - t361 * t245 + t246 * t364 + t356 * t344) * t283 + (t229 * t270 + t231 * t271 - t245 * t362 + t246 * t363 + t344 * t355) * t282) * t283 / 0.2e1 + (((t262 * t310 + t263 * t308) * t303 + t360 * t266 - t359 * t265 + t354 * t305) * t293 + ((t230 * t310 + t232 * t308) * t303 + t364 * t266 - t361 * t265 + t356 * t305) * t283 + ((t229 * t310 + t231 * t308) * t303 + t363 * t266 - t362 * t265 + t355 * t305) * t282) * t293 / 0.2e1 + ((Icges(2,5) * t304 - Icges(2,6) * t302 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t302 + Icges(2,6) * t304 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
