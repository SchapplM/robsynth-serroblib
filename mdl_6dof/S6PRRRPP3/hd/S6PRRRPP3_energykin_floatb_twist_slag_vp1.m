% Calculate kinetic energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:22
% EndTime: 2019-03-08 22:54:25
% DurationCPUTime: 3.49s
% Computational Cost: add. (2790->334), mult. (6605->472), div. (0->0), fcn. (8019->10), ass. (0->149)
t359 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t358 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t357 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t356 = -Icges(6,4) + Icges(5,5) + Icges(7,5);
t355 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t354 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t353 = rSges(7,1) + pkin(5);
t352 = rSges(7,3) + qJ(6);
t298 = sin(pkin(10));
t300 = cos(pkin(10));
t304 = cos(qJ(2));
t301 = cos(pkin(6));
t303 = sin(qJ(2));
t328 = t301 * t303;
t264 = t298 * t304 + t300 * t328;
t299 = sin(pkin(6));
t302 = sin(qJ(3));
t330 = t299 * t302;
t337 = cos(qJ(3));
t247 = t264 * t337 - t300 * t330;
t327 = t301 * t304;
t263 = t298 * t303 - t300 * t327;
t335 = sin(qJ(4));
t336 = cos(qJ(4));
t219 = t247 * t335 - t263 * t336;
t220 = t247 * t336 + t263 * t335;
t317 = t299 * t337;
t246 = t264 * t302 + t300 * t317;
t349 = t357 * t219 + t359 * t220 + t356 * t246;
t266 = -t298 * t328 + t300 * t304;
t249 = t266 * t337 + t298 * t330;
t265 = t298 * t327 + t300 * t303;
t221 = t249 * t335 - t265 * t336;
t222 = t249 * t336 + t265 * t335;
t248 = t266 * t302 - t298 * t317;
t348 = t357 * t221 + t359 * t222 + t356 * t248;
t347 = t354 * t219 + t357 * t220 + t355 * t246;
t346 = t354 * t221 + t357 * t222 + t355 * t248;
t345 = t355 * t219 + t356 * t220 + t358 * t246;
t344 = t355 * t221 + t356 * t222 + t358 * t248;
t271 = t301 * t302 + t303 * t317;
t329 = t299 * t304;
t250 = t271 * t335 + t329 * t336;
t251 = t271 * t336 - t329 * t335;
t270 = -t301 * t337 + t303 * t330;
t343 = t357 * t250 + t359 * t251 + t356 * t270;
t342 = t354 * t250 + t357 * t251 + t355 * t270;
t341 = t355 * t250 + t356 * t251 + t358 * t270;
t334 = pkin(7) * t301;
t333 = Icges(2,4) * t298;
t332 = t298 * t299;
t331 = t299 * t300;
t326 = rSges(7,2) * t219 + t352 * t220 + t353 * t246;
t325 = rSges(7,2) * t221 + t352 * t222 + t353 * t248;
t324 = rSges(7,2) * t250 + t352 * t251 + t353 * t270;
t323 = qJD(2) * t299;
t322 = V_base(5) * qJ(1) + V_base(1);
t318 = qJD(1) + V_base(3);
t279 = t298 * t323 + V_base(4);
t291 = qJD(2) * t301 + V_base(6);
t245 = qJD(3) * t265 + t279;
t278 = -t300 * t323 + V_base(5);
t244 = qJD(3) * t263 + t278;
t267 = -qJD(3) * t329 + t291;
t273 = pkin(1) * t298 - pkin(7) * t331;
t316 = -t273 * V_base(6) + V_base(5) * t334 + t322;
t274 = pkin(1) * t300 + pkin(7) * t332;
t315 = V_base(4) * t273 - t274 * V_base(5) + t318;
t314 = V_base(6) * t274 + V_base(2) + (-qJ(1) - t334) * V_base(4);
t238 = pkin(2) * t264 + pkin(8) * t263;
t272 = (pkin(2) * t303 - pkin(8) * t304) * t299;
t313 = -t238 * t291 + t278 * t272 + t316;
t239 = pkin(2) * t266 + pkin(8) * t265;
t312 = t279 * t238 - t239 * t278 + t315;
t311 = t291 * t239 - t272 * t279 + t314;
t212 = pkin(3) * t247 + pkin(9) * t246;
t240 = pkin(3) * t271 + pkin(9) * t270;
t310 = -t212 * t267 + t244 * t240 + t313;
t213 = pkin(3) * t249 + pkin(9) * t248;
t309 = t245 * t212 - t213 * t244 + t312;
t308 = t267 * t213 - t240 * t245 + t311;
t214 = pkin(4) * t251 + qJ(5) * t250;
t215 = qJD(4) * t246 + t244;
t307 = qJD(5) * t221 + t215 * t214 + t310;
t185 = pkin(4) * t220 + qJ(5) * t219;
t216 = qJD(4) * t248 + t245;
t306 = qJD(5) * t250 + t216 * t185 + t309;
t186 = pkin(4) * t222 + qJ(5) * t221;
t241 = qJD(4) * t270 + t267;
t305 = qJD(5) * t219 + t241 * t186 + t308;
t296 = Icges(2,4) * t300;
t287 = rSges(2,1) * t300 - rSges(2,2) * t298;
t286 = rSges(2,1) * t298 + rSges(2,2) * t300;
t285 = Icges(2,1) * t300 - t333;
t284 = Icges(2,1) * t298 + t296;
t283 = -Icges(2,2) * t298 + t296;
t282 = Icges(2,2) * t300 + t333;
t277 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t276 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t275 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t258 = t301 * rSges(3,3) + (rSges(3,1) * t303 + rSges(3,2) * t304) * t299;
t257 = Icges(3,5) * t301 + (Icges(3,1) * t303 + Icges(3,4) * t304) * t299;
t256 = Icges(3,6) * t301 + (Icges(3,4) * t303 + Icges(3,2) * t304) * t299;
t255 = Icges(3,3) * t301 + (Icges(3,5) * t303 + Icges(3,6) * t304) * t299;
t254 = V_base(5) * rSges(2,3) - t286 * V_base(6) + t322;
t253 = t287 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t242 = t286 * V_base(4) - t287 * V_base(5) + t318;
t237 = t271 * rSges(4,1) - t270 * rSges(4,2) - rSges(4,3) * t329;
t236 = Icges(4,1) * t271 - Icges(4,4) * t270 - Icges(4,5) * t329;
t235 = Icges(4,4) * t271 - Icges(4,2) * t270 - Icges(4,6) * t329;
t234 = Icges(4,5) * t271 - Icges(4,6) * t270 - Icges(4,3) * t329;
t233 = rSges(3,1) * t266 - rSges(3,2) * t265 + rSges(3,3) * t332;
t232 = rSges(3,1) * t264 - rSges(3,2) * t263 - rSges(3,3) * t331;
t231 = Icges(3,1) * t266 - Icges(3,4) * t265 + Icges(3,5) * t332;
t230 = Icges(3,1) * t264 - Icges(3,4) * t263 - Icges(3,5) * t331;
t229 = Icges(3,4) * t266 - Icges(3,2) * t265 + Icges(3,6) * t332;
t228 = Icges(3,4) * t264 - Icges(3,2) * t263 - Icges(3,6) * t331;
t227 = Icges(3,5) * t266 - Icges(3,6) * t265 + Icges(3,3) * t332;
t226 = Icges(3,5) * t264 - Icges(3,6) * t263 - Icges(3,3) * t331;
t210 = rSges(5,1) * t251 - rSges(5,2) * t250 + rSges(5,3) * t270;
t209 = rSges(6,1) * t270 - rSges(6,2) * t251 + rSges(6,3) * t250;
t197 = rSges(4,1) * t249 - rSges(4,2) * t248 + rSges(4,3) * t265;
t196 = rSges(4,1) * t247 - rSges(4,2) * t246 + rSges(4,3) * t263;
t195 = Icges(4,1) * t249 - Icges(4,4) * t248 + Icges(4,5) * t265;
t194 = Icges(4,1) * t247 - Icges(4,4) * t246 + Icges(4,5) * t263;
t193 = Icges(4,4) * t249 - Icges(4,2) * t248 + Icges(4,6) * t265;
t192 = Icges(4,4) * t247 - Icges(4,2) * t246 + Icges(4,6) * t263;
t191 = Icges(4,5) * t249 - Icges(4,6) * t248 + Icges(4,3) * t265;
t190 = Icges(4,5) * t247 - Icges(4,6) * t246 + Icges(4,3) * t263;
t184 = -t232 * t291 + t258 * t278 + t316;
t183 = t233 * t291 - t258 * t279 + t314;
t181 = rSges(5,1) * t222 - rSges(5,2) * t221 + rSges(5,3) * t248;
t180 = rSges(5,1) * t220 - rSges(5,2) * t219 + rSges(5,3) * t246;
t179 = rSges(6,1) * t248 - rSges(6,2) * t222 + rSges(6,3) * t221;
t177 = rSges(6,1) * t246 - rSges(6,2) * t220 + rSges(6,3) * t219;
t156 = t232 * t279 - t233 * t278 + t315;
t154 = -t196 * t267 + t237 * t244 + t313;
t153 = t197 * t267 - t237 * t245 + t311;
t152 = t196 * t245 - t197 * t244 + t312;
t151 = -t180 * t241 + t210 * t215 + t310;
t150 = t181 * t241 - t210 * t216 + t308;
t149 = t180 * t216 - t181 * t215 + t309;
t148 = t209 * t215 + (-t177 - t185) * t241 + t307;
t147 = t179 * t241 + (-t209 - t214) * t216 + t305;
t146 = t177 * t216 + (-t179 - t186) * t215 + t306;
t145 = qJD(6) * t222 + t324 * t215 + (-t185 - t326) * t241 + t307;
t144 = qJD(6) * t220 + t325 * t241 + (-t214 - t324) * t216 + t305;
t143 = qJD(6) * t251 + t326 * t216 + (-t186 - t325) * t215 + t306;
t1 = m(7) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(6) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(5) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(4) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(3) * (t156 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + t278 * ((-t227 * t331 - t229 * t263 + t231 * t264) * t279 + (-t226 * t331 - t228 * t263 + t230 * t264) * t278 + (-t255 * t331 - t256 * t263 + t257 * t264) * t291) / 0.2e1 + t279 * ((t227 * t332 - t229 * t265 + t231 * t266) * t279 + (t226 * t332 - t228 * t265 + t230 * t266) * t278 + (t255 * t332 - t256 * t265 + t257 * t266) * t291) / 0.2e1 + t267 * ((-t191 * t329 - t270 * t193 + t271 * t195) * t245 + (-t190 * t329 - t270 * t192 + t271 * t194) * t244 + (-t234 * t329 - t270 * t235 + t271 * t236) * t267) / 0.2e1 + t291 * ((t226 * t278 + t227 * t279 + t255 * t291) * t301 + ((t229 * t304 + t231 * t303) * t279 + (t228 * t304 + t230 * t303) * t278 + (t256 * t304 + t257 * t303) * t291) * t299) / 0.2e1 + m(2) * (t242 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + t244 * ((t191 * t263 - t193 * t246 + t195 * t247) * t245 + (t190 * t263 - t192 * t246 + t194 * t247) * t244 + (t234 * t263 - t235 * t246 + t236 * t247) * t267) / 0.2e1 + t245 * ((t191 * t265 - t193 * t248 + t195 * t249) * t245 + (t190 * t265 - t192 * t248 + t194 * t249) * t244 + (t234 * t265 - t235 * t248 + t236 * t249) * t267) / 0.2e1 + m(1) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + ((-t282 * t298 + t284 * t300 + Icges(1,4)) * V_base(5) + (-t283 * t298 + t285 * t300 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t282 * t300 + t284 * t298 + Icges(1,2)) * V_base(5) + (t283 * t300 + t285 * t298 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t219 * t342 + t220 * t343 + t246 * t341) * t241 + (t219 * t346 + t220 * t348 + t246 * t344) * t216 + (t347 * t219 + t349 * t220 + t345 * t246) * t215) * t215 / 0.2e1 + ((t221 * t342 + t222 * t343 + t248 * t341) * t241 + (t346 * t221 + t348 * t222 + t344 * t248) * t216 + (t221 * t347 + t222 * t349 + t248 * t345) * t215) * t216 / 0.2e1 + ((t250 * t342 + t251 * t343 + t270 * t341) * t241 + (t250 * t346 + t251 * t348 + t270 * t344) * t216 + (t250 * t347 + t251 * t349 + t270 * t345) * t215) * t241 / 0.2e1 + ((Icges(2,5) * t300 - Icges(2,6) * t298 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t298 + Icges(2,6) * t300 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
