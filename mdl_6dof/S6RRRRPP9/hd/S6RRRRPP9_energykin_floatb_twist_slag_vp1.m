% Calculate kinetic energy for
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:42
% EndTime: 2019-03-09 21:41:45
% DurationCPUTime: 3.60s
% Computational Cost: add. (2850->334), mult. (6605->474), div. (0->0), fcn. (8019->10), ass. (0->150)
t357 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t356 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t355 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t354 = -Icges(6,4) + Icges(5,5) + Icges(7,5);
t353 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t352 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t351 = rSges(7,1) + pkin(5);
t350 = rSges(7,3) + qJ(6);
t300 = cos(pkin(6));
t303 = sin(qJ(1));
t304 = cos(qJ(2));
t327 = t303 * t304;
t302 = sin(qJ(2));
t305 = cos(qJ(1));
t328 = t302 * t305;
t268 = t300 * t328 + t327;
t301 = sin(qJ(3));
t299 = sin(pkin(6));
t330 = t299 * t305;
t337 = cos(qJ(3));
t248 = t268 * t337 - t301 * t330;
t326 = t304 * t305;
t329 = t302 * t303;
t267 = -t300 * t326 + t329;
t335 = sin(qJ(4));
t336 = cos(qJ(4));
t219 = t248 * t335 - t267 * t336;
t220 = t248 * t336 + t267 * t335;
t318 = t299 * t337;
t247 = t268 * t301 + t305 * t318;
t349 = t355 * t219 + t357 * t220 + t354 * t247;
t270 = -t300 * t329 + t326;
t332 = t299 * t303;
t250 = t270 * t337 + t301 * t332;
t269 = t300 * t327 + t328;
t221 = t250 * t335 - t269 * t336;
t222 = t250 * t336 + t269 * t335;
t249 = t270 * t301 - t303 * t318;
t348 = t355 * t221 + t357 * t222 + t354 * t249;
t347 = t352 * t219 + t355 * t220 + t353 * t247;
t346 = t352 * t221 + t355 * t222 + t353 * t249;
t345 = t353 * t219 + t354 * t220 + t356 * t247;
t344 = t353 * t221 + t354 * t222 + t356 * t249;
t266 = t300 * t301 + t302 * t318;
t331 = t299 * t304;
t243 = t266 * t335 + t331 * t336;
t244 = t266 * t336 - t331 * t335;
t265 = t299 * t301 * t302 - t300 * t337;
t343 = t355 * t243 + t357 * t244 + t354 * t265;
t342 = t352 * t243 + t355 * t244 + t353 * t265;
t341 = t353 * t243 + t354 * t244 + t356 * t265;
t334 = pkin(8) * t300;
t333 = Icges(2,4) * t303;
t325 = rSges(7,2) * t219 + t350 * t220 + t247 * t351;
t324 = rSges(7,2) * t221 + t350 * t222 + t249 * t351;
t323 = rSges(7,2) * t243 + t350 * t244 + t265 * t351;
t322 = qJD(2) * t299;
t321 = V_base(5) * pkin(7) + V_base(1);
t279 = t303 * t322 + V_base(4);
t296 = V_base(6) + qJD(1);
t246 = qJD(3) * t269 + t279;
t280 = qJD(2) * t300 + t296;
t278 = -t305 * t322 + V_base(5);
t273 = t303 * pkin(1) - pkin(8) * t330;
t317 = -t273 * t296 + V_base(5) * t334 + t321;
t274 = pkin(1) * t305 + pkin(8) * t332;
t316 = V_base(4) * t273 - t274 * V_base(5) + V_base(3);
t245 = qJD(3) * t267 + t278;
t263 = -qJD(3) * t331 + t280;
t315 = t296 * t274 + V_base(2) + (-pkin(7) - t334) * V_base(4);
t240 = pkin(2) * t268 + pkin(9) * t267;
t272 = (pkin(2) * t302 - pkin(9) * t304) * t299;
t314 = -t240 * t280 + t278 * t272 + t317;
t241 = pkin(2) * t270 + pkin(9) * t269;
t313 = t279 * t240 - t241 * t278 + t316;
t312 = t280 * t241 - t272 * t279 + t315;
t213 = pkin(3) * t248 + pkin(10) * t247;
t238 = pkin(3) * t266 + pkin(10) * t265;
t311 = -t213 * t263 + t245 * t238 + t314;
t214 = pkin(3) * t250 + pkin(10) * t249;
t310 = t246 * t213 - t214 * t245 + t313;
t212 = pkin(4) * t244 + qJ(5) * t243;
t215 = qJD(4) * t247 + t245;
t309 = qJD(5) * t221 + t215 * t212 + t311;
t185 = pkin(4) * t220 + qJ(5) * t219;
t216 = qJD(4) * t249 + t246;
t308 = qJD(5) * t243 + t216 * t185 + t310;
t307 = t263 * t214 - t238 * t246 + t312;
t186 = pkin(4) * t222 + qJ(5) * t221;
t239 = qJD(4) * t265 + t263;
t306 = qJD(5) * t219 + t239 * t186 + t307;
t297 = Icges(2,4) * t305;
t288 = rSges(2,1) * t305 - t303 * rSges(2,2);
t287 = t303 * rSges(2,1) + rSges(2,2) * t305;
t286 = Icges(2,1) * t305 - t333;
t285 = Icges(2,1) * t303 + t297;
t284 = -Icges(2,2) * t303 + t297;
t283 = Icges(2,2) * t305 + t333;
t277 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t276 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t275 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t258 = rSges(3,3) * t300 + (rSges(3,1) * t302 + rSges(3,2) * t304) * t299;
t257 = Icges(3,5) * t300 + (Icges(3,1) * t302 + Icges(3,4) * t304) * t299;
t256 = Icges(3,6) * t300 + (Icges(3,4) * t302 + Icges(3,2) * t304) * t299;
t255 = Icges(3,3) * t300 + (Icges(3,5) * t302 + Icges(3,6) * t304) * t299;
t254 = V_base(5) * rSges(2,3) - t287 * t296 + t321;
t253 = t288 * t296 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t251 = t287 * V_base(4) - t288 * V_base(5) + V_base(3);
t237 = rSges(3,1) * t270 - rSges(3,2) * t269 + rSges(3,3) * t332;
t236 = t268 * rSges(3,1) - t267 * rSges(3,2) - rSges(3,3) * t330;
t235 = Icges(3,1) * t270 - Icges(3,4) * t269 + Icges(3,5) * t332;
t234 = Icges(3,1) * t268 - Icges(3,4) * t267 - Icges(3,5) * t330;
t233 = Icges(3,4) * t270 - Icges(3,2) * t269 + Icges(3,6) * t332;
t232 = Icges(3,4) * t268 - Icges(3,2) * t267 - Icges(3,6) * t330;
t231 = Icges(3,5) * t270 - Icges(3,6) * t269 + Icges(3,3) * t332;
t230 = Icges(3,5) * t268 - Icges(3,6) * t267 - Icges(3,3) * t330;
t229 = rSges(4,1) * t266 - rSges(4,2) * t265 - rSges(4,3) * t331;
t228 = Icges(4,1) * t266 - Icges(4,4) * t265 - Icges(4,5) * t331;
t227 = Icges(4,4) * t266 - Icges(4,2) * t265 - Icges(4,6) * t331;
t226 = Icges(4,5) * t266 - Icges(4,6) * t265 - Icges(4,3) * t331;
t210 = rSges(4,1) * t250 - rSges(4,2) * t249 + rSges(4,3) * t269;
t209 = rSges(4,1) * t248 - rSges(4,2) * t247 + rSges(4,3) * t267;
t207 = Icges(4,1) * t250 - Icges(4,4) * t249 + Icges(4,5) * t269;
t206 = Icges(4,1) * t248 - Icges(4,4) * t247 + Icges(4,5) * t267;
t205 = Icges(4,4) * t250 - Icges(4,2) * t249 + Icges(4,6) * t269;
t204 = Icges(4,4) * t248 - Icges(4,2) * t247 + Icges(4,6) * t267;
t203 = Icges(4,5) * t250 - Icges(4,6) * t249 + Icges(4,3) * t269;
t202 = Icges(4,5) * t248 - Icges(4,6) * t247 + Icges(4,3) * t267;
t201 = rSges(5,1) * t244 - rSges(5,2) * t243 + rSges(5,3) * t265;
t200 = rSges(6,1) * t265 - rSges(6,2) * t244 + rSges(6,3) * t243;
t184 = -t236 * t280 + t258 * t278 + t317;
t183 = t237 * t280 - t258 * t279 + t315;
t181 = rSges(5,1) * t222 - rSges(5,2) * t221 + rSges(5,3) * t249;
t180 = rSges(5,1) * t220 - rSges(5,2) * t219 + rSges(5,3) * t247;
t179 = rSges(6,1) * t249 - rSges(6,2) * t222 + rSges(6,3) * t221;
t177 = rSges(6,1) * t247 - rSges(6,2) * t220 + rSges(6,3) * t219;
t157 = t236 * t279 - t237 * t278 + t316;
t154 = -t209 * t263 + t229 * t245 + t314;
t153 = t210 * t263 - t229 * t246 + t312;
t152 = t209 * t246 - t210 * t245 + t313;
t151 = -t180 * t239 + t201 * t215 + t311;
t150 = t181 * t239 - t201 * t216 + t307;
t149 = t180 * t216 - t181 * t215 + t310;
t148 = t200 * t215 + (-t177 - t185) * t239 + t309;
t147 = t179 * t239 + (-t200 - t212) * t216 + t306;
t146 = t177 * t216 + (-t179 - t186) * t215 + t308;
t145 = qJD(6) * t222 + t323 * t215 + (-t185 - t325) * t239 + t309;
t144 = qJD(6) * t220 + t324 * t239 + (-t212 - t323) * t216 + t306;
t143 = qJD(6) * t244 + t325 * t216 + (-t186 - t324) * t215 + t308;
t1 = m(1) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + t245 * ((t203 * t267 - t205 * t247 + t207 * t248) * t246 + (t202 * t267 - t204 * t247 + t206 * t248) * t245 + (t226 * t267 - t227 * t247 + t228 * t248) * t263) / 0.2e1 + t246 * ((t203 * t269 - t205 * t249 + t207 * t250) * t246 + (t202 * t269 - t204 * t249 + t206 * t250) * t245 + (t226 * t269 - t227 * t249 + t228 * t250) * t263) / 0.2e1 + m(2) * (t251 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(3) * (t157 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + t280 * ((t230 * t278 + t231 * t279 + t255 * t280) * t300 + ((t233 * t304 + t235 * t302) * t279 + (t232 * t304 + t234 * t302) * t278 + (t256 * t304 + t257 * t302) * t280) * t299) / 0.2e1 + m(4) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(5) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(6) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(7) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + t278 * ((-t231 * t330 - t267 * t233 + t268 * t235) * t279 + (-t230 * t330 - t267 * t232 + t268 * t234) * t278 + (-t255 * t330 - t267 * t256 + t268 * t257) * t280) / 0.2e1 + t263 * ((-t203 * t331 - t205 * t265 + t207 * t266) * t246 + (-t202 * t331 - t204 * t265 + t206 * t266) * t245 + (-t226 * t331 - t227 * t265 + t228 * t266) * t263) / 0.2e1 + t279 * ((t231 * t332 - t233 * t269 + t235 * t270) * t279 + (t230 * t332 - t232 * t269 + t234 * t270) * t278 + (t255 * t332 - t256 * t269 + t257 * t270) * t280) / 0.2e1 + ((-t303 * t283 + t285 * t305 + Icges(1,4)) * V_base(5) + (-t303 * t284 + t286 * t305 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t283 * t305 + t303 * t285 + Icges(1,2)) * V_base(5) + (t284 * t305 + t303 * t286 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t219 * t342 + t220 * t343 + t247 * t341) * t239 + (t219 * t346 + t220 * t348 + t247 * t344) * t216 + (t347 * t219 + t349 * t220 + t345 * t247) * t215) * t215 / 0.2e1 + ((t221 * t342 + t222 * t343 + t249 * t341) * t239 + (t346 * t221 + t348 * t222 + t344 * t249) * t216 + (t347 * t221 + t222 * t349 + t345 * t249) * t215) * t216 / 0.2e1 + ((t342 * t243 + t343 * t244 + t341 * t265) * t239 + (t243 * t346 + t244 * t348 + t265 * t344) * t216 + (t347 * t243 + t244 * t349 + t345 * t265) * t215) * t239 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t303 + Icges(2,6) * t305) * V_base(5) + (Icges(2,5) * t305 - Icges(2,6) * t303) * V_base(4) + Icges(2,3) * t296 / 0.2e1) * t296;
T  = t1;
