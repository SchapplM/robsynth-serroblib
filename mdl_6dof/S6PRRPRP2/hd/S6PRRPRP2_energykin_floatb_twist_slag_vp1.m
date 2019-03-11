% Calculate kinetic energy for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:17
% EndTime: 2019-03-08 21:28:22
% DurationCPUTime: 4.37s
% Computational Cost: add. (3399->381), mult. (5929->533), div. (0->0), fcn. (7055->12), ass. (0->171)
t390 = Icges(6,1) + Icges(7,1);
t389 = -Icges(6,4) + Icges(7,5);
t388 = Icges(7,4) + Icges(6,5);
t387 = Icges(6,2) + Icges(7,3);
t386 = Icges(7,2) + Icges(6,3);
t385 = -Icges(6,6) + Icges(7,6);
t384 = Icges(4,3) + Icges(5,3);
t383 = rSges(7,1) + pkin(5);
t382 = rSges(7,3) + qJ(6);
t312 = sin(pkin(10));
t314 = cos(pkin(10));
t321 = cos(qJ(2));
t315 = cos(pkin(6));
t319 = sin(qJ(2));
t349 = t315 * t319;
t277 = t312 * t321 + t314 * t349;
t343 = qJ(3) + pkin(11);
t309 = sin(t343);
t335 = cos(t343);
t313 = sin(pkin(6));
t355 = t313 * t314;
t249 = t277 * t335 - t309 * t355;
t348 = t315 * t321;
t276 = t312 * t319 - t314 * t348;
t317 = sin(qJ(5));
t361 = cos(qJ(5));
t221 = t249 * t317 - t276 * t361;
t222 = t249 * t361 + t276 * t317;
t334 = t313 * t335;
t248 = t277 * t309 + t314 * t334;
t380 = t387 * t221 + t389 * t222 + t385 * t248;
t279 = -t312 * t349 + t314 * t321;
t356 = t312 * t313;
t251 = t279 * t335 + t309 * t356;
t278 = t312 * t348 + t314 * t319;
t223 = t251 * t317 - t278 * t361;
t224 = t251 * t361 + t278 * t317;
t250 = t279 * t309 - t312 * t334;
t379 = t387 * t223 + t389 * t224 + t385 * t250;
t378 = t385 * t221 + t388 * t222 + t386 * t248;
t377 = t385 * t223 + t388 * t224 + t386 * t250;
t376 = t389 * t221 + t390 * t222 + t388 * t248;
t375 = t389 * t223 + t390 * t224 + t388 * t250;
t268 = t315 * t309 + t319 * t334;
t351 = t313 * t321;
t252 = t268 * t317 + t351 * t361;
t253 = t268 * t361 - t317 * t351;
t353 = t313 * t319;
t267 = t309 * t353 - t315 * t335;
t374 = t387 * t252 + t389 * t253 + t385 * t267;
t373 = t385 * t252 + t388 * t253 + t386 * t267;
t372 = t389 * t252 + t390 * t253 + t388 * t267;
t318 = sin(qJ(3));
t320 = cos(qJ(3));
t352 = t313 * t320;
t257 = -t277 * t318 - t314 * t352;
t354 = t313 * t318;
t336 = t314 * t354;
t258 = t277 * t320 - t336;
t371 = Icges(4,5) * t258 + Icges(5,5) * t249 + Icges(4,6) * t257 - Icges(5,6) * t248 + t384 * t276;
t259 = -t279 * t318 + t312 * t352;
t337 = t312 * t354;
t260 = t279 * t320 + t337;
t370 = Icges(4,5) * t260 + Icges(5,5) * t251 + Icges(4,6) * t259 - Icges(5,6) * t250 + t384 * t278;
t283 = t315 * t320 - t318 * t353;
t350 = t315 * t318;
t284 = t319 * t352 + t350;
t369 = Icges(4,5) * t284 + Icges(5,5) * t268 + Icges(4,6) * t283 - Icges(5,6) * t267 - t384 * t351;
t360 = pkin(7) * t315;
t359 = pkin(3) * t320;
t357 = Icges(2,4) * t312;
t347 = rSges(7,2) * t248 + t382 * t221 + t383 * t222;
t346 = rSges(7,2) * t250 + t382 * t223 + t383 * t224;
t345 = rSges(7,2) * t267 + t382 * t252 + t383 * t253;
t344 = qJD(2) * t313;
t342 = V_base(5) * qJ(1) + V_base(1);
t338 = qJD(1) + V_base(3);
t293 = t312 * t344 + V_base(4);
t304 = qJD(2) * t315 + V_base(6);
t256 = qJD(3) * t278 + t293;
t292 = -t314 * t344 + V_base(5);
t255 = qJD(3) * t276 + t292;
t280 = -qJD(3) * t351 + t304;
t286 = pkin(1) * t312 - pkin(7) * t355;
t333 = -t286 * V_base(6) + V_base(5) * t360 + t342;
t287 = pkin(1) * t314 + pkin(7) * t356;
t332 = V_base(4) * t286 - V_base(5) * t287 + t338;
t331 = V_base(6) * t287 + V_base(2) + (-qJ(1) - t360) * V_base(4);
t246 = pkin(2) * t277 + pkin(8) * t276;
t285 = (pkin(2) * t319 - pkin(8) * t321) * t313;
t330 = -t246 * t304 + t292 * t285 + t333;
t247 = pkin(2) * t279 + pkin(8) * t278;
t329 = t293 * t246 - t292 * t247 + t332;
t328 = t304 * t247 - t285 * t293 + t331;
t243 = pkin(3) * t350 + (-qJ(4) * t321 + t319 * t359) * t313;
t327 = qJD(4) * t278 + t255 * t243 + t330;
t204 = pkin(3) * t337 + qJ(4) * t278 + t279 * t359;
t326 = qJD(4) * t276 + t280 * t204 + t328;
t203 = -pkin(3) * t336 + qJ(4) * t276 + t277 * t359;
t325 = -qJD(4) * t351 + t256 * t203 + t329;
t216 = pkin(4) * t249 + pkin(9) * t248;
t237 = pkin(4) * t268 + pkin(9) * t267;
t324 = t255 * t237 + (-t203 - t216) * t280 + t327;
t217 = pkin(4) * t251 + pkin(9) * t250;
t323 = t280 * t217 + (-t237 - t243) * t256 + t326;
t322 = t256 * t216 + (-t204 - t217) * t255 + t325;
t310 = Icges(2,4) * t314;
t301 = rSges(2,1) * t314 - rSges(2,2) * t312;
t300 = rSges(2,1) * t312 + rSges(2,2) * t314;
t299 = Icges(2,1) * t314 - t357;
t298 = Icges(2,1) * t312 + t310;
t297 = -Icges(2,2) * t312 + t310;
t296 = Icges(2,2) * t314 + t357;
t291 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t290 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t289 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t269 = t315 * rSges(3,3) + (rSges(3,1) * t319 + rSges(3,2) * t321) * t313;
t266 = Icges(3,5) * t315 + (Icges(3,1) * t319 + Icges(3,4) * t321) * t313;
t265 = Icges(3,6) * t315 + (Icges(3,4) * t319 + Icges(3,2) * t321) * t313;
t264 = Icges(3,3) * t315 + (Icges(3,5) * t319 + Icges(3,6) * t321) * t313;
t263 = V_base(5) * rSges(2,3) - t300 * V_base(6) + t342;
t262 = t301 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t254 = t300 * V_base(4) - t301 * V_base(5) + t338;
t245 = qJD(5) * t267 + t280;
t244 = t284 * rSges(4,1) + t283 * rSges(4,2) - rSges(4,3) * t351;
t242 = Icges(4,1) * t284 + Icges(4,4) * t283 - Icges(4,5) * t351;
t241 = Icges(4,4) * t284 + Icges(4,2) * t283 - Icges(4,6) * t351;
t239 = rSges(3,1) * t279 - rSges(3,2) * t278 + rSges(3,3) * t356;
t238 = rSges(3,1) * t277 - rSges(3,2) * t276 - rSges(3,3) * t355;
t236 = Icges(3,1) * t279 - Icges(3,4) * t278 + Icges(3,5) * t356;
t235 = Icges(3,1) * t277 - Icges(3,4) * t276 - Icges(3,5) * t355;
t234 = Icges(3,4) * t279 - Icges(3,2) * t278 + Icges(3,6) * t356;
t233 = Icges(3,4) * t277 - Icges(3,2) * t276 - Icges(3,6) * t355;
t232 = Icges(3,5) * t279 - Icges(3,6) * t278 + Icges(3,3) * t356;
t231 = Icges(3,5) * t277 - Icges(3,6) * t276 - Icges(3,3) * t355;
t228 = t268 * rSges(5,1) - t267 * rSges(5,2) - rSges(5,3) * t351;
t227 = Icges(5,1) * t268 - Icges(5,4) * t267 - Icges(5,5) * t351;
t226 = Icges(5,4) * t268 - Icges(5,2) * t267 - Icges(5,6) * t351;
t220 = qJD(5) * t250 + t256;
t219 = qJD(5) * t248 + t255;
t213 = rSges(4,1) * t260 + rSges(4,2) * t259 + rSges(4,3) * t278;
t212 = rSges(4,1) * t258 + rSges(4,2) * t257 + rSges(4,3) * t276;
t211 = Icges(4,1) * t260 + Icges(4,4) * t259 + Icges(4,5) * t278;
t210 = Icges(4,1) * t258 + Icges(4,4) * t257 + Icges(4,5) * t276;
t209 = Icges(4,4) * t260 + Icges(4,2) * t259 + Icges(4,6) * t278;
t208 = Icges(4,4) * t258 + Icges(4,2) * t257 + Icges(4,6) * t276;
t202 = rSges(5,1) * t251 - rSges(5,2) * t250 + rSges(5,3) * t278;
t201 = rSges(5,1) * t249 - rSges(5,2) * t248 + rSges(5,3) * t276;
t200 = Icges(5,1) * t251 - Icges(5,4) * t250 + Icges(5,5) * t278;
t199 = Icges(5,1) * t249 - Icges(5,4) * t248 + Icges(5,5) * t276;
t198 = Icges(5,4) * t251 - Icges(5,2) * t250 + Icges(5,6) * t278;
t197 = Icges(5,4) * t249 - Icges(5,2) * t248 + Icges(5,6) * t276;
t194 = rSges(6,1) * t253 - rSges(6,2) * t252 + rSges(6,3) * t267;
t184 = -t238 * t304 + t269 * t292 + t333;
t183 = t239 * t304 - t269 * t293 + t331;
t179 = t238 * t293 - t239 * t292 + t332;
t178 = rSges(6,1) * t224 - rSges(6,2) * t223 + rSges(6,3) * t250;
t176 = rSges(6,1) * t222 - rSges(6,2) * t221 + rSges(6,3) * t248;
t162 = -t212 * t280 + t244 * t255 + t330;
t161 = t213 * t280 - t244 * t256 + t328;
t160 = t212 * t256 - t213 * t255 + t329;
t159 = t228 * t255 + (-t201 - t203) * t280 + t327;
t158 = t202 * t280 + (-t228 - t243) * t256 + t326;
t157 = t256 * t201 + (-t202 - t204) * t255 + t325;
t156 = -t176 * t245 + t194 * t219 + t324;
t155 = t178 * t245 - t194 * t220 + t323;
t154 = t220 * t176 - t219 * t178 + t322;
t153 = qJD(6) * t223 + t219 * t345 - t245 * t347 + t324;
t152 = qJD(6) * t221 - t220 * t345 + t245 * t346 + t323;
t151 = qJD(6) * t252 - t219 * t346 + t220 * t347 + t322;
t1 = t292 * ((-t232 * t355 - t234 * t276 + t236 * t277) * t293 + (-t231 * t355 - t233 * t276 + t235 * t277) * t292 + (-t264 * t355 - t265 * t276 + t266 * t277) * t304) / 0.2e1 + t293 * ((t232 * t356 - t234 * t278 + t236 * t279) * t293 + (t231 * t356 - t233 * t278 + t235 * t279) * t292 + (t264 * t356 - t265 * t278 + t266 * t279) * t304) / 0.2e1 + m(2) * (t254 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(1) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(3) * (t179 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + t304 * ((t231 * t292 + t232 * t293 + t264 * t304) * t315 + ((t234 * t321 + t236 * t319) * t293 + (t233 * t321 + t235 * t319) * t292 + (t265 * t321 + t266 * t319) * t304) * t313) / 0.2e1 + ((t221 * t374 + t222 * t372 + t248 * t373) * t245 + (t221 * t379 + t222 * t375 + t248 * t377) * t220 + (t380 * t221 + t376 * t222 + t378 * t248) * t219) * t219 / 0.2e1 + ((t223 * t374 + t224 * t372 + t250 * t373) * t245 + (t379 * t223 + t375 * t224 + t377 * t250) * t220 + (t223 * t380 + t224 * t376 + t250 * t378) * t219) * t220 / 0.2e1 + ((t374 * t252 + t372 * t253 + t373 * t267) * t245 + (t252 * t379 + t253 * t375 + t267 * t377) * t220 + (t252 * t380 + t253 * t376 + t267 * t378) * t219) * t245 / 0.2e1 + ((-t226 * t248 + t227 * t249 + t241 * t257 + t242 * t258 + t276 * t369) * t280 + (-t198 * t248 + t200 * t249 + t209 * t257 + t211 * t258 + t276 * t370) * t256 + (-t197 * t248 + t199 * t249 + t208 * t257 + t210 * t258 + t371 * t276) * t255) * t255 / 0.2e1 + ((-t226 * t250 + t227 * t251 + t241 * t259 + t242 * t260 + t278 * t369) * t280 + (-t198 * t250 + t200 * t251 + t209 * t259 + t211 * t260 + t370 * t278) * t256 + (-t197 * t250 + t199 * t251 + t208 * t259 + t210 * t260 + t278 * t371) * t255) * t256 / 0.2e1 + ((-t267 * t226 + t268 * t227 + t283 * t241 + t284 * t242 - t369 * t351) * t280 + (-t267 * t198 + t268 * t200 + t283 * t209 + t284 * t211 - t351 * t370) * t256 + (-t267 * t197 + t268 * t199 + t283 * t208 + t284 * t210 - t351 * t371) * t255) * t280 / 0.2e1 + ((-t296 * t312 + t298 * t314 + Icges(1,4)) * V_base(5) + (-t297 * t312 + t299 * t314 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t314 + t298 * t312 + Icges(1,2)) * V_base(5) + (t297 * t314 + t299 * t312 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t314 - Icges(2,6) * t312 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t312 + Icges(2,6) * t314 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
