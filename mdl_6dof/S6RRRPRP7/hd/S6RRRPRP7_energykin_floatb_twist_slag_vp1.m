% Calculate kinetic energy for
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:45
% EndTime: 2019-03-09 17:03:49
% DurationCPUTime: 4.40s
% Computational Cost: add. (3459->381), mult. (5929->534), div. (0->0), fcn. (7055->12), ass. (0->172)
t388 = Icges(6,1) + Icges(7,1);
t387 = -Icges(6,4) + Icges(7,5);
t386 = Icges(7,4) + Icges(6,5);
t385 = Icges(6,2) + Icges(7,3);
t384 = Icges(7,2) + Icges(6,3);
t383 = -Icges(6,6) + Icges(7,6);
t382 = Icges(4,3) + Icges(5,3);
t381 = rSges(7,1) + pkin(5);
t380 = rSges(7,3) + qJ(6);
t314 = cos(pkin(6));
t319 = sin(qJ(1));
t321 = cos(qJ(2));
t348 = t319 * t321;
t318 = sin(qJ(2));
t322 = cos(qJ(1));
t349 = t318 * t322;
t281 = t314 * t349 + t348;
t342 = qJ(3) + pkin(11);
t309 = sin(t342);
t336 = cos(t342);
t313 = sin(pkin(6));
t352 = t313 * t322;
t251 = t281 * t336 - t309 * t352;
t347 = t321 * t322;
t350 = t318 * t319;
t280 = -t314 * t347 + t350;
t316 = sin(qJ(5));
t361 = cos(qJ(5));
t221 = t251 * t316 - t280 * t361;
t222 = t251 * t361 + t280 * t316;
t335 = t313 * t336;
t250 = t281 * t309 + t322 * t335;
t379 = t385 * t221 + t387 * t222 + t383 * t250;
t283 = -t314 * t350 + t347;
t355 = t313 * t319;
t253 = t283 * t336 + t309 * t355;
t282 = t314 * t348 + t349;
t223 = t253 * t316 - t282 * t361;
t224 = t253 * t361 + t282 * t316;
t252 = t283 * t309 - t319 * t335;
t378 = t385 * t223 + t387 * t224 + t383 * t252;
t377 = t383 * t221 + t386 * t222 + t384 * t250;
t376 = t383 * t223 + t386 * t224 + t384 * t252;
t375 = t387 * t221 + t388 * t222 + t386 * t250;
t374 = t387 * t223 + t388 * t224 + t386 * t252;
t268 = t314 * t309 + t318 * t335;
t353 = t313 * t321;
t248 = t268 * t316 + t353 * t361;
t249 = t268 * t361 - t316 * t353;
t356 = t313 * t318;
t267 = t309 * t356 - t314 * t336;
t373 = t385 * t248 + t387 * t249 + t383 * t267;
t372 = t383 * t248 + t386 * t249 + t384 * t267;
t371 = t387 * t248 + t388 * t249 + t386 * t267;
t317 = sin(qJ(3));
t320 = cos(qJ(3));
t256 = -t281 * t317 - t320 * t352;
t337 = t317 * t352;
t257 = t281 * t320 - t337;
t370 = Icges(4,5) * t257 + Icges(5,5) * t251 + Icges(4,6) * t256 - Icges(5,6) * t250 + t382 * t280;
t354 = t313 * t320;
t258 = -t283 * t317 + t319 * t354;
t338 = t317 * t355;
t259 = t283 * t320 + t338;
t369 = Icges(4,5) * t259 + Icges(5,5) * t253 + Icges(4,6) * t258 - Icges(5,6) * t252 + t382 * t282;
t278 = t314 * t320 - t317 * t356;
t351 = t314 * t317;
t279 = t318 * t354 + t351;
t368 = Icges(4,5) * t279 + Icges(5,5) * t268 + Icges(4,6) * t278 - Icges(5,6) * t267 - t382 * t353;
t360 = pkin(8) * t314;
t359 = pkin(3) * t320;
t357 = Icges(2,4) * t319;
t346 = rSges(7,2) * t250 + t380 * t221 + t222 * t381;
t345 = rSges(7,2) * t252 + t380 * t223 + t224 * t381;
t344 = rSges(7,2) * t267 + t380 * t248 + t249 * t381;
t343 = qJD(2) * t313;
t341 = V_base(5) * pkin(7) + V_base(1);
t293 = t319 * t343 + V_base(4);
t310 = V_base(6) + qJD(1);
t255 = qJD(3) * t282 + t293;
t294 = qJD(2) * t314 + t310;
t292 = -t322 * t343 + V_base(5);
t286 = t319 * pkin(1) - pkin(8) * t352;
t334 = -t286 * t310 + V_base(5) * t360 + t341;
t287 = pkin(1) * t322 + pkin(8) * t355;
t333 = V_base(4) * t286 - t287 * V_base(5) + V_base(3);
t254 = qJD(3) * t280 + t292;
t276 = -qJD(3) * t353 + t294;
t332 = t310 * t287 + V_base(2) + (-pkin(7) - t360) * V_base(4);
t246 = t281 * pkin(2) + t280 * pkin(9);
t285 = (pkin(2) * t318 - pkin(9) * t321) * t313;
t331 = -t246 * t294 + t292 * t285 + t334;
t247 = pkin(2) * t283 + pkin(9) * t282;
t330 = t293 * t246 - t247 * t292 + t333;
t235 = pkin(3) * t351 + (-qJ(4) * t321 + t318 * t359) * t313;
t329 = qJD(4) * t282 + t254 * t235 + t331;
t328 = t294 * t247 - t285 * t293 + t332;
t205 = pkin(3) * t338 + qJ(4) * t282 + t283 * t359;
t327 = qJD(4) * t280 + t276 * t205 + t328;
t204 = -pkin(3) * t337 + qJ(4) * t280 + t281 * t359;
t326 = -qJD(4) * t353 + t255 * t204 + t330;
t217 = pkin(4) * t251 + pkin(10) * t250;
t231 = pkin(4) * t268 + pkin(10) * t267;
t325 = t254 * t231 + (-t204 - t217) * t276 + t329;
t218 = pkin(4) * t253 + pkin(10) * t252;
t324 = t276 * t218 + (-t231 - t235) * t255 + t327;
t323 = t255 * t217 + (-t205 - t218) * t254 + t326;
t311 = Icges(2,4) * t322;
t302 = rSges(2,1) * t322 - t319 * rSges(2,2);
t301 = t319 * rSges(2,1) + rSges(2,2) * t322;
t300 = Icges(2,1) * t322 - t357;
t299 = Icges(2,1) * t319 + t311;
t298 = -Icges(2,2) * t319 + t311;
t297 = Icges(2,2) * t322 + t357;
t290 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t289 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t288 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t269 = rSges(3,3) * t314 + (rSges(3,1) * t318 + rSges(3,2) * t321) * t313;
t266 = Icges(3,5) * t314 + (Icges(3,1) * t318 + Icges(3,4) * t321) * t313;
t265 = Icges(3,6) * t314 + (Icges(3,4) * t318 + Icges(3,2) * t321) * t313;
t264 = Icges(3,3) * t314 + (Icges(3,5) * t318 + Icges(3,6) * t321) * t313;
t263 = V_base(5) * rSges(2,3) - t301 * t310 + t341;
t262 = t302 * t310 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t260 = t301 * V_base(4) - t302 * V_base(5) + V_base(3);
t245 = qJD(5) * t267 + t276;
t244 = rSges(3,1) * t283 - rSges(3,2) * t282 + rSges(3,3) * t355;
t243 = t281 * rSges(3,1) - t280 * rSges(3,2) - rSges(3,3) * t352;
t242 = Icges(3,1) * t283 - Icges(3,4) * t282 + Icges(3,5) * t355;
t241 = Icges(3,1) * t281 - Icges(3,4) * t280 - Icges(3,5) * t352;
t240 = Icges(3,4) * t283 - Icges(3,2) * t282 + Icges(3,6) * t355;
t239 = Icges(3,4) * t281 - Icges(3,2) * t280 - Icges(3,6) * t352;
t238 = Icges(3,5) * t283 - Icges(3,6) * t282 + Icges(3,3) * t355;
t237 = Icges(3,5) * t281 - Icges(3,6) * t280 - Icges(3,3) * t352;
t236 = rSges(4,1) * t279 + rSges(4,2) * t278 - rSges(4,3) * t353;
t234 = Icges(4,1) * t279 + Icges(4,4) * t278 - Icges(4,5) * t353;
t233 = Icges(4,4) * t279 + Icges(4,2) * t278 - Icges(4,6) * t353;
t228 = rSges(5,1) * t268 - rSges(5,2) * t267 - rSges(5,3) * t353;
t227 = Icges(5,1) * t268 - Icges(5,4) * t267 - Icges(5,5) * t353;
t226 = Icges(5,4) * t268 - Icges(5,2) * t267 - Icges(5,6) * t353;
t220 = qJD(5) * t252 + t255;
t219 = qJD(5) * t250 + t254;
t213 = rSges(4,1) * t259 + rSges(4,2) * t258 + rSges(4,3) * t282;
t212 = rSges(4,1) * t257 + rSges(4,2) * t256 + rSges(4,3) * t280;
t211 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t282;
t210 = Icges(4,1) * t257 + Icges(4,4) * t256 + Icges(4,5) * t280;
t209 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t282;
t208 = Icges(4,4) * t257 + Icges(4,2) * t256 + Icges(4,6) * t280;
t203 = rSges(5,1) * t253 - rSges(5,2) * t252 + rSges(5,3) * t282;
t202 = rSges(5,1) * t251 - rSges(5,2) * t250 + rSges(5,3) * t280;
t200 = Icges(5,1) * t253 - Icges(5,4) * t252 + Icges(5,5) * t282;
t199 = Icges(5,1) * t251 - Icges(5,4) * t250 + Icges(5,5) * t280;
t198 = Icges(5,4) * t253 - Icges(5,2) * t252 + Icges(5,6) * t282;
t197 = Icges(5,4) * t251 - Icges(5,2) * t250 + Icges(5,6) * t280;
t194 = rSges(6,1) * t249 - rSges(6,2) * t248 + rSges(6,3) * t267;
t181 = -t243 * t294 + t269 * t292 + t334;
t180 = t244 * t294 - t269 * t293 + t332;
t179 = t243 * t293 - t244 * t292 + t333;
t178 = rSges(6,1) * t224 - rSges(6,2) * t223 + rSges(6,3) * t252;
t176 = rSges(6,1) * t222 - rSges(6,2) * t221 + rSges(6,3) * t250;
t162 = -t212 * t276 + t236 * t254 + t331;
t161 = t213 * t276 - t236 * t255 + t328;
t160 = t212 * t255 - t213 * t254 + t330;
t159 = t228 * t254 + (-t202 - t204) * t276 + t329;
t158 = t203 * t276 + (-t228 - t235) * t255 + t327;
t157 = t202 * t255 + (-t203 - t205) * t254 + t326;
t156 = -t176 * t245 + t194 * t219 + t325;
t155 = t178 * t245 - t194 * t220 + t324;
t154 = t176 * t220 - t178 * t219 + t323;
t153 = qJD(6) * t223 + t219 * t344 - t245 * t346 + t325;
t152 = qJD(6) * t221 - t220 * t344 + t245 * t345 + t324;
t151 = qJD(6) * t248 - t219 * t345 + t220 * t346 + t323;
t1 = m(1) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + m(2) * (t260 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t293 * ((t238 * t355 - t240 * t282 + t242 * t283) * t293 + (t237 * t355 - t239 * t282 + t241 * t283) * t292 + (t264 * t355 - t265 * t282 + t266 * t283) * t294) / 0.2e1 + t292 * ((-t238 * t352 - t280 * t240 + t281 * t242) * t293 + (-t237 * t352 - t280 * t239 + t281 * t241) * t292 + (-t264 * t352 - t280 * t265 + t281 * t266) * t294) / 0.2e1 + m(3) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + t294 * ((t237 * t292 + t238 * t293 + t264 * t294) * t314 + ((t240 * t321 + t242 * t318) * t293 + (t239 * t321 + t241 * t318) * t292 + (t265 * t321 + t266 * t318) * t294) * t313) / 0.2e1 + ((t221 * t373 + t222 * t371 + t250 * t372) * t245 + (t221 * t378 + t222 * t374 + t250 * t376) * t220 + (t379 * t221 + t375 * t222 + t377 * t250) * t219) * t219 / 0.2e1 + ((t223 * t373 + t224 * t371 + t252 * t372) * t245 + (t378 * t223 + t374 * t224 + t376 * t252) * t220 + (t223 * t379 + t375 * t224 + t377 * t252) * t219) * t220 / 0.2e1 + ((t373 * t248 + t371 * t249 + t372 * t267) * t245 + (t248 * t378 + t249 * t374 + t267 * t376) * t220 + (t248 * t379 + t375 * t249 + t377 * t267) * t219) * t245 / 0.2e1 + ((-t226 * t250 + t227 * t251 + t233 * t256 + t234 * t257 + t280 * t368) * t276 + (-t198 * t250 + t200 * t251 + t209 * t256 + t211 * t257 + t280 * t369) * t255 + (-t197 * t250 + t199 * t251 + t208 * t256 + t210 * t257 + t370 * t280) * t254) * t254 / 0.2e1 + ((-t226 * t252 + t227 * t253 + t233 * t258 + t234 * t259 + t282 * t368) * t276 + (-t198 * t252 + t200 * t253 + t209 * t258 + t211 * t259 + t369 * t282) * t255 + (-t197 * t252 + t199 * t253 + t208 * t258 + t210 * t259 + t282 * t370) * t254) * t255 / 0.2e1 + ((-t226 * t267 + t227 * t268 + t233 * t278 + t234 * t279 - t368 * t353) * t276 + (-t198 * t267 + t200 * t268 + t209 * t278 + t211 * t279 - t353 * t369) * t255 + (-t197 * t267 + t199 * t268 + t208 * t278 + t210 * t279 - t353 * t370) * t254) * t276 / 0.2e1 + ((-t319 * t297 + t299 * t322 + Icges(1,4)) * V_base(5) + (-t319 * t298 + t300 * t322 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t297 * t322 + t319 * t299 + Icges(1,2)) * V_base(5) + (t298 * t322 + t319 * t300 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t319 + Icges(2,6) * t322) * V_base(5) + (Icges(2,5) * t322 - Icges(2,6) * t319) * V_base(4) + Icges(2,3) * t310 / 0.2e1) * t310;
T  = t1;
