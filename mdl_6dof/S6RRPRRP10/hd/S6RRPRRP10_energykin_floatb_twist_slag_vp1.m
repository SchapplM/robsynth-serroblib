% Calculate kinetic energy for
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:26
% EndTime: 2019-03-09 12:36:31
% DurationCPUTime: 4.35s
% Computational Cost: add. (3404->386), mult. (5819->542), div. (0->0), fcn. (6923->12), ass. (0->176)
t386 = Icges(6,1) + Icges(7,1);
t385 = -Icges(6,4) + Icges(7,5);
t384 = Icges(7,4) + Icges(6,5);
t383 = Icges(6,2) + Icges(7,3);
t382 = Icges(7,2) + Icges(6,3);
t381 = -Icges(6,6) + Icges(7,6);
t380 = rSges(7,1) + pkin(5);
t379 = rSges(7,3) + qJ(6);
t316 = cos(pkin(6));
t320 = sin(qJ(1));
t321 = cos(qJ(2));
t349 = t320 * t321;
t319 = sin(qJ(2));
t322 = cos(qJ(1));
t350 = t319 * t322;
t281 = t316 * t350 + t349;
t342 = pkin(11) + qJ(4);
t309 = sin(t342);
t336 = cos(t342);
t314 = sin(pkin(6));
t352 = t314 * t322;
t251 = t281 * t336 - t309 * t352;
t348 = t321 * t322;
t351 = t319 * t320;
t280 = -t316 * t348 + t351;
t318 = sin(qJ(5));
t360 = cos(qJ(5));
t221 = t251 * t318 - t280 * t360;
t222 = t251 * t360 + t280 * t318;
t335 = t314 * t336;
t250 = t281 * t309 + t322 * t335;
t378 = t383 * t221 + t385 * t222 + t381 * t250;
t283 = -t316 * t351 + t348;
t354 = t314 * t320;
t253 = t283 * t336 + t309 * t354;
t282 = t316 * t349 + t350;
t223 = t253 * t318 - t282 * t360;
t224 = t253 * t360 + t282 * t318;
t252 = t283 * t309 - t320 * t335;
t377 = t383 * t223 + t385 * t224 + t381 * t252;
t376 = t381 * t221 + t384 * t222 + t382 * t250;
t375 = t381 * t223 + t384 * t224 + t382 * t252;
t374 = t385 * t221 + t386 * t222 + t384 * t250;
t373 = t385 * t223 + t386 * t224 + t384 * t252;
t268 = t316 * t309 + t319 * t335;
t353 = t314 * t321;
t248 = t268 * t318 + t353 * t360;
t249 = t268 * t360 - t318 * t353;
t355 = t314 * t319;
t267 = t309 * t355 - t316 * t336;
t372 = t383 * t248 + t385 * t249 + t381 * t267;
t371 = t381 * t248 + t384 * t249 + t382 * t267;
t370 = t385 * t248 + t386 * t249 + t384 * t267;
t313 = sin(pkin(11));
t315 = cos(pkin(11));
t254 = -t281 * t313 - t315 * t352;
t337 = t313 * t352;
t255 = t281 * t315 - t337;
t206 = Icges(4,5) * t255 + Icges(4,6) * t254 + Icges(4,3) * t280;
t239 = Icges(3,4) * t281 - Icges(3,2) * t280 - Icges(3,6) * t352;
t369 = t206 - t239;
t256 = -t283 * t313 + t315 * t354;
t338 = t313 * t354;
t257 = t283 * t315 + t338;
t207 = Icges(4,5) * t257 + Icges(4,6) * t256 + Icges(4,3) * t282;
t240 = Icges(3,4) * t283 - Icges(3,2) * t282 + Icges(3,6) * t354;
t368 = t207 - t240;
t278 = -t313 * t355 + t315 * t316;
t356 = t313 * t316;
t279 = t315 * t355 + t356;
t231 = Icges(4,5) * t279 + Icges(4,6) * t278 - Icges(4,3) * t353;
t265 = Icges(3,6) * t316 + (Icges(3,4) * t319 + Icges(3,2) * t321) * t314;
t367 = t231 - t265;
t359 = pkin(8) * t316;
t358 = pkin(3) * t315;
t357 = Icges(2,4) * t320;
t346 = rSges(7,2) * t250 + t379 * t221 + t222 * t380;
t345 = rSges(7,2) * t252 + t379 * t223 + t224 * t380;
t344 = rSges(7,2) * t267 + t379 * t248 + t249 * t380;
t343 = qJD(2) * t314;
t341 = V_base(5) * pkin(7) + V_base(1);
t293 = t320 * t343 + V_base(4);
t310 = V_base(6) + qJD(1);
t259 = qJD(4) * t282 + t293;
t294 = qJD(2) * t316 + t310;
t292 = -t322 * t343 + V_base(5);
t286 = t320 * pkin(1) - pkin(8) * t352;
t334 = -t286 * t310 + V_base(5) * t359 + t341;
t287 = pkin(1) * t322 + pkin(8) * t354;
t333 = V_base(4) * t286 - t287 * V_base(5) + V_base(3);
t258 = qJD(4) * t280 + t292;
t276 = -qJD(4) * t353 + t294;
t284 = (pkin(2) * t319 - qJ(3) * t321) * t314;
t332 = qJD(3) * t282 + t292 * t284 + t334;
t331 = t310 * t287 + V_base(2) + (-pkin(7) - t359) * V_base(4);
t247 = pkin(2) * t283 + qJ(3) * t282;
t330 = qJD(3) * t280 + t294 * t247 + t331;
t246 = t281 * pkin(2) + t280 * qJ(3);
t329 = -qJD(3) * t353 + t293 * t246 + t333;
t204 = -pkin(3) * t337 + pkin(9) * t280 + t281 * t358;
t236 = pkin(3) * t356 + (-pkin(9) * t321 + t319 * t358) * t314;
t328 = t292 * t236 + (-t204 - t246) * t294 + t332;
t205 = pkin(3) * t338 + pkin(9) * t282 + t283 * t358;
t327 = t294 * t205 + (-t236 - t284) * t293 + t330;
t326 = t293 * t204 + (-t205 - t247) * t292 + t329;
t216 = pkin(4) * t251 + pkin(10) * t250;
t235 = pkin(4) * t268 + pkin(10) * t267;
t325 = -t216 * t276 + t258 * t235 + t328;
t217 = pkin(4) * t253 + pkin(10) * t252;
t324 = t276 * t217 - t235 * t259 + t327;
t323 = t259 * t216 - t217 * t258 + t326;
t311 = Icges(2,4) * t322;
t302 = rSges(2,1) * t322 - t320 * rSges(2,2);
t301 = t320 * rSges(2,1) + rSges(2,2) * t322;
t300 = Icges(2,1) * t322 - t357;
t299 = Icges(2,1) * t320 + t311;
t298 = -Icges(2,2) * t320 + t311;
t297 = Icges(2,2) * t322 + t357;
t290 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t289 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t288 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t269 = rSges(3,3) * t316 + (rSges(3,1) * t319 + rSges(3,2) * t321) * t314;
t266 = Icges(3,5) * t316 + (Icges(3,1) * t319 + Icges(3,4) * t321) * t314;
t264 = Icges(3,3) * t316 + (Icges(3,5) * t319 + Icges(3,6) * t321) * t314;
t263 = V_base(5) * rSges(2,3) - t301 * t310 + t341;
t262 = t302 * t310 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t260 = t301 * V_base(4) - t302 * V_base(5) + V_base(3);
t245 = qJD(5) * t267 + t276;
t244 = rSges(3,1) * t283 - rSges(3,2) * t282 + rSges(3,3) * t354;
t243 = t281 * rSges(3,1) - t280 * rSges(3,2) - rSges(3,3) * t352;
t242 = Icges(3,1) * t283 - Icges(3,4) * t282 + Icges(3,5) * t354;
t241 = Icges(3,1) * t281 - Icges(3,4) * t280 - Icges(3,5) * t352;
t238 = Icges(3,5) * t283 - Icges(3,6) * t282 + Icges(3,3) * t354;
t237 = Icges(3,5) * t281 - Icges(3,6) * t280 - Icges(3,3) * t352;
t234 = rSges(4,1) * t279 + rSges(4,2) * t278 - rSges(4,3) * t353;
t233 = Icges(4,1) * t279 + Icges(4,4) * t278 - Icges(4,5) * t353;
t232 = Icges(4,4) * t279 + Icges(4,2) * t278 - Icges(4,6) * t353;
t228 = rSges(5,1) * t268 - rSges(5,2) * t267 - rSges(5,3) * t353;
t227 = Icges(5,1) * t268 - Icges(5,4) * t267 - Icges(5,5) * t353;
t226 = Icges(5,4) * t268 - Icges(5,2) * t267 - Icges(5,6) * t353;
t225 = Icges(5,5) * t268 - Icges(5,6) * t267 - Icges(5,3) * t353;
t219 = qJD(5) * t252 + t259;
t218 = qJD(5) * t250 + t258;
t213 = rSges(4,1) * t257 + rSges(4,2) * t256 + rSges(4,3) * t282;
t212 = rSges(4,1) * t255 + rSges(4,2) * t254 + rSges(4,3) * t280;
t211 = Icges(4,1) * t257 + Icges(4,4) * t256 + Icges(4,5) * t282;
t210 = Icges(4,1) * t255 + Icges(4,4) * t254 + Icges(4,5) * t280;
t209 = Icges(4,4) * t257 + Icges(4,2) * t256 + Icges(4,6) * t282;
t208 = Icges(4,4) * t255 + Icges(4,2) * t254 + Icges(4,6) * t280;
t203 = rSges(5,1) * t253 - rSges(5,2) * t252 + rSges(5,3) * t282;
t202 = rSges(5,1) * t251 - rSges(5,2) * t250 + rSges(5,3) * t280;
t200 = Icges(5,1) * t253 - Icges(5,4) * t252 + Icges(5,5) * t282;
t199 = Icges(5,1) * t251 - Icges(5,4) * t250 + Icges(5,5) * t280;
t198 = Icges(5,4) * t253 - Icges(5,2) * t252 + Icges(5,6) * t282;
t197 = Icges(5,4) * t251 - Icges(5,2) * t250 + Icges(5,6) * t280;
t196 = Icges(5,5) * t253 - Icges(5,6) * t252 + Icges(5,3) * t282;
t195 = Icges(5,5) * t251 - Icges(5,6) * t250 + Icges(5,3) * t280;
t194 = rSges(6,1) * t249 - rSges(6,2) * t248 + rSges(6,3) * t267;
t181 = -t243 * t294 + t269 * t292 + t334;
t180 = t244 * t294 - t269 * t293 + t331;
t179 = t243 * t293 - t244 * t292 + t333;
t178 = rSges(6,1) * t224 - rSges(6,2) * t223 + rSges(6,3) * t252;
t176 = rSges(6,1) * t222 - rSges(6,2) * t221 + rSges(6,3) * t250;
t162 = t234 * t292 + (-t212 - t246) * t294 + t332;
t161 = t213 * t294 + (-t234 - t284) * t293 + t330;
t160 = t212 * t293 + (-t213 - t247) * t292 + t329;
t159 = -t202 * t276 + t228 * t258 + t328;
t158 = t203 * t276 - t228 * t259 + t327;
t157 = t202 * t259 - t203 * t258 + t326;
t156 = -t176 * t245 + t194 * t218 + t325;
t155 = t178 * t245 - t194 * t219 + t324;
t154 = t176 * t219 - t178 * t218 + t323;
t153 = qJD(6) * t223 + t218 * t344 - t245 * t346 + t325;
t152 = qJD(6) * t221 - t219 * t344 + t245 * t345 + t324;
t151 = qJD(6) * t248 - t218 * t345 + t219 * t346 + t323;
t1 = m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + t276 * ((-t196 * t353 - t198 * t267 + t200 * t268) * t259 + (-t195 * t353 - t197 * t267 + t199 * t268) * t258 + (-t225 * t353 - t226 * t267 + t227 * t268) * t276) / 0.2e1 + m(1) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + t259 * ((t196 * t282 - t198 * t252 + t200 * t253) * t259 + (t195 * t282 - t197 * t252 + t199 * t253) * t258 + (t225 * t282 - t226 * t252 + t227 * t253) * t276) / 0.2e1 + t258 * ((t196 * t280 - t198 * t250 + t200 * t251) * t259 + (t195 * t280 - t197 * t250 + t199 * t251) * t258 + (t225 * t280 - t226 * t250 + t227 * t251) * t276) / 0.2e1 + m(2) * (t260 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(3) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + ((t221 * t372 + t222 * t370 + t250 * t371) * t245 + (t221 * t377 + t222 * t373 + t250 * t375) * t219 + (t378 * t221 + t374 * t222 + t376 * t250) * t218) * t218 / 0.2e1 + ((t223 * t372 + t224 * t370 + t252 * t371) * t245 + (t377 * t223 + t373 * t224 + t375 * t252) * t219 + (t223 * t378 + t374 * t224 + t376 * t252) * t218) * t219 / 0.2e1 + ((t372 * t248 + t370 * t249 + t371 * t267) * t245 + (t248 * t377 + t249 * t373 + t267 * t375) * t219 + (t248 * t378 + t374 * t249 + t376 * t267) * t218) * t245 / 0.2e1 + ((t232 * t254 + t233 * t255 - t264 * t352 + t281 * t266 + t280 * t367) * t294 + (t209 * t254 + t211 * t255 - t238 * t352 + t281 * t242 + t280 * t368) * t293 + (t208 * t254 + t210 * t255 - t237 * t352 + t281 * t241 + t369 * t280) * t292) * t292 / 0.2e1 + ((t232 * t256 + t233 * t257 + t264 * t354 + t266 * t283 + t282 * t367) * t294 + (t209 * t256 + t211 * t257 + t238 * t354 + t242 * t283 + t368 * t282) * t293 + (t208 * t256 + t210 * t257 + t237 * t354 + t241 * t283 + t282 * t369) * t292) * t293 / 0.2e1 + ((t237 * t292 + t238 * t293 + t264 * t294) * t316 + ((t240 * t321 + t242 * t319) * t293 + (t239 * t321 + t241 * t319) * t292 + (t265 * t321 + t266 * t319) * t294) * t314 + (-t207 * t353 + t209 * t278 + t211 * t279) * t293 + (-t206 * t353 + t208 * t278 + t210 * t279) * t292 + (-t231 * t353 + t232 * t278 + t233 * t279) * t294) * t294 / 0.2e1 + ((-t320 * t297 + t299 * t322 + Icges(1,4)) * V_base(5) + (-t320 * t298 + t300 * t322 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t297 * t322 + t320 * t299 + Icges(1,2)) * V_base(5) + (t298 * t322 + t320 * t300 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t320 + Icges(2,6) * t322) * V_base(5) + (Icges(2,5) * t322 - Icges(2,6) * t320) * V_base(4) + Icges(2,3) * t310 / 0.2e1) * t310;
T  = t1;
