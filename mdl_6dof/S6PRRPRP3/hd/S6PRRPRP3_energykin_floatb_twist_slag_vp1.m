% Calculate kinetic energy for
% S6PRRPRP3
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:06
% EndTime: 2019-03-08 21:34:10
% DurationCPUTime: 4.26s
% Computational Cost: add. (3152->381), mult. (6547->530), div. (0->0), fcn. (7924->12), ass. (0->169)
t389 = Icges(6,1) + Icges(7,1);
t388 = -Icges(6,4) + Icges(7,5);
t387 = Icges(7,4) + Icges(6,5);
t386 = Icges(4,2) + Icges(5,3);
t385 = Icges(6,2) + Icges(7,3);
t384 = Icges(7,2) + Icges(6,3);
t383 = -Icges(6,6) + Icges(7,6);
t382 = rSges(7,1) + pkin(5);
t381 = rSges(7,3) + qJ(6);
t314 = sin(pkin(10));
t317 = cos(pkin(10));
t322 = cos(qJ(2));
t318 = cos(pkin(6));
t321 = sin(qJ(2));
t350 = t318 * t321;
t278 = t314 * t322 + t317 * t350;
t315 = sin(pkin(6));
t320 = sin(qJ(3));
t352 = t315 * t320;
t360 = cos(qJ(3));
t262 = t278 * t360 - t317 * t352;
t349 = t318 * t322;
t277 = t314 * t321 - t317 * t349;
t343 = pkin(11) + qJ(5);
t310 = sin(t343);
t335 = cos(t343);
t226 = t262 * t310 - t277 * t335;
t227 = t262 * t335 + t277 * t310;
t336 = t315 * t360;
t261 = t278 * t320 + t317 * t336;
t380 = t226 * t385 + t227 * t388 + t261 * t383;
t280 = -t314 * t350 + t317 * t322;
t264 = t280 * t360 + t314 * t352;
t279 = t314 * t349 + t317 * t321;
t228 = t264 * t310 - t279 * t335;
t229 = t264 * t335 + t279 * t310;
t263 = t280 * t320 - t314 * t336;
t379 = t228 * t385 + t229 * t388 + t263 * t383;
t378 = t226 * t383 + t227 * t387 + t261 * t384;
t377 = t228 * t383 + t229 * t387 + t263 * t384;
t376 = t388 * t226 + t227 * t389 + t387 * t261;
t375 = t388 * t228 + t229 * t389 + t387 * t263;
t313 = sin(pkin(11));
t316 = cos(pkin(11));
t230 = -t262 * t313 + t277 * t316;
t356 = t277 * t313;
t231 = t262 * t316 + t356;
t374 = -Icges(4,4) * t262 + Icges(5,5) * t231 - Icges(4,6) * t277 + Icges(5,6) * t230 + t261 * t386;
t232 = -t264 * t313 + t279 * t316;
t355 = t279 * t313;
t233 = t264 * t316 + t355;
t373 = -Icges(4,4) * t264 + Icges(5,5) * t233 - Icges(4,6) * t279 + Icges(5,6) * t232 + t263 * t386;
t285 = t318 * t320 + t321 * t336;
t351 = t315 * t322;
t252 = t285 * t310 + t335 * t351;
t253 = t285 * t335 - t310 * t351;
t284 = -t318 * t360 + t321 * t352;
t372 = t252 * t385 + t253 * t388 + t284 * t383;
t371 = t252 * t383 + t253 * t387 + t284 * t384;
t370 = t388 * t252 + t253 * t389 + t387 * t284;
t259 = -t285 * t313 - t316 * t351;
t337 = t313 * t351;
t260 = t285 * t316 - t337;
t369 = -Icges(4,4) * t285 + Icges(5,5) * t260 + Icges(4,6) * t351 + Icges(5,6) * t259 + t284 * t386;
t359 = pkin(7) * t318;
t358 = pkin(4) * t316;
t357 = Icges(2,4) * t314;
t354 = t314 * t315;
t353 = t315 * t317;
t347 = rSges(7,2) * t261 + t381 * t226 + t227 * t382;
t346 = rSges(7,2) * t263 + t381 * t228 + t229 * t382;
t345 = rSges(7,2) * t284 + t381 * t252 + t253 * t382;
t344 = qJD(2) * t315;
t342 = V_base(5) * qJ(1) + V_base(1);
t338 = qJD(1) + V_base(3);
t293 = t314 * t344 + V_base(4);
t304 = qJD(2) * t318 + V_base(6);
t258 = qJD(3) * t279 + t293;
t292 = -t317 * t344 + V_base(5);
t257 = qJD(3) * t277 + t292;
t281 = -qJD(3) * t351 + t304;
t287 = pkin(1) * t314 - pkin(7) * t353;
t334 = -t287 * V_base(6) + V_base(5) * t359 + t342;
t288 = pkin(1) * t317 + pkin(7) * t354;
t333 = V_base(4) * t287 - t288 * V_base(5) + t338;
t332 = V_base(6) * t288 + V_base(2) + (-qJ(1) - t359) * V_base(4);
t248 = pkin(2) * t278 + pkin(8) * t277;
t286 = (pkin(2) * t321 - pkin(8) * t322) * t315;
t331 = -t248 * t304 + t292 * t286 + t334;
t249 = pkin(2) * t280 + pkin(8) * t279;
t330 = t293 * t248 - t249 * t292 + t333;
t329 = t304 * t249 - t286 * t293 + t332;
t250 = t285 * pkin(3) + t284 * qJ(4);
t328 = qJD(4) * t263 + t257 * t250 + t331;
t222 = pkin(3) * t262 + qJ(4) * t261;
t327 = qJD(4) * t284 + t258 * t222 + t330;
t223 = pkin(3) * t264 + qJ(4) * t263;
t326 = qJD(4) * t261 + t281 * t223 + t329;
t165 = pkin(4) * t356 + pkin(9) * t261 + t262 * t358;
t206 = -pkin(4) * t337 + pkin(9) * t284 + t285 * t358;
t325 = t257 * t206 + (-t165 - t222) * t281 + t328;
t166 = pkin(4) * t355 + pkin(9) * t263 + t264 * t358;
t324 = t258 * t165 + (-t166 - t223) * t257 + t327;
t323 = t281 * t166 + (-t206 - t250) * t258 + t326;
t311 = Icges(2,4) * t317;
t301 = rSges(2,1) * t317 - rSges(2,2) * t314;
t300 = rSges(2,1) * t314 + rSges(2,2) * t317;
t299 = Icges(2,1) * t317 - t357;
t298 = Icges(2,1) * t314 + t311;
t297 = -Icges(2,2) * t314 + t311;
t296 = Icges(2,2) * t317 + t357;
t291 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t290 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t289 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t273 = t318 * rSges(3,3) + (rSges(3,1) * t321 + rSges(3,2) * t322) * t315;
t272 = Icges(3,5) * t318 + (Icges(3,1) * t321 + Icges(3,4) * t322) * t315;
t271 = Icges(3,6) * t318 + (Icges(3,4) * t321 + Icges(3,2) * t322) * t315;
t270 = Icges(3,3) * t318 + (Icges(3,5) * t321 + Icges(3,6) * t322) * t315;
t267 = V_base(5) * rSges(2,3) - t300 * V_base(6) + t342;
t266 = t301 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t256 = t300 * V_base(4) - t301 * V_base(5) + t338;
t251 = qJD(5) * t284 + t281;
t247 = t285 * rSges(4,1) - t284 * rSges(4,2) - rSges(4,3) * t351;
t246 = Icges(4,1) * t285 - Icges(4,4) * t284 - Icges(4,5) * t351;
t244 = Icges(4,5) * t285 - Icges(4,6) * t284 - Icges(4,3) * t351;
t243 = rSges(3,1) * t280 - rSges(3,2) * t279 + rSges(3,3) * t354;
t242 = rSges(3,1) * t278 - rSges(3,2) * t277 - rSges(3,3) * t353;
t241 = Icges(3,1) * t280 - Icges(3,4) * t279 + Icges(3,5) * t354;
t240 = Icges(3,1) * t278 - Icges(3,4) * t277 - Icges(3,5) * t353;
t239 = Icges(3,4) * t280 - Icges(3,2) * t279 + Icges(3,6) * t354;
t238 = Icges(3,4) * t278 - Icges(3,2) * t277 - Icges(3,6) * t353;
t237 = Icges(3,5) * t280 - Icges(3,6) * t279 + Icges(3,3) * t354;
t236 = Icges(3,5) * t278 - Icges(3,6) * t277 - Icges(3,3) * t353;
t225 = qJD(5) * t263 + t258;
t224 = qJD(5) * t261 + t257;
t218 = rSges(5,1) * t260 + rSges(5,2) * t259 + rSges(5,3) * t284;
t217 = rSges(4,1) * t264 - rSges(4,2) * t263 + rSges(4,3) * t279;
t216 = rSges(4,1) * t262 - rSges(4,2) * t261 + rSges(4,3) * t277;
t215 = Icges(5,1) * t260 + Icges(5,4) * t259 + Icges(5,5) * t284;
t214 = Icges(5,4) * t260 + Icges(5,2) * t259 + Icges(5,6) * t284;
t212 = Icges(4,1) * t264 - Icges(4,4) * t263 + Icges(4,5) * t279;
t211 = Icges(4,1) * t262 - Icges(4,4) * t261 + Icges(4,5) * t277;
t208 = Icges(4,5) * t264 - Icges(4,6) * t263 + Icges(4,3) * t279;
t207 = Icges(4,5) * t262 - Icges(4,6) * t261 + Icges(4,3) * t277;
t205 = rSges(6,1) * t253 - rSges(6,2) * t252 + rSges(6,3) * t284;
t196 = -t242 * t304 + t273 * t292 + t334;
t195 = t243 * t304 - t273 * t293 + t332;
t191 = rSges(5,1) * t233 + rSges(5,2) * t232 + rSges(5,3) * t263;
t190 = rSges(5,1) * t231 + rSges(5,2) * t230 + rSges(5,3) * t261;
t189 = Icges(5,1) * t233 + Icges(5,4) * t232 + Icges(5,5) * t263;
t188 = Icges(5,1) * t231 + Icges(5,4) * t230 + Icges(5,5) * t261;
t187 = Icges(5,4) * t233 + Icges(5,2) * t232 + Icges(5,6) * t263;
t186 = Icges(5,4) * t231 + Icges(5,2) * t230 + Icges(5,6) * t261;
t183 = t242 * t293 - t243 * t292 + t333;
t182 = rSges(6,1) * t229 - rSges(6,2) * t228 + rSges(6,3) * t263;
t180 = rSges(6,1) * t227 - rSges(6,2) * t226 + rSges(6,3) * t261;
t162 = -t216 * t281 + t247 * t257 + t331;
t161 = t217 * t281 - t247 * t258 + t329;
t160 = t216 * t258 - t217 * t257 + t330;
t159 = t218 * t257 + (-t190 - t222) * t281 + t328;
t158 = t191 * t281 + (-t218 - t250) * t258 + t326;
t157 = t190 * t258 + (-t191 - t223) * t257 + t327;
t156 = -t180 * t251 + t205 * t224 + t325;
t155 = t182 * t251 - t205 * t225 + t323;
t154 = t180 * t225 - t182 * t224 + t324;
t153 = qJD(6) * t228 + t224 * t345 - t251 * t347 + t325;
t152 = qJD(6) * t226 - t225 * t345 + t251 * t346 + t323;
t151 = qJD(6) * t252 - t224 * t346 + t225 * t347 + t324;
t1 = m(1) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + t304 * ((t236 * t292 + t237 * t293 + t270 * t304) * t318 + ((t239 * t322 + t241 * t321) * t293 + (t238 * t322 + t240 * t321) * t292 + (t271 * t322 + t272 * t321) * t304) * t315) / 0.2e1 + m(3) * (t183 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + t293 * ((t237 * t354 - t239 * t279 + t241 * t280) * t293 + (t236 * t354 - t238 * t279 + t240 * t280) * t292 + (t270 * t354 - t271 * t279 + t272 * t280) * t304) / 0.2e1 + t292 * ((-t237 * t353 - t239 * t277 + t241 * t278) * t293 + (-t236 * t353 - t238 * t277 + t240 * t278) * t292 + (-t270 * t353 - t271 * t277 + t272 * t278) * t304) / 0.2e1 + m(2) * (t256 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + ((t226 * t372 + t227 * t370 + t261 * t371) * t251 + (t226 * t379 + t227 * t375 + t261 * t377) * t225 + (t380 * t226 + t376 * t227 + t378 * t261) * t224) * t224 / 0.2e1 + ((t228 * t372 + t229 * t370 + t263 * t371) * t251 + (t379 * t228 + t375 * t229 + t377 * t263) * t225 + (t228 * t380 + t376 * t229 + t378 * t263) * t224) * t225 / 0.2e1 + ((t372 * t252 + t370 * t253 + t371 * t284) * t251 + (t252 * t379 + t253 * t375 + t284 * t377) * t225 + (t252 * t380 + t376 * t253 + t378 * t284) * t224) * t251 / 0.2e1 + ((t214 * t230 + t215 * t231 + t244 * t277 + t246 * t262 + t261 * t369) * t281 + (t187 * t230 + t189 * t231 + t208 * t277 + t212 * t262 + t261 * t373) * t258 + (t186 * t230 + t188 * t231 + t207 * t277 + t211 * t262 + t374 * t261) * t257) * t257 / 0.2e1 + ((t214 * t232 + t215 * t233 + t244 * t279 + t246 * t264 + t263 * t369) * t281 + (t187 * t232 + t189 * t233 + t208 * t279 + t212 * t264 + t373 * t263) * t258 + (t186 * t232 + t188 * t233 + t207 * t279 + t211 * t264 + t263 * t374) * t257) * t258 / 0.2e1 + ((t214 * t259 + t215 * t260 - t244 * t351 + t285 * t246 + t369 * t284) * t281 + (t187 * t259 + t189 * t260 - t208 * t351 + t285 * t212 + t284 * t373) * t258 + (t186 * t259 + t188 * t260 - t207 * t351 + t285 * t211 + t284 * t374) * t257) * t281 / 0.2e1 + ((-t296 * t314 + t298 * t317 + Icges(1,4)) * V_base(5) + (-t297 * t314 + t299 * t317 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t317 + t298 * t314 + Icges(1,2)) * V_base(5) + (t297 * t317 + t299 * t314 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t317 - Icges(2,6) * t314 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t314 + Icges(2,6) * t317 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
