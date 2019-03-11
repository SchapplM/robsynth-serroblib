% Calculate kinetic energy for
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:02
% EndTime: 2019-03-09 21:17:06
% DurationCPUTime: 3.96s
% Computational Cost: add. (3284->380), mult. (6709->532), div. (0->0), fcn. (8131->12), ass. (0->169)
t384 = Icges(6,1) + Icges(7,1);
t383 = -Icges(6,4) + Icges(7,5);
t382 = Icges(7,4) + Icges(6,5);
t381 = Icges(6,2) + Icges(7,3);
t380 = -Icges(6,6) + Icges(7,6);
t379 = Icges(7,2) + Icges(5,3) + Icges(6,3);
t378 = rSges(7,1) + pkin(5);
t377 = rSges(7,3) + qJ(6);
t319 = cos(pkin(6));
t324 = sin(qJ(1));
t326 = cos(qJ(2));
t352 = t324 * t326;
t323 = sin(qJ(2));
t327 = cos(qJ(1));
t353 = t323 * t327;
t286 = t319 * t353 + t352;
t322 = sin(qJ(3));
t318 = sin(pkin(6));
t355 = t318 * t327;
t364 = cos(qJ(3));
t265 = t286 * t364 - t322 * t355;
t351 = t326 * t327;
t354 = t323 * t324;
t285 = -t319 * t351 + t354;
t346 = qJ(4) + pkin(11);
t314 = sin(t346);
t340 = cos(t346);
t230 = t265 * t314 - t285 * t340;
t231 = t265 * t340 + t285 * t314;
t341 = t318 * t364;
t264 = t286 * t322 + t327 * t341;
t376 = t381 * t230 + t383 * t231 + t380 * t264;
t288 = -t319 * t354 + t351;
t357 = t318 * t324;
t267 = t288 * t364 + t322 * t357;
t287 = t319 * t352 + t353;
t232 = t267 * t314 - t287 * t340;
t233 = t267 * t340 + t287 * t314;
t266 = t288 * t322 - t324 * t341;
t375 = t381 * t232 + t383 * t233 + t380 * t266;
t374 = t383 * t230 + t384 * t231 + t382 * t264;
t373 = t383 * t232 + t384 * t233 + t382 * t266;
t284 = t319 * t322 + t323 * t341;
t356 = t318 * t326;
t256 = t284 * t314 + t340 * t356;
t257 = t284 * t340 - t314 * t356;
t283 = t318 * t322 * t323 - t319 * t364;
t372 = t381 * t256 + t383 * t257 + t380 * t283;
t371 = t383 * t256 + t384 * t257 + t382 * t283;
t321 = sin(qJ(4));
t325 = cos(qJ(4));
t234 = -t265 * t321 + t285 * t325;
t359 = t285 * t321;
t235 = t265 * t325 + t359;
t370 = Icges(5,5) * t235 + Icges(5,6) * t234 + t380 * t230 + t382 * t231 + t379 * t264;
t236 = -t267 * t321 + t287 * t325;
t358 = t287 * t321;
t237 = t267 * t325 + t358;
t369 = Icges(5,5) * t237 + Icges(5,6) * t236 + t380 * t232 + t382 * t233 + t379 * t266;
t260 = -t284 * t321 - t325 * t356;
t342 = t321 * t356;
t261 = t284 * t325 - t342;
t368 = Icges(5,5) * t261 + Icges(5,6) * t260 + t380 * t256 + t382 * t257 + t379 * t283;
t363 = pkin(8) * t319;
t362 = pkin(4) * t325;
t360 = Icges(2,4) * t324;
t350 = rSges(7,2) * t264 + t377 * t230 + t231 * t378;
t349 = rSges(7,2) * t266 + t377 * t232 + t233 * t378;
t348 = rSges(7,2) * t283 + t377 * t256 + t257 * t378;
t347 = qJD(2) * t318;
t345 = V_base(5) * pkin(7) + V_base(1);
t297 = t324 * t347 + V_base(4);
t315 = V_base(6) + qJD(1);
t263 = qJD(3) * t287 + t297;
t298 = qJD(2) * t319 + t315;
t296 = -t327 * t347 + V_base(5);
t291 = t324 * pkin(1) - pkin(8) * t355;
t339 = -t291 * t315 + V_base(5) * t363 + t345;
t292 = pkin(1) * t327 + pkin(8) * t357;
t338 = V_base(4) * t291 - t292 * V_base(5) + V_base(3);
t262 = qJD(3) * t285 + t296;
t281 = -qJD(3) * t356 + t298;
t337 = t315 * t292 + V_base(2) + (-pkin(7) - t363) * V_base(4);
t254 = pkin(2) * t286 + pkin(9) * t285;
t290 = (pkin(2) * t323 - pkin(9) * t326) * t318;
t336 = -t254 * t298 + t296 * t290 + t339;
t255 = pkin(2) * t288 + pkin(9) * t287;
t335 = t297 * t254 - t255 * t296 + t338;
t334 = t298 * t255 - t290 * t297 + t337;
t226 = pkin(3) * t265 + pkin(10) * t264;
t252 = pkin(3) * t284 + pkin(10) * t283;
t333 = -t226 * t281 + t262 * t252 + t336;
t227 = pkin(3) * t267 + pkin(10) * t266;
t332 = t263 * t226 - t227 * t262 + t335;
t210 = -pkin(4) * t342 + qJ(5) * t283 + t284 * t362;
t228 = qJD(4) * t264 + t262;
t331 = qJD(5) * t266 + t228 * t210 + t333;
t169 = pkin(4) * t359 + qJ(5) * t264 + t265 * t362;
t229 = qJD(4) * t266 + t263;
t330 = qJD(5) * t283 + t229 * t169 + t332;
t329 = t281 * t227 - t252 * t263 + t334;
t170 = pkin(4) * t358 + qJ(5) * t266 + t267 * t362;
t253 = qJD(4) * t283 + t281;
t328 = qJD(5) * t264 + t253 * t170 + t329;
t316 = Icges(2,4) * t327;
t306 = rSges(2,1) * t327 - t324 * rSges(2,2);
t305 = t324 * rSges(2,1) + rSges(2,2) * t327;
t304 = Icges(2,1) * t327 - t360;
t303 = Icges(2,1) * t324 + t316;
t302 = -Icges(2,2) * t324 + t316;
t301 = Icges(2,2) * t327 + t360;
t295 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t294 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t293 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t277 = rSges(3,3) * t319 + (rSges(3,1) * t323 + rSges(3,2) * t326) * t318;
t276 = Icges(3,5) * t319 + (Icges(3,1) * t323 + Icges(3,4) * t326) * t318;
t275 = Icges(3,6) * t319 + (Icges(3,4) * t323 + Icges(3,2) * t326) * t318;
t274 = Icges(3,3) * t319 + (Icges(3,5) * t323 + Icges(3,6) * t326) * t318;
t271 = V_base(5) * rSges(2,3) - t305 * t315 + t345;
t270 = t306 * t315 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t268 = t305 * V_base(4) - t306 * V_base(5) + V_base(3);
t251 = rSges(3,1) * t288 - rSges(3,2) * t287 + rSges(3,3) * t357;
t250 = t286 * rSges(3,1) - t285 * rSges(3,2) - rSges(3,3) * t355;
t249 = Icges(3,1) * t288 - Icges(3,4) * t287 + Icges(3,5) * t357;
t248 = Icges(3,1) * t286 - Icges(3,4) * t285 - Icges(3,5) * t355;
t247 = Icges(3,4) * t288 - Icges(3,2) * t287 + Icges(3,6) * t357;
t246 = Icges(3,4) * t286 - Icges(3,2) * t285 - Icges(3,6) * t355;
t245 = Icges(3,5) * t288 - Icges(3,6) * t287 + Icges(3,3) * t357;
t244 = Icges(3,5) * t286 - Icges(3,6) * t285 - Icges(3,3) * t355;
t243 = rSges(4,1) * t284 - rSges(4,2) * t283 - rSges(4,3) * t356;
t242 = Icges(4,1) * t284 - Icges(4,4) * t283 - Icges(4,5) * t356;
t241 = Icges(4,4) * t284 - Icges(4,2) * t283 - Icges(4,6) * t356;
t240 = Icges(4,5) * t284 - Icges(4,6) * t283 - Icges(4,3) * t356;
t223 = rSges(4,1) * t267 - rSges(4,2) * t266 + rSges(4,3) * t287;
t222 = rSges(4,1) * t265 - rSges(4,2) * t264 + rSges(4,3) * t285;
t220 = Icges(4,1) * t267 - Icges(4,4) * t266 + Icges(4,5) * t287;
t219 = Icges(4,1) * t265 - Icges(4,4) * t264 + Icges(4,5) * t285;
t218 = Icges(4,4) * t267 - Icges(4,2) * t266 + Icges(4,6) * t287;
t217 = Icges(4,4) * t265 - Icges(4,2) * t264 + Icges(4,6) * t285;
t216 = Icges(4,5) * t267 - Icges(4,6) * t266 + Icges(4,3) * t287;
t215 = Icges(4,5) * t265 - Icges(4,6) * t264 + Icges(4,3) * t285;
t214 = rSges(5,1) * t261 + rSges(5,2) * t260 + rSges(5,3) * t283;
t213 = Icges(5,1) * t261 + Icges(5,4) * t260 + Icges(5,5) * t283;
t212 = Icges(5,4) * t261 + Icges(5,2) * t260 + Icges(5,6) * t283;
t209 = rSges(6,1) * t257 - rSges(6,2) * t256 + rSges(6,3) * t283;
t198 = -t250 * t298 + t277 * t296 + t339;
t197 = t251 * t298 - t277 * t297 + t337;
t196 = rSges(5,1) * t237 + rSges(5,2) * t236 + rSges(5,3) * t266;
t195 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t264;
t194 = Icges(5,1) * t237 + Icges(5,4) * t236 + Icges(5,5) * t266;
t193 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t264;
t192 = Icges(5,4) * t237 + Icges(5,2) * t236 + Icges(5,6) * t266;
t191 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t264;
t188 = t250 * t297 - t251 * t296 + t338;
t187 = rSges(6,1) * t233 - rSges(6,2) * t232 + rSges(6,3) * t266;
t185 = rSges(6,1) * t231 - rSges(6,2) * t230 + rSges(6,3) * t264;
t166 = -t222 * t281 + t243 * t262 + t336;
t165 = t223 * t281 - t243 * t263 + t334;
t164 = t222 * t263 - t223 * t262 + t335;
t163 = -t195 * t253 + t214 * t228 + t333;
t162 = t196 * t253 - t214 * t229 + t329;
t161 = t195 * t229 - t196 * t228 + t332;
t160 = t209 * t228 + (-t169 - t185) * t253 + t331;
t159 = t187 * t253 + (-t209 - t210) * t229 + t328;
t158 = t185 * t229 + (-t170 - t187) * t228 + t330;
t157 = qJD(6) * t232 + t348 * t228 + (-t169 - t350) * t253 + t331;
t156 = qJD(6) * t230 + t349 * t253 + (-t210 - t348) * t229 + t328;
t155 = qJD(6) * t256 + t350 * t229 + (-t170 - t349) * t228 + t330;
t1 = m(7) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + t281 * ((-t216 * t356 - t218 * t283 + t220 * t284) * t263 + (-t215 * t356 - t217 * t283 + t219 * t284) * t262 + (-t240 * t356 - t241 * t283 + t242 * t284) * t281) / 0.2e1 + t297 * ((t245 * t357 - t247 * t287 + t249 * t288) * t297 + (t244 * t357 - t246 * t287 + t248 * t288) * t296 + (t274 * t357 - t275 * t287 + t276 * t288) * t298) / 0.2e1 + t296 * ((-t245 * t355 - t285 * t247 + t286 * t249) * t297 + (-t244 * t355 - t285 * t246 + t286 * t248) * t296 + (-t274 * t355 - t285 * t275 + t286 * t276) * t298) / 0.2e1 + t262 * ((t216 * t285 - t218 * t264 + t220 * t265) * t263 + (t215 * t285 - t217 * t264 + t219 * t265) * t262 + (t240 * t285 - t241 * t264 + t242 * t265) * t281) / 0.2e1 + m(2) * (t268 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(4) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(5) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + t298 * ((t244 * t296 + t245 * t297 + t274 * t298) * t319 + ((t247 * t326 + t249 * t323) * t297 + (t246 * t326 + t248 * t323) * t296 + (t275 * t326 + t276 * t323) * t298) * t318) / 0.2e1 + m(6) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(1) * (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) / 0.2e1 + t263 * ((t216 * t287 - t218 * t266 + t220 * t267) * t263 + (t215 * t287 - t217 * t266 + t219 * t267) * t262 + (t240 * t287 - t241 * t266 + t242 * t267) * t281) / 0.2e1 + m(3) * (t188 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + ((-t324 * t301 + t303 * t327 + Icges(1,4)) * V_base(5) + (-t324 * t302 + t304 * t327 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t301 * t327 + t324 * t303 + Icges(1,2)) * V_base(5) + (t302 * t327 + t324 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t212 * t234 + t213 * t235 + t230 * t372 + t231 * t371 + t264 * t368) * t253 + (t192 * t234 + t194 * t235 + t230 * t375 + t231 * t373 + t264 * t369) * t229 + (t191 * t234 + t193 * t235 + t376 * t230 + t374 * t231 + t370 * t264) * t228) * t228 / 0.2e1 + ((t212 * t236 + t213 * t237 + t232 * t372 + t233 * t371 + t266 * t368) * t253 + (t192 * t236 + t194 * t237 + t375 * t232 + t373 * t233 + t369 * t266) * t229 + (t191 * t236 + t193 * t237 + t232 * t376 + t374 * t233 + t370 * t266) * t228) * t229 / 0.2e1 + ((t212 * t260 + t213 * t261 + t372 * t256 + t371 * t257 + t368 * t283) * t253 + (t192 * t260 + t194 * t261 + t256 * t375 + t257 * t373 + t283 * t369) * t229 + (t191 * t260 + t193 * t261 + t256 * t376 + t374 * t257 + t370 * t283) * t228) * t253 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t324 + Icges(2,6) * t327) * V_base(5) + (Icges(2,5) * t327 - Icges(2,6) * t324) * V_base(4) + Icges(2,3) * t315 / 0.2e1) * t315;
T  = t1;
