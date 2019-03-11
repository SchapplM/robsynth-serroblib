% Calculate kinetic energy for
% S6RRRPRP10
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:22
% EndTime: 2019-03-09 17:29:27
% DurationCPUTime: 4.34s
% Computational Cost: add. (3212->381), mult. (6547->532), div. (0->0), fcn. (7924->12), ass. (0->170)
t391 = Icges(6,1) + Icges(7,1);
t390 = -Icges(6,4) + Icges(7,5);
t389 = Icges(7,4) + Icges(6,5);
t388 = Icges(4,2) + Icges(5,3);
t387 = Icges(6,2) + Icges(7,3);
t386 = Icges(7,2) + Icges(6,3);
t385 = -Icges(6,6) + Icges(7,6);
t384 = rSges(7,1) + pkin(5);
t383 = rSges(7,3) + qJ(6);
t321 = cos(pkin(6));
t325 = sin(qJ(1));
t326 = cos(qJ(2));
t353 = t325 * t326;
t324 = sin(qJ(2));
t327 = cos(qJ(1));
t354 = t324 * t327;
t286 = t321 * t354 + t353;
t323 = sin(qJ(3));
t319 = sin(pkin(6));
t356 = t319 * t327;
t364 = cos(qJ(3));
t265 = t286 * t364 - t323 * t356;
t352 = t326 * t327;
t355 = t324 * t325;
t285 = -t321 * t352 + t355;
t346 = pkin(11) + qJ(5);
t314 = sin(t346);
t340 = cos(t346);
t230 = t265 * t314 - t285 * t340;
t231 = t265 * t340 + t285 * t314;
t341 = t319 * t364;
t264 = t286 * t323 + t327 * t341;
t382 = t387 * t230 + t390 * t231 + t385 * t264;
t288 = -t321 * t355 + t352;
t358 = t319 * t325;
t267 = t288 * t364 + t323 * t358;
t287 = t321 * t353 + t354;
t232 = t267 * t314 - t287 * t340;
t233 = t267 * t340 + t287 * t314;
t266 = t288 * t323 - t325 * t341;
t381 = t387 * t232 + t390 * t233 + t385 * t266;
t380 = t385 * t230 + t389 * t231 + t386 * t264;
t379 = t385 * t232 + t389 * t233 + t386 * t266;
t378 = t390 * t230 + t391 * t231 + t389 * t264;
t377 = t390 * t232 + t391 * t233 + t389 * t266;
t318 = sin(pkin(11));
t320 = cos(pkin(11));
t234 = -t265 * t318 + t285 * t320;
t360 = t285 * t318;
t235 = t265 * t320 + t360;
t376 = -Icges(4,4) * t265 + Icges(5,5) * t235 - Icges(4,6) * t285 + Icges(5,6) * t234 + t388 * t264;
t236 = -t267 * t318 + t287 * t320;
t359 = t287 * t318;
t237 = t267 * t320 + t359;
t375 = -Icges(4,4) * t267 + Icges(5,5) * t237 - Icges(4,6) * t287 + Icges(5,6) * t236 + t388 * t266;
t284 = t321 * t323 + t324 * t341;
t357 = t319 * t326;
t256 = t284 * t314 + t340 * t357;
t257 = t284 * t340 - t314 * t357;
t283 = t319 * t323 * t324 - t321 * t364;
t374 = t387 * t256 + t390 * t257 + t385 * t283;
t373 = t385 * t256 + t389 * t257 + t386 * t283;
t372 = t390 * t256 + t391 * t257 + t389 * t283;
t260 = -t284 * t318 - t320 * t357;
t342 = t318 * t357;
t261 = t284 * t320 - t342;
t371 = -Icges(4,4) * t284 + Icges(5,5) * t261 + Icges(4,6) * t357 + Icges(5,6) * t260 + t388 * t283;
t363 = pkin(8) * t321;
t362 = pkin(4) * t320;
t361 = Icges(2,4) * t325;
t350 = rSges(7,2) * t264 + t383 * t230 + t231 * t384;
t349 = rSges(7,2) * t266 + t383 * t232 + t233 * t384;
t348 = rSges(7,2) * t283 + t383 * t256 + t257 * t384;
t347 = qJD(2) * t319;
t345 = V_base(5) * pkin(7) + V_base(1);
t297 = t325 * t347 + V_base(4);
t315 = V_base(6) + qJD(1);
t263 = qJD(3) * t287 + t297;
t298 = qJD(2) * t321 + t315;
t296 = -t327 * t347 + V_base(5);
t291 = t325 * pkin(1) - pkin(8) * t356;
t339 = -t291 * t315 + V_base(5) * t363 + t345;
t292 = pkin(1) * t327 + pkin(8) * t358;
t338 = V_base(4) * t291 - t292 * V_base(5) + V_base(3);
t262 = qJD(3) * t285 + t296;
t281 = -qJD(3) * t357 + t298;
t337 = t315 * t292 + V_base(2) + (-pkin(7) - t363) * V_base(4);
t254 = pkin(2) * t286 + pkin(9) * t285;
t290 = (pkin(2) * t324 - pkin(9) * t326) * t319;
t336 = -t254 * t298 + t296 * t290 + t339;
t255 = pkin(2) * t288 + pkin(9) * t287;
t335 = t297 * t254 - t255 * t296 + t338;
t252 = pkin(3) * t284 + qJ(4) * t283;
t334 = qJD(4) * t266 + t262 * t252 + t336;
t226 = pkin(3) * t265 + qJ(4) * t264;
t333 = qJD(4) * t283 + t263 * t226 + t335;
t332 = t298 * t255 - t290 * t297 + t337;
t227 = pkin(3) * t267 + qJ(4) * t266;
t331 = qJD(4) * t264 + t281 * t227 + t332;
t169 = pkin(4) * t360 + pkin(10) * t264 + t265 * t362;
t210 = -pkin(4) * t342 + pkin(10) * t283 + t284 * t362;
t330 = t262 * t210 + (-t169 - t226) * t281 + t334;
t170 = pkin(4) * t359 + pkin(10) * t266 + t267 * t362;
t329 = t263 * t169 + (-t170 - t227) * t262 + t333;
t328 = t281 * t170 + (-t210 - t252) * t263 + t331;
t316 = Icges(2,4) * t327;
t306 = rSges(2,1) * t327 - t325 * rSges(2,2);
t305 = t325 * rSges(2,1) + rSges(2,2) * t327;
t304 = Icges(2,1) * t327 - t361;
t303 = Icges(2,1) * t325 + t316;
t302 = -Icges(2,2) * t325 + t316;
t301 = Icges(2,2) * t327 + t361;
t295 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t294 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t293 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t277 = rSges(3,3) * t321 + (rSges(3,1) * t324 + rSges(3,2) * t326) * t319;
t276 = Icges(3,5) * t321 + (Icges(3,1) * t324 + Icges(3,4) * t326) * t319;
t275 = Icges(3,6) * t321 + (Icges(3,4) * t324 + Icges(3,2) * t326) * t319;
t274 = Icges(3,3) * t321 + (Icges(3,5) * t324 + Icges(3,6) * t326) * t319;
t271 = V_base(5) * rSges(2,3) - t305 * t315 + t345;
t270 = t306 * t315 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t268 = t305 * V_base(4) - t306 * V_base(5) + V_base(3);
t253 = qJD(5) * t283 + t281;
t251 = rSges(3,1) * t288 - rSges(3,2) * t287 + rSges(3,3) * t358;
t250 = t286 * rSges(3,1) - t285 * rSges(3,2) - rSges(3,3) * t356;
t249 = Icges(3,1) * t288 - Icges(3,4) * t287 + Icges(3,5) * t358;
t248 = Icges(3,1) * t286 - Icges(3,4) * t285 - Icges(3,5) * t356;
t247 = Icges(3,4) * t288 - Icges(3,2) * t287 + Icges(3,6) * t358;
t246 = Icges(3,4) * t286 - Icges(3,2) * t285 - Icges(3,6) * t356;
t245 = Icges(3,5) * t288 - Icges(3,6) * t287 + Icges(3,3) * t358;
t244 = Icges(3,5) * t286 - Icges(3,6) * t285 - Icges(3,3) * t356;
t243 = rSges(4,1) * t284 - rSges(4,2) * t283 - rSges(4,3) * t357;
t242 = Icges(4,1) * t284 - Icges(4,4) * t283 - Icges(4,5) * t357;
t240 = Icges(4,5) * t284 - Icges(4,6) * t283 - Icges(4,3) * t357;
t229 = qJD(5) * t266 + t263;
t228 = qJD(5) * t264 + t262;
t223 = rSges(4,1) * t267 - rSges(4,2) * t266 + rSges(4,3) * t287;
t222 = rSges(4,1) * t265 - rSges(4,2) * t264 + rSges(4,3) * t285;
t221 = Icges(4,1) * t267 - Icges(4,4) * t266 + Icges(4,5) * t287;
t220 = Icges(4,1) * t265 - Icges(4,4) * t264 + Icges(4,5) * t285;
t217 = Icges(4,5) * t267 - Icges(4,6) * t266 + Icges(4,3) * t287;
t216 = Icges(4,5) * t265 - Icges(4,6) * t264 + Icges(4,3) * t285;
t214 = rSges(5,1) * t261 + rSges(5,2) * t260 + rSges(5,3) * t283;
t213 = Icges(5,1) * t261 + Icges(5,4) * t260 + Icges(5,5) * t283;
t212 = Icges(5,4) * t261 + Icges(5,2) * t260 + Icges(5,6) * t283;
t209 = rSges(6,1) * t257 - rSges(6,2) * t256 + rSges(6,3) * t283;
t198 = -t250 * t298 + t277 * t296 + t339;
t197 = t251 * t298 - t277 * t297 + t337;
t195 = rSges(5,1) * t237 + rSges(5,2) * t236 + rSges(5,3) * t266;
t194 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t264;
t193 = Icges(5,1) * t237 + Icges(5,4) * t236 + Icges(5,5) * t266;
t192 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t264;
t191 = Icges(5,4) * t237 + Icges(5,2) * t236 + Icges(5,6) * t266;
t190 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t264;
t187 = t250 * t297 - t251 * t296 + t338;
t186 = rSges(6,1) * t233 - rSges(6,2) * t232 + rSges(6,3) * t266;
t184 = rSges(6,1) * t231 - rSges(6,2) * t230 + rSges(6,3) * t264;
t166 = -t222 * t281 + t243 * t262 + t336;
t165 = t223 * t281 - t243 * t263 + t332;
t164 = t222 * t263 - t223 * t262 + t335;
t163 = t214 * t262 + (-t194 - t226) * t281 + t334;
t162 = t195 * t281 + (-t214 - t252) * t263 + t331;
t161 = t194 * t263 + (-t195 - t227) * t262 + t333;
t160 = -t184 * t253 + t209 * t228 + t330;
t159 = t186 * t253 - t209 * t229 + t328;
t158 = t184 * t229 - t186 * t228 + t329;
t157 = qJD(6) * t232 + t228 * t348 - t253 * t350 + t330;
t156 = qJD(6) * t230 - t229 * t348 + t253 * t349 + t328;
t155 = qJD(6) * t256 - t228 * t349 + t229 * t350 + t329;
t1 = m(4) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(5) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(6) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(7) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + t297 * ((t245 * t358 - t247 * t287 + t249 * t288) * t297 + (t244 * t358 - t246 * t287 + t248 * t288) * t296 + (t274 * t358 - t275 * t287 + t276 * t288) * t298) / 0.2e1 + t296 * ((-t245 * t356 - t285 * t247 + t286 * t249) * t297 + (-t244 * t356 - t285 * t246 + t286 * t248) * t296 + (-t274 * t356 - t285 * t275 + t286 * t276) * t298) / 0.2e1 + t298 * ((t244 * t296 + t245 * t297 + t274 * t298) * t321 + ((t247 * t326 + t249 * t324) * t297 + (t246 * t326 + t248 * t324) * t296 + (t275 * t326 + t276 * t324) * t298) * t319) / 0.2e1 + m(1) * (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) / 0.2e1 + m(2) * (t268 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(3) * (t187 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + ((t230 * t374 + t231 * t372 + t264 * t373) * t253 + (t230 * t381 + t231 * t377 + t264 * t379) * t229 + (t382 * t230 + t378 * t231 + t380 * t264) * t228) * t228 / 0.2e1 + ((t232 * t374 + t233 * t372 + t266 * t373) * t253 + (t381 * t232 + t377 * t233 + t379 * t266) * t229 + (t232 * t382 + t378 * t233 + t380 * t266) * t228) * t229 / 0.2e1 + ((t374 * t256 + t372 * t257 + t373 * t283) * t253 + (t256 * t381 + t257 * t377 + t283 * t379) * t229 + (t256 * t382 + t378 * t257 + t380 * t283) * t228) * t253 / 0.2e1 + ((t212 * t234 + t213 * t235 + t240 * t285 + t242 * t265 + t264 * t371) * t281 + (t191 * t234 + t193 * t235 + t217 * t285 + t221 * t265 + t264 * t375) * t263 + (t190 * t234 + t192 * t235 + t216 * t285 + t220 * t265 + t376 * t264) * t262) * t262 / 0.2e1 + ((t212 * t236 + t213 * t237 + t240 * t287 + t242 * t267 + t266 * t371) * t281 + (t191 * t236 + t193 * t237 + t217 * t287 + t221 * t267 + t375 * t266) * t263 + (t190 * t236 + t192 * t237 + t216 * t287 + t220 * t267 + t266 * t376) * t262) * t263 / 0.2e1 + ((t212 * t260 + t213 * t261 - t240 * t357 + t242 * t284 + t371 * t283) * t281 + (t191 * t260 + t193 * t261 - t217 * t357 + t221 * t284 + t283 * t375) * t263 + (t190 * t260 + t192 * t261 - t216 * t357 + t220 * t284 + t283 * t376) * t262) * t281 / 0.2e1 + ((-t325 * t301 + t303 * t327 + Icges(1,4)) * V_base(5) + (-t325 * t302 + t304 * t327 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t301 * t327 + t325 * t303 + Icges(1,2)) * V_base(5) + (t302 * t327 + t325 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t325 + Icges(2,6) * t327) * V_base(5) + (Icges(2,5) * t327 - Icges(2,6) * t325) * V_base(4) + Icges(2,3) * t315 / 0.2e1) * t315;
T  = t1;
