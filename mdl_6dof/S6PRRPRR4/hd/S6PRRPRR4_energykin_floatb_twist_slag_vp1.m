% Calculate kinetic energy for
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:10
% EndTime: 2019-03-08 22:10:16
% DurationCPUTime: 6.18s
% Computational Cost: add. (2915->392), mult. (6948->565), div. (0->0), fcn. (8499->12), ass. (0->170)
t376 = Icges(4,1) + Icges(5,1);
t375 = -Icges(4,4) + Icges(5,5);
t374 = Icges(5,4) + Icges(4,5);
t373 = Icges(4,2) + Icges(5,3);
t372 = Icges(5,2) + Icges(4,3);
t371 = Icges(4,6) - Icges(5,6);
t318 = sin(pkin(11));
t320 = cos(pkin(11));
t326 = cos(qJ(2));
t324 = sin(qJ(2));
t354 = cos(pkin(6));
t340 = t324 * t354;
t284 = t318 * t326 + t320 * t340;
t323 = sin(qJ(3));
t319 = sin(pkin(6));
t356 = cos(qJ(3));
t342 = t319 * t356;
t264 = t284 * t323 + t320 * t342;
t350 = t319 * t323;
t265 = t284 * t356 - t320 * t350;
t339 = t326 * t354;
t283 = t318 * t324 - t320 * t339;
t368 = t373 * t264 + t375 * t265 - t371 * t283;
t286 = -t318 * t340 + t320 * t326;
t266 = t286 * t323 - t318 * t342;
t267 = t286 * t356 + t318 * t350;
t285 = t318 * t339 + t320 * t324;
t367 = t373 * t266 + t375 * t267 - t371 * t285;
t366 = -t371 * t264 + t374 * t265 + t372 * t283;
t365 = -t371 * t266 + t374 * t267 + t372 * t285;
t364 = t375 * t264 + t376 * t265 + t374 * t283;
t363 = t375 * t266 + t376 * t267 + t374 * t285;
t290 = t324 * t350 - t354 * t356;
t291 = t323 * t354 + t324 * t342;
t349 = t319 * t326;
t362 = t373 * t290 + t375 * t291 + t371 * t349;
t361 = -t371 * t290 + t374 * t291 - t372 * t349;
t360 = t375 * t290 + t376 * t291 - t374 * t349;
t355 = cos(qJ(5));
t353 = Icges(2,4) * t318;
t352 = t318 * t319;
t351 = t319 * t320;
t348 = qJD(2) * t319;
t347 = V_base(5) * qJ(1) + V_base(1);
t343 = qJD(1) + V_base(3);
t341 = t354 * pkin(7);
t299 = t318 * t348 + V_base(4);
t310 = qJD(2) * t354 + V_base(6);
t263 = qJD(3) * t285 + t299;
t234 = -qJD(5) * t285 + t263;
t298 = -t320 * t348 + V_base(5);
t262 = qJD(3) * t283 + t298;
t287 = -qJD(3) * t349 + t310;
t293 = pkin(1) * t318 - pkin(7) * t351;
t338 = -t293 * V_base(6) + V_base(5) * t341 + t347;
t294 = pkin(1) * t320 + pkin(7) * t352;
t337 = V_base(4) * t293 - t294 * V_base(5) + t343;
t233 = -qJD(5) * t283 + t262;
t272 = qJD(5) * t349 + t287;
t336 = V_base(6) * t294 + V_base(2) + (-t341 - qJ(1)) * V_base(4);
t254 = pkin(2) * t284 + pkin(8) * t283;
t292 = (pkin(2) * t324 - pkin(8) * t326) * t319;
t335 = -t254 * t310 + t298 * t292 + t338;
t255 = pkin(2) * t286 + pkin(8) * t285;
t334 = t299 * t254 - t255 * t298 + t337;
t256 = pkin(3) * t291 + qJ(4) * t290;
t333 = qJD(4) * t266 + t262 * t256 + t335;
t225 = pkin(3) * t265 + qJ(4) * t264;
t332 = qJD(4) * t290 + t263 * t225 + t334;
t331 = t310 * t255 - t299 * t292 + t336;
t226 = pkin(3) * t267 + qJ(4) * t266;
t330 = qJD(4) * t264 + t287 * t226 + t331;
t230 = pkin(4) * t265 - pkin(9) * t283;
t269 = t291 * pkin(4) + pkin(9) * t349;
t329 = t262 * t269 + (-t225 - t230) * t287 + t333;
t231 = pkin(4) * t267 - pkin(9) * t285;
t328 = t263 * t230 + (-t226 - t231) * t262 + t332;
t327 = t287 * t231 + (-t256 - t269) * t263 + t330;
t325 = cos(qJ(6));
t322 = sin(qJ(5));
t321 = sin(qJ(6));
t316 = Icges(2,4) * t320;
t307 = rSges(2,1) * t320 - rSges(2,2) * t318;
t306 = rSges(2,1) * t318 + rSges(2,2) * t320;
t305 = Icges(2,1) * t320 - t353;
t304 = Icges(2,1) * t318 + t316;
t303 = -Icges(2,2) * t318 + t316;
t302 = Icges(2,2) * t320 + t353;
t297 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t296 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t295 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t276 = t354 * rSges(3,3) + (rSges(3,1) * t324 + rSges(3,2) * t326) * t319;
t275 = Icges(3,5) * t354 + (Icges(3,1) * t324 + Icges(3,4) * t326) * t319;
t274 = Icges(3,6) * t354 + (Icges(3,4) * t324 + Icges(3,2) * t326) * t319;
t273 = Icges(3,3) * t354 + (Icges(3,5) * t324 + Icges(3,6) * t326) * t319;
t271 = V_base(5) * rSges(2,3) - t306 * V_base(6) + t347;
t270 = t307 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t261 = t306 * V_base(4) - t307 * V_base(5) + t343;
t253 = t290 * t322 + t291 * t355;
t252 = -t290 * t355 + t291 * t322;
t251 = t291 * rSges(4,1) - t290 * rSges(4,2) - rSges(4,3) * t349;
t250 = t291 * rSges(5,1) - rSges(5,2) * t349 + t290 * rSges(5,3);
t243 = rSges(3,1) * t286 - rSges(3,2) * t285 + rSges(3,3) * t352;
t242 = rSges(3,1) * t284 - rSges(3,2) * t283 - rSges(3,3) * t351;
t241 = Icges(3,1) * t286 - Icges(3,4) * t285 + Icges(3,5) * t352;
t240 = Icges(3,1) * t284 - Icges(3,4) * t283 - Icges(3,5) * t351;
t239 = Icges(3,4) * t286 - Icges(3,2) * t285 + Icges(3,6) * t352;
t238 = Icges(3,4) * t284 - Icges(3,2) * t283 - Icges(3,6) * t351;
t237 = Icges(3,5) * t286 - Icges(3,6) * t285 + Icges(3,3) * t352;
t236 = Icges(3,5) * t284 - Icges(3,6) * t283 - Icges(3,3) * t351;
t229 = t253 * t325 + t321 * t349;
t228 = -t253 * t321 + t325 * t349;
t224 = qJD(6) * t252 + t272;
t221 = t266 * t322 + t267 * t355;
t220 = -t266 * t355 + t267 * t322;
t219 = t264 * t322 + t265 * t355;
t218 = -t264 * t355 + t265 * t322;
t215 = rSges(4,1) * t267 - rSges(4,2) * t266 + rSges(4,3) * t285;
t214 = rSges(5,1) * t267 + rSges(5,2) * t285 + rSges(5,3) * t266;
t213 = rSges(4,1) * t265 - rSges(4,2) * t264 + rSges(4,3) * t283;
t212 = rSges(5,1) * t265 + rSges(5,2) * t283 + rSges(5,3) * t264;
t199 = pkin(5) * t253 + pkin(10) * t252;
t198 = t253 * rSges(6,1) - t252 * rSges(6,2) + rSges(6,3) * t349;
t197 = t221 * t325 - t285 * t321;
t196 = -t221 * t321 - t285 * t325;
t195 = t219 * t325 - t283 * t321;
t194 = -t219 * t321 - t283 * t325;
t193 = Icges(6,1) * t253 - Icges(6,4) * t252 + Icges(6,5) * t349;
t192 = Icges(6,4) * t253 - Icges(6,2) * t252 + Icges(6,6) * t349;
t191 = Icges(6,5) * t253 - Icges(6,6) * t252 + Icges(6,3) * t349;
t189 = -t242 * t310 + t276 * t298 + t338;
t188 = t310 * t243 - t299 * t276 + t336;
t187 = qJD(6) * t220 + t234;
t186 = qJD(6) * t218 + t233;
t185 = pkin(5) * t221 + pkin(10) * t220;
t184 = pkin(5) * t219 + pkin(10) * t218;
t183 = t242 * t299 - t243 * t298 + t337;
t182 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t252;
t181 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t252;
t180 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t252;
t179 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t252;
t178 = rSges(6,1) * t221 - rSges(6,2) * t220 - rSges(6,3) * t285;
t177 = rSges(6,1) * t219 - rSges(6,2) * t218 - rSges(6,3) * t283;
t176 = Icges(6,1) * t221 - Icges(6,4) * t220 - Icges(6,5) * t285;
t175 = Icges(6,1) * t219 - Icges(6,4) * t218 - Icges(6,5) * t283;
t174 = Icges(6,4) * t221 - Icges(6,2) * t220 - Icges(6,6) * t285;
t173 = Icges(6,4) * t219 - Icges(6,2) * t218 - Icges(6,6) * t283;
t172 = Icges(6,5) * t221 - Icges(6,6) * t220 - Icges(6,3) * t285;
t171 = Icges(6,5) * t219 - Icges(6,6) * t218 - Icges(6,3) * t283;
t170 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t220;
t169 = rSges(7,1) * t195 + rSges(7,2) * t194 + rSges(7,3) * t218;
t168 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t220;
t167 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t218;
t166 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t220;
t165 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t218;
t164 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t220;
t163 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t218;
t162 = -t213 * t287 + t251 * t262 + t335;
t161 = t287 * t215 - t263 * t251 + t331;
t160 = t213 * t263 - t215 * t262 + t334;
t159 = t250 * t262 + (-t212 - t225) * t287 + t333;
t158 = t287 * t214 + (-t250 - t256) * t263 + t330;
t157 = t212 * t263 + (-t214 - t226) * t262 + t332;
t156 = -t177 * t272 + t198 * t233 + t329;
t155 = t272 * t178 - t234 * t198 + t327;
t154 = t177 * t234 - t178 * t233 + t328;
t153 = -t169 * t224 + t182 * t186 - t184 * t272 + t199 * t233 + t329;
t152 = t224 * t170 - t187 * t182 + t272 * t185 - t234 * t199 + t327;
t151 = t169 * t187 - t170 * t186 + t184 * t234 - t185 * t233 + t328;
t1 = t299 * ((t237 * t352 - t285 * t239 + t286 * t241) * t299 + (t236 * t352 - t238 * t285 + t240 * t286) * t298 + (t273 * t352 - t274 * t285 + t275 * t286) * t310) / 0.2e1 + t272 * ((t172 * t349 - t252 * t174 + t253 * t176) * t234 + (t171 * t349 - t252 * t173 + t253 * t175) * t233 + (t191 * t349 - t252 * t192 + t253 * t193) * t272) / 0.2e1 + t298 * ((-t237 * t351 - t239 * t283 + t241 * t284) * t299 + (-t236 * t351 - t283 * t238 + t284 * t240) * t298 + (-t273 * t351 - t274 * t283 + t275 * t284) * t310) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(3) * (t183 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + t186 * ((t164 * t218 + t166 * t194 + t168 * t195) * t187 + (t218 * t163 + t194 * t165 + t195 * t167) * t186 + (t179 * t218 + t180 * t194 + t181 * t195) * t224) / 0.2e1 + t187 * ((t220 * t164 + t196 * t166 + t197 * t168) * t187 + (t163 * t220 + t165 * t196 + t167 * t197) * t186 + (t179 * t220 + t180 * t196 + t181 * t197) * t224) / 0.2e1 + t224 * ((t164 * t252 + t166 * t228 + t168 * t229) * t187 + (t163 * t252 + t165 * t228 + t167 * t229) * t186 + (t252 * t179 + t228 * t180 + t229 * t181) * t224) / 0.2e1 + m(2) * (t261 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + t233 * ((-t172 * t283 - t174 * t218 + t176 * t219) * t234 + (-t283 * t171 - t218 * t173 + t219 * t175) * t233 + (-t191 * t283 - t192 * t218 + t193 * t219) * t272) / 0.2e1 + t234 * ((-t285 * t172 - t220 * t174 + t176 * t221) * t234 + (-t171 * t285 - t173 * t220 + t175 * t221) * t233 + (-t191 * t285 - t192 * t220 + t193 * t221) * t272) / 0.2e1 + m(1) * (t295 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + t310 * (((t239 * t326 + t241 * t324) * t299 + (t238 * t326 + t240 * t324) * t298 + (t274 * t326 + t275 * t324) * t310) * t319 + (t236 * t298 + t237 * t299 + t273 * t310) * t354) / 0.2e1 + ((t362 * t264 + t360 * t265 + t361 * t283) * t287 + (t367 * t264 + t363 * t265 + t365 * t283) * t263 + (t368 * t264 + t364 * t265 + t366 * t283) * t262) * t262 / 0.2e1 + ((t362 * t266 + t360 * t267 + t361 * t285) * t287 + (t367 * t266 + t363 * t267 + t365 * t285) * t263 + (t368 * t266 + t364 * t267 + t366 * t285) * t262) * t263 / 0.2e1 + ((t362 * t290 + t360 * t291 - t361 * t349) * t287 + (t367 * t290 + t363 * t291 - t365 * t349) * t263 + (t368 * t290 + t364 * t291 - t366 * t349) * t262) * t287 / 0.2e1 + ((-t302 * t318 + t304 * t320 + Icges(1,4)) * V_base(5) + (-t318 * t303 + t320 * t305 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t320 * t302 + t318 * t304 + Icges(1,2)) * V_base(5) + (t303 * t320 + t305 * t318 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t318 + Icges(2,6) * t320 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t320 - Icges(2,6) * t318 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
