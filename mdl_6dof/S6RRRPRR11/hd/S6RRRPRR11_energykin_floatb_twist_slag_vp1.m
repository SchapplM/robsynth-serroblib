% Calculate kinetic energy for
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:44
% EndTime: 2019-03-09 19:24:49
% DurationCPUTime: 5.92s
% Computational Cost: add. (2975->392), mult. (6948->569), div. (0->0), fcn. (8499->12), ass. (0->169)
t372 = Icges(4,1) + Icges(5,1);
t371 = -Icges(4,4) + Icges(5,5);
t370 = Icges(5,4) + Icges(4,5);
t369 = Icges(4,2) + Icges(5,3);
t368 = Icges(5,2) + Icges(4,3);
t367 = Icges(4,6) - Icges(5,6);
t323 = sin(qJ(2));
t324 = sin(qJ(1));
t326 = cos(qJ(2));
t327 = cos(qJ(1));
t352 = cos(pkin(6));
t340 = t327 * t352;
t288 = t323 * t340 + t324 * t326;
t322 = sin(qJ(3));
t319 = sin(pkin(6));
t354 = cos(qJ(3));
t343 = t319 * t354;
t263 = t288 * t322 + t327 * t343;
t348 = t319 * t327;
t264 = t288 * t354 - t322 * t348;
t287 = t323 * t324 - t326 * t340;
t366 = t369 * t263 + t371 * t264 - t367 * t287;
t341 = t324 * t352;
t290 = -t323 * t341 + t327 * t326;
t265 = t290 * t322 - t324 * t343;
t350 = t319 * t324;
t266 = t290 * t354 + t322 * t350;
t289 = t327 * t323 + t326 * t341;
t365 = t369 * t265 + t371 * t266 - t367 * t289;
t364 = -t367 * t263 + t370 * t264 + t368 * t287;
t363 = -t367 * t265 + t370 * t266 + t368 * t289;
t362 = t371 * t263 + t372 * t264 + t370 * t287;
t361 = t371 * t265 + t372 * t266 + t370 * t289;
t285 = t319 * t322 * t323 - t352 * t354;
t286 = t322 * t352 + t323 * t343;
t349 = t319 * t326;
t360 = t369 * t285 + t371 * t286 + t367 * t349;
t359 = -t367 * t285 + t370 * t286 - t368 * t349;
t358 = t371 * t285 + t372 * t286 - t370 * t349;
t353 = cos(qJ(5));
t351 = Icges(2,4) * t324;
t347 = qJD(2) * t319;
t346 = V_base(5) * pkin(7) + V_base(1);
t342 = t352 * pkin(8);
t299 = t324 * t347 + V_base(4);
t316 = V_base(6) + qJD(1);
t262 = qJD(3) * t289 + t299;
t300 = qJD(2) * t352 + t316;
t234 = -qJD(5) * t289 + t262;
t298 = -t327 * t347 + V_base(5);
t293 = t324 * pkin(1) - pkin(8) * t348;
t339 = -t293 * t316 + V_base(5) * t342 + t346;
t294 = pkin(1) * t327 + pkin(8) * t350;
t338 = V_base(4) * t293 - t294 * V_base(5) + V_base(3);
t261 = qJD(3) * t287 + t298;
t233 = -qJD(5) * t287 + t261;
t283 = -qJD(3) * t349 + t300;
t272 = qJD(5) * t349 + t283;
t255 = pkin(2) * t288 + pkin(9) * t287;
t292 = (pkin(2) * t323 - pkin(9) * t326) * t319;
t337 = -t255 * t300 + t298 * t292 + t339;
t256 = pkin(2) * t290 + pkin(9) * t289;
t336 = t299 * t255 - t256 * t298 + t338;
t335 = t316 * t294 + V_base(2) + (-t342 - pkin(7)) * V_base(4);
t254 = pkin(3) * t286 + qJ(4) * t285;
t334 = qJD(4) * t265 + t261 * t254 + t337;
t225 = pkin(3) * t264 + qJ(4) * t263;
t333 = qJD(4) * t285 + t262 * t225 + t336;
t332 = t300 * t256 - t299 * t292 + t335;
t226 = pkin(3) * t266 + qJ(4) * t265;
t331 = qJD(4) * t263 + t283 * t226 + t332;
t230 = pkin(4) * t264 - pkin(10) * t287;
t271 = pkin(4) * t286 + pkin(10) * t349;
t330 = t261 * t271 + (-t225 - t230) * t283 + t334;
t231 = pkin(4) * t266 - pkin(10) * t289;
t329 = t262 * t230 + (-t226 - t231) * t261 + t333;
t328 = t283 * t231 + (-t254 - t271) * t262 + t331;
t325 = cos(qJ(6));
t321 = sin(qJ(5));
t320 = sin(qJ(6));
t317 = Icges(2,4) * t327;
t308 = rSges(2,1) * t327 - t324 * rSges(2,2);
t307 = t324 * rSges(2,1) + rSges(2,2) * t327;
t306 = Icges(2,1) * t327 - t351;
t305 = Icges(2,1) * t324 + t317;
t304 = -Icges(2,2) * t324 + t317;
t303 = Icges(2,2) * t327 + t351;
t297 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t296 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t295 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t276 = t352 * rSges(3,3) + (rSges(3,1) * t323 + rSges(3,2) * t326) * t319;
t275 = Icges(3,5) * t352 + (Icges(3,1) * t323 + Icges(3,4) * t326) * t319;
t274 = Icges(3,6) * t352 + (Icges(3,4) * t323 + Icges(3,2) * t326) * t319;
t273 = Icges(3,3) * t352 + (Icges(3,5) * t323 + Icges(3,6) * t326) * t319;
t270 = V_base(5) * rSges(2,3) - t307 * t316 + t346;
t269 = t308 * t316 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t267 = t307 * V_base(4) - t308 * V_base(5) + V_base(3);
t253 = t285 * t321 + t286 * t353;
t252 = -t285 * t353 + t286 * t321;
t251 = rSges(3,1) * t290 - rSges(3,2) * t289 + rSges(3,3) * t350;
t250 = t288 * rSges(3,1) - t287 * rSges(3,2) - rSges(3,3) * t348;
t249 = Icges(3,1) * t290 - Icges(3,4) * t289 + Icges(3,5) * t350;
t248 = Icges(3,1) * t288 - Icges(3,4) * t287 - Icges(3,5) * t348;
t247 = Icges(3,4) * t290 - Icges(3,2) * t289 + Icges(3,6) * t350;
t246 = Icges(3,4) * t288 - Icges(3,2) * t287 - Icges(3,6) * t348;
t245 = Icges(3,5) * t290 - Icges(3,6) * t289 + Icges(3,3) * t350;
t244 = Icges(3,5) * t288 - Icges(3,6) * t287 - Icges(3,3) * t348;
t243 = rSges(4,1) * t286 - rSges(4,2) * t285 - rSges(4,3) * t349;
t242 = rSges(5,1) * t286 - rSges(5,2) * t349 + rSges(5,3) * t285;
t229 = t253 * t325 + t320 * t349;
t228 = -t253 * t320 + t325 * t349;
t224 = t265 * t321 + t266 * t353;
t223 = -t265 * t353 + t266 * t321;
t222 = t263 * t321 + t264 * t353;
t221 = -t263 * t353 + t264 * t321;
t218 = qJD(6) * t252 + t272;
t216 = rSges(4,1) * t266 - rSges(4,2) * t265 + rSges(4,3) * t289;
t215 = rSges(5,1) * t266 + rSges(5,2) * t289 + rSges(5,3) * t265;
t214 = rSges(4,1) * t264 - rSges(4,2) * t263 + rSges(4,3) * t287;
t213 = rSges(5,1) * t264 + rSges(5,2) * t287 + rSges(5,3) * t263;
t199 = pkin(5) * t253 + pkin(11) * t252;
t198 = t224 * t325 - t289 * t320;
t197 = -t224 * t320 - t289 * t325;
t196 = t222 * t325 - t287 * t320;
t195 = -t222 * t320 - t287 * t325;
t194 = rSges(6,1) * t253 - rSges(6,2) * t252 + rSges(6,3) * t349;
t193 = Icges(6,1) * t253 - Icges(6,4) * t252 + Icges(6,5) * t349;
t192 = Icges(6,4) * t253 - Icges(6,2) * t252 + Icges(6,6) * t349;
t191 = Icges(6,5) * t253 - Icges(6,6) * t252 + Icges(6,3) * t349;
t189 = qJD(6) * t223 + t234;
t188 = qJD(6) * t221 + t233;
t187 = -t250 * t300 + t276 * t298 + t339;
t186 = t300 * t251 - t299 * t276 + t335;
t185 = pkin(5) * t224 + pkin(11) * t223;
t184 = pkin(5) * t222 + pkin(11) * t221;
t183 = t250 * t299 - t251 * t298 + t338;
t182 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t252;
t181 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t252;
t180 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t252;
t179 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t252;
t178 = rSges(6,1) * t224 - rSges(6,2) * t223 - rSges(6,3) * t289;
t177 = rSges(6,1) * t222 - rSges(6,2) * t221 - rSges(6,3) * t287;
t176 = Icges(6,1) * t224 - Icges(6,4) * t223 - Icges(6,5) * t289;
t175 = Icges(6,1) * t222 - Icges(6,4) * t221 - Icges(6,5) * t287;
t174 = Icges(6,4) * t224 - Icges(6,2) * t223 - Icges(6,6) * t289;
t173 = Icges(6,4) * t222 - Icges(6,2) * t221 - Icges(6,6) * t287;
t172 = Icges(6,5) * t224 - Icges(6,6) * t223 - Icges(6,3) * t289;
t171 = Icges(6,5) * t222 - Icges(6,6) * t221 - Icges(6,3) * t287;
t170 = rSges(7,1) * t198 + rSges(7,2) * t197 + rSges(7,3) * t223;
t169 = rSges(7,1) * t196 + rSges(7,2) * t195 + rSges(7,3) * t221;
t168 = Icges(7,1) * t198 + Icges(7,4) * t197 + Icges(7,5) * t223;
t167 = Icges(7,1) * t196 + Icges(7,4) * t195 + Icges(7,5) * t221;
t166 = Icges(7,4) * t198 + Icges(7,2) * t197 + Icges(7,6) * t223;
t165 = Icges(7,4) * t196 + Icges(7,2) * t195 + Icges(7,6) * t221;
t164 = Icges(7,5) * t198 + Icges(7,6) * t197 + Icges(7,3) * t223;
t163 = Icges(7,5) * t196 + Icges(7,6) * t195 + Icges(7,3) * t221;
t162 = -t214 * t283 + t243 * t261 + t337;
t161 = t283 * t216 - t262 * t243 + t332;
t160 = t214 * t262 - t216 * t261 + t336;
t159 = t242 * t261 + (-t213 - t225) * t283 + t334;
t158 = t283 * t215 + (-t242 - t254) * t262 + t331;
t157 = t213 * t262 + (-t215 - t226) * t261 + t333;
t156 = -t177 * t272 + t194 * t233 + t330;
t155 = t272 * t178 - t234 * t194 + t328;
t154 = t177 * t234 - t178 * t233 + t329;
t153 = -t169 * t218 + t182 * t188 - t184 * t272 + t199 * t233 + t330;
t152 = t218 * t170 - t189 * t182 + t272 * t185 - t234 * t199 + t328;
t151 = t169 * t189 - t170 * t188 + t184 * t234 - t185 * t233 + t329;
t1 = t300 * (((t247 * t326 + t249 * t323) * t299 + (t246 * t326 + t248 * t323) * t298 + (t274 * t326 + t275 * t323) * t300) * t319 + (t244 * t298 + t245 * t299 + t273 * t300) * t352) / 0.2e1 + t272 * ((t172 * t349 - t174 * t252 + t176 * t253) * t234 + (t171 * t349 - t173 * t252 + t175 * t253) * t233 + (t191 * t349 - t192 * t252 + t193 * t253) * t272) / 0.2e1 + t299 * ((t245 * t350 - t247 * t289 + t249 * t290) * t299 + (t244 * t350 - t246 * t289 + t248 * t290) * t298 + (t273 * t350 - t274 * t289 + t275 * t290) * t300) / 0.2e1 + t298 * ((-t245 * t348 - t287 * t247 + t288 * t249) * t299 + (-t244 * t348 - t287 * t246 + t288 * t248) * t298 + (-t273 * t348 - t287 * t274 + t288 * t275) * t300) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(3) * (t183 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + t188 * ((t164 * t221 + t166 * t195 + t168 * t196) * t189 + (t221 * t163 + t195 * t165 + t196 * t167) * t188 + (t179 * t221 + t180 * t195 + t181 * t196) * t218) / 0.2e1 + t189 * ((t164 * t223 + t166 * t197 + t168 * t198) * t189 + (t163 * t223 + t165 * t197 + t167 * t198) * t188 + (t179 * t223 + t180 * t197 + t181 * t198) * t218) / 0.2e1 + t218 * ((t164 * t252 + t166 * t228 + t168 * t229) * t189 + (t163 * t252 + t165 * t228 + t167 * t229) * t188 + (t179 * t252 + t180 * t228 + t181 * t229) * t218) / 0.2e1 + m(2) * (t267 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + t233 * ((-t172 * t287 - t174 * t221 + t176 * t222) * t234 + (-t171 * t287 - t173 * t221 + t175 * t222) * t233 + (-t191 * t287 - t192 * t221 + t193 * t222) * t272) / 0.2e1 + t234 * ((-t172 * t289 - t174 * t223 + t176 * t224) * t234 + (-t171 * t289 - t173 * t223 + t175 * t224) * t233 + (-t191 * t289 - t192 * t223 + t193 * t224) * t272) / 0.2e1 + m(1) * (t295 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + ((t263 * t360 + t264 * t358 + t287 * t359) * t283 + (t263 * t365 + t264 * t361 + t287 * t363) * t262 + (t366 * t263 + t362 * t264 + t364 * t287) * t261) * t261 / 0.2e1 + ((t265 * t360 + t266 * t358 + t289 * t359) * t283 + (t365 * t265 + t361 * t266 + t363 * t289) * t262 + (t265 * t366 + t362 * t266 + t364 * t289) * t261) * t262 / 0.2e1 + ((t360 * t285 + t358 * t286 - t359 * t349) * t283 + (t285 * t365 + t286 * t361 - t349 * t363) * t262 + (t285 * t366 + t362 * t286 - t364 * t349) * t261) * t283 / 0.2e1 + ((-t324 * t303 + t305 * t327 + Icges(1,4)) * V_base(5) + (-t324 * t304 + t306 * t327 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t303 * t327 + t324 * t305 + Icges(1,2)) * V_base(5) + (t304 * t327 + t324 * t306 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t324 + Icges(2,6) * t327) * V_base(5) + (Icges(2,5) * t327 - Icges(2,6) * t324) * V_base(4) + Icges(2,3) * t316 / 0.2e1) * t316;
T  = t1;
