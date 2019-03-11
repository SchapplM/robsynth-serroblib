% Calculate kinetic energy for
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP14_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP14_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:40
% EndTime: 2019-03-09 13:03:43
% DurationCPUTime: 4.01s
% Computational Cost: add. (2332->340), mult. (5309->475), div. (0->0), fcn. (6265->10), ass. (0->157)
t380 = Icges(3,1) + Icges(4,2);
t379 = Icges(6,1) + Icges(7,1);
t378 = Icges(3,4) + Icges(4,6);
t377 = -Icges(6,4) + Icges(7,5);
t376 = Icges(7,4) + Icges(6,5);
t375 = Icges(3,5) - Icges(4,4);
t374 = Icges(3,2) + Icges(4,3);
t373 = Icges(6,2) + Icges(7,3);
t372 = Icges(7,2) + Icges(6,3);
t371 = Icges(3,6) - Icges(4,5);
t370 = Icges(6,6) - Icges(7,6);
t369 = Icges(3,3) + Icges(4,1);
t368 = rSges(7,1) + pkin(5);
t367 = rSges(7,3) + qJ(6);
t304 = cos(pkin(6));
t308 = sin(qJ(1));
t309 = cos(qJ(2));
t332 = t308 * t309;
t307 = sin(qJ(2));
t310 = cos(qJ(1));
t333 = t307 * t310;
t272 = t304 * t332 + t333;
t306 = sin(qJ(4));
t303 = sin(pkin(6));
t342 = cos(qJ(4));
t323 = t303 * t342;
t241 = t272 * t306 + t308 * t323;
t331 = t309 * t310;
t334 = t307 * t308;
t273 = -t304 * t334 + t331;
t305 = sin(qJ(5));
t341 = cos(qJ(5));
t204 = t241 * t305 - t273 * t341;
t205 = t241 * t341 + t273 * t305;
t337 = t303 * t308;
t240 = -t272 * t342 + t306 * t337;
t366 = t204 * t373 + t205 * t377 - t240 * t370;
t270 = -t304 * t331 + t334;
t243 = t270 * t306 - t310 * t323;
t271 = t304 * t333 + t332;
t206 = t243 * t305 - t271 * t341;
t207 = t243 * t341 + t271 * t305;
t335 = t303 * t310;
t242 = t270 * t342 + t306 * t335;
t365 = t206 * t373 + t207 * t377 + t242 * t370;
t364 = -t204 * t370 + t205 * t376 + t240 * t372;
t363 = -t206 * t370 + t207 * t376 - t242 * t372;
t362 = t204 * t377 + t205 * t379 + t240 * t376;
t361 = t206 * t377 + t207 * t379 - t242 * t376;
t336 = t303 * t309;
t269 = t304 * t342 - t306 * t336;
t338 = t303 * t307;
t236 = t269 * t305 - t338 * t341;
t237 = t269 * t341 + t305 * t338;
t268 = t304 * t306 + t309 * t323;
t360 = t236 * t373 + t237 * t377 - t268 * t370;
t359 = -t236 * t370 + t237 * t376 + t268 * t372;
t358 = t236 * t377 + t237 * t379 + t268 * t376;
t357 = t272 * t374 - t273 * t378 - t337 * t371;
t356 = t270 * t374 - t271 * t378 + t335 * t371;
t355 = -t378 * t272 + t273 * t380 + t375 * t337;
t354 = -t378 * t270 + t271 * t380 - t375 * t335;
t353 = -t272 * t371 + t273 * t375 + t337 * t369;
t352 = -t270 * t371 + t271 * t375 - t335 * t369;
t351 = t369 * t304 + (t307 * t375 + t309 * t371) * t303;
t350 = t371 * t304 + (t307 * t378 + t309 * t374) * t303;
t349 = t375 * t304 + (t307 * t380 + t378 * t309) * t303;
t340 = pkin(8) * t304;
t339 = Icges(2,4) * t308;
t330 = rSges(7,2) * t240 + t367 * t204 + t368 * t205;
t329 = -rSges(7,2) * t242 + t367 * t206 + t368 * t207;
t328 = rSges(7,2) * t268 + t367 * t236 + t368 * t237;
t327 = qJD(2) * t303;
t326 = V_base(5) * pkin(7) + V_base(1);
t283 = t308 * t327 + V_base(4);
t300 = V_base(6) + qJD(1);
t239 = qJD(4) * t273 + t283;
t284 = qJD(2) * t304 + t300;
t266 = qJD(4) * t338 + t284;
t282 = -t310 * t327 + V_base(5);
t277 = t308 * pkin(1) - pkin(8) * t335;
t322 = -t277 * t300 + V_base(5) * t340 + t326;
t278 = pkin(1) * t310 + pkin(8) * t337;
t321 = V_base(4) * t277 - t278 * V_base(5) + V_base(3);
t238 = qJD(4) * t271 + t282;
t274 = (pkin(2) * t307 - qJ(3) * t309) * t303;
t320 = qJD(3) * t272 + t282 * t274 + t322;
t319 = t300 * t278 + V_base(2) + (-pkin(7) - t340) * V_base(4);
t235 = pkin(2) * t273 + qJ(3) * t272;
t318 = qJD(3) * t270 + t284 * t235 + t319;
t234 = pkin(2) * t271 + qJ(3) * t270;
t317 = -qJD(3) * t336 + t283 * t234 + t321;
t249 = -pkin(3) * t335 + t271 * pkin(9);
t276 = pkin(3) * t304 + pkin(9) * t338;
t316 = t282 * t276 + (-t234 - t249) * t284 + t320;
t248 = pkin(3) * t337 + pkin(9) * t273;
t315 = t284 * t248 + (-t274 - t276) * t283 + t318;
t314 = t283 * t249 + (-t235 - t248) * t282 + t317;
t201 = pkin(4) * t243 - pkin(10) * t242;
t231 = pkin(4) * t269 + pkin(10) * t268;
t313 = -t201 * t266 + t238 * t231 + t316;
t200 = pkin(4) * t241 + pkin(10) * t240;
t312 = t266 * t200 - t231 * t239 + t315;
t311 = -t200 * t238 + t239 * t201 + t314;
t301 = Icges(2,4) * t310;
t292 = rSges(2,1) * t310 - t308 * rSges(2,2);
t291 = t308 * rSges(2,1) + rSges(2,2) * t310;
t290 = Icges(2,1) * t310 - t339;
t289 = Icges(2,1) * t308 + t301;
t288 = -Icges(2,2) * t308 + t301;
t287 = Icges(2,2) * t310 + t339;
t281 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t280 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t279 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t258 = rSges(4,1) * t304 + (-rSges(4,2) * t307 - rSges(4,3) * t309) * t303;
t257 = rSges(3,3) * t304 + (rSges(3,1) * t307 + rSges(3,2) * t309) * t303;
t247 = V_base(5) * rSges(2,3) - t291 * t300 + t326;
t246 = t292 * t300 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t244 = t291 * V_base(4) - t292 * V_base(5) + V_base(3);
t232 = qJD(5) * t268 + t266;
t229 = rSges(3,1) * t273 - rSges(3,2) * t272 + rSges(3,3) * t337;
t228 = t271 * rSges(3,1) - t270 * rSges(3,2) - rSges(3,3) * t335;
t227 = -rSges(4,1) * t335 - t271 * rSges(4,2) + t270 * rSges(4,3);
t226 = rSges(4,1) * t337 - rSges(4,2) * t273 + rSges(4,3) * t272;
t213 = rSges(5,1) * t269 - rSges(5,2) * t268 + rSges(5,3) * t338;
t212 = Icges(5,1) * t269 - Icges(5,4) * t268 + Icges(5,5) * t338;
t211 = Icges(5,4) * t269 - Icges(5,2) * t268 + Icges(5,6) * t338;
t210 = Icges(5,5) * t269 - Icges(5,6) * t268 + Icges(5,3) * t338;
t203 = qJD(5) * t240 + t239;
t202 = -qJD(5) * t242 + t238;
t197 = rSges(5,1) * t243 + rSges(5,2) * t242 + rSges(5,3) * t271;
t196 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t273;
t194 = Icges(5,1) * t243 + Icges(5,4) * t242 + Icges(5,5) * t271;
t193 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t273;
t192 = Icges(5,4) * t243 + Icges(5,2) * t242 + Icges(5,6) * t271;
t191 = Icges(5,4) * t241 - Icges(5,2) * t240 + Icges(5,6) * t273;
t190 = Icges(5,5) * t243 + Icges(5,6) * t242 + Icges(5,3) * t271;
t189 = Icges(5,5) * t241 - Icges(5,6) * t240 + Icges(5,3) * t273;
t188 = rSges(6,1) * t237 - rSges(6,2) * t236 + rSges(6,3) * t268;
t177 = -t228 * t284 + t257 * t282 + t322;
t176 = t229 * t284 - t257 * t283 + t319;
t175 = rSges(6,1) * t207 - rSges(6,2) * t206 - rSges(6,3) * t242;
t173 = rSges(6,1) * t205 - rSges(6,2) * t204 + rSges(6,3) * t240;
t159 = t228 * t283 - t229 * t282 + t321;
t158 = t258 * t282 + (-t227 - t234) * t284 + t320;
t157 = t226 * t284 + (-t258 - t274) * t283 + t318;
t156 = t227 * t283 + (-t226 - t235) * t282 + t317;
t155 = -t197 * t266 + t213 * t238 + t316;
t154 = t196 * t266 - t213 * t239 + t315;
t153 = -t196 * t238 + t197 * t239 + t314;
t152 = -t175 * t232 + t188 * t202 + t313;
t151 = t173 * t232 - t188 * t203 + t312;
t150 = -t173 * t202 + t175 * t203 + t311;
t149 = qJD(6) * t204 + t202 * t328 - t232 * t329 + t313;
t148 = qJD(6) * t206 - t203 * t328 + t232 * t330 + t312;
t147 = qJD(6) * t236 - t202 * t330 + t203 * t329 + t311;
t1 = t266 * ((t189 * t338 - t191 * t268 + t193 * t269) * t239 + (t190 * t338 - t192 * t268 + t194 * t269) * t238 + (t210 * t338 - t211 * t268 + t212 * t269) * t266) / 0.2e1 + m(7) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(5) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(4) * (t156 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(3) * (t159 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(2) * (t244 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + t238 * ((t189 * t271 + t191 * t242 + t193 * t243) * t239 + (t190 * t271 + t192 * t242 + t194 * t243) * t238 + (t210 * t271 + t211 * t242 + t212 * t243) * t266) / 0.2e1 + t239 * ((t189 * t273 - t191 * t240 + t193 * t241) * t239 + (t190 * t273 - t192 * t240 + t194 * t241) * t238 + (t210 * t273 - t211 * t240 + t212 * t241) * t266) / 0.2e1 + m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + ((t206 * t360 + t207 * t358 - t242 * t359) * t232 + (t206 * t366 + t362 * t207 - t364 * t242) * t203 + (t365 * t206 + t361 * t207 - t363 * t242) * t202) * t202 / 0.2e1 + ((t204 * t360 + t205 * t358 + t240 * t359) * t232 + (t366 * t204 + t362 * t205 + t364 * t240) * t203 + (t204 * t365 + t205 * t361 + t240 * t363) * t202) * t203 / 0.2e1 + ((t360 * t236 + t358 * t237 + t359 * t268) * t232 + (t236 * t366 + t362 * t237 + t364 * t268) * t203 + (t236 * t365 + t237 * t361 + t268 * t363) * t202) * t232 / 0.2e1 + ((-t270 * t350 + t271 * t349 - t335 * t351) * t284 + (t270 * t357 + t271 * t355 - t335 * t353) * t283 + (t356 * t270 + t354 * t271 - t352 * t335) * t282) * t282 / 0.2e1 + ((-t272 * t350 + t273 * t349 + t337 * t351) * t284 + (t357 * t272 + t355 * t273 + t353 * t337) * t283 + (t272 * t356 + t273 * t354 + t337 * t352) * t282) * t283 / 0.2e1 + ((t282 * t352 + t283 * t353 + t284 * t351) * t304 + ((t307 * t349 + t309 * t350) * t284 + (t307 * t355 - t309 * t357) * t283 + (t307 * t354 - t309 * t356) * t282) * t303) * t284 / 0.2e1 + ((-t308 * t287 + t289 * t310 + Icges(1,4)) * V_base(5) + (-t308 * t288 + t290 * t310 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t287 * t310 + t308 * t289 + Icges(1,2)) * V_base(5) + (t288 * t310 + t308 * t290 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t308 + Icges(2,6) * t310) * V_base(5) + (Icges(2,5) * t310 - Icges(2,6) * t308) * V_base(4) + Icges(2,3) * t300 / 0.2e1) * t300;
T  = t1;
