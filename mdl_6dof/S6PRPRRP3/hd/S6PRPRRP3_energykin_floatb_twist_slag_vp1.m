% Calculate kinetic energy for
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:36
% EndTime: 2019-03-08 20:04:41
% DurationCPUTime: 4.44s
% Computational Cost: add. (3397->391), mult. (5858->547), div. (0->0), fcn. (6961->12), ass. (0->177)
t377 = Icges(6,1) + Icges(7,1);
t376 = Icges(6,4) + Icges(7,4);
t375 = Icges(6,5) + Icges(7,5);
t374 = Icges(6,2) + Icges(7,2);
t373 = Icges(6,6) + Icges(7,6);
t372 = Icges(6,3) + Icges(7,3);
t371 = rSges(7,3) + qJ(6);
t298 = sin(pkin(10));
t301 = cos(pkin(10));
t308 = cos(qJ(2));
t302 = cos(pkin(6));
t306 = sin(qJ(2));
t338 = t302 * t306;
t263 = t298 * t308 + t301 * t338;
t331 = pkin(11) + qJ(4);
t294 = sin(t331);
t322 = cos(t331);
t299 = sin(pkin(6));
t341 = t299 * t301;
t235 = t263 * t322 - t294 * t341;
t337 = t302 * t308;
t262 = t298 * t306 - t301 * t337;
t305 = sin(qJ(5));
t307 = cos(qJ(5));
t206 = -t235 * t305 + t262 * t307;
t345 = t262 * t305;
t207 = t235 * t307 + t345;
t321 = t299 * t322;
t234 = t263 * t294 + t301 * t321;
t369 = t373 * t206 + t375 * t207 + t372 * t234;
t265 = -t298 * t338 + t301 * t308;
t342 = t298 * t299;
t237 = t265 * t322 + t294 * t342;
t264 = t298 * t337 + t301 * t306;
t208 = -t237 * t305 + t264 * t307;
t344 = t264 * t305;
t209 = t237 * t307 + t344;
t236 = t265 * t294 - t298 * t321;
t368 = t373 * t208 + t375 * t209 + t372 * t236;
t367 = t374 * t206 + t376 * t207 + t373 * t234;
t366 = t374 * t208 + t376 * t209 + t373 * t236;
t365 = t376 * t206 + t377 * t207 + t375 * t234;
t364 = t376 * t208 + t377 * t209 + t375 * t236;
t254 = t302 * t294 + t306 * t321;
t339 = t299 * t308;
t238 = -t254 * t305 - t307 * t339;
t323 = t305 * t339;
t239 = t254 * t307 - t323;
t340 = t299 * t306;
t253 = t294 * t340 - t302 * t322;
t363 = t373 * t238 + t375 * t239 + t372 * t253;
t362 = t374 * t238 + t376 * t239 + t373 * t253;
t361 = t376 * t238 + t377 * t239 + t375 * t253;
t297 = sin(pkin(11));
t300 = cos(pkin(11));
t241 = -t263 * t297 - t300 * t341;
t324 = t297 * t341;
t242 = t263 * t300 - t324;
t193 = Icges(4,5) * t242 + Icges(4,6) * t241 + Icges(4,3) * t262;
t222 = Icges(3,4) * t263 - Icges(3,2) * t262 - Icges(3,6) * t341;
t360 = t193 - t222;
t243 = -t265 * t297 + t300 * t342;
t325 = t297 * t342;
t244 = t265 * t300 + t325;
t194 = Icges(4,5) * t244 + Icges(4,6) * t243 + Icges(4,3) * t264;
t223 = Icges(3,4) * t265 - Icges(3,2) * t264 + Icges(3,6) * t342;
t359 = t194 - t223;
t260 = -t297 * t340 + t300 * t302;
t343 = t297 * t302;
t261 = t300 * t340 + t343;
t217 = Icges(4,5) * t261 + Icges(4,6) * t260 - Icges(4,3) * t339;
t251 = Icges(3,6) * t302 + (Icges(3,4) * t306 + Icges(3,2) * t308) * t299;
t358 = t217 - t251;
t350 = pkin(7) * t302;
t349 = pkin(3) * t300;
t348 = pkin(5) * t307;
t346 = Icges(2,4) * t298;
t335 = rSges(7,1) * t207 + rSges(7,2) * t206 + pkin(5) * t345 + t234 * t371 + t235 * t348;
t334 = rSges(7,1) * t209 + rSges(7,2) * t208 + pkin(5) * t344 + t236 * t371 + t237 * t348;
t333 = rSges(7,1) * t239 + rSges(7,2) * t238 - pkin(5) * t323 + t253 * t371 + t254 * t348;
t332 = qJD(2) * t299;
t330 = V_base(5) * qJ(1) + V_base(1);
t326 = qJD(1) + V_base(3);
t277 = t298 * t332 + V_base(4);
t288 = qJD(2) * t302 + V_base(6);
t246 = qJD(4) * t264 + t277;
t276 = -t301 * t332 + V_base(5);
t245 = qJD(4) * t262 + t276;
t266 = -qJD(4) * t339 + t288;
t270 = pkin(1) * t298 - pkin(7) * t341;
t320 = -t270 * V_base(6) + V_base(5) * t350 + t330;
t271 = pkin(1) * t301 + pkin(7) * t342;
t319 = V_base(4) * t270 - V_base(5) * t271 + t326;
t318 = V_base(6) * t271 + V_base(2) + (-qJ(1) - t350) * V_base(4);
t269 = (pkin(2) * t306 - qJ(3) * t308) * t299;
t317 = qJD(3) * t264 + t276 * t269 + t320;
t233 = pkin(2) * t265 + qJ(3) * t264;
t316 = qJD(3) * t262 + t288 * t233 + t318;
t232 = pkin(2) * t263 + qJ(3) * t262;
t315 = -qJD(3) * t339 + t277 * t232 + t319;
t190 = -pkin(3) * t324 + pkin(8) * t262 + t263 * t349;
t228 = pkin(3) * t343 + (-pkin(8) * t308 + t306 * t349) * t299;
t314 = t276 * t228 + (-t190 - t232) * t288 + t317;
t191 = pkin(3) * t325 + pkin(8) * t264 + t265 * t349;
t313 = t288 * t191 + (-t228 - t269) * t277 + t316;
t312 = t277 * t190 + (-t191 - t233) * t276 + t315;
t202 = pkin(4) * t235 + pkin(9) * t234;
t227 = t254 * pkin(4) + t253 * pkin(9);
t311 = -t202 * t266 + t245 * t227 + t314;
t203 = pkin(4) * t237 + pkin(9) * t236;
t310 = t266 * t203 - t227 * t246 + t313;
t309 = t246 * t202 - t245 * t203 + t312;
t295 = Icges(2,4) * t301;
t285 = rSges(2,1) * t301 - rSges(2,2) * t298;
t284 = rSges(2,1) * t298 + rSges(2,2) * t301;
t283 = Icges(2,1) * t301 - t346;
t282 = Icges(2,1) * t298 + t295;
t281 = -Icges(2,2) * t298 + t295;
t280 = Icges(2,2) * t301 + t346;
t275 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t274 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t273 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t255 = t302 * rSges(3,3) + (rSges(3,1) * t306 + rSges(3,2) * t308) * t299;
t252 = Icges(3,5) * t302 + (Icges(3,1) * t306 + Icges(3,4) * t308) * t299;
t250 = Icges(3,3) * t302 + (Icges(3,5) * t306 + Icges(3,6) * t308) * t299;
t249 = V_base(5) * rSges(2,3) - t284 * V_base(6) + t330;
t248 = t285 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t240 = t284 * V_base(4) - t285 * V_base(5) + t326;
t231 = qJD(5) * t253 + t266;
t230 = rSges(3,1) * t265 - rSges(3,2) * t264 + rSges(3,3) * t342;
t229 = rSges(3,1) * t263 - rSges(3,2) * t262 - rSges(3,3) * t341;
t226 = t261 * rSges(4,1) + t260 * rSges(4,2) - rSges(4,3) * t339;
t225 = Icges(3,1) * t265 - Icges(3,4) * t264 + Icges(3,5) * t342;
t224 = Icges(3,1) * t263 - Icges(3,4) * t262 - Icges(3,5) * t341;
t221 = Icges(3,5) * t265 - Icges(3,6) * t264 + Icges(3,3) * t342;
t220 = Icges(3,5) * t263 - Icges(3,6) * t262 - Icges(3,3) * t341;
t219 = Icges(4,1) * t261 + Icges(4,4) * t260 - Icges(4,5) * t339;
t218 = Icges(4,4) * t261 + Icges(4,2) * t260 - Icges(4,6) * t339;
t214 = t254 * rSges(5,1) - t253 * rSges(5,2) - rSges(5,3) * t339;
t213 = Icges(5,1) * t254 - Icges(5,4) * t253 - Icges(5,5) * t339;
t212 = Icges(5,4) * t254 - Icges(5,2) * t253 - Icges(5,6) * t339;
t211 = Icges(5,5) * t254 - Icges(5,6) * t253 - Icges(5,3) * t339;
t205 = qJD(5) * t236 + t246;
t204 = qJD(5) * t234 + t245;
t200 = rSges(4,1) * t244 + rSges(4,2) * t243 + rSges(4,3) * t264;
t199 = rSges(4,1) * t242 + rSges(4,2) * t241 + rSges(4,3) * t262;
t198 = Icges(4,1) * t244 + Icges(4,4) * t243 + Icges(4,5) * t264;
t197 = Icges(4,1) * t242 + Icges(4,4) * t241 + Icges(4,5) * t262;
t196 = Icges(4,4) * t244 + Icges(4,2) * t243 + Icges(4,6) * t264;
t195 = Icges(4,4) * t242 + Icges(4,2) * t241 + Icges(4,6) * t262;
t189 = rSges(5,1) * t237 - rSges(5,2) * t236 + rSges(5,3) * t264;
t188 = rSges(5,1) * t235 - rSges(5,2) * t234 + rSges(5,3) * t262;
t187 = Icges(5,1) * t237 - Icges(5,4) * t236 + Icges(5,5) * t264;
t186 = Icges(5,1) * t235 - Icges(5,4) * t234 + Icges(5,5) * t262;
t185 = Icges(5,4) * t237 - Icges(5,2) * t236 + Icges(5,6) * t264;
t184 = Icges(5,4) * t235 - Icges(5,2) * t234 + Icges(5,6) * t262;
t183 = Icges(5,5) * t237 - Icges(5,6) * t236 + Icges(5,3) * t264;
t182 = Icges(5,5) * t235 - Icges(5,6) * t234 + Icges(5,3) * t262;
t181 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t253;
t169 = -t229 * t288 + t255 * t276 + t320;
t168 = t230 * t288 - t255 * t277 + t318;
t167 = t229 * t277 - t230 * t276 + t319;
t166 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t236;
t164 = rSges(6,1) * t207 + rSges(6,2) * t206 + rSges(6,3) * t234;
t148 = t226 * t276 + (-t199 - t232) * t288 + t317;
t147 = t200 * t288 + (-t226 - t269) * t277 + t316;
t146 = t277 * t199 + (-t200 - t233) * t276 + t315;
t145 = -t188 * t266 + t214 * t245 + t314;
t144 = t189 * t266 - t214 * t246 + t313;
t143 = t246 * t188 - t245 * t189 + t312;
t142 = -t164 * t231 + t181 * t204 + t311;
t141 = t166 * t231 - t181 * t205 + t310;
t140 = t205 * t164 - t204 * t166 + t309;
t139 = qJD(6) * t236 + t204 * t333 - t231 * t335 + t311;
t138 = qJD(6) * t234 - t205 * t333 + t231 * t334 + t310;
t137 = qJD(6) * t253 - t204 * t334 + t205 * t335 + t309;
t1 = t266 * ((-t183 * t339 - t253 * t185 + t254 * t187) * t246 + (-t182 * t339 - t253 * t184 + t254 * t186) * t245 + (-t211 * t339 - t253 * t212 + t254 * t213) * t266) / 0.2e1 + m(3) * (t167 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(5) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(4) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(6) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(7) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(2) * (t240 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t245 * ((t183 * t262 - t185 * t234 + t187 * t235) * t246 + (t182 * t262 - t184 * t234 + t186 * t235) * t245 + (t211 * t262 - t212 * t234 + t213 * t235) * t266) / 0.2e1 + t246 * ((t183 * t264 - t185 * t236 + t187 * t237) * t246 + (t182 * t264 - t184 * t236 + t186 * t237) * t245 + (t211 * t264 - t212 * t236 + t213 * t237) * t266) / 0.2e1 + m(1) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + ((t206 * t362 + t207 * t361 + t234 * t363) * t231 + (t206 * t366 + t207 * t364 + t234 * t368) * t205 + (t367 * t206 + t365 * t207 + t369 * t234) * t204) * t204 / 0.2e1 + ((t208 * t362 + t209 * t361 + t236 * t363) * t231 + (t366 * t208 + t364 * t209 + t368 * t236) * t205 + (t208 * t367 + t209 * t365 + t236 * t369) * t204) * t205 / 0.2e1 + ((t362 * t238 + t361 * t239 + t363 * t253) * t231 + (t238 * t366 + t239 * t364 + t253 * t368) * t205 + (t238 * t367 + t239 * t365 + t253 * t369) * t204) * t231 / 0.2e1 + ((t218 * t241 + t219 * t242 - t250 * t341 + t252 * t263 + t262 * t358) * t288 + (t196 * t241 + t198 * t242 - t221 * t341 + t225 * t263 + t262 * t359) * t277 + (t195 * t241 + t197 * t242 - t220 * t341 + t224 * t263 + t360 * t262) * t276) * t276 / 0.2e1 + ((t218 * t243 + t219 * t244 + t250 * t342 + t252 * t265 + t264 * t358) * t288 + (t196 * t243 + t198 * t244 + t221 * t342 + t225 * t265 + t359 * t264) * t277 + (t195 * t243 + t197 * t244 + t220 * t342 + t224 * t265 + t264 * t360) * t276) * t277 / 0.2e1 + ((t220 * t276 + t221 * t277 + t250 * t288) * t302 + ((t223 * t308 + t225 * t306) * t277 + (t222 * t308 + t224 * t306) * t276 + (t251 * t308 + t252 * t306) * t288) * t299 + (-t194 * t339 + t260 * t196 + t261 * t198) * t277 + (-t193 * t339 + t260 * t195 + t261 * t197) * t276 + (-t217 * t339 + t260 * t218 + t261 * t219) * t288) * t288 / 0.2e1 + ((-t280 * t298 + t282 * t301 + Icges(1,4)) * V_base(5) + (-t281 * t298 + t283 * t301 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t280 * t301 + t282 * t298 + Icges(1,2)) * V_base(5) + (t281 * t301 + t283 * t298 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t301 - Icges(2,6) * t298 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t298 + Icges(2,6) * t301 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
