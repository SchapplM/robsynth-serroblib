% Calculate kinetic energy for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:10
% EndTime: 2019-03-09 17:14:14
% DurationCPUTime: 3.52s
% Computational Cost: add. (1630->321), mult. (3375->457), div. (0->0), fcn. (3655->8), ass. (0->152)
t345 = Icges(4,1) + Icges(5,1);
t344 = Icges(6,1) + Icges(7,1);
t343 = -Icges(4,4) + Icges(5,5);
t342 = Icges(5,4) + Icges(4,5);
t341 = Icges(6,4) + Icges(7,4);
t340 = Icges(6,5) + Icges(7,5);
t339 = Icges(4,2) + Icges(5,3);
t338 = Icges(6,2) + Icges(7,2);
t337 = -Icges(5,6) + Icges(4,6);
t336 = Icges(6,6) + Icges(7,6);
t335 = -Icges(4,3) - Icges(5,2);
t334 = Icges(6,3) + Icges(7,3);
t262 = sin(qJ(3));
t266 = cos(qJ(3));
t268 = cos(qJ(1));
t264 = sin(qJ(1));
t267 = cos(qJ(2));
t298 = t264 * t267;
t219 = t262 * t298 + t266 * t268;
t220 = -t262 * t268 + t266 * t298;
t261 = sin(qJ(5));
t265 = cos(qJ(5));
t178 = t219 * t265 - t220 * t261;
t303 = t219 * t261;
t179 = t220 * t265 + t303;
t263 = sin(qJ(2));
t300 = t263 * t264;
t333 = t336 * t178 + t340 * t179 - t334 * t300;
t297 = t267 * t268;
t221 = t262 * t297 - t264 * t266;
t222 = t264 * t262 + t266 * t297;
t180 = t221 * t265 - t222 * t261;
t302 = t221 * t261;
t181 = t222 * t265 + t302;
t299 = t263 * t268;
t332 = t336 * t180 + t340 * t181 - t334 * t299;
t331 = t338 * t178 + t341 * t179 - t336 * t300;
t330 = t338 * t180 + t341 * t181 - t336 * t299;
t329 = t341 * t178 + t344 * t179 - t340 * t300;
t328 = t341 * t180 + t344 * t181 - t340 * t299;
t212 = (-t261 * t266 + t262 * t265) * t263;
t301 = t261 * t262;
t213 = (t265 * t266 + t301) * t263;
t327 = t336 * t212 + t340 * t213 + t334 * t267;
t326 = t338 * t212 + t341 * t213 + t336 * t267;
t325 = t341 * t212 + t344 * t213 + t340 * t267;
t324 = t339 * t219 + t343 * t220 - t337 * t300;
t323 = t339 * t221 + t343 * t222 - t337 * t299;
t322 = -t337 * t219 + t342 * t220 - t335 * t300;
t321 = -t337 * t221 + t342 * t222 - t335 * t299;
t320 = t343 * t219 + t345 * t220 + t342 * t300;
t319 = t343 * t221 + t345 * t222 + t342 * t299;
t318 = t337 * t267 + (t339 * t262 + t343 * t266) * t263;
t317 = t335 * t267 + (-t337 * t262 + t342 * t266) * t263;
t316 = -t342 * t267 + (t343 * t262 + t345 * t266) * t263;
t308 = pkin(5) * t265;
t306 = Icges(2,4) * t264;
t305 = Icges(3,4) * t263;
t304 = Icges(3,4) * t267;
t287 = t263 * qJ(6);
t296 = rSges(7,1) * t179 + rSges(7,2) * t178 - rSges(7,3) * t300 + pkin(5) * t303 + t220 * t308 - t264 * t287;
t295 = t181 * rSges(7,1) + t180 * rSges(7,2) - rSges(7,3) * t299 + pkin(5) * t302 + t222 * t308 - t268 * t287;
t294 = rSges(7,1) * t213 + rSges(7,2) * t212 + (pkin(5) * t301 + t266 * t308) * t263 + (rSges(7,3) + qJ(6)) * t267;
t293 = qJD(3) * t263;
t292 = qJD(5) * t263;
t291 = qJD(6) * t263;
t290 = V_base(5) * pkin(6) + V_base(1);
t250 = qJD(2) * t264 + V_base(4);
t256 = V_base(6) + qJD(1);
t218 = t268 * t293 + t250;
t286 = pkin(2) * t267 + pkin(8) * t263;
t249 = -qJD(2) * t268 + V_base(5);
t285 = rSges(3,1) * t267 - rSges(3,2) * t263;
t284 = Icges(3,1) * t267 - t305;
t283 = -Icges(3,2) * t263 + t304;
t282 = Icges(3,5) * t267 - Icges(3,6) * t263;
t217 = t264 * t293 + t249;
t248 = pkin(1) * t268 + t264 * pkin(7);
t281 = -V_base(4) * pkin(6) + t256 * t248 + V_base(2);
t247 = t264 * pkin(1) - pkin(7) * t268;
t280 = V_base(4) * t247 - t248 * V_base(5) + V_base(3);
t225 = t286 * t264;
t246 = pkin(2) * t263 - pkin(8) * t267;
t279 = t249 * t246 + (-t225 - t247) * t256 + t290;
t278 = (-Icges(3,3) * t268 + t264 * t282) * t249 + (Icges(3,3) * t264 + t268 * t282) * t250 + (Icges(3,5) * t263 + Icges(3,6) * t267) * t256;
t226 = t286 * t268;
t277 = t256 * t226 - t246 * t250 + t281;
t223 = (pkin(3) * t266 + qJ(4) * t262) * t263;
t276 = qJD(4) * t221 + t217 * t223 + t279;
t275 = t250 * t225 - t226 * t249 + t280;
t183 = pkin(3) * t222 + qJ(4) * t221;
t242 = -qJD(3) * t267 + t256;
t274 = qJD(4) * t219 + t242 * t183 + t277;
t182 = pkin(3) * t220 + qJ(4) * t219;
t273 = qJD(4) * t263 * t262 + t218 * t182 + t275;
t190 = pkin(4) * t220 - pkin(9) * t300;
t229 = pkin(4) * t263 * t266 + pkin(9) * t267;
t272 = t217 * t229 + (-t182 - t190) * t242 + t276;
t191 = t222 * pkin(4) - pkin(9) * t299;
t271 = t242 * t191 + (-t223 - t229) * t218 + t274;
t270 = t218 * t190 + (-t183 - t191) * t217 + t273;
t202 = -Icges(3,6) * t268 + t264 * t283;
t203 = Icges(3,6) * t264 + t268 * t283;
t206 = -Icges(3,5) * t268 + t264 * t284;
t207 = Icges(3,5) * t264 + t268 * t284;
t236 = Icges(3,2) * t267 + t305;
t239 = Icges(3,1) * t263 + t304;
t269 = (-t203 * t263 + t207 * t267) * t250 + (-t202 * t263 + t206 * t267) * t249 + (-t236 * t263 + t239 * t267) * t256;
t258 = Icges(2,4) * t268;
t245 = rSges(2,1) * t268 - t264 * rSges(2,2);
t244 = t264 * rSges(2,1) + rSges(2,2) * t268;
t243 = rSges(3,1) * t263 + rSges(3,2) * t267;
t241 = Icges(2,1) * t268 - t306;
t240 = Icges(2,1) * t264 + t258;
t238 = -Icges(2,2) * t264 + t258;
t237 = Icges(2,2) * t268 + t306;
t232 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t231 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t230 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t227 = (-qJD(3) + qJD(5)) * t267 + t256;
t211 = t264 * rSges(3,3) + t268 * t285;
t210 = -rSges(3,3) * t268 + t264 * t285;
t209 = -rSges(4,3) * t267 + (rSges(4,1) * t266 - rSges(4,2) * t262) * t263;
t208 = -rSges(5,2) * t267 + (rSges(5,1) * t266 + rSges(5,3) * t262) * t263;
t194 = -t268 * t292 + t218;
t193 = -t264 * t292 + t217;
t189 = V_base(5) * rSges(2,3) - t244 * t256 + t290;
t188 = t245 * t256 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t186 = t244 * V_base(4) - t245 * V_base(5) + V_base(3);
t175 = t222 * rSges(4,1) - t221 * rSges(4,2) + rSges(4,3) * t299;
t174 = t222 * rSges(5,1) + rSges(5,2) * t299 + t221 * rSges(5,3);
t173 = rSges(4,1) * t220 - rSges(4,2) * t219 + rSges(4,3) * t300;
t172 = rSges(5,1) * t220 + rSges(5,2) * t300 + rSges(5,3) * t219;
t159 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t267;
t149 = t243 * t249 + (-t210 - t247) * t256 + t290;
t148 = t211 * t256 - t243 * t250 + t281;
t145 = t210 * t250 - t211 * t249 + t280;
t144 = t181 * rSges(6,1) + t180 * rSges(6,2) - rSges(6,3) * t299;
t142 = rSges(6,1) * t179 + rSges(6,2) * t178 - rSges(6,3) * t300;
t128 = -t173 * t242 + t209 * t217 + t279;
t127 = t175 * t242 - t209 * t218 + t277;
t126 = t173 * t218 - t175 * t217 + t275;
t125 = t208 * t217 + (-t172 - t182) * t242 + t276;
t124 = t174 * t242 + (-t208 - t223) * t218 + t274;
t123 = t172 * t218 + (-t174 - t183) * t217 + t273;
t122 = -t142 * t227 + t159 * t193 + t272;
t121 = t144 * t227 - t159 * t194 + t271;
t120 = t142 * t194 - t144 * t193 + t270;
t119 = t193 * t294 - t227 * t296 - t268 * t291 + t272;
t118 = -t194 * t294 + t227 * t295 - t264 * t291 + t271;
t117 = qJD(6) * t267 - t193 * t295 + t194 * t296 + t270;
t1 = t250 * (t264 * t278 + t268 * t269) / 0.2e1 + t249 * (t264 * t269 - t278 * t268) / 0.2e1 + m(1) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + m(2) * (t186 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(3) * (t145 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(5) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(4) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(7) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + ((t178 * t326 + t179 * t325 - t300 * t327) * t227 + (t178 * t330 + t179 * t328 - t300 * t332) * t194 + (t331 * t178 + t329 * t179 - t333 * t300) * t193) * t193 / 0.2e1 + ((t180 * t326 + t181 * t325 - t299 * t327) * t227 + (t330 * t180 + t328 * t181 - t332 * t299) * t194 + (t331 * t180 + t329 * t181 - t299 * t333) * t193) * t194 / 0.2e1 + ((t219 * t318 + t220 * t316 + t300 * t317) * t242 + (t219 * t323 + t220 * t319 + t300 * t321) * t218 + (t324 * t219 + t320 * t220 + t322 * t300) * t217) * t217 / 0.2e1 + ((t221 * t318 + t222 * t316 + t299 * t317) * t242 + (t323 * t221 + t319 * t222 + t321 * t299) * t218 + (t221 * t324 + t222 * t320 + t299 * t322) * t217) * t218 / 0.2e1 + ((t326 * t212 + t325 * t213 + t327 * t267) * t227 + (t212 * t330 + t213 * t328 + t332 * t267) * t194 + (t331 * t212 + t329 * t213 + t333 * t267) * t193) * t227 / 0.2e1 + ((-t217 * t322 - t218 * t321 - t242 * t317) * t267 + ((t262 * t318 + t266 * t316) * t242 + (t262 * t323 + t266 * t319) * t218 + (t262 * t324 + t266 * t320) * t217) * t263) * t242 / 0.2e1 + ((t203 * t267 + t207 * t263) * t250 + (t202 * t267 + t206 * t263) * t249 + (t267 * t236 + t263 * t239 + Icges(2,3)) * t256) * t256 / 0.2e1 + ((-t264 * t237 + t240 * t268 + Icges(1,4)) * V_base(5) + (-t264 * t238 + t268 * t241 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t268 * t237 + t264 * t240 + Icges(1,2)) * V_base(5) + (t238 * t268 + t264 * t241 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t256 * (Icges(2,5) * t268 - Icges(2,6) * t264) + V_base(5) * t256 * (Icges(2,5) * t264 + Icges(2,6) * t268) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
