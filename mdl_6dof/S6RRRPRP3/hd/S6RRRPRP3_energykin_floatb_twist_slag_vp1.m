% Calculate kinetic energy for
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:02
% EndTime: 2019-03-09 16:39:05
% DurationCPUTime: 3.28s
% Computational Cost: add. (2192->335), mult. (2279->484), div. (0->0), fcn. (2183->10), ass. (0->167)
t343 = Icges(6,1) + Icges(7,1);
t342 = -Icges(6,4) + Icges(7,5);
t341 = Icges(7,4) + Icges(6,5);
t340 = Icges(6,2) + Icges(7,3);
t339 = -Icges(7,6) + Icges(6,6);
t338 = -Icges(6,3) - Icges(7,2);
t337 = rSges(7,1) + pkin(5);
t336 = rSges(7,3) + qJ(6);
t257 = pkin(10) + qJ(5);
t249 = cos(t257);
t258 = qJ(2) + qJ(3);
t254 = cos(t258);
t265 = cos(qJ(1));
t248 = sin(t257);
t263 = sin(qJ(1));
t311 = t248 * t263;
t195 = t249 * t265 + t254 * t311;
t303 = t263 * t249;
t196 = -t248 * t265 + t254 * t303;
t253 = sin(t258);
t310 = t253 * t263;
t335 = t195 * t340 + t196 * t342 - t310 * t339;
t308 = t254 * t265;
t197 = t248 * t308 - t303;
t198 = t249 * t308 + t311;
t309 = t253 * t265;
t334 = t197 * t340 + t198 * t342 - t309 * t339;
t333 = -t195 * t339 + t196 * t341 - t310 * t338;
t332 = -t197 * t339 + t198 * t341 - t309 * t338;
t331 = t342 * t195 + t196 * t343 + t341 * t310;
t330 = t342 * t197 + t198 * t343 + t341 * t309;
t329 = t339 * t254 + (t248 * t340 + t249 * t342) * t253;
t328 = t338 * t254 + (-t248 * t339 + t249 * t341) * t253;
t327 = -t341 * t254 + (t342 * t248 + t249 * t343) * t253;
t262 = sin(qJ(2));
t320 = pkin(2) * t262;
t264 = cos(qJ(2));
t319 = pkin(2) * t264;
t260 = cos(pkin(10));
t318 = pkin(4) * t260;
t316 = Icges(2,4) * t263;
t315 = Icges(3,4) * t262;
t314 = Icges(3,4) * t264;
t313 = Icges(4,4) * t253;
t312 = Icges(4,4) * t254;
t259 = sin(pkin(10));
t307 = t259 * t263;
t306 = t259 * t265;
t305 = t260 * t263;
t304 = t260 * t265;
t301 = rSges(7,2) * t310 + t336 * t195 + t337 * t196;
t300 = rSges(7,2) * t309 + t336 * t197 + t337 * t198;
t299 = -rSges(7,2) * t254 + (t336 * t248 + t337 * t249) * t253;
t181 = -pkin(8) * t265 + t263 * t319;
t241 = t263 * pkin(1) - t265 * pkin(7);
t298 = -t181 - t241;
t297 = qJD(4) * t253;
t296 = qJD(5) * t253;
t295 = V_base(5) * pkin(6) + V_base(1);
t288 = pkin(3) * t254 + qJ(4) * t253;
t207 = t288 * t263;
t292 = -t207 + t298;
t245 = qJD(2) * t263 + V_base(4);
t250 = V_base(6) + qJD(1);
t244 = -qJD(2) * t265 + V_base(5);
t291 = t244 * t320 + t295;
t220 = qJD(3) * t263 + t245;
t290 = rSges(3,1) * t264 - rSges(3,2) * t262;
t289 = rSges(4,1) * t254 - rSges(4,2) * t253;
t287 = Icges(3,1) * t264 - t315;
t286 = Icges(4,1) * t254 - t313;
t285 = -Icges(3,2) * t262 + t314;
t284 = -Icges(4,2) * t253 + t312;
t283 = Icges(3,5) * t264 - Icges(3,6) * t262;
t282 = Icges(4,5) * t254 - Icges(4,6) * t253;
t217 = pkin(3) * t253 - qJ(4) * t254;
t219 = V_base(5) + (-qJD(2) - qJD(3)) * t265;
t281 = t219 * t217 + t265 * t297 + t291;
t242 = t265 * pkin(1) + t263 * pkin(7);
t280 = -V_base(4) * pkin(6) + t250 * t242 + V_base(2);
t279 = V_base(4) * t241 - t242 * V_base(5) + V_base(3);
t278 = (-Icges(4,3) * t265 + t263 * t282) * t219 + (Icges(4,3) * t263 + t265 * t282) * t220 + (Icges(4,5) * t253 + Icges(4,6) * t254) * t250;
t277 = (-Icges(3,3) * t265 + t263 * t283) * t244 + (Icges(3,3) * t263 + t265 * t283) * t245 + (Icges(3,5) * t262 + Icges(3,6) * t264) * t250;
t276 = pkin(9) * t253 + t254 * t318;
t182 = pkin(8) * t263 + t265 * t319;
t275 = t245 * t181 - t182 * t244 + t279;
t274 = t250 * t182 - t245 * t320 + t280;
t208 = t288 * t265;
t273 = t250 * t208 + t263 * t297 + t274;
t148 = -pkin(4) * t306 + t263 * t276;
t161 = -pkin(9) * t254 + t253 * t318;
t272 = t219 * t161 + (-t148 + t292) * t250 + t281;
t271 = -qJD(4) * t254 + t220 * t207 + t275;
t149 = pkin(4) * t307 + t265 * t276;
t270 = t250 * t149 + (-t161 - t217) * t220 + t273;
t269 = t220 * t148 + (-t149 - t208) * t219 + t271;
t186 = -Icges(4,6) * t265 + t263 * t284;
t187 = Icges(4,6) * t263 + t265 * t284;
t188 = -Icges(4,5) * t265 + t263 * t286;
t189 = Icges(4,5) * t263 + t265 * t286;
t215 = Icges(4,2) * t254 + t313;
t216 = Icges(4,1) * t253 + t312;
t268 = (-t187 * t253 + t189 * t254) * t220 + (-t186 * t253 + t188 * t254) * t219 + (-t215 * t253 + t216 * t254) * t250;
t201 = -Icges(3,6) * t265 + t263 * t285;
t202 = Icges(3,6) * t263 + t265 * t285;
t203 = -Icges(3,5) * t265 + t263 * t287;
t204 = Icges(3,5) * t263 + t265 * t287;
t230 = Icges(3,2) * t264 + t315;
t233 = Icges(3,1) * t262 + t314;
t267 = (-t202 * t262 + t204 * t264) * t245 + (-t201 * t262 + t203 * t264) * t244 + (-t230 * t262 + t233 * t264) * t250;
t255 = Icges(2,4) * t265;
t238 = rSges(2,1) * t265 - rSges(2,2) * t263;
t237 = rSges(2,1) * t263 + rSges(2,2) * t265;
t236 = rSges(3,1) * t262 + rSges(3,2) * t264;
t235 = Icges(2,1) * t265 - t316;
t234 = Icges(2,1) * t263 + t255;
t232 = -Icges(2,2) * t263 + t255;
t231 = Icges(2,2) * t265 + t316;
t226 = -qJD(5) * t254 + t250;
t225 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t224 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t223 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t218 = rSges(4,1) * t253 + rSges(4,2) * t254;
t212 = t254 * t304 + t307;
t211 = -t254 * t306 + t305;
t210 = t254 * t305 - t306;
t209 = -t254 * t307 - t304;
t206 = rSges(3,3) * t263 + t265 * t290;
t205 = -rSges(3,3) * t265 + t263 * t290;
t194 = t265 * t296 + t220;
t193 = t263 * t296 + t219;
t192 = rSges(4,3) * t263 + t265 * t289;
t191 = -rSges(4,3) * t265 + t263 * t289;
t180 = -rSges(5,3) * t254 + (rSges(5,1) * t260 - rSges(5,2) * t259) * t253;
t179 = -Icges(5,5) * t254 + (Icges(5,1) * t260 - Icges(5,4) * t259) * t253;
t178 = -Icges(5,6) * t254 + (Icges(5,4) * t260 - Icges(5,2) * t259) * t253;
t177 = -Icges(5,3) * t254 + (Icges(5,5) * t260 - Icges(5,6) * t259) * t253;
t176 = V_base(5) * rSges(2,3) - t237 * t250 + t295;
t175 = t238 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = t237 * V_base(4) - t238 * V_base(5) + V_base(3);
t171 = -rSges(6,3) * t254 + (rSges(6,1) * t249 - rSges(6,2) * t248) * t253;
t157 = rSges(5,1) * t212 + rSges(5,2) * t211 + rSges(5,3) * t309;
t156 = rSges(5,1) * t210 + rSges(5,2) * t209 + rSges(5,3) * t310;
t155 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t309;
t154 = Icges(5,1) * t210 + Icges(5,4) * t209 + Icges(5,5) * t310;
t153 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t309;
t152 = Icges(5,4) * t210 + Icges(5,2) * t209 + Icges(5,6) * t310;
t151 = Icges(5,5) * t212 + Icges(5,6) * t211 + Icges(5,3) * t309;
t150 = Icges(5,5) * t210 + Icges(5,6) * t209 + Icges(5,3) * t310;
t146 = rSges(6,1) * t198 - rSges(6,2) * t197 + rSges(6,3) * t309;
t144 = rSges(6,1) * t196 - rSges(6,2) * t195 + rSges(6,3) * t310;
t129 = t236 * t244 + (-t205 - t241) * t250 + t295;
t128 = t206 * t250 - t236 * t245 + t280;
t127 = t205 * t245 - t206 * t244 + t279;
t126 = t218 * t219 + (-t191 + t298) * t250 + t291;
t125 = t192 * t250 - t218 * t220 + t274;
t124 = t191 * t220 - t192 * t219 + t275;
t123 = t180 * t219 + (-t156 + t292) * t250 + t281;
t122 = t157 * t250 + (-t180 - t217) * t220 + t273;
t121 = t156 * t220 + (-t157 - t208) * t219 + t271;
t120 = -t144 * t226 + t171 * t193 + t272;
t119 = t146 * t226 - t171 * t194 + t270;
t118 = t144 * t194 - t146 * t193 + t269;
t117 = qJD(6) * t197 + t193 * t299 - t226 * t301 + t272;
t116 = qJD(6) * t195 - t194 * t299 + t226 * t300 + t270;
t115 = qJD(6) * t248 * t253 - t193 * t300 + t194 * t301 + t269;
t1 = t245 * (t277 * t263 + t267 * t265) / 0.2e1 + t244 * (t267 * t263 - t277 * t265) / 0.2e1 + m(1) * (t223 ^ 2 + t224 ^ 2 + t225 ^ 2) / 0.2e1 + m(2) * (t173 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(5) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(4) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(7) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + ((t195 * t329 + t196 * t327 + t310 * t328) * t226 + (t195 * t334 + t196 * t330 + t310 * t332) * t194 + (t335 * t195 + t331 * t196 + t333 * t310) * t193) * t193 / 0.2e1 + ((t197 * t329 + t198 * t327 + t309 * t328) * t226 + (t334 * t197 + t330 * t198 + t332 * t309) * t194 + (t197 * t335 + t331 * t198 + t333 * t309) * t193) * t194 / 0.2e1 + ((t151 * t310 + t153 * t209 + t155 * t210) * t220 + (t150 * t310 + t209 * t152 + t210 * t154) * t219 + (t177 * t310 + t178 * t209 + t179 * t210) * t250 + t268 * t263 - t278 * t265) * t219 / 0.2e1 + ((t151 * t309 + t211 * t153 + t212 * t155) * t220 + (t150 * t309 + t152 * t211 + t154 * t212) * t219 + (t177 * t309 + t178 * t211 + t179 * t212) * t250 + t278 * t263 + t268 * t265) * t220 / 0.2e1 + ((-t193 * t333 - t194 * t332 - t226 * t328) * t254 + ((t248 * t329 + t249 * t327) * t226 + (t248 * t334 + t249 * t330) * t194 + (t248 * t335 + t331 * t249) * t193) * t253) * t226 / 0.2e1 + ((-t231 * t263 + t234 * t265 + Icges(1,4)) * V_base(5) + (-t263 * t232 + t265 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t265 * t231 + t263 * t234 + Icges(1,2)) * V_base(5) + (t232 * t265 + t235 * t263 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t150 * t219 - t151 * t220) * t254 + ((-t153 * t259 + t155 * t260) * t220 + (-t152 * t259 + t154 * t260) * t219) * t253 + (t202 * t264 + t204 * t262) * t245 + (t201 * t264 + t203 * t262) * t244 + (t187 * t254 + t189 * t253) * t220 + (t186 * t254 + t188 * t253) * t219 + (t264 * t230 + t262 * t233 + Icges(2,3) + (-t177 + t215) * t254 + (-t178 * t259 + t179 * t260 + t216) * t253) * t250) * t250 / 0.2e1 + t250 * V_base(4) * (Icges(2,5) * t265 - Icges(2,6) * t263) + t250 * V_base(5) * (Icges(2,5) * t263 + Icges(2,6) * t265) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
