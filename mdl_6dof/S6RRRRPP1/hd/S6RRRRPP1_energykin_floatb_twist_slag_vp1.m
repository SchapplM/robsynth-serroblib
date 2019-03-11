% Calculate kinetic energy for
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:43:58
% EndTime: 2019-03-09 20:44:00
% DurationCPUTime: 2.88s
% Computational Cost: add. (2246->326), mult. (2324->462), div. (0->0), fcn. (2228->10), ass. (0->163)
t340 = Icges(6,1) + Icges(7,1);
t339 = -Icges(6,4) + Icges(7,5);
t338 = Icges(7,4) + Icges(6,5);
t337 = Icges(6,2) + Icges(7,3);
t336 = -Icges(7,6) + Icges(6,6);
t335 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t334 = rSges(7,1) + pkin(5);
t333 = rSges(7,3) + qJ(6);
t257 = qJ(4) + pkin(10);
t249 = cos(t257);
t258 = qJ(2) + qJ(3);
t254 = cos(t258);
t265 = cos(qJ(1));
t248 = sin(t257);
t262 = sin(qJ(1));
t309 = t248 * t262;
t195 = t249 * t265 + t254 * t309;
t303 = t262 * t249;
t196 = -t248 * t265 + t254 * t303;
t253 = sin(t258);
t308 = t253 * t262;
t332 = t337 * t195 + t339 * t196 - t336 * t308;
t306 = t254 * t265;
t197 = t248 * t306 - t303;
t198 = t249 * t306 + t309;
t307 = t253 * t265;
t331 = t337 * t197 + t339 * t198 - t336 * t307;
t330 = t339 * t195 + t340 * t196 + t338 * t308;
t329 = t339 * t197 + t340 * t198 + t338 * t307;
t328 = t336 * t254 + (t337 * t248 + t339 * t249) * t253;
t327 = -t338 * t254 + (t339 * t248 + t340 * t249) * t253;
t263 = cos(qJ(4));
t301 = t263 * t265;
t260 = sin(qJ(4));
t305 = t260 * t262;
t209 = -t254 * t305 - t301;
t302 = t262 * t263;
t304 = t260 * t265;
t210 = t254 * t302 - t304;
t326 = Icges(5,5) * t210 + Icges(5,6) * t209 - t336 * t195 + t338 * t196 - t335 * t308;
t211 = -t254 * t304 + t302;
t212 = t254 * t301 + t305;
t325 = Icges(5,5) * t212 + Icges(5,6) * t211 - t336 * t197 + t338 * t198 - t335 * t307;
t324 = t335 * t254 + (Icges(5,5) * t263 - Icges(5,6) * t260 - t336 * t248 + t338 * t249) * t253;
t261 = sin(qJ(2));
t319 = pkin(2) * t261;
t264 = cos(qJ(2));
t318 = pkin(2) * t264;
t317 = pkin(4) * t263;
t314 = Icges(2,4) * t262;
t313 = Icges(3,4) * t261;
t312 = Icges(3,4) * t264;
t311 = Icges(4,4) * t253;
t310 = Icges(4,4) * t254;
t300 = rSges(7,2) * t308 + t333 * t195 + t334 * t196;
t299 = rSges(7,2) * t307 + t333 * t197 + t334 * t198;
t298 = -rSges(7,2) * t254 + (t333 * t248 + t334 * t249) * t253;
t181 = -pkin(8) * t265 + t262 * t318;
t241 = t262 * pkin(1) - t265 * pkin(7);
t297 = -t181 - t241;
t296 = qJD(4) * t253;
t295 = qJD(5) * t253;
t294 = V_base(5) * pkin(6) + V_base(1);
t245 = qJD(2) * t262 + V_base(4);
t250 = V_base(6) + qJD(1);
t244 = -qJD(2) * t265 + V_base(5);
t291 = t244 * t319 + t294;
t220 = qJD(3) * t262 + t245;
t290 = pkin(3) * t254 + pkin(9) * t253;
t289 = rSges(3,1) * t264 - rSges(3,2) * t261;
t288 = rSges(4,1) * t254 - rSges(4,2) * t253;
t287 = Icges(3,1) * t264 - t313;
t286 = Icges(4,1) * t254 - t311;
t285 = -Icges(3,2) * t261 + t312;
t284 = -Icges(4,2) * t253 + t310;
t283 = Icges(3,5) * t264 - Icges(3,6) * t261;
t282 = Icges(4,5) * t254 - Icges(4,6) * t253;
t242 = t265 * pkin(1) + t262 * pkin(7);
t281 = -V_base(4) * pkin(6) + t250 * t242 + V_base(2);
t280 = V_base(4) * t241 - t242 * V_base(5) + V_base(3);
t219 = V_base(5) + (-qJD(2) - qJD(3)) * t265;
t279 = qJ(5) * t253 + t254 * t317;
t278 = (-Icges(4,3) * t265 + t262 * t282) * t219 + (Icges(4,3) * t262 + t265 * t282) * t220 + (Icges(4,5) * t253 + Icges(4,6) * t254) * t250;
t277 = (-Icges(3,3) * t265 + t262 * t283) * t244 + (Icges(3,3) * t262 + t265 * t283) * t245 + (Icges(3,5) * t261 + Icges(3,6) * t264) * t250;
t207 = t290 * t262;
t218 = pkin(3) * t253 - pkin(9) * t254;
t276 = t219 * t218 + (-t207 + t297) * t250 + t291;
t182 = pkin(8) * t262 + t265 * t318;
t275 = t245 * t181 - t182 * t244 + t280;
t274 = t250 * t182 - t245 * t319 + t281;
t162 = -qJ(5) * t254 + t253 * t317;
t193 = t262 * t296 + t219;
t273 = t193 * t162 + t265 * t295 + t276;
t208 = t290 * t265;
t272 = t220 * t207 - t208 * t219 + t275;
t271 = t250 * t208 - t218 * t220 + t274;
t150 = pkin(4) * t305 + t265 * t279;
t226 = -qJD(4) * t254 + t250;
t270 = t226 * t150 + t262 * t295 + t271;
t149 = -pkin(4) * t304 + t262 * t279;
t194 = t265 * t296 + t220;
t269 = -qJD(5) * t254 + t194 * t149 + t272;
t186 = -Icges(4,6) * t265 + t262 * t284;
t187 = Icges(4,6) * t262 + t265 * t284;
t188 = -Icges(4,5) * t265 + t262 * t286;
t189 = Icges(4,5) * t262 + t265 * t286;
t215 = Icges(4,2) * t254 + t311;
t216 = Icges(4,1) * t253 + t310;
t268 = (-t187 * t253 + t189 * t254) * t220 + (-t186 * t253 + t188 * t254) * t219 + (-t215 * t253 + t216 * t254) * t250;
t201 = -Icges(3,6) * t265 + t262 * t285;
t202 = Icges(3,6) * t262 + t265 * t285;
t203 = -Icges(3,5) * t265 + t262 * t287;
t204 = Icges(3,5) * t262 + t265 * t287;
t230 = Icges(3,2) * t264 + t313;
t233 = Icges(3,1) * t261 + t312;
t267 = (-t202 * t261 + t204 * t264) * t245 + (-t201 * t261 + t203 * t264) * t244 + (-t230 * t261 + t233 * t264) * t250;
t255 = Icges(2,4) * t265;
t238 = rSges(2,1) * t265 - rSges(2,2) * t262;
t237 = rSges(2,1) * t262 + rSges(2,2) * t265;
t236 = rSges(3,1) * t261 + rSges(3,2) * t264;
t235 = Icges(2,1) * t265 - t314;
t234 = Icges(2,1) * t262 + t255;
t232 = -Icges(2,2) * t262 + t255;
t231 = Icges(2,2) * t265 + t314;
t225 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t224 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t223 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t217 = rSges(4,1) * t253 + rSges(4,2) * t254;
t206 = rSges(3,3) * t262 + t265 * t289;
t205 = -rSges(3,3) * t265 + t262 * t289;
t192 = rSges(4,3) * t262 + t265 * t288;
t191 = -rSges(4,3) * t265 + t262 * t288;
t180 = -rSges(5,3) * t254 + (rSges(5,1) * t263 - rSges(5,2) * t260) * t253;
t179 = -Icges(5,5) * t254 + (Icges(5,1) * t263 - Icges(5,4) * t260) * t253;
t178 = -Icges(5,6) * t254 + (Icges(5,4) * t263 - Icges(5,2) * t260) * t253;
t176 = V_base(5) * rSges(2,3) - t237 * t250 + t294;
t175 = t238 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = t237 * V_base(4) - t238 * V_base(5) + V_base(3);
t171 = -rSges(6,3) * t254 + (rSges(6,1) * t249 - rSges(6,2) * t248) * t253;
t160 = rSges(5,1) * t212 + rSges(5,2) * t211 + rSges(5,3) * t307;
t159 = rSges(5,1) * t210 + rSges(5,2) * t209 + rSges(5,3) * t308;
t156 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t307;
t155 = Icges(5,1) * t210 + Icges(5,4) * t209 + Icges(5,5) * t308;
t154 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t307;
t153 = Icges(5,4) * t210 + Icges(5,2) * t209 + Icges(5,6) * t308;
t147 = rSges(6,1) * t198 - rSges(6,2) * t197 + rSges(6,3) * t307;
t145 = rSges(6,1) * t196 - rSges(6,2) * t195 + rSges(6,3) * t308;
t130 = t236 * t244 + (-t205 - t241) * t250 + t294;
t129 = t206 * t250 - t236 * t245 + t281;
t127 = t205 * t245 - t206 * t244 + t280;
t126 = t217 * t219 + (-t191 + t297) * t250 + t291;
t125 = t192 * t250 - t217 * t220 + t274;
t124 = t191 * t220 - t192 * t219 + t275;
t123 = -t159 * t226 + t180 * t193 + t276;
t122 = t160 * t226 - t180 * t194 + t271;
t121 = t159 * t194 - t160 * t193 + t272;
t120 = t171 * t193 + (-t145 - t149) * t226 + t273;
t119 = t147 * t226 + (-t162 - t171) * t194 + t270;
t118 = t145 * t194 + (-t147 - t150) * t193 + t269;
t117 = qJD(6) * t197 + t298 * t193 + (-t149 - t300) * t226 + t273;
t116 = qJD(6) * t195 + t299 * t226 + (-t162 - t298) * t194 + t270;
t115 = qJD(6) * t248 * t253 + t300 * t194 + (-t150 - t299) * t193 + t269;
t1 = t245 * (t277 * t262 + t267 * t265) / 0.2e1 + t244 * (t267 * t262 - t277 * t265) / 0.2e1 + t220 * (t278 * t262 + t268 * t265) / 0.2e1 + t219 * (t268 * t262 - t278 * t265) / 0.2e1 + m(1) * (t223 ^ 2 + t224 ^ 2 + t225 ^ 2) / 0.2e1 + m(2) * (t173 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(6) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(4) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(7) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + ((-t231 * t262 + t234 * t265 + Icges(1,4)) * V_base(5) + (-t262 * t232 + t265 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t265 * t231 + t262 * t234 + Icges(1,2)) * V_base(5) + (t232 * t265 + t235 * t262 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t178 * t209 + t210 * t179 + t195 * t328 + t196 * t327 + t308 * t324) * t226 + (t209 * t154 + t210 * t156 + t195 * t331 + t196 * t329 + t308 * t325) * t194 + (t209 * t153 + t210 * t155 + t332 * t195 + t330 * t196 + t326 * t308) * t193) * t193 / 0.2e1 + ((t211 * t178 + t212 * t179 + t197 * t328 + t198 * t327 + t307 * t324) * t226 + (t211 * t154 + t212 * t156 + t331 * t197 + t329 * t198 + t325 * t307) * t194 + (t211 * t153 + t155 * t212 + t197 * t332 + t330 * t198 + t326 * t307) * t193) * t194 / 0.2e1 + ((-t193 * t326 - t194 * t325 - t226 * t324) * t254 + ((-t178 * t260 + t179 * t263 + t248 * t328 + t249 * t327) * t226 + (-t154 * t260 + t156 * t263 + t248 * t331 + t249 * t329) * t194 + (-t153 * t260 + t155 * t263 + t248 * t332 + t330 * t249) * t193) * t253) * t226 / 0.2e1 + ((t264 * t202 + t261 * t204) * t245 + (t264 * t201 + t261 * t203) * t244 + (t254 * t187 + t253 * t189) * t220 + (t254 * t186 + t253 * t188) * t219 + (t254 * t215 + t253 * t216 + t264 * t230 + t261 * t233 + Icges(2,3)) * t250) * t250 / 0.2e1 + t250 * V_base(4) * (Icges(2,5) * t265 - Icges(2,6) * t262) + t250 * V_base(5) * (Icges(2,5) * t262 + Icges(2,6) * t265) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
