% Calculate kinetic energy for
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:03
% EndTime: 2019-03-10 01:18:07
% DurationCPUTime: 4.10s
% Computational Cost: add. (2389->362), mult. (2881->545), div. (0->0), fcn. (2831->10), ass. (0->171)
t337 = Icges(6,1) + Icges(7,1);
t336 = Icges(6,4) + Icges(7,4);
t335 = -Icges(7,5) - Icges(6,5);
t334 = Icges(6,2) + Icges(7,2);
t333 = -Icges(7,6) - Icges(6,6);
t332 = -Icges(7,3) - Icges(6,3);
t262 = qJ(3) + qJ(4);
t257 = qJ(5) + t262;
t248 = sin(t257);
t249 = cos(t257);
t268 = cos(qJ(1));
t265 = sin(qJ(1));
t267 = cos(qJ(2));
t309 = t265 * t267;
t183 = -t248 * t309 - t249 * t268;
t184 = -t248 * t268 + t249 * t309;
t264 = sin(qJ(2));
t311 = t264 * t265;
t331 = -t183 * t333 - t184 * t335 - t311 * t332;
t308 = t267 * t268;
t185 = -t248 * t308 + t249 * t265;
t186 = t248 * t265 + t249 * t308;
t310 = t264 * t268;
t330 = -t185 * t333 - t186 * t335 - t310 * t332;
t329 = t183 * t334 + t184 * t336 - t311 * t333;
t328 = t185 * t334 + t186 * t336 - t310 * t333;
t327 = t336 * t183 + t184 * t337 - t335 * t311;
t326 = t336 * t185 + t186 * t337 - t335 * t310;
t325 = t332 * t267 + (t248 * t333 - t249 * t335) * t264;
t324 = t333 * t267 + (-t248 * t334 + t249 * t336) * t264;
t323 = t335 * t267 + (-t336 * t248 + t249 * t337) * t264;
t266 = cos(qJ(3));
t318 = t266 * pkin(3);
t316 = Icges(2,4) * t265;
t315 = Icges(3,4) * t264;
t314 = Icges(3,4) * t267;
t263 = sin(qJ(3));
t313 = t263 * t265;
t312 = t263 * t268;
t304 = pkin(5) * t249;
t280 = qJ(6) * t264 + t267 * t304;
t303 = pkin(5) * t248;
t307 = rSges(7,1) * t184 + rSges(7,2) * t183 + rSges(7,3) * t311 + t265 * t280 - t268 * t303;
t306 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t310 + t265 * t303 + t268 * t280;
t305 = (-qJ(6) - rSges(7,3)) * t267 + (rSges(7,1) * t249 - rSges(7,2) * t248 + t304) * t264;
t254 = cos(t262);
t302 = pkin(4) * t254;
t299 = qJD(3) * t264;
t298 = qJD(4) * t264;
t297 = qJD(5) * t264;
t296 = qJD(6) * t264;
t295 = -qJD(3) - qJD(4);
t294 = V_base(5) * pkin(6) + V_base(1);
t241 = qJD(2) * t265 + V_base(4);
t251 = V_base(6) + qJD(1);
t253 = sin(t262);
t291 = pkin(4) * t253;
t209 = t268 * t299 + t241;
t290 = pkin(2) * t267 + pkin(8) * t264;
t240 = -qJD(2) * t268 + V_base(5);
t289 = rSges(3,1) * t267 - rSges(3,2) * t264;
t182 = t268 * t298 + t209;
t288 = Icges(3,1) * t267 - t315;
t287 = -Icges(3,2) * t264 + t314;
t286 = Icges(3,5) * t267 - Icges(3,6) * t264;
t208 = t265 * t299 + t240;
t239 = pkin(1) * t268 + pkin(7) * t265;
t285 = -V_base(4) * pkin(6) + t251 * t239 + V_base(2);
t238 = pkin(1) * t265 - pkin(7) * t268;
t284 = V_base(4) * t238 - t239 * V_base(5) + V_base(3);
t181 = t265 * t298 + t208;
t283 = pkin(9) * t264 + t267 * t318;
t215 = t290 * t265;
t237 = t264 * pkin(2) - t267 * pkin(8);
t282 = t240 * t237 + (-t215 - t238) * t251 + t294;
t281 = (-Icges(3,3) * t268 + t265 * t286) * t240 + (Icges(3,3) * t265 + t268 * t286) * t241 + (Icges(3,5) * t264 + Icges(3,6) * t267) * t251;
t279 = pkin(10) * t264 + t267 * t302;
t216 = t290 * t268;
t278 = t251 * t216 - t237 * t241 + t285;
t277 = t241 * t215 - t216 * t240 + t284;
t158 = -pkin(3) * t312 + t265 * t283;
t175 = -pkin(9) * t267 + t264 * t318;
t233 = -qJD(3) * t267 + t251;
t276 = -t158 * t233 + t208 * t175 + t282;
t159 = pkin(3) * t313 + t268 * t283;
t275 = t233 * t159 - t175 * t209 + t278;
t274 = t209 * t158 - t159 * t208 + t277;
t118 = t265 * t279 - t268 * t291;
t161 = -pkin(10) * t267 + t264 * t302;
t217 = t267 * t295 + t251;
t273 = -t118 * t217 + t181 * t161 + t276;
t119 = t265 * t291 + t268 * t279;
t272 = t217 * t119 - t161 * t182 + t275;
t271 = t182 * t118 - t119 * t181 + t274;
t192 = -Icges(3,6) * t268 + t265 * t287;
t193 = Icges(3,6) * t265 + t268 * t287;
t195 = -Icges(3,5) * t268 + t265 * t288;
t196 = Icges(3,5) * t265 + t268 * t288;
t227 = Icges(3,2) * t267 + t315;
t230 = Icges(3,1) * t264 + t314;
t270 = (-t193 * t264 + t196 * t267) * t241 + (-t192 * t264 + t195 * t267) * t240 + (-t227 * t264 + t230 * t267) * t251;
t256 = Icges(2,4) * t268;
t236 = rSges(2,1) * t268 - rSges(2,2) * t265;
t235 = rSges(2,1) * t265 + rSges(2,2) * t268;
t234 = rSges(3,1) * t264 + rSges(3,2) * t267;
t232 = Icges(2,1) * t268 - t316;
t231 = Icges(2,1) * t265 + t256;
t229 = -Icges(2,2) * t265 + t256;
t228 = Icges(2,2) * t268 + t316;
t222 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t221 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t220 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t213 = t266 * t308 + t313;
t212 = -t263 * t308 + t265 * t266;
t211 = t266 * t309 - t312;
t210 = -t263 * t309 - t266 * t268;
t205 = (-qJD(5) + t295) * t267 + t251;
t204 = t253 * t265 + t254 * t308;
t203 = -t253 * t308 + t254 * t265;
t202 = -t253 * t268 + t254 * t309;
t201 = -t253 * t309 - t254 * t268;
t200 = rSges(3,3) * t265 + t268 * t289;
t199 = -rSges(3,3) * t268 + t265 * t289;
t198 = -rSges(4,3) * t267 + (rSges(4,1) * t266 - rSges(4,2) * t263) * t264;
t194 = -Icges(4,5) * t267 + (Icges(4,1) * t266 - Icges(4,4) * t263) * t264;
t191 = -Icges(4,6) * t267 + (Icges(4,4) * t266 - Icges(4,2) * t263) * t264;
t188 = -Icges(4,3) * t267 + (Icges(4,5) * t266 - Icges(4,6) * t263) * t264;
t180 = -rSges(5,3) * t267 + (rSges(5,1) * t254 - rSges(5,2) * t253) * t264;
t178 = -Icges(5,5) * t267 + (Icges(5,1) * t254 - Icges(5,4) * t253) * t264;
t177 = -Icges(5,6) * t267 + (Icges(5,4) * t254 - Icges(5,2) * t253) * t264;
t176 = -Icges(5,3) * t267 + (Icges(5,5) * t254 - Icges(5,6) * t253) * t264;
t174 = -rSges(6,3) * t267 + (rSges(6,1) * t249 - rSges(6,2) * t248) * t264;
t166 = V_base(5) * rSges(2,3) - t235 * t251 + t294;
t165 = t236 * t251 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t164 = t235 * V_base(4) - t236 * V_base(5) + V_base(3);
t163 = t268 * t297 + t182;
t162 = t265 * t297 + t181;
t157 = rSges(4,1) * t213 + rSges(4,2) * t212 + rSges(4,3) * t310;
t156 = rSges(4,1) * t211 + rSges(4,2) * t210 + rSges(4,3) * t311;
t155 = Icges(4,1) * t213 + Icges(4,4) * t212 + Icges(4,5) * t310;
t154 = Icges(4,1) * t211 + Icges(4,4) * t210 + Icges(4,5) * t311;
t153 = Icges(4,4) * t213 + Icges(4,2) * t212 + Icges(4,6) * t310;
t152 = Icges(4,4) * t211 + Icges(4,2) * t210 + Icges(4,6) * t311;
t151 = Icges(4,5) * t213 + Icges(4,6) * t212 + Icges(4,3) * t310;
t150 = Icges(4,5) * t211 + Icges(4,6) * t210 + Icges(4,3) * t311;
t149 = rSges(5,1) * t204 + rSges(5,2) * t203 + rSges(5,3) * t310;
t148 = rSges(5,1) * t202 + rSges(5,2) * t201 + rSges(5,3) * t311;
t147 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t310;
t146 = Icges(5,1) * t202 + Icges(5,4) * t201 + Icges(5,5) * t311;
t145 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t310;
t144 = Icges(5,4) * t202 + Icges(5,2) * t201 + Icges(5,6) * t311;
t143 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t310;
t142 = Icges(5,5) * t202 + Icges(5,6) * t201 + Icges(5,3) * t311;
t139 = rSges(6,1) * t186 + rSges(6,2) * t185 + rSges(6,3) * t310;
t137 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t311;
t121 = t234 * t240 + (-t199 - t238) * t251 + t294;
t120 = t200 * t251 - t234 * t241 + t285;
t117 = t199 * t241 - t200 * t240 + t284;
t112 = -t156 * t233 + t198 * t208 + t282;
t111 = t157 * t233 - t198 * t209 + t278;
t110 = t156 * t209 - t157 * t208 + t277;
t109 = -t148 * t217 + t180 * t181 + t276;
t108 = t149 * t217 - t180 * t182 + t275;
t107 = t148 * t182 - t149 * t181 + t274;
t106 = -t137 * t205 + t162 * t174 + t273;
t105 = t139 * t205 - t163 * t174 + t272;
t104 = t137 * t163 - t139 * t162 + t271;
t103 = t162 * t305 - t205 * t307 + t268 * t296 + t273;
t102 = -t163 * t305 + t205 * t306 + t265 * t296 + t272;
t101 = -qJD(6) * t267 - t162 * t306 + t163 * t307 + t271;
t1 = t209 * ((t151 * t310 + t212 * t153 + t213 * t155) * t209 + (t150 * t310 + t152 * t212 + t154 * t213) * t208 + (t188 * t310 + t191 * t212 + t194 * t213) * t233) / 0.2e1 + t182 * ((t143 * t310 + t203 * t145 + t204 * t147) * t182 + (t142 * t310 + t144 * t203 + t146 * t204) * t181 + (t176 * t310 + t177 * t203 + t178 * t204) * t217) / 0.2e1 + t208 * ((t151 * t311 + t153 * t210 + t155 * t211) * t209 + (t150 * t311 + t210 * t152 + t211 * t154) * t208 + (t188 * t311 + t191 * t210 + t194 * t211) * t233) / 0.2e1 + t181 * ((t143 * t311 + t145 * t201 + t147 * t202) * t182 + (t142 * t311 + t201 * t144 + t202 * t146) * t181 + (t176 * t311 + t177 * t201 + t178 * t202) * t217) / 0.2e1 + t241 * (t281 * t265 + t270 * t268) / 0.2e1 + t240 * (t270 * t265 - t281 * t268) / 0.2e1 + m(1) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + m(2) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(5) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + t217 * ((-t142 * t181 - t143 * t182 - t176 * t217) * t267 + ((-t145 * t253 + t147 * t254) * t182 + (-t144 * t253 + t146 * t254) * t181 + (-t177 * t253 + t178 * t254) * t217) * t264) / 0.2e1 + t233 * ((-t150 * t208 - t151 * t209 - t188 * t233) * t267 + ((-t153 * t263 + t155 * t266) * t209 + (-t152 * t263 + t154 * t266) * t208 + (-t191 * t263 + t194 * t266) * t233) * t264) / 0.2e1 + ((t183 * t324 + t184 * t323 + t325 * t311) * t205 + (t183 * t328 + t184 * t326 + t311 * t330) * t163 + (t329 * t183 + t327 * t184 + t331 * t311) * t162) * t162 / 0.2e1 + ((t185 * t324 + t186 * t323 + t310 * t325) * t205 + (t328 * t185 + t326 * t186 + t330 * t310) * t163 + (t329 * t185 + t327 * t186 + t310 * t331) * t162) * t163 / 0.2e1 + ((-t162 * t331 - t330 * t163 - t325 * t205) * t267 + ((-t248 * t324 + t249 * t323) * t205 + (-t248 * t328 + t249 * t326) * t163 + (-t248 * t329 + t249 * t327) * t162) * t264) * t205 / 0.2e1 + ((t193 * t267 + t196 * t264) * t241 + (t192 * t267 + t195 * t264) * t240 + (t267 * t227 + t264 * t230 + Icges(2,3)) * t251) * t251 / 0.2e1 + ((-t228 * t265 + t231 * t268 + Icges(1,4)) * V_base(5) + (-t265 * t229 + t268 * t232 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t268 * t228 + t265 * t231 + Icges(1,2)) * V_base(5) + (t229 * t268 + t232 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t251 * (Icges(2,5) * t268 - Icges(2,6) * t265) + V_base(5) * t251 * (Icges(2,5) * t265 + Icges(2,6) * t268) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
