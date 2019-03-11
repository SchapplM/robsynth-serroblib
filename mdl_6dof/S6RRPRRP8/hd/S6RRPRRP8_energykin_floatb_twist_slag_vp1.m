% Calculate kinetic energy for
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:21:52
% EndTime: 2019-03-09 12:21:56
% DurationCPUTime: 3.96s
% Computational Cost: add. (2221->359), mult. (2668->526), div. (0->0), fcn. (2628->10), ass. (0->166)
t335 = Icges(6,1) + Icges(7,1);
t334 = -Icges(6,4) + Icges(7,5);
t333 = Icges(7,4) + Icges(6,5);
t332 = Icges(6,2) + Icges(7,3);
t331 = -Icges(7,6) + Icges(6,6);
t330 = -Icges(6,3) - Icges(7,2);
t329 = rSges(7,1) + pkin(5);
t328 = rSges(7,3) + qJ(6);
t260 = pkin(10) + qJ(4);
t253 = qJ(5) + t260;
t248 = sin(t253);
t249 = cos(t253);
t267 = cos(qJ(1));
t265 = sin(qJ(1));
t266 = cos(qJ(2));
t303 = t265 * t266;
t186 = t248 * t303 + t249 * t267;
t187 = -t248 * t267 + t249 * t303;
t264 = sin(qJ(2));
t306 = t264 * t265;
t327 = t332 * t186 + t334 * t187 - t331 * t306;
t302 = t266 * t267;
t188 = t248 * t302 - t265 * t249;
t189 = t265 * t248 + t249 * t302;
t305 = t264 * t267;
t326 = t332 * t188 + t334 * t189 - t331 * t305;
t325 = -t331 * t186 + t333 * t187 - t330 * t306;
t324 = -t331 * t188 + t333 * t189 - t330 * t305;
t323 = t334 * t186 + t335 * t187 + t333 * t306;
t322 = t334 * t188 + t335 * t189 + t333 * t305;
t321 = t331 * t266 + (t332 * t248 + t334 * t249) * t264;
t320 = t330 * t266 + (-t331 * t248 + t333 * t249) * t264;
t319 = -t333 * t266 + (t334 * t248 + t335 * t249) * t264;
t262 = cos(pkin(10));
t311 = t262 * pkin(3);
t310 = Icges(2,4) * t265;
t309 = Icges(3,4) * t264;
t308 = Icges(3,4) * t266;
t261 = sin(pkin(10));
t307 = t261 * t267;
t304 = t265 * t261;
t300 = rSges(7,2) * t306 + t328 * t186 + t329 * t187;
t299 = rSges(7,2) * t305 + t328 * t188 + t329 * t189;
t298 = -rSges(7,2) * t266 + (t328 * t248 + t329 * t249) * t264;
t285 = pkin(2) * t266 + qJ(3) * t264;
t215 = t285 * t265;
t240 = t265 * pkin(1) - pkin(7) * t267;
t297 = -t215 - t240;
t252 = cos(t260);
t296 = pkin(4) * t252;
t294 = qJD(3) * t264;
t293 = qJD(4) * t264;
t292 = qJD(5) * t264;
t291 = V_base(5) * pkin(6) + V_base(1);
t243 = qJD(2) * t265 + V_base(4);
t254 = V_base(6) + qJD(1);
t251 = sin(t260);
t288 = pkin(4) * t251;
t214 = t267 * t293 + t243;
t235 = pkin(2) * t264 - qJ(3) * t266;
t242 = -qJD(2) * t267 + V_base(5);
t287 = t242 * t235 + t267 * t294 + t291;
t286 = rSges(3,1) * t266 - rSges(3,2) * t264;
t284 = Icges(3,1) * t266 - t309;
t283 = -Icges(3,2) * t264 + t308;
t282 = Icges(3,5) * t266 - Icges(3,6) * t264;
t213 = t265 * t293 + t242;
t241 = pkin(1) * t267 + t265 * pkin(7);
t281 = -V_base(4) * pkin(6) + t254 * t241 + V_base(2);
t280 = V_base(4) * t240 - t241 * V_base(5) + V_base(3);
t279 = (-Icges(3,3) * t267 + t265 * t282) * t242 + (Icges(3,3) * t265 + t267 * t282) * t243 + (Icges(3,5) * t264 + Icges(3,6) * t266) * t254;
t278 = pkin(8) * t264 + t266 * t311;
t216 = t285 * t267;
t277 = t254 * t216 + t265 * t294 + t281;
t276 = pkin(9) * t264 + t266 * t296;
t275 = -qJD(3) * t266 + t243 * t215 + t280;
t163 = -pkin(3) * t307 + t265 * t278;
t178 = -pkin(8) * t266 + t264 * t311;
t274 = t242 * t178 + (-t163 + t297) * t254 + t287;
t164 = pkin(3) * t304 + t267 * t278;
t273 = t254 * t164 + (-t178 - t235) * t243 + t277;
t122 = t265 * t276 - t267 * t288;
t165 = -pkin(9) * t266 + t264 * t296;
t234 = -qJD(4) * t266 + t254;
t272 = -t122 * t234 + t213 * t165 + t274;
t271 = t243 * t163 + (-t164 - t216) * t242 + t275;
t123 = t265 * t288 + t267 * t276;
t270 = t234 * t123 - t165 * t214 + t273;
t269 = t214 * t122 - t123 * t213 + t271;
t202 = -Icges(3,6) * t267 + t265 * t283;
t203 = Icges(3,6) * t265 + t267 * t283;
t204 = -Icges(3,5) * t267 + t265 * t284;
t205 = Icges(3,5) * t265 + t267 * t284;
t228 = Icges(3,2) * t266 + t309;
t231 = Icges(3,1) * t264 + t308;
t268 = (-t203 * t264 + t205 * t266) * t243 + (-t202 * t264 + t204 * t266) * t242 + (-t228 * t264 + t231 * t266) * t254;
t257 = Icges(2,4) * t267;
t238 = rSges(2,1) * t267 - t265 * rSges(2,2);
t237 = t265 * rSges(2,1) + rSges(2,2) * t267;
t236 = rSges(3,1) * t264 + rSges(3,2) * t266;
t233 = Icges(2,1) * t267 - t310;
t232 = Icges(2,1) * t265 + t257;
t230 = -Icges(2,2) * t265 + t257;
t229 = Icges(2,2) * t267 + t310;
t224 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t223 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t222 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t218 = (-qJD(4) - qJD(5)) * t266 + t254;
t212 = t262 * t302 + t304;
t211 = -t261 * t302 + t265 * t262;
t210 = t262 * t303 - t307;
t209 = -t261 * t303 - t262 * t267;
t207 = t265 * rSges(3,3) + t267 * t286;
t206 = -rSges(3,3) * t267 + t265 * t286;
t199 = t265 * t251 + t252 * t302;
t198 = -t251 * t302 + t265 * t252;
t197 = -t251 * t267 + t252 * t303;
t196 = -t251 * t303 - t252 * t267;
t195 = -rSges(4,3) * t266 + (rSges(4,1) * t262 - rSges(4,2) * t261) * t264;
t193 = -Icges(4,5) * t266 + (Icges(4,1) * t262 - Icges(4,4) * t261) * t264;
t192 = -Icges(4,6) * t266 + (Icges(4,4) * t262 - Icges(4,2) * t261) * t264;
t191 = -Icges(4,3) * t266 + (Icges(4,5) * t262 - Icges(4,6) * t261) * t264;
t185 = t267 * t292 + t214;
t184 = t265 * t292 + t213;
t182 = -rSges(5,3) * t266 + (rSges(5,1) * t252 - rSges(5,2) * t251) * t264;
t181 = -Icges(5,5) * t266 + (Icges(5,1) * t252 - Icges(5,4) * t251) * t264;
t180 = -Icges(5,6) * t266 + (Icges(5,4) * t252 - Icges(5,2) * t251) * t264;
t179 = -Icges(5,3) * t266 + (Icges(5,5) * t252 - Icges(5,6) * t251) * t264;
t177 = -rSges(6,3) * t266 + (rSges(6,1) * t249 - rSges(6,2) * t248) * t264;
t175 = V_base(5) * rSges(2,3) - t237 * t254 + t291;
t174 = t238 * t254 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t167 = t237 * V_base(4) - t238 * V_base(5) + V_base(3);
t162 = t212 * rSges(4,1) + t211 * rSges(4,2) + rSges(4,3) * t305;
t161 = rSges(4,1) * t210 + rSges(4,2) * t209 + rSges(4,3) * t306;
t160 = Icges(4,1) * t212 + Icges(4,4) * t211 + Icges(4,5) * t305;
t159 = Icges(4,1) * t210 + Icges(4,4) * t209 + Icges(4,5) * t306;
t158 = Icges(4,4) * t212 + Icges(4,2) * t211 + Icges(4,6) * t305;
t157 = Icges(4,4) * t210 + Icges(4,2) * t209 + Icges(4,6) * t306;
t156 = Icges(4,5) * t212 + Icges(4,6) * t211 + Icges(4,3) * t305;
t155 = Icges(4,5) * t210 + Icges(4,6) * t209 + Icges(4,3) * t306;
t151 = t199 * rSges(5,1) + t198 * rSges(5,2) + rSges(5,3) * t305;
t150 = rSges(5,1) * t197 + rSges(5,2) * t196 + rSges(5,3) * t306;
t149 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t305;
t148 = Icges(5,1) * t197 + Icges(5,4) * t196 + Icges(5,5) * t306;
t147 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t305;
t146 = Icges(5,4) * t197 + Icges(5,2) * t196 + Icges(5,6) * t306;
t145 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t305;
t144 = Icges(5,5) * t197 + Icges(5,6) * t196 + Icges(5,3) * t306;
t141 = t189 * rSges(6,1) - t188 * rSges(6,2) + rSges(6,3) * t305;
t139 = rSges(6,1) * t187 - rSges(6,2) * t186 + rSges(6,3) * t306;
t125 = t236 * t242 + (-t206 - t240) * t254 + t291;
t124 = t207 * t254 - t236 * t243 + t281;
t120 = t206 * t243 - t207 * t242 + t280;
t118 = t195 * t242 + (-t161 + t297) * t254 + t287;
t117 = t162 * t254 + (-t195 - t235) * t243 + t277;
t116 = t161 * t243 + (-t162 - t216) * t242 + t275;
t115 = -t150 * t234 + t182 * t213 + t274;
t114 = t151 * t234 - t182 * t214 + t273;
t113 = t150 * t214 - t151 * t213 + t271;
t112 = -t139 * t218 + t177 * t184 + t272;
t111 = t141 * t218 - t177 * t185 + t270;
t110 = t139 * t185 - t141 * t184 + t269;
t109 = qJD(6) * t188 + t184 * t298 - t218 * t300 + t272;
t108 = qJD(6) * t186 - t185 * t298 + t218 * t299 + t270;
t107 = qJD(6) * t248 * t264 - t184 * t299 + t185 * t300 + t269;
t1 = t214 * ((t145 * t305 + t198 * t147 + t199 * t149) * t214 + (t144 * t305 + t198 * t146 + t199 * t148) * t213 + (t179 * t305 + t198 * t180 + t199 * t181) * t234) / 0.2e1 + t213 * ((t145 * t306 + t147 * t196 + t149 * t197) * t214 + (t144 * t306 + t196 * t146 + t197 * t148) * t213 + (t179 * t306 + t180 * t196 + t181 * t197) * t234) / 0.2e1 + m(1) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(2) * (t167 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + t234 * ((-t144 * t213 - t145 * t214 - t179 * t234) * t266 + ((-t147 * t251 + t149 * t252) * t214 + (-t146 * t251 + t148 * t252) * t213 + (-t180 * t251 + t181 * t252) * t234) * t264) / 0.2e1 + ((t186 * t321 + t187 * t319 + t306 * t320) * t218 + (t186 * t326 + t187 * t322 + t306 * t324) * t185 + (t327 * t186 + t323 * t187 + t325 * t306) * t184) * t184 / 0.2e1 + ((t188 * t321 + t189 * t319 + t305 * t320) * t218 + (t326 * t188 + t322 * t189 + t324 * t305) * t185 + (t188 * t327 + t323 * t189 + t325 * t305) * t184) * t185 / 0.2e1 + ((-t184 * t325 - t185 * t324 - t218 * t320) * t266 + ((t248 * t321 + t249 * t319) * t218 + (t248 * t326 + t249 * t322) * t185 + (t248 * t327 + t323 * t249) * t184) * t264) * t218 / 0.2e1 + ((t156 * t306 + t158 * t209 + t160 * t210) * t243 + (t155 * t306 + t209 * t157 + t210 * t159) * t242 + (t191 * t306 + t192 * t209 + t193 * t210) * t254 + t265 * t268 - t279 * t267) * t242 / 0.2e1 + ((t156 * t305 + t211 * t158 + t212 * t160) * t243 + (t155 * t305 + t211 * t157 + t212 * t159) * t242 + (t191 * t305 + t211 * t192 + t212 * t193) * t254 + t265 * t279 + t267 * t268) * t243 / 0.2e1 + ((-t265 * t229 + t232 * t267 + Icges(1,4)) * V_base(5) + (-t265 * t230 + t267 * t233 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t229 + t265 * t232 + Icges(1,2)) * V_base(5) + (t230 * t267 + t265 * t233 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t155 * t242 - t156 * t243) * t266 + ((-t158 * t261 + t160 * t262) * t243 + (-t157 * t261 + t159 * t262) * t242) * t264 + (t203 * t266 + t205 * t264) * t243 + (t202 * t266 + t204 * t264) * t242 + (Icges(2,3) + (-t191 + t228) * t266 + (-t192 * t261 + t193 * t262 + t231) * t264) * t254) * t254 / 0.2e1 + t254 * V_base(4) * (Icges(2,5) * t267 - Icges(2,6) * t265) + t254 * V_base(5) * (Icges(2,5) * t265 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
