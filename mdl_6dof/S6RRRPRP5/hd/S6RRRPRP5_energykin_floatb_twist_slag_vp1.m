% Calculate kinetic energy for
% S6RRRPRP5
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:39
% EndTime: 2019-03-09 16:47:43
% DurationCPUTime: 3.86s
% Computational Cost: add. (2254->351), mult. (2723->504), div. (0->0), fcn. (2683->10), ass. (0->163)
t339 = Icges(6,1) + Icges(7,1);
t338 = -Icges(6,4) + Icges(7,5);
t337 = Icges(7,4) + Icges(6,5);
t336 = Icges(6,2) + Icges(7,3);
t335 = -Icges(7,6) + Icges(6,6);
t334 = -Icges(5,3) - Icges(4,3);
t333 = -Icges(6,3) - Icges(7,2);
t332 = rSges(7,1) + pkin(5);
t331 = rSges(7,3) + qJ(6);
t260 = qJ(3) + pkin(10);
t253 = qJ(5) + t260;
t248 = sin(t253);
t249 = cos(t253);
t267 = cos(qJ(1));
t264 = sin(qJ(1));
t266 = cos(qJ(2));
t301 = t264 * t266;
t186 = t248 * t301 + t249 * t267;
t187 = -t248 * t267 + t249 * t301;
t263 = sin(qJ(2));
t304 = t263 * t264;
t330 = t336 * t186 + t338 * t187 - t335 * t304;
t300 = t266 * t267;
t188 = t248 * t300 - t264 * t249;
t189 = t264 * t248 + t249 * t300;
t303 = t263 * t267;
t329 = t336 * t188 + t338 * t189 - t335 * t303;
t328 = -t335 * t186 + t337 * t187 - t333 * t304;
t327 = -t335 * t188 + t337 * t189 - t333 * t303;
t326 = t338 * t186 + t339 * t187 + t337 * t304;
t325 = t338 * t188 + t339 * t189 + t337 * t303;
t251 = sin(t260);
t252 = cos(t260);
t191 = -t251 * t301 - t252 * t267;
t192 = -t251 * t267 + t252 * t301;
t262 = sin(qJ(3));
t265 = cos(qJ(3));
t211 = -t262 * t301 - t265 * t267;
t305 = t262 * t267;
t212 = t265 * t301 - t305;
t324 = Icges(4,5) * t212 + Icges(5,5) * t192 + Icges(4,6) * t211 + Icges(5,6) * t191 - t334 * t304;
t193 = -t251 * t300 + t264 * t252;
t194 = t264 * t251 + t252 * t300;
t213 = -t262 * t300 + t264 * t265;
t302 = t264 * t262;
t214 = t265 * t300 + t302;
t323 = Icges(4,5) * t214 + Icges(5,5) * t194 + Icges(4,6) * t213 + Icges(5,6) * t193 - t334 * t303;
t322 = t335 * t266 + (t336 * t248 + t338 * t249) * t263;
t321 = t333 * t266 + (-t335 * t248 + t337 * t249) * t263;
t320 = -t337 * t266 + (t338 * t248 + t339 * t249) * t263;
t319 = t334 * t266 + (Icges(4,5) * t265 + Icges(5,5) * t252 - Icges(4,6) * t262 - Icges(5,6) * t251) * t263;
t310 = t265 * pkin(3);
t308 = Icges(2,4) * t264;
t307 = Icges(3,4) * t263;
t306 = Icges(3,4) * t266;
t299 = rSges(7,2) * t304 + t331 * t186 + t332 * t187;
t298 = rSges(7,2) * t303 + t331 * t188 + t332 * t189;
t297 = -rSges(7,2) * t266 + (t331 * t248 + t332 * t249) * t263;
t296 = pkin(4) * t252;
t294 = qJD(3) * t263;
t293 = qJD(4) * t263;
t292 = qJD(5) * t263;
t291 = V_base(5) * pkin(6) + V_base(1);
t243 = qJD(2) * t264 + V_base(4);
t254 = V_base(6) + qJD(1);
t288 = pkin(4) * t251;
t210 = t267 * t294 + t243;
t287 = pkin(2) * t266 + pkin(8) * t263;
t242 = -qJD(2) * t267 + V_base(5);
t286 = rSges(3,1) * t266 - rSges(3,2) * t263;
t285 = Icges(3,1) * t266 - t307;
t284 = -Icges(3,2) * t263 + t306;
t283 = Icges(3,5) * t266 - Icges(3,6) * t263;
t209 = t264 * t294 + t242;
t241 = pkin(1) * t267 + t264 * pkin(7);
t282 = -V_base(4) * pkin(6) + t254 * t241 + V_base(2);
t240 = t264 * pkin(1) - pkin(7) * t267;
t281 = V_base(4) * t240 - t241 * V_base(5) + V_base(3);
t280 = qJ(4) * t263 + t266 * t310;
t216 = t287 * t264;
t239 = pkin(2) * t263 - pkin(8) * t266;
t279 = t242 * t239 + (-t216 - t240) * t254 + t291;
t278 = (-Icges(3,3) * t267 + t264 * t283) * t242 + (Icges(3,3) * t264 + t267 * t283) * t243 + (Icges(3,5) * t263 + Icges(3,6) * t266) * t254;
t277 = pkin(9) * t263 + t266 * t296;
t217 = t287 * t267;
t276 = t254 * t217 - t239 * t243 + t282;
t178 = -qJ(4) * t266 + t263 * t310;
t275 = t209 * t178 + t267 * t293 + t279;
t274 = t243 * t216 - t217 * t242 + t281;
t162 = pkin(3) * t302 + t267 * t280;
t234 = -qJD(3) * t266 + t254;
t273 = t234 * t162 + t264 * t293 + t276;
t161 = -pkin(3) * t305 + t264 * t280;
t272 = -qJD(4) * t266 + t210 * t161 + t274;
t122 = t264 * t277 - t267 * t288;
t166 = -pkin(9) * t266 + t263 * t296;
t271 = t209 * t166 + (-t122 - t161) * t234 + t275;
t123 = t264 * t288 + t267 * t277;
t270 = t234 * t123 + (-t166 - t178) * t210 + t273;
t269 = t210 * t122 + (-t123 - t162) * t209 + t272;
t200 = -Icges(3,6) * t267 + t264 * t284;
t201 = Icges(3,6) * t264 + t267 * t284;
t203 = -Icges(3,5) * t267 + t264 * t285;
t204 = Icges(3,5) * t264 + t267 * t285;
t228 = Icges(3,2) * t266 + t307;
t231 = Icges(3,1) * t263 + t306;
t268 = (-t201 * t263 + t204 * t266) * t243 + (-t200 * t263 + t203 * t266) * t242 + (-t228 * t263 + t231 * t266) * t254;
t256 = Icges(2,4) * t267;
t237 = rSges(2,1) * t267 - t264 * rSges(2,2);
t236 = t264 * rSges(2,1) + rSges(2,2) * t267;
t235 = rSges(3,1) * t263 + rSges(3,2) * t266;
t233 = Icges(2,1) * t267 - t308;
t232 = Icges(2,1) * t264 + t256;
t230 = -Icges(2,2) * t264 + t256;
t229 = Icges(2,2) * t267 + t308;
t223 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t222 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t221 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t218 = (-qJD(3) - qJD(5)) * t266 + t254;
t207 = t264 * rSges(3,3) + t267 * t286;
t206 = -rSges(3,3) * t267 + t264 * t286;
t205 = -rSges(4,3) * t266 + (rSges(4,1) * t265 - rSges(4,2) * t262) * t263;
t202 = -Icges(4,5) * t266 + (Icges(4,1) * t265 - Icges(4,4) * t262) * t263;
t199 = -Icges(4,6) * t266 + (Icges(4,4) * t265 - Icges(4,2) * t262) * t263;
t185 = t267 * t292 + t210;
t184 = t264 * t292 + t209;
t182 = -rSges(5,3) * t266 + (rSges(5,1) * t252 - rSges(5,2) * t251) * t263;
t181 = -Icges(5,5) * t266 + (Icges(5,1) * t252 - Icges(5,4) * t251) * t263;
t180 = -Icges(5,6) * t266 + (Icges(5,4) * t252 - Icges(5,2) * t251) * t263;
t177 = -rSges(6,3) * t266 + (rSges(6,1) * t249 - rSges(6,2) * t248) * t263;
t175 = V_base(5) * rSges(2,3) - t236 * t254 + t291;
t174 = t237 * t254 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t167 = t236 * V_base(4) - t237 * V_base(5) + V_base(3);
t164 = t214 * rSges(4,1) + t213 * rSges(4,2) + rSges(4,3) * t303;
t163 = rSges(4,1) * t212 + rSges(4,2) * t211 + rSges(4,3) * t304;
t160 = Icges(4,1) * t214 + Icges(4,4) * t213 + Icges(4,5) * t303;
t159 = Icges(4,1) * t212 + Icges(4,4) * t211 + Icges(4,5) * t304;
t158 = Icges(4,4) * t214 + Icges(4,2) * t213 + Icges(4,6) * t303;
t157 = Icges(4,4) * t212 + Icges(4,2) * t211 + Icges(4,6) * t304;
t152 = t194 * rSges(5,1) + t193 * rSges(5,2) + rSges(5,3) * t303;
t151 = rSges(5,1) * t192 + rSges(5,2) * t191 + rSges(5,3) * t304;
t150 = Icges(5,1) * t194 + Icges(5,4) * t193 + Icges(5,5) * t303;
t149 = Icges(5,1) * t192 + Icges(5,4) * t191 + Icges(5,5) * t304;
t148 = Icges(5,4) * t194 + Icges(5,2) * t193 + Icges(5,6) * t303;
t147 = Icges(5,4) * t192 + Icges(5,2) * t191 + Icges(5,6) * t304;
t142 = t189 * rSges(6,1) - t188 * rSges(6,2) + rSges(6,3) * t303;
t140 = rSges(6,1) * t187 - rSges(6,2) * t186 + rSges(6,3) * t304;
t125 = t235 * t242 + (-t206 - t240) * t254 + t291;
t124 = t207 * t254 - t235 * t243 + t282;
t120 = t206 * t243 - t207 * t242 + t281;
t118 = -t163 * t234 + t205 * t209 + t279;
t117 = t164 * t234 - t205 * t210 + t276;
t116 = t163 * t210 - t164 * t209 + t274;
t115 = t182 * t209 + (-t151 - t161) * t234 + t275;
t114 = t152 * t234 + (-t178 - t182) * t210 + t273;
t113 = t151 * t210 + (-t152 - t162) * t209 + t272;
t112 = -t140 * t218 + t177 * t184 + t271;
t111 = t142 * t218 - t177 * t185 + t270;
t110 = t140 * t185 - t142 * t184 + t269;
t109 = qJD(6) * t188 + t184 * t297 - t218 * t299 + t271;
t108 = qJD(6) * t186 - t185 * t297 + t218 * t298 + t270;
t107 = qJD(6) * t248 * t263 - t184 * t298 + t185 * t299 + t269;
t1 = m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(2) * (t167 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t243 * (t264 * t278 + t267 * t268) / 0.2e1 + t242 * (t264 * t268 - t278 * t267) / 0.2e1 + ((t186 * t322 + t187 * t320 + t304 * t321) * t218 + (t186 * t329 + t187 * t325 + t304 * t327) * t185 + (t330 * t186 + t326 * t187 + t328 * t304) * t184) * t184 / 0.2e1 + ((t188 * t322 + t189 * t320 + t303 * t321) * t218 + (t329 * t188 + t325 * t189 + t327 * t303) * t185 + (t188 * t330 + t326 * t189 + t328 * t303) * t184) * t185 / 0.2e1 + ((t180 * t191 + t181 * t192 + t199 * t211 + t202 * t212 + t304 * t319) * t234 + (t148 * t191 + t150 * t192 + t158 * t211 + t160 * t212 + t304 * t323) * t210 + (t191 * t147 + t192 * t149 + t211 * t157 + t212 * t159 + t324 * t304) * t209) * t209 / 0.2e1 + ((t193 * t180 + t194 * t181 + t213 * t199 + t214 * t202 + t303 * t319) * t234 + (t193 * t148 + t194 * t150 + t213 * t158 + t214 * t160 + t323 * t303) * t210 + (t193 * t147 + t194 * t149 + t213 * t157 + t214 * t159 + t303 * t324) * t209) * t210 / 0.2e1 + ((-t184 * t328 - t185 * t327 - t218 * t321) * t266 + ((t248 * t322 + t249 * t320) * t218 + (t248 * t329 + t249 * t325) * t185 + (t248 * t330 + t326 * t249) * t184) * t263) * t218 / 0.2e1 + ((-t209 * t324 - t210 * t323 - t234 * t319) * t266 + ((-t180 * t251 + t181 * t252 - t199 * t262 + t202 * t265) * t234 + (-t148 * t251 + t150 * t252 - t158 * t262 + t160 * t265) * t210 + (-t147 * t251 + t149 * t252 - t157 * t262 + t159 * t265) * t209) * t263) * t234 / 0.2e1 + ((t201 * t266 + t204 * t263) * t243 + (t200 * t266 + t203 * t263) * t242 + (t266 * t228 + t263 * t231 + Icges(2,3)) * t254) * t254 / 0.2e1 + ((-t264 * t229 + t232 * t267 + Icges(1,4)) * V_base(5) + (-t264 * t230 + t267 * t233 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t229 + t264 * t232 + Icges(1,2)) * V_base(5) + (t230 * t267 + t264 * t233 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t254 * (Icges(2,5) * t267 - Icges(2,6) * t264) + V_base(5) * t254 * (Icges(2,5) * t264 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
