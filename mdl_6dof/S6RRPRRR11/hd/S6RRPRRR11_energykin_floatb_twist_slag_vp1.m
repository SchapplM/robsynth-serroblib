% Calculate kinetic energy for
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:07
% EndTime: 2019-03-09 14:29:12
% DurationCPUTime: 5.30s
% Computational Cost: add. (1705->362), mult. (2352->544), div. (0->0), fcn. (2246->10), ass. (0->173)
t342 = Icges(3,4) + Icges(4,6);
t341 = Icges(3,1) + Icges(4,2);
t340 = -Icges(3,2) - Icges(4,3);
t266 = cos(qJ(2));
t339 = t342 * t266;
t263 = sin(qJ(2));
t338 = t342 * t263;
t337 = Icges(4,4) - Icges(3,5);
t336 = Icges(4,5) - Icges(3,6);
t335 = t263 * t340 + t339;
t334 = t266 * t341 - t338;
t333 = Icges(4,1) + Icges(3,3);
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t332 = t335 * t264 + t336 * t267;
t331 = -t336 * t264 + t335 * t267;
t330 = t334 * t264 + t337 * t267;
t329 = -t337 * t264 + t334 * t267;
t328 = t266 * t340 - t338;
t327 = t263 * t341 + t339;
t326 = t336 * t263 - t337 * t266;
t239 = -qJD(2) * t267 + V_base(5);
t240 = qJD(2) * t264 + V_base(4);
t250 = V_base(6) + qJD(1);
t325 = (t263 * t328 + t266 * t327) * t250 + (-t263 * t331 + t266 * t329) * t240 + (-t263 * t332 + t266 * t330) * t239;
t324 = (-t337 * t263 - t336 * t266) * t250 + (t264 * t333 + t326 * t267) * t240 + (t326 * t264 - t333 * t267) * t239;
t262 = sin(qJ(4));
t320 = pkin(4) * t262;
t319 = t263 * pkin(8);
t265 = cos(qJ(4));
t318 = t265 * pkin(4);
t316 = Icges(2,4) * t264;
t311 = t263 * t264;
t310 = t263 * t267;
t309 = t264 * t265;
t308 = t264 * t266;
t307 = t265 * t267;
t306 = t266 * t267;
t261 = qJ(4) + qJ(5);
t291 = pkin(2) * t266 + qJ(3) * t263;
t206 = t291 * t264;
t237 = pkin(1) * t264 - pkin(7) * t267;
t305 = -t206 - t237;
t255 = cos(t261);
t304 = pkin(5) * t255;
t302 = qJD(3) * t263;
t301 = qJD(4) * t266;
t300 = qJD(5) * t266;
t299 = qJD(6) * t266;
t298 = V_base(5) * pkin(6) + V_base(1);
t254 = sin(t261);
t295 = pkin(5) * t254;
t201 = t267 * t301 + t240;
t231 = qJD(4) * t263 + t250;
t232 = pkin(2) * t263 - qJ(3) * t266;
t294 = t239 * t232 + t267 * t302 + t298;
t293 = rSges(3,1) * t266 - rSges(3,2) * t263;
t292 = -rSges(4,2) * t266 + rSges(4,3) * t263;
t165 = t267 * t300 + t201;
t209 = qJD(5) * t263 + t231;
t200 = t264 * t301 + t239;
t238 = pkin(1) * t267 + pkin(7) * t264;
t284 = -V_base(4) * pkin(6) + t250 * t238 + V_base(2);
t283 = V_base(4) * t237 - t238 * V_base(5) + V_base(3);
t164 = t264 * t300 + t200;
t282 = pkin(9) * t266 + t263 * t320;
t207 = t291 * t267;
t279 = t250 * t207 + t264 * t302 + t284;
t278 = pkin(10) * t266 + t263 * t295;
t277 = -qJD(3) * t266 + t240 * t206 + t283;
t214 = -t267 * pkin(3) + pkin(8) * t308;
t276 = t239 * t319 + (-t214 + t305) * t250 + t294;
t147 = t264 * t282 - t267 * t318;
t198 = pkin(9) * t263 - t266 * t320;
t275 = -t147 * t231 + t200 * t198 + t276;
t213 = t264 * pkin(3) + pkin(8) * t306;
t274 = t250 * t213 + (-t232 - t319) * t240 + t279;
t273 = t240 * t214 + (-t207 - t213) * t239 + t277;
t146 = t264 * t318 + t267 * t282;
t272 = t231 * t146 - t198 * t201 + t274;
t271 = -t146 * t200 + t201 * t147 + t273;
t257 = qJ(6) + t261;
t256 = Icges(2,4) * t267;
t248 = cos(t257);
t247 = sin(t257);
t236 = rSges(2,1) * t267 - rSges(2,2) * t264;
t235 = rSges(2,1) * t264 + rSges(2,2) * t267;
t234 = rSges(3,1) * t263 + rSges(3,2) * t266;
t233 = -rSges(4,2) * t263 - rSges(4,3) * t266;
t230 = Icges(2,1) * t267 - t316;
t229 = Icges(2,1) * t264 + t256;
t227 = -Icges(2,2) * t264 + t256;
t226 = Icges(2,2) * t267 + t316;
t217 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t216 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t215 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t205 = t262 * t311 - t307;
t204 = t262 * t267 + t263 * t309;
t203 = t262 * t310 + t309;
t202 = -t262 * t264 + t263 * t307;
t196 = qJD(6) * t263 + t209;
t195 = t254 * t311 - t255 * t267;
t194 = t254 * t267 + t255 * t311;
t193 = t254 * t310 + t255 * t264;
t192 = -t254 * t264 + t255 * t310;
t191 = -rSges(4,1) * t267 + t264 * t292;
t190 = rSges(4,1) * t264 + t267 * t292;
t189 = rSges(3,3) * t264 + t267 * t293;
t188 = rSges(5,3) * t263 + (-rSges(5,1) * t262 - rSges(5,2) * t265) * t266;
t187 = -rSges(3,3) * t267 + t264 * t293;
t178 = Icges(5,5) * t263 + (-Icges(5,1) * t262 - Icges(5,4) * t265) * t266;
t175 = Icges(5,6) * t263 + (-Icges(5,4) * t262 - Icges(5,2) * t265) * t266;
t172 = Icges(5,3) * t263 + (-Icges(5,5) * t262 - Icges(5,6) * t265) * t266;
t169 = t247 * t311 - t248 * t267;
t168 = t247 * t267 + t248 * t311;
t167 = t247 * t310 + t248 * t264;
t166 = -t247 * t264 + t248 * t310;
t163 = rSges(6,3) * t263 + (-rSges(6,1) * t254 - rSges(6,2) * t255) * t266;
t161 = Icges(6,5) * t263 + (-Icges(6,1) * t254 - Icges(6,4) * t255) * t266;
t160 = Icges(6,6) * t263 + (-Icges(6,4) * t254 - Icges(6,2) * t255) * t266;
t159 = Icges(6,3) * t263 + (-Icges(6,5) * t254 - Icges(6,6) * t255) * t266;
t158 = rSges(7,3) * t263 + (-rSges(7,1) * t247 - rSges(7,2) * t248) * t266;
t157 = Icges(7,5) * t263 + (-Icges(7,1) * t247 - Icges(7,4) * t248) * t266;
t156 = Icges(7,6) * t263 + (-Icges(7,4) * t247 - Icges(7,2) * t248) * t266;
t155 = Icges(7,3) * t263 + (-Icges(7,5) * t247 - Icges(7,6) * t248) * t266;
t154 = V_base(5) * rSges(2,3) - t235 * t250 + t298;
t153 = t236 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t152 = t235 * V_base(4) - t236 * V_base(5) + V_base(3);
t151 = t267 * t299 + t165;
t150 = t264 * t299 + t164;
t148 = pkin(10) * t263 - t266 * t295;
t145 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t308;
t144 = rSges(5,1) * t203 + rSges(5,2) * t202 + rSges(5,3) * t306;
t143 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t308;
t142 = Icges(5,1) * t203 + Icges(5,4) * t202 + Icges(5,5) * t306;
t141 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t308;
t140 = Icges(5,4) * t203 + Icges(5,2) * t202 + Icges(5,6) * t306;
t139 = Icges(5,5) * t205 + Icges(5,6) * t204 + Icges(5,3) * t308;
t138 = Icges(5,5) * t203 + Icges(5,6) * t202 + Icges(5,3) * t306;
t136 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t308;
t135 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t306;
t134 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t308;
t133 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t306;
t132 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t308;
t131 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t306;
t130 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t308;
t129 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t306;
t128 = rSges(7,1) * t169 + rSges(7,2) * t168 + rSges(7,3) * t308;
t127 = rSges(7,1) * t167 + rSges(7,2) * t166 + rSges(7,3) * t306;
t125 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t308;
t124 = Icges(7,1) * t167 + Icges(7,4) * t166 + Icges(7,5) * t306;
t123 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t308;
t122 = Icges(7,4) * t167 + Icges(7,2) * t166 + Icges(7,6) * t306;
t121 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t308;
t120 = Icges(7,5) * t167 + Icges(7,6) * t166 + Icges(7,3) * t306;
t119 = t234 * t239 + (-t187 - t237) * t250 + t298;
t118 = t189 * t250 - t234 * t240 + t284;
t117 = t264 * t278 - t267 * t304;
t116 = t264 * t304 + t267 * t278;
t115 = t187 * t240 - t189 * t239 + t283;
t114 = t233 * t239 + (-t191 + t305) * t250 + t294;
t113 = t190 * t250 + (-t232 - t233) * t240 + t279;
t112 = t191 * t240 + (-t190 - t207) * t239 + t277;
t111 = -t145 * t231 + t188 * t200 + t276;
t110 = t144 * t231 - t188 * t201 + t274;
t109 = -t144 * t200 + t145 * t201 + t273;
t108 = -t136 * t209 + t163 * t164 + t275;
t107 = t135 * t209 - t163 * t165 + t272;
t106 = -t135 * t164 + t136 * t165 + t271;
t105 = -t117 * t209 - t128 * t196 + t148 * t164 + t150 * t158 + t275;
t104 = t116 * t209 + t127 * t196 - t148 * t165 - t151 * t158 + t272;
t103 = -t116 * t164 + t117 * t165 - t127 * t150 + t128 * t151 + t271;
t1 = t231 * ((t138 * t201 + t139 * t200 + t172 * t231) * t263 + ((-t140 * t265 - t142 * t262) * t201 + (-t141 * t265 - t143 * t262) * t200 + (-t175 * t265 - t178 * t262) * t231) * t266) / 0.2e1 + t209 * ((t129 * t165 + t130 * t164 + t159 * t209) * t263 + ((-t131 * t255 - t133 * t254) * t165 + (-t132 * t255 - t134 * t254) * t164 + (-t160 * t255 - t161 * t254) * t209) * t266) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t201 * ((t138 * t306 + t202 * t140 + t203 * t142) * t201 + (t139 * t306 + t141 * t202 + t143 * t203) * t200 + (t172 * t306 + t175 * t202 + t178 * t203) * t231) / 0.2e1 + t165 * ((t129 * t306 + t192 * t131 + t193 * t133) * t165 + (t130 * t306 + t132 * t192 + t134 * t193) * t164 + (t159 * t306 + t160 * t192 + t161 * t193) * t209) / 0.2e1 + t151 * ((t120 * t306 + t166 * t122 + t167 * t124) * t151 + (t121 * t306 + t123 * t166 + t125 * t167) * t150 + (t155 * t306 + t156 * t166 + t157 * t167) * t196) / 0.2e1 + m(3) * (t115 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t200 * ((t138 * t308 + t140 * t204 + t142 * t205) * t201 + (t139 * t308 + t204 * t141 + t205 * t143) * t200 + (t172 * t308 + t175 * t204 + t178 * t205) * t231) / 0.2e1 + t164 * ((t129 * t308 + t131 * t194 + t133 * t195) * t165 + (t130 * t308 + t194 * t132 + t195 * t134) * t164 + (t159 * t308 + t160 * t194 + t161 * t195) * t209) / 0.2e1 + t150 * ((t120 * t308 + t122 * t168 + t124 * t169) * t151 + (t121 * t308 + t168 * t123 + t169 * t125) * t150 + (t155 * t308 + t156 * t168 + t157 * t169) * t196) / 0.2e1 + m(2) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(1) * (t215 ^ 2 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + t196 * ((t120 * t151 + t121 * t150 + t155 * t196) * t263 + ((-t122 * t248 - t124 * t247) * t151 + (-t123 * t248 - t125 * t247) * t150 + (-t156 * t248 - t157 * t247) * t196) * t266) / 0.2e1 + (t325 * t264 - t324 * t267) * t239 / 0.2e1 + (t324 * t264 + t325 * t267) * t240 / 0.2e1 + ((-t226 * t264 + t229 * t267 + Icges(1,4)) * V_base(5) + (-t264 * t227 + t267 * t230 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t226 + t264 * t229 + Icges(1,2)) * V_base(5) + (t227 * t267 + t230 * t264 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t329 * t263 + t331 * t266) * t240 + (t330 * t263 + t332 * t266) * t239 + (t327 * t263 - t328 * t266 + Icges(2,3)) * t250) * t250 / 0.2e1 + t250 * V_base(4) * (Icges(2,5) * t267 - Icges(2,6) * t264) + t250 * V_base(5) * (Icges(2,5) * t264 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
