% Calculate kinetic energy for
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:45
% EndTime: 2019-03-09 10:15:48
% DurationCPUTime: 3.84s
% Computational Cost: add. (2399->377), mult. (2318->538), div. (0->0), fcn. (2212->12), ass. (0->180)
t336 = Icges(3,3) + Icges(4,3);
t335 = -Icges(6,3) - Icges(5,3);
t256 = qJ(2) + pkin(10);
t244 = sin(t256);
t246 = cos(t256);
t260 = sin(qJ(2));
t263 = cos(qJ(2));
t334 = Icges(3,5) * t263 + Icges(4,5) * t246 - Icges(3,6) * t260 - Icges(4,6) * t244;
t255 = qJ(4) + pkin(11);
t245 = cos(t255);
t264 = cos(qJ(1));
t243 = sin(t255);
t261 = sin(qJ(1));
t305 = t261 * t243;
t184 = -t245 * t264 - t246 * t305;
t304 = t261 * t245;
t185 = -t243 * t264 + t246 * t304;
t262 = cos(qJ(4));
t301 = t262 * t264;
t259 = sin(qJ(4));
t303 = t261 * t259;
t201 = -t246 * t303 - t301;
t302 = t261 * t262;
t308 = t259 * t264;
t202 = t246 * t302 - t308;
t311 = t244 * t261;
t333 = Icges(5,5) * t202 + Icges(6,5) * t185 + Icges(5,6) * t201 + Icges(6,6) * t184 - t335 * t311;
t309 = t246 * t264;
t186 = -t243 * t309 + t304;
t187 = t245 * t309 + t305;
t203 = -t246 * t308 + t302;
t204 = t246 * t301 + t303;
t310 = t244 * t264;
t332 = Icges(5,5) * t204 + Icges(6,5) * t187 + Icges(5,6) * t203 + Icges(6,6) * t186 - t335 * t310;
t331 = t335 * t246 + (Icges(5,5) * t262 + Icges(6,5) * t245 - Icges(5,6) * t259 - Icges(6,6) * t243) * t244;
t312 = Icges(4,4) * t246;
t283 = -Icges(4,2) * t244 + t312;
t177 = -Icges(4,6) * t264 + t261 * t283;
t178 = Icges(4,6) * t261 + t264 * t283;
t313 = Icges(4,4) * t244;
t285 = Icges(4,1) * t246 - t313;
t179 = -Icges(4,5) * t264 + t261 * t285;
t180 = Icges(4,5) * t261 + t264 * t285;
t314 = Icges(3,4) * t263;
t284 = -Icges(3,2) * t260 + t314;
t190 = -Icges(3,6) * t264 + t261 * t284;
t191 = Icges(3,6) * t261 + t264 * t284;
t315 = Icges(3,4) * t260;
t286 = Icges(3,1) * t263 - t315;
t192 = -Icges(3,5) * t264 + t261 * t286;
t193 = Icges(3,5) * t261 + t264 * t286;
t207 = Icges(4,2) * t246 + t313;
t208 = Icges(4,1) * t244 + t312;
t222 = Icges(3,2) * t263 + t315;
t225 = Icges(3,1) * t260 + t314;
t237 = -qJD(2) * t264 + V_base(5);
t238 = qJD(2) * t261 + V_base(4);
t248 = V_base(6) + qJD(1);
t330 = (-t207 * t244 + t208 * t246 - t222 * t260 + t225 * t263) * t248 + (-t178 * t244 + t180 * t246 - t191 * t260 + t193 * t263) * t238 + (-t177 * t244 + t179 * t246 - t190 * t260 + t192 * t263) * t237;
t329 = (Icges(3,5) * t260 + Icges(4,5) * t244 + Icges(3,6) * t263 + Icges(4,6) * t246) * t248 + (t336 * t261 + t334 * t264) * t238 + (t334 * t261 - t336 * t264) * t237;
t322 = pkin(2) * t260;
t320 = pkin(2) * t263;
t319 = t262 * pkin(4);
t316 = Icges(2,4) * t261;
t247 = qJ(6) + t255;
t239 = sin(t247);
t307 = t261 * t239;
t240 = cos(t247);
t306 = t261 * t240;
t172 = -qJ(3) * t264 + t261 * t320;
t235 = t261 * pkin(1) - pkin(7) * t264;
t300 = -t172 - t235;
t299 = pkin(5) * t245;
t297 = qJD(4) * t244;
t296 = qJD(5) * t244;
t295 = qJD(6) * t244;
t294 = V_base(5) * pkin(6) + V_base(1);
t291 = pkin(5) * t243;
t200 = t264 * t297 + t238;
t290 = qJD(3) * t261 + t237 * t322 + t294;
t289 = pkin(3) * t246 + pkin(8) * t244;
t288 = rSges(3,1) * t263 - rSges(3,2) * t260;
t287 = rSges(4,1) * t246 - rSges(4,2) * t244;
t199 = t261 * t297 + t237;
t236 = pkin(1) * t264 + t261 * pkin(7);
t280 = -V_base(4) * pkin(6) + t248 * t236 + V_base(2);
t279 = V_base(4) * t235 - t236 * V_base(5) + V_base(3);
t278 = t238 * t172 + t279;
t277 = qJ(5) * t244 + t246 * t319;
t274 = pkin(9) * t244 + t246 * t299;
t173 = qJ(3) * t261 + t264 * t320;
t273 = -qJD(3) * t264 + t248 * t173 + t280;
t196 = t289 * t261;
t210 = pkin(3) * t244 - pkin(8) * t246;
t272 = t237 * t210 + (-t196 + t300) * t248 + t290;
t197 = t289 * t264;
t271 = t238 * t196 + (-t173 - t197) * t237 + t278;
t147 = -qJ(5) * t246 + t244 * t319;
t270 = t199 * t147 + t264 * t296 + t272;
t269 = t248 * t197 + (-t210 - t322) * t238 + t273;
t136 = -pkin(4) * t308 + t261 * t277;
t268 = -qJD(5) * t246 + t200 * t136 + t271;
t137 = pkin(4) * t303 + t264 * t277;
t214 = -qJD(4) * t246 + t248;
t267 = t214 * t137 + t261 * t296 + t269;
t251 = Icges(2,4) * t264;
t234 = rSges(2,1) * t264 - t261 * rSges(2,2);
t233 = t261 * rSges(2,1) + rSges(2,2) * t264;
t232 = rSges(3,1) * t260 + rSges(3,2) * t263;
t227 = Icges(2,1) * t264 - t316;
t226 = Icges(2,1) * t261 + t251;
t224 = -Icges(2,2) * t261 + t251;
t223 = Icges(2,2) * t264 + t316;
t217 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t216 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t215 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t209 = rSges(4,1) * t244 + rSges(4,2) * t246;
t198 = (-qJD(4) - qJD(6)) * t246 + t248;
t195 = t261 * rSges(3,3) + t264 * t288;
t194 = -rSges(3,3) * t264 + t261 * t288;
t182 = t261 * rSges(4,3) + t264 * t287;
t181 = -rSges(4,3) * t264 + t261 * t287;
t171 = t240 * t309 + t307;
t170 = -t239 * t309 + t306;
t169 = -t239 * t264 + t246 * t306;
t168 = -t240 * t264 - t246 * t307;
t167 = -rSges(5,3) * t246 + (rSges(5,1) * t262 - rSges(5,2) * t259) * t244;
t166 = V_base(5) * rSges(2,3) - t233 * t248 + t294;
t165 = t234 * t248 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t164 = -Icges(5,5) * t246 + (Icges(5,1) * t262 - Icges(5,4) * t259) * t244;
t163 = -Icges(5,6) * t246 + (Icges(5,4) * t262 - Icges(5,2) * t259) * t244;
t161 = t264 * t295 + t200;
t160 = t261 * t295 + t199;
t158 = t233 * V_base(4) - t234 * V_base(5) + V_base(3);
t156 = -rSges(6,3) * t246 + (rSges(6,1) * t245 - rSges(6,2) * t243) * t244;
t155 = -Icges(6,5) * t246 + (Icges(6,1) * t245 - Icges(6,4) * t243) * t244;
t154 = -Icges(6,6) * t246 + (Icges(6,4) * t245 - Icges(6,2) * t243) * t244;
t152 = -rSges(7,3) * t246 + (rSges(7,1) * t240 - rSges(7,2) * t239) * t244;
t150 = -Icges(7,5) * t246 + (Icges(7,1) * t240 - Icges(7,4) * t239) * t244;
t149 = -Icges(7,6) * t246 + (Icges(7,4) * t240 - Icges(7,2) * t239) * t244;
t148 = -Icges(7,3) * t246 + (Icges(7,5) * t240 - Icges(7,6) * t239) * t244;
t146 = -pkin(9) * t246 + t244 * t299;
t145 = t204 * rSges(5,1) + t203 * rSges(5,2) + rSges(5,3) * t310;
t144 = rSges(5,1) * t202 + rSges(5,2) * t201 + rSges(5,3) * t311;
t143 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t310;
t142 = Icges(5,1) * t202 + Icges(5,4) * t201 + Icges(5,5) * t311;
t141 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t310;
t140 = Icges(5,4) * t202 + Icges(5,2) * t201 + Icges(5,6) * t311;
t134 = t187 * rSges(6,1) + t186 * rSges(6,2) + rSges(6,3) * t310;
t133 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t311;
t132 = Icges(6,1) * t187 + Icges(6,4) * t186 + Icges(6,5) * t310;
t131 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t311;
t130 = Icges(6,4) * t187 + Icges(6,2) * t186 + Icges(6,6) * t310;
t129 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t311;
t126 = t171 * rSges(7,1) + t170 * rSges(7,2) + rSges(7,3) * t310;
t125 = rSges(7,1) * t169 + rSges(7,2) * t168 + rSges(7,3) * t311;
t124 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t310;
t123 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t311;
t122 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t310;
t121 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t311;
t120 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t310;
t119 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t311;
t117 = t232 * t237 + (-t194 - t235) * t248 + t294;
t116 = t195 * t248 - t232 * t238 + t280;
t114 = t261 * t291 + t264 * t274;
t113 = t261 * t274 - t264 * t291;
t112 = t194 * t238 - t195 * t237 + t279;
t111 = t209 * t237 + (-t181 + t300) * t248 + t290;
t110 = t248 * t182 + (-t209 - t322) * t238 + t273;
t109 = t181 * t238 + (-t173 - t182) * t237 + t278;
t108 = -t144 * t214 + t167 * t199 + t272;
t107 = t214 * t145 - t200 * t167 + t269;
t106 = t144 * t200 - t145 * t199 + t271;
t105 = t156 * t199 + (-t133 - t136) * t214 + t270;
t104 = t214 * t134 + (-t147 - t156) * t200 + t267;
t103 = t133 * t200 + (-t134 - t137) * t199 + t268;
t102 = -t125 * t198 + t146 * t199 + t152 * t160 + (-t113 - t136) * t214 + t270;
t101 = t214 * t114 + t198 * t126 - t161 * t152 + (-t146 - t147) * t200 + t267;
t100 = t113 * t200 + t125 * t161 - t126 * t160 + (-t114 - t137) * t199 + t268;
t1 = m(1) * (t215 ^ 2 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + m(3) * (t112 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(4) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(2) * (t158 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + t161 * ((t120 * t310 + t170 * t122 + t171 * t124) * t161 + (t119 * t310 + t170 * t121 + t171 * t123) * t160 + (t148 * t310 + t170 * t149 + t171 * t150) * t198) / 0.2e1 + t160 * ((t120 * t311 + t122 * t168 + t124 * t169) * t161 + (t119 * t311 + t168 * t121 + t169 * t123) * t160 + (t148 * t311 + t149 * t168 + t150 * t169) * t198) / 0.2e1 + t198 * ((-t119 * t160 - t120 * t161 - t148 * t198) * t246 + ((-t122 * t239 + t124 * t240) * t161 + (-t121 * t239 + t123 * t240) * t160 + (-t149 * t239 + t150 * t240) * t198) * t244) / 0.2e1 + ((t154 * t184 + t155 * t185 + t163 * t201 + t164 * t202 + t311 * t331) * t214 + (t130 * t184 + t132 * t185 + t141 * t201 + t143 * t202 + t311 * t332) * t200 + (t184 * t129 + t185 * t131 + t201 * t140 + t202 * t142 + t333 * t311) * t199) * t199 / 0.2e1 + ((t186 * t154 + t187 * t155 + t203 * t163 + t204 * t164 + t310 * t331) * t214 + (t186 * t130 + t187 * t132 + t203 * t141 + t204 * t143 + t332 * t310) * t200 + (t186 * t129 + t187 * t131 + t203 * t140 + t204 * t142 + t310 * t333) * t199) * t200 / 0.2e1 + ((-t199 * t333 - t332 * t200 - t331 * t214) * t246 + ((-t154 * t243 + t155 * t245 - t163 * t259 + t164 * t262) * t214 + (-t130 * t243 + t132 * t245 - t141 * t259 + t143 * t262) * t200 + (-t129 * t243 + t131 * t245 - t140 * t259 + t142 * t262) * t199) * t244) * t214 / 0.2e1 + (t330 * t261 - t329 * t264) * t237 / 0.2e1 + (t329 * t261 + t330 * t264) * t238 / 0.2e1 + ((-t261 * t223 + t226 * t264 + Icges(1,4)) * V_base(5) + (-t261 * t224 + t227 * t264 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t223 * t264 + t261 * t226 + Icges(1,2)) * V_base(5) + (t224 * t264 + t261 * t227 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t178 * t246 + t180 * t244 + t191 * t263 + t193 * t260) * t238 + (t177 * t246 + t179 * t244 + t190 * t263 + t192 * t260) * t237 + (t207 * t246 + t208 * t244 + t222 * t263 + t225 * t260 + Icges(2,3)) * t248) * t248 / 0.2e1 + t248 * V_base(4) * (Icges(2,5) * t264 - Icges(2,6) * t261) + t248 * V_base(5) * (Icges(2,5) * t261 + Icges(2,6) * t264) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
