% Calculate kinetic energy for
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:20
% EndTime: 2019-03-09 03:44:24
% DurationCPUTime: 3.95s
% Computational Cost: add. (1779->323), mult. (1743->461), div. (0->0), fcn. (1579->10), ass. (0->158)
t324 = Icges(4,4) + Icges(5,6);
t323 = Icges(4,1) + Icges(5,2);
t322 = -Icges(4,2) - Icges(5,3);
t247 = cos(qJ(3));
t321 = t324 * t247;
t244 = sin(qJ(3));
t320 = t324 * t244;
t319 = Icges(5,4) - Icges(4,5);
t318 = Icges(5,5) - Icges(4,6);
t317 = t244 * t322 + t321;
t316 = t247 * t323 - t320;
t315 = Icges(5,1) + Icges(4,3);
t241 = qJ(1) + pkin(10);
t232 = sin(t241);
t233 = cos(t241);
t314 = t317 * t232 + t318 * t233;
t313 = -t318 * t232 + t317 * t233;
t312 = t316 * t232 + t319 * t233;
t311 = -t319 * t232 + t316 * t233;
t310 = t247 * t322 - t320;
t309 = t244 * t323 + t321;
t308 = t318 * t244 - t319 * t247;
t203 = -qJD(3) * t233 + V_base(5);
t204 = qJD(3) * t232 + V_base(4);
t234 = V_base(6) + qJD(1);
t305 = (t310 * t244 + t309 * t247) * t234 + (-t313 * t244 + t311 * t247) * t204 + (-t314 * t244 + t312 * t247) * t203;
t304 = (-t319 * t244 - t318 * t247) * t234 + (t315 * t232 + t308 * t233) * t204 + (t308 * t232 - t315 * t233) * t203;
t245 = sin(qJ(1));
t300 = pkin(1) * t245;
t248 = cos(qJ(1));
t299 = pkin(1) * t248;
t298 = t244 * pkin(8);
t246 = cos(qJ(5));
t297 = pkin(5) * t246;
t296 = -pkin(6) - qJ(2);
t294 = Icges(2,4) * t245;
t293 = Icges(3,4) * t232;
t288 = t232 * t247;
t287 = t233 * t247;
t242 = qJ(5) + qJ(6);
t236 = sin(t242);
t286 = t236 * t244;
t237 = cos(t242);
t285 = t237 * t244;
t243 = sin(qJ(5));
t284 = t243 * t244;
t283 = t244 * t246;
t282 = qJD(4) * t244;
t281 = qJD(5) * t247;
t280 = qJD(6) * t247;
t279 = t234 * t299 + V_base(2);
t278 = V_base(5) * pkin(6) + V_base(1);
t198 = pkin(2) * t232 - pkin(7) * t233;
t275 = -t198 - t300;
t274 = V_base(5) * qJ(2) + t278;
t273 = V_base(4) * t300 + qJD(2) + V_base(3);
t172 = t233 * t281 + t204;
t221 = qJD(5) * t244 + t234;
t269 = pkin(3) * t247 + qJ(4) * t244;
t180 = t269 * t232;
t272 = -t180 + t275;
t271 = rSges(4,1) * t247 - rSges(4,2) * t244;
t270 = -rSges(5,2) * t247 + rSges(5,3) * t244;
t222 = pkin(3) * t244 - qJ(4) * t247;
t262 = t203 * t222 + t233 * t282 + t274;
t171 = t232 * t281 + t203;
t261 = pkin(5) * t284 + pkin(9) * t247;
t199 = pkin(2) * t233 + pkin(7) * t232;
t258 = t234 * t199 + t296 * V_base(4) + t279;
t181 = t269 * t233;
t257 = t234 * t181 + t232 * t282 + t258;
t256 = V_base(4) * t198 + (-t199 - t299) * V_base(5) + t273;
t187 = -t233 * pkin(4) + pkin(8) * t288;
t255 = t203 * t298 + (-t187 + t272) * t234 + t262;
t254 = -qJD(4) * t247 + t204 * t180 + t256;
t186 = t232 * pkin(4) + pkin(8) * t287;
t253 = t234 * t186 + (-t222 - t298) * t204 + t257;
t252 = t204 * t187 + (-t181 - t186) * t203 + t254;
t239 = Icges(2,4) * t248;
t230 = Icges(3,4) * t233;
t226 = rSges(2,1) * t248 - rSges(2,2) * t245;
t225 = rSges(2,1) * t245 + rSges(2,2) * t248;
t224 = rSges(4,1) * t244 + rSges(4,2) * t247;
t223 = -rSges(5,2) * t244 - rSges(5,3) * t247;
t216 = Icges(2,1) * t248 - t294;
t215 = Icges(2,1) * t245 + t239;
t213 = -Icges(2,2) * t245 + t239;
t212 = Icges(2,2) * t248 + t294;
t202 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t201 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t200 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t197 = rSges(3,1) * t233 - rSges(3,2) * t232;
t196 = rSges(3,1) * t232 + rSges(3,2) * t233;
t195 = qJD(6) * t244 + t221;
t194 = Icges(3,1) * t233 - t293;
t193 = Icges(3,1) * t232 + t230;
t192 = -Icges(3,2) * t232 + t230;
t191 = Icges(3,2) * t233 + t293;
t184 = -t247 * t243 * pkin(5) + pkin(9) * t244;
t182 = rSges(6,3) * t244 + (-rSges(6,1) * t243 - rSges(6,2) * t246) * t247;
t179 = Icges(6,5) * t244 + (-Icges(6,1) * t243 - Icges(6,4) * t246) * t247;
t178 = Icges(6,6) * t244 + (-Icges(6,4) * t243 - Icges(6,2) * t246) * t247;
t177 = Icges(6,3) * t244 + (-Icges(6,5) * t243 - Icges(6,6) * t246) * t247;
t176 = t232 * t284 - t233 * t246;
t175 = t232 * t283 + t233 * t243;
t174 = t232 * t246 + t233 * t284;
t173 = -t232 * t243 + t233 * t283;
t168 = rSges(7,3) * t244 + (-rSges(7,1) * t236 - rSges(7,2) * t237) * t247;
t167 = Icges(7,5) * t244 + (-Icges(7,1) * t236 - Icges(7,4) * t237) * t247;
t166 = Icges(7,6) * t244 + (-Icges(7,4) * t236 - Icges(7,2) * t237) * t247;
t165 = Icges(7,3) * t244 + (-Icges(7,5) * t236 - Icges(7,6) * t237) * t247;
t164 = t232 * t286 - t233 * t237;
t163 = t232 * t285 + t233 * t236;
t162 = t232 * t237 + t233 * t286;
t161 = -t232 * t236 + t233 * t285;
t159 = -rSges(5,1) * t233 + t232 * t270;
t158 = rSges(5,1) * t232 + t233 * t270;
t157 = rSges(4,3) * t232 + t233 * t271;
t156 = -rSges(4,3) * t233 + t232 * t271;
t155 = V_base(5) * rSges(2,3) - t225 * t234 + t278;
t154 = t226 * t234 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t140 = t225 * V_base(4) - t226 * V_base(5) + V_base(3);
t139 = t233 * t280 + t172;
t138 = t232 * t280 + t171;
t136 = V_base(5) * rSges(3,3) + (-t196 - t300) * t234 + t274;
t135 = t197 * t234 + (-rSges(3,3) + t296) * V_base(4) + t279;
t134 = t232 * t261 - t233 * t297;
t133 = t232 * t297 + t233 * t261;
t132 = t196 * V_base(4) + (-t197 - t299) * V_base(5) + t273;
t131 = rSges(6,1) * t176 + rSges(6,2) * t175 + rSges(6,3) * t288;
t130 = rSges(6,1) * t174 + rSges(6,2) * t173 + rSges(6,3) * t287;
t129 = Icges(6,1) * t176 + Icges(6,4) * t175 + Icges(6,5) * t288;
t128 = Icges(6,1) * t174 + Icges(6,4) * t173 + Icges(6,5) * t287;
t127 = Icges(6,4) * t176 + Icges(6,2) * t175 + Icges(6,6) * t288;
t126 = Icges(6,4) * t174 + Icges(6,2) * t173 + Icges(6,6) * t287;
t125 = Icges(6,5) * t176 + Icges(6,6) * t175 + Icges(6,3) * t288;
t124 = Icges(6,5) * t174 + Icges(6,6) * t173 + Icges(6,3) * t287;
t123 = rSges(7,1) * t164 + rSges(7,2) * t163 + rSges(7,3) * t288;
t122 = rSges(7,1) * t162 + rSges(7,2) * t161 + rSges(7,3) * t287;
t121 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t288;
t120 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t287;
t119 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t288;
t118 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t287;
t117 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t288;
t116 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t287;
t115 = t203 * t224 + (-t156 + t275) * t234 + t274;
t114 = t157 * t234 - t204 * t224 + t258;
t113 = t156 * t204 - t157 * t203 + t256;
t112 = t203 * t223 + (-t159 + t272) * t234 + t262;
t111 = t158 * t234 + (-t222 - t223) * t204 + t257;
t110 = t159 * t204 + (-t158 - t181) * t203 + t254;
t109 = -t131 * t221 + t171 * t182 + t255;
t108 = t130 * t221 - t172 * t182 + t253;
t107 = -t130 * t171 + t131 * t172 + t252;
t106 = -t123 * t195 - t134 * t221 + t138 * t168 + t171 * t184 + t255;
t105 = t122 * t195 + t133 * t221 - t139 * t168 - t172 * t184 + t253;
t104 = -t122 * t138 + t123 * t139 - t133 * t171 + t134 * t172 + t252;
t1 = t172 * ((t124 * t287 + t173 * t126 + t174 * t128) * t172 + (t125 * t287 + t127 * t173 + t129 * t174) * t171 + (t173 * t178 + t174 * t179 + t177 * t287) * t221) / 0.2e1 + t139 * ((t116 * t287 + t161 * t118 + t162 * t120) * t139 + (t117 * t287 + t119 * t161 + t121 * t162) * t138 + (t161 * t166 + t162 * t167 + t165 * t287) * t195) / 0.2e1 + t171 * ((t124 * t288 + t126 * t175 + t128 * t176) * t172 + (t125 * t288 + t175 * t127 + t176 * t129) * t171 + (t175 * t178 + t176 * t179 + t177 * t288) * t221) / 0.2e1 + t138 * ((t116 * t288 + t118 * t163 + t120 * t164) * t139 + (t117 * t288 + t163 * t119 + t164 * t121) * t138 + (t163 * t166 + t164 * t167 + t165 * t288) * t195) / 0.2e1 + m(3) * (t132 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + t221 * ((t124 * t172 + t125 * t171 + t177 * t221) * t244 + ((-t126 * t246 - t128 * t243) * t172 + (-t127 * t246 - t129 * t243) * t171 + (-t178 * t246 - t179 * t243) * t221) * t247) / 0.2e1 + t195 * ((t116 * t139 + t117 * t138 + t165 * t195) * t244 + ((-t118 * t237 - t120 * t236) * t139 + (-t119 * t237 - t121 * t236) * t138 + (-t166 * t237 - t167 * t236) * t195) * t247) / 0.2e1 + m(1) * (t200 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(2) * (t140 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + (t305 * t232 - t304 * t233) * t203 / 0.2e1 + (t304 * t232 + t305 * t233) * t204 / 0.2e1 + ((-t191 * t232 + t193 * t233 - t212 * t245 + t215 * t248 + Icges(1,4)) * V_base(5) + (-t232 * t192 + t233 * t194 - t245 * t213 + t248 * t216 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t233 * t191 + t232 * t193 + t248 * t212 + t245 * t215 + Icges(1,2)) * V_base(5) + (t192 * t233 + t194 * t232 + t213 * t248 + t216 * t245 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t311 * t244 + t313 * t247) * t204 + (t312 * t244 + t314 * t247) * t203 + (t309 * t244 - t310 * t247 + Icges(2,3) + Icges(3,3)) * t234) * t234 / 0.2e1 + t234 * V_base(5) * (Icges(2,5) * t245 + Icges(3,5) * t232 + Icges(2,6) * t248 + Icges(3,6) * t233) + t234 * V_base(4) * (Icges(2,5) * t248 + Icges(3,5) * t233 - Icges(2,6) * t245 - Icges(3,6) * t232) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
