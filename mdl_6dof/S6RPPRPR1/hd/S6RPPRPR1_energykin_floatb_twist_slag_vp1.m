% Calculate kinetic energy for
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:38:57
% EndTime: 2019-03-09 01:39:00
% DurationCPUTime: 3.08s
% Computational Cost: add. (2148->351), mult. (1676->490), div. (0->0), fcn. (1512->12), ass. (0->172)
t235 = qJ(1) + pkin(9);
t225 = sin(t235);
t228 = cos(t235);
t242 = sin(qJ(1));
t243 = cos(qJ(1));
t305 = Icges(2,5) * t242 + Icges(3,5) * t225 + Icges(2,6) * t243 + Icges(3,6) * t228;
t304 = Icges(2,5) * t243 + Icges(3,5) * t228 - Icges(2,6) * t242 - Icges(3,6) * t225;
t300 = pkin(1) * t242;
t299 = pkin(1) * t243;
t237 = sin(pkin(10));
t298 = pkin(3) * t237;
t239 = cos(pkin(10));
t297 = pkin(3) * t239;
t238 = cos(pkin(11));
t296 = pkin(5) * t238;
t295 = -pkin(6) - qJ(2);
t294 = Icges(2,4) * t242;
t293 = Icges(3,4) * t225;
t292 = Icges(4,4) * t237;
t291 = Icges(4,4) * t239;
t234 = pkin(10) + qJ(4);
t224 = sin(t234);
t290 = Icges(5,4) * t224;
t227 = cos(t234);
t289 = Icges(5,4) * t227;
t288 = t224 * t225;
t287 = t224 * t228;
t286 = t225 * t227;
t236 = sin(pkin(11));
t285 = t225 * t236;
t284 = t225 * t238;
t283 = t227 * t228;
t282 = t228 * t236;
t281 = t228 * t238;
t278 = qJD(5) * t224;
t277 = qJD(6) * t224;
t229 = V_base(6) + qJD(1);
t276 = t229 * t299 + V_base(2);
t275 = V_base(5) * pkin(6) + V_base(1);
t205 = qJD(4) * t225 + V_base(4);
t191 = pkin(2) * t225 - qJ(3) * t228;
t272 = -t191 - t300;
t193 = pkin(2) * t228 + qJ(3) * t225;
t271 = -t193 - t299;
t270 = V_base(5) * qJ(2) + t275;
t269 = V_base(4) * t300 + qJD(2) + V_base(3);
t135 = -pkin(7) * t228 + t225 * t297;
t268 = -t135 + t272;
t267 = qJD(3) * t225 + t270;
t266 = V_base(4) * t191 + t269;
t204 = -qJD(4) * t228 + V_base(5);
t265 = rSges(4,1) * t239 - rSges(4,2) * t237;
t264 = rSges(5,1) * t227 - rSges(5,2) * t224;
t263 = pkin(4) * t227 + qJ(5) * t224;
t262 = Icges(4,1) * t239 - t292;
t261 = Icges(5,1) * t227 - t290;
t260 = -Icges(4,2) * t237 + t291;
t259 = -Icges(5,2) * t224 + t289;
t258 = Icges(4,5) * t239 - Icges(4,6) * t237;
t257 = Icges(5,5) * t227 - Icges(5,6) * t224;
t170 = t263 * t225;
t256 = -t170 + t268;
t255 = V_base(5) * t298 + t267;
t254 = -qJD(3) * t228 + t229 * t193 + t276;
t189 = pkin(4) * t224 - qJ(5) * t227;
t253 = t204 * t189 + t228 * t278 + t255;
t252 = (-Icges(5,3) * t228 + t225 * t257) * t204 + (Icges(5,3) * t225 + t228 * t257) * t205 + (Icges(5,5) * t224 + Icges(5,6) * t227) * t229;
t251 = pkin(8) * t224 + t227 * t296;
t250 = (-Icges(4,3) * t228 + t225 * t258) * V_base(5) + (Icges(4,3) * t225 + t228 * t258) * V_base(4) + (Icges(4,5) * t237 + Icges(4,6) * t239) * t229;
t136 = pkin(7) * t225 + t228 * t297;
t249 = V_base(4) * t135 + (-t136 + t271) * V_base(5) + t266;
t248 = t229 * t136 + (t295 - t298) * V_base(4) + t254;
t171 = t263 * t228;
t247 = t229 * t171 + t225 * t278 + t248;
t246 = -qJD(5) * t227 + t205 * t170 + t249;
t141 = -Icges(5,6) * t228 + t225 * t259;
t142 = Icges(5,6) * t225 + t228 * t259;
t144 = -Icges(5,5) * t228 + t225 * t261;
t145 = Icges(5,5) * t225 + t228 * t261;
t183 = Icges(5,2) * t227 + t290;
t186 = Icges(5,1) * t224 + t289;
t245 = (-t142 * t224 + t145 * t227) * t205 + (-t141 * t224 + t144 * t227) * t204 + (-t183 * t224 + t186 * t227) * t229;
t155 = -Icges(4,6) * t228 + t225 * t260;
t156 = Icges(4,6) * t225 + t228 * t260;
t158 = -Icges(4,5) * t228 + t225 * t262;
t159 = Icges(4,5) * t225 + t228 * t262;
t202 = Icges(4,2) * t239 + t292;
t203 = Icges(4,1) * t237 + t291;
t244 = (-t156 * t237 + t159 * t239) * V_base(4) + (-t155 * t237 + t158 * t239) * V_base(5) + (-t202 * t237 + t203 * t239) * t229;
t233 = pkin(11) + qJ(6);
t231 = Icges(2,4) * t243;
t226 = cos(t233);
t223 = sin(t233);
t220 = Icges(3,4) * t228;
t214 = rSges(2,1) * t243 - t242 * rSges(2,2);
t213 = t242 * rSges(2,1) + rSges(2,2) * t243;
t212 = Icges(2,1) * t243 - t294;
t211 = Icges(2,1) * t242 + t231;
t210 = -Icges(2,2) * t242 + t231;
t209 = Icges(2,2) * t243 + t294;
t206 = rSges(4,1) * t237 + rSges(4,2) * t239;
t200 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t199 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t198 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t197 = -qJD(6) * t227 + t229;
t194 = rSges(3,1) * t228 - rSges(3,2) * t225;
t192 = rSges(3,1) * t225 + rSges(3,2) * t228;
t190 = rSges(5,1) * t224 + rSges(5,2) * t227;
t188 = Icges(3,1) * t228 - t293;
t187 = Icges(3,1) * t225 + t220;
t185 = -Icges(3,2) * t225 + t220;
t184 = Icges(3,2) * t228 + t293;
t177 = t228 * t277 + t205;
t176 = t225 * t277 + t204;
t175 = t227 * t281 + t285;
t174 = -t227 * t282 + t284;
t173 = t227 * t284 - t282;
t172 = -t227 * t285 - t281;
t168 = V_base(5) * rSges(2,3) - t213 * t229 + t275;
t167 = t214 * t229 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t166 = t223 * t225 + t226 * t283;
t165 = -t223 * t283 + t225 * t226;
t164 = -t223 * t228 + t226 * t286;
t163 = -t223 * t286 - t226 * t228;
t162 = rSges(4,3) * t225 + t228 * t265;
t161 = -rSges(4,3) * t228 + t225 * t265;
t160 = -rSges(6,3) * t227 + (rSges(6,1) * t238 - rSges(6,2) * t236) * t224;
t157 = -Icges(6,5) * t227 + (Icges(6,1) * t238 - Icges(6,4) * t236) * t224;
t154 = -Icges(6,6) * t227 + (Icges(6,4) * t238 - Icges(6,2) * t236) * t224;
t151 = -Icges(6,3) * t227 + (Icges(6,5) * t238 - Icges(6,6) * t236) * t224;
t150 = t213 * V_base(4) - t214 * V_base(5) + V_base(3);
t148 = rSges(5,3) * t225 + t228 * t264;
t147 = -rSges(5,3) * t228 + t225 * t264;
t146 = -rSges(7,3) * t227 + (rSges(7,1) * t226 - rSges(7,2) * t223) * t224;
t143 = -Icges(7,5) * t227 + (Icges(7,1) * t226 - Icges(7,4) * t223) * t224;
t140 = -Icges(7,6) * t227 + (Icges(7,4) * t226 - Icges(7,2) * t223) * t224;
t137 = -Icges(7,3) * t227 + (Icges(7,5) * t226 - Icges(7,6) * t223) * t224;
t134 = -pkin(8) * t227 + t224 * t296;
t130 = V_base(5) * rSges(3,3) + (-t192 - t300) * t229 + t270;
t129 = t194 * t229 + (-rSges(3,3) + t295) * V_base(4) + t276;
t128 = V_base(4) * t192 + (-t194 - t299) * V_base(5) + t269;
t127 = rSges(6,1) * t175 + rSges(6,2) * t174 + rSges(6,3) * t287;
t126 = rSges(6,1) * t173 + rSges(6,2) * t172 + rSges(6,3) * t288;
t125 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t287;
t124 = Icges(6,1) * t173 + Icges(6,4) * t172 + Icges(6,5) * t288;
t123 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t287;
t122 = Icges(6,4) * t173 + Icges(6,2) * t172 + Icges(6,6) * t288;
t121 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t287;
t120 = Icges(6,5) * t173 + Icges(6,6) * t172 + Icges(6,3) * t288;
t119 = pkin(5) * t285 + t228 * t251;
t118 = -pkin(5) * t282 + t225 * t251;
t117 = rSges(7,1) * t166 + rSges(7,2) * t165 + rSges(7,3) * t287;
t116 = rSges(7,1) * t164 + rSges(7,2) * t163 + rSges(7,3) * t288;
t115 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t287;
t114 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t288;
t113 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t287;
t112 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t288;
t111 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t287;
t110 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t288;
t109 = t206 * V_base(5) + (-t161 + t272) * t229 + t267;
t108 = t162 * t229 + (-t206 + t295) * V_base(4) + t254;
t107 = V_base(4) * t161 + (-t162 + t271) * V_base(5) + t266;
t106 = t190 * t204 + (-t147 + t268) * t229 + t255;
t105 = t148 * t229 - t190 * t205 + t248;
t104 = t205 * t147 - t204 * t148 + t249;
t103 = t160 * t204 + (-t126 + t256) * t229 + t253;
t102 = t127 * t229 + (-t160 - t189) * t205 + t247;
t101 = t205 * t126 + (-t127 - t171) * t204 + t246;
t100 = -t116 * t197 + t134 * t204 + t146 * t176 + (-t118 + t256) * t229 + t253;
t99 = t117 * t197 + t119 * t229 - t146 * t177 + (-t134 - t189) * t205 + t247;
t98 = t177 * t116 - t176 * t117 + t205 * t118 + (-t119 - t171) * t204 + t246;
t1 = t197 * ((-t110 * t176 - t111 * t177 - t137 * t197) * t227 + ((-t113 * t223 + t115 * t226) * t177 + (-t112 * t223 + t114 * t226) * t176 + (-t140 * t223 + t143 * t226) * t197) * t224) / 0.2e1 + m(1) * (t198 ^ 2 + t199 ^ 2 + t200 ^ 2) / 0.2e1 + m(2) * (t150 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(3) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + t177 * ((t111 * t287 + t113 * t165 + t115 * t166) * t177 + (t110 * t287 + t112 * t165 + t114 * t166) * t176 + (t137 * t287 + t140 * t165 + t143 * t166) * t197) / 0.2e1 + t176 * ((t111 * t288 + t113 * t163 + t115 * t164) * t177 + (t110 * t288 + t112 * t163 + t114 * t164) * t176 + (t137 * t288 + t140 * t163 + t143 * t164) * t197) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + ((t121 * t288 + t123 * t172 + t125 * t173) * t205 + (t120 * t288 + t122 * t172 + t124 * t173) * t204 + (t151 * t288 + t154 * t172 + t157 * t173) * t229 + t225 * t245 - t228 * t252) * t204 / 0.2e1 + ((t121 * t287 + t123 * t174 + t125 * t175) * t205 + (t120 * t287 + t122 * t174 + t124 * t175) * t204 + (t151 * t287 + t154 * t174 + t157 * t175) * t229 + t225 * t252 + t228 * t245) * t205 / 0.2e1 + (t225 * t250 + t228 * t244 + t304 * t229 + (-t184 * t225 + t187 * t228 - t242 * t209 + t211 * t243 + Icges(1,4)) * V_base(5) + (-t185 * t225 + t188 * t228 - t242 * t210 + t212 * t243 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t225 * t244 - t228 * t250 + t305 * t229 + (t184 * t228 + t187 * t225 + t209 * t243 + t242 * t211 + Icges(1,2)) * V_base(5) + (t185 * t228 + t188 * t225 + t210 * t243 + t242 * t212 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t120 * t204 - t121 * t205) * t227 + ((-t123 * t236 + t125 * t238) * t205 + (-t122 * t236 + t124 * t238) * t204) * t224 + (t142 * t227 + t145 * t224) * t205 + (t141 * t227 + t144 * t224) * t204 + (t155 * t239 + t158 * t237 + t305) * V_base(5) + (t156 * t239 + t159 * t237 + t304) * V_base(4) + (t202 * t239 + t203 * t237 + Icges(2,3) + Icges(3,3) + (-t151 + t183) * t227 + (-t154 * t236 + t157 * t238 + t186) * t224) * t229) * t229 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
