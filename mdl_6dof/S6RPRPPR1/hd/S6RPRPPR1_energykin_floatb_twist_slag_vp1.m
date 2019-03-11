% Calculate kinetic energy for
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:37:55
% EndTime: 2019-03-09 02:37:57
% DurationCPUTime: 2.77s
% Computational Cost: add. (2192->345), mult. (1698->470), div. (0->0), fcn. (1534->12), ass. (0->168)
t307 = Icges(4,3) + Icges(5,3);
t232 = qJ(3) + pkin(10);
t222 = sin(t232);
t225 = cos(t232);
t238 = sin(qJ(3));
t240 = cos(qJ(3));
t306 = Icges(4,5) * t240 + Icges(5,5) * t225 - Icges(4,6) * t238 - Icges(5,6) * t222;
t233 = qJ(1) + pkin(9);
t223 = sin(t233);
t226 = cos(t233);
t286 = Icges(5,4) * t225;
t258 = -Icges(5,2) * t222 + t286;
t139 = -Icges(5,6) * t226 + t223 * t258;
t140 = Icges(5,6) * t223 + t226 * t258;
t287 = Icges(5,4) * t222;
t260 = Icges(5,1) * t225 - t287;
t142 = -Icges(5,5) * t226 + t223 * t260;
t143 = Icges(5,5) * t223 + t226 * t260;
t288 = Icges(4,4) * t240;
t259 = -Icges(4,2) * t238 + t288;
t159 = -Icges(4,6) * t226 + t223 * t259;
t160 = Icges(4,6) * t223 + t226 * t259;
t289 = Icges(4,4) * t238;
t261 = Icges(4,1) * t240 - t289;
t161 = -Icges(4,5) * t226 + t223 * t261;
t162 = Icges(4,5) * t223 + t226 * t261;
t182 = Icges(5,2) * t225 + t287;
t185 = Icges(5,1) * t222 + t286;
t200 = -qJD(3) * t226 + V_base(5);
t201 = qJD(3) * t223 + V_base(4);
t205 = Icges(4,2) * t240 + t289;
t208 = Icges(4,1) * t238 + t288;
t227 = V_base(6) + qJD(1);
t303 = (-t182 * t222 + t185 * t225 - t205 * t238 + t208 * t240) * t227 + (-t140 * t222 + t143 * t225 - t160 * t238 + t162 * t240) * t201 + (-t139 * t222 + t142 * t225 - t159 * t238 + t161 * t240) * t200;
t302 = (Icges(4,5) * t238 + Icges(5,5) * t222 + Icges(4,6) * t240 + Icges(5,6) * t225) * t227 + (t307 * t223 + t306 * t226) * t201 + (t306 * t223 - t307 * t226) * t200;
t239 = sin(qJ(1));
t298 = pkin(1) * t239;
t241 = cos(qJ(1));
t297 = pkin(1) * t241;
t296 = pkin(3) * t238;
t295 = pkin(3) * t240;
t235 = cos(pkin(11));
t294 = pkin(5) * t235;
t293 = -pkin(6) - qJ(2);
t291 = Icges(2,4) * t239;
t290 = Icges(3,4) * t223;
t285 = t222 * t223;
t284 = t222 * t226;
t283 = t223 * t225;
t234 = sin(pkin(11));
t282 = t223 * t234;
t281 = t223 * t235;
t280 = t225 * t226;
t279 = t226 * t234;
t278 = t226 * t235;
t134 = qJ(4) * t223 + t226 * t295;
t262 = pkin(4) * t225 + qJ(5) * t222;
t169 = t262 * t226;
t276 = -t134 - t169;
t275 = qJD(5) * t222;
t274 = qJD(6) * t222;
t273 = t227 * t297 + V_base(2);
t272 = V_base(5) * pkin(6) + V_base(1);
t192 = pkin(2) * t223 - pkin(7) * t226;
t269 = -t192 - t298;
t188 = pkin(4) * t222 - qJ(5) * t225;
t268 = -t188 - t296;
t267 = V_base(5) * qJ(2) + t272;
t266 = V_base(4) * t298 + qJD(2) + V_base(3);
t133 = -qJ(4) * t226 + t223 * t295;
t265 = -t133 + t269;
t264 = rSges(4,1) * t240 - rSges(4,2) * t238;
t263 = rSges(5,1) * t225 - rSges(5,2) * t222;
t168 = t262 * t223;
t255 = -t168 + t265;
t254 = qJD(4) * t223 + t200 * t296 + t267;
t253 = t200 * t188 + t226 * t275 + t254;
t250 = pkin(8) * t222 + t225 * t294;
t193 = pkin(2) * t226 + pkin(7) * t223;
t249 = t227 * t193 + t293 * V_base(4) + t273;
t248 = V_base(4) * t192 + (-t193 - t297) * V_base(5) + t266;
t247 = t201 * t133 + t248;
t246 = -qJD(4) * t226 + t227 * t134 + t249;
t245 = t227 * t169 + t223 * t275 + t246;
t244 = -qJD(5) * t225 + t201 * t168 + t247;
t231 = pkin(11) + qJ(6);
t229 = Icges(2,4) * t241;
t224 = cos(t231);
t221 = sin(t231);
t218 = Icges(3,4) * t226;
t213 = rSges(2,1) * t241 - t239 * rSges(2,2);
t212 = t239 * rSges(2,1) + rSges(2,2) * t241;
t211 = rSges(4,1) * t238 + rSges(4,2) * t240;
t210 = Icges(2,1) * t241 - t291;
t209 = Icges(2,1) * t239 + t229;
t207 = -Icges(2,2) * t239 + t229;
t206 = Icges(2,2) * t241 + t291;
t199 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t198 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t197 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t196 = -qJD(6) * t225 + t227;
t191 = rSges(3,1) * t226 - rSges(3,2) * t223;
t190 = rSges(3,1) * t223 + rSges(3,2) * t226;
t189 = rSges(5,1) * t222 + rSges(5,2) * t225;
t187 = Icges(3,1) * t226 - t290;
t186 = Icges(3,1) * t223 + t218;
t184 = -Icges(3,2) * t223 + t218;
t183 = Icges(3,2) * t226 + t290;
t175 = t226 * t274 + t201;
t174 = t223 * t274 + t200;
t173 = t225 * t278 + t282;
t172 = -t225 * t279 + t281;
t171 = t225 * t281 - t279;
t170 = -t225 * t282 - t278;
t166 = rSges(4,3) * t223 + t226 * t264;
t165 = -rSges(4,3) * t226 + t223 * t264;
t164 = V_base(5) * rSges(2,3) - t212 * t227 + t272;
t163 = t213 * t227 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t156 = t221 * t223 + t224 * t280;
t155 = -t221 * t280 + t223 * t224;
t154 = -t221 * t226 + t224 * t283;
t153 = -t221 * t283 - t224 * t226;
t152 = -rSges(6,3) * t225 + (rSges(6,1) * t235 - rSges(6,2) * t234) * t222;
t151 = -Icges(6,5) * t225 + (Icges(6,1) * t235 - Icges(6,4) * t234) * t222;
t150 = -Icges(6,6) * t225 + (Icges(6,4) * t235 - Icges(6,2) * t234) * t222;
t149 = -Icges(6,3) * t225 + (Icges(6,5) * t235 - Icges(6,6) * t234) * t222;
t148 = t212 * V_base(4) - t213 * V_base(5) + V_base(3);
t146 = rSges(5,3) * t223 + t226 * t263;
t145 = -rSges(5,3) * t226 + t223 * t263;
t144 = -rSges(7,3) * t225 + (rSges(7,1) * t224 - rSges(7,2) * t221) * t222;
t141 = -Icges(7,5) * t225 + (Icges(7,1) * t224 - Icges(7,4) * t221) * t222;
t138 = -Icges(7,6) * t225 + (Icges(7,4) * t224 - Icges(7,2) * t221) * t222;
t135 = -Icges(7,3) * t225 + (Icges(7,5) * t224 - Icges(7,6) * t221) * t222;
t132 = -pkin(8) * t225 + t222 * t294;
t129 = V_base(5) * rSges(3,3) + (-t190 - t298) * t227 + t267;
t128 = t191 * t227 + (-rSges(3,3) + t293) * V_base(4) + t273;
t126 = V_base(4) * t190 + (-t191 - t297) * V_base(5) + t266;
t125 = rSges(6,1) * t173 + rSges(6,2) * t172 + rSges(6,3) * t284;
t124 = rSges(6,1) * t171 + rSges(6,2) * t170 + rSges(6,3) * t285;
t123 = Icges(6,1) * t173 + Icges(6,4) * t172 + Icges(6,5) * t284;
t122 = Icges(6,1) * t171 + Icges(6,4) * t170 + Icges(6,5) * t285;
t121 = Icges(6,4) * t173 + Icges(6,2) * t172 + Icges(6,6) * t284;
t120 = Icges(6,4) * t171 + Icges(6,2) * t170 + Icges(6,6) * t285;
t119 = Icges(6,5) * t173 + Icges(6,6) * t172 + Icges(6,3) * t284;
t118 = Icges(6,5) * t171 + Icges(6,6) * t170 + Icges(6,3) * t285;
t117 = pkin(5) * t282 + t226 * t250;
t116 = -pkin(5) * t279 + t223 * t250;
t115 = rSges(7,1) * t156 + rSges(7,2) * t155 + rSges(7,3) * t284;
t114 = rSges(7,1) * t154 + rSges(7,2) * t153 + rSges(7,3) * t285;
t113 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t284;
t112 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t285;
t111 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t284;
t110 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t285;
t109 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t284;
t108 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t285;
t107 = t200 * t211 + (-t165 + t269) * t227 + t267;
t106 = t166 * t227 - t201 * t211 + t249;
t105 = t201 * t165 - t200 * t166 + t248;
t104 = t189 * t200 + (-t145 + t265) * t227 + t254;
t103 = t146 * t227 + (-t189 - t296) * t201 + t246;
t102 = t201 * t145 + (-t134 - t146) * t200 + t247;
t101 = t152 * t200 + (-t124 + t255) * t227 + t253;
t100 = t125 * t227 + (-t152 + t268) * t201 + t245;
t99 = t201 * t124 + (-t125 + t276) * t200 + t244;
t98 = -t114 * t196 + t132 * t200 + t144 * t174 + (-t116 + t255) * t227 + t253;
t97 = t115 * t196 + t117 * t227 - t144 * t175 + (-t132 + t268) * t201 + t245;
t96 = t175 * t114 - t174 * t115 + t201 * t116 + (-t117 + t276) * t200 + t244;
t1 = t175 * ((t109 * t284 + t155 * t111 + t156 * t113) * t175 + (t108 * t284 + t110 * t155 + t112 * t156) * t174 + (t135 * t284 + t138 * t155 + t141 * t156) * t196) / 0.2e1 + t174 * ((t109 * t285 + t111 * t153 + t113 * t154) * t175 + (t108 * t285 + t153 * t110 + t154 * t112) * t174 + (t135 * t285 + t138 * t153 + t141 * t154) * t196) / 0.2e1 + t196 * ((-t108 * t174 - t109 * t175 - t135 * t196) * t225 + ((-t111 * t221 + t113 * t224) * t175 + (-t110 * t221 + t112 * t224) * t174 + (-t138 * t221 + t141 * t224) * t196) * t222) / 0.2e1 + m(1) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + m(2) * (t148 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(3) * (t126 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + ((t119 * t285 + t121 * t170 + t123 * t171) * t201 + (t118 * t285 + t170 * t120 + t171 * t122) * t200 + (t149 * t285 + t150 * t170 + t151 * t171) * t227 - t302 * t226 + t303 * t223) * t200 / 0.2e1 + ((t119 * t284 + t172 * t121 + t173 * t123) * t201 + (t118 * t284 + t120 * t172 + t122 * t173) * t200 + (t149 * t284 + t150 * t172 + t151 * t173) * t227 + t303 * t226 + t302 * t223) * t201 / 0.2e1 + ((-t183 * t223 + t186 * t226 - t239 * t206 + t209 * t241 + Icges(1,4)) * V_base(5) + (-t223 * t184 + t226 * t187 - t239 * t207 + t241 * t210 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t226 * t183 + t223 * t186 + t241 * t206 + t239 * t209 + Icges(1,2)) * V_base(5) + (t184 * t226 + t187 * t223 + t207 * t241 + t239 * t210 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t160 * t240 + t162 * t238 + (-t119 + t140) * t225 + (-t121 * t234 + t123 * t235 + t143) * t222) * t201 + (t159 * t240 + t161 * t238 + (-t118 + t139) * t225 + (-t120 * t234 + t122 * t235 + t142) * t222) * t200 + (t240 * t205 + t238 * t208 + Icges(2,3) + Icges(3,3) + (-t149 + t182) * t225 + (-t150 * t234 + t151 * t235 + t185) * t222) * t227) * t227 / 0.2e1 + t227 * V_base(5) * (Icges(2,5) * t239 + Icges(3,5) * t223 + Icges(2,6) * t241 + Icges(3,6) * t226) + t227 * V_base(4) * (Icges(2,5) * t241 + Icges(3,5) * t226 - Icges(2,6) * t239 - Icges(3,6) * t223) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
