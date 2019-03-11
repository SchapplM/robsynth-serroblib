% Calculate kinetic energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:18
% EndTime: 2019-03-09 05:34:21
% DurationCPUTime: 3.51s
% Computational Cost: add. (1206->316), mult. (2402->441), div. (0->0), fcn. (2460->8), ass. (0->149)
t309 = Icges(2,4) + Icges(3,6);
t308 = Icges(2,1) + Icges(3,2);
t307 = Icges(5,1) + Icges(6,1);
t306 = -Icges(3,4) + Icges(2,5);
t305 = Icges(5,4) - Icges(6,5);
t304 = Icges(6,4) + Icges(5,5);
t303 = Icges(3,5) - Icges(2,6);
t302 = Icges(2,2) + Icges(3,3);
t301 = Icges(5,2) + Icges(6,3);
t300 = Icges(6,6) - Icges(5,6);
t299 = Icges(5,3) + Icges(6,2);
t237 = cos(qJ(1));
t298 = t309 * t237;
t233 = sin(qJ(1));
t297 = t309 * t233;
t232 = sin(qJ(3));
t235 = cos(qJ(4));
t264 = t237 * t235;
t231 = sin(qJ(4));
t268 = t233 * t231;
t177 = t232 * t268 - t264;
t267 = t233 * t235;
t269 = t231 * t237;
t178 = t232 * t267 + t269;
t236 = cos(qJ(3));
t266 = t233 * t236;
t296 = t177 * t301 - t178 * t305 - t266 * t300;
t179 = t232 * t269 + t267;
t180 = -t232 * t264 + t268;
t265 = t236 * t237;
t295 = -t179 * t301 - t180 * t305 + t265 * t300;
t294 = t177 * t300 + t178 * t304 - t266 * t299;
t293 = -t179 * t300 + t180 * t304 + t265 * t299;
t292 = -t177 * t305 + t178 * t307 - t266 * t304;
t291 = t179 * t305 + t180 * t307 + t265 * t304;
t290 = (t231 * t301 - t235 * t305) * t236 + t300 * t232;
t289 = (t231 * t300 + t235 * t304) * t236 + t299 * t232;
t288 = (-t231 * t305 + t235 * t307) * t236 + t304 * t232;
t287 = -t237 * t302 - t297;
t286 = t233 * t302 - t298;
t285 = t233 * t308 + t298;
t284 = t237 * t308 - t297;
t273 = Icges(4,4) * t232;
t252 = Icges(4,2) * t236 + t273;
t160 = Icges(4,6) * t237 + t233 * t252;
t161 = Icges(4,6) * t233 - t237 * t252;
t272 = Icges(4,4) * t236;
t253 = Icges(4,1) * t232 + t272;
t164 = Icges(4,5) * t237 + t233 * t253;
t165 = Icges(4,5) * t233 - t237 * t253;
t198 = -Icges(4,2) * t232 + t272;
t203 = Icges(4,1) * t236 - t273;
t216 = qJD(3) * t233 + V_base(5);
t217 = qJD(3) * t237 + V_base(4);
t222 = V_base(6) + qJD(1);
t281 = (t160 * t236 + t164 * t232) * t217 + (t161 * t236 + t165 * t232) * t216 + (t198 * t236 + t203 * t232) * t222;
t276 = pkin(7) * t233;
t275 = pkin(7) * t237;
t263 = qJD(4) * t236;
t207 = t233 * pkin(1) - qJ(2) * t237;
t262 = V_base(4) * t207 + V_base(3);
t261 = V_base(5) * pkin(6) + V_base(1);
t258 = -t207 - t276;
t257 = qJD(2) * t233 + t261;
t175 = t237 * t263 + t216;
t206 = qJD(4) * t232 + t222;
t256 = V_base(5) * pkin(2) + t257;
t255 = pkin(3) * t232 - pkin(8) * t236;
t254 = rSges(4,1) * t232 + rSges(4,2) * t236;
t251 = Icges(4,5) * t232 + Icges(4,6) * t236;
t211 = pkin(1) * t237 + t233 * qJ(2);
t247 = -qJD(2) * t237 + t222 * t211 + V_base(2);
t246 = (Icges(4,3) * t237 + t233 * t251) * t217 + (Icges(4,3) * t233 - t237 * t251) * t216 + (Icges(4,5) * t236 - Icges(4,6) * t232) * t222;
t245 = V_base(4) * t276 + (-t211 - t275) * V_base(5) + t262;
t244 = t222 * t275 + (-pkin(2) - pkin(6)) * V_base(4) + t247;
t184 = t255 * t237;
t214 = pkin(3) * t236 + pkin(8) * t232;
t243 = t216 * t214 + (t184 + t258) * t222 + t256;
t183 = t255 * t233;
t242 = -t216 * t183 - t217 * t184 + t245;
t182 = (pkin(4) * t235 + qJ(5) * t231) * t236;
t241 = qJD(5) * t177 + t175 * t182 + t243;
t240 = t222 * t183 - t217 * t214 + t244;
t143 = pkin(4) * t180 - qJ(5) * t179;
t176 = -t233 * t263 + t217;
t239 = qJD(5) * t236 * t231 + t176 * t143 + t242;
t142 = pkin(4) * t178 + qJ(5) * t177;
t238 = -qJD(5) * t179 + t206 * t142 + t240;
t234 = cos(qJ(6));
t230 = sin(qJ(6));
t213 = rSges(2,1) * t237 - t233 * rSges(2,2);
t212 = -rSges(3,2) * t237 + t233 * rSges(3,3);
t210 = rSges(4,1) * t236 - rSges(4,2) * t232;
t209 = t233 * rSges(2,1) + rSges(2,2) * t237;
t208 = -t233 * rSges(3,2) - rSges(3,3) * t237;
t190 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t189 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t188 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t187 = pkin(5) * t235 * t236 - pkin(9) * t232;
t185 = -qJD(6) * t232 + t206;
t171 = (t230 * t231 + t234 * t235) * t236;
t170 = (-t230 * t235 + t231 * t234) * t236;
t169 = t233 * rSges(4,3) - t237 * t254;
t168 = rSges(5,3) * t232 + (rSges(5,1) * t235 - rSges(5,2) * t231) * t236;
t167 = rSges(6,2) * t232 + (rSges(6,1) * t235 + rSges(6,3) * t231) * t236;
t166 = rSges(4,3) * t237 + t233 * t254;
t152 = (-qJD(4) + qJD(6)) * t266 + t217;
t151 = -qJD(6) * t265 + t175;
t149 = t180 * pkin(5) - pkin(9) * t265;
t148 = pkin(5) * t178 + pkin(9) * t266;
t147 = V_base(5) * rSges(2,3) - t209 * t222 + t261;
t146 = t213 * t222 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = t209 * V_base(4) - t213 * V_base(5) + V_base(3);
t141 = -t179 * t230 + t180 * t234;
t140 = -t179 * t234 - t180 * t230;
t139 = t177 * t230 + t178 * t234;
t138 = t177 * t234 - t178 * t230;
t137 = t180 * rSges(5,1) + t179 * rSges(5,2) + rSges(5,3) * t265;
t136 = t180 * rSges(6,1) + rSges(6,2) * t265 - t179 * rSges(6,3);
t135 = rSges(5,1) * t178 - rSges(5,2) * t177 - rSges(5,3) * t266;
t134 = rSges(6,1) * t178 - rSges(6,2) * t266 + rSges(6,3) * t177;
t121 = rSges(7,1) * t171 + rSges(7,2) * t170 - rSges(7,3) * t232;
t119 = Icges(7,1) * t171 + Icges(7,4) * t170 - Icges(7,5) * t232;
t118 = Icges(7,4) * t171 + Icges(7,2) * t170 - Icges(7,6) * t232;
t117 = Icges(7,5) * t171 + Icges(7,6) * t170 - Icges(7,3) * t232;
t116 = V_base(5) * rSges(3,1) + (-t207 - t208) * t222 + t257;
t115 = t222 * t212 + (-rSges(3,1) - pkin(6)) * V_base(4) + t247;
t113 = t208 * V_base(4) + (-t211 - t212) * V_base(5) + t262;
t112 = t210 * t216 + (-t169 + t258) * t222 + t256;
t111 = t222 * t166 - t217 * t210 + t244;
t110 = t141 * rSges(7,1) + t140 * rSges(7,2) - rSges(7,3) * t265;
t109 = rSges(7,1) * t139 + rSges(7,2) * t138 + rSges(7,3) * t266;
t108 = Icges(7,1) * t141 + Icges(7,4) * t140 - Icges(7,5) * t265;
t107 = Icges(7,1) * t139 + Icges(7,4) * t138 + Icges(7,5) * t266;
t106 = Icges(7,4) * t141 + Icges(7,2) * t140 - Icges(7,6) * t265;
t105 = Icges(7,4) * t139 + Icges(7,2) * t138 + Icges(7,6) * t266;
t104 = Icges(7,5) * t141 + Icges(7,6) * t140 - Icges(7,3) * t265;
t103 = Icges(7,5) * t139 + Icges(7,6) * t138 + Icges(7,3) * t266;
t102 = -t216 * t166 + t217 * t169 + t245;
t101 = -t137 * t206 + t168 * t175 + t243;
t100 = t206 * t135 - t176 * t168 + t240;
t99 = -t175 * t135 + t176 * t137 + t242;
t98 = t167 * t175 + (-t136 - t143) * t206 + t241;
t97 = t206 * t134 + (-t167 - t182) * t176 + t238;
t96 = t176 * t136 + (-t134 - t142) * t175 + t239;
t95 = -t110 * t185 + t121 * t151 + t175 * t187 + (-t143 - t149) * t206 + t241;
t94 = t185 * t109 - t152 * t121 + t206 * t148 + (-t182 - t187) * t176 + t238;
t93 = -t151 * t109 + t152 * t110 + t176 * t149 + (-t142 - t148) * t175 + t239;
t1 = t185 * ((-t103 * t232 + t105 * t170 + t107 * t171) * t152 + (-t104 * t232 + t106 * t170 + t108 * t171) * t151 + (-t232 * t117 + t170 * t118 + t171 * t119) * t185) / 0.2e1 + m(1) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(2) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(3) * (t113 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(6) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + t152 * ((t103 * t266 + t138 * t105 + t139 * t107) * t152 + (t104 * t266 + t106 * t138 + t108 * t139) * t151 + (t117 * t266 + t118 * t138 + t119 * t139) * t185) / 0.2e1 + t151 * ((-t103 * t265 + t140 * t105 + t141 * t107) * t152 + (-t104 * t265 + t140 * t106 + t141 * t108) * t151 + (-t117 * t265 + t140 * t118 + t141 * t119) * t185) / 0.2e1 + t217 * (t281 * t233 + t246 * t237) / 0.2e1 + t216 * (t246 * t233 - t281 * t237) / 0.2e1 + ((-t179 * t290 + t180 * t288 + t265 * t289) * t206 + (-t179 * t296 + t292 * t180 + t294 * t265) * t176 + (-t295 * t179 + t291 * t180 + t293 * t265) * t175) * t175 / 0.2e1 + ((t177 * t290 + t178 * t288 - t266 * t289) * t206 + (t296 * t177 + t292 * t178 - t294 * t266) * t176 + (t177 * t295 + t178 * t291 - t266 * t293) * t175) * t176 / 0.2e1 + (((t231 * t290 + t235 * t288) * t206 + (t231 * t296 + t292 * t235) * t176 + (t231 * t295 + t235 * t291) * t175) * t236 + (t175 * t293 + t176 * t294 + t206 * t289) * t232) * t206 / 0.2e1 + ((-t160 * t232 + t164 * t236) * t217 + (-t161 * t232 + t165 * t236) * t216 + (-t232 * t198 + t236 * t203 + Icges(3,1) + Icges(2,3)) * t222) * t222 / 0.2e1 + ((t233 * t287 + t237 * t285 + Icges(1,4)) * V_base(5) + (t286 * t233 + t284 * t237 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t285 * t233 - t287 * t237 + Icges(1,2)) * V_base(5) + (t233 * t284 - t237 * t286 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t222 * (t233 * t306 - t237 * t303) + V_base(4) * t222 * (t233 * t303 + t306 * t237) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
