% Calculate kinetic energy for
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:10
% EndTime: 2019-03-09 02:41:13
% DurationCPUTime: 3.16s
% Computational Cost: add. (1807->296), mult. (1464->394), div. (0->0), fcn. (1244->10), ass. (0->154)
t310 = Icges(5,4) + Icges(6,6);
t309 = Icges(5,1) + Icges(6,2);
t308 = -Icges(5,2) - Icges(6,3);
t217 = qJ(3) + pkin(10);
t211 = cos(t217);
t307 = t310 * t211;
t209 = sin(t217);
t306 = t310 * t209;
t305 = Icges(6,4) - Icges(5,5);
t304 = Icges(6,5) - Icges(5,6);
t303 = t308 * t209 + t307;
t302 = t309 * t211 - t306;
t218 = qJ(1) + pkin(9);
t210 = sin(t218);
t212 = cos(t218);
t301 = t210 * t303 + t212 * t304;
t300 = -t210 * t304 + t212 * t303;
t299 = t302 * t210 + t212 * t305;
t298 = -t210 * t305 + t302 * t212;
t297 = t308 * t211 - t306;
t296 = t309 * t209 + t307;
t295 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t221 = sin(qJ(3));
t224 = cos(qJ(3));
t294 = Icges(4,5) * t224 - Icges(4,6) * t221 + t304 * t209 - t211 * t305;
t276 = Icges(4,4) * t224;
t247 = -Icges(4,2) * t221 + t276;
t140 = -Icges(4,6) * t212 + t210 * t247;
t141 = Icges(4,6) * t210 + t212 * t247;
t277 = Icges(4,4) * t221;
t249 = Icges(4,1) * t224 - t277;
t143 = -Icges(4,5) * t212 + t210 * t249;
t144 = Icges(4,5) * t210 + t212 * t249;
t189 = -qJD(3) * t212 + V_base(5);
t190 = qJD(3) * t210 + V_base(4);
t194 = Icges(4,2) * t224 + t277;
t197 = Icges(4,1) * t221 + t276;
t213 = V_base(6) + qJD(1);
t291 = (-t194 * t221 + t197 * t224 + t297 * t209 + t296 * t211) * t213 + (-t141 * t221 + t144 * t224 - t300 * t209 + t298 * t211) * t190 + (-t140 * t221 + t143 * t224 - t301 * t209 + t299 * t211) * t189;
t290 = (Icges(4,5) * t221 + Icges(4,6) * t224 - t209 * t305 - t304 * t211) * t213 + (t295 * t210 + t294 * t212) * t190 + (t294 * t210 - t295 * t212) * t189;
t222 = sin(qJ(1));
t286 = pkin(1) * t222;
t225 = cos(qJ(1));
t285 = pkin(1) * t225;
t284 = pkin(3) * t221;
t283 = pkin(8) * t209;
t282 = pkin(3) * t224;
t281 = -pkin(6) - qJ(2);
t279 = Icges(2,4) * t222;
t278 = Icges(3,4) * t210;
t271 = t210 * t211;
t220 = sin(qJ(6));
t270 = t210 * t220;
t223 = cos(qJ(6));
t269 = t210 * t223;
t268 = t211 * t212;
t267 = t212 * t220;
t266 = t212 * t223;
t117 = qJ(4) * t210 + t212 * t282;
t250 = pkin(4) * t211 + qJ(5) * t209;
t152 = t250 * t212;
t265 = -t117 - t152;
t264 = qJD(5) * t209;
t263 = qJD(6) * t211;
t262 = t213 * t285 + V_base(2);
t261 = V_base(5) * pkin(6) + V_base(1);
t181 = pkin(2) * t210 - pkin(7) * t212;
t258 = -t181 - t286;
t176 = pkin(4) * t209 - qJ(5) * t211;
t257 = -t176 - t284;
t256 = V_base(5) * qJ(2) + t261;
t255 = V_base(4) * t286 + qJD(2) + V_base(3);
t116 = -qJ(4) * t212 + t210 * t282;
t254 = -t116 + t258;
t253 = rSges(4,1) * t224 - rSges(4,2) * t221;
t252 = rSges(5,1) * t211 - rSges(5,2) * t209;
t251 = -rSges(6,2) * t211 + rSges(6,3) * t209;
t151 = t250 * t210;
t240 = -t151 + t254;
t239 = qJD(4) * t210 + t189 * t284 + t256;
t238 = t189 * t176 + t212 * t264 + t239;
t182 = pkin(2) * t212 + pkin(7) * t210;
t234 = t213 * t182 + t281 * V_base(4) + t262;
t233 = V_base(4) * t181 + (-t182 - t285) * V_base(5) + t255;
t232 = t190 * t116 + t233;
t231 = -qJD(4) * t212 + t213 * t117 + t234;
t230 = t213 * t152 + t210 * t264 + t231;
t229 = -qJD(5) * t211 + t190 * t151 + t232;
t215 = Icges(2,4) * t225;
t207 = Icges(3,4) * t212;
t202 = rSges(2,1) * t225 - t222 * rSges(2,2);
t201 = t222 * rSges(2,1) + rSges(2,2) * t225;
t200 = rSges(4,1) * t221 + rSges(4,2) * t224;
t199 = Icges(2,1) * t225 - t279;
t198 = Icges(2,1) * t222 + t215;
t196 = -Icges(2,2) * t222 + t215;
t195 = Icges(2,2) * t225 + t279;
t188 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t187 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t186 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t185 = qJD(6) * t209 + t213;
t180 = rSges(3,1) * t212 - rSges(3,2) * t210;
t179 = rSges(3,1) * t210 + rSges(3,2) * t212;
t178 = rSges(5,1) * t209 + rSges(5,2) * t211;
t177 = -rSges(6,2) * t209 - rSges(6,3) * t211;
t175 = Icges(3,1) * t212 - t278;
t174 = Icges(3,1) * t210 + t207;
t172 = -Icges(3,2) * t210 + t207;
t171 = Icges(3,2) * t212 + t278;
t161 = -pkin(5) * t212 + pkin(8) * t271;
t160 = pkin(5) * t210 + pkin(8) * t268;
t158 = t209 * t270 - t266;
t157 = t209 * t269 + t267;
t156 = t209 * t267 + t269;
t155 = t209 * t266 - t270;
t154 = t212 * t263 + t190;
t153 = t210 * t263 + t189;
t149 = rSges(4,3) * t210 + t212 * t253;
t148 = rSges(7,3) * t209 + (-rSges(7,1) * t220 - rSges(7,2) * t223) * t211;
t147 = -rSges(4,3) * t212 + t210 * t253;
t146 = V_base(5) * rSges(2,3) - t201 * t213 + t261;
t145 = t202 * t213 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t142 = Icges(7,5) * t209 + (-Icges(7,1) * t220 - Icges(7,4) * t223) * t211;
t139 = Icges(7,6) * t209 + (-Icges(7,4) * t220 - Icges(7,2) * t223) * t211;
t136 = Icges(7,3) * t209 + (-Icges(7,5) * t220 - Icges(7,6) * t223) * t211;
t135 = t201 * V_base(4) - t202 * V_base(5) + V_base(3);
t133 = -rSges(6,1) * t212 + t210 * t251;
t132 = rSges(6,1) * t210 + t212 * t251;
t131 = rSges(5,3) * t210 + t212 * t252;
t130 = -rSges(5,3) * t212 + t210 * t252;
t113 = V_base(5) * rSges(3,3) + (-t179 - t286) * t213 + t256;
t112 = t180 * t213 + (-rSges(3,3) + t281) * V_base(4) + t262;
t110 = V_base(4) * t179 + (-t180 - t285) * V_base(5) + t255;
t109 = rSges(7,1) * t158 + rSges(7,2) * t157 + rSges(7,3) * t271;
t108 = rSges(7,1) * t156 + rSges(7,2) * t155 + rSges(7,3) * t268;
t107 = Icges(7,1) * t158 + Icges(7,4) * t157 + Icges(7,5) * t271;
t106 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t268;
t105 = Icges(7,4) * t158 + Icges(7,2) * t157 + Icges(7,6) * t271;
t104 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t268;
t103 = Icges(7,5) * t158 + Icges(7,6) * t157 + Icges(7,3) * t271;
t102 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t268;
t101 = t189 * t200 + (-t147 + t258) * t213 + t256;
t100 = t149 * t213 - t190 * t200 + t234;
t99 = t190 * t147 - t189 * t149 + t233;
t98 = t178 * t189 + (-t130 + t254) * t213 + t239;
t97 = t131 * t213 + (-t178 - t284) * t190 + t231;
t96 = t190 * t130 + (-t117 - t131) * t189 + t232;
t95 = t177 * t189 + (-t133 + t240) * t213 + t238;
t94 = t132 * t213 + (-t177 + t257) * t190 + t230;
t93 = t190 * t133 + (-t132 + t265) * t189 + t229;
t92 = t189 * t283 - t109 * t185 + t148 * t153 + (-t161 + t240) * t213 + t238;
t91 = t108 * t185 - t148 * t154 + t160 * t213 + (t257 - t283) * t190 + t230;
t90 = -t153 * t108 + t154 * t109 + t190 * t161 + (-t160 + t265) * t189 + t229;
t1 = m(1) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(2) * (t135 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t110 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + t185 * ((t102 * t154 + t103 * t153 + t136 * t185) * t209 + ((-t104 * t223 - t106 * t220) * t154 + (-t105 * t223 - t107 * t220) * t153 + (-t139 * t223 - t142 * t220) * t185) * t211) / 0.2e1 + t154 * ((t102 * t268 + t155 * t104 + t156 * t106) * t154 + (t103 * t268 + t105 * t155 + t107 * t156) * t153 + (t136 * t268 + t139 * t155 + t142 * t156) * t185) / 0.2e1 + t153 * ((t102 * t271 + t104 * t157 + t106 * t158) * t154 + (t103 * t271 + t157 * t105 + t158 * t107) * t153 + (t136 * t271 + t139 * t157 + t142 * t158) * t185) / 0.2e1 + (t291 * t210 - t290 * t212) * t189 / 0.2e1 + (t290 * t210 + t291 * t212) * t190 / 0.2e1 + ((-t171 * t210 + t174 * t212 - t222 * t195 + t198 * t225 + Icges(1,4)) * V_base(5) + (-t172 * t210 + t175 * t212 - t222 * t196 + t199 * t225 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t171 * t212 + t174 * t210 + t195 * t225 + t222 * t198 + Icges(1,2)) * V_base(5) + (t172 * t212 + t175 * t210 + t196 * t225 + t222 * t199 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t141 * t224 + t144 * t221 + t298 * t209 + t300 * t211) * t190 + (t140 * t224 + t143 * t221 + t299 * t209 + t301 * t211) * t189 + (t194 * t224 + t197 * t221 + t296 * t209 - t297 * t211 + Icges(2,3) + Icges(3,3)) * t213) * t213 / 0.2e1 + t213 * V_base(5) * (Icges(2,5) * t222 + Icges(3,5) * t210 + Icges(2,6) * t225 + Icges(3,6) * t212) + t213 * V_base(4) * (Icges(2,5) * t225 + Icges(3,5) * t212 - Icges(2,6) * t222 - Icges(3,6) * t210) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
