% Calculate kinetic energy for
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:11
% EndTime: 2019-03-09 01:46:14
% DurationCPUTime: 2.78s
% Computational Cost: add. (1467->294), mult. (2214->393), div. (0->0), fcn. (2466->10), ass. (0->147)
t292 = Icges(2,4) - Icges(3,5);
t291 = Icges(2,1) + Icges(3,1);
t290 = Icges(3,4) + Icges(2,5);
t289 = Icges(2,2) + Icges(3,3);
t288 = Icges(2,6) - Icges(3,6);
t287 = Icges(5,3) + Icges(6,3);
t215 = qJ(4) + pkin(10);
t207 = sin(t215);
t208 = cos(t215);
t218 = sin(qJ(4));
t220 = cos(qJ(4));
t286 = -Icges(5,5) * t220 - Icges(6,5) * t208 + Icges(5,6) * t218 + Icges(6,6) * t207;
t271 = sin(qJ(1));
t285 = t292 * t271;
t272 = cos(qJ(1));
t284 = t292 * t272;
t283 = -t272 * t289 - t285;
t282 = t271 * t289 - t284;
t281 = t271 * t291 + t284;
t280 = t272 * t291 - t285;
t266 = sin(pkin(9));
t267 = cos(pkin(9));
t172 = -t266 * t271 - t267 * t272;
t173 = t266 * t272 - t267 * t271;
t261 = Icges(6,4) * t208;
t237 = Icges(6,2) * t207 - t261;
t114 = -Icges(6,6) * t172 + t173 * t237;
t115 = Icges(6,6) * t173 + t172 * t237;
t262 = Icges(6,4) * t207;
t239 = -Icges(6,1) * t208 + t262;
t116 = -Icges(6,5) * t172 + t173 * t239;
t117 = Icges(6,5) * t173 + t172 * t239;
t263 = Icges(5,4) * t220;
t238 = Icges(5,2) * t218 - t263;
t123 = -Icges(5,6) * t172 + t173 * t238;
t124 = Icges(5,6) * t173 + t172 * t238;
t264 = Icges(5,4) * t218;
t240 = -Icges(5,1) * t220 + t264;
t125 = -Icges(5,5) * t172 + t173 * t240;
t126 = Icges(5,5) * t173 + t172 * t240;
t159 = -qJD(4) * t172 + V_base(5);
t160 = qJD(4) * t173 + V_base(4);
t168 = -Icges(6,2) * t208 - t262;
t169 = -Icges(6,1) * t207 - t261;
t184 = -Icges(5,2) * t220 - t264;
t189 = -Icges(5,1) * t218 - t263;
t209 = V_base(6) + qJD(1);
t279 = (t168 * t207 - t169 * t208 + t184 * t218 - t189 * t220) * t209 + (t115 * t207 - t117 * t208 + t124 * t218 - t126 * t220) * t160 + (t114 * t207 - t116 * t208 + t123 * t218 - t125 * t220) * t159;
t278 = (-Icges(5,5) * t218 - Icges(6,5) * t207 - Icges(5,6) * t220 - Icges(6,6) * t208) * t209 + (t172 * t286 + t173 * t287) * t160 + (-t172 * t287 + t173 * t286) * t159;
t270 = pkin(4) * t218;
t269 = pkin(4) * t220;
t265 = Icges(4,4) * t172;
t260 = t172 * t207;
t259 = t173 * t207;
t217 = sin(qJ(6));
t258 = t208 * t217;
t219 = cos(qJ(6));
t257 = t208 * t219;
t256 = qJD(6) * t207;
t195 = pkin(1) * t271 - qJ(2) * t272;
t255 = V_base(4) * t195 + V_base(3);
t254 = V_base(5) * pkin(6) + V_base(1);
t251 = t272 * pkin(2);
t250 = t271 * pkin(2);
t247 = qJD(2) * t271 + t254;
t246 = -t195 - t250;
t198 = pkin(1) * t272 + qJ(2) * t271;
t245 = -t198 - t251;
t244 = -pkin(5) * t208 - pkin(8) * t207;
t243 = V_base(4) * t250 - qJD(3) + t255;
t242 = -rSges(5,1) * t220 + rSges(5,2) * t218;
t241 = -rSges(6,1) * t208 + rSges(6,2) * t207;
t149 = -pkin(3) * t173 - pkin(7) * t172;
t234 = -t149 + t246;
t233 = -qJD(2) * t272 + t209 * t198 + V_base(2);
t110 = -qJ(5) * t172 - t173 * t269;
t232 = -t110 + t234;
t231 = -V_base(5) * qJ(3) + t247;
t228 = qJD(5) * t173 + t231;
t227 = V_base(4) * qJ(3) + t209 * t251 + t233;
t150 = -pkin(3) * t172 + pkin(7) * t173;
t226 = -V_base(4) * pkin(6) + t209 * t150 + t227;
t225 = V_base(4) * t149 + (-t150 + t245) * V_base(5) + t243;
t224 = t160 * t110 + t225;
t111 = qJ(5) * t173 - t172 * t269;
t223 = -qJD(5) * t172 + t209 * t111 + t160 * t270 + t226;
t200 = rSges(2,1) * t272 - rSges(2,2) * t271;
t199 = rSges(3,1) * t272 + rSges(3,3) * t271;
t197 = rSges(2,1) * t271 + rSges(2,2) * t272;
t196 = rSges(3,1) * t271 - rSges(3,3) * t272;
t194 = -t218 * rSges(5,1) - rSges(5,2) * t220;
t178 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t177 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t176 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t175 = qJD(6) * t208 + t209;
t171 = -pkin(5) * t207 + pkin(8) * t208;
t170 = -rSges(6,1) * t207 - rSges(6,2) * t208;
t165 = Icges(4,4) * t173;
t157 = rSges(7,3) * t208 + (-rSges(7,1) * t219 + rSges(7,2) * t217) * t207;
t156 = V_base(5) * rSges(2,3) - t197 * t209 + t254;
t155 = t200 * t209 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t154 = Icges(7,5) * t208 + (-Icges(7,1) * t219 + Icges(7,4) * t217) * t207;
t153 = Icges(7,6) * t208 + (-Icges(7,4) * t219 + Icges(7,2) * t217) * t207;
t152 = Icges(7,3) * t208 + (-Icges(7,5) * t219 + Icges(7,6) * t217) * t207;
t151 = t197 * V_base(4) - t200 * V_base(5) + V_base(3);
t148 = -rSges(4,1) * t172 - rSges(4,2) * t173;
t147 = -rSges(4,1) * t173 + rSges(4,2) * t172;
t146 = -Icges(4,1) * t172 - t165;
t145 = -Icges(4,1) * t173 + t265;
t144 = -Icges(4,2) * t173 - t265;
t143 = Icges(4,2) * t172 - t165;
t138 = -t172 * t257 + t173 * t217;
t137 = t172 * t258 + t173 * t219;
t136 = -t172 * t217 - t173 * t257;
t135 = -t172 * t219 + t173 * t258;
t134 = -t172 * t256 + t160;
t133 = -t173 * t256 + t159;
t132 = t244 * t172;
t131 = t244 * t173;
t130 = V_base(5) * rSges(3,2) + (-t195 - t196) * t209 + t247;
t129 = t209 * t199 + (-rSges(3,2) - pkin(6)) * V_base(4) + t233;
t128 = t173 * rSges(5,3) + t172 * t242;
t127 = -t172 * rSges(5,3) + t173 * t242;
t120 = t196 * V_base(4) + (-t198 - t199) * V_base(5) + t255;
t119 = rSges(6,3) * t173 + t172 * t241;
t118 = -rSges(6,3) * t172 + t173 * t241;
t108 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t147 + t246) * t209 + t247;
t107 = t209 * t148 + (rSges(4,3) - pkin(6)) * V_base(4) + t227;
t105 = rSges(7,1) * t138 + rSges(7,2) * t137 - rSges(7,3) * t260;
t104 = rSges(7,1) * t136 + rSges(7,2) * t135 - rSges(7,3) * t259;
t103 = Icges(7,1) * t138 + Icges(7,4) * t137 - Icges(7,5) * t260;
t102 = Icges(7,1) * t136 + Icges(7,4) * t135 - Icges(7,5) * t259;
t101 = Icges(7,4) * t138 + Icges(7,2) * t137 - Icges(7,6) * t260;
t100 = Icges(7,4) * t136 + Icges(7,2) * t135 - Icges(7,6) * t259;
t99 = Icges(7,5) * t138 + Icges(7,6) * t137 - Icges(7,3) * t260;
t98 = Icges(7,5) * t136 + Icges(7,6) * t135 - Icges(7,3) * t259;
t97 = V_base(4) * t147 + (-t148 + t245) * V_base(5) + t243;
t96 = t159 * t194 + (-t127 + t234) * t209 + t231;
t95 = t209 * t128 - t160 * t194 + t226;
t94 = t160 * t127 - t159 * t128 + t225;
t93 = (t170 - t270) * t159 + (-t118 + t232) * t209 + t228;
t92 = t209 * t119 - t160 * t170 + t223;
t91 = t160 * t118 + (-t111 - t119) * t159 + t224;
t90 = -t175 * t104 + t133 * t157 + (t171 - t270) * t159 + (-t131 + t232) * t209 + t228;
t89 = t175 * t105 + t209 * t132 - t134 * t157 - t160 * t171 + t223;
t88 = t134 * t104 - t133 * t105 + t160 * t131 + (-t111 - t132) * t159 + t224;
t1 = m(3) * (t120 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t133 * ((t101 * t135 + t103 * t136 - t259 * t99) * t134 + (t135 * t100 + t136 * t102 - t98 * t259) * t133 + (t135 * t153 + t136 * t154 - t152 * t259) * t175) / 0.2e1 + t134 * ((t137 * t101 + t138 * t103 - t99 * t260) * t134 + (t100 * t137 + t102 * t138 - t260 * t98) * t133 + (t137 * t153 + t138 * t154 - t152 * t260) * t175) / 0.2e1 + m(2) * (t151 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(1) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + t175 * ((t133 * t98 + t134 * t99 + t152 * t175) * t208 + ((t101 * t217 - t103 * t219) * t134 + (t100 * t217 - t102 * t219) * t133 + (t153 * t217 - t154 * t219) * t175) * t207) / 0.2e1 + (-t278 * t172 + t279 * t173) * t159 / 0.2e1 + (t279 * t172 + t278 * t173) * t160 / 0.2e1 + ((-t143 * t173 - t145 * t172 + t271 * t283 + t281 * t272 + Icges(1,4)) * V_base(5) + (-t144 * t173 - t146 * t172 + t282 * t271 + t280 * t272 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t143 * t172 - t145 * t173 + t281 * t271 - t283 * t272 + Icges(1,2)) * V_base(5) + (t144 * t172 - t146 * t173 + t271 * t280 - t272 * t282 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t115 * t208 - t117 * t207 - t124 * t220 - t218 * t126) * t160 + (-t114 * t208 - t116 * t207 - t123 * t220 - t218 * t125) * t159 + (-t168 * t208 - t169 * t207 - t184 * t220 - t218 * t189 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t209) * t209 / 0.2e1 + t209 * V_base(5) * (Icges(4,5) * t173 - Icges(4,6) * t172 + t271 * t290 + t272 * t288) + t209 * V_base(4) * (Icges(4,5) * t172 + Icges(4,6) * t173 - t271 * t288 + t272 * t290) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
