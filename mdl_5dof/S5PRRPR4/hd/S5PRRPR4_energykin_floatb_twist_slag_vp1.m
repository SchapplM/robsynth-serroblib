% Calculate kinetic energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:23
% EndTime: 2019-12-05 16:21:26
% DurationCPUTime: 3.19s
% Computational Cost: add. (1445->325), mult. (1974->476), div. (0->0), fcn. (1924->10), ass. (0->148)
t271 = -Icges(5,3) - Icges(4,3);
t213 = qJ(3) + pkin(9);
t206 = sin(t213);
t207 = cos(t213);
t215 = cos(pkin(8));
t214 = sin(pkin(8));
t220 = cos(qJ(2));
t254 = t214 * t220;
t152 = -t206 * t254 - t215 * t207;
t153 = -t215 * t206 + t207 * t254;
t219 = cos(qJ(3));
t217 = sin(qJ(3));
t250 = t217 * t220;
t166 = -t214 * t250 - t215 * t219;
t249 = t219 * t220;
t253 = t215 * t217;
t167 = t214 * t249 - t253;
t218 = sin(qJ(2));
t255 = t214 * t218;
t268 = Icges(4,5) * t167 + Icges(5,5) * t153 + Icges(4,6) * t166 + Icges(5,6) * t152 - t255 * t271;
t251 = t215 * t220;
t154 = -t206 * t251 + t214 * t207;
t155 = t214 * t206 + t207 * t251;
t168 = t214 * t219 - t215 * t250;
t256 = t214 * t217;
t169 = t215 * t249 + t256;
t252 = t215 * t218;
t267 = Icges(4,5) * t169 + Icges(5,5) * t155 + Icges(4,6) * t168 + Icges(5,6) * t154 - t252 * t271;
t266 = t271 * t220 + (Icges(4,5) * t219 + Icges(5,5) * t207 - Icges(4,6) * t217 - Icges(5,6) * t206) * t218;
t261 = t219 * pkin(3);
t259 = Icges(2,4) * t214;
t258 = Icges(3,4) * t218;
t257 = Icges(3,4) * t220;
t248 = pkin(4) * t207;
t246 = qJD(3) * t218;
t245 = qJD(4) * t218;
t244 = qJD(5) * t218;
t243 = V_base(5) * qJ(1) + V_base(1);
t239 = qJD(1) + V_base(3);
t196 = qJD(2) * t214 + V_base(4);
t238 = pkin(4) * t206;
t165 = t215 * t246 + t196;
t237 = pkin(2) * t220 + pkin(6) * t218;
t195 = -qJD(2) * t215 + V_base(5);
t236 = rSges(3,1) * t220 - rSges(3,2) * t218;
t235 = Icges(3,1) * t220 - t258;
t234 = -Icges(3,2) * t218 + t257;
t233 = Icges(3,5) * t220 - Icges(3,6) * t218;
t164 = t214 * t246 + t195;
t189 = pkin(1) * t215 + pkin(5) * t214;
t232 = -V_base(4) * qJ(1) + V_base(6) * t189 + V_base(2);
t188 = pkin(1) * t214 - pkin(5) * t215;
t231 = V_base(4) * t188 - V_base(5) * t189 + t239;
t230 = qJ(4) * t218 + t220 * t261;
t229 = pkin(7) * t218 + t220 * t248;
t170 = t237 * t214;
t194 = t218 * pkin(2) - pkin(6) * t220;
t228 = t195 * t194 + (-t170 - t188) * V_base(6) + t243;
t227 = (-Icges(3,3) * t215 + t214 * t233) * t195 + (Icges(3,3) * t214 + t215 * t233) * t196 + (Icges(3,5) * t218 + Icges(3,6) * t220) * V_base(6);
t171 = t237 * t215;
t226 = V_base(6) * t171 - t194 * t196 + t232;
t134 = -qJ(4) * t220 + t218 * t261;
t225 = t164 * t134 + t215 * t245 + t228;
t224 = t196 * t170 - t195 * t171 + t231;
t122 = pkin(3) * t256 + t215 * t230;
t197 = -qJD(3) * t220 + V_base(6);
t223 = t197 * t122 + t214 * t245 + t226;
t121 = -pkin(3) * t253 + t214 * t230;
t222 = -qJD(4) * t220 + t165 * t121 + t224;
t148 = -Icges(3,6) * t215 + t214 * t234;
t149 = Icges(3,6) * t214 + t215 * t234;
t150 = -Icges(3,5) * t215 + t214 * t235;
t151 = Icges(3,5) * t214 + t215 * t235;
t191 = Icges(3,2) * t220 + t258;
t192 = Icges(3,1) * t218 + t257;
t221 = (-t149 * t218 + t151 * t220) * t196 + (-t148 * t218 + t150 * t220) * t195 + (-t191 * t218 + t192 * t220) * V_base(6);
t209 = qJ(5) + t213;
t208 = Icges(2,4) * t215;
t203 = cos(t209);
t202 = sin(t209);
t193 = t218 * rSges(3,1) + rSges(3,2) * t220;
t187 = rSges(2,1) * t215 - rSges(2,2) * t214;
t186 = rSges(2,1) * t214 + rSges(2,2) * t215;
t185 = Icges(2,1) * t215 - t259;
t184 = Icges(2,1) * t214 + t208;
t183 = -Icges(2,2) * t214 + t208;
t182 = Icges(2,2) * t215 + t259;
t178 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t177 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t176 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t174 = V_base(6) + (-qJD(3) - qJD(5)) * t220;
t161 = -rSges(4,3) * t220 + (rSges(4,1) * t219 - rSges(4,2) * t217) * t218;
t160 = -Icges(4,5) * t220 + (Icges(4,1) * t219 - Icges(4,4) * t217) * t218;
t159 = -Icges(4,6) * t220 + (Icges(4,4) * t219 - Icges(4,2) * t217) * t218;
t157 = t214 * rSges(3,3) + t215 * t236;
t156 = -t215 * rSges(3,3) + t214 * t236;
t145 = t214 * t202 + t203 * t251;
t144 = -t202 * t251 + t214 * t203;
t143 = -t215 * t202 + t203 * t254;
t142 = -t202 * t254 - t215 * t203;
t141 = -rSges(5,3) * t220 + (rSges(5,1) * t207 - rSges(5,2) * t206) * t218;
t140 = t215 * t244 + t165;
t139 = t214 * t244 + t164;
t137 = -Icges(5,5) * t220 + (Icges(5,1) * t207 - Icges(5,4) * t206) * t218;
t136 = -Icges(5,6) * t220 + (Icges(5,4) * t207 - Icges(5,2) * t206) * t218;
t133 = V_base(5) * rSges(2,3) - t186 * V_base(6) + t243;
t132 = t187 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t131 = -rSges(6,3) * t220 + (rSges(6,1) * t203 - rSges(6,2) * t202) * t218;
t130 = -Icges(6,5) * t220 + (Icges(6,1) * t203 - Icges(6,4) * t202) * t218;
t129 = -Icges(6,6) * t220 + (Icges(6,4) * t203 - Icges(6,2) * t202) * t218;
t128 = -Icges(6,3) * t220 + (Icges(6,5) * t203 - Icges(6,6) * t202) * t218;
t127 = t186 * V_base(4) - t187 * V_base(5) + t239;
t126 = -pkin(7) * t220 + t218 * t248;
t124 = rSges(4,1) * t169 + rSges(4,2) * t168 + rSges(4,3) * t252;
t123 = rSges(4,1) * t167 + rSges(4,2) * t166 + rSges(4,3) * t255;
t120 = Icges(4,1) * t169 + Icges(4,4) * t168 + Icges(4,5) * t252;
t119 = Icges(4,1) * t167 + Icges(4,4) * t166 + Icges(4,5) * t255;
t118 = Icges(4,4) * t169 + Icges(4,2) * t168 + Icges(4,6) * t252;
t117 = Icges(4,4) * t167 + Icges(4,2) * t166 + Icges(4,6) * t255;
t113 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t252;
t112 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t255;
t111 = Icges(5,1) * t155 + Icges(5,4) * t154 + Icges(5,5) * t252;
t110 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t255;
t109 = Icges(5,4) * t155 + Icges(5,2) * t154 + Icges(5,6) * t252;
t108 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t255;
t105 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t252;
t104 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t255;
t103 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t252;
t102 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t255;
t101 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t252;
t100 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t255;
t99 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t252;
t98 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t255;
t96 = t193 * t195 + (-t156 - t188) * V_base(6) + t243;
t95 = t157 * V_base(6) - t193 * t196 + t232;
t94 = t214 * t238 + t215 * t229;
t93 = t214 * t229 - t215 * t238;
t92 = t156 * t196 - t157 * t195 + t231;
t91 = -t123 * t197 + t161 * t164 + t228;
t90 = t124 * t197 - t161 * t165 + t226;
t89 = t123 * t165 - t124 * t164 + t224;
t88 = t141 * t164 + (-t112 - t121) * t197 + t225;
t87 = t113 * t197 + (-t134 - t141) * t165 + t223;
t86 = t165 * t112 + (-t113 - t122) * t164 + t222;
t85 = -t104 * t174 + t126 * t164 + t131 * t139 + (-t121 - t93) * t197 + t225;
t84 = t105 * t174 - t131 * t140 + t197 * t94 + (-t126 - t134) * t165 + t223;
t83 = t140 * t104 - t139 * t105 + t165 * t93 + (-t122 - t94) * t164 + t222;
t1 = m(1) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(2) * (t127 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t196 * (t227 * t214 + t221 * t215) / 0.2e1 + t195 * (t221 * t214 - t227 * t215) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t140 * ((t144 * t101 + t145 * t103 + t99 * t252) * t140 + (t100 * t144 + t102 * t145 + t252 * t98) * t139 + (t128 * t252 + t129 * t144 + t130 * t145) * t174) / 0.2e1 + t139 * ((t101 * t142 + t103 * t143 + t255 * t99) * t140 + (t142 * t100 + t143 * t102 + t98 * t255) * t139 + (t128 * t255 + t129 * t142 + t130 * t143) * t174) / 0.2e1 + t174 * ((-t128 * t174 - t98 * t139 - t99 * t140) * t220 + ((-t101 * t202 + t103 * t203) * t140 + (-t100 * t202 + t102 * t203) * t139 + (-t129 * t202 + t130 * t203) * t174) * t218) / 0.2e1 + ((t136 * t152 + t137 * t153 + t159 * t166 + t160 * t167 + t266 * t255) * t197 + (t109 * t152 + t111 * t153 + t118 * t166 + t120 * t167 + t267 * t255) * t165 + (t152 * t108 + t153 * t110 + t166 * t117 + t167 * t119 + t268 * t255) * t164) * t164 / 0.2e1 + ((t136 * t154 + t137 * t155 + t159 * t168 + t160 * t169 + t266 * t252) * t197 + (t154 * t109 + t155 * t111 + t168 * t118 + t169 * t120 + t267 * t252) * t165 + (t108 * t154 + t110 * t155 + t117 * t168 + t119 * t169 + t268 * t252) * t164) * t165 / 0.2e1 + ((-t268 * t164 - t267 * t165 - t266 * t197) * t220 + ((-t136 * t206 + t137 * t207 - t159 * t217 + t160 * t219) * t197 + (-t109 * t206 + t111 * t207 - t118 * t217 + t120 * t219) * t165 + (-t108 * t206 + t110 * t207 - t117 * t217 + t119 * t219) * t164) * t218) * t197 / 0.2e1 + ((-t182 * t214 + t184 * t215 + Icges(1,4)) * V_base(5) + (-t214 * t183 + t215 * t185 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t215 * t182 + t214 * t184 + Icges(1,2)) * V_base(5) + (t183 * t215 + t185 * t214 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t149 * t220 + t218 * t151) * t196 + (t148 * t220 + t218 * t150) * t195 + (t220 * t191 + t218 * t192 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t215 - Icges(2,6) * t214 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t214 + Icges(2,6) * t215 + Icges(1,6));
T = t1;
