% Calculate kinetic energy for
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:30
% EndTime: 2019-12-05 17:09:32
% DurationCPUTime: 2.26s
% Computational Cost: add. (1484->307), mult. (1682->466), div. (0->0), fcn. (1576->10), ass. (0->151)
t208 = sin(qJ(2));
t260 = pkin(2) * t208;
t210 = cos(qJ(2));
t259 = pkin(2) * t210;
t209 = cos(qJ(4));
t258 = pkin(4) * t209;
t205 = sin(pkin(9));
t255 = Icges(2,4) * t205;
t254 = Icges(3,4) * t208;
t253 = Icges(3,4) * t210;
t204 = qJ(2) + qJ(3);
t199 = sin(t204);
t252 = Icges(4,4) * t199;
t201 = cos(t204);
t251 = Icges(4,4) * t201;
t250 = t199 * t205;
t206 = cos(pkin(9));
t249 = t199 * t206;
t248 = t201 * t205;
t247 = t201 * t206;
t207 = sin(qJ(4));
t246 = t205 * t207;
t245 = t205 * t209;
t244 = t206 * t207;
t243 = t206 * t209;
t123 = -pkin(6) * t206 + t205 * t259;
t183 = t205 * pkin(1) - t206 * pkin(5);
t242 = -t123 - t183;
t241 = qJD(4) * t199;
t240 = qJD(5) * t199;
t239 = V_base(5) * qJ(1) + V_base(1);
t235 = qJD(1) + V_base(3);
t192 = qJD(2) * t205 + V_base(4);
t191 = -qJD(2) * t206 + V_base(5);
t234 = t191 * t260 + t239;
t165 = qJD(3) * t205 + t192;
t233 = pkin(3) * t201 + pkin(7) * t199;
t232 = rSges(3,1) * t210 - rSges(3,2) * t208;
t231 = rSges(4,1) * t201 - rSges(4,2) * t199;
t140 = t206 * t241 + t165;
t230 = Icges(3,1) * t210 - t254;
t229 = Icges(4,1) * t201 - t252;
t228 = -Icges(3,2) * t208 + t253;
t227 = -Icges(4,2) * t199 + t251;
t226 = Icges(3,5) * t210 - Icges(3,6) * t208;
t225 = Icges(4,5) * t201 - Icges(4,6) * t199;
t184 = t206 * pkin(1) + t205 * pkin(5);
t224 = -V_base(4) * qJ(1) + V_base(6) * t184 + V_base(2);
t164 = V_base(5) + (-qJD(2) - qJD(3)) * t206;
t223 = V_base(4) * t183 - t184 * V_base(5) + t235;
t139 = t205 * t241 + t164;
t222 = pkin(8) * t199 + t201 * t258;
t221 = (-Icges(4,3) * t206 + t205 * t225) * t164 + (Icges(4,3) * t205 + t206 * t225) * t165 + (Icges(4,5) * t199 + Icges(4,6) * t201) * V_base(6);
t220 = (-Icges(3,3) * t206 + t205 * t226) * t191 + (Icges(3,3) * t205 + t206 * t226) * t192 + (Icges(3,5) * t208 + Icges(3,6) * t210) * V_base(6);
t124 = pkin(6) * t205 + t206 * t259;
t219 = V_base(6) * t124 - t192 * t260 + t224;
t154 = t233 * t205;
t169 = pkin(3) * t199 - pkin(7) * t201;
t218 = t164 * t169 + (-t154 + t242) * V_base(6) + t234;
t217 = t192 * t123 - t124 * t191 + t223;
t155 = t233 * t206;
t216 = V_base(6) * t155 - t165 * t169 + t219;
t215 = t165 * t154 - t155 * t164 + t217;
t133 = -Icges(4,6) * t206 + t205 * t227;
t134 = Icges(4,6) * t205 + t206 * t227;
t135 = -Icges(4,5) * t206 + t205 * t229;
t136 = Icges(4,5) * t205 + t206 * t229;
t162 = Icges(4,2) * t201 + t252;
t163 = Icges(4,1) * t199 + t251;
t214 = (-t134 * t199 + t136 * t201) * t165 + (-t133 * t199 + t135 * t201) * t164 + (-t162 * t199 + t163 * t201) * V_base(6);
t148 = -Icges(3,6) * t206 + t205 * t228;
t149 = Icges(3,6) * t205 + t206 * t228;
t150 = -Icges(3,5) * t206 + t205 * t230;
t151 = Icges(3,5) * t205 + t206 * t230;
t186 = Icges(3,2) * t210 + t254;
t187 = Icges(3,1) * t208 + t253;
t213 = (-t149 * t208 + t151 * t210) * t192 + (-t148 * t208 + t150 * t210) * t191 + (-t186 * t208 + t187 * t210) * V_base(6);
t203 = qJ(4) + qJ(5);
t200 = cos(t203);
t198 = sin(t203);
t197 = Icges(2,4) * t206;
t190 = rSges(3,1) * t208 + rSges(3,2) * t210;
t182 = -qJD(4) * t201 + V_base(6);
t181 = rSges(2,1) * t206 - rSges(2,2) * t205;
t180 = rSges(2,1) * t205 + rSges(2,2) * t206;
t179 = Icges(2,1) * t206 - t255;
t178 = Icges(2,1) * t205 + t197;
t177 = -Icges(2,2) * t205 + t197;
t176 = Icges(2,2) * t206 + t255;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t166 = rSges(4,1) * t199 + rSges(4,2) * t201;
t160 = V_base(6) + (-qJD(4) - qJD(5)) * t201;
t159 = t201 * t243 + t246;
t158 = -t201 * t244 + t245;
t157 = t201 * t245 - t244;
t156 = -t201 * t246 - t243;
t153 = rSges(3,3) * t205 + t206 * t232;
t152 = -rSges(3,3) * t206 + t205 * t232;
t145 = t198 * t205 + t200 * t247;
t144 = -t198 * t247 + t200 * t205;
t143 = -t198 * t206 + t200 * t248;
t142 = -t198 * t248 - t200 * t206;
t138 = rSges(4,3) * t205 + t206 * t231;
t137 = -rSges(4,3) * t206 + t205 * t231;
t130 = -rSges(5,3) * t201 + (rSges(5,1) * t209 - rSges(5,2) * t207) * t199;
t129 = V_base(5) * rSges(2,3) - t180 * V_base(6) + t239;
t128 = t181 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t127 = -Icges(5,5) * t201 + (Icges(5,1) * t209 - Icges(5,4) * t207) * t199;
t126 = -Icges(5,6) * t201 + (Icges(5,4) * t209 - Icges(5,2) * t207) * t199;
t125 = -Icges(5,3) * t201 + (Icges(5,5) * t209 - Icges(5,6) * t207) * t199;
t120 = -rSges(6,3) * t201 + (rSges(6,1) * t200 - rSges(6,2) * t198) * t199;
t119 = -Icges(6,5) * t201 + (Icges(6,1) * t200 - Icges(6,4) * t198) * t199;
t118 = -Icges(6,6) * t201 + (Icges(6,4) * t200 - Icges(6,2) * t198) * t199;
t117 = -Icges(6,3) * t201 + (Icges(6,5) * t200 - Icges(6,6) * t198) * t199;
t116 = t180 * V_base(4) - t181 * V_base(5) + t235;
t115 = -pkin(8) * t201 + t199 * t258;
t114 = t206 * t240 + t140;
t113 = t205 * t240 + t139;
t110 = rSges(5,1) * t159 + rSges(5,2) * t158 + rSges(5,3) * t249;
t109 = rSges(5,1) * t157 + rSges(5,2) * t156 + rSges(5,3) * t250;
t108 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t249;
t107 = Icges(5,1) * t157 + Icges(5,4) * t156 + Icges(5,5) * t250;
t106 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t249;
t105 = Icges(5,4) * t157 + Icges(5,2) * t156 + Icges(5,6) * t250;
t104 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t249;
t103 = Icges(5,5) * t157 + Icges(5,6) * t156 + Icges(5,3) * t250;
t102 = pkin(4) * t246 + t206 * t222;
t101 = -pkin(4) * t244 + t205 * t222;
t100 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t249;
t99 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t250;
t98 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t249;
t97 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t250;
t96 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t249;
t95 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t250;
t94 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t249;
t93 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t250;
t92 = t190 * t191 + (-t152 - t183) * V_base(6) + t239;
t91 = t153 * V_base(6) - t190 * t192 + t224;
t90 = t152 * t192 - t153 * t191 + t223;
t89 = t164 * t166 + (-t137 + t242) * V_base(6) + t234;
t88 = t138 * V_base(6) - t165 * t166 + t219;
t87 = t137 * t165 - t138 * t164 + t217;
t86 = -t109 * t182 + t130 * t139 + t218;
t85 = t110 * t182 - t130 * t140 + t216;
t84 = t109 * t140 - t110 * t139 + t215;
t83 = -t101 * t182 + t113 * t120 + t115 * t139 - t160 * t99 + t218;
t82 = t100 * t160 + t102 * t182 - t114 * t120 - t115 * t140 + t216;
t81 = -t100 * t113 + t101 * t140 - t102 * t139 + t114 * t99 + t215;
t1 = m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(3) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + t192 * (t205 * t220 + t206 * t213) / 0.2e1 + t191 * (t205 * t213 - t206 * t220) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t165 * (t205 * t221 + t206 * t214) / 0.2e1 + t164 * (t205 * t214 - t206 * t221) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t140 * ((t104 * t249 + t158 * t106 + t159 * t108) * t140 + (t103 * t249 + t105 * t158 + t107 * t159) * t139 + (t125 * t249 + t126 * t158 + t127 * t159) * t182) / 0.2e1 + t139 * ((t104 * t250 + t106 * t156 + t108 * t157) * t140 + (t103 * t250 + t156 * t105 + t157 * t107) * t139 + (t125 * t250 + t126 * t156 + t127 * t157) * t182) / 0.2e1 + t182 * ((-t103 * t139 - t104 * t140 - t125 * t182) * t201 + ((-t106 * t207 + t108 * t209) * t140 + (-t105 * t207 + t107 * t209) * t139 + (-t126 * t207 + t127 * t209) * t182) * t199) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t114 * ((t144 * t96 + t145 * t98 + t249 * t94) * t114 + (t144 * t95 + t145 * t97 + t249 * t93) * t113 + (t117 * t249 + t118 * t144 + t119 * t145) * t160) / 0.2e1 + t113 * ((t142 * t96 + t143 * t98 + t250 * t94) * t114 + (t142 * t95 + t143 * t97 + t250 * t93) * t113 + (t117 * t250 + t118 * t142 + t119 * t143) * t160) / 0.2e1 + t160 * ((-t93 * t113 - t94 * t114 - t117 * t160) * t201 + ((-t198 * t96 + t200 * t98) * t114 + (-t198 * t95 + t200 * t97) * t113 + (-t118 * t198 + t119 * t200) * t160) * t199) / 0.2e1 + ((-t176 * t205 + t178 * t206 + Icges(1,4)) * V_base(5) + (-t177 * t205 + t179 * t206 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t176 * t206 + t178 * t205 + Icges(1,2)) * V_base(5) + (t177 * t206 + t179 * t205 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t149 * t210 + t151 * t208) * t192 + (t148 * t210 + t150 * t208) * t191 + (t134 * t201 + t136 * t199) * t165 + (t133 * t201 + t135 * t199) * t164 + (t162 * t201 + t163 * t199 + t186 * t210 + t187 * t208 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t206 - Icges(2,6) * t205 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t205 + Icges(2,6) * t206 + Icges(1,6));
T = t1;
