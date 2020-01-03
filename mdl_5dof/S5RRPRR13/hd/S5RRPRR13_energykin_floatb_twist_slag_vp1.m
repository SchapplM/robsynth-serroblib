% Calculate kinetic energy for
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:04
% EndTime: 2019-12-31 20:32:06
% DurationCPUTime: 2.82s
% Computational Cost: add. (1463->333), mult. (1929->501), div. (0->0), fcn. (1879->10), ass. (0->149)
t210 = cos(pkin(9));
t253 = t210 * pkin(3);
t213 = sin(qJ(1));
t252 = Icges(2,4) * t213;
t212 = sin(qJ(2));
t251 = Icges(3,4) * t212;
t214 = cos(qJ(2));
t250 = Icges(3,4) * t214;
t209 = sin(pkin(9));
t215 = cos(qJ(1));
t249 = t209 * t215;
t248 = t212 * t213;
t247 = t212 * t215;
t246 = t213 * t209;
t245 = t213 * t214;
t244 = t214 * t215;
t230 = pkin(2) * t214 + qJ(3) * t212;
t164 = t230 * t213;
t188 = t213 * pkin(1) - pkin(6) * t215;
t242 = -t164 - t188;
t208 = pkin(9) + qJ(4);
t200 = cos(t208);
t241 = pkin(4) * t200;
t239 = qJD(3) * t212;
t238 = qJD(4) * t212;
t237 = qJD(5) * t212;
t236 = V_base(5) * pkin(5) + V_base(1);
t191 = qJD(2) * t213 + V_base(4);
t202 = V_base(6) + qJD(1);
t199 = sin(t208);
t233 = pkin(4) * t199;
t163 = t215 * t238 + t191;
t184 = pkin(2) * t212 - qJ(3) * t214;
t190 = -qJD(2) * t215 + V_base(5);
t232 = t190 * t184 + t215 * t239 + t236;
t231 = rSges(3,1) * t214 - rSges(3,2) * t212;
t229 = Icges(3,1) * t214 - t251;
t228 = -Icges(3,2) * t212 + t250;
t227 = Icges(3,5) * t214 - Icges(3,6) * t212;
t162 = t213 * t238 + t190;
t189 = pkin(1) * t215 + t213 * pkin(6);
t226 = -V_base(4) * pkin(5) + t202 * t189 + V_base(2);
t225 = V_base(4) * t188 - t189 * V_base(5) + V_base(3);
t224 = (-Icges(3,3) * t215 + t213 * t227) * t190 + (Icges(3,3) * t213 + t215 * t227) * t191 + (Icges(3,5) * t212 + Icges(3,6) * t214) * t202;
t223 = pkin(7) * t212 + t214 * t253;
t165 = t230 * t215;
t222 = t202 * t165 + t213 * t239 + t226;
t221 = pkin(8) * t212 + t214 * t241;
t220 = -qJD(3) * t214 + t191 * t164 + t225;
t117 = -pkin(3) * t249 + t213 * t223;
t128 = -pkin(7) * t214 + t212 * t253;
t219 = t190 * t128 + (-t117 + t242) * t202 + t232;
t118 = pkin(3) * t246 + t215 * t223;
t218 = t202 * t118 + (-t128 - t184) * t191 + t222;
t217 = t191 * t117 + (-t118 - t165) * t190 + t220;
t151 = -Icges(3,6) * t215 + t213 * t228;
t152 = Icges(3,6) * t213 + t215 * t228;
t153 = -Icges(3,5) * t215 + t213 * t229;
t154 = Icges(3,5) * t213 + t215 * t229;
t177 = Icges(3,2) * t214 + t251;
t180 = Icges(3,1) * t212 + t250;
t216 = (-t152 * t212 + t154 * t214) * t191 + (-t151 * t212 + t153 * t214) * t190 + (-t177 * t212 + t180 * t214) * t202;
t205 = Icges(2,4) * t215;
t201 = qJ(5) + t208;
t197 = cos(t201);
t196 = sin(t201);
t187 = rSges(2,1) * t215 - t213 * rSges(2,2);
t186 = t213 * rSges(2,1) + rSges(2,2) * t215;
t185 = rSges(3,1) * t212 + rSges(3,2) * t214;
t183 = -qJD(4) * t214 + t202;
t182 = Icges(2,1) * t215 - t252;
t181 = Icges(2,1) * t213 + t205;
t179 = -Icges(2,2) * t213 + t205;
t178 = Icges(2,2) * t215 + t252;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t167 = (-qJD(4) - qJD(5)) * t214 + t202;
t161 = t210 * t244 + t246;
t160 = -t209 * t244 + t213 * t210;
t159 = t210 * t245 - t249;
t158 = -t209 * t245 - t210 * t215;
t156 = t213 * rSges(3,3) + t215 * t231;
t155 = -rSges(3,3) * t215 + t213 * t231;
t148 = t213 * t199 + t200 * t244;
t147 = -t199 * t244 + t213 * t200;
t146 = -t199 * t215 + t200 * t245;
t145 = -t199 * t245 - t200 * t215;
t144 = -rSges(4,3) * t214 + (rSges(4,1) * t210 - rSges(4,2) * t209) * t212;
t142 = -Icges(4,5) * t214 + (Icges(4,1) * t210 - Icges(4,4) * t209) * t212;
t141 = -Icges(4,6) * t214 + (Icges(4,4) * t210 - Icges(4,2) * t209) * t212;
t140 = -Icges(4,3) * t214 + (Icges(4,5) * t210 - Icges(4,6) * t209) * t212;
t139 = t213 * t196 + t197 * t244;
t138 = -t196 * t244 + t213 * t197;
t137 = -t196 * t215 + t197 * t245;
t136 = -t196 * t245 - t197 * t215;
t135 = t215 * t237 + t163;
t134 = t213 * t237 + t162;
t132 = -rSges(5,3) * t214 + (rSges(5,1) * t200 - rSges(5,2) * t199) * t212;
t131 = -Icges(5,5) * t214 + (Icges(5,1) * t200 - Icges(5,4) * t199) * t212;
t130 = -Icges(5,6) * t214 + (Icges(5,4) * t200 - Icges(5,2) * t199) * t212;
t129 = -Icges(5,3) * t214 + (Icges(5,5) * t200 - Icges(5,6) * t199) * t212;
t127 = -rSges(6,3) * t214 + (rSges(6,1) * t197 - rSges(6,2) * t196) * t212;
t126 = V_base(5) * rSges(2,3) - t186 * t202 + t236;
t125 = t187 * t202 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t124 = -Icges(6,5) * t214 + (Icges(6,1) * t197 - Icges(6,4) * t196) * t212;
t123 = -Icges(6,6) * t214 + (Icges(6,4) * t197 - Icges(6,2) * t196) * t212;
t122 = -Icges(6,3) * t214 + (Icges(6,5) * t197 - Icges(6,6) * t196) * t212;
t121 = t186 * V_base(4) - t187 * V_base(5) + V_base(3);
t119 = -pkin(8) * t214 + t212 * t241;
t116 = t161 * rSges(4,1) + t160 * rSges(4,2) + rSges(4,3) * t247;
t115 = rSges(4,1) * t159 + rSges(4,2) * t158 + rSges(4,3) * t248;
t114 = Icges(4,1) * t161 + Icges(4,4) * t160 + Icges(4,5) * t247;
t113 = Icges(4,1) * t159 + Icges(4,4) * t158 + Icges(4,5) * t248;
t112 = Icges(4,4) * t161 + Icges(4,2) * t160 + Icges(4,6) * t247;
t111 = Icges(4,4) * t159 + Icges(4,2) * t158 + Icges(4,6) * t248;
t110 = Icges(4,5) * t161 + Icges(4,6) * t160 + Icges(4,3) * t247;
t109 = Icges(4,5) * t159 + Icges(4,6) * t158 + Icges(4,3) * t248;
t107 = t148 * rSges(5,1) + t147 * rSges(5,2) + rSges(5,3) * t247;
t106 = rSges(5,1) * t146 + rSges(5,2) * t145 + rSges(5,3) * t248;
t105 = Icges(5,1) * t148 + Icges(5,4) * t147 + Icges(5,5) * t247;
t104 = Icges(5,1) * t146 + Icges(5,4) * t145 + Icges(5,5) * t248;
t103 = Icges(5,4) * t148 + Icges(5,2) * t147 + Icges(5,6) * t247;
t102 = Icges(5,4) * t146 + Icges(5,2) * t145 + Icges(5,6) * t248;
t101 = Icges(5,5) * t148 + Icges(5,6) * t147 + Icges(5,3) * t247;
t100 = Icges(5,5) * t146 + Icges(5,6) * t145 + Icges(5,3) * t248;
t98 = t139 * rSges(6,1) + t138 * rSges(6,2) + rSges(6,3) * t247;
t97 = rSges(6,1) * t137 + rSges(6,2) * t136 + rSges(6,3) * t248;
t96 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t247;
t95 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t248;
t94 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t247;
t93 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t248;
t92 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t247;
t91 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t248;
t90 = t185 * t190 + (-t155 - t188) * t202 + t236;
t89 = t156 * t202 - t185 * t191 + t226;
t88 = t213 * t233 + t215 * t221;
t87 = t213 * t221 - t215 * t233;
t86 = t155 * t191 - t156 * t190 + t225;
t85 = t144 * t190 + (-t115 + t242) * t202 + t232;
t84 = t116 * t202 + (-t144 - t184) * t191 + t222;
t83 = t115 * t191 + (-t116 - t165) * t190 + t220;
t82 = -t106 * t183 + t132 * t162 + t219;
t81 = t107 * t183 - t132 * t163 + t218;
t80 = t106 * t163 - t107 * t162 + t217;
t79 = t119 * t162 + t127 * t134 - t167 * t97 - t183 * t87 + t219;
t78 = -t119 * t163 - t127 * t135 + t167 * t98 + t183 * t88 + t218;
t77 = -t134 * t98 + t135 * t97 - t162 * t88 + t163 * t87 + t217;
t1 = m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t121 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t86 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t163 * ((t101 * t247 + t147 * t103 + t148 * t105) * t163 + (t100 * t247 + t147 * t102 + t148 * t104) * t162 + (t129 * t247 + t147 * t130 + t148 * t131) * t183) / 0.2e1 + t162 * ((t101 * t248 + t103 * t145 + t105 * t146) * t163 + (t100 * t248 + t145 * t102 + t146 * t104) * t162 + (t129 * t248 + t130 * t145 + t131 * t146) * t183) / 0.2e1 + t183 * ((-t100 * t162 - t101 * t163 - t129 * t183) * t214 + ((-t103 * t199 + t105 * t200) * t163 + (-t102 * t199 + t104 * t200) * t162 + (-t130 * t199 + t131 * t200) * t183) * t212) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + t135 * ((t138 * t94 + t139 * t96 + t92 * t247) * t135 + (t138 * t93 + t139 * t95 + t247 * t91) * t134 + (t122 * t247 + t138 * t123 + t139 * t124) * t167) / 0.2e1 + t134 * ((t136 * t94 + t137 * t96 + t248 * t92) * t135 + (t136 * t93 + t137 * t95 + t91 * t248) * t134 + (t122 * t248 + t123 * t136 + t124 * t137) * t167) / 0.2e1 + t167 * ((-t122 * t167 - t91 * t134 - t92 * t135) * t214 + ((-t196 * t94 + t197 * t96) * t135 + (-t196 * t93 + t197 * t95) * t134 + (-t123 * t196 + t124 * t197) * t167) * t212) / 0.2e1 + (t216 * t213 - t224 * t215 + (t110 * t248 + t112 * t158 + t114 * t159) * t191 + (t109 * t248 + t111 * t158 + t113 * t159) * t190 + (t140 * t248 + t141 * t158 + t142 * t159) * t202) * t190 / 0.2e1 + (t224 * t213 + t216 * t215 + (t110 * t247 + t160 * t112 + t161 * t114) * t191 + (t109 * t247 + t160 * t111 + t161 * t113) * t190 + (t140 * t247 + t160 * t141 + t161 * t142) * t202) * t191 / 0.2e1 + ((-t213 * t178 + t181 * t215 + Icges(1,4)) * V_base(5) + (-t213 * t179 + t182 * t215 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t215 * t178 + t213 * t181 + Icges(1,2)) * V_base(5) + (t179 * t215 + t213 * t182 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t152 * t214 + t154 * t212) * t191 + (t151 * t214 + t153 * t212) * t190 + (-t109 * t190 - t110 * t191) * t214 + ((-t112 * t209 + t114 * t210) * t191 + (-t111 * t209 + t113 * t210) * t190) * t212 + (Icges(2,3) + (t177 - t140) * t214 + (-t141 * t209 + t142 * t210 + t180) * t212) * t202) * t202 / 0.2e1 + t202 * V_base(4) * (Icges(2,5) * t215 - Icges(2,6) * t213) + t202 * V_base(5) * (Icges(2,5) * t213 + Icges(2,6) * t215) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
