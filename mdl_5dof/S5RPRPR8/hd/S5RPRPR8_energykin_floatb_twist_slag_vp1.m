% Calculate kinetic energy for
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:14
% EndTime: 2019-12-31 18:21:16
% DurationCPUTime: 2.03s
% Computational Cost: add. (1402->293), mult. (1378->418), div. (0->0), fcn. (1270->10), ass. (0->137)
t198 = sin(qJ(1));
t240 = pkin(1) * t198;
t200 = cos(qJ(1));
t239 = pkin(1) * t200;
t195 = cos(pkin(9));
t238 = pkin(4) * t195;
t237 = -pkin(5) - qJ(2);
t236 = Icges(2,4) * t198;
t193 = qJ(1) + pkin(8);
t185 = sin(t193);
t235 = Icges(3,4) * t185;
t197 = sin(qJ(3));
t234 = Icges(4,4) * t197;
t199 = cos(qJ(3));
t233 = Icges(4,4) * t199;
t194 = sin(pkin(9));
t232 = t185 * t194;
t231 = t185 * t197;
t230 = t185 * t199;
t187 = cos(t193);
t229 = t187 * t194;
t228 = t187 * t197;
t227 = t187 * t199;
t226 = t194 * t199;
t225 = t195 * t199;
t223 = qJD(4) * t197;
t222 = qJD(5) * t197;
t188 = V_base(6) + qJD(1);
t221 = t188 * t239 + V_base(2);
t220 = V_base(5) * pkin(5) + V_base(1);
t162 = qJD(3) * t185 + V_base(4);
t156 = pkin(2) * t185 - pkin(6) * t187;
t217 = -t156 - t240;
t216 = V_base(5) * qJ(2) + t220;
t215 = V_base(4) * t240 + qJD(2) + V_base(3);
t212 = pkin(3) * t199 + qJ(4) * t197;
t144 = t212 * t185;
t214 = -t144 + t217;
t161 = -qJD(3) * t187 + V_base(5);
t213 = rSges(4,1) * t199 - rSges(4,2) * t197;
t211 = Icges(4,1) * t199 - t234;
t210 = -Icges(4,2) * t197 + t233;
t209 = Icges(4,5) * t199 - Icges(4,6) * t197;
t175 = pkin(3) * t197 - qJ(4) * t199;
t208 = t161 * t175 + t187 * t223 + t216;
t207 = (-Icges(4,3) * t187 + t185 * t209) * t161 + (Icges(4,3) * t185 + t187 * t209) * t162 + (Icges(4,5) * t197 + Icges(4,6) * t199) * t188;
t206 = pkin(7) * t197 + t199 * t238;
t157 = pkin(2) * t187 + pkin(6) * t185;
t205 = t188 * t157 + t237 * V_base(4) + t221;
t145 = t212 * t187;
t204 = t188 * t145 + t185 * t223 + t205;
t203 = V_base(4) * t156 + (-t157 - t239) * V_base(5) + t215;
t202 = -qJD(4) * t199 + t162 * t144 + t203;
t115 = -Icges(4,6) * t187 + t185 * t210;
t116 = Icges(4,6) * t185 + t187 * t210;
t117 = -Icges(4,5) * t187 + t185 * t211;
t118 = Icges(4,5) * t185 + t187 * t211;
t166 = Icges(4,2) * t199 + t234;
t169 = Icges(4,1) * t197 + t233;
t201 = (-t116 * t197 + t118 * t199) * t162 + (-t115 * t197 + t117 * t199) * t161 + (-t166 * t197 + t169 * t199) * t188;
t192 = pkin(9) + qJ(5);
t190 = Icges(2,4) * t200;
t186 = cos(t192);
t184 = sin(t192);
t182 = Icges(3,4) * t187;
t178 = rSges(2,1) * t200 - t198 * rSges(2,2);
t177 = t198 * rSges(2,1) + rSges(2,2) * t200;
t176 = rSges(4,1) * t197 + rSges(4,2) * t199;
t174 = -qJD(5) * t199 + t188;
t171 = Icges(2,1) * t200 - t236;
t170 = Icges(2,1) * t198 + t190;
t168 = -Icges(2,2) * t198 + t190;
t167 = Icges(2,2) * t200 + t236;
t160 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t159 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t158 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t155 = rSges(3,1) * t187 - rSges(3,2) * t185;
t154 = rSges(3,1) * t185 + rSges(3,2) * t187;
t153 = Icges(3,1) * t187 - t235;
t152 = Icges(3,1) * t185 + t182;
t151 = -Icges(3,2) * t185 + t182;
t150 = Icges(3,2) * t187 + t235;
t143 = -rSges(5,3) * t199 + (rSges(5,1) * t195 - rSges(5,2) * t194) * t197;
t142 = t187 * t222 + t162;
t141 = t185 * t222 + t161;
t140 = -Icges(5,5) * t199 + (Icges(5,1) * t195 - Icges(5,4) * t194) * t197;
t139 = -Icges(5,6) * t199 + (Icges(5,4) * t195 - Icges(5,2) * t194) * t197;
t138 = -Icges(5,3) * t199 + (Icges(5,5) * t195 - Icges(5,6) * t194) * t197;
t136 = t187 * t225 + t232;
t135 = t185 * t195 - t187 * t226;
t134 = t185 * t225 - t229;
t133 = -t185 * t226 - t187 * t195;
t132 = -rSges(6,3) * t199 + (rSges(6,1) * t186 - rSges(6,2) * t184) * t197;
t131 = -Icges(6,5) * t199 + (Icges(6,1) * t186 - Icges(6,4) * t184) * t197;
t130 = -Icges(6,6) * t199 + (Icges(6,4) * t186 - Icges(6,2) * t184) * t197;
t129 = -Icges(6,3) * t199 + (Icges(6,5) * t186 - Icges(6,6) * t184) * t197;
t128 = t184 * t185 + t186 * t227;
t127 = -t184 * t227 + t185 * t186;
t126 = -t184 * t187 + t186 * t230;
t125 = -t184 * t230 - t186 * t187;
t123 = -pkin(7) * t199 + t197 * t238;
t122 = rSges(4,3) * t185 + t187 * t213;
t121 = -rSges(4,3) * t187 + t185 * t213;
t120 = V_base(5) * rSges(2,3) - t177 * t188 + t220;
t119 = t178 * t188 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t112 = t177 * V_base(4) - t178 * V_base(5) + V_base(3);
t110 = V_base(5) * rSges(3,3) + (-t154 - t240) * t188 + t216;
t109 = t155 * t188 + (-rSges(3,3) + t237) * V_base(4) + t221;
t108 = V_base(4) * t154 + (-t155 - t239) * V_base(5) + t215;
t107 = pkin(4) * t232 + t187 * t206;
t106 = -pkin(4) * t229 + t185 * t206;
t105 = rSges(5,1) * t136 + rSges(5,2) * t135 + rSges(5,3) * t228;
t104 = rSges(5,1) * t134 + rSges(5,2) * t133 + rSges(5,3) * t231;
t103 = Icges(5,1) * t136 + Icges(5,4) * t135 + Icges(5,5) * t228;
t102 = Icges(5,1) * t134 + Icges(5,4) * t133 + Icges(5,5) * t231;
t101 = Icges(5,4) * t136 + Icges(5,2) * t135 + Icges(5,6) * t228;
t100 = Icges(5,4) * t134 + Icges(5,2) * t133 + Icges(5,6) * t231;
t99 = Icges(5,5) * t136 + Icges(5,6) * t135 + Icges(5,3) * t228;
t98 = Icges(5,5) * t134 + Icges(5,6) * t133 + Icges(5,3) * t231;
t97 = rSges(6,1) * t128 + rSges(6,2) * t127 + rSges(6,3) * t228;
t96 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t231;
t95 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t228;
t94 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t231;
t93 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t228;
t92 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t231;
t91 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t228;
t90 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t231;
t89 = t161 * t176 + (-t121 + t217) * t188 + t216;
t88 = t122 * t188 - t162 * t176 + t205;
t87 = t162 * t121 - t161 * t122 + t203;
t86 = t143 * t161 + (-t104 + t214) * t188 + t208;
t85 = t105 * t188 + (-t143 - t175) * t162 + t204;
t84 = t162 * t104 + (-t105 - t145) * t161 + t202;
t83 = t123 * t161 + t132 * t141 - t174 * t96 + (-t106 + t214) * t188 + t208;
t82 = t107 * t188 - t132 * t142 + t174 * t97 + (-t123 - t175) * t162 + t204;
t81 = t162 * t106 - t141 * t97 + t142 * t96 + (-t107 - t145) * t161 + t202;
t1 = m(1) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(2) * (t112 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t142 * ((t127 * t93 + t128 * t95 + t91 * t228) * t142 + (t127 * t92 + t128 * t94 + t228 * t90) * t141 + (t127 * t130 + t128 * t131 + t129 * t228) * t174) / 0.2e1 + t141 * ((t125 * t93 + t126 * t95 + t231 * t91) * t142 + (t125 * t92 + t126 * t94 + t90 * t231) * t141 + (t125 * t130 + t126 * t131 + t129 * t231) * t174) / 0.2e1 + t174 * ((-t129 * t174 - t90 * t141 - t91 * t142) * t199 + ((-t184 * t93 + t186 * t95) * t142 + (-t184 * t92 + t186 * t94) * t141 + (-t130 * t184 + t131 * t186) * t174) * t197) / 0.2e1 + (t201 * t185 - t207 * t187 + (t101 * t133 + t103 * t134 + t231 * t99) * t162 + (t133 * t100 + t134 * t102 + t231 * t98) * t161 + (t133 * t139 + t134 * t140 + t138 * t231) * t188) * t161 / 0.2e1 + (t207 * t185 + t201 * t187 + (t135 * t101 + t136 * t103 + t228 * t99) * t162 + (t100 * t135 + t102 * t136 + t228 * t98) * t161 + (t135 * t139 + t136 * t140 + t138 * t228) * t188) * t162 / 0.2e1 + ((-t150 * t185 + t152 * t187 - t198 * t167 + t170 * t200 + Icges(1,4)) * V_base(5) + (-t185 * t151 + t187 * t153 - t198 * t168 + t200 * t171 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t187 * t150 + t185 * t152 + t200 * t167 + t198 * t170 + Icges(1,2)) * V_base(5) + (t151 * t187 + t153 * t185 + t168 * t200 + t198 * t171 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t116 * t199 + t118 * t197) * t162 + (t115 * t199 + t117 * t197) * t161 + (-t98 * t161 - t99 * t162) * t199 + ((-t101 * t194 + t103 * t195) * t162 + (-t100 * t194 + t102 * t195) * t161) * t197 + (Icges(2,3) + Icges(3,3) + (t166 - t138) * t199 + (-t139 * t194 + t140 * t195 + t169) * t197) * t188) * t188 / 0.2e1 + t188 * V_base(5) * (Icges(2,5) * t198 + Icges(3,5) * t185 + Icges(2,6) * t200 + Icges(3,6) * t187) + t188 * V_base(4) * (Icges(2,5) * t200 + Icges(3,5) * t187 - Icges(2,6) * t198 - Icges(3,6) * t185) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
