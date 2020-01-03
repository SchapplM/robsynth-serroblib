% Calculate kinetic energy for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:08
% EndTime: 2019-12-31 19:03:10
% DurationCPUTime: 2.21s
% Computational Cost: add. (1462->291), mult. (1438->429), div. (0->0), fcn. (1330->10), ass. (0->139)
t200 = sin(qJ(1));
t243 = pkin(1) * t200;
t203 = cos(qJ(1));
t242 = pkin(1) * t203;
t201 = cos(qJ(4));
t241 = pkin(4) * t201;
t240 = -pkin(5) - qJ(2);
t238 = Icges(2,4) * t200;
t196 = qJ(1) + pkin(9);
t188 = sin(t196);
t237 = Icges(3,4) * t188;
t199 = sin(qJ(3));
t236 = Icges(4,4) * t199;
t202 = cos(qJ(3));
t235 = Icges(4,4) * t202;
t198 = sin(qJ(4));
t234 = t188 * t198;
t233 = t188 * t199;
t189 = cos(t196);
t232 = t189 * t198;
t231 = t189 * t199;
t197 = qJ(4) + qJ(5);
t191 = sin(t197);
t230 = t191 * t202;
t192 = cos(t197);
t229 = t192 * t202;
t228 = t198 * t202;
t227 = t201 * t202;
t226 = qJD(4) * t199;
t225 = qJD(5) * t199;
t190 = V_base(6) + qJD(1);
t224 = t190 * t242 + V_base(2);
t223 = V_base(5) * pkin(5) + V_base(1);
t166 = qJD(3) * t188 + V_base(4);
t160 = pkin(2) * t188 - pkin(6) * t189;
t220 = -t160 - t243;
t219 = V_base(5) * qJ(2) + t223;
t218 = V_base(4) * t243 + qJD(2) + V_base(3);
t138 = t189 * t226 + t166;
t217 = pkin(3) * t202 + pkin(7) * t199;
t165 = -qJD(3) * t189 + V_base(5);
t216 = rSges(4,1) * t202 - rSges(4,2) * t199;
t215 = Icges(4,1) * t202 - t236;
t214 = -Icges(4,2) * t199 + t235;
t213 = Icges(4,5) * t202 - Icges(4,6) * t199;
t137 = t188 * t226 + t165;
t212 = pkin(8) * t199 + t202 * t241;
t211 = (-Icges(4,3) * t189 + t188 * t213) * t165 + (Icges(4,3) * t188 + t189 * t213) * t166 + (Icges(4,5) * t199 + Icges(4,6) * t202) * t190;
t161 = pkin(2) * t189 + pkin(6) * t188;
t210 = t190 * t161 + t240 * V_base(4) + t224;
t148 = t217 * t188;
t182 = t199 * pkin(3) - t202 * pkin(7);
t209 = t165 * t182 + (-t148 + t220) * t190 + t219;
t208 = V_base(4) * t160 + (-t161 - t242) * V_base(5) + t218;
t149 = t217 * t189;
t207 = t190 * t149 - t166 * t182 + t210;
t206 = t166 * t148 - t149 * t165 + t208;
t118 = -Icges(4,6) * t189 + t188 * t214;
t119 = Icges(4,6) * t188 + t189 * t214;
t120 = -Icges(4,5) * t189 + t188 * t215;
t121 = Icges(4,5) * t188 + t189 * t215;
t170 = Icges(4,2) * t202 + t236;
t173 = Icges(4,1) * t199 + t235;
t205 = (-t119 * t199 + t121 * t202) * t166 + (-t118 * t199 + t120 * t202) * t165 + (-t170 * t199 + t173 * t202) * t190;
t194 = Icges(2,4) * t203;
t186 = Icges(3,4) * t189;
t181 = rSges(2,1) * t203 - rSges(2,2) * t200;
t180 = rSges(2,1) * t200 + rSges(2,2) * t203;
t179 = rSges(4,1) * t199 + rSges(4,2) * t202;
t178 = -qJD(4) * t202 + t190;
t175 = Icges(2,1) * t203 - t238;
t174 = Icges(2,1) * t200 + t194;
t172 = -Icges(2,2) * t200 + t194;
t171 = Icges(2,2) * t203 + t238;
t164 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t163 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t162 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t159 = rSges(3,1) * t189 - rSges(3,2) * t188;
t158 = rSges(3,1) * t188 + rSges(3,2) * t189;
t157 = (-qJD(4) - qJD(5)) * t202 + t190;
t156 = Icges(3,1) * t189 - t237;
t155 = Icges(3,1) * t188 + t186;
t154 = -Icges(3,2) * t188 + t186;
t153 = Icges(3,2) * t189 + t237;
t146 = -rSges(5,3) * t202 + (rSges(5,1) * t201 - rSges(5,2) * t198) * t199;
t145 = -Icges(5,5) * t202 + (Icges(5,1) * t201 - Icges(5,4) * t198) * t199;
t144 = -Icges(5,6) * t202 + (Icges(5,4) * t201 - Icges(5,2) * t198) * t199;
t143 = -Icges(5,3) * t202 + (Icges(5,5) * t201 - Icges(5,6) * t198) * t199;
t142 = t189 * t227 + t234;
t141 = t188 * t201 - t189 * t228;
t140 = t188 * t227 - t232;
t139 = -t188 * t228 - t189 * t201;
t135 = -rSges(6,3) * t202 + (rSges(6,1) * t192 - rSges(6,2) * t191) * t199;
t134 = -Icges(6,5) * t202 + (Icges(6,1) * t192 - Icges(6,4) * t191) * t199;
t133 = -Icges(6,6) * t202 + (Icges(6,4) * t192 - Icges(6,2) * t191) * t199;
t132 = -Icges(6,3) * t202 + (Icges(6,5) * t192 - Icges(6,6) * t191) * t199;
t131 = t188 * t191 + t189 * t229;
t130 = t188 * t192 - t189 * t230;
t129 = t188 * t229 - t189 * t191;
t128 = -t188 * t230 - t189 * t192;
t127 = -pkin(8) * t202 + t199 * t241;
t125 = rSges(4,3) * t188 + t189 * t216;
t124 = -rSges(4,3) * t189 + t188 * t216;
t123 = V_base(5) * rSges(2,3) - t180 * t190 + t223;
t122 = t181 * t190 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t115 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t114 = t189 * t225 + t138;
t113 = t188 * t225 + t137;
t111 = V_base(5) * rSges(3,3) + (-t158 - t243) * t190 + t219;
t110 = t159 * t190 + (-rSges(3,3) + t240) * V_base(4) + t224;
t109 = t158 * V_base(4) + (-t159 - t242) * V_base(5) + t218;
t108 = pkin(4) * t234 + t189 * t212;
t107 = -pkin(4) * t232 + t188 * t212;
t106 = rSges(5,1) * t142 + rSges(5,2) * t141 + rSges(5,3) * t231;
t105 = rSges(5,1) * t140 + rSges(5,2) * t139 + rSges(5,3) * t233;
t104 = Icges(5,1) * t142 + Icges(5,4) * t141 + Icges(5,5) * t231;
t103 = Icges(5,1) * t140 + Icges(5,4) * t139 + Icges(5,5) * t233;
t102 = Icges(5,4) * t142 + Icges(5,2) * t141 + Icges(5,6) * t231;
t101 = Icges(5,4) * t140 + Icges(5,2) * t139 + Icges(5,6) * t233;
t100 = Icges(5,5) * t142 + Icges(5,6) * t141 + Icges(5,3) * t231;
t99 = Icges(5,5) * t140 + Icges(5,6) * t139 + Icges(5,3) * t233;
t98 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t231;
t97 = rSges(6,1) * t129 + rSges(6,2) * t128 + rSges(6,3) * t233;
t96 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t231;
t95 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t233;
t94 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t231;
t93 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t233;
t92 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t231;
t91 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t233;
t90 = t165 * t179 + (-t124 + t220) * t190 + t219;
t89 = t125 * t190 - t166 * t179 + t210;
t88 = t124 * t166 - t125 * t165 + t208;
t87 = -t105 * t178 + t137 * t146 + t209;
t86 = t106 * t178 - t138 * t146 + t207;
t85 = t105 * t138 - t106 * t137 + t206;
t84 = -t107 * t178 + t113 * t135 + t127 * t137 - t157 * t97 + t209;
t83 = t108 * t178 - t114 * t135 - t127 * t138 + t157 * t98 + t207;
t82 = t107 * t138 - t108 * t137 - t113 * t98 + t114 * t97 + t206;
t1 = m(1) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(2) * (t115 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(3) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + t166 * (t211 * t188 + t205 * t189) / 0.2e1 + t165 * (t205 * t188 - t211 * t189) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t138 * ((t100 * t231 + t141 * t102 + t142 * t104) * t138 + (t101 * t141 + t103 * t142 + t231 * t99) * t137 + (t141 * t144 + t142 * t145 + t143 * t231) * t178) / 0.2e1 + t137 * ((t100 * t233 + t102 * t139 + t104 * t140) * t138 + (t139 * t101 + t140 * t103 + t99 * t233) * t137 + (t139 * t144 + t140 * t145 + t143 * t233) * t178) / 0.2e1 + t178 * ((-t100 * t138 - t137 * t99 - t143 * t178) * t202 + ((-t102 * t198 + t104 * t201) * t138 + (-t101 * t198 + t103 * t201) * t137 + (-t144 * t198 + t145 * t201) * t178) * t199) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + t114 * ((t130 * t94 + t131 * t96 + t92 * t231) * t114 + (t130 * t93 + t131 * t95 + t231 * t91) * t113 + (t130 * t133 + t131 * t134 + t132 * t231) * t157) / 0.2e1 + t113 * ((t128 * t94 + t129 * t96 + t233 * t92) * t114 + (t128 * t93 + t129 * t95 + t91 * t233) * t113 + (t128 * t133 + t129 * t134 + t132 * t233) * t157) / 0.2e1 + t157 * ((-t113 * t91 - t114 * t92 - t132 * t157) * t202 + ((-t191 * t94 + t192 * t96) * t114 + (-t191 * t93 + t192 * t95) * t113 + (-t133 * t191 + t134 * t192) * t157) * t199) / 0.2e1 + ((t119 * t202 + t121 * t199) * t166 + (t118 * t202 + t120 * t199) * t165 + (t202 * t170 + t199 * t173 + Icges(2,3) + Icges(3,3)) * t190) * t190 / 0.2e1 + ((-t153 * t188 + t155 * t189 - t171 * t200 + t174 * t203 + Icges(1,4)) * V_base(5) + (-t188 * t154 + t189 * t156 - t200 * t172 + t203 * t175 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t189 * t153 + t188 * t155 + t203 * t171 + t200 * t174 + Icges(1,2)) * V_base(5) + (t154 * t189 + t156 * t188 + t172 * t203 + t175 * t200 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t190 * (Icges(2,5) * t200 + Icges(3,5) * t188 + Icges(2,6) * t203 + Icges(3,6) * t189) + V_base(4) * t190 * (Icges(2,5) * t203 + Icges(3,5) * t189 - Icges(2,6) * t200 - Icges(3,6) * t188) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
