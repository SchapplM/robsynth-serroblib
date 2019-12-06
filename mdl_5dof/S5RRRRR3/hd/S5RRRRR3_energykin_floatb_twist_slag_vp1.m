% Calculate kinetic energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:02
% EndTime: 2019-12-05 18:55:04
% DurationCPUTime: 2.30s
% Computational Cost: add. (1399->296), mult. (1578->459), div. (0->0), fcn. (1488->10), ass. (0->148)
t189 = sin(qJ(2));
t238 = pkin(1) * t189;
t192 = cos(qJ(2));
t237 = pkin(1) * t192;
t191 = cos(qJ(4));
t236 = pkin(3) * t191;
t190 = sin(qJ(1));
t235 = Icges(2,4) * t190;
t234 = Icges(3,4) * t189;
t233 = Icges(3,4) * t192;
t187 = qJ(2) + qJ(3);
t181 = sin(t187);
t232 = Icges(4,4) * t181;
t183 = cos(t187);
t231 = Icges(4,4) * t183;
t230 = t181 * t190;
t193 = cos(qJ(1));
t229 = t181 * t193;
t228 = t183 * t193;
t188 = sin(qJ(4));
t227 = t188 * t193;
t186 = qJ(4) + qJ(5);
t180 = sin(t186);
t226 = t190 * t180;
t182 = cos(t186);
t225 = t190 * t182;
t224 = t190 * t188;
t223 = t190 * t191;
t222 = t191 * t193;
t221 = qJD(4) * t181;
t220 = qJD(5) * t181;
t219 = V_base(5) * pkin(4) + V_base(1);
t218 = t190 * t237;
t217 = t193 * t237;
t175 = qJD(2) * t190 + V_base(4);
t214 = t183 * t236;
t177 = V_base(6) + qJD(1);
t174 = -qJD(2) * t193 + V_base(5);
t213 = t174 * t238 + t219;
t153 = qJD(3) * t190 + t175;
t212 = pkin(2) * t183 + pkin(5) * t181;
t211 = rSges(3,1) * t192 - rSges(3,2) * t189;
t210 = rSges(4,1) * t183 - rSges(4,2) * t181;
t209 = -V_base(4) * pkin(4) + V_base(2);
t125 = t193 * t221 + t153;
t208 = Icges(3,1) * t192 - t234;
t207 = Icges(4,1) * t183 - t232;
t206 = -Icges(3,2) * t189 + t233;
t205 = -Icges(4,2) * t181 + t231;
t204 = Icges(3,5) * t192 - Icges(3,6) * t189;
t203 = Icges(4,5) * t183 - Icges(4,6) * t181;
t152 = V_base(5) + (-qJD(2) - qJD(3)) * t193;
t202 = -t174 * t217 + t175 * t218 + V_base(3);
t124 = t190 * t221 + t152;
t201 = (-Icges(4,3) * t193 + t190 * t203) * t152 + (Icges(4,3) * t190 + t193 * t203) * t153 + (Icges(4,5) * t181 + Icges(4,6) * t183) * t177;
t200 = (-Icges(3,3) * t193 + t190 * t204) * t174 + (Icges(3,3) * t190 + t193 * t204) * t175 + (Icges(3,5) * t189 + Icges(3,6) * t192) * t177;
t199 = -t175 * t238 + t177 * t217 + t209;
t138 = t212 * t190;
t139 = t212 * t193;
t198 = t153 * t138 - t152 * t139 + t202;
t151 = pkin(2) * t181 - pkin(5) * t183;
t197 = t152 * t151 + (-t138 - t218) * t177 + t213;
t196 = t177 * t139 - t151 * t153 + t199;
t118 = -Icges(4,6) * t193 + t190 * t205;
t119 = Icges(4,6) * t190 + t193 * t205;
t120 = -Icges(4,5) * t193 + t190 * t207;
t121 = Icges(4,5) * t190 + t193 * t207;
t148 = Icges(4,2) * t183 + t232;
t149 = Icges(4,1) * t181 + t231;
t195 = (-t119 * t181 + t121 * t183) * t153 + (-t118 * t181 + t120 * t183) * t152 + (-t148 * t181 + t149 * t183) * t177;
t132 = -Icges(3,6) * t193 + t190 * t206;
t133 = Icges(3,6) * t190 + t193 * t206;
t134 = -Icges(3,5) * t193 + t190 * t208;
t135 = Icges(3,5) * t190 + t193 * t208;
t163 = Icges(3,2) * t192 + t234;
t166 = Icges(3,1) * t189 + t233;
t194 = (-t133 * t189 + t135 * t192) * t175 + (-t132 * t189 + t134 * t192) * t174 + (-t163 * t189 + t166 * t192) * t177;
t184 = Icges(2,4) * t193;
t171 = rSges(2,1) * t193 - t190 * rSges(2,2);
t170 = t190 * rSges(2,1) + rSges(2,2) * t193;
t169 = rSges(3,1) * t189 + rSges(3,2) * t192;
t168 = Icges(2,1) * t193 - t235;
t167 = Icges(2,1) * t190 + t184;
t165 = -Icges(2,2) * t190 + t184;
t164 = Icges(2,2) * t193 + t235;
t159 = -qJD(4) * t183 + t177;
t158 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t157 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t156 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t150 = rSges(4,1) * t181 + rSges(4,2) * t183;
t145 = t183 * t222 + t224;
t144 = -t183 * t227 + t223;
t143 = t183 * t223 - t227;
t142 = -t183 * t224 - t222;
t141 = t236 * t181;
t140 = (-qJD(4) - qJD(5)) * t183 + t177;
t137 = t190 * rSges(3,3) + t193 * t211;
t136 = -rSges(3,3) * t193 + t190 * t211;
t129 = t182 * t228 + t226;
t128 = -t180 * t228 + t225;
t127 = -t180 * t193 + t183 * t225;
t126 = -t182 * t193 - t183 * t226;
t123 = t190 * rSges(4,3) + t193 * t210;
t122 = -rSges(4,3) * t193 + t190 * t210;
t114 = -rSges(5,3) * t183 + (rSges(5,1) * t191 - rSges(5,2) * t188) * t181;
t113 = -Icges(5,5) * t183 + (Icges(5,1) * t191 - Icges(5,4) * t188) * t181;
t112 = -Icges(5,6) * t183 + (Icges(5,4) * t191 - Icges(5,2) * t188) * t181;
t111 = -Icges(5,3) * t183 + (Icges(5,5) * t191 - Icges(5,6) * t188) * t181;
t110 = V_base(5) * rSges(2,3) - t170 * t177 + t219;
t109 = t171 * t177 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t107 = t170 * V_base(4) - t171 * V_base(5) + V_base(3);
t106 = -rSges(6,3) * t183 + (rSges(6,1) * t182 - rSges(6,2) * t180) * t181;
t105 = -Icges(6,5) * t183 + (Icges(6,1) * t182 - Icges(6,4) * t180) * t181;
t104 = -Icges(6,6) * t183 + (Icges(6,4) * t182 - Icges(6,2) * t180) * t181;
t103 = -Icges(6,3) * t183 + (Icges(6,5) * t182 - Icges(6,6) * t180) * t181;
t102 = pkin(3) * t224 + t193 * t214;
t101 = -pkin(3) * t227 + t190 * t214;
t100 = t193 * t220 + t125;
t99 = t190 * t220 + t124;
t97 = t145 * rSges(5,1) + t144 * rSges(5,2) + rSges(5,3) * t229;
t96 = rSges(5,1) * t143 + rSges(5,2) * t142 + rSges(5,3) * t230;
t95 = Icges(5,1) * t145 + Icges(5,4) * t144 + Icges(5,5) * t229;
t94 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t230;
t93 = Icges(5,4) * t145 + Icges(5,2) * t144 + Icges(5,6) * t229;
t92 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t230;
t91 = Icges(5,5) * t145 + Icges(5,6) * t144 + Icges(5,3) * t229;
t90 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t230;
t89 = -t136 * t177 + t169 * t174 + t219;
t88 = t137 * t177 - t169 * t175 + t209;
t87 = t129 * rSges(6,1) + t128 * rSges(6,2) + rSges(6,3) * t229;
t86 = rSges(6,1) * t127 + rSges(6,2) * t126 + rSges(6,3) * t230;
t85 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t229;
t84 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t230;
t83 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t229;
t82 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t230;
t81 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t229;
t80 = Icges(6,5) * t127 + Icges(6,6) * t126 + Icges(6,3) * t230;
t79 = t136 * t175 - t137 * t174 + V_base(3);
t78 = t150 * t152 + (-t122 - t218) * t177 + t213;
t77 = t123 * t177 - t150 * t153 + t199;
t76 = t153 * t122 - t152 * t123 + t202;
t75 = t114 * t124 - t159 * t96 + t197;
t74 = -t114 * t125 + t159 * t97 + t196;
t73 = -t124 * t97 + t125 * t96 + t198;
t72 = -t101 * t159 + t106 * t99 + t124 * t141 - t140 * t86 + t197;
t71 = -t100 * t106 + t102 * t159 - t125 * t141 + t140 * t87 + t196;
t70 = t100 * t86 + t125 * t101 - t124 * t102 - t99 * t87 + t198;
t1 = m(1) * (t156 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(2) * (t107 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(3) * (t79 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t175 * (t200 * t190 + t194 * t193) / 0.2e1 + t174 * (t194 * t190 - t200 * t193) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + t153 * (t201 * t190 + t195 * t193) / 0.2e1 + t152 * (t195 * t190 - t201 * t193) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + t125 * ((t144 * t93 + t145 * t95 + t91 * t229) * t125 + (t144 * t92 + t145 * t94 + t229 * t90) * t124 + (t111 * t229 + t144 * t112 + t145 * t113) * t159) / 0.2e1 + t124 * ((t142 * t93 + t143 * t95 + t230 * t91) * t125 + (t142 * t92 + t143 * t94 + t90 * t230) * t124 + (t111 * t230 + t112 * t142 + t113 * t143) * t159) / 0.2e1 + t159 * ((-t111 * t159 - t124 * t90 - t125 * t91) * t183 + ((-t188 * t93 + t191 * t95) * t125 + (-t188 * t92 + t191 * t94) * t124 + (-t112 * t188 + t113 * t191) * t159) * t181) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + t100 * ((t128 * t83 + t129 * t85 + t81 * t229) * t100 + (t128 * t82 + t129 * t84 + t229 * t80) * t99 + (t103 * t229 + t128 * t104 + t129 * t105) * t140) / 0.2e1 + t99 * ((t126 * t83 + t127 * t85 + t230 * t81) * t100 + (t126 * t82 + t127 * t84 + t80 * t230) * t99 + (t103 * t230 + t104 * t126 + t105 * t127) * t140) / 0.2e1 + t140 * ((-t100 * t81 - t103 * t140 - t80 * t99) * t183 + ((-t180 * t83 + t182 * t85) * t100 + (-t180 * t82 + t182 * t84) * t99 + (-t104 * t180 + t105 * t182) * t140) * t181) / 0.2e1 + ((-t190 * t164 + t167 * t193 + Icges(1,4)) * V_base(5) + (-t190 * t165 + t168 * t193 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t164 * t193 + t190 * t167 + Icges(1,2)) * V_base(5) + (t165 * t193 + t190 * t168 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t133 * t192 + t135 * t189) * t175 + (t132 * t192 + t134 * t189) * t174 + (t119 * t183 + t121 * t181) * t153 + (t118 * t183 + t120 * t181) * t152 + (t148 * t183 + t149 * t181 + t163 * t192 + t166 * t189 + Icges(2,3)) * t177) * t177 / 0.2e1 + t177 * V_base(4) * (Icges(2,5) * t193 - Icges(2,6) * t190) + V_base(5) * t177 * (Icges(2,5) * t190 + Icges(2,6) * t193) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
