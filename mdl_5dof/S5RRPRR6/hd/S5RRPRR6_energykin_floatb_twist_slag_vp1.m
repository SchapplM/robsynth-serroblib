% Calculate kinetic energy for
% S5RRPRR6
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:26
% EndTime: 2020-01-03 12:05:28
% DurationCPUTime: 2.53s
% Computational Cost: add. (1410->296), mult. (1396->426), div. (0->0), fcn. (1288->10), ass. (0->141)
t182 = qJ(1) + qJ(2);
t175 = sin(t182);
t177 = cos(t182);
t183 = sin(pkin(9));
t184 = cos(pkin(9));
t222 = Icges(4,4) * t184;
t200 = -Icges(4,2) * t183 + t222;
t105 = -Icges(4,6) * t177 + t175 * t200;
t223 = Icges(4,4) * t183;
t201 = Icges(4,1) * t184 - t223;
t107 = -Icges(4,5) * t177 + t175 * t201;
t224 = Icges(3,4) * t177;
t234 = Icges(3,1) * t175 - t105 * t183 + t107 * t184 + t224;
t106 = -Icges(4,6) * t175 - t177 * t200;
t108 = -Icges(4,5) * t175 - t177 * t201;
t170 = Icges(3,4) * t175;
t233 = -Icges(3,1) * t177 - t106 * t183 + t108 * t184 + t170;
t187 = cos(qJ(4));
t227 = pkin(4) * t187;
t232 = pkin(8) * t183 + t184 * t227;
t231 = -pkin(5) - pkin(6);
t186 = sin(qJ(1));
t229 = pkin(1) * t186;
t188 = cos(qJ(1));
t228 = pkin(1) * t188;
t225 = Icges(2,4) * t188;
t221 = t175 * t183;
t220 = t175 * t184;
t185 = sin(qJ(4));
t219 = t175 * t185;
t218 = t177 * t183;
t217 = t177 * t184;
t216 = t177 * t185;
t215 = t184 * t185;
t214 = t184 * t187;
t213 = qJD(4) * t183;
t212 = -qJD(4) - qJD(5);
t145 = -pkin(2) * t177 - qJ(3) * t175;
t211 = V_base(5) * t145 + V_base(1);
t173 = V_base(4) + qJD(1);
t210 = t173 * t229 + V_base(3);
t209 = V_base(6) * pkin(5) + V_base(2);
t206 = t188 * V_base(5);
t149 = t175 * t213 + V_base(5);
t143 = pkin(2) * t175 - qJ(3) * t177;
t205 = -t143 - t229;
t171 = qJD(2) + t173;
t204 = V_base(6) * pkin(6) + t173 * t228 + t209;
t203 = pkin(3) * t184 + pkin(7) * t183;
t202 = rSges(4,1) * t184 - rSges(4,2) * t183;
t199 = Icges(4,5) * t184 - Icges(4,6) * t183;
t155 = Icges(4,2) * t184 + t223;
t156 = Icges(4,1) * t183 + t222;
t196 = t155 * t183 - t156 * t184;
t195 = -qJD(3) * t175 + t171 * t143 + t210;
t194 = -qJD(3) * t177 + t204;
t193 = -(-Icges(4,3) * t177 + t175 * t199) * V_base(5) - (-Icges(4,3) * t175 - t177 * t199) * V_base(6) - (Icges(4,5) * t183 + Icges(4,6) * t184) * t171;
t134 = t203 * t177;
t158 = t183 * pkin(3) - t184 * pkin(7);
t192 = V_base(6) * t158 + (t134 - t145) * t171 + t194;
t133 = t203 * t175;
t191 = t171 * t133 + (-t158 + t231) * V_base(5) + t195;
t190 = -V_base(5) * t134 + (-t133 + t205) * V_base(6) - pkin(1) * t206 + t211;
t181 = qJ(4) + qJ(5);
t178 = Icges(2,4) * t186;
t176 = cos(t181);
t174 = sin(t181);
t167 = -rSges(2,1) * t188 + rSges(2,2) * t186;
t166 = rSges(2,1) * t186 + rSges(2,2) * t188;
t164 = -Icges(2,1) * t188 + t178;
t163 = Icges(2,1) * t186 + t225;
t162 = Icges(2,2) * t186 - t225;
t161 = Icges(2,2) * t188 + t178;
t157 = rSges(4,1) * t183 + rSges(4,2) * t184;
t153 = -qJD(4) * t184 + t171;
t152 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t151 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t150 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t148 = -t177 * t213 + V_base(6);
t146 = -rSges(3,1) * t177 + rSges(3,2) * t175;
t144 = rSges(3,1) * t175 + rSges(3,2) * t177;
t140 = Icges(3,2) * t175 - t224;
t139 = Icges(3,2) * t177 + t170;
t138 = -Icges(3,5) * t177 + Icges(3,6) * t175;
t137 = Icges(3,5) * t175 + Icges(3,6) * t177;
t136 = t184 * t212 + t171;
t132 = -t177 * t214 - t219;
t131 = -t175 * t187 + t177 * t215;
t130 = t175 * t214 - t216;
t129 = -t175 * t215 - t177 * t187;
t127 = -rSges(5,3) * t184 + (rSges(5,1) * t187 - rSges(5,2) * t185) * t183;
t126 = -Icges(5,5) * t184 + (Icges(5,1) * t187 - Icges(5,4) * t185) * t183;
t125 = -Icges(5,6) * t184 + (Icges(5,4) * t187 - Icges(5,2) * t185) * t183;
t124 = -Icges(5,3) * t184 + (Icges(5,5) * t187 - Icges(5,6) * t185) * t183;
t122 = qJD(5) * t221 + t149;
t121 = t212 * t218 + V_base(6);
t120 = -t174 * t175 - t176 * t217;
t119 = t174 * t217 - t175 * t176;
t118 = -t174 * t177 + t176 * t220;
t117 = -t174 * t220 - t176 * t177;
t116 = -rSges(6,3) * t184 + (rSges(6,1) * t176 - rSges(6,2) * t174) * t183;
t115 = -Icges(6,5) * t184 + (Icges(6,1) * t176 - Icges(6,4) * t174) * t183;
t114 = -Icges(6,6) * t184 + (Icges(6,4) * t176 - Icges(6,2) * t174) * t183;
t113 = -Icges(6,3) * t184 + (Icges(6,5) * t176 - Icges(6,6) * t174) * t183;
t112 = -pkin(8) * t184 + t183 * t227;
t111 = -rSges(4,3) * t175 - t177 * t202;
t110 = -rSges(4,3) * t177 + t175 * t202;
t102 = V_base(6) * rSges(2,3) - t167 * t173 + t209;
t101 = t166 * t173 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t100 = -t166 * V_base(6) + t167 * V_base(5) + V_base(1);
t99 = V_base(6) * rSges(3,3) - t146 * t171 + t204;
t98 = t144 * t171 + (-rSges(3,3) + t231) * V_base(5) + t210;
t97 = -t144 * V_base(6) + t146 * V_base(5) + V_base(1) + (-t186 * V_base(6) - t206) * pkin(1);
t96 = rSges(5,1) * t132 + rSges(5,2) * t131 - rSges(5,3) * t218;
t95 = rSges(5,1) * t130 + rSges(5,2) * t129 + rSges(5,3) * t221;
t94 = -pkin(4) * t219 - t177 * t232;
t93 = -pkin(4) * t216 + t175 * t232;
t92 = Icges(5,1) * t132 + Icges(5,4) * t131 - Icges(5,5) * t218;
t91 = Icges(5,1) * t130 + Icges(5,4) * t129 + Icges(5,5) * t221;
t90 = Icges(5,4) * t132 + Icges(5,2) * t131 - Icges(5,6) * t218;
t89 = Icges(5,4) * t130 + Icges(5,2) * t129 + Icges(5,6) * t221;
t88 = Icges(5,5) * t132 + Icges(5,6) * t131 - Icges(5,3) * t218;
t87 = Icges(5,5) * t130 + Icges(5,6) * t129 + Icges(5,3) * t221;
t86 = rSges(6,1) * t120 + rSges(6,2) * t119 - rSges(6,3) * t218;
t85 = rSges(6,1) * t118 + rSges(6,2) * t117 + rSges(6,3) * t221;
t84 = Icges(6,1) * t120 + Icges(6,4) * t119 - Icges(6,5) * t218;
t83 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t221;
t82 = Icges(6,4) * t120 + Icges(6,2) * t119 - Icges(6,6) * t218;
t81 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t221;
t80 = Icges(6,5) * t120 + Icges(6,6) * t119 - Icges(6,3) * t218;
t79 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t221;
t78 = t157 * V_base(6) + (-t111 - t145) * t171 + t194;
t77 = t110 * t171 + (-t157 + t231) * V_base(5) + t195;
t76 = (t111 - t228) * V_base(5) + (-t110 + t205) * V_base(6) + t211;
t75 = t127 * t148 - t153 * t96 + t192;
t74 = -t127 * t149 + t153 * t95 + t191;
t73 = -t148 * t95 + t149 * t96 + t190;
t72 = t112 * t148 + t116 * t121 - t136 * t86 - t153 * t94 + t192;
t71 = -t112 * t149 - t116 * t122 + t136 * t85 + t153 * t93 + t191;
t70 = -t121 * t85 + t122 * t86 - t148 * t93 + t149 * t94 + t190;
t1 = m(1) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(3) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + t153 * ((-t124 * t153 - t88 * t148 - t87 * t149) * t184 + ((-t125 * t185 + t126 * t187) * t153 + (-t185 * t89 + t187 * t91) * t149 + (-t185 * t90 + t187 * t92) * t148) * t183) / 0.2e1 + t149 * ((t124 * t221 + t125 * t129 + t126 * t130) * t153 + (t129 * t89 + t130 * t91 + t87 * t221) * t149 + (t129 * t90 + t130 * t92 + t221 * t88) * t148) / 0.2e1 + t148 * ((-t124 * t218 + t125 * t131 + t126 * t132) * t153 + (t131 * t89 + t132 * t91 - t218 * t87) * t149 + (t131 * t90 + t132 * t92 - t88 * t218) * t148) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + t136 * ((-t113 * t136 - t80 * t121 - t79 * t122) * t184 + ((-t114 * t174 + t115 * t176) * t136 + (-t174 * t81 + t176 * t83) * t122 + (-t174 * t82 + t176 * t84) * t121) * t183) / 0.2e1 + t122 * ((t113 * t221 + t114 * t117 + t115 * t118) * t136 + (t117 * t81 + t118 * t83 + t79 * t221) * t122 + (t117 * t82 + t118 * t84 + t221 * t80) * t121) / 0.2e1 + t121 * ((-t113 * t218 + t114 * t119 + t115 * t120) * t136 + (t119 * t81 + t120 * t83 - t218 * t79) * t122 + (t119 * t82 + t120 * t84 - t80 * t218) * t121) / 0.2e1 + ((t106 * t184 + t108 * t183 + t138) * V_base(6) + (t105 * t184 + t107 * t183 + t137) * V_base(5) + (t155 * t184 + t156 * t183 + Icges(3,3)) * t171) * t171 / 0.2e1 + (t193 * t177 + (-t196 * t175 + t137) * t171 + (t140 * t177 + t162 * t188 + t164 * t186 + t175 * t233 + Icges(1,6)) * V_base(6) + (t139 * t177 + t161 * t188 + t163 * t186 + t234 * t175 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t193 * t175 + (t196 * t177 + t138) * t171 + (t140 * t175 + t162 * t186 - t164 * t188 - t233 * t177 + Icges(1,3)) * V_base(6) + (t139 * t175 + t161 * t186 - t163 * t188 - t177 * t234 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t186 + Icges(2,6) * t188) * V_base(5) + (-Icges(2,5) * t188 + Icges(2,6) * t186) * V_base(6) + Icges(2,3) * t173 / 0.2e1) * t173;
T = t1;
