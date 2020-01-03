% Calculate kinetic energy for
% S5RRPRR5
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:27
% EndTime: 2020-01-03 12:03:29
% DurationCPUTime: 2.77s
% Computational Cost: add. (1258->250), mult. (928->343), div. (0->0), fcn. (708->10), ass. (0->133)
t168 = qJ(1) + qJ(2);
t160 = sin(t168);
t161 = cos(t168);
t169 = sin(pkin(9));
t170 = cos(pkin(9));
t219 = Icges(4,4) * t170;
t192 = -Icges(4,2) * t169 + t219;
t102 = -Icges(4,6) * t161 + t160 * t192;
t220 = Icges(4,4) * t169;
t195 = Icges(4,1) * t170 - t220;
t104 = -Icges(4,5) * t161 + t160 * t195;
t221 = Icges(3,4) * t161;
t234 = Icges(3,1) * t160 - t102 * t169 + t104 * t170 + t221;
t103 = -Icges(4,6) * t160 - t161 * t192;
t105 = -Icges(4,5) * t160 - t161 * t195;
t151 = Icges(3,4) * t160;
t233 = -Icges(3,1) * t161 - t103 * t169 + t105 * t170 + t151;
t167 = pkin(9) + qJ(4);
t158 = qJ(5) + t167;
t153 = cos(t158);
t152 = sin(t158);
t216 = Icges(6,4) * t152;
t110 = Icges(6,2) * t153 + t216;
t215 = Icges(6,4) * t153;
t111 = Icges(6,1) * t152 + t215;
t211 = -qJD(4) - qJD(5);
t112 = t160 * t211 + V_base(6);
t113 = t161 * t211 + V_base(5);
t159 = V_base(4) + qJD(1);
t154 = qJD(2) + t159;
t190 = -Icges(6,2) * t152 + t215;
t83 = -Icges(6,6) * t161 + t160 * t190;
t84 = -Icges(6,6) * t160 - t161 * t190;
t193 = Icges(6,1) * t153 - t216;
t85 = -Icges(6,5) * t161 + t160 * t193;
t86 = -Icges(6,5) * t160 - t161 * t193;
t232 = (t110 * t152 - t111 * t153) * t154 + (t152 * t83 - t153 * t85) * t113 + (t152 * t84 - t153 * t86) * t112;
t157 = cos(t167);
t156 = sin(t167);
t218 = Icges(5,4) * t156;
t117 = Icges(5,2) * t157 + t218;
t217 = Icges(5,4) * t157;
t118 = Icges(5,1) * t156 + t217;
t138 = -qJD(4) * t160 + V_base(6);
t139 = -qJD(4) * t161 + V_base(5);
t191 = -Icges(5,2) * t156 + t217;
t91 = -Icges(5,6) * t161 + t160 * t191;
t92 = -Icges(5,6) * t160 - t161 * t191;
t194 = Icges(5,1) * t157 - t218;
t93 = -Icges(5,5) * t161 + t160 * t194;
t94 = -Icges(5,5) * t160 - t161 * t194;
t231 = (t117 * t156 - t118 * t157) * t154 + (t156 * t91 - t157 * t93) * t139 + (t156 * t92 - t157 * t94) * t138;
t230 = -pkin(5) - pkin(6);
t172 = sin(qJ(1));
t228 = pkin(1) * t172;
t173 = cos(qJ(1));
t227 = pkin(1) * t173;
t226 = pkin(3) * t169;
t225 = pkin(4) * t156;
t224 = t170 * pkin(3);
t128 = -pkin(2) * t161 - qJ(3) * t160;
t80 = -pkin(7) * t160 - t161 * t224;
t223 = -t128 - t80;
t222 = Icges(2,4) * t173;
t213 = pkin(4) * t157;
t210 = V_base(5) * t128 + V_base(1);
t209 = t159 * t228 + V_base(3);
t208 = V_base(6) * pkin(5) + V_base(2);
t205 = t173 * V_base(5);
t126 = pkin(2) * t160 - qJ(3) * t161;
t204 = -t126 - t228;
t203 = V_base(6) * pkin(6) + t159 * t227 + t208;
t202 = rSges(4,1) * t170 - rSges(4,2) * t169;
t201 = rSges(5,1) * t157 - rSges(5,2) * t156;
t200 = rSges(6,1) * t153 - rSges(6,2) * t152;
t189 = Icges(4,5) * t170 - Icges(4,6) * t169;
t188 = Icges(5,5) * t157 - Icges(5,6) * t156;
t187 = Icges(6,5) * t153 - Icges(6,6) * t152;
t135 = Icges(4,2) * t170 + t220;
t136 = Icges(4,1) * t169 + t219;
t182 = t135 * t169 - t136 * t170;
t181 = -qJD(3) * t160 + t154 * t126 + t209;
t180 = -qJD(3) * t161 + t203;
t179 = -(Icges(6,5) * t152 + Icges(6,6) * t153) * t154 - (-Icges(6,3) * t160 - t161 * t187) * t112 - (-Icges(6,3) * t161 + t160 * t187) * t113;
t178 = -(Icges(5,5) * t156 + Icges(5,6) * t157) * t154 - (-Icges(5,3) * t160 - t161 * t188) * t138 - (-Icges(5,3) * t161 + t160 * t188) * t139;
t177 = V_base(6) * t226 + t180;
t176 = -(-Icges(4,3) * t161 + t160 * t189) * V_base(5) - (-Icges(4,3) * t160 - t161 * t189) * V_base(6) - (Icges(4,5) * t169 + Icges(4,6) * t170) * t154;
t79 = -pkin(7) * t161 + t160 * t224;
t175 = t154 * t79 + (-t226 + t230) * V_base(5) + t181;
t174 = V_base(5) * t80 + (t204 - t79) * V_base(6) - pkin(1) * t205 + t210;
t163 = Icges(2,4) * t172;
t147 = -rSges(2,1) * t173 + t172 * rSges(2,2);
t146 = t172 * rSges(2,1) + rSges(2,2) * t173;
t145 = -Icges(2,1) * t173 + t163;
t144 = Icges(2,1) * t172 + t222;
t143 = Icges(2,2) * t172 - t222;
t142 = Icges(2,2) * t173 + t163;
t137 = rSges(4,1) * t169 + rSges(4,2) * t170;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = -rSges(3,1) * t161 + rSges(3,2) * t160;
t127 = rSges(3,1) * t160 + rSges(3,2) * t161;
t123 = Icges(3,2) * t160 - t221;
t122 = Icges(3,2) * t161 + t151;
t121 = -Icges(3,5) * t161 + Icges(3,6) * t160;
t120 = Icges(3,5) * t160 + Icges(3,6) * t161;
t119 = rSges(5,1) * t156 + rSges(5,2) * t157;
t115 = rSges(6,1) * t152 + rSges(6,2) * t153;
t107 = -rSges(4,3) * t160 - t161 * t202;
t106 = -rSges(4,3) * t161 + t160 * t202;
t99 = V_base(6) * rSges(2,3) - t147 * t159 + t208;
t98 = t146 * t159 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t146 * V_base(6) + t147 * V_base(5) + V_base(1);
t96 = -rSges(5,3) * t160 - t161 * t201;
t95 = -rSges(5,3) * t161 + t160 * t201;
t88 = -rSges(6,3) * t160 - t161 * t200;
t87 = -rSges(6,3) * t161 + t160 * t200;
t76 = V_base(6) * rSges(3,3) - t129 * t154 + t203;
t75 = t127 * t154 + (-rSges(3,3) + t230) * V_base(5) + t209;
t74 = -pkin(8) * t160 - t161 * t213;
t73 = -pkin(8) * t161 + t160 * t213;
t72 = -V_base(6) * t127 + V_base(5) * t129 + V_base(1) + (-t172 * V_base(6) - t205) * pkin(1);
t71 = t137 * V_base(6) + (-t107 - t128) * t154 + t180;
t70 = t106 * t154 + (-t137 + t230) * V_base(5) + t181;
t69 = (t107 - t227) * V_base(5) + (-t106 + t204) * V_base(6) + t210;
t68 = t119 * t138 + (-t96 + t223) * t154 + t177;
t67 = -t119 * t139 + t154 * t95 + t175;
t66 = -t138 * t95 + t139 * t96 + t174;
t65 = t138 * t225 + t112 * t115 + (-t74 - t88 + t223) * t154 + t177;
t64 = -t139 * t225 - t113 * t115 + (t73 + t87) * t154 + t175;
t63 = -t112 * t87 + t113 * t88 - t138 * t73 + t139 * t74 + t174;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t139 * (-t231 * t160 + t178 * t161) / 0.2e1 + t138 * (t178 * t160 + t231 * t161) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t113 * (-t232 * t160 + t179 * t161) / 0.2e1 + t112 * (t179 * t160 + t232 * t161) / 0.2e1 + ((t156 * t93 + t157 * t91) * t139 + (t156 * t94 + t157 * t92) * t138 + (t152 * t85 + t153 * t83) * t113 + (t152 * t86 + t153 * t84) * t112 + (t103 * t170 + t105 * t169 + t121) * V_base(6) + (t102 * t170 + t104 * t169 + t120) * V_base(5) + (t153 * t110 + t152 * t111 + t157 * t117 + t156 * t118 + t170 * t135 + t169 * t136 + Icges(3,3)) * t154) * t154 / 0.2e1 + (t176 * t161 + (-t182 * t160 + t120) * t154 + (t123 * t161 + t143 * t173 + t172 * t145 + t233 * t160 + Icges(1,6)) * V_base(6) + (t161 * t122 + t173 * t142 + t172 * t144 + t234 * t160 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t176 * t160 + (t182 * t161 + t121) * t154 + (t160 * t123 + t172 * t143 - t173 * t145 - t233 * t161 + Icges(1,3)) * V_base(6) + (t122 * t160 + t172 * t142 - t144 * t173 - t234 * t161 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t172 + Icges(2,6) * t173) * V_base(5) + (-Icges(2,5) * t173 + Icges(2,6) * t172) * V_base(6) + Icges(2,3) * t159 / 0.2e1) * t159;
T = t1;
