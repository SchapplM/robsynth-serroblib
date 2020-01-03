% Calculate kinetic energy for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:14:44
% EndTime: 2020-01-03 12:14:45
% DurationCPUTime: 1.56s
% Computational Cost: add. (1342->244), mult. (970->352), div. (0->0), fcn. (750->10), ass. (0->131)
t170 = qJ(1) + qJ(2);
t160 = sin(t170);
t212 = -qJD(3) - qJD(4);
t207 = -qJD(5) + t212;
t108 = t207 * t160 + V_base(6);
t162 = cos(t170);
t109 = t207 * t162 + V_base(5);
t169 = qJ(3) + qJ(4);
t164 = qJ(5) + t169;
t156 = cos(t164);
t155 = sin(t164);
t216 = Icges(6,4) * t155;
t115 = Icges(6,2) * t156 + t216;
t215 = Icges(6,4) * t156;
t116 = Icges(6,1) * t155 + t215;
t158 = V_base(4) + qJD(1);
t154 = qJD(2) + t158;
t193 = -Icges(6,2) * t155 + t215;
t83 = -Icges(6,6) * t162 + t193 * t160;
t84 = -Icges(6,6) * t160 - t193 * t162;
t196 = Icges(6,1) * t156 - t216;
t85 = -Icges(6,5) * t162 + t196 * t160;
t86 = -Icges(6,5) * t160 - t196 * t162;
t233 = (t115 * t155 - t116 * t156) * t154 + (t155 * t83 - t156 * t85) * t109 + (t155 * t84 - t156 * t86) * t108;
t111 = t212 * t160 + V_base(6);
t112 = t212 * t162 + V_base(5);
t161 = cos(t169);
t159 = sin(t169);
t218 = Icges(5,4) * t159;
t122 = Icges(5,2) * t161 + t218;
t217 = Icges(5,4) * t161;
t125 = Icges(5,1) * t159 + t217;
t194 = -Icges(5,2) * t159 + t217;
t91 = -Icges(5,6) * t162 + t194 * t160;
t92 = -Icges(5,6) * t160 - t194 * t162;
t197 = Icges(5,1) * t161 - t218;
t93 = -Icges(5,5) * t162 + t197 * t160;
t94 = -Icges(5,5) * t160 - t197 * t162;
t232 = (t122 * t159 - t125 * t161) * t154 + (t159 * t91 - t161 * t93) * t112 + (t159 * t92 - t161 * t94) * t111;
t171 = sin(qJ(3));
t173 = cos(qJ(3));
t219 = Icges(4,4) * t173;
t195 = -Icges(4,2) * t171 + t219;
t102 = -Icges(4,6) * t162 + t195 * t160;
t103 = -Icges(4,6) * t160 - t195 * t162;
t220 = Icges(4,4) * t171;
t198 = Icges(4,1) * t173 - t220;
t104 = -Icges(4,5) * t162 + t198 * t160;
t105 = -Icges(4,5) * t160 - t198 * t162;
t137 = -qJD(3) * t160 + V_base(6);
t138 = -qJD(3) * t162 + V_base(5);
t142 = Icges(4,2) * t173 + t220;
t145 = Icges(4,1) * t171 + t219;
t231 = (t142 * t171 - t145 * t173) * t154 + (t102 * t171 - t104 * t173) * t138 + (t103 * t171 - t105 * t173) * t137;
t230 = -pkin(5) - pkin(6);
t228 = pkin(1) * t158;
t227 = pkin(3) * t171;
t226 = pkin(4) * t159;
t225 = t173 * pkin(3);
t132 = -t162 * pkin(2) - t160 * pkin(7);
t80 = -pkin(8) * t160 - t225 * t162;
t223 = -t132 - t80;
t174 = cos(qJ(1));
t222 = Icges(2,4) * t174;
t221 = Icges(3,4) * t162;
t214 = pkin(4) * t161;
t172 = sin(qJ(1));
t211 = t172 * t228 + V_base(3);
t210 = V_base(6) * pkin(5) + V_base(2);
t206 = V_base(6) * pkin(6) + t174 * t228 + t210;
t205 = rSges(4,1) * t173 - rSges(4,2) * t171;
t204 = rSges(5,1) * t161 - rSges(5,2) * t159;
t203 = rSges(6,1) * t156 - rSges(6,2) * t155;
t192 = Icges(4,5) * t173 - Icges(4,6) * t171;
t191 = Icges(5,5) * t161 - Icges(5,6) * t159;
t190 = Icges(6,5) * t156 - Icges(6,6) * t155;
t184 = t137 * t227 + t206;
t183 = -(-Icges(6,3) * t160 - t190 * t162) * t108 - (-Icges(6,3) * t162 + t190 * t160) * t109 - (Icges(6,5) * t155 + Icges(6,6) * t156) * t154;
t182 = -(-Icges(5,3) * t160 - t191 * t162) * t111 - (-Icges(5,3) * t162 + t191 * t160) * t112 - (Icges(5,5) * t159 + Icges(5,6) * t161) * t154;
t181 = -(-Icges(4,3) * t162 + t192 * t160) * t138 - (-Icges(4,3) * t160 - t192 * t162) * t137 - (Icges(4,5) * t171 + Icges(4,6) * t173) * t154;
t131 = t160 * pkin(2) - t162 * pkin(7);
t180 = t154 * t131 + t230 * V_base(5) + t211;
t179 = V_base(1) + (-V_base(6) * t172 - t174 * V_base(5)) * pkin(1);
t79 = -pkin(8) * t162 + t225 * t160;
t178 = -t138 * t227 + t154 * t79 + t180;
t177 = -t131 * V_base(6) + V_base(5) * t132 + t179;
t176 = -t137 * t79 + t138 * t80 + t177;
t163 = Icges(2,4) * t172;
t153 = Icges(3,4) * t160;
t150 = -rSges(2,1) * t174 + rSges(2,2) * t172;
t149 = rSges(2,1) * t172 + rSges(2,2) * t174;
t148 = rSges(4,1) * t171 + rSges(4,2) * t173;
t147 = -Icges(2,1) * t174 + t163;
t146 = Icges(2,1) * t172 + t222;
t144 = Icges(2,2) * t172 - t222;
t143 = Icges(2,2) * t174 + t163;
t136 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t135 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t134 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = -rSges(3,1) * t162 + rSges(3,2) * t160;
t129 = rSges(3,1) * t160 + rSges(3,2) * t162;
t128 = rSges(5,1) * t159 + rSges(5,2) * t161;
t127 = -Icges(3,1) * t162 + t153;
t126 = Icges(3,1) * t160 + t221;
t124 = Icges(3,2) * t160 - t221;
t123 = Icges(3,2) * t162 + t153;
t117 = rSges(6,1) * t155 + rSges(6,2) * t156;
t107 = -rSges(4,3) * t160 - t205 * t162;
t106 = -rSges(4,3) * t162 + t205 * t160;
t99 = V_base(6) * rSges(2,3) - t150 * t158 + t210;
t98 = t149 * t158 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t149 * V_base(6) + t150 * V_base(5) + V_base(1);
t96 = -rSges(5,3) * t160 - t204 * t162;
t95 = -rSges(5,3) * t162 + t204 * t160;
t88 = -rSges(6,3) * t160 - t203 * t162;
t87 = -rSges(6,3) * t162 + t203 * t160;
t76 = V_base(6) * rSges(3,3) - t130 * t154 + t206;
t75 = t129 * t154 + (-rSges(3,3) + t230) * V_base(5) + t211;
t74 = -pkin(9) * t160 - t214 * t162;
t73 = -pkin(9) * t162 + t214 * t160;
t72 = -t129 * V_base(6) + t130 * V_base(5) + t179;
t71 = t137 * t148 + (-t107 - t132) * t154 + t206;
t70 = t106 * t154 - t138 * t148 + t180;
t69 = -t106 * t137 + t107 * t138 + t177;
t68 = t111 * t128 + (-t96 + t223) * t154 + t184;
t67 = -t112 * t128 + t154 * t95 + t178;
t66 = -t111 * t95 + t112 * t96 + t176;
t65 = t111 * t226 + t108 * t117 + (-t74 - t88 + t223) * t154 + t184;
t64 = -t112 * t226 - t109 * t117 + (t73 + t87) * t154 + t178;
t63 = -t108 * t87 + t109 * t88 - t111 * t73 + t112 * t74 + t176;
t1 = m(1) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(2) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t138 * (-t231 * t160 + t181 * t162) / 0.2e1 + t137 * (t181 * t160 + t231 * t162) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t112 * (-t232 * t160 + t182 * t162) / 0.2e1 + t111 * (t182 * t160 + t232 * t162) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t109 * (-t233 * t160 + t183 * t162) / 0.2e1 + t108 * (t183 * t160 + t233 * t162) / 0.2e1 + ((t124 * t162 + t127 * t160 + t144 * t174 + t147 * t172 + Icges(1,6)) * V_base(6) + (t162 * t123 + t160 * t126 + t174 * t143 + t172 * t146 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t160 * t124 - t162 * t127 + t172 * t144 - t174 * t147 + Icges(1,3)) * V_base(6) + (t123 * t160 - t126 * t162 + t143 * t172 - t146 * t174 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t102 * t173 + t104 * t171) * t138 + (t103 * t173 + t105 * t171) * t137 + (t159 * t93 + t161 * t91) * t112 + (t159 * t94 + t161 * t92) * t111 + (t155 * t85 + t156 * t83) * t109 + (t155 * t86 + t156 * t84) * t108 + (t156 * t115 + t155 * t116 + t161 * t122 + t159 * t125 + t173 * t142 + t171 * t145 + Icges(3,3)) * t154) * t154 / 0.2e1 + t154 * V_base(6) * (-Icges(3,5) * t162 + Icges(3,6) * t160) + t154 * V_base(5) * (Icges(3,5) * t160 + Icges(3,6) * t162) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t172 + Icges(2,6) * t174) * V_base(5) + (-Icges(2,5) * t174 + Icges(2,6) * t172) * V_base(6) + Icges(2,3) * t158 / 0.2e1) * t158;
T = t1;
