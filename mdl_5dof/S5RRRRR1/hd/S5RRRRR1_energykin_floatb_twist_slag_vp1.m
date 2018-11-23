% Calculate kinetic energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S5RRRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:51:07
% EndTime: 2018-11-16 14:51:09
% DurationCPUTime: 1.98s
% Computational Cost: add. (1360->277), mult. (1376->423), div. (0->0), fcn. (1214->10), ass. (0->149)
t182 = sin(qJ(1));
t238 = pkin(1) * t182;
t185 = cos(qJ(1));
t237 = pkin(1) * t185;
t181 = sin(qJ(2));
t236 = pkin(2) * t181;
t179 = qJ(2) + qJ(3);
t173 = sin(t179);
t235 = pkin(3) * t173;
t184 = cos(qJ(2));
t234 = t184 * pkin(2);
t233 = Icges(2,4) * t182;
t232 = Icges(3,4) * t181;
t231 = Icges(3,4) * t184;
t230 = Icges(4,4) * t173;
t174 = cos(t179);
t229 = Icges(4,4) * t174;
t176 = qJ(4) + t179;
t166 = sin(t176);
t228 = Icges(5,4) * t166;
t167 = cos(t176);
t227 = Icges(5,4) * t167;
t226 = t166 * t182;
t225 = t166 * t185;
t180 = sin(qJ(5));
t224 = t180 * t185;
t223 = t182 * t180;
t183 = cos(qJ(5));
t222 = t182 * t183;
t221 = t183 * t185;
t220 = pkin(3) * t174;
t219 = qJD(5) * t166;
t218 = -qJD(2) - qJD(3);
t217 = V_base(5) * pkin(5) + V_base(1);
t163 = qJD(2) * t185 + V_base(5);
t169 = V_base(6) + qJD(1);
t136 = t234 * t182;
t214 = -t136 - t238;
t142 = qJD(3) * t185 + t163;
t102 = t220 * t182;
t213 = -t102 + t214;
t212 = pkin(4) * t167 + pkin(6) * t166;
t211 = rSges(3,1) * t184 - rSges(3,2) * t181;
t210 = rSges(4,1) * t174 - rSges(4,2) * t173;
t209 = rSges(5,1) * t167 - rSges(5,2) * t166;
t129 = qJD(4) * t185 + t142;
t208 = Icges(3,1) * t184 - t232;
t207 = Icges(4,1) * t174 - t230;
t206 = Icges(5,1) * t167 - t228;
t205 = -Icges(3,2) * t181 + t231;
t204 = -Icges(4,2) * t173 + t229;
t203 = -Icges(5,2) * t166 + t227;
t202 = Icges(3,5) * t184 - Icges(3,6) * t181;
t201 = Icges(4,5) * t174 - Icges(4,6) * t173;
t200 = Icges(5,5) * t167 - Icges(5,6) * t166;
t199 = -t163 * t236 + t217;
t198 = -V_base(4) * pkin(5) + t169 * t237 + V_base(2);
t197 = -t237 * V_base(5) + V_base(4) * t238 + V_base(3);
t130 = V_base(4) + (-qJD(4) + t218) * t182;
t196 = (Icges(5,3) * t185 + t182 * t200) * t129 + (-Icges(5,3) * t182 + t185 * t200) * t130 + (-Icges(5,5) * t166 - Icges(5,6) * t167) * t169;
t143 = t182 * t218 + V_base(4);
t195 = (Icges(4,3) * t185 + t182 * t201) * t142 + (-Icges(4,3) * t182 + t185 * t201) * t143 + (-Icges(4,5) * t173 - Icges(4,6) * t174) * t169;
t164 = -qJD(2) * t182 + V_base(4);
t194 = (Icges(3,3) * t185 + t182 * t202) * t163 + (-Icges(3,3) * t182 + t185 * t202) * t164 + (-Icges(3,5) * t181 - Icges(3,6) * t184) * t169;
t137 = t234 * t185;
t193 = t169 * t137 + t164 * t236 + t198;
t192 = -t142 * t235 + t199;
t103 = t220 * t185;
t191 = t169 * t103 + t143 * t235 + t193;
t190 = t164 * t136 - t163 * t137 + t197;
t189 = t143 * t102 - t142 * t103 + t190;
t132 = -Icges(5,2) * t167 - t228;
t133 = -Icges(5,1) * t166 - t227;
t96 = Icges(5,6) * t185 + t182 * t203;
t97 = -Icges(5,6) * t182 + t185 * t203;
t98 = Icges(5,5) * t185 + t182 * t206;
t99 = -Icges(5,5) * t182 + t185 * t206;
t188 = (-t166 * t97 + t167 * t99) * t130 + (-t166 * t96 + t167 * t98) * t129 + (-t132 * t166 + t133 * t167) * t169;
t106 = Icges(4,6) * t185 + t182 * t204;
t107 = -Icges(4,6) * t182 + t185 * t204;
t108 = Icges(4,5) * t185 + t182 * t207;
t109 = -Icges(4,5) * t182 + t185 * t207;
t139 = -Icges(4,2) * t174 - t230;
t140 = -Icges(4,1) * t173 - t229;
t187 = (-t107 * t173 + t109 * t174) * t143 + (-t106 * t173 + t108 * t174) * t142 + (-t139 * t173 + t140 * t174) * t169;
t117 = Icges(3,6) * t185 + t182 * t205;
t118 = -Icges(3,6) * t182 + t185 * t205;
t119 = Icges(3,5) * t185 + t182 * t208;
t120 = -Icges(3,5) * t182 + t185 * t208;
t153 = -Icges(3,2) * t184 - t232;
t156 = -Icges(3,1) * t181 - t231;
t186 = (-t118 * t181 + t120 * t184) * t164 + (-t117 * t181 + t119 * t184) * t163 + (-t153 * t181 + t156 * t184) * t169;
t175 = Icges(2,4) * t185;
t161 = rSges(2,1) * t185 - t182 * rSges(2,2);
t160 = t182 * rSges(2,1) + rSges(2,2) * t185;
t159 = -rSges(3,1) * t181 - rSges(3,2) * t184;
t158 = Icges(2,1) * t185 - t233;
t157 = Icges(2,1) * t182 + t175;
t155 = -Icges(2,2) * t182 + t175;
t154 = Icges(2,2) * t185 + t233;
t149 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t148 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t147 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t144 = qJD(5) * t167 + t169;
t141 = -rSges(4,1) * t173 - rSges(4,2) * t174;
t135 = -pkin(4) * t166 + pkin(6) * t167;
t134 = -rSges(5,1) * t166 - rSges(5,2) * t167;
t128 = t167 * t221 - t223;
t127 = -t167 * t224 - t222;
t126 = t167 * t222 + t224;
t125 = -t167 * t223 + t221;
t123 = -t182 * rSges(3,3) + t185 * t211;
t122 = rSges(3,3) * t185 + t182 * t211;
t114 = t212 * t185;
t113 = t212 * t182;
t111 = -t182 * rSges(4,3) + t185 * t210;
t110 = rSges(4,3) * t185 + t182 * t210;
t101 = -t182 * rSges(5,3) + t185 * t209;
t100 = rSges(5,3) * t185 + t182 * t209;
t93 = V_base(5) * rSges(2,3) - t160 * t169 + t217;
t92 = t161 * t169 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t91 = t185 * t219 + t130;
t90 = t182 * t219 + t129;
t89 = t160 * V_base(4) - t161 * V_base(5) + V_base(3);
t88 = rSges(6,3) * t167 + (-rSges(6,1) * t183 + rSges(6,2) * t180) * t166;
t87 = Icges(6,5) * t167 + (-Icges(6,1) * t183 + Icges(6,4) * t180) * t166;
t86 = Icges(6,6) * t167 + (-Icges(6,4) * t183 + Icges(6,2) * t180) * t166;
t85 = Icges(6,3) * t167 + (-Icges(6,5) * t183 + Icges(6,6) * t180) * t166;
t82 = t128 * rSges(6,1) + t127 * rSges(6,2) + rSges(6,3) * t225;
t81 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t226;
t80 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t225;
t79 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t226;
t78 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t225;
t77 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t226;
t76 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t225;
t75 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t226;
t74 = t159 * t163 + (-t122 - t238) * t169 + t217;
t73 = t123 * t169 - t159 * t164 + t198;
t72 = t164 * t122 - t163 * t123 + t197;
t71 = t141 * t142 + (-t110 + t214) * t169 + t199;
t70 = t111 * t169 - t141 * t143 + t193;
t69 = t143 * t110 - t142 * t111 + t190;
t68 = t129 * t134 + (-t100 + t213) * t169 + t192;
t67 = t101 * t169 - t130 * t134 + t191;
t66 = t130 * t100 - t129 * t101 + t189;
t65 = t129 * t135 - t144 * t81 + t88 * t90 + (-t113 + t213) * t169 + t192;
t64 = t114 * t169 - t130 * t135 + t144 * t82 - t88 * t91 + t191;
t63 = t130 * t113 - t129 * t114 + t91 * t81 - t90 * t82 + t189;
t1 = m(1) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(2) * (t89 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + t164 * (-t182 * t194 + t186 * t185) / 0.2e1 + t163 * (t186 * t182 + t185 * t194) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t143 * (-t182 * t195 + t187 * t185) / 0.2e1 + t142 * (t187 * t182 + t185 * t195) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t130 * (-t182 * t196 + t185 * t188) / 0.2e1 + t129 * (t182 * t188 + t185 * t196) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t91 * ((t127 * t78 + t128 * t80 + t225 * t76) * t91 + (t127 * t77 + t128 * t79 + t225 * t75) * t90 + (t127 * t86 + t128 * t87 + t225 * t85) * t144) / 0.2e1 + t90 * ((t125 * t78 + t126 * t80 + t226 * t76) * t91 + (t125 * t77 + t126 * t79 + t226 * t75) * t90 + (t125 * t86 + t126 * t87 + t226 * t85) * t144) / 0.2e1 + t144 * ((t85 * t144 + t75 * t90 + t76 * t91) * t167 + ((t180 * t78 - t183 * t80) * t91 + (t180 * t77 - t183 * t79) * t90 + (t180 * t86 - t183 * t87) * t144) * t166) / 0.2e1 + ((-t182 * t154 + t157 * t185 + Icges(1,4)) * V_base(5) + (-t182 * t155 + t185 * t158 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t154 * t185 + t182 * t157 + Icges(1,2)) * V_base(5) + (t155 * t185 + t182 * t158 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t118 * t184 - t120 * t181) * t164 + (-t117 * t184 - t119 * t181) * t163 + (-t107 * t174 - t109 * t173) * t143 + (-t106 * t174 - t108 * t173) * t142 + (-t166 * t99 - t167 * t97) * t130 + (-t166 * t98 - t167 * t96) * t129 + (-t167 * t132 - t166 * t133 - t174 * t139 - t173 * t140 - t184 * t153 - t181 * t156 + Icges(2,3)) * t169) * t169 / 0.2e1 + V_base(4) * t169 * (Icges(2,5) * t185 - Icges(2,6) * t182) + V_base(5) * t169 * (Icges(2,5) * t182 + Icges(2,6) * t185) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
