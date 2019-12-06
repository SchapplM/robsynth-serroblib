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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 19:00:00
% EndTime: 2019-12-05 19:00:01
% DurationCPUTime: 1.52s
% Computational Cost: add. (1342->242), mult. (970->352), div. (0->0), fcn. (750->10), ass. (0->130)
t172 = qJ(1) + qJ(2);
t163 = sin(t172);
t137 = qJD(3) * t163 + V_base(6);
t111 = qJD(4) * t163 + t137;
t108 = qJD(5) * t163 + t111;
t165 = cos(t172);
t138 = qJD(3) * t165 + V_base(5);
t112 = qJD(4) * t165 + t138;
t109 = qJD(5) * t165 + t112;
t171 = qJ(3) + qJ(4);
t166 = qJ(5) + t171;
t159 = cos(t166);
t158 = sin(t166);
t216 = Icges(6,4) * t158;
t115 = Icges(6,2) * t159 + t216;
t215 = Icges(6,4) * t159;
t116 = Icges(6,1) * t158 + t215;
t161 = V_base(4) + qJD(1);
t157 = qJD(2) + t161;
t195 = -Icges(6,2) * t158 + t215;
t83 = Icges(6,6) * t165 - t195 * t163;
t84 = Icges(6,6) * t163 + t195 * t165;
t198 = Icges(6,1) * t159 - t216;
t85 = Icges(6,5) * t165 - t198 * t163;
t86 = Icges(6,5) * t163 + t198 * t165;
t236 = (t115 * t158 - t116 * t159) * t157 + (t158 * t83 - t159 * t85) * t109 + (t158 * t84 - t159 * t86) * t108;
t164 = cos(t171);
t162 = sin(t171);
t218 = Icges(5,4) * t162;
t122 = Icges(5,2) * t164 + t218;
t217 = Icges(5,4) * t164;
t125 = Icges(5,1) * t162 + t217;
t196 = -Icges(5,2) * t162 + t217;
t91 = Icges(5,6) * t165 - t196 * t163;
t92 = Icges(5,6) * t163 + t196 * t165;
t199 = Icges(5,1) * t164 - t218;
t93 = Icges(5,5) * t165 - t199 * t163;
t94 = Icges(5,5) * t163 + t199 * t165;
t235 = (t122 * t162 - t125 * t164) * t157 + (t162 * t91 - t164 * t93) * t112 + (t162 * t92 - t164 * t94) * t111;
t173 = sin(qJ(3));
t175 = cos(qJ(3));
t219 = Icges(4,4) * t175;
t197 = -Icges(4,2) * t173 + t219;
t102 = Icges(4,6) * t165 - t197 * t163;
t103 = Icges(4,6) * t163 + t197 * t165;
t220 = Icges(4,4) * t173;
t200 = Icges(4,1) * t175 - t220;
t104 = Icges(4,5) * t165 - t200 * t163;
t105 = Icges(4,5) * t163 + t200 * t165;
t142 = Icges(4,2) * t175 + t220;
t145 = Icges(4,1) * t173 + t219;
t234 = (t142 * t173 - t145 * t175) * t157 + (t102 * t173 - t104 * t175) * t138 + (t103 * t173 - t105 * t175) * t137;
t233 = -pkin(5) - pkin(6);
t174 = sin(qJ(1));
t231 = pkin(1) * t174;
t176 = cos(qJ(1));
t230 = pkin(1) * t176;
t229 = pkin(3) * t173;
t228 = pkin(4) * t162;
t227 = t175 * pkin(3);
t132 = t165 * pkin(2) + t163 * pkin(7);
t80 = pkin(8) * t163 + t227 * t165;
t225 = -t132 - t80;
t224 = Icges(2,4) * t174;
t223 = Icges(2,4) * t176;
t222 = Icges(3,4) * t163;
t221 = Icges(3,4) * t165;
t214 = pkin(4) * t164;
t212 = V_base(6) * pkin(5) + V_base(2);
t209 = V_base(5) * t230 + V_base(6) * t231 + V_base(1);
t208 = rSges(4,1) * t175 - rSges(4,2) * t173;
t207 = rSges(5,1) * t164 - rSges(5,2) * t162;
t206 = rSges(6,1) * t159 - rSges(6,2) * t158;
t205 = -t161 * t231 + V_base(3);
t194 = Icges(4,5) * t175 - Icges(4,6) * t173;
t193 = Icges(5,5) * t164 - Icges(5,6) * t162;
t192 = Icges(6,5) * t159 - Icges(6,6) * t158;
t186 = V_base(6) * pkin(6) - t161 * t230 + t212;
t185 = (Icges(6,3) * t163 + t192 * t165) * t108 + (Icges(6,3) * t165 - t192 * t163) * t109 + (Icges(6,5) * t158 + Icges(6,6) * t159) * t157;
t184 = (Icges(5,3) * t163 + t193 * t165) * t111 + (Icges(5,3) * t165 - t193 * t163) * t112 + (Icges(5,5) * t162 + Icges(5,6) * t164) * t157;
t183 = (Icges(4,3) * t165 - t194 * t163) * t138 + (Icges(4,3) * t163 + t194 * t165) * t137 + (Icges(4,5) * t173 + Icges(4,6) * t175) * t157;
t182 = t137 * t229 + t186;
t131 = -t163 * pkin(2) + t165 * pkin(7);
t181 = -t131 * V_base(6) + V_base(5) * t132 + t209;
t180 = t157 * t131 + t233 * V_base(5) + t205;
t79 = pkin(8) * t165 - t227 * t163;
t179 = -t137 * t79 + t138 * t80 + t181;
t178 = -t138 * t229 + t157 * t79 + t180;
t150 = rSges(2,1) * t176 - rSges(2,2) * t174;
t149 = -rSges(2,1) * t174 - rSges(2,2) * t176;
t148 = rSges(4,1) * t173 + rSges(4,2) * t175;
t147 = Icges(2,1) * t176 - t224;
t146 = -Icges(2,1) * t174 - t223;
t144 = -Icges(2,2) * t174 + t223;
t143 = -Icges(2,2) * t176 - t224;
t136 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t135 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t134 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = rSges(3,1) * t165 - rSges(3,2) * t163;
t129 = -rSges(3,1) * t163 - rSges(3,2) * t165;
t128 = rSges(5,1) * t162 + rSges(5,2) * t164;
t127 = Icges(3,1) * t165 - t222;
t126 = -Icges(3,1) * t163 - t221;
t124 = -Icges(3,2) * t163 + t221;
t123 = -Icges(3,2) * t165 - t222;
t117 = rSges(6,1) * t158 + rSges(6,2) * t159;
t107 = rSges(4,3) * t163 + t208 * t165;
t106 = rSges(4,3) * t165 - t208 * t163;
t99 = V_base(6) * rSges(2,3) - t150 * t161 + t212;
t98 = t149 * t161 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t149 * V_base(6) + t150 * V_base(5) + V_base(1);
t96 = rSges(5,3) * t163 + t207 * t165;
t95 = rSges(5,3) * t165 - t207 * t163;
t88 = rSges(6,3) * t163 + t206 * t165;
t87 = rSges(6,3) * t165 - t206 * t163;
t76 = V_base(6) * rSges(3,3) - t130 * t157 + t186;
t75 = t129 * t157 + (-rSges(3,3) + t233) * V_base(5) + t205;
t74 = pkin(9) * t163 + t214 * t165;
t73 = pkin(9) * t165 - t214 * t163;
t72 = -t129 * V_base(6) + t130 * V_base(5) + t209;
t71 = t137 * t148 + (-t107 - t132) * t157 + t186;
t70 = t106 * t157 - t138 * t148 + t180;
t69 = -t106 * t137 + t107 * t138 + t181;
t68 = t111 * t128 + (-t96 + t225) * t157 + t182;
t67 = -t112 * t128 + t157 * t95 + t178;
t66 = -t111 * t95 + t112 * t96 + t179;
t65 = t111 * t228 + t108 * t117 + (-t74 - t88 + t225) * t157 + t182;
t64 = -t112 * t228 - t109 * t117 + (t73 + t87) * t157 + t178;
t63 = -t108 * t87 + t109 * t88 - t111 * t73 + t112 * t74 + t179;
t1 = m(1) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(2) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t138 * (t234 * t163 + t183 * t165) / 0.2e1 + t137 * (t183 * t163 - t234 * t165) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t112 * (t235 * t163 + t184 * t165) / 0.2e1 + t111 * (t184 * t163 - t235 * t165) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t109 * (t236 * t163 + t185 * t165) / 0.2e1 + t108 * (t185 * t163 - t236 * t165) / 0.2e1 + ((-t124 * t165 - t127 * t163 - t144 * t176 - t147 * t174 + Icges(1,6)) * V_base(6) + (-t165 * t123 - t163 * t126 - t176 * t143 - t174 * t146 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((-t163 * t124 + t165 * t127 - t174 * t144 + t176 * t147 + Icges(1,3)) * V_base(6) + (-t123 * t163 + t126 * t165 - t143 * t174 + t146 * t176 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t102 * t175 + t104 * t173) * t138 + (t103 * t175 + t105 * t173) * t137 + (t162 * t93 + t164 * t91) * t112 + (t162 * t94 + t164 * t92) * t111 + (t158 * t85 + t159 * t83) * t109 + (t158 * t86 + t159 * t84) * t108 + (t159 * t115 + t158 * t116 + t164 * t122 + t162 * t125 + t175 * t142 + t173 * t145 + Icges(3,3)) * t157) * t157 / 0.2e1 + t157 * V_base(6) * (Icges(3,5) * t165 - Icges(3,6) * t163) + t157 * V_base(5) * (-Icges(3,5) * t163 - Icges(3,6) * t165) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t174 - Icges(2,6) * t176) * V_base(5) + (Icges(2,5) * t176 - Icges(2,6) * t174) * V_base(6) + Icges(2,3) * t161 / 0.2e1) * t161;
T = t1;
