% Calculate kinetic energy for
% S5RPRPR4
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
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:21
% EndTime: 2022-01-23 09:22:22
% DurationCPUTime: 1.55s
% Computational Cost: add. (1262->243), mult. (946->323), div. (0->0), fcn. (726->10), ass. (0->126)
t232 = Icges(4,3) + Icges(5,3);
t169 = qJ(3) + pkin(9);
t158 = sin(t169);
t160 = cos(t169);
t172 = sin(qJ(3));
t174 = cos(qJ(3));
t231 = Icges(4,5) * t174 + Icges(5,5) * t160 - Icges(4,6) * t172 - Icges(5,6) * t158;
t170 = qJ(1) + pkin(8);
t159 = sin(t170);
t161 = cos(t170);
t213 = Icges(4,4) * t174;
t192 = -Icges(4,2) * t172 + t213;
t101 = -Icges(4,6) * t161 + t159 * t192;
t102 = Icges(4,6) * t159 + t161 * t192;
t214 = Icges(4,4) * t172;
t195 = Icges(4,1) * t174 - t214;
t103 = -Icges(4,5) * t161 + t159 * t195;
t104 = Icges(4,5) * t159 + t161 * t195;
t212 = Icges(5,4) * t158;
t121 = Icges(5,2) * t160 + t212;
t211 = Icges(5,4) * t160;
t124 = Icges(5,1) * t158 + t211;
t136 = -qJD(3) * t161 + V_base(5);
t137 = qJD(3) * t159 + V_base(4);
t141 = Icges(4,2) * t174 + t214;
t144 = Icges(4,1) * t172 + t213;
t163 = V_base(6) + qJD(1);
t191 = -Icges(5,2) * t158 + t211;
t92 = -Icges(5,6) * t161 + t159 * t191;
t93 = Icges(5,6) * t159 + t161 * t191;
t194 = Icges(5,1) * t160 - t212;
t94 = -Icges(5,5) * t161 + t159 * t194;
t95 = Icges(5,5) * t159 + t161 * t194;
t228 = (-t121 * t158 + t124 * t160 - t141 * t172 + t144 * t174) * t163 + (-t102 * t172 + t104 * t174 - t158 * t93 + t160 * t95) * t137 + (-t101 * t172 + t103 * t174 - t158 * t92 + t160 * t94) * t136;
t227 = (Icges(4,5) * t172 + Icges(5,5) * t158 + Icges(4,6) * t174 + Icges(5,6) * t160) * t163 + (t232 * t159 + t231 * t161) * t137 + (t231 * t159 - t232 * t161) * t136;
t173 = sin(qJ(1));
t223 = pkin(1) * t173;
t175 = cos(qJ(1));
t222 = pkin(1) * t175;
t221 = pkin(3) * t172;
t220 = pkin(4) * t158;
t219 = t174 * pkin(3);
t218 = -pkin(5) - qJ(2);
t216 = Icges(2,4) * t173;
t215 = Icges(3,4) * t159;
t162 = qJ(5) + t169;
t155 = sin(t162);
t210 = Icges(6,4) * t155;
t156 = cos(t162);
t209 = Icges(6,4) * t156;
t208 = pkin(4) * t160;
t206 = t163 * t222 + V_base(2);
t205 = V_base(5) * pkin(5) + V_base(1);
t130 = pkin(2) * t159 - pkin(6) * t161;
t202 = -t130 - t223;
t201 = V_base(5) * qJ(2) + t205;
t200 = V_base(4) * t223 + qJD(2) + V_base(3);
t80 = -qJ(4) * t161 + t159 * t219;
t199 = t202 - t80;
t198 = rSges(4,1) * t174 - rSges(4,2) * t172;
t197 = rSges(5,1) * t160 - rSges(5,2) * t158;
t196 = rSges(6,1) * t156 - rSges(6,2) * t155;
t193 = Icges(6,1) * t156 - t210;
t190 = -Icges(6,2) * t155 + t209;
t187 = Icges(6,5) * t156 - Icges(6,6) * t155;
t186 = qJD(4) * t159 + t136 * t221 + t201;
t110 = V_base(5) + (-qJD(3) - qJD(5)) * t161;
t111 = qJD(5) * t159 + t137;
t185 = (-Icges(6,3) * t161 + t159 * t187) * t110 + (Icges(6,3) * t159 + t161 * t187) * t111 + (Icges(6,5) * t155 + Icges(6,6) * t156) * t163;
t131 = pkin(2) * t161 + pkin(6) * t159;
t182 = t163 * t131 + t218 * V_base(4) + t206;
t181 = V_base(4) * t130 + (-t131 - t222) * V_base(5) + t200;
t180 = t137 * t80 + t181;
t81 = qJ(4) * t159 + t161 * t219;
t179 = -qJD(4) * t161 + t163 * t81 + t182;
t114 = Icges(6,2) * t156 + t210;
t115 = Icges(6,1) * t155 + t209;
t84 = -Icges(6,6) * t161 + t159 * t190;
t85 = Icges(6,6) * t159 + t161 * t190;
t86 = -Icges(6,5) * t161 + t159 * t193;
t87 = Icges(6,5) * t159 + t161 * t193;
t178 = (-t155 * t85 + t156 * t87) * t111 + (-t155 * t84 + t156 * t86) * t110 + (-t114 * t155 + t115 * t156) * t163;
t165 = Icges(2,4) * t175;
t154 = Icges(3,4) * t161;
t149 = rSges(2,1) * t175 - t173 * rSges(2,2);
t148 = t173 * rSges(2,1) + rSges(2,2) * t175;
t147 = rSges(4,1) * t172 + rSges(4,2) * t174;
t146 = Icges(2,1) * t175 - t216;
t145 = Icges(2,1) * t173 + t165;
t143 = -Icges(2,2) * t173 + t165;
t142 = Icges(2,2) * t175 + t216;
t135 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t134 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t133 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t161 - rSges(3,2) * t159;
t128 = rSges(3,1) * t159 + rSges(3,2) * t161;
t127 = rSges(5,1) * t158 + rSges(5,2) * t160;
t126 = Icges(3,1) * t161 - t215;
t125 = Icges(3,1) * t159 + t154;
t123 = -Icges(3,2) * t159 + t154;
t122 = Icges(3,2) * t161 + t215;
t116 = rSges(6,1) * t155 + rSges(6,2) * t156;
t108 = rSges(4,3) * t159 + t161 * t198;
t107 = -rSges(4,3) * t161 + t159 * t198;
t106 = V_base(5) * rSges(2,3) - t148 * t163 + t205;
t105 = t149 * t163 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t98 = t148 * V_base(4) - t149 * V_base(5) + V_base(3);
t97 = rSges(5,3) * t159 + t161 * t197;
t96 = -rSges(5,3) * t161 + t159 * t197;
t89 = rSges(6,3) * t159 + t161 * t196;
t88 = -rSges(6,3) * t161 + t159 * t196;
t78 = V_base(5) * rSges(3,3) + (-t128 - t223) * t163 + t201;
t77 = t129 * t163 + (-rSges(3,3) + t218) * V_base(4) + t206;
t75 = pkin(7) * t159 + t161 * t208;
t74 = -pkin(7) * t161 + t159 * t208;
t73 = V_base(4) * t128 + (-t129 - t222) * V_base(5) + t200;
t72 = t136 * t147 + (-t107 + t202) * t163 + t201;
t71 = t108 * t163 - t137 * t147 + t182;
t70 = t137 * t107 - t136 * t108 + t181;
t69 = t127 * t136 + (t199 - t96) * t163 + t186;
t68 = t163 * t97 + (-t127 - t221) * t137 + t179;
t67 = t137 * t96 + (-t81 - t97) * t136 + t180;
t66 = t136 * t220 + t110 * t116 + (t199 - t74 - t88) * t163 + t186;
t65 = -t111 * t116 + (t75 + t89) * t163 + (-t220 - t221) * t137 + t179;
t64 = -t110 * t89 + t111 * t88 + t137 * t74 + (-t75 - t81) * t136 + t180;
t1 = m(1) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t106 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t111 * (t185 * t159 + t178 * t161) / 0.2e1 + t110 * (t178 * t159 - t185 * t161) / 0.2e1 + (t228 * t159 - t227 * t161) * t136 / 0.2e1 + (t227 * t159 + t228 * t161) * t137 / 0.2e1 + ((-t122 * t159 + t125 * t161 - t173 * t142 + t145 * t175 + Icges(1,4)) * V_base(5) + (-t123 * t159 + t126 * t161 - t173 * t143 + t146 * t175 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t122 * t161 + t125 * t159 + t142 * t175 + t173 * t145 + Icges(1,2)) * V_base(5) + (t123 * t161 + t126 * t159 + t143 * t175 + t173 * t146 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t155 * t87 + t156 * t85) * t111 + (t155 * t86 + t156 * t84) * t110 + (t102 * t174 + t104 * t172 + t158 * t95 + t160 * t93) * t137 + (t101 * t174 + t103 * t172 + t158 * t94 + t160 * t92) * t136 + (t114 * t156 + t115 * t155 + t121 * t160 + t124 * t158 + t141 * t174 + t144 * t172 + Icges(2,3) + Icges(3,3)) * t163) * t163 / 0.2e1 + t163 * V_base(5) * (Icges(2,5) * t173 + Icges(3,5) * t159 + Icges(2,6) * t175 + Icges(3,6) * t161) + t163 * V_base(4) * (Icges(2,5) * t175 + Icges(3,5) * t161 - Icges(2,6) * t173 - Icges(3,6) * t159) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
