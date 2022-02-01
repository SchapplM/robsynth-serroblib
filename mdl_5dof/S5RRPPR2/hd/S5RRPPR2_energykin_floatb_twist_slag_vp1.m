% Calculate kinetic energy for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:23
% EndTime: 2022-01-20 10:05:25
% DurationCPUTime: 1.62s
% Computational Cost: add. (1211->254), mult. (946->341), div. (0->0), fcn. (780->10), ass. (0->125)
t166 = qJ(1) + qJ(2);
t158 = pkin(8) + t166;
t155 = sin(t158);
t156 = cos(t158);
t160 = sin(t166);
t161 = cos(t166);
t213 = Icges(3,5) * t160 + Icges(4,5) * t155 + Icges(3,6) * t161 + Icges(4,6) * t156;
t212 = Icges(3,5) * t161 + Icges(4,5) * t156 - Icges(3,6) * t160 - Icges(4,6) * t155;
t211 = -pkin(5) - pkin(6);
t170 = sin(qJ(1));
t209 = pkin(1) * t170;
t172 = cos(qJ(1));
t208 = pkin(1) * t172;
t207 = pkin(2) * t160;
t206 = pkin(2) * t161;
t205 = Icges(2,4) * t170;
t204 = Icges(3,4) * t160;
t203 = Icges(4,4) * t155;
t167 = sin(pkin(9));
t202 = Icges(5,4) * t167;
t168 = cos(pkin(9));
t201 = Icges(5,4) * t168;
t200 = t155 * t167;
t199 = t156 * t167;
t169 = sin(qJ(5));
t198 = t168 * t169;
t171 = cos(qJ(5));
t197 = t168 * t171;
t196 = qJD(5) * t167;
t195 = -qJ(3) + t211;
t159 = V_base(6) + qJD(1);
t194 = t159 * t208 + V_base(2);
t193 = V_base(4) * t209 + V_base(3);
t192 = V_base(5) * pkin(5) + V_base(1);
t117 = pkin(3) * t155 - qJ(4) * t156;
t189 = -t117 - t207;
t157 = qJD(2) + t159;
t188 = t157 * t206 + t194;
t187 = -t206 - t208;
t186 = pkin(4) * t168 + pkin(7) * t167;
t185 = V_base(4) * t207 + qJD(3) + t193;
t184 = rSges(5,1) * t168 - rSges(5,2) * t167;
t183 = Icges(5,1) * t168 - t202;
t182 = -Icges(5,2) * t167 + t201;
t181 = Icges(5,5) * t168 - Icges(5,6) * t167;
t119 = pkin(3) * t156 + qJ(4) * t155;
t180 = -t119 + t187;
t179 = V_base(4) * t117 + t185;
t178 = V_base(5) * pkin(6) - t159 * t209 + t192;
t177 = -qJD(4) * t156 + t157 * t119 + t188;
t176 = V_base(5) * qJ(3) + t178;
t175 = qJD(4) * t155 + t176;
t174 = (Icges(5,5) * t167 + Icges(5,6) * t168) * t157 + (-Icges(5,3) * t156 + t181 * t155) * V_base(5) + (Icges(5,3) * t155 + t181 * t156) * V_base(4);
t137 = Icges(5,2) * t168 + t202;
t138 = Icges(5,1) * t167 + t201;
t90 = -Icges(5,6) * t156 + t182 * t155;
t91 = Icges(5,6) * t155 + t182 * t156;
t92 = -Icges(5,5) * t156 + t183 * t155;
t93 = Icges(5,5) * t155 + t183 * t156;
t173 = (-t167 * t91 + t168 * t93) * V_base(4) + (-t167 * t90 + t168 * t92) * V_base(5) + (-t137 * t167 + t138 * t168) * t157;
t163 = Icges(2,4) * t172;
t154 = Icges(3,4) * t161;
t152 = Icges(4,4) * t156;
t148 = rSges(2,1) * t172 - t170 * rSges(2,2);
t147 = t170 * rSges(2,1) + rSges(2,2) * t172;
t146 = Icges(2,1) * t172 - t205;
t145 = Icges(2,1) * t170 + t163;
t144 = -Icges(2,2) * t170 + t163;
t143 = Icges(2,2) * t172 + t205;
t140 = pkin(4) * t167 - pkin(7) * t168;
t139 = rSges(5,1) * t167 + rSges(5,2) * t168;
t135 = -qJD(5) * t168 + t157;
t134 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t133 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t132 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = rSges(3,1) * t161 - rSges(3,2) * t160;
t129 = rSges(3,1) * t160 + rSges(3,2) * t161;
t128 = Icges(3,1) * t161 - t204;
t127 = Icges(3,1) * t160 + t154;
t126 = -Icges(3,2) * t160 + t154;
t125 = Icges(3,2) * t161 + t204;
t122 = t156 * t196 + V_base(4);
t121 = t155 * t196 + V_base(5);
t120 = rSges(4,1) * t156 - rSges(4,2) * t155;
t118 = rSges(4,1) * t155 + rSges(4,2) * t156;
t116 = Icges(4,1) * t156 - t203;
t115 = Icges(4,1) * t155 + t152;
t114 = -Icges(4,2) * t155 + t152;
t113 = Icges(4,2) * t156 + t203;
t109 = -rSges(6,3) * t168 + (rSges(6,1) * t171 - rSges(6,2) * t169) * t167;
t108 = -Icges(6,5) * t168 + (Icges(6,1) * t171 - Icges(6,4) * t169) * t167;
t107 = -Icges(6,6) * t168 + (Icges(6,4) * t171 - Icges(6,2) * t169) * t167;
t106 = -Icges(6,3) * t168 + (Icges(6,5) * t171 - Icges(6,6) * t169) * t167;
t105 = t186 * t156;
t104 = t186 * t155;
t103 = t155 * t169 + t156 * t197;
t102 = t155 * t171 - t156 * t198;
t101 = t155 * t197 - t156 * t169;
t100 = -t155 * t198 - t156 * t171;
t98 = V_base(5) * rSges(2,3) - t147 * t159 + t192;
t97 = t148 * t159 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t96 = t147 * V_base(4) - t148 * V_base(5) + V_base(3);
t95 = rSges(5,3) * t155 + t184 * t156;
t94 = -rSges(5,3) * t156 + t184 * t155;
t87 = V_base(5) * rSges(3,3) - t129 * t157 + t178;
t86 = t130 * t157 + (-rSges(3,3) + t211) * V_base(4) + t194;
t85 = V_base(4) * t129 + (-t130 - t208) * V_base(5) + t193;
t84 = rSges(6,1) * t103 + rSges(6,2) * t102 + rSges(6,3) * t199;
t83 = rSges(6,1) * t101 + rSges(6,2) * t100 + rSges(6,3) * t200;
t82 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t199;
t81 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t200;
t80 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t199;
t79 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t200;
t78 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t199;
t77 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t200;
t76 = V_base(5) * rSges(4,3) + (-t118 - t207) * t157 + t176;
t75 = t120 * t157 + (-rSges(4,3) + t195) * V_base(4) + t188;
t74 = V_base(4) * t118 + (-t120 + t187) * V_base(5) + t185;
t73 = t139 * V_base(5) + (t189 - t94) * t157 + t175;
t72 = t157 * t95 + (-t139 + t195) * V_base(4) + t177;
t71 = V_base(4) * t94 + (t180 - t95) * V_base(5) + t179;
t70 = t109 * t121 - t135 * t83 + t140 * V_base(5) + (-t104 + t189) * t157 + t175;
t69 = t105 * t157 - t109 * t122 + t135 * t84 + (-t140 + t195) * V_base(4) + t177;
t68 = V_base(4) * t104 - t121 * t84 + t122 * t83 + (-t105 + t180) * V_base(5) + t179;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + t122 * ((t102 * t80 + t103 * t82 + t78 * t199) * t122 + (t102 * t79 + t103 * t81 + t77 * t199) * t121 + (t102 * t107 + t103 * t108 + t106 * t199) * t135) / 0.2e1 + t121 * ((t100 * t80 + t101 * t82 + t78 * t200) * t122 + (t100 * t79 + t101 * t81 + t77 * t200) * t121 + (t100 * t107 + t101 * t108 + t106 * t200) * t135) / 0.2e1 + t135 * ((-t106 * t135 - t121 * t77 - t122 * t78) * t168 + ((-t169 * t80 + t171 * t82) * t122 + (-t169 * t79 + t171 * t81) * t121 + (-t107 * t169 + t108 * t171) * t135) * t167) / 0.2e1 + ((t167 * t92 + t168 * t90 + t213) * V_base(5) + (t167 * t93 + t168 * t91 + t212) * V_base(4) + (t168 * t137 + t167 * t138 + Icges(3,3) + Icges(4,3)) * t157) * t157 / 0.2e1 + (t174 * t155 + t173 * t156 + t212 * t157 + (-t113 * t155 + t115 * t156 - t125 * t160 + t127 * t161 - t170 * t143 + t145 * t172 + Icges(1,4)) * V_base(5) + (-t155 * t114 + t156 * t116 - t160 * t126 + t161 * t128 - t170 * t144 + t172 * t146 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t173 * t155 - t174 * t156 + t213 * t157 + (t156 * t113 + t155 * t115 + t161 * t125 + t160 * t127 + t172 * t143 + t170 * t145 + Icges(1,2)) * V_base(5) + (t114 * t156 + t116 * t155 + t126 * t161 + t128 * t160 + t144 * t172 + t170 * t146 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t170 + Icges(2,6) * t172) * V_base(5) + (Icges(2,5) * t172 - Icges(2,6) * t170) * V_base(4) + Icges(2,3) * t159 / 0.2e1) * t159;
T = t1;
