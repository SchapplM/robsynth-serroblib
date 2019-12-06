% Calculate kinetic energy for
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:16
% EndTime: 2019-12-05 16:17:17
% DurationCPUTime: 1.81s
% Computational Cost: add. (1187->253), mult. (946->345), div. (0->0), fcn. (780->10), ass. (0->124)
t170 = cos(pkin(8));
t210 = pkin(1) * t170;
t162 = V_base(6) + qJD(2);
t209 = pkin(2) * t162;
t208 = -pkin(5) - qJ(1);
t168 = sin(pkin(8));
t207 = Icges(2,4) * t168;
t166 = pkin(8) + qJ(2);
t158 = sin(t166);
t206 = Icges(3,4) * t158;
t161 = qJ(3) + t166;
t155 = sin(t161);
t205 = Icges(4,4) * t155;
t167 = sin(pkin(9));
t204 = Icges(5,4) * t167;
t169 = cos(pkin(9));
t203 = Icges(5,4) * t169;
t202 = t155 * t167;
t156 = cos(t161);
t201 = t156 * t167;
t171 = sin(qJ(5));
t200 = t169 * t171;
t172 = cos(qJ(5));
t199 = t169 * t172;
t198 = qJD(5) * t167;
t197 = -pkin(6) + t208;
t190 = pkin(1) * V_base(6);
t196 = t170 * t190 + V_base(2);
t195 = V_base(5) * qJ(1) + V_base(1);
t191 = qJD(1) + V_base(3);
t159 = cos(t166);
t189 = t159 * t209 + t196;
t188 = V_base(4) * t168 * pkin(1) + t191;
t157 = qJD(3) + t162;
t187 = -pkin(2) * t159 - t210;
t186 = pkin(4) * t169 + pkin(7) * t167;
t185 = V_base(4) * pkin(2) * t158 + t188;
t184 = rSges(5,1) * t169 - rSges(5,2) * t167;
t183 = Icges(5,1) * t169 - t204;
t182 = -Icges(5,2) * t167 + t203;
t181 = Icges(5,5) * t169 - Icges(5,6) * t167;
t119 = pkin(3) * t156 + qJ(4) * t155;
t180 = -t119 + t187;
t117 = pkin(3) * t155 - qJ(4) * t156;
t179 = V_base(4) * t117 + t185;
t178 = -qJD(4) * t156 + t157 * t119 + t189;
t177 = V_base(5) * pkin(5) - t168 * t190 + t195;
t176 = (Icges(5,5) * t167 + Icges(5,6) * t169) * t157 + (-Icges(5,3) * t156 + t155 * t181) * V_base(5) + (Icges(5,3) * t155 + t156 * t181) * V_base(4);
t175 = V_base(5) * pkin(6) - t158 * t209 + t177;
t174 = qJD(4) * t155 + t175;
t139 = Icges(5,2) * t169 + t204;
t142 = Icges(5,1) * t167 + t203;
t91 = -Icges(5,6) * t156 + t155 * t182;
t92 = Icges(5,6) * t155 + t156 * t182;
t93 = -Icges(5,5) * t156 + t155 * t183;
t94 = Icges(5,5) * t155 + t156 * t183;
t173 = (-t167 * t92 + t169 * t94) * V_base(4) + (-t167 * t91 + t169 * t93) * V_base(5) + (-t139 * t167 + t142 * t169) * t157;
t160 = Icges(2,4) * t170;
t154 = Icges(3,4) * t159;
t151 = Icges(4,4) * t156;
t148 = pkin(4) * t167 - pkin(7) * t169;
t147 = rSges(2,1) * t170 - rSges(2,2) * t168;
t146 = rSges(2,1) * t168 + rSges(2,2) * t170;
t145 = rSges(5,1) * t167 + rSges(5,2) * t169;
t144 = Icges(2,1) * t170 - t207;
t143 = Icges(2,1) * t168 + t160;
t141 = -Icges(2,2) * t168 + t160;
t140 = Icges(2,2) * t170 + t207;
t134 = -qJD(5) * t169 + t157;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = rSges(3,1) * t159 - rSges(3,2) * t158;
t129 = rSges(3,1) * t158 + rSges(3,2) * t159;
t128 = Icges(3,1) * t159 - t206;
t127 = Icges(3,1) * t158 + t154;
t126 = -Icges(3,2) * t158 + t154;
t125 = Icges(3,2) * t159 + t206;
t122 = t156 * t198 + V_base(4);
t121 = t155 * t198 + V_base(5);
t120 = rSges(4,1) * t156 - rSges(4,2) * t155;
t118 = rSges(4,1) * t155 + rSges(4,2) * t156;
t116 = Icges(4,1) * t156 - t205;
t115 = Icges(4,1) * t155 + t151;
t114 = -Icges(4,2) * t155 + t151;
t113 = Icges(4,2) * t156 + t205;
t112 = Icges(4,5) * t156 - Icges(4,6) * t155;
t111 = Icges(4,5) * t155 + Icges(4,6) * t156;
t109 = -t169 * rSges(6,3) + (rSges(6,1) * t172 - rSges(6,2) * t171) * t167;
t108 = -Icges(6,5) * t169 + (Icges(6,1) * t172 - Icges(6,4) * t171) * t167;
t107 = -Icges(6,6) * t169 + (Icges(6,4) * t172 - Icges(6,2) * t171) * t167;
t106 = -Icges(6,3) * t169 + (Icges(6,5) * t172 - Icges(6,6) * t171) * t167;
t105 = t186 * t156;
t104 = t186 * t155;
t103 = t155 * t171 + t156 * t199;
t102 = t155 * t172 - t156 * t200;
t101 = t155 * t199 - t156 * t171;
t100 = -t155 * t200 - t156 * t172;
t98 = V_base(5) * rSges(2,3) - t146 * V_base(6) + t195;
t97 = t147 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t96 = rSges(5,3) * t155 + t156 * t184;
t95 = -rSges(5,3) * t156 + t155 * t184;
t88 = t146 * V_base(4) - t147 * V_base(5) + t191;
t87 = V_base(5) * rSges(3,3) - t129 * t162 + t177;
t86 = t130 * t162 + (-rSges(3,3) + t208) * V_base(4) + t196;
t85 = t129 * V_base(4) + (-t130 - t210) * V_base(5) + t188;
t84 = rSges(6,1) * t103 + rSges(6,2) * t102 + rSges(6,3) * t201;
t83 = rSges(6,1) * t101 + rSges(6,2) * t100 + rSges(6,3) * t202;
t82 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t201;
t81 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t202;
t80 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t201;
t79 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t202;
t78 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t201;
t77 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t202;
t76 = V_base(5) * rSges(4,3) - t118 * t157 + t175;
t75 = t120 * t157 + (-rSges(4,3) + t197) * V_base(4) + t189;
t74 = t118 * V_base(4) + (-t120 + t187) * V_base(5) + t185;
t73 = t145 * V_base(5) + (-t117 - t95) * t157 + t174;
t72 = t157 * t96 + (-t145 + t197) * V_base(4) + t178;
t71 = t95 * V_base(4) + (t180 - t96) * V_base(5) + t179;
t70 = t109 * t121 - t134 * t83 + t148 * V_base(5) + (-t104 - t117) * t157 + t174;
t69 = t105 * t157 - t109 * t122 + t134 * t84 + (-t148 + t197) * V_base(4) + t178;
t68 = t104 * V_base(4) - t121 * t84 + t122 * t83 + (-t105 + t180) * V_base(5) + t179;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t88 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + t122 * ((t102 * t80 + t103 * t82 + t201 * t78) * t122 + (t102 * t79 + t103 * t81 + t201 * t77) * t121 + (t102 * t107 + t103 * t108 + t106 * t201) * t134) / 0.2e1 + t121 * ((t100 * t80 + t101 * t82 + t202 * t78) * t122 + (t100 * t79 + t101 * t81 + t202 * t77) * t121 + (t100 * t107 + t101 * t108 + t106 * t202) * t134) / 0.2e1 + t134 * ((-t106 * t134 - t77 * t121 - t78 * t122) * t169 + ((-t171 * t80 + t172 * t82) * t122 + (-t171 * t79 + t172 * t81) * t121 + (-t107 * t171 + t108 * t172) * t134) * t167) / 0.2e1 + ((t167 * t93 + t169 * t91 + t111) * V_base(5) + (t167 * t94 + t169 * t92 + t112) * V_base(4) + (t169 * t139 + t167 * t142 + Icges(4,3)) * t157) * t157 / 0.2e1 + (t112 * t157 + t155 * t176 + t156 * t173 + (-t113 * t155 + t115 * t156 - t125 * t158 + t127 * t159 - t140 * t168 + t143 * t170 + Icges(1,4)) * V_base(5) + (-t155 * t114 + t156 * t116 - t158 * t126 + t159 * t128 - t168 * t141 + t170 * t144 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t111 * t157 + t155 * t173 - t176 * t156 + (t156 * t113 + t155 * t115 + t159 * t125 + t158 * t127 + t170 * t140 + t168 * t143 + Icges(1,2)) * V_base(5) + (t114 * t156 + t116 * t155 + t126 * t159 + t128 * t158 + t141 * t170 + t144 * t168 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t168 + Icges(2,6) * t170 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t170 - Icges(2,6) * t168 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,5) * t158 + Icges(3,6) * t159) * V_base(5) + (Icges(3,5) * t159 - Icges(3,6) * t158) * V_base(4) + Icges(3,3) * t162 / 0.2e1) * t162;
T = t1;
