% Calculate kinetic energy for
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:35
% EndTime: 2019-12-05 15:42:37
% DurationCPUTime: 1.55s
% Computational Cost: add. (1213->248), mult. (928->340), div. (0->0), fcn. (708->10), ass. (0->132)
t174 = cos(pkin(8));
t226 = pkin(1) * t174;
t171 = sin(pkin(9));
t225 = pkin(3) * t171;
t169 = pkin(9) + qJ(4);
t158 = sin(t169);
t224 = pkin(4) * t158;
t173 = cos(pkin(9));
t223 = t173 * pkin(3);
t222 = -pkin(5) - qJ(1);
t170 = pkin(8) + qJ(2);
t159 = sin(t170);
t161 = cos(t170);
t127 = t159 * pkin(2) - t161 * qJ(3);
t80 = -pkin(6) * t161 + t159 * t223;
t221 = -t127 - t80;
t172 = sin(pkin(8));
t220 = Icges(2,4) * t172;
t219 = Icges(3,4) * t159;
t218 = Icges(4,4) * t171;
t217 = Icges(4,4) * t173;
t216 = Icges(5,4) * t158;
t160 = cos(t169);
t215 = Icges(5,4) * t160;
t163 = qJ(5) + t169;
t155 = sin(t163);
t214 = Icges(6,4) * t155;
t156 = cos(t163);
t213 = Icges(6,4) * t156;
t211 = pkin(4) * t160;
t203 = pkin(1) * V_base(6);
t209 = t174 * t203 + V_base(2);
t208 = V_base(5) * qJ(1) + V_base(1);
t204 = qJD(1) + V_base(3);
t145 = qJD(4) * t159 + V_base(4);
t129 = t161 * pkin(2) + t159 * qJ(3);
t202 = -t129 - t226;
t201 = V_base(4) * t172 * pkin(1) + t204;
t200 = V_base(4) * t127 + t201;
t199 = rSges(4,1) * t173 - rSges(4,2) * t171;
t198 = rSges(5,1) * t160 - rSges(5,2) * t158;
t197 = rSges(6,1) * t156 - rSges(6,2) * t155;
t196 = Icges(4,1) * t173 - t218;
t195 = Icges(5,1) * t160 - t216;
t194 = Icges(6,1) * t156 - t214;
t193 = -Icges(4,2) * t171 + t217;
t192 = -Icges(5,2) * t158 + t215;
t191 = -Icges(6,2) * t155 + t213;
t190 = Icges(4,5) * t173 - Icges(4,6) * t171;
t189 = Icges(5,5) * t160 - Icges(5,6) * t158;
t188 = Icges(6,5) * t156 - Icges(6,6) * t155;
t164 = V_base(6) + qJD(2);
t187 = -qJD(3) * t161 + t164 * t129 + t209;
t110 = V_base(5) + (-qJD(4) - qJD(5)) * t161;
t111 = qJD(5) * t159 + t145;
t186 = t110 * (-Icges(6,3) * t161 + t159 * t188) + t111 * (Icges(6,3) * t159 + t161 * t188) + (Icges(6,5) * t155 + Icges(6,6) * t156) * t164;
t144 = -qJD(4) * t161 + V_base(5);
t185 = (Icges(5,5) * t158 + Icges(5,6) * t160) * t164 + t144 * (-Icges(5,3) * t161 + t159 * t189) + t145 * (Icges(5,3) * t159 + t161 * t189);
t184 = V_base(5) * pkin(5) - t172 * t203 + t208;
t183 = qJD(3) * t159 + t184;
t182 = (Icges(4,3) * t159 + t161 * t190) * V_base(4) + (Icges(4,5) * t171 + Icges(4,6) * t173) * t164 + (-Icges(4,3) * t161 + t159 * t190) * V_base(5);
t181 = V_base(5) * t225 + t183;
t81 = pkin(6) * t159 + t161 * t223;
t180 = V_base(4) * t80 + (t202 - t81) * V_base(5) + t200;
t179 = t164 * t81 + (t222 - t225) * V_base(4) + t187;
t114 = Icges(6,2) * t156 + t214;
t115 = Icges(6,1) * t155 + t213;
t84 = -Icges(6,6) * t161 + t159 * t191;
t85 = Icges(6,6) * t159 + t161 * t191;
t86 = -Icges(6,5) * t161 + t159 * t194;
t87 = Icges(6,5) * t159 + t161 * t194;
t178 = (-t155 * t85 + t156 * t87) * t111 + (-t155 * t84 + t156 * t86) * t110 + (-t114 * t155 + t115 * t156) * t164;
t120 = Icges(5,2) * t160 + t216;
t123 = Icges(5,1) * t158 + t215;
t93 = -Icges(5,6) * t161 + t159 * t192;
t94 = Icges(5,6) * t159 + t161 * t192;
t95 = -Icges(5,5) * t161 + t159 * t195;
t96 = Icges(5,5) * t159 + t161 * t195;
t177 = (-t158 * t94 + t160 * t96) * t145 + (-t158 * t93 + t160 * t95) * t144 + (-t120 * t158 + t123 * t160) * t164;
t101 = -Icges(4,6) * t161 + t159 * t193;
t102 = Icges(4,6) * t159 + t161 * t193;
t103 = -Icges(4,5) * t161 + t159 * t196;
t104 = Icges(4,5) * t159 + t161 * t196;
t138 = Icges(4,2) * t173 + t218;
t141 = Icges(4,1) * t171 + t217;
t176 = (-t102 * t171 + t104 * t173) * V_base(4) + (-t101 * t171 + t103 * t173) * V_base(5) + (-t138 * t171 + t141 * t173) * t164;
t162 = Icges(2,4) * t174;
t154 = Icges(3,4) * t161;
t148 = rSges(2,1) * t174 - rSges(2,2) * t172;
t147 = rSges(2,1) * t172 + rSges(2,2) * t174;
t146 = rSges(4,1) * t171 + rSges(4,2) * t173;
t143 = Icges(2,1) * t174 - t220;
t142 = Icges(2,1) * t172 + t162;
t140 = -Icges(2,2) * t172 + t162;
t139 = Icges(2,2) * t174 + t220;
t134 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t133 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t132 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = rSges(3,1) * t161 - rSges(3,2) * t159;
t128 = rSges(3,1) * t159 + rSges(3,2) * t161;
t126 = rSges(5,1) * t158 + rSges(5,2) * t160;
t125 = Icges(3,1) * t161 - t219;
t124 = Icges(3,1) * t159 + t154;
t122 = -Icges(3,2) * t159 + t154;
t121 = Icges(3,2) * t161 + t219;
t119 = Icges(3,5) * t161 - Icges(3,6) * t159;
t118 = Icges(3,5) * t159 + Icges(3,6) * t161;
t116 = rSges(6,1) * t155 + rSges(6,2) * t156;
t108 = V_base(5) * rSges(2,3) - t147 * V_base(6) + t208;
t107 = t148 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t106 = rSges(4,3) * t159 + t161 * t199;
t105 = -rSges(4,3) * t161 + t159 * t199;
t98 = rSges(5,3) * t159 + t161 * t198;
t97 = -rSges(5,3) * t161 + t159 * t198;
t90 = t147 * V_base(4) - t148 * V_base(5) + t204;
t89 = rSges(6,3) * t159 + t161 * t197;
t88 = -rSges(6,3) * t161 + t159 * t197;
t77 = V_base(5) * rSges(3,3) - t128 * t164 + t184;
t76 = t130 * t164 + (-rSges(3,3) + t222) * V_base(4) + t209;
t75 = pkin(7) * t159 + t161 * t211;
t74 = -pkin(7) * t161 + t159 * t211;
t73 = t128 * V_base(4) + (-t130 - t226) * V_base(5) + t201;
t72 = t146 * V_base(5) + (-t105 - t127) * t164 + t183;
t71 = t106 * t164 + (-t146 + t222) * V_base(4) + t187;
t70 = t105 * V_base(4) + (-t106 + t202) * V_base(5) + t200;
t69 = t126 * t144 + (-t97 + t221) * t164 + t181;
t68 = -t126 * t145 + t164 * t98 + t179;
t67 = -t144 * t98 + t145 * t97 + t180;
t66 = t144 * t224 + t110 * t116 + (-t74 - t88 + t221) * t164 + t181;
t65 = -t145 * t224 - t111 * t116 + (t75 + t89) * t164 + t179;
t64 = -t110 * t89 + t111 * t88 - t144 * t75 + t145 * t74 + t180;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t107 ^ 2 + t108 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t145 * (t185 * t159 + t177 * t161) / 0.2e1 + t144 * (t177 * t159 - t185 * t161) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t111 * (t186 * t159 + t178 * t161) / 0.2e1 + t110 * (t178 * t159 - t186 * t161) / 0.2e1 + ((t158 * t96 + t160 * t94) * t145 + (t158 * t95 + t160 * t93) * t144 + (t155 * t87 + t156 * t85) * t111 + (t155 * t86 + t156 * t84) * t110 + (t101 * t173 + t103 * t171 + t118) * V_base(5) + (t102 * t173 + t104 * t171 + t119) * V_base(4) + (t114 * t156 + t115 * t155 + t120 * t160 + t123 * t158 + t138 * t173 + t141 * t171 + Icges(3,3)) * t164) * t164 / 0.2e1 + (t119 * t164 + t182 * t159 + t176 * t161 + (-t121 * t159 + t124 * t161 - t139 * t172 + t142 * t174 + Icges(1,4)) * V_base(5) + (-t122 * t159 + t125 * t161 - t140 * t172 + t143 * t174 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t118 * t164 + t176 * t159 - t182 * t161 + (t121 * t161 + t124 * t159 + t139 * t174 + t142 * t172 + Icges(1,2)) * V_base(5) + (t122 * t161 + t125 * t159 + t140 * t174 + t143 * t172 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t172 + Icges(2,6) * t174 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t174 - Icges(2,6) * t172 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
