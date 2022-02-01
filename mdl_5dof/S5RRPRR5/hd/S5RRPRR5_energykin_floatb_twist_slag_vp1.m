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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:02:08
% EndTime: 2022-01-20 11:02:10
% DurationCPUTime: 1.80s
% Computational Cost: add. (1258->248), mult. (928->343), div. (0->0), fcn. (708->10), ass. (0->132)
t226 = -pkin(5) - pkin(6);
t175 = sin(qJ(1));
t224 = pkin(1) * t175;
t176 = cos(qJ(1));
t223 = pkin(1) * t176;
t172 = sin(pkin(9));
t222 = pkin(3) * t172;
t170 = pkin(9) + qJ(4);
t159 = sin(t170);
t221 = pkin(4) * t159;
t173 = cos(pkin(9));
t220 = t173 * pkin(3);
t171 = qJ(1) + qJ(2);
t163 = sin(t171);
t164 = cos(t171);
t127 = pkin(2) * t163 - qJ(3) * t164;
t80 = -pkin(7) * t164 + t163 * t220;
t219 = -t127 - t80;
t218 = Icges(2,4) * t175;
t217 = Icges(3,4) * t163;
t216 = Icges(4,4) * t172;
t215 = Icges(4,4) * t173;
t214 = Icges(5,4) * t159;
t160 = cos(t170);
t213 = Icges(5,4) * t160;
t161 = qJ(5) + t170;
t155 = sin(t161);
t212 = Icges(6,4) * t155;
t156 = cos(t161);
t211 = Icges(6,4) * t156;
t209 = pkin(4) * t160;
t162 = V_base(6) + qJD(1);
t207 = t162 * t223 + V_base(2);
t206 = V_base(4) * t224 + V_base(3);
t205 = V_base(5) * pkin(5) + V_base(1);
t140 = qJD(4) * t163 + V_base(4);
t129 = pkin(2) * t164 + qJ(3) * t163;
t202 = -t129 - t223;
t201 = t127 * V_base(4) + t206;
t200 = rSges(4,1) * t173 - rSges(4,2) * t172;
t199 = rSges(5,1) * t160 - rSges(5,2) * t159;
t198 = rSges(6,1) * t156 - rSges(6,2) * t155;
t197 = Icges(4,1) * t173 - t216;
t196 = Icges(5,1) * t160 - t214;
t195 = Icges(6,1) * t156 - t212;
t194 = -Icges(4,2) * t172 + t215;
t193 = -Icges(5,2) * t159 + t213;
t192 = -Icges(6,2) * t155 + t211;
t191 = Icges(4,5) * t173 - Icges(4,6) * t172;
t190 = Icges(5,5) * t160 - Icges(5,6) * t159;
t189 = Icges(6,5) * t156 - Icges(6,6) * t155;
t157 = qJD(2) + t162;
t188 = -qJD(3) * t164 + t129 * t157 + t207;
t187 = V_base(5) * pkin(6) - t162 * t224 + t205;
t113 = V_base(5) + (-qJD(4) - qJD(5)) * t164;
t114 = qJD(5) * t163 + t140;
t186 = (Icges(6,5) * t155 + Icges(6,6) * t156) * t157 + (-Icges(6,3) * t164 + t163 * t189) * t113 + (Icges(6,3) * t163 + t164 * t189) * t114;
t139 = -qJD(4) * t164 + V_base(5);
t185 = (Icges(5,5) * t159 + Icges(5,6) * t160) * t157 + (-Icges(5,3) * t164 + t163 * t190) * t139 + (Icges(5,3) * t163 + t164 * t190) * t140;
t184 = qJD(3) * t163 + t187;
t183 = V_base(5) * t222 + t184;
t182 = (-Icges(4,3) * t164 + t163 * t191) * V_base(5) + (Icges(4,3) * t163 + t164 * t191) * V_base(4) + (Icges(4,5) * t172 + Icges(4,6) * t173) * t157;
t81 = pkin(7) * t163 + t164 * t220;
t181 = V_base(4) * t80 + (t202 - t81) * V_base(5) + t201;
t180 = t157 * t81 + (-t222 + t226) * V_base(4) + t188;
t111 = Icges(6,2) * t156 + t212;
t112 = Icges(6,1) * t155 + t211;
t84 = -Icges(6,6) * t164 + t163 * t192;
t85 = Icges(6,6) * t163 + t164 * t192;
t86 = -Icges(6,5) * t164 + t163 * t195;
t87 = Icges(6,5) * t163 + t164 * t195;
t179 = (-t155 * t85 + t156 * t87) * t114 + (-t155 * t84 + t156 * t86) * t113 + (-t111 * t155 + t112 * t156) * t157;
t118 = Icges(5,2) * t160 + t214;
t119 = Icges(5,1) * t159 + t213;
t92 = -Icges(5,6) * t164 + t163 * t193;
t93 = Icges(5,6) * t163 + t164 * t193;
t94 = -Icges(5,5) * t164 + t163 * t196;
t95 = Icges(5,5) * t163 + t164 * t196;
t178 = (-t159 * t93 + t160 * t95) * t140 + (-t159 * t92 + t160 * t94) * t139 + (-t118 * t159 + t119 * t160) * t157;
t103 = -Icges(4,6) * t164 + t163 * t194;
t104 = Icges(4,6) * t163 + t164 * t194;
t105 = -Icges(4,5) * t164 + t163 * t197;
t106 = Icges(4,5) * t163 + t164 * t197;
t136 = Icges(4,2) * t173 + t216;
t137 = Icges(4,1) * t172 + t215;
t177 = (-t104 * t172 + t106 * t173) * V_base(4) + (-t103 * t172 + t105 * t173) * V_base(5) + (-t136 * t172 + t137 * t173) * t157;
t166 = Icges(2,4) * t176;
t154 = Icges(3,4) * t164;
t148 = rSges(2,1) * t176 - rSges(2,2) * t175;
t147 = rSges(2,1) * t175 + rSges(2,2) * t176;
t146 = Icges(2,1) * t176 - t218;
t145 = Icges(2,1) * t175 + t166;
t144 = -Icges(2,2) * t175 + t166;
t143 = Icges(2,2) * t176 + t218;
t138 = rSges(4,1) * t172 + rSges(4,2) * t173;
t134 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t133 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t132 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t130 = rSges(3,1) * t164 - rSges(3,2) * t163;
t128 = rSges(3,1) * t163 + rSges(3,2) * t164;
t126 = Icges(3,1) * t164 - t217;
t125 = Icges(3,1) * t163 + t154;
t124 = -Icges(3,2) * t163 + t154;
t123 = Icges(3,2) * t164 + t217;
t122 = Icges(3,5) * t164 - Icges(3,6) * t163;
t121 = Icges(3,5) * t163 + Icges(3,6) * t164;
t120 = rSges(5,1) * t159 + rSges(5,2) * t160;
t116 = rSges(6,1) * t155 + rSges(6,2) * t156;
t108 = rSges(4,3) * t163 + t164 * t200;
t107 = -rSges(4,3) * t164 + t163 * t200;
t100 = V_base(5) * rSges(2,3) - t147 * t162 + t205;
t99 = t148 * t162 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t98 = t147 * V_base(4) - t148 * V_base(5) + V_base(3);
t97 = rSges(5,3) * t163 + t164 * t199;
t96 = -rSges(5,3) * t164 + t163 * t199;
t89 = rSges(6,3) * t163 + t164 * t198;
t88 = -rSges(6,3) * t164 + t163 * t198;
t77 = V_base(5) * rSges(3,3) - t128 * t157 + t187;
t76 = t130 * t157 + (-rSges(3,3) + t226) * V_base(4) + t207;
t75 = pkin(8) * t163 + t164 * t209;
t74 = -pkin(8) * t164 + t163 * t209;
t73 = V_base(4) * t128 + (-t130 - t223) * V_base(5) + t206;
t72 = t138 * V_base(5) + (-t107 - t127) * t157 + t184;
t71 = t108 * t157 + (-t138 + t226) * V_base(4) + t188;
t70 = V_base(4) * t107 + (-t108 + t202) * V_base(5) + t201;
t69 = t120 * t139 + (-t96 + t219) * t157 + t183;
t68 = -t120 * t140 + t157 * t97 + t180;
t67 = -t139 * t97 + t140 * t96 + t181;
t66 = t139 * t221 + t113 * t116 + (-t74 - t88 + t219) * t157 + t183;
t65 = -t140 * t221 - t114 * t116 + (t75 + t89) * t157 + t180;
t64 = -t113 * t89 + t114 * t88 - t139 * t75 + t140 * t74 + t181;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t140 * (t185 * t163 + t178 * t164) / 0.2e1 + t139 * (t178 * t163 - t185 * t164) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t114 * (t186 * t163 + t179 * t164) / 0.2e1 + t113 * (t179 * t163 - t186 * t164) / 0.2e1 + ((t159 * t95 + t160 * t93) * t140 + (t159 * t94 + t160 * t92) * t139 + (t155 * t87 + t156 * t85) * t114 + (t155 * t86 + t156 * t84) * t113 + (t103 * t173 + t105 * t172 + t121) * V_base(5) + (t104 * t173 + t106 * t172 + t122) * V_base(4) + (t111 * t156 + t112 * t155 + t118 * t160 + t119 * t159 + t136 * t173 + t137 * t172 + Icges(3,3)) * t157) * t157 / 0.2e1 + (t122 * t157 + t182 * t163 + t177 * t164 + (-t123 * t163 + t125 * t164 - t143 * t175 + t145 * t176 + Icges(1,4)) * V_base(5) + (-t124 * t163 + t126 * t164 - t175 * t144 + t146 * t176 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t121 * t157 + t177 * t163 - t182 * t164 + (t123 * t164 + t125 * t163 + t143 * t176 + t175 * t145 + Icges(1,2)) * V_base(5) + (t124 * t164 + t126 * t163 + t144 * t176 + t146 * t175 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t175 + Icges(2,6) * t176) * V_base(5) + (Icges(2,5) * t176 - Icges(2,6) * t175) * V_base(4) + Icges(2,3) * t162 / 0.2e1) * t162;
T = t1;
