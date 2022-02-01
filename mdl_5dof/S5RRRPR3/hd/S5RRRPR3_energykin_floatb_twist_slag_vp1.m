% Calculate kinetic energy for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:23
% EndTime: 2022-01-20 11:42:25
% DurationCPUTime: 1.61s
% Computational Cost: add. (1294->242), mult. (946->329), div. (0->0), fcn. (726->10), ass. (0->126)
t230 = Icges(4,3) + Icges(5,3);
t170 = qJ(3) + pkin(9);
t159 = sin(t170);
t160 = cos(t170);
t173 = sin(qJ(3));
t175 = cos(qJ(3));
t229 = Icges(4,5) * t175 + Icges(5,5) * t160 - Icges(4,6) * t173 - Icges(5,6) * t159;
t171 = qJ(1) + qJ(2);
t163 = sin(t171);
t164 = cos(t171);
t212 = Icges(4,4) * t175;
t194 = -Icges(4,2) * t173 + t212;
t103 = -Icges(4,6) * t164 + t163 * t194;
t104 = Icges(4,6) * t163 + t164 * t194;
t213 = Icges(4,4) * t173;
t197 = Icges(4,1) * t175 - t213;
t105 = -Icges(4,5) * t164 + t163 * t197;
t106 = Icges(4,5) * t163 + t164 * t197;
t211 = Icges(5,4) * t159;
t118 = Icges(5,2) * t160 + t211;
t210 = Icges(5,4) * t160;
t119 = Icges(5,1) * t159 + t210;
t136 = -qJD(3) * t164 + V_base(5);
t137 = qJD(3) * t163 + V_base(4);
t141 = Icges(4,2) * t175 + t213;
t144 = Icges(4,1) * t173 + t212;
t162 = V_base(6) + qJD(1);
t157 = qJD(2) + t162;
t193 = -Icges(5,2) * t159 + t210;
t92 = -Icges(5,6) * t164 + t163 * t193;
t93 = Icges(5,6) * t163 + t164 * t193;
t196 = Icges(5,1) * t160 - t211;
t94 = -Icges(5,5) * t164 + t163 * t196;
t95 = Icges(5,5) * t163 + t164 * t196;
t228 = (-t118 * t159 + t119 * t160 - t141 * t173 + t144 * t175) * t157 + (-t104 * t173 + t106 * t175 - t159 * t93 + t160 * t95) * t137 + (-t103 * t173 + t105 * t175 - t159 * t92 + t160 * t94) * t136;
t227 = (Icges(4,5) * t173 + Icges(5,5) * t159 + Icges(4,6) * t175 + Icges(5,6) * t160) * t157 + (t230 * t163 + t229 * t164) * t137 + (t229 * t163 - t230 * t164) * t136;
t226 = -pkin(5) - pkin(6);
t174 = sin(qJ(1));
t222 = pkin(1) * t174;
t176 = cos(qJ(1));
t221 = pkin(1) * t176;
t220 = pkin(3) * t173;
t219 = pkin(4) * t159;
t218 = t175 * pkin(3);
t130 = pkin(2) * t163 - pkin(7) * t164;
t86 = -qJ(4) * t164 + t163 * t218;
t216 = -t130 - t86;
t215 = Icges(2,4) * t174;
t214 = Icges(3,4) * t163;
t161 = qJ(5) + t170;
t155 = sin(t161);
t209 = Icges(6,4) * t155;
t156 = cos(t161);
t208 = Icges(6,4) * t156;
t207 = pkin(4) * t160;
t205 = t162 * t221 + V_base(2);
t204 = V_base(4) * t222 + V_base(3);
t203 = V_base(5) * pkin(5) + V_base(1);
t200 = rSges(4,1) * t175 - rSges(4,2) * t173;
t199 = rSges(5,1) * t160 - rSges(5,2) * t159;
t198 = rSges(6,1) * t156 - rSges(6,2) * t155;
t195 = Icges(6,1) * t156 - t209;
t192 = -Icges(6,2) * t155 + t208;
t189 = Icges(6,5) * t156 - Icges(6,6) * t155;
t188 = V_base(5) * pkin(6) - t162 * t222 + t203;
t113 = V_base(5) + (-qJD(3) - qJD(5)) * t164;
t114 = qJD(5) * t163 + t137;
t187 = (Icges(6,5) * t155 + Icges(6,6) * t156) * t157 + (-Icges(6,3) * t164 + t163 * t189) * t113 + (Icges(6,3) * t163 + t164 * t189) * t114;
t131 = pkin(2) * t164 + pkin(7) * t163;
t184 = t157 * t131 + t226 * V_base(4) + t205;
t183 = qJD(4) * t163 + t136 * t220 + t188;
t182 = V_base(4) * t130 + (-t131 - t221) * V_base(5) + t204;
t181 = t137 * t86 + t182;
t87 = qJ(4) * t163 + t164 * t218;
t180 = -qJD(4) * t164 + t157 * t87 + t184;
t111 = Icges(6,2) * t156 + t209;
t112 = Icges(6,1) * t155 + t208;
t82 = -Icges(6,6) * t164 + t163 * t192;
t83 = Icges(6,6) * t163 + t164 * t192;
t84 = -Icges(6,5) * t164 + t163 * t195;
t85 = Icges(6,5) * t163 + t164 * t195;
t179 = (-t155 * t83 + t156 * t85) * t114 + (-t155 * t82 + t156 * t84) * t113 + (-t111 * t155 + t112 * t156) * t157;
t165 = Icges(2,4) * t176;
t154 = Icges(3,4) * t164;
t149 = rSges(2,1) * t176 - t174 * rSges(2,2);
t148 = t174 * rSges(2,1) + rSges(2,2) * t176;
t147 = rSges(4,1) * t173 + rSges(4,2) * t175;
t146 = Icges(2,1) * t176 - t215;
t145 = Icges(2,1) * t174 + t165;
t143 = -Icges(2,2) * t174 + t165;
t142 = Icges(2,2) * t176 + t215;
t135 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t134 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t133 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t164 - rSges(3,2) * t163;
t128 = rSges(3,1) * t163 + rSges(3,2) * t164;
t127 = Icges(3,1) * t164 - t214;
t126 = Icges(3,1) * t163 + t154;
t125 = -Icges(3,2) * t163 + t154;
t124 = Icges(3,2) * t164 + t214;
t121 = rSges(5,1) * t159 + rSges(5,2) * t160;
t115 = rSges(6,1) * t155 + rSges(6,2) * t156;
t108 = rSges(4,3) * t163 + t164 * t200;
t107 = -rSges(4,3) * t164 + t163 * t200;
t100 = V_base(5) * rSges(2,3) - t148 * t162 + t203;
t99 = t149 * t162 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t98 = t148 * V_base(4) - t149 * V_base(5) + V_base(3);
t97 = rSges(5,3) * t163 + t164 * t199;
t96 = -rSges(5,3) * t164 + t163 * t199;
t89 = rSges(6,3) * t163 + t164 * t198;
t88 = -rSges(6,3) * t164 + t163 * t198;
t77 = V_base(5) * rSges(3,3) - t128 * t157 + t188;
t76 = t129 * t157 + (-rSges(3,3) + t226) * V_base(4) + t205;
t75 = pkin(8) * t163 + t164 * t207;
t74 = -pkin(8) * t164 + t163 * t207;
t73 = V_base(4) * t128 + (-t129 - t221) * V_base(5) + t204;
t72 = t136 * t147 + (-t107 - t130) * t157 + t188;
t71 = t108 * t157 - t137 * t147 + t184;
t70 = t137 * t107 - t136 * t108 + t182;
t69 = t121 * t136 + (-t96 + t216) * t157 + t183;
t68 = t157 * t97 + (-t121 - t220) * t137 + t180;
t67 = t137 * t96 + (-t87 - t97) * t136 + t181;
t66 = t136 * t219 + t113 * t115 + (-t74 - t88 + t216) * t157 + t183;
t65 = -t114 * t115 + (t75 + t89) * t157 + (-t219 - t220) * t137 + t180;
t64 = -t113 * t89 + t114 * t88 + t137 * t74 + (-t75 - t87) * t136 + t181;
t1 = m(1) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t114 * (t163 * t187 + t164 * t179) / 0.2e1 + t113 * (t163 * t179 - t164 * t187) / 0.2e1 + (t228 * t163 - t227 * t164) * t136 / 0.2e1 + (t227 * t163 + t228 * t164) * t137 / 0.2e1 + ((-t124 * t163 + t126 * t164 - t174 * t142 + t145 * t176 + Icges(1,4)) * V_base(5) + (-t125 * t163 + t127 * t164 - t174 * t143 + t146 * t176 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t124 * t164 + t126 * t163 + t142 * t176 + t174 * t145 + Icges(1,2)) * V_base(5) + (t125 * t164 + t127 * t163 + t143 * t176 + t174 * t146 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t155 * t85 + t156 * t83) * t114 + (t155 * t84 + t156 * t82) * t113 + (t104 * t175 + t106 * t173 + t159 * t95 + t160 * t93) * t137 + (t103 * t175 + t105 * t173 + t159 * t94 + t160 * t92) * t136 + (t111 * t156 + t112 * t155 + t118 * t160 + t119 * t159 + t141 * t175 + t144 * t173 + Icges(3,3)) * t157) * t157 / 0.2e1 + V_base(4) * t157 * (Icges(3,5) * t164 - Icges(3,6) * t163) + V_base(5) * t157 * (Icges(3,5) * t163 + Icges(3,6) * t164) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t174 + Icges(2,6) * t176) * V_base(5) + (Icges(2,5) * t176 - Icges(2,6) * t174) * V_base(4) + Icges(2,3) * t162 / 0.2e1) * t162;
T = t1;
