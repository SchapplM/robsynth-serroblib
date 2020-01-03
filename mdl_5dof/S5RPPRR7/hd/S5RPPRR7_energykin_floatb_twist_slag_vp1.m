% Calculate kinetic energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:28
% EndTime: 2019-12-31 17:59:30
% DurationCPUTime: 1.89s
% Computational Cost: add. (944->240), mult. (978->327), div. (0->0), fcn. (814->8), ass. (0->119)
t220 = Icges(3,4) + Icges(4,6);
t219 = Icges(3,1) + Icges(4,2);
t218 = -Icges(4,4) + Icges(3,5);
t217 = Icges(4,5) - Icges(3,6);
t216 = Icges(3,2) + Icges(4,3);
t159 = qJ(1) + pkin(8);
t153 = cos(t159);
t215 = t220 * t153;
t152 = sin(t159);
t214 = t220 * t152;
t213 = -t216 * t153 - t214;
t212 = t216 * t152 - t215;
t211 = t219 * t152 + t215;
t210 = t219 * t153 - t214;
t129 = qJD(4) * t152 + V_base(5);
t130 = qJD(4) * t153 + V_base(4);
t161 = sin(qJ(4));
t164 = cos(qJ(4));
t197 = Icges(5,4) * t164;
t134 = -Icges(5,2) * t161 + t197;
t198 = Icges(5,4) * t161;
t137 = Icges(5,1) * t164 - t198;
t154 = V_base(6) + qJD(1);
t174 = Icges(5,2) * t164 + t198;
t85 = Icges(5,6) * t153 + t174 * t152;
t86 = Icges(5,6) * t152 - t174 * t153;
t175 = Icges(5,1) * t161 + t197;
t87 = Icges(5,5) * t153 + t175 * t152;
t88 = Icges(5,5) * t152 - t175 * t153;
t207 = (t161 * t87 + t164 * t85) * t130 + (t161 * t88 + t164 * t86) * t129 + (t134 * t164 + t137 * t161) * t154;
t162 = sin(qJ(1));
t205 = pkin(1) * t162;
t165 = cos(qJ(1));
t204 = pkin(1) * t165;
t203 = pkin(6) * t152;
t202 = pkin(6) * t153;
t201 = -pkin(5) - qJ(2);
t200 = Icges(2,4) * t162;
t194 = t152 * t164;
t193 = t153 * t164;
t160 = sin(qJ(5));
t192 = t160 * t161;
t163 = cos(qJ(5));
t191 = t161 * t163;
t190 = qJD(5) * t164;
t189 = t154 * t204 + V_base(2);
t188 = V_base(5) * pkin(5) + V_base(1);
t119 = pkin(2) * t152 - qJ(3) * t153;
t185 = -t119 - t205;
t122 = pkin(2) * t153 + qJ(3) * t152;
t184 = -t122 - t204;
t183 = V_base(5) * qJ(2) + t188;
t182 = V_base(4) * t205 + qJD(2) + V_base(3);
t181 = qJD(3) * t152 + t183;
t180 = pkin(4) * t161 - pkin(7) * t164;
t179 = V_base(4) * t119 + t182;
t178 = rSges(5,1) * t161 + rSges(5,2) * t164;
t173 = Icges(5,5) * t161 + Icges(5,6) * t164;
t171 = V_base(5) * pkin(3) + t181;
t170 = t185 - t203;
t169 = -qJD(3) * t153 + t154 * t122 + t189;
t168 = (Icges(5,3) * t152 - t173 * t153) * t129 + (Icges(5,3) * t153 + t173 * t152) * t130 + (Icges(5,5) * t164 - Icges(5,6) * t161) * t154;
t167 = t154 * t202 + (-pkin(3) + t201) * V_base(4) + t169;
t166 = V_base(4) * t203 + (t184 - t202) * V_base(5) + t179;
t156 = Icges(2,4) * t165;
t144 = pkin(4) * t164 + pkin(7) * t161;
t143 = rSges(2,1) * t165 - t162 * rSges(2,2);
t142 = rSges(5,1) * t164 - rSges(5,2) * t161;
t141 = t162 * rSges(2,1) + rSges(2,2) * t165;
t140 = qJD(5) * t161 + t154;
t139 = Icges(2,1) * t165 - t200;
t138 = Icges(2,1) * t162 + t156;
t136 = -Icges(2,2) * t162 + t156;
t135 = Icges(2,2) * t165 + t200;
t127 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t126 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t125 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t124 = rSges(3,1) * t153 - rSges(3,2) * t152;
t123 = -rSges(4,2) * t153 + rSges(4,3) * t152;
t121 = rSges(3,1) * t152 + rSges(3,2) * t153;
t120 = -rSges(4,2) * t152 - rSges(4,3) * t153;
t105 = t180 * t153;
t104 = t180 * t152;
t103 = rSges(6,3) * t161 + (rSges(6,1) * t163 - rSges(6,2) * t160) * t164;
t101 = Icges(6,5) * t161 + (Icges(6,1) * t163 - Icges(6,4) * t160) * t164;
t100 = Icges(6,6) * t161 + (Icges(6,4) * t163 - Icges(6,2) * t160) * t164;
t99 = Icges(6,3) * t161 + (Icges(6,5) * t163 - Icges(6,6) * t160) * t164;
t98 = t152 * t160 - t153 * t191;
t97 = t152 * t163 + t153 * t192;
t96 = t152 * t191 + t153 * t160;
t95 = -t152 * t192 + t153 * t163;
t94 = -t152 * t190 + t130;
t93 = t153 * t190 + t129;
t92 = rSges(5,3) * t152 - t178 * t153;
t91 = rSges(5,3) * t153 + t178 * t152;
t90 = V_base(5) * rSges(2,3) - t141 * t154 + t188;
t89 = t143 * t154 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t82 = t141 * V_base(4) - t143 * V_base(5) + V_base(3);
t81 = V_base(5) * rSges(3,3) + (-t121 - t205) * t154 + t183;
t80 = t124 * t154 + (-rSges(3,3) + t201) * V_base(4) + t189;
t79 = V_base(4) * t121 + (-t124 - t204) * V_base(5) + t182;
t78 = rSges(6,1) * t98 + rSges(6,2) * t97 + rSges(6,3) * t193;
t77 = rSges(6,1) * t96 + rSges(6,2) * t95 - rSges(6,3) * t194;
t76 = Icges(6,1) * t98 + Icges(6,4) * t97 + Icges(6,5) * t193;
t75 = Icges(6,1) * t96 + Icges(6,4) * t95 - Icges(6,5) * t194;
t74 = Icges(6,4) * t98 + Icges(6,2) * t97 + Icges(6,6) * t193;
t73 = Icges(6,4) * t96 + Icges(6,2) * t95 - Icges(6,6) * t194;
t72 = Icges(6,5) * t98 + Icges(6,6) * t97 + Icges(6,3) * t193;
t71 = Icges(6,5) * t96 + Icges(6,6) * t95 - Icges(6,3) * t194;
t70 = V_base(5) * rSges(4,1) + (-t120 + t185) * t154 + t181;
t69 = t123 * t154 + (-rSges(4,1) + t201) * V_base(4) + t169;
t68 = V_base(4) * t120 + (-t123 + t184) * V_base(5) + t179;
t67 = t129 * t142 + (t170 - t92) * t154 + t171;
t66 = -t130 * t142 + t154 * t91 + t167;
t65 = -t129 * t91 + t130 * t92 + t166;
t64 = t103 * t93 + t129 * t144 - t140 * t78 + (t105 + t170) * t154 + t171;
t63 = -t103 * t94 + t104 * t154 - t130 * t144 + t140 * t77 + t167;
t62 = -t129 * t104 - t130 * t105 - t93 * t77 + t94 * t78 + t166;
t1 = m(1) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(2) * (t82 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t130 * (t207 * t152 + t168 * t153) / 0.2e1 + t129 * (t168 * t152 - t207 * t153) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + t94 * ((-t71 * t194 + t95 * t73 + t96 * t75) * t94 + (-t72 * t194 + t74 * t95 + t76 * t96) * t93 + (t100 * t95 + t101 * t96 - t99 * t194) * t140) / 0.2e1 + t93 * ((t71 * t193 + t73 * t97 + t75 * t98) * t94 + (t72 * t193 + t97 * t74 + t98 * t76) * t93 + (t100 * t97 + t101 * t98 + t99 * t193) * t140) / 0.2e1 + t140 * ((t99 * t140 + t71 * t94 + t72 * t93) * t161 + ((-t160 * t73 + t163 * t75) * t94 + (-t160 * t74 + t163 * t76) * t93 + (-t100 * t160 + t101 * t163) * t140) * t164) / 0.2e1 + ((-t161 * t85 + t164 * t87) * t130 + (-t161 * t86 + t164 * t88) * t129 + (-t161 * t134 + t164 * t137 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t154) * t154 / 0.2e1 + ((-t162 * t135 + t138 * t165 + t213 * t152 + t211 * t153 + Icges(1,4)) * V_base(5) + (-t162 * t136 + t165 * t139 + t212 * t152 + t210 * t153 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t165 * t135 + t162 * t138 + t211 * t152 - t213 * t153 + Icges(1,2)) * V_base(5) + (t136 * t165 + t162 * t139 + t210 * t152 - t212 * t153 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t154 * (Icges(2,5) * t162 + Icges(2,6) * t165 + t218 * t152 - t217 * t153) + V_base(4) * t154 * (Icges(2,5) * t165 - Icges(2,6) * t162 + t217 * t152 + t218 * t153) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
