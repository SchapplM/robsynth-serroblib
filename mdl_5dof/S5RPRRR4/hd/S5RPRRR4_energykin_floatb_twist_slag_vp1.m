% Calculate kinetic energy for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:25
% EndTime: 2022-01-23 09:34:27
% DurationCPUTime: 1.39s
% Computational Cost: add. (990->212), mult. (596->271), div. (0->0), fcn. (372->10), ass. (0->108)
t153 = sin(qJ(1));
t189 = pkin(1) * t153;
t155 = cos(qJ(1));
t188 = pkin(1) * t155;
t151 = qJ(1) + pkin(9);
t142 = sin(t151);
t187 = pkin(2) * t142;
t143 = cos(t151);
t186 = pkin(2) * t143;
t145 = V_base(6) + qJD(1);
t140 = qJD(3) + t145;
t185 = pkin(3) * t140;
t184 = -pkin(5) - qJ(2);
t183 = Icges(2,4) * t153;
t182 = Icges(3,4) * t142;
t144 = qJ(3) + t151;
t138 = sin(t144);
t181 = Icges(4,4) * t138;
t141 = qJ(4) + t144;
t134 = sin(t141);
t180 = Icges(5,4) * t134;
t152 = sin(qJ(5));
t179 = Icges(6,4) * t152;
t154 = cos(qJ(5));
t178 = Icges(6,4) * t154;
t177 = -pkin(6) + t184;
t176 = t145 * t188 + V_base(2);
t175 = V_base(5) * pkin(5) + V_base(1);
t172 = -pkin(7) + t177;
t171 = t145 * t186 + t176;
t170 = V_base(5) * qJ(2) + t175;
t169 = V_base(4) * t189 + qJD(2) + V_base(3);
t139 = cos(t144);
t168 = t139 * t185 + t171;
t167 = -t186 - t188;
t166 = V_base(4) * t187 + t169;
t165 = rSges(6,1) * t154 - rSges(6,2) * t152;
t164 = Icges(6,1) * t154 - t179;
t163 = -Icges(6,2) * t152 + t178;
t162 = Icges(6,5) * t154 - Icges(6,6) * t152;
t161 = V_base(4) * pkin(3) * t138 + t166;
t160 = -pkin(3) * t139 + t167;
t135 = cos(t141);
t110 = -qJD(5) * t135 + V_base(5);
t111 = qJD(5) * t134 + V_base(4);
t133 = qJD(4) + t140;
t159 = t110 * (-Icges(6,3) * t135 + t134 * t162) + t111 * (Icges(6,3) * t134 + t135 * t162) + (Icges(6,5) * t152 + Icges(6,6) * t154) * t133;
t158 = V_base(5) * pkin(6) + (-t187 - t189) * t145 + t170;
t157 = V_base(5) * pkin(7) - t138 * t185 + t158;
t120 = Icges(6,2) * t154 + t179;
t123 = Icges(6,1) * t152 + t178;
t74 = -Icges(6,6) * t135 + t134 * t163;
t75 = Icges(6,6) * t134 + t135 * t163;
t76 = -Icges(6,5) * t135 + t134 * t164;
t77 = Icges(6,5) * t134 + t135 * t164;
t156 = (-t152 * t75 + t154 * t77) * t111 + (-t152 * t74 + t154 * t76) * t110 + (-t120 * t152 + t123 * t154) * t133;
t147 = Icges(2,4) * t155;
t137 = Icges(3,4) * t143;
t132 = Icges(4,4) * t139;
t129 = Icges(5,4) * t135;
t128 = rSges(2,1) * t155 - t153 * rSges(2,2);
t127 = t153 * rSges(2,1) + rSges(2,2) * t155;
t126 = rSges(6,1) * t152 + rSges(6,2) * t154;
t125 = Icges(2,1) * t155 - t183;
t124 = Icges(2,1) * t153 + t147;
t122 = -Icges(2,2) * t153 + t147;
t121 = Icges(2,2) * t155 + t183;
t114 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t113 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t112 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t108 = rSges(3,1) * t143 - rSges(3,2) * t142;
t107 = rSges(3,1) * t142 + rSges(3,2) * t143;
t106 = Icges(3,1) * t143 - t182;
t105 = Icges(3,1) * t142 + t137;
t104 = -Icges(3,2) * t142 + t137;
t103 = Icges(3,2) * t143 + t182;
t100 = rSges(4,1) * t139 - rSges(4,2) * t138;
t99 = rSges(4,1) * t138 + rSges(4,2) * t139;
t98 = Icges(4,1) * t139 - t181;
t97 = Icges(4,1) * t138 + t132;
t96 = -Icges(4,2) * t138 + t132;
t95 = Icges(4,2) * t139 + t181;
t92 = pkin(4) * t135 + pkin(8) * t134;
t91 = pkin(4) * t134 - pkin(8) * t135;
t90 = rSges(5,1) * t135 - rSges(5,2) * t134;
t89 = rSges(5,1) * t134 + rSges(5,2) * t135;
t88 = Icges(5,1) * t135 - t180;
t87 = Icges(5,1) * t134 + t129;
t86 = -Icges(5,2) * t134 + t129;
t85 = Icges(5,2) * t135 + t180;
t82 = V_base(5) * rSges(2,3) - t127 * t145 + t175;
t81 = t128 * t145 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t80 = t127 * V_base(4) - t128 * V_base(5) + V_base(3);
t79 = rSges(6,3) * t134 + t135 * t165;
t78 = -rSges(6,3) * t135 + t134 * t165;
t71 = V_base(5) * rSges(3,3) + (-t107 - t189) * t145 + t170;
t70 = t108 * t145 + (-rSges(3,3) + t184) * V_base(4) + t176;
t69 = V_base(4) * t107 + (-t108 - t188) * V_base(5) + t169;
t68 = V_base(5) * rSges(4,3) - t140 * t99 + t158;
t67 = t100 * t140 + (-rSges(4,3) + t177) * V_base(4) + t171;
t66 = V_base(4) * t99 + (-t100 + t167) * V_base(5) + t166;
t65 = V_base(5) * rSges(5,3) - t133 * t89 + t157;
t64 = t133 * t90 + (-rSges(5,3) + t172) * V_base(4) + t168;
t63 = V_base(4) * t89 + (t160 - t90) * V_base(5) + t161;
t62 = t110 * t126 + (-t78 - t91) * t133 + t157;
t61 = -t111 * t126 + (t79 + t92) * t133 + t172 * V_base(4) + t168;
t60 = -t110 * t79 + t111 * t78 + V_base(4) * t91 + (t160 - t92) * V_base(5) + t161;
t1 = m(1) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(2) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t111 * (t159 * t134 + t156 * t135) / 0.2e1 + t110 * (t156 * t134 - t159 * t135) / 0.2e1 + ((t152 * t77 + t154 * t75) * t111 + (t152 * t76 + t154 * t74) * t110 + (t154 * t120 + t152 * t123 + Icges(5,3)) * t133) * t133 / 0.2e1 + ((-t103 * t142 + t105 * t143 - t153 * t121 + t124 * t155 - t134 * t85 + t135 * t87 - t138 * t95 + t139 * t97 + Icges(1,4)) * V_base(5) + (-t142 * t104 + t143 * t106 - t153 * t122 + t155 * t125 - t134 * t86 + t135 * t88 - t138 * t96 + t139 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t143 * t103 + t142 * t105 + t155 * t121 + t153 * t124 + t134 * t87 + t135 * t85 + t138 * t97 + t139 * t95 + Icges(1,2)) * V_base(5) + (t104 * t143 + t106 * t142 + t122 * t155 + t153 * t125 + t134 * t88 + t135 * t86 + t138 * t98 + t139 * t96 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t133 * (Icges(5,5) * t135 - Icges(5,6) * t134) + V_base(5) * t133 * (Icges(5,5) * t134 + Icges(5,6) * t135) + ((Icges(2,5) * t155 + Icges(3,5) * t143 - Icges(2,6) * t153 - Icges(3,6) * t142) * V_base(4) + (Icges(2,5) * t153 + Icges(3,5) * t142 + Icges(2,6) * t155 + Icges(3,6) * t143) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t145) * t145 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(4,5) * t138 + Icges(4,6) * t139) * V_base(5) + (Icges(4,5) * t139 - Icges(4,6) * t138) * V_base(4) + Icges(4,3) * t140 / 0.2e1) * t140;
T = t1;
