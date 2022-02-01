% Calculate kinetic energy for
% S5RRRPR2
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:30
% EndTime: 2022-01-20 11:30:31
% DurationCPUTime: 1.31s
% Computational Cost: add. (1010->212), mult. (596->269), div. (0->0), fcn. (372->10), ass. (0->109)
t191 = -pkin(5) - pkin(6);
t153 = sin(qJ(1));
t189 = pkin(1) * t153;
t155 = cos(qJ(1));
t188 = pkin(1) * t155;
t151 = qJ(1) + qJ(2);
t143 = sin(t151);
t187 = pkin(2) * t143;
t144 = cos(t151);
t186 = pkin(2) * t144;
t147 = qJ(3) + t151;
t140 = sin(t147);
t185 = pkin(3) * t140;
t141 = cos(t147);
t184 = pkin(3) * t141;
t183 = Icges(2,4) * t153;
t182 = Icges(3,4) * t143;
t181 = Icges(4,4) * t140;
t139 = pkin(9) + t147;
t134 = sin(t139);
t180 = Icges(5,4) * t134;
t152 = sin(qJ(5));
t179 = Icges(6,4) * t152;
t154 = cos(qJ(5));
t178 = Icges(6,4) * t154;
t177 = -pkin(7) + t191;
t142 = V_base(6) + qJD(1);
t176 = t142 * t188 + V_base(2);
t175 = V_base(4) * t189 + V_base(3);
t174 = V_base(5) * pkin(5) + V_base(1);
t171 = -qJ(4) + t177;
t138 = qJD(2) + t142;
t170 = t138 * t186 + t176;
t169 = V_base(4) * t187 + t175;
t133 = qJD(3) + t138;
t168 = t133 * t184 + t170;
t167 = -t186 - t188;
t166 = rSges(6,1) * t154 - rSges(6,2) * t152;
t165 = Icges(6,1) * t154 - t179;
t164 = -Icges(6,2) * t152 + t178;
t163 = Icges(6,5) * t154 - Icges(6,6) * t152;
t162 = V_base(4) * t185 + qJD(4) + t169;
t161 = t167 - t184;
t160 = V_base(5) * pkin(6) - t142 * t189 + t174;
t135 = cos(t139);
t110 = -qJD(5) * t135 + V_base(5);
t111 = qJD(5) * t134 + V_base(4);
t159 = (-Icges(6,3) * t135 + t134 * t163) * t110 + (Icges(6,3) * t134 + t135 * t163) * t111 + (Icges(6,5) * t152 + Icges(6,6) * t154) * t133;
t158 = V_base(5) * pkin(7) - t138 * t187 + t160;
t157 = V_base(5) * qJ(4) + t158;
t120 = Icges(6,2) * t154 + t179;
t123 = Icges(6,1) * t152 + t178;
t74 = -Icges(6,6) * t135 + t134 * t164;
t75 = Icges(6,6) * t134 + t135 * t164;
t76 = -Icges(6,5) * t135 + t134 * t165;
t77 = Icges(6,5) * t134 + t135 * t165;
t156 = (-t152 * t75 + t154 * t77) * t111 + (-t152 * t74 + t154 * t76) * t110 + (-t120 * t152 + t123 * t154) * t133;
t146 = Icges(2,4) * t155;
t137 = Icges(3,4) * t144;
t132 = Icges(4,4) * t141;
t129 = Icges(5,4) * t135;
t128 = rSges(2,1) * t155 - t153 * rSges(2,2);
t127 = t153 * rSges(2,1) + rSges(2,2) * t155;
t126 = rSges(6,1) * t152 + rSges(6,2) * t154;
t125 = Icges(2,1) * t155 - t183;
t124 = Icges(2,1) * t153 + t146;
t122 = -Icges(2,2) * t153 + t146;
t121 = Icges(2,2) * t155 + t183;
t115 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t114 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t113 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t109 = rSges(3,1) * t144 - rSges(3,2) * t143;
t108 = rSges(3,1) * t143 + rSges(3,2) * t144;
t107 = Icges(3,1) * t144 - t182;
t106 = Icges(3,1) * t143 + t137;
t105 = -Icges(3,2) * t143 + t137;
t104 = Icges(3,2) * t144 + t182;
t100 = rSges(4,1) * t141 - rSges(4,2) * t140;
t99 = rSges(4,1) * t140 + rSges(4,2) * t141;
t98 = Icges(4,1) * t141 - t181;
t97 = Icges(4,1) * t140 + t132;
t96 = -Icges(4,2) * t140 + t132;
t95 = Icges(4,2) * t141 + t181;
t92 = pkin(4) * t135 + pkin(8) * t134;
t91 = pkin(4) * t134 - pkin(8) * t135;
t90 = rSges(5,1) * t135 - rSges(5,2) * t134;
t89 = rSges(5,1) * t134 + rSges(5,2) * t135;
t88 = Icges(5,1) * t135 - t180;
t87 = Icges(5,1) * t134 + t129;
t86 = -Icges(5,2) * t134 + t129;
t85 = Icges(5,2) * t135 + t180;
t82 = V_base(5) * rSges(2,3) - t127 * t142 + t174;
t81 = t128 * t142 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t80 = t127 * V_base(4) - t128 * V_base(5) + V_base(3);
t79 = rSges(6,3) * t134 + t135 * t166;
t78 = -rSges(6,3) * t135 + t134 * t166;
t71 = V_base(5) * rSges(3,3) - t108 * t138 + t160;
t70 = t109 * t138 + (-rSges(3,3) + t191) * V_base(4) + t176;
t69 = V_base(4) * t108 + (-t109 - t188) * V_base(5) + t175;
t68 = V_base(5) * rSges(4,3) - t133 * t99 + t158;
t67 = t100 * t133 + (-rSges(4,3) + t177) * V_base(4) + t170;
t66 = V_base(4) * t99 + (-t100 + t167) * V_base(5) + t169;
t65 = V_base(5) * rSges(5,3) + (-t89 - t185) * t133 + t157;
t64 = t133 * t90 + (-rSges(5,3) + t171) * V_base(4) + t168;
t63 = V_base(4) * t89 + (t161 - t90) * V_base(5) + t162;
t62 = t110 * t126 + (-t78 - t91 - t185) * t133 + t157;
t61 = -t111 * t126 + (t79 + t92) * t133 + t171 * V_base(4) + t168;
t60 = -t110 * t79 + t111 * t78 + V_base(4) * t91 + (t161 - t92) * V_base(5) + t162;
t1 = m(1) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(2) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t111 * (t159 * t134 + t156 * t135) / 0.2e1 + t110 * (t156 * t134 - t159 * t135) / 0.2e1 + ((t152 * t77 + t154 * t75) * t111 + (t152 * t76 + t154 * t74) * t110 + (t154 * t120 + t152 * t123 + Icges(4,3) + Icges(5,3)) * t133) * t133 / 0.2e1 + ((-t104 * t143 + t106 * t144 - t153 * t121 + t124 * t155 - t134 * t85 + t135 * t87 - t140 * t95 + t141 * t97 + Icges(1,4)) * V_base(5) + (-t143 * t105 + t144 * t107 - t153 * t122 + t155 * t125 - t134 * t86 + t135 * t88 - t140 * t96 + t141 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t144 * t104 + t143 * t106 + t155 * t121 + t153 * t124 + t134 * t87 + t135 * t85 + t140 * t97 + t141 * t95 + Icges(1,2)) * V_base(5) + (t105 * t144 + t107 * t143 + t122 * t155 + t153 * t125 + t134 * t88 + t135 * t86 + t140 * t98 + t141 * t96 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t133 * (Icges(4,5) * t140 + Icges(5,5) * t134 + Icges(4,6) * t141 + Icges(5,6) * t135) + V_base(4) * t133 * (Icges(4,5) * t141 + Icges(5,5) * t135 - Icges(4,6) * t140 - Icges(5,6) * t134) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t153 + Icges(2,6) * t155) * V_base(5) + (Icges(2,5) * t155 - Icges(2,6) * t153) * V_base(4) + Icges(2,3) * t142 / 0.2e1) * t142 + ((Icges(3,5) * t143 + Icges(3,6) * t144) * V_base(5) + (Icges(3,5) * t144 - Icges(3,6) * t143) * V_base(4) + Icges(3,3) * t138 / 0.2e1) * t138;
T = t1;
