% Calculate kinetic energy for
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:02
% EndTime: 2019-12-05 17:06:03
% DurationCPUTime: 1.44s
% Computational Cost: add. (977->211), mult. (596->274), div. (0->0), fcn. (372->10), ass. (0->108)
t153 = cos(pkin(9));
t190 = pkin(1) * t153;
t146 = V_base(6) + qJD(2);
t189 = pkin(2) * t146;
t140 = qJD(3) + t146;
t188 = pkin(3) * t140;
t187 = -pkin(5) - qJ(1);
t152 = sin(pkin(9));
t186 = Icges(2,4) * t152;
t151 = pkin(9) + qJ(2);
t142 = sin(t151);
t185 = Icges(3,4) * t142;
t145 = qJ(3) + t151;
t138 = sin(t145);
t184 = Icges(4,4) * t138;
t141 = qJ(4) + t145;
t133 = sin(t141);
t183 = Icges(5,4) * t133;
t154 = sin(qJ(5));
t182 = Icges(6,4) * t154;
t155 = cos(qJ(5));
t181 = Icges(6,4) * t155;
t180 = -pkin(6) + t187;
t173 = pkin(1) * V_base(6);
t179 = t153 * t173 + V_base(2);
t178 = V_base(5) * qJ(1) + V_base(1);
t174 = qJD(1) + V_base(3);
t172 = -pkin(7) + t180;
t143 = cos(t151);
t171 = t143 * t189 + t179;
t170 = V_base(4) * t152 * pkin(1) + t174;
t139 = cos(t145);
t169 = t139 * t188 + t171;
t168 = -pkin(2) * t143 - t190;
t167 = V_base(4) * pkin(2) * t142 + t170;
t166 = rSges(6,1) * t155 - rSges(6,2) * t154;
t165 = Icges(6,1) * t155 - t182;
t164 = -Icges(6,2) * t154 + t181;
t163 = Icges(6,5) * t155 - Icges(6,6) * t154;
t162 = V_base(4) * pkin(3) * t138 + t167;
t161 = -pkin(3) * t139 + t168;
t134 = cos(t141);
t110 = -qJD(5) * t134 + V_base(5);
t111 = qJD(5) * t133 + V_base(4);
t132 = qJD(4) + t140;
t160 = (-Icges(6,3) * t134 + t133 * t163) * t110 + (Icges(6,3) * t133 + t134 * t163) * t111 + (Icges(6,5) * t154 + Icges(6,6) * t155) * t132;
t159 = V_base(5) * pkin(5) - t152 * t173 + t178;
t158 = V_base(5) * pkin(6) - t142 * t189 + t159;
t157 = V_base(5) * pkin(7) - t138 * t188 + t158;
t126 = Icges(6,2) * t155 + t182;
t127 = Icges(6,1) * t154 + t181;
t74 = -Icges(6,6) * t134 + t133 * t164;
t75 = Icges(6,6) * t133 + t134 * t164;
t76 = -Icges(6,5) * t134 + t133 * t165;
t77 = Icges(6,5) * t133 + t134 * t165;
t156 = (-t154 * t75 + t155 * t77) * t111 + (-t154 * t74 + t155 * t76) * t110 + (-t126 * t154 + t127 * t155) * t132;
t144 = Icges(2,4) * t153;
t137 = Icges(3,4) * t143;
t131 = Icges(4,4) * t139;
t129 = Icges(5,4) * t134;
t128 = t154 * rSges(6,1) + rSges(6,2) * t155;
t124 = rSges(2,1) * t153 - rSges(2,2) * t152;
t123 = rSges(2,1) * t152 + rSges(2,2) * t153;
t121 = Icges(2,1) * t153 - t186;
t120 = Icges(2,1) * t152 + t144;
t119 = -Icges(2,2) * t152 + t144;
t118 = Icges(2,2) * t153 + t186;
t114 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t113 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t112 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t108 = rSges(3,1) * t143 - rSges(3,2) * t142;
t107 = rSges(3,1) * t142 + rSges(3,2) * t143;
t106 = Icges(3,1) * t143 - t185;
t105 = Icges(3,1) * t142 + t137;
t104 = -Icges(3,2) * t142 + t137;
t103 = Icges(3,2) * t143 + t185;
t100 = rSges(4,1) * t139 - rSges(4,2) * t138;
t99 = rSges(4,1) * t138 + rSges(4,2) * t139;
t98 = Icges(4,1) * t139 - t184;
t97 = Icges(4,1) * t138 + t131;
t96 = -Icges(4,2) * t138 + t131;
t95 = Icges(4,2) * t139 + t184;
t92 = pkin(4) * t134 + pkin(8) * t133;
t91 = pkin(4) * t133 - pkin(8) * t134;
t90 = rSges(5,1) * t134 - rSges(5,2) * t133;
t89 = rSges(5,1) * t133 + rSges(5,2) * t134;
t88 = Icges(5,1) * t134 - t183;
t87 = Icges(5,1) * t133 + t129;
t86 = -Icges(5,2) * t133 + t129;
t85 = Icges(5,2) * t134 + t183;
t82 = V_base(5) * rSges(2,3) - t123 * V_base(6) + t178;
t81 = t124 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t80 = t123 * V_base(4) - t124 * V_base(5) + t174;
t79 = t133 * rSges(6,3) + t134 * t166;
t78 = -t134 * rSges(6,3) + t133 * t166;
t71 = V_base(5) * rSges(3,3) - t107 * t146 + t159;
t70 = t108 * t146 + (-rSges(3,3) + t187) * V_base(4) + t179;
t69 = t107 * V_base(4) + (-t108 - t190) * V_base(5) + t170;
t68 = V_base(5) * rSges(4,3) - t140 * t99 + t158;
t67 = t100 * t140 + (-rSges(4,3) + t180) * V_base(4) + t171;
t66 = t99 * V_base(4) + (-t100 + t168) * V_base(5) + t167;
t65 = V_base(5) * rSges(5,3) - t132 * t89 + t157;
t64 = t132 * t90 + (-rSges(5,3) + t172) * V_base(4) + t169;
t63 = t89 * V_base(4) + (t161 - t90) * V_base(5) + t162;
t62 = t110 * t128 + (-t78 - t91) * t132 + t157;
t61 = -t111 * t128 + (t79 + t92) * t132 + t172 * V_base(4) + t169;
t60 = -t110 * t79 + t111 * t78 + t91 * V_base(4) + (t161 - t92) * V_base(5) + t162;
t1 = m(1) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(2) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t111 * (t160 * t133 + t156 * t134) / 0.2e1 + t110 * (t156 * t133 - t160 * t134) / 0.2e1 + ((t154 * t77 + t155 * t75) * t111 + (t154 * t76 + t155 * t74) * t110 + (t155 * t126 + t154 * t127 + Icges(5,3)) * t132) * t132 / 0.2e1 + ((-t103 * t142 + t105 * t143 - t118 * t152 + t120 * t153 - t133 * t85 + t134 * t87 - t138 * t95 + t139 * t97 + Icges(1,4)) * V_base(5) + (-t142 * t104 + t143 * t106 - t152 * t119 + t153 * t121 - t133 * t86 + t134 * t88 - t138 * t96 + t139 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t143 * t103 + t142 * t105 + t153 * t118 + t152 * t120 + t133 * t87 + t134 * t85 + t138 * t97 + t139 * t95 + Icges(1,2)) * V_base(5) + (t104 * t143 + t106 * t142 + t119 * t153 + t121 * t152 + t133 * t88 + t134 * t86 + t138 * t98 + t139 * t96 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t132 * (Icges(5,5) * t134 - Icges(5,6) * t133) + V_base(5) * t132 * (Icges(5,5) * t133 + Icges(5,6) * t134) + ((Icges(2,5) * t152 + Icges(2,6) * t153 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t153 - Icges(2,6) * t152 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,5) * t142 + Icges(3,6) * t143) * V_base(5) + (Icges(3,5) * t143 - Icges(3,6) * t142) * V_base(4) + Icges(3,3) * t146 / 0.2e1) * t146 + ((Icges(4,5) * t138 + Icges(4,6) * t139) * V_base(5) + (Icges(4,5) * t139 - Icges(4,6) * t138) * V_base(4) + Icges(4,3) * t140 / 0.2e1) * t140;
T = t1;
