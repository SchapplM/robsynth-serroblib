% Calculate kinetic energy for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:38
% EndTime: 2019-12-31 17:25:40
% DurationCPUTime: 1.29s
% Computational Cost: add. (839->228), mult. (1064->347), div. (0->0), fcn. (952->8), ass. (0->115)
t149 = sin(qJ(2));
t191 = pkin(2) * t149;
t152 = cos(qJ(2));
t190 = pkin(2) * t152;
t150 = sin(qJ(1));
t153 = cos(qJ(1));
t135 = t150 * pkin(1) - t153 * pkin(5);
t83 = -pkin(6) * t153 + t190 * t150;
t188 = -t135 - t83;
t187 = Icges(2,4) * t150;
t186 = Icges(3,4) * t149;
t185 = Icges(3,4) * t152;
t147 = qJ(2) + qJ(3);
t143 = sin(t147);
t184 = Icges(4,4) * t143;
t144 = cos(t147);
t183 = Icges(4,4) * t144;
t182 = t143 * t150;
t181 = t143 * t153;
t148 = sin(qJ(4));
t180 = t148 * t150;
t179 = t148 * t153;
t151 = cos(qJ(4));
t178 = t150 * t151;
t177 = t151 * t153;
t176 = qJD(4) * t143;
t175 = V_base(5) * pkin(4) + V_base(1);
t138 = qJD(2) * t150 + V_base(4);
t140 = V_base(6) + qJD(1);
t137 = -qJD(2) * t153 + V_base(5);
t172 = t137 * t191 + t175;
t116 = qJD(3) * t150 + t138;
t171 = pkin(3) * t144 + pkin(7) * t143;
t170 = rSges(3,1) * t152 - rSges(3,2) * t149;
t169 = rSges(4,1) * t144 - rSges(4,2) * t143;
t168 = Icges(3,1) * t152 - t186;
t167 = Icges(4,1) * t144 - t184;
t166 = -Icges(3,2) * t149 + t185;
t165 = -Icges(4,2) * t143 + t183;
t164 = Icges(3,5) * t152 - Icges(3,6) * t149;
t163 = Icges(4,5) * t144 - Icges(4,6) * t143;
t136 = t153 * pkin(1) + t150 * pkin(5);
t162 = -V_base(4) * pkin(4) + t140 * t136 + V_base(2);
t161 = V_base(4) * t135 - t136 * V_base(5) + V_base(3);
t115 = V_base(5) + (-qJD(2) - qJD(3)) * t153;
t160 = (Icges(4,5) * t143 + Icges(4,6) * t144) * t140 + (-Icges(4,3) * t153 + t163 * t150) * t115 + (Icges(4,3) * t150 + t163 * t153) * t116;
t159 = (Icges(3,5) * t149 + Icges(3,6) * t152) * t140 + (-Icges(3,3) * t153 + t164 * t150) * t137 + (Icges(3,3) * t150 + t164 * t153) * t138;
t84 = pkin(6) * t150 + t190 * t153;
t158 = -t137 * t84 + t138 * t83 + t161;
t157 = -t138 * t191 + t140 * t84 + t162;
t111 = Icges(4,2) * t144 + t184;
t112 = Icges(4,1) * t143 + t183;
t87 = -Icges(4,6) * t153 + t165 * t150;
t88 = Icges(4,6) * t150 + t165 * t153;
t89 = -Icges(4,5) * t153 + t167 * t150;
t90 = Icges(4,5) * t150 + t167 * t153;
t156 = (-t143 * t88 + t144 * t90) * t116 + (-t143 * t87 + t144 * t89) * t115 + (-t111 * t143 + t112 * t144) * t140;
t100 = Icges(3,5) * t150 + t168 * t153;
t126 = Icges(3,2) * t152 + t186;
t129 = Icges(3,1) * t149 + t185;
t97 = -Icges(3,6) * t153 + t166 * t150;
t98 = Icges(3,6) * t150 + t166 * t153;
t99 = -Icges(3,5) * t153 + t168 * t150;
t155 = (t100 * t152 - t149 * t98) * t138 + (-t149 * t97 + t152 * t99) * t137 + (-t126 * t149 + t129 * t152) * t140;
t145 = Icges(2,4) * t153;
t134 = rSges(2,1) * t153 - rSges(2,2) * t150;
t133 = rSges(2,1) * t150 + rSges(2,2) * t153;
t132 = rSges(3,1) * t149 + rSges(3,2) * t152;
t131 = Icges(2,1) * t153 - t187;
t130 = Icges(2,1) * t150 + t145;
t128 = -Icges(2,2) * t150 + t145;
t127 = Icges(2,2) * t153 + t187;
t122 = -qJD(4) * t144 + t140;
t121 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t120 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t119 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t114 = pkin(3) * t143 - pkin(7) * t144;
t113 = rSges(4,1) * t143 + rSges(4,2) * t144;
t108 = t144 * t177 + t180;
t107 = -t144 * t179 + t178;
t106 = t144 * t178 - t179;
t105 = -t144 * t180 - t177;
t104 = t171 * t153;
t103 = t171 * t150;
t102 = rSges(3,3) * t150 + t170 * t153;
t101 = -rSges(3,3) * t153 + t170 * t150;
t94 = t153 * t176 + t116;
t93 = t150 * t176 + t115;
t92 = rSges(4,3) * t150 + t169 * t153;
t91 = -rSges(4,3) * t153 + t169 * t150;
t82 = -rSges(5,3) * t144 + (rSges(5,1) * t151 - rSges(5,2) * t148) * t143;
t81 = -Icges(5,5) * t144 + (Icges(5,1) * t151 - Icges(5,4) * t148) * t143;
t80 = -Icges(5,6) * t144 + (Icges(5,4) * t151 - Icges(5,2) * t148) * t143;
t79 = -Icges(5,3) * t144 + (Icges(5,5) * t151 - Icges(5,6) * t148) * t143;
t78 = V_base(5) * rSges(2,3) - t133 * t140 + t175;
t77 = t134 * t140 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t76 = t133 * V_base(4) - t134 * V_base(5) + V_base(3);
t73 = rSges(5,1) * t108 + rSges(5,2) * t107 + rSges(5,3) * t181;
t72 = rSges(5,1) * t106 + rSges(5,2) * t105 + rSges(5,3) * t182;
t71 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t181;
t70 = Icges(5,1) * t106 + Icges(5,4) * t105 + Icges(5,5) * t182;
t69 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t181;
t68 = Icges(5,4) * t106 + Icges(5,2) * t105 + Icges(5,6) * t182;
t67 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t181;
t66 = Icges(5,5) * t106 + Icges(5,6) * t105 + Icges(5,3) * t182;
t65 = t132 * t137 + (-t101 - t135) * t140 + t175;
t64 = t102 * t140 - t132 * t138 + t162;
t63 = t101 * t138 - t102 * t137 + t161;
t62 = t113 * t115 + (-t91 + t188) * t140 + t172;
t61 = -t113 * t116 + t140 * t92 + t157;
t60 = -t115 * t92 + t116 * t91 + t158;
t59 = t114 * t115 - t122 * t72 + t82 * t93 + (-t103 + t188) * t140 + t172;
t58 = t104 * t140 - t114 * t116 + t122 * t73 - t82 * t94 + t157;
t57 = t103 * t116 - t104 * t115 + t72 * t94 - t73 * t93 + t158;
t1 = m(1) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(2) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t138 * (t159 * t150 + t155 * t153) / 0.2e1 + t137 * (t155 * t150 - t159 * t153) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t116 * (t160 * t150 + t156 * t153) / 0.2e1 + t115 * (t156 * t150 - t160 * t153) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t94 * ((t107 * t69 + t108 * t71 + t67 * t181) * t94 + (t107 * t68 + t108 * t70 + t66 * t181) * t93 + (t107 * t80 + t108 * t81 + t79 * t181) * t122) / 0.2e1 + t93 * ((t105 * t69 + t106 * t71 + t67 * t182) * t94 + (t105 * t68 + t106 * t70 + t66 * t182) * t93 + (t105 * t80 + t106 * t81 + t79 * t182) * t122) / 0.2e1 + t122 * ((-t79 * t122 - t66 * t93 - t67 * t94) * t144 + ((-t148 * t69 + t151 * t71) * t94 + (-t148 * t68 + t151 * t70) * t93 + (-t148 * t80 + t151 * t81) * t122) * t143) / 0.2e1 + ((-t127 * t150 + t130 * t153 + Icges(1,4)) * V_base(5) + (-t150 * t128 + t153 * t131 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t153 * t127 + t150 * t130 + Icges(1,2)) * V_base(5) + (t128 * t153 + t131 * t150 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t149 + t152 * t98) * t138 + (t149 * t99 + t152 * t97) * t137 + (t143 * t90 + t144 * t88) * t116 + (t143 * t89 + t144 * t87) * t115 + (t144 * t111 + t143 * t112 + t152 * t126 + t149 * t129 + Icges(2,3)) * t140) * t140 / 0.2e1 + V_base(4) * t140 * (Icges(2,5) * t153 - Icges(2,6) * t150) + t140 * V_base(5) * (Icges(2,5) * t150 + Icges(2,6) * t153) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
