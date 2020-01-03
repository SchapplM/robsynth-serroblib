% Calculate kinetic energy for
% S4PRRR5
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:33
% EndTime: 2019-12-31 16:33:34
% DurationCPUTime: 1.23s
% Computational Cost: add. (807->228), mult. (1064->342), div. (0->0), fcn. (952->8), ass. (0->115)
t150 = sin(qJ(2));
t192 = pkin(2) * t150;
t152 = cos(qJ(2));
t191 = pkin(2) * t152;
t147 = sin(pkin(7));
t148 = cos(pkin(7));
t131 = t147 * pkin(1) - t148 * pkin(4);
t77 = -pkin(5) * t148 + t191 * t147;
t189 = -t131 - t77;
t188 = Icges(2,4) * t147;
t187 = Icges(3,4) * t150;
t186 = Icges(3,4) * t152;
t146 = qJ(2) + qJ(3);
t143 = sin(t146);
t185 = Icges(4,4) * t143;
t144 = cos(t146);
t184 = Icges(4,4) * t144;
t183 = t143 * t147;
t182 = t143 * t148;
t149 = sin(qJ(4));
t181 = t147 * t149;
t151 = cos(qJ(4));
t180 = t147 * t151;
t179 = t148 * t149;
t178 = t148 * t151;
t177 = qJD(4) * t143;
t176 = V_base(5) * qJ(1) + V_base(1);
t172 = qJD(1) + V_base(3);
t138 = qJD(2) * t147 + V_base(4);
t137 = -qJD(2) * t148 + V_base(5);
t171 = t137 * t192 + t176;
t113 = qJD(3) * t147 + t138;
t170 = pkin(3) * t144 + pkin(6) * t143;
t169 = rSges(3,1) * t152 - rSges(3,2) * t150;
t168 = rSges(4,1) * t144 - rSges(4,2) * t143;
t167 = Icges(3,1) * t152 - t187;
t166 = Icges(4,1) * t144 - t185;
t165 = -Icges(3,2) * t150 + t186;
t164 = -Icges(4,2) * t143 + t184;
t163 = Icges(3,5) * t152 - Icges(3,6) * t150;
t162 = Icges(4,5) * t144 - Icges(4,6) * t143;
t132 = t148 * pkin(1) + t147 * pkin(4);
t161 = -V_base(4) * qJ(1) + V_base(6) * t132 + V_base(2);
t112 = V_base(5) + (-qJD(2) - qJD(3)) * t148;
t160 = V_base(4) * t131 - t132 * V_base(5) + t172;
t159 = (Icges(4,5) * t143 + Icges(4,6) * t144) * V_base(6) + (-Icges(4,3) * t148 + t162 * t147) * t112 + (Icges(4,3) * t147 + t162 * t148) * t113;
t158 = (Icges(3,5) * t150 + Icges(3,6) * t152) * V_base(6) + (-Icges(3,3) * t148 + t163 * t147) * t137 + (Icges(3,3) * t147 + t163 * t148) * t138;
t78 = pkin(5) * t147 + t191 * t148;
t157 = -t138 * t192 + V_base(6) * t78 + t161;
t156 = -t137 * t78 + t138 * t77 + t160;
t110 = Icges(4,2) * t144 + t185;
t111 = Icges(4,1) * t143 + t184;
t87 = -Icges(4,6) * t148 + t164 * t147;
t88 = Icges(4,6) * t147 + t164 * t148;
t89 = -Icges(4,5) * t148 + t166 * t147;
t90 = Icges(4,5) * t147 + t166 * t148;
t155 = (-t143 * t88 + t144 * t90) * t113 + (-t143 * t87 + t144 * t89) * t112 + (-t110 * t143 + t111 * t144) * V_base(6);
t100 = Icges(3,5) * t147 + t167 * t148;
t134 = Icges(3,2) * t152 + t187;
t135 = Icges(3,1) * t150 + t186;
t97 = -Icges(3,6) * t148 + t165 * t147;
t98 = Icges(3,6) * t147 + t165 * t148;
t99 = -Icges(3,5) * t148 + t167 * t147;
t154 = (t100 * t152 - t150 * t98) * t138 + (-t150 * t97 + t152 * t99) * t137 + (-t134 * t150 + t135 * t152) * V_base(6);
t142 = Icges(2,4) * t148;
t136 = rSges(3,1) * t150 + rSges(3,2) * t152;
t130 = -qJD(4) * t144 + V_base(6);
t129 = rSges(2,1) * t148 - rSges(2,2) * t147;
t128 = rSges(2,1) * t147 + rSges(2,2) * t148;
t127 = Icges(2,1) * t148 - t188;
t126 = Icges(2,1) * t147 + t142;
t125 = -Icges(2,2) * t147 + t142;
t124 = Icges(2,2) * t148 + t188;
t121 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t120 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t119 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t117 = pkin(3) * t143 - pkin(6) * t144;
t114 = rSges(4,1) * t143 + rSges(4,2) * t144;
t108 = t144 * t178 + t181;
t107 = -t144 * t179 + t180;
t106 = t144 * t180 - t179;
t105 = -t144 * t181 - t178;
t104 = t170 * t148;
t103 = t170 * t147;
t102 = rSges(3,3) * t147 + t169 * t148;
t101 = -rSges(3,3) * t148 + t169 * t147;
t94 = t148 * t177 + t113;
t93 = t147 * t177 + t112;
t92 = rSges(4,3) * t147 + t168 * t148;
t91 = -rSges(4,3) * t148 + t168 * t147;
t84 = -rSges(5,3) * t144 + (rSges(5,1) * t151 - rSges(5,2) * t149) * t143;
t83 = V_base(5) * rSges(2,3) - t128 * V_base(6) + t176;
t82 = t129 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t81 = -Icges(5,5) * t144 + (Icges(5,1) * t151 - Icges(5,4) * t149) * t143;
t80 = -Icges(5,6) * t144 + (Icges(5,4) * t151 - Icges(5,2) * t149) * t143;
t79 = -Icges(5,3) * t144 + (Icges(5,5) * t151 - Icges(5,6) * t149) * t143;
t75 = t128 * V_base(4) - t129 * V_base(5) + t172;
t73 = rSges(5,1) * t108 + rSges(5,2) * t107 + rSges(5,3) * t182;
t72 = rSges(5,1) * t106 + rSges(5,2) * t105 + rSges(5,3) * t183;
t71 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t182;
t70 = Icges(5,1) * t106 + Icges(5,4) * t105 + Icges(5,5) * t183;
t69 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t182;
t68 = Icges(5,4) * t106 + Icges(5,2) * t105 + Icges(5,6) * t183;
t67 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t182;
t66 = Icges(5,5) * t106 + Icges(5,6) * t105 + Icges(5,3) * t183;
t65 = t136 * t137 + (-t101 - t131) * V_base(6) + t176;
t64 = t102 * V_base(6) - t136 * t138 + t161;
t63 = t101 * t138 - t102 * t137 + t160;
t62 = t112 * t114 + (-t91 + t189) * V_base(6) + t171;
t61 = -t113 * t114 + t92 * V_base(6) + t157;
t60 = -t112 * t92 + t113 * t91 + t156;
t59 = t112 * t117 - t130 * t72 + t84 * t93 + (-t103 + t189) * V_base(6) + t171;
t58 = t104 * V_base(6) - t113 * t117 + t130 * t73 - t84 * t94 + t157;
t57 = t103 * t113 - t104 * t112 + t72 * t94 - t73 * t93 + t156;
t1 = m(1) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(2) * (t75 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t138 * (t158 * t147 + t154 * t148) / 0.2e1 + t137 * (t154 * t147 - t158 * t148) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t113 * (t159 * t147 + t155 * t148) / 0.2e1 + t112 * (t155 * t147 - t159 * t148) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t94 * ((t107 * t69 + t108 * t71 + t67 * t182) * t94 + (t107 * t68 + t108 * t70 + t66 * t182) * t93 + (t107 * t80 + t108 * t81 + t79 * t182) * t130) / 0.2e1 + t93 * ((t105 * t69 + t106 * t71 + t67 * t183) * t94 + (t105 * t68 + t106 * t70 + t66 * t183) * t93 + (t105 * t80 + t106 * t81 + t79 * t183) * t130) / 0.2e1 + t130 * ((-t79 * t130 - t66 * t93 - t67 * t94) * t144 + ((-t149 * t69 + t151 * t71) * t94 + (-t149 * t68 + t151 * t70) * t93 + (-t149 * t80 + t151 * t81) * t130) * t143) / 0.2e1 + ((-t124 * t147 + t126 * t148 + Icges(1,4)) * V_base(5) + (-t147 * t125 + t148 * t127 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t148 * t124 + t147 * t126 + Icges(1,2)) * V_base(5) + (t125 * t148 + t127 * t147 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t150 + t152 * t98) * t138 + (t150 * t99 + t152 * t97) * t137 + (t143 * t90 + t144 * t88) * t113 + (t143 * t89 + t144 * t87) * t112 + (t144 * t110 + t143 * t111 + t152 * t134 + t150 * t135 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t148 - Icges(2,6) * t147 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t147 + Icges(2,6) * t148 + Icges(1,6));
T = t1;
