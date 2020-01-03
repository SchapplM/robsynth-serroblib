% Calculate kinetic energy for
% S4RRPR7
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:54
% EndTime: 2019-12-31 17:05:56
% DurationCPUTime: 1.49s
% Computational Cost: add. (815->225), mult. (1040->325), div. (0->0), fcn. (928->8), ass. (0->111)
t196 = Icges(3,3) + Icges(4,3);
t145 = qJ(2) + pkin(7);
t138 = sin(t145);
t139 = cos(t145);
t148 = sin(qJ(2));
t151 = cos(qJ(2));
t195 = Icges(3,5) * t151 + Icges(4,5) * t139 - Icges(3,6) * t148 - Icges(4,6) * t138;
t182 = Icges(4,4) * t138;
t111 = Icges(4,2) * t139 + t182;
t181 = Icges(4,4) * t139;
t112 = Icges(4,1) * t138 + t181;
t184 = Icges(3,4) * t148;
t124 = Icges(3,2) * t151 + t184;
t183 = Icges(3,4) * t151;
t127 = Icges(3,1) * t148 + t183;
t152 = cos(qJ(1));
t135 = -qJD(2) * t152 + V_base(5);
t149 = sin(qJ(1));
t136 = qJD(2) * t149 + V_base(4);
t140 = V_base(6) + qJD(1);
t163 = -Icges(4,2) * t138 + t181;
t87 = -Icges(4,6) * t152 + t163 * t149;
t88 = Icges(4,6) * t149 + t163 * t152;
t165 = Icges(4,1) * t139 - t182;
t89 = -Icges(4,5) * t152 + t165 * t149;
t90 = Icges(4,5) * t149 + t165 * t152;
t164 = -Icges(3,2) * t148 + t183;
t95 = -Icges(3,6) * t152 + t164 * t149;
t96 = Icges(3,6) * t149 + t164 * t152;
t166 = Icges(3,1) * t151 - t184;
t97 = -Icges(3,5) * t152 + t166 * t149;
t98 = Icges(3,5) * t149 + t166 * t152;
t194 = (-t111 * t138 + t112 * t139 - t124 * t148 + t127 * t151) * t140 + (-t138 * t88 + t139 * t90 - t148 * t96 + t151 * t98) * t136 + (-t138 * t87 + t139 * t89 - t148 * t95 + t151 * t97) * t135;
t193 = (Icges(3,5) * t148 + Icges(4,5) * t138 + Icges(3,6) * t151 + Icges(4,6) * t139) * t140 + (t196 * t149 + t195 * t152) * t136 + (t195 * t149 - t196 * t152) * t135;
t189 = pkin(2) * t148;
t188 = pkin(2) * t151;
t133 = t149 * pkin(1) - pkin(5) * t152;
t83 = -qJ(3) * t152 + t188 * t149;
t186 = -t133 - t83;
t185 = Icges(2,4) * t149;
t180 = t138 * t149;
t179 = t138 * t152;
t147 = sin(qJ(4));
t178 = t147 * t152;
t177 = t149 * t147;
t150 = cos(qJ(4));
t176 = t149 * t150;
t175 = t150 * t152;
t174 = qJD(4) * t138;
t173 = V_base(5) * pkin(4) + V_base(1);
t170 = qJD(3) * t149 + t135 * t189 + t173;
t169 = pkin(3) * t139 + pkin(6) * t138;
t168 = rSges(3,1) * t151 - rSges(3,2) * t148;
t167 = rSges(4,1) * t139 - rSges(4,2) * t138;
t134 = pkin(1) * t152 + t149 * pkin(5);
t160 = -V_base(4) * pkin(4) + t140 * t134 + V_base(2);
t159 = V_base(4) * t133 - t134 * V_base(5) + V_base(3);
t158 = t136 * t83 + t159;
t84 = qJ(3) * t149 + t188 * t152;
t155 = -qJD(3) * t152 + t140 * t84 + t160;
t143 = Icges(2,4) * t152;
t132 = rSges(2,1) * t152 - t149 * rSges(2,2);
t131 = t149 * rSges(2,1) + rSges(2,2) * t152;
t130 = rSges(3,1) * t148 + rSges(3,2) * t151;
t129 = Icges(2,1) * t152 - t185;
t128 = Icges(2,1) * t149 + t143;
t126 = -Icges(2,2) * t149 + t143;
t125 = Icges(2,2) * t152 + t185;
t120 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t119 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t118 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t117 = -qJD(4) * t139 + t140;
t114 = pkin(3) * t138 - pkin(6) * t139;
t113 = rSges(4,1) * t138 + rSges(4,2) * t139;
t108 = t139 * t175 + t177;
t107 = -t139 * t178 + t176;
t106 = t139 * t176 - t178;
t105 = -t139 * t177 - t175;
t104 = t152 * t174 + t136;
t103 = t149 * t174 + t135;
t102 = t169 * t152;
t101 = t169 * t149;
t100 = t149 * rSges(3,3) + t168 * t152;
t99 = -rSges(3,3) * t152 + t168 * t149;
t92 = t149 * rSges(4,3) + t167 * t152;
t91 = -rSges(4,3) * t152 + t167 * t149;
t82 = -rSges(5,3) * t139 + (rSges(5,1) * t150 - rSges(5,2) * t147) * t138;
t81 = V_base(5) * rSges(2,3) - t131 * t140 + t173;
t80 = t132 * t140 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t79 = -Icges(5,5) * t139 + (Icges(5,1) * t150 - Icges(5,4) * t147) * t138;
t78 = -Icges(5,6) * t139 + (Icges(5,4) * t150 - Icges(5,2) * t147) * t138;
t77 = -Icges(5,3) * t139 + (Icges(5,5) * t150 - Icges(5,6) * t147) * t138;
t76 = t131 * V_base(4) - t132 * V_base(5) + V_base(3);
t73 = t108 * rSges(5,1) + t107 * rSges(5,2) + rSges(5,3) * t179;
t72 = rSges(5,1) * t106 + rSges(5,2) * t105 + rSges(5,3) * t180;
t71 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t179;
t70 = Icges(5,1) * t106 + Icges(5,4) * t105 + Icges(5,5) * t180;
t69 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t179;
t68 = Icges(5,4) * t106 + Icges(5,2) * t105 + Icges(5,6) * t180;
t67 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t179;
t66 = Icges(5,5) * t106 + Icges(5,6) * t105 + Icges(5,3) * t180;
t65 = t130 * t135 + (-t133 - t99) * t140 + t173;
t64 = t100 * t140 - t130 * t136 + t160;
t63 = -t100 * t135 + t136 * t99 + t159;
t62 = t113 * t135 + (-t91 + t186) * t140 + t170;
t61 = t140 * t92 + (-t113 - t189) * t136 + t155;
t60 = t136 * t91 + (-t84 - t92) * t135 + t158;
t59 = t103 * t82 + t114 * t135 - t117 * t72 + (-t101 + t186) * t140 + t170;
t58 = t140 * t102 - t104 * t82 + t117 * t73 + (-t114 - t189) * t136 + t155;
t57 = t101 * t136 - t103 * t73 + t104 * t72 + (-t102 - t84) * t135 + t158;
t1 = m(1) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(2) * (t76 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t104 * ((t107 * t69 + t108 * t71 + t67 * t179) * t104 + (t107 * t68 + t108 * t70 + t66 * t179) * t103 + (t107 * t78 + t108 * t79 + t77 * t179) * t117) / 0.2e1 + t103 * ((t105 * t69 + t106 * t71 + t67 * t180) * t104 + (t105 * t68 + t106 * t70 + t66 * t180) * t103 + (t105 * t78 + t106 * t79 + t77 * t180) * t117) / 0.2e1 + t117 * ((-t66 * t103 - t67 * t104 - t77 * t117) * t139 + ((-t147 * t69 + t150 * t71) * t104 + (-t147 * t68 + t150 * t70) * t103 + (-t147 * t78 + t150 * t79) * t117) * t138) / 0.2e1 + (t194 * t149 - t193 * t152) * t135 / 0.2e1 + (t193 * t149 + t194 * t152) * t136 / 0.2e1 + ((-t149 * t125 + t128 * t152 + Icges(1,4)) * V_base(5) + (-t149 * t126 + t152 * t129 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t152 * t125 + t149 * t128 + Icges(1,2)) * V_base(5) + (t126 * t152 + t149 * t129 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t138 * t90 + t139 * t88 + t148 * t98 + t151 * t96) * t136 + (t138 * t89 + t139 * t87 + t148 * t97 + t151 * t95) * t135 + (t139 * t111 + t138 * t112 + t151 * t124 + t148 * t127 + Icges(2,3)) * t140) * t140 / 0.2e1 + t140 * V_base(4) * (Icges(2,5) * t152 - Icges(2,6) * t149) + V_base(5) * t140 * (Icges(2,5) * t149 + Icges(2,6) * t152) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
