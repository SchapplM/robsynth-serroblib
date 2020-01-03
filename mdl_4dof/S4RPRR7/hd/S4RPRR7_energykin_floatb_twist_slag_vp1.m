% Calculate kinetic energy for
% S4RPRR7
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR7_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:40
% EndTime: 2019-12-31 16:53:41
% DurationCPUTime: 1.41s
% Computational Cost: add. (797->231), mult. (1022->339), div. (0->0), fcn. (910->8), ass. (0->116)
t146 = sin(pkin(7));
t190 = pkin(2) * t146;
t147 = cos(pkin(7));
t189 = pkin(2) * t147;
t150 = sin(qJ(1));
t152 = cos(qJ(1));
t130 = t150 * pkin(1) - qJ(2) * t152;
t83 = -pkin(5) * t152 + t150 * t189;
t188 = -t130 - t83;
t187 = Icges(2,4) * t150;
t186 = Icges(3,4) * t146;
t185 = Icges(3,4) * t147;
t145 = pkin(7) + qJ(3);
t138 = sin(t145);
t184 = Icges(4,4) * t138;
t139 = cos(t145);
t183 = Icges(4,4) * t139;
t182 = t138 * t150;
t181 = t138 * t152;
t149 = sin(qJ(4));
t180 = t149 * t152;
t179 = t150 * t149;
t151 = cos(qJ(4));
t178 = t150 * t151;
t177 = t151 * t152;
t175 = qJD(4) * t138;
t174 = V_base(4) * t130 + V_base(3);
t173 = V_base(5) * pkin(4) + V_base(1);
t135 = qJD(3) * t150 + V_base(4);
t140 = V_base(6) + qJD(1);
t170 = qJD(2) * t150 + t173;
t169 = V_base(5) * t190 + t170;
t168 = pkin(3) * t139 + pkin(6) * t138;
t134 = -qJD(3) * t152 + V_base(5);
t167 = rSges(3,1) * t147 - rSges(3,2) * t146;
t166 = rSges(4,1) * t139 - rSges(4,2) * t138;
t165 = Icges(3,1) * t147 - t186;
t164 = Icges(4,1) * t139 - t184;
t163 = -Icges(3,2) * t146 + t185;
t162 = -Icges(4,2) * t138 + t183;
t161 = Icges(3,5) * t147 - Icges(3,6) * t146;
t160 = Icges(4,5) * t139 - Icges(4,6) * t138;
t132 = pkin(1) * t152 + t150 * qJ(2);
t159 = -qJD(2) * t152 + t140 * t132 + V_base(2);
t158 = (Icges(4,5) * t138 + Icges(4,6) * t139) * t140 + (-Icges(4,3) * t152 + t150 * t160) * t134 + (Icges(4,3) * t150 + t152 * t160) * t135;
t84 = pkin(5) * t150 + t152 * t189;
t157 = V_base(4) * t83 + (-t132 - t84) * V_base(5) + t174;
t156 = (Icges(3,5) * t146 + Icges(3,6) * t147) * t140 + (-Icges(3,3) * t152 + t150 * t161) * V_base(5) + (Icges(3,3) * t150 + t152 * t161) * V_base(4);
t155 = t140 * t84 + (-pkin(4) - t190) * V_base(4) + t159;
t111 = Icges(4,2) * t139 + t184;
t112 = Icges(4,1) * t138 + t183;
t87 = -Icges(4,6) * t152 + t150 * t162;
t88 = Icges(4,6) * t150 + t152 * t162;
t89 = -Icges(4,5) * t152 + t150 * t164;
t90 = Icges(4,5) * t150 + t152 * t164;
t154 = (-t138 * t88 + t139 * t90) * t135 + (-t138 * t87 + t139 * t89) * t134 + (-t111 * t138 + t112 * t139) * t140;
t121 = Icges(3,2) * t147 + t186;
t122 = Icges(3,1) * t146 + t185;
t95 = -Icges(3,6) * t152 + t150 * t163;
t96 = Icges(3,6) * t150 + t152 * t163;
t97 = -Icges(3,5) * t152 + t150 * t165;
t98 = Icges(3,5) * t150 + t152 * t165;
t153 = (-t146 * t96 + t147 * t98) * V_base(4) + (-t146 * t95 + t147 * t97) * V_base(5) + (-t121 * t146 + t122 * t147) * t140;
t143 = Icges(2,4) * t152;
t133 = rSges(2,1) * t152 - t150 * rSges(2,2);
t131 = t150 * rSges(2,1) + rSges(2,2) * t152;
t129 = Icges(2,1) * t152 - t187;
t128 = Icges(2,1) * t150 + t143;
t127 = -Icges(2,2) * t150 + t143;
t126 = Icges(2,2) * t152 + t187;
t125 = Icges(2,5) * t152 - Icges(2,6) * t150;
t124 = Icges(2,5) * t150 + Icges(2,6) * t152;
t123 = rSges(3,1) * t146 + rSges(3,2) * t147;
t119 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t118 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t117 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t116 = -qJD(4) * t139 + t140;
t114 = pkin(3) * t138 - pkin(6) * t139;
t113 = rSges(4,1) * t138 + rSges(4,2) * t139;
t108 = t139 * t177 + t179;
t107 = -t139 * t180 + t178;
t106 = t139 * t178 - t180;
t105 = -t139 * t179 - t177;
t104 = t152 * t175 + t135;
t103 = t150 * t175 + t134;
t102 = t168 * t152;
t101 = t168 * t150;
t100 = t150 * rSges(3,3) + t152 * t167;
t99 = -rSges(3,3) * t152 + t150 * t167;
t92 = t150 * rSges(4,3) + t152 * t166;
t91 = -rSges(4,3) * t152 + t150 * t166;
t82 = -rSges(5,3) * t139 + (rSges(5,1) * t151 - rSges(5,2) * t149) * t138;
t81 = V_base(5) * rSges(2,3) - t131 * t140 + t173;
t80 = t133 * t140 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t79 = -Icges(5,5) * t139 + (Icges(5,1) * t151 - Icges(5,4) * t149) * t138;
t78 = -Icges(5,6) * t139 + (Icges(5,4) * t151 - Icges(5,2) * t149) * t138;
t77 = -Icges(5,3) * t139 + (Icges(5,5) * t151 - Icges(5,6) * t149) * t138;
t76 = t131 * V_base(4) - t133 * V_base(5) + V_base(3);
t73 = t108 * rSges(5,1) + t107 * rSges(5,2) + rSges(5,3) * t181;
t72 = rSges(5,1) * t106 + rSges(5,2) * t105 + rSges(5,3) * t182;
t71 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t181;
t70 = Icges(5,1) * t106 + Icges(5,4) * t105 + Icges(5,5) * t182;
t69 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t181;
t68 = Icges(5,4) * t106 + Icges(5,2) * t105 + Icges(5,6) * t182;
t67 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t181;
t66 = Icges(5,5) * t106 + Icges(5,6) * t105 + Icges(5,3) * t182;
t65 = t123 * V_base(5) + (-t130 - t99) * t140 + t170;
t64 = t140 * t100 + (-pkin(4) - t123) * V_base(4) + t159;
t63 = t99 * V_base(4) + (-t100 - t132) * V_base(5) + t174;
t62 = t113 * t134 + (-t91 + t188) * t140 + t169;
t61 = -t135 * t113 + t140 * t92 + t155;
t60 = -t134 * t92 + t135 * t91 + t157;
t59 = t103 * t82 + t114 * t134 - t116 * t72 + (-t101 + t188) * t140 + t169;
t58 = t140 * t102 - t104 * t82 - t135 * t114 + t116 * t73 + t155;
t57 = t101 * t135 - t102 * t134 - t103 * t73 + t104 * t72 + t157;
t1 = m(1) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(2) * (t76 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t135 * (t158 * t150 + t154 * t152) / 0.2e1 + t134 * (t154 * t150 - t158 * t152) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t104 * ((t107 * t69 + t108 * t71 + t67 * t181) * t104 + (t107 * t68 + t108 * t70 + t181 * t66) * t103 + (t107 * t78 + t108 * t79 + t181 * t77) * t116) / 0.2e1 + t103 * ((t105 * t69 + t106 * t71 + t182 * t67) * t104 + (t105 * t68 + t106 * t70 + t66 * t182) * t103 + (t105 * t78 + t106 * t79 + t182 * t77) * t116) / 0.2e1 + t116 * ((-t66 * t103 - t67 * t104 - t77 * t116) * t139 + ((-t149 * t69 + t151 * t71) * t104 + (-t149 * t68 + t151 * t70) * t103 + (-t149 * t78 + t151 * t79) * t116) * t138) / 0.2e1 + ((t138 * t90 + t139 * t88) * t135 + (t138 * t89 + t139 * t87) * t134 + (t146 * t97 + t147 * t95 + t124) * V_base(5) + (t146 * t98 + t147 * t96 + t125) * V_base(4) + (t139 * t111 + t138 * t112 + t147 * t121 + t146 * t122 + Icges(2,3)) * t140) * t140 / 0.2e1 + (t125 * t140 + t156 * t150 + t153 * t152 + (-t150 * t126 + t128 * t152 + Icges(1,4)) * V_base(5) + (-t150 * t127 + t152 * t129 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t124 * t140 + t153 * t150 - t156 * t152 + (t152 * t126 + t150 * t128 + Icges(1,2)) * V_base(5) + (t127 * t152 + t150 * t129 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
