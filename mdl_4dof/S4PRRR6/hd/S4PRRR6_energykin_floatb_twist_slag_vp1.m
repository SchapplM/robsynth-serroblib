% Calculate kinetic energy for
% S4PRRR6
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:36
% EndTime: 2019-12-31 16:34:38
% DurationCPUTime: 1.70s
% Computational Cost: add. (802->252), mult. (1330->386), div. (0->0), fcn. (1274->8), ass. (0->118)
t167 = cos(qJ(3));
t202 = pkin(3) * t167;
t163 = sin(pkin(7));
t200 = Icges(2,4) * t163;
t166 = sin(qJ(2));
t199 = Icges(3,4) * t166;
t168 = cos(qJ(2));
t198 = Icges(3,4) * t168;
t165 = sin(qJ(3));
t197 = t163 * t165;
t196 = t163 * t166;
t195 = t163 * t168;
t164 = cos(pkin(7));
t194 = t164 * t165;
t193 = t164 * t166;
t192 = t164 * t168;
t191 = t165 * t168;
t190 = t167 * t168;
t189 = qJD(3) * t166;
t188 = qJD(4) * t166;
t187 = V_base(5) * qJ(1) + V_base(1);
t183 = qJD(1) + V_base(3);
t152 = qJD(2) * t163 + V_base(4);
t123 = t164 * t189 + t152;
t182 = pkin(2) * t168 + pkin(5) * t166;
t151 = -qJD(2) * t164 + V_base(5);
t181 = rSges(3,1) * t168 - rSges(3,2) * t166;
t180 = Icges(3,1) * t168 - t199;
t179 = -Icges(3,2) * t166 + t198;
t178 = Icges(3,5) * t168 - Icges(3,6) * t166;
t122 = t163 * t189 + t151;
t145 = pkin(1) * t164 + pkin(4) * t163;
t177 = -V_base(4) * qJ(1) + V_base(6) * t145 + V_base(2);
t144 = pkin(1) * t163 - pkin(4) * t164;
t176 = V_base(4) * t144 - t145 * V_base(5) + t183;
t175 = pkin(6) * t166 + t168 * t202;
t128 = t182 * t163;
t150 = t166 * pkin(2) - t168 * pkin(5);
t174 = t151 * t150 + (-t128 - t144) * V_base(6) + t187;
t173 = (-Icges(3,3) * t164 + t163 * t178) * t151 + (Icges(3,3) * t163 + t164 * t178) * t152 + (Icges(3,5) * t166 + Icges(3,6) * t168) * V_base(6);
t129 = t182 * t164;
t172 = V_base(6) * t129 - t150 * t152 + t177;
t171 = t152 * t128 - t129 * t151 + t176;
t106 = -Icges(3,6) * t164 + t163 * t179;
t107 = Icges(3,6) * t163 + t164 * t179;
t108 = -Icges(3,5) * t164 + t163 * t180;
t109 = Icges(3,5) * t163 + t164 * t180;
t147 = Icges(3,2) * t168 + t199;
t148 = Icges(3,1) * t166 + t198;
t170 = (-t107 * t166 + t109 * t168) * t152 + (-t106 * t166 + t108 * t168) * t151 + (-t147 * t166 + t148 * t168) * V_base(6);
t162 = qJ(3) + qJ(4);
t160 = cos(t162);
t159 = sin(t162);
t158 = Icges(2,4) * t164;
t153 = -qJD(3) * t168 + V_base(6);
t149 = rSges(3,1) * t166 + rSges(3,2) * t168;
t143 = rSges(2,1) * t164 - rSges(2,2) * t163;
t142 = rSges(2,1) * t163 + rSges(2,2) * t164;
t141 = Icges(2,1) * t164 - t200;
t140 = Icges(2,1) * t163 + t158;
t139 = -Icges(2,2) * t163 + t158;
t138 = Icges(2,2) * t164 + t200;
t135 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t134 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t133 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t132 = V_base(6) + (-qJD(3) - qJD(4)) * t168;
t127 = t164 * t190 + t197;
t126 = t163 * t167 - t164 * t191;
t125 = t163 * t190 - t194;
t124 = -t163 * t191 - t164 * t167;
t119 = -rSges(4,3) * t168 + (rSges(4,1) * t167 - rSges(4,2) * t165) * t166;
t118 = -Icges(4,5) * t168 + (Icges(4,1) * t167 - Icges(4,4) * t165) * t166;
t117 = -Icges(4,6) * t168 + (Icges(4,4) * t167 - Icges(4,2) * t165) * t166;
t116 = -Icges(4,3) * t168 + (Icges(4,5) * t167 - Icges(4,6) * t165) * t166;
t115 = t159 * t163 + t160 * t192;
t114 = -t159 * t192 + t160 * t163;
t113 = -t159 * t164 + t160 * t195;
t112 = -t159 * t195 - t160 * t164;
t111 = rSges(3,3) * t163 + t164 * t181;
t110 = -rSges(3,3) * t164 + t163 * t181;
t103 = -rSges(5,3) * t168 + (rSges(5,1) * t160 - rSges(5,2) * t159) * t166;
t102 = -Icges(5,5) * t168 + (Icges(5,1) * t160 - Icges(5,4) * t159) * t166;
t101 = -Icges(5,6) * t168 + (Icges(5,4) * t160 - Icges(5,2) * t159) * t166;
t100 = -Icges(5,3) * t168 + (Icges(5,5) * t160 - Icges(5,6) * t159) * t166;
t99 = t164 * t188 + t123;
t98 = t163 * t188 + t122;
t97 = -pkin(6) * t168 + t166 * t202;
t95 = V_base(5) * rSges(2,3) - t142 * V_base(6) + t187;
t94 = t143 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t93 = t142 * V_base(4) - t143 * V_base(5) + t183;
t92 = pkin(3) * t197 + t164 * t175;
t91 = -pkin(3) * t194 + t163 * t175;
t90 = rSges(4,1) * t127 + rSges(4,2) * t126 + rSges(4,3) * t193;
t89 = rSges(4,1) * t125 + rSges(4,2) * t124 + rSges(4,3) * t196;
t88 = Icges(4,1) * t127 + Icges(4,4) * t126 + Icges(4,5) * t193;
t87 = Icges(4,1) * t125 + Icges(4,4) * t124 + Icges(4,5) * t196;
t86 = Icges(4,4) * t127 + Icges(4,2) * t126 + Icges(4,6) * t193;
t85 = Icges(4,4) * t125 + Icges(4,2) * t124 + Icges(4,6) * t196;
t84 = Icges(4,5) * t127 + Icges(4,6) * t126 + Icges(4,3) * t193;
t83 = Icges(4,5) * t125 + Icges(4,6) * t124 + Icges(4,3) * t196;
t82 = rSges(5,1) * t115 + rSges(5,2) * t114 + rSges(5,3) * t193;
t81 = rSges(5,1) * t113 + rSges(5,2) * t112 + rSges(5,3) * t196;
t80 = Icges(5,1) * t115 + Icges(5,4) * t114 + Icges(5,5) * t193;
t79 = Icges(5,1) * t113 + Icges(5,4) * t112 + Icges(5,5) * t196;
t78 = Icges(5,4) * t115 + Icges(5,2) * t114 + Icges(5,6) * t193;
t77 = Icges(5,4) * t113 + Icges(5,2) * t112 + Icges(5,6) * t196;
t76 = Icges(5,5) * t115 + Icges(5,6) * t114 + Icges(5,3) * t193;
t75 = Icges(5,5) * t113 + Icges(5,6) * t112 + Icges(5,3) * t196;
t74 = t149 * t151 + (-t110 - t144) * V_base(6) + t187;
t73 = t111 * V_base(6) - t149 * t152 + t177;
t72 = t110 * t152 - t111 * t151 + t176;
t71 = t119 * t122 - t153 * t89 + t174;
t70 = -t119 * t123 + t153 * t90 + t172;
t69 = -t122 * t90 + t123 * t89 + t171;
t68 = t103 * t98 + t122 * t97 - t132 * t81 - t153 * t91 + t174;
t67 = -t103 * t99 - t123 * t97 + t132 * t82 + t153 * t92 + t172;
t66 = -t122 * t92 + t123 * t91 + t81 * t99 - t82 * t98 + t171;
t1 = m(1) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(2) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + t152 * (t173 * t163 + t170 * t164) / 0.2e1 + t151 * (t170 * t163 - t173 * t164) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t123 * ((t126 * t86 + t127 * t88 + t84 * t193) * t123 + (t126 * t85 + t127 * t87 + t193 * t83) * t122 + (t116 * t193 + t117 * t126 + t118 * t127) * t153) / 0.2e1 + t122 * ((t124 * t86 + t125 * t88 + t196 * t84) * t123 + (t124 * t85 + t125 * t87 + t83 * t196) * t122 + (t116 * t196 + t117 * t124 + t118 * t125) * t153) / 0.2e1 + t153 * ((-t116 * t153 - t83 * t122 - t84 * t123) * t168 + ((-t165 * t86 + t167 * t88) * t123 + (-t165 * t85 + t167 * t87) * t122 + (-t117 * t165 + t118 * t167) * t153) * t166) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t99 * ((t114 * t78 + t115 * t80 + t76 * t193) * t99 + (t114 * t77 + t115 * t79 + t193 * t75) * t98 + (t100 * t193 + t101 * t114 + t102 * t115) * t132) / 0.2e1 + t98 * ((t112 * t78 + t113 * t80 + t196 * t76) * t99 + (t112 * t77 + t113 * t79 + t75 * t196) * t98 + (t100 * t196 + t101 * t112 + t102 * t113) * t132) / 0.2e1 + t132 * ((-t100 * t132 - t75 * t98 - t76 * t99) * t168 + ((-t159 * t78 + t160 * t80) * t99 + (-t159 * t77 + t160 * t79) * t98 + (-t101 * t159 + t102 * t160) * t132) * t166) / 0.2e1 + ((-t138 * t163 + t140 * t164 + Icges(1,4)) * V_base(5) + (-t139 * t163 + t141 * t164 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t138 * t164 + t140 * t163 + Icges(1,2)) * V_base(5) + (t139 * t164 + t141 * t163 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t107 * t168 + t109 * t166) * t152 + (t106 * t168 + t108 * t166) * t151 + (t147 * t168 + t148 * t166 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t164 - Icges(2,6) * t163 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t163 + Icges(2,6) * t164 + Icges(1,6));
T = t1;
