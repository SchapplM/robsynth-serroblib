% Calculate kinetic energy for
% S4RRRR5
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:19
% EndTime: 2019-12-31 17:27:21
% DurationCPUTime: 1.79s
% Computational Cost: add. (834->252), mult. (1330->389), div. (0->0), fcn. (1274->8), ass. (0->116)
t160 = cos(qJ(3));
t192 = pkin(3) * t160;
t159 = sin(qJ(1));
t190 = Icges(2,4) * t159;
t158 = sin(qJ(2));
t189 = Icges(3,4) * t158;
t161 = cos(qJ(2));
t188 = Icges(3,4) * t161;
t157 = sin(qJ(3));
t187 = t157 * t159;
t162 = cos(qJ(1));
t186 = t157 * t162;
t185 = t158 * t159;
t184 = t158 * t162;
t183 = t159 * t161;
t182 = t161 * t162;
t181 = qJD(3) * t158;
t180 = qJD(4) * t158;
t179 = V_base(5) * pkin(4) + V_base(1);
t146 = qJD(2) * t159 + V_base(4);
t150 = V_base(6) + qJD(1);
t116 = t162 * t181 + t146;
t176 = pkin(2) * t161 + pkin(6) * t158;
t145 = -qJD(2) * t162 + V_base(5);
t175 = rSges(3,1) * t161 - rSges(3,2) * t158;
t174 = Icges(3,1) * t161 - t189;
t173 = -Icges(3,2) * t158 + t188;
t172 = Icges(3,5) * t161 - Icges(3,6) * t158;
t115 = t159 * t181 + t145;
t144 = pkin(1) * t162 + pkin(5) * t159;
t171 = -V_base(4) * pkin(4) + t150 * t144 + V_base(2);
t143 = pkin(1) * t159 - pkin(5) * t162;
t170 = V_base(4) * t143 - t144 * V_base(5) + V_base(3);
t169 = (Icges(3,3) * t159 + t162 * t172) * t146 + (Icges(3,5) * t158 + Icges(3,6) * t161) * t150 + (-Icges(3,3) * t162 + t159 * t172) * t145;
t168 = pkin(7) * t158 + t161 * t192;
t122 = t176 * t159;
t142 = t158 * pkin(2) - t161 * pkin(6);
t167 = t145 * t142 + (-t122 - t143) * t150 + t179;
t123 = t176 * t162;
t166 = t150 * t123 - t142 * t146 + t171;
t165 = t146 * t122 - t123 * t145 + t170;
t102 = -Icges(3,6) * t162 + t159 * t173;
t103 = Icges(3,6) * t159 + t162 * t173;
t105 = -Icges(3,5) * t162 + t159 * t174;
t106 = Icges(3,5) * t159 + t162 * t174;
t132 = Icges(3,2) * t161 + t189;
t135 = Icges(3,1) * t158 + t188;
t164 = (-t103 * t158 + t106 * t161) * t146 + (-t102 * t158 + t105 * t161) * t145 + (-t132 * t158 + t135 * t161) * t150;
t156 = qJ(3) + qJ(4);
t154 = Icges(2,4) * t162;
t153 = cos(t156);
t152 = sin(t156);
t141 = rSges(2,1) * t162 - rSges(2,2) * t159;
t140 = rSges(2,1) * t159 + rSges(2,2) * t162;
t139 = rSges(3,1) * t158 + rSges(3,2) * t161;
t138 = -qJD(3) * t161 + t150;
t137 = Icges(2,1) * t162 - t190;
t136 = Icges(2,1) * t159 + t154;
t134 = -Icges(2,2) * t159 + t154;
t133 = Icges(2,2) * t162 + t190;
t128 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t127 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t126 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t124 = (-qJD(3) - qJD(4)) * t161 + t150;
t120 = t160 * t182 + t187;
t119 = -t157 * t182 + t159 * t160;
t118 = t160 * t183 - t186;
t117 = -t157 * t183 - t160 * t162;
t113 = t152 * t159 + t153 * t182;
t112 = -t152 * t182 + t153 * t159;
t111 = -t152 * t162 + t153 * t183;
t110 = -t152 * t183 - t153 * t162;
t109 = rSges(3,3) * t159 + t162 * t175;
t108 = -rSges(3,3) * t162 + t159 * t175;
t107 = -rSges(4,3) * t161 + (rSges(4,1) * t160 - rSges(4,2) * t157) * t158;
t104 = -Icges(4,5) * t161 + (Icges(4,1) * t160 - Icges(4,4) * t157) * t158;
t101 = -Icges(4,6) * t161 + (Icges(4,4) * t160 - Icges(4,2) * t157) * t158;
t98 = -Icges(4,3) * t161 + (Icges(4,5) * t160 - Icges(4,6) * t157) * t158;
t96 = t162 * t180 + t116;
t95 = t159 * t180 + t115;
t94 = -rSges(5,3) * t161 + (rSges(5,1) * t153 - rSges(5,2) * t152) * t158;
t92 = -Icges(5,5) * t161 + (Icges(5,1) * t153 - Icges(5,4) * t152) * t158;
t91 = -Icges(5,6) * t161 + (Icges(5,4) * t153 - Icges(5,2) * t152) * t158;
t90 = -Icges(5,3) * t161 + (Icges(5,5) * t153 - Icges(5,6) * t152) * t158;
t89 = -pkin(7) * t161 + t158 * t192;
t88 = V_base(5) * rSges(2,3) - t140 * t150 + t179;
t87 = t141 * t150 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t86 = t140 * V_base(4) - t141 * V_base(5) + V_base(3);
t85 = pkin(3) * t187 + t162 * t168;
t84 = -pkin(3) * t186 + t159 * t168;
t83 = rSges(4,1) * t120 + rSges(4,2) * t119 + rSges(4,3) * t184;
t82 = rSges(4,1) * t118 + rSges(4,2) * t117 + rSges(4,3) * t185;
t81 = Icges(4,1) * t120 + Icges(4,4) * t119 + Icges(4,5) * t184;
t80 = Icges(4,1) * t118 + Icges(4,4) * t117 + Icges(4,5) * t185;
t79 = Icges(4,4) * t120 + Icges(4,2) * t119 + Icges(4,6) * t184;
t78 = Icges(4,4) * t118 + Icges(4,2) * t117 + Icges(4,6) * t185;
t77 = Icges(4,5) * t120 + Icges(4,6) * t119 + Icges(4,3) * t184;
t76 = Icges(4,5) * t118 + Icges(4,6) * t117 + Icges(4,3) * t185;
t75 = rSges(5,1) * t113 + rSges(5,2) * t112 + rSges(5,3) * t184;
t74 = rSges(5,1) * t111 + rSges(5,2) * t110 + rSges(5,3) * t185;
t73 = Icges(5,1) * t113 + Icges(5,4) * t112 + Icges(5,5) * t184;
t72 = Icges(5,1) * t111 + Icges(5,4) * t110 + Icges(5,5) * t185;
t71 = Icges(5,4) * t113 + Icges(5,2) * t112 + Icges(5,6) * t184;
t70 = Icges(5,4) * t111 + Icges(5,2) * t110 + Icges(5,6) * t185;
t69 = Icges(5,5) * t113 + Icges(5,6) * t112 + Icges(5,3) * t184;
t68 = Icges(5,5) * t111 + Icges(5,6) * t110 + Icges(5,3) * t185;
t67 = t139 * t145 + (-t108 - t143) * t150 + t179;
t66 = t109 * t150 - t139 * t146 + t171;
t65 = t108 * t146 - t109 * t145 + t170;
t64 = t107 * t115 - t138 * t82 + t167;
t63 = -t107 * t116 + t138 * t83 + t166;
t62 = -t115 * t83 + t116 * t82 + t165;
t61 = t115 * t89 - t124 * t74 - t138 * t84 + t94 * t95 + t167;
t60 = -t116 * t89 + t124 * t75 + t138 * t85 - t94 * t96 + t166;
t59 = -t115 * t85 + t116 * t84 + t74 * t96 - t75 * t95 + t165;
t1 = m(1) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(2) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t146 * (t169 * t159 + t164 * t162) / 0.2e1 + t145 * (t164 * t159 - t169 * t162) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + t116 * ((t119 * t79 + t120 * t81 + t77 * t184) * t116 + (t119 * t78 + t120 * t80 + t184 * t76) * t115 + (t101 * t119 + t104 * t120 + t184 * t98) * t138) / 0.2e1 + t115 * ((t117 * t79 + t118 * t81 + t185 * t77) * t116 + (t117 * t78 + t118 * t80 + t76 * t185) * t115 + (t101 * t117 + t104 * t118 + t185 * t98) * t138) / 0.2e1 + t138 * ((-t76 * t115 - t77 * t116 - t98 * t138) * t161 + ((-t157 * t79 + t160 * t81) * t116 + (-t157 * t78 + t160 * t80) * t115 + (-t101 * t157 + t104 * t160) * t138) * t158) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + t96 * ((t112 * t71 + t113 * t73 + t69 * t184) * t96 + (t112 * t70 + t113 * t72 + t184 * t68) * t95 + (t112 * t91 + t113 * t92 + t184 * t90) * t124) / 0.2e1 + t95 * ((t110 * t71 + t111 * t73 + t185 * t69) * t96 + (t110 * t70 + t111 * t72 + t68 * t185) * t95 + (t110 * t91 + t111 * t92 + t185 * t90) * t124) / 0.2e1 + t124 * ((-t90 * t124 - t68 * t95 - t69 * t96) * t161 + ((-t152 * t71 + t153 * t73) * t96 + (-t152 * t70 + t153 * t72) * t95 + (-t152 * t91 + t153 * t92) * t124) * t158) / 0.2e1 + ((t103 * t161 + t106 * t158) * t146 + (t102 * t161 + t105 * t158) * t145 + (t132 * t161 + t135 * t158 + Icges(2,3)) * t150) * t150 / 0.2e1 + ((-t133 * t159 + t136 * t162 + Icges(1,4)) * V_base(5) + (-t134 * t159 + t137 * t162 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t133 * t162 + t136 * t159 + Icges(1,2)) * V_base(5) + (t134 * t162 + t137 * t159 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t150 * (Icges(2,5) * t162 - Icges(2,6) * t159) + V_base(5) * t150 * (Icges(2,5) * t159 + Icges(2,6) * t162) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
