% Calculate kinetic energy for
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:55
% EndTime: 2019-12-31 16:22:57
% DurationCPUTime: 1.44s
% Computational Cost: add. (783->225), mult. (1040->320), div. (0->0), fcn. (928->8), ass. (0->111)
t198 = Icges(3,3) + Icges(4,3);
t144 = qJ(2) + pkin(7);
t140 = sin(t144);
t141 = cos(t144);
t149 = sin(qJ(2));
t151 = cos(qJ(2));
t197 = Icges(3,5) * t151 + Icges(4,5) * t141 - Icges(3,6) * t149 - Icges(4,6) * t140;
t183 = Icges(4,4) * t140;
t110 = Icges(4,2) * t141 + t183;
t182 = Icges(4,4) * t141;
t111 = Icges(4,1) * t140 + t182;
t185 = Icges(3,4) * t149;
t132 = Icges(3,2) * t151 + t185;
t184 = Icges(3,4) * t151;
t133 = Icges(3,1) * t149 + t184;
t146 = cos(pkin(6));
t135 = -qJD(2) * t146 + V_base(5);
t145 = sin(pkin(6));
t136 = qJD(2) * t145 + V_base(4);
t162 = -Icges(4,2) * t140 + t182;
t85 = -Icges(4,6) * t146 + t162 * t145;
t86 = Icges(4,6) * t145 + t162 * t146;
t164 = Icges(4,1) * t141 - t183;
t87 = -Icges(4,5) * t146 + t164 * t145;
t88 = Icges(4,5) * t145 + t164 * t146;
t163 = -Icges(3,2) * t149 + t184;
t95 = -Icges(3,6) * t146 + t163 * t145;
t96 = Icges(3,6) * t145 + t163 * t146;
t165 = Icges(3,1) * t151 - t185;
t97 = -Icges(3,5) * t146 + t165 * t145;
t98 = Icges(3,5) * t145 + t165 * t146;
t194 = (-t110 * t140 + t111 * t141 - t132 * t149 + t133 * t151) * V_base(6) + (-t140 * t86 + t141 * t88 - t149 * t96 + t151 * t98) * t136 + (-t140 * t85 + t141 * t87 - t149 * t95 + t151 * t97) * t135;
t193 = (Icges(3,5) * t149 + Icges(4,5) * t140 + Icges(3,6) * t151 + Icges(4,6) * t141) * V_base(6) + (t198 * t145 + t197 * t146) * t136 + (t197 * t145 - t198 * t146) * t135;
t190 = pkin(2) * t149;
t189 = pkin(2) * t151;
t129 = pkin(1) * t145 - pkin(4) * t146;
t81 = -qJ(3) * t146 + t189 * t145;
t187 = -t129 - t81;
t186 = Icges(2,4) * t145;
t181 = t140 * t145;
t180 = t140 * t146;
t148 = sin(qJ(4));
t179 = t145 * t148;
t150 = cos(qJ(4));
t178 = t145 * t150;
t177 = t146 * t148;
t176 = t146 * t150;
t175 = qJD(4) * t140;
t174 = V_base(5) * qJ(1) + V_base(1);
t170 = qJD(1) + V_base(3);
t169 = qJD(3) * t145 + t135 * t190 + t174;
t168 = pkin(3) * t141 + pkin(5) * t140;
t167 = rSges(3,1) * t151 - rSges(3,2) * t149;
t166 = rSges(4,1) * t141 - rSges(4,2) * t140;
t130 = pkin(1) * t146 + pkin(4) * t145;
t159 = -V_base(4) * qJ(1) + V_base(6) * t130 + V_base(2);
t158 = V_base(4) * t129 - t130 * V_base(5) + t170;
t157 = t136 * t81 + t158;
t82 = qJ(3) * t145 + t189 * t146;
t154 = -qJD(3) * t146 + V_base(6) * t82 + t159;
t142 = Icges(2,4) * t146;
t134 = t149 * rSges(3,1) + rSges(3,2) * t151;
t128 = rSges(2,1) * t146 - rSges(2,2) * t145;
t127 = rSges(2,1) * t145 + rSges(2,2) * t146;
t126 = -qJD(4) * t141 + V_base(6);
t125 = Icges(2,1) * t146 - t186;
t124 = Icges(2,1) * t145 + t142;
t123 = -Icges(2,2) * t145 + t142;
t122 = Icges(2,2) * t146 + t186;
t119 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t118 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t117 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t113 = pkin(3) * t140 - pkin(5) * t141;
t112 = rSges(4,1) * t140 + rSges(4,2) * t141;
t108 = t141 * t176 + t179;
t107 = -t141 * t177 + t178;
t106 = t141 * t178 - t177;
t105 = -t141 * t179 - t176;
t104 = t146 * t175 + t136;
t103 = t145 * t175 + t135;
t102 = t168 * t146;
t101 = t168 * t145;
t100 = t145 * rSges(3,3) + t167 * t146;
t99 = -t146 * rSges(3,3) + t167 * t145;
t92 = rSges(4,3) * t145 + t166 * t146;
t91 = -rSges(4,3) * t146 + t166 * t145;
t90 = V_base(5) * rSges(2,3) - t127 * V_base(6) + t174;
t89 = t128 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t80 = -rSges(5,3) * t141 + (rSges(5,1) * t150 - rSges(5,2) * t148) * t140;
t79 = -Icges(5,5) * t141 + (Icges(5,1) * t150 - Icges(5,4) * t148) * t140;
t78 = -Icges(5,6) * t141 + (Icges(5,4) * t150 - Icges(5,2) * t148) * t140;
t77 = -Icges(5,3) * t141 + (Icges(5,5) * t150 - Icges(5,6) * t148) * t140;
t75 = t127 * V_base(4) - t128 * V_base(5) + t170;
t73 = rSges(5,1) * t108 + rSges(5,2) * t107 + rSges(5,3) * t180;
t72 = rSges(5,1) * t106 + rSges(5,2) * t105 + rSges(5,3) * t181;
t71 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t180;
t70 = Icges(5,1) * t106 + Icges(5,4) * t105 + Icges(5,5) * t181;
t69 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t180;
t68 = Icges(5,4) * t106 + Icges(5,2) * t105 + Icges(5,6) * t181;
t67 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t180;
t66 = Icges(5,5) * t106 + Icges(5,6) * t105 + Icges(5,3) * t181;
t65 = t134 * t135 + (-t129 - t99) * V_base(6) + t174;
t64 = t100 * V_base(6) - t134 * t136 + t159;
t63 = -t100 * t135 + t136 * t99 + t158;
t62 = t112 * t135 + (-t91 + t187) * V_base(6) + t169;
t61 = t92 * V_base(6) + (-t112 - t190) * t136 + t154;
t60 = t136 * t91 + (-t82 - t92) * t135 + t157;
t59 = t103 * t80 + t113 * t135 - t126 * t72 + (-t101 + t187) * V_base(6) + t169;
t58 = t102 * V_base(6) - t104 * t80 + t126 * t73 + (-t113 - t190) * t136 + t154;
t57 = t101 * t136 - t103 * t73 + t104 * t72 + (-t102 - t82) * t135 + t157;
t1 = m(1) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(2) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t104 * ((t107 * t69 + t108 * t71 + t67 * t180) * t104 + (t107 * t68 + t108 * t70 + t66 * t180) * t103 + (t107 * t78 + t108 * t79 + t77 * t180) * t126) / 0.2e1 + t103 * ((t105 * t69 + t106 * t71 + t67 * t181) * t104 + (t105 * t68 + t106 * t70 + t66 * t181) * t103 + (t105 * t78 + t106 * t79 + t77 * t181) * t126) / 0.2e1 + t126 * ((-t66 * t103 - t67 * t104 - t77 * t126) * t141 + ((-t148 * t69 + t150 * t71) * t104 + (-t148 * t68 + t150 * t70) * t103 + (-t148 * t78 + t150 * t79) * t126) * t140) / 0.2e1 + (t194 * t145 - t193 * t146) * t135 / 0.2e1 + (t193 * t145 + t194 * t146) * t136 / 0.2e1 + ((-t122 * t145 + t124 * t146 + Icges(1,4)) * V_base(5) + (-t145 * t123 + t146 * t125 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t146 * t122 + t145 * t124 + Icges(1,2)) * V_base(5) + (t123 * t146 + t125 * t145 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t140 * t88 + t141 * t86 + t149 * t98 + t151 * t96) * t136 + (t140 * t87 + t141 * t85 + t149 * t97 + t151 * t95) * t135 + (t141 * t110 + t140 * t111 + t151 * t132 + t149 * t133 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t146 - Icges(2,6) * t145 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t145 + Icges(2,6) * t146 + Icges(1,6));
T = t1;
