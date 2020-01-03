% Calculate kinetic energy for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:23
% EndTime: 2019-12-31 16:18:24
% DurationCPUTime: 1.14s
% Computational Cost: add. (765->231), mult. (1022->334), div. (0->0), fcn. (910->8), ass. (0->116)
t146 = sin(pkin(6));
t148 = cos(pkin(6));
t193 = Icges(2,5) * t148 - Icges(2,6) * t146 + Icges(1,5);
t192 = Icges(2,5) * t146 + Icges(2,6) * t148 + Icges(1,6);
t145 = sin(pkin(7));
t191 = pkin(2) * t145;
t147 = cos(pkin(7));
t190 = pkin(2) * t147;
t130 = pkin(1) * t146 - qJ(2) * t148;
t77 = -pkin(4) * t148 + t190 * t146;
t189 = -t130 - t77;
t188 = Icges(2,4) * t146;
t187 = Icges(3,4) * t145;
t186 = Icges(3,4) * t147;
t144 = pkin(7) + qJ(3);
t140 = sin(t144);
t185 = Icges(4,4) * t140;
t141 = cos(t144);
t184 = Icges(4,4) * t141;
t183 = t140 * t146;
t182 = t140 * t148;
t150 = sin(qJ(4));
t181 = t146 * t150;
t151 = cos(qJ(4));
t180 = t146 * t151;
t179 = t148 * t150;
t178 = t148 * t151;
t176 = qJD(4) * t140;
t175 = V_base(5) * qJ(1) + V_base(1);
t171 = qJD(1) + V_base(3);
t135 = qJD(3) * t146 + V_base(4);
t170 = qJD(2) * t146 + t175;
t169 = V_base(4) * t130 + t171;
t168 = V_base(5) * t191 + t170;
t167 = pkin(3) * t141 + pkin(5) * t140;
t134 = -qJD(3) * t148 + V_base(5);
t166 = rSges(3,1) * t147 - rSges(3,2) * t145;
t165 = rSges(4,1) * t141 - rSges(4,2) * t140;
t164 = Icges(3,1) * t147 - t187;
t163 = Icges(4,1) * t141 - t185;
t162 = -Icges(3,2) * t145 + t186;
t161 = -Icges(4,2) * t140 + t184;
t160 = Icges(3,5) * t147 - Icges(3,6) * t145;
t159 = Icges(4,5) * t141 - Icges(4,6) * t140;
t132 = pkin(1) * t148 + qJ(2) * t146;
t158 = -qJD(2) * t148 + V_base(6) * t132 + V_base(2);
t157 = (Icges(4,5) * t140 + Icges(4,6) * t141) * V_base(6) + (-Icges(4,3) * t148 + t159 * t146) * t134 + (Icges(4,3) * t146 + t159 * t148) * t135;
t78 = pkin(4) * t146 + t190 * t148;
t156 = V_base(4) * t77 + (-t132 - t78) * V_base(5) + t169;
t155 = (Icges(3,5) * t145 + Icges(3,6) * t147) * V_base(6) + (-Icges(3,3) * t148 + t160 * t146) * V_base(5) + (Icges(3,3) * t146 + t160 * t148) * V_base(4);
t154 = V_base(6) * t78 + (-qJ(1) - t191) * V_base(4) + t158;
t110 = Icges(4,2) * t141 + t185;
t111 = Icges(4,1) * t140 + t184;
t85 = -Icges(4,6) * t148 + t161 * t146;
t86 = Icges(4,6) * t146 + t161 * t148;
t87 = -Icges(4,5) * t148 + t163 * t146;
t88 = Icges(4,5) * t146 + t163 * t148;
t153 = (-t140 * t86 + t141 * t88) * t135 + (-t140 * t85 + t141 * t87) * t134 + (-t110 * t140 + t111 * t141) * V_base(6);
t122 = Icges(3,2) * t147 + t187;
t125 = Icges(3,1) * t145 + t186;
t95 = -Icges(3,6) * t148 + t162 * t146;
t96 = Icges(3,6) * t146 + t162 * t148;
t97 = -Icges(3,5) * t148 + t164 * t146;
t98 = Icges(3,5) * t146 + t164 * t148;
t152 = (-t145 * t96 + t147 * t98) * V_base(4) + (-t145 * t95 + t147 * t97) * V_base(5) + (-t122 * t145 + t125 * t147) * V_base(6);
t142 = Icges(2,4) * t148;
t133 = rSges(2,1) * t148 - rSges(2,2) * t146;
t131 = rSges(2,1) * t146 + rSges(2,2) * t148;
t129 = rSges(3,1) * t145 + rSges(3,2) * t147;
t128 = -qJD(4) * t141 + V_base(6);
t127 = Icges(2,1) * t148 - t188;
t126 = Icges(2,1) * t146 + t142;
t124 = -Icges(2,2) * t146 + t142;
t123 = Icges(2,2) * t148 + t188;
t118 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t117 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t116 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t113 = pkin(3) * t140 - pkin(5) * t141;
t112 = rSges(4,1) * t140 + rSges(4,2) * t141;
t108 = t141 * t178 + t181;
t107 = -t141 * t179 + t180;
t106 = t141 * t180 - t179;
t105 = -t141 * t181 - t178;
t104 = t148 * t176 + t135;
t103 = t146 * t176 + t134;
t102 = t167 * t148;
t101 = t167 * t146;
t100 = rSges(3,3) * t146 + t166 * t148;
t99 = -rSges(3,3) * t148 + t166 * t146;
t92 = rSges(4,3) * t146 + t165 * t148;
t91 = -rSges(4,3) * t148 + t165 * t146;
t90 = V_base(5) * rSges(2,3) - t131 * V_base(6) + t175;
t89 = t133 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t82 = -t141 * rSges(5,3) + (rSges(5,1) * t151 - rSges(5,2) * t150) * t140;
t81 = -Icges(5,5) * t141 + (Icges(5,1) * t151 - Icges(5,4) * t150) * t140;
t80 = -Icges(5,6) * t141 + (Icges(5,4) * t151 - Icges(5,2) * t150) * t140;
t79 = -Icges(5,3) * t141 + (Icges(5,5) * t151 - Icges(5,6) * t150) * t140;
t74 = t131 * V_base(4) - t133 * V_base(5) + t171;
t73 = rSges(5,1) * t108 + rSges(5,2) * t107 + rSges(5,3) * t182;
t72 = rSges(5,1) * t106 + rSges(5,2) * t105 + rSges(5,3) * t183;
t71 = Icges(5,1) * t108 + Icges(5,4) * t107 + Icges(5,5) * t182;
t70 = Icges(5,1) * t106 + Icges(5,4) * t105 + Icges(5,5) * t183;
t69 = Icges(5,4) * t108 + Icges(5,2) * t107 + Icges(5,6) * t182;
t68 = Icges(5,4) * t106 + Icges(5,2) * t105 + Icges(5,6) * t183;
t67 = Icges(5,5) * t108 + Icges(5,6) * t107 + Icges(5,3) * t182;
t66 = Icges(5,5) * t106 + Icges(5,6) * t105 + Icges(5,3) * t183;
t65 = t129 * V_base(5) + (-t130 - t99) * V_base(6) + t170;
t64 = t100 * V_base(6) + (-qJ(1) - t129) * V_base(4) + t158;
t63 = t99 * V_base(4) + (-t100 - t132) * V_base(5) + t169;
t62 = t112 * t134 + (-t91 + t189) * V_base(6) + t168;
t61 = -t112 * t135 + t92 * V_base(6) + t154;
t60 = -t134 * t92 + t135 * t91 + t156;
t59 = t103 * t82 + t113 * t134 - t128 * t72 + (-t101 + t189) * V_base(6) + t168;
t58 = t102 * V_base(6) - t104 * t82 - t113 * t135 + t128 * t73 + t154;
t57 = t101 * t135 - t102 * t134 - t103 * t73 + t104 * t72 + t156;
t1 = m(1) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t135 * (t157 * t146 + t153 * t148) / 0.2e1 + t134 * (t153 * t146 - t157 * t148) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t104 * ((t107 * t69 + t108 * t71 + t67 * t182) * t104 + (t107 * t68 + t108 * t70 + t66 * t182) * t103 + (t107 * t80 + t108 * t81 + t79 * t182) * t128) / 0.2e1 + t103 * ((t105 * t69 + t106 * t71 + t67 * t183) * t104 + (t105 * t68 + t106 * t70 + t66 * t183) * t103 + (t105 * t80 + t106 * t81 + t79 * t183) * t128) / 0.2e1 + t128 * ((-t66 * t103 - t67 * t104 - t79 * t128) * t141 + ((-t150 * t69 + t151 * t71) * t104 + (-t150 * t68 + t151 * t70) * t103 + (-t150 * t80 + t151 * t81) * t128) * t140) / 0.2e1 + (t155 * t146 + t152 * t148 + t193 * V_base(6) + (-t123 * t146 + t126 * t148 + Icges(1,4)) * V_base(5) + (-t146 * t124 + t148 * t127 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t152 * t146 - t155 * t148 + t192 * V_base(6) + (t148 * t123 + t146 * t126 + Icges(1,2)) * V_base(5) + (t124 * t148 + t127 * t146 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t140 * t88 + t141 * t86) * t135 + (t140 * t87 + t141 * t85) * t134 + (t145 * t97 + t147 * t95 + t192) * V_base(5) + (t145 * t98 + t147 * t96 + t193) * V_base(4) + (t141 * t110 + t140 * t111 + t147 * t122 + t145 * t125 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1;
T = t1;
