% Calculate kinetic energy for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:43
% EndTime: 2019-12-31 16:44:45
% DurationCPUTime: 1.56s
% Computational Cost: add. (644->180), mult. (803->241), div. (0->0), fcn. (635->6), ass. (0->97)
t202 = Icges(4,4) - Icges(5,5);
t201 = Icges(4,1) + Icges(5,1);
t200 = Icges(4,2) + Icges(5,3);
t129 = pkin(6) + qJ(3);
t123 = sin(t129);
t199 = t202 * t123;
t124 = cos(t129);
t198 = t202 * t124;
t197 = Icges(5,4) + Icges(4,5);
t196 = Icges(4,6) - Icges(5,6);
t195 = t200 * t123 - t198;
t194 = t201 * t124 - t199;
t193 = rSges(5,1) + pkin(3);
t192 = rSges(5,3) + qJ(4);
t133 = sin(qJ(1));
t134 = cos(qJ(1));
t191 = t195 * t133 + t196 * t134;
t190 = -t196 * t133 + t195 * t134;
t189 = t194 * t133 - t197 * t134;
t188 = t197 * t133 + t194 * t134;
t187 = -t200 * t124 - t199;
t186 = t201 * t123 + t198;
t185 = Icges(5,2) + Icges(4,3);
t184 = -t196 * t123 + t197 * t124;
t183 = t192 * t123 + t193 * t124;
t119 = -qJD(3) * t134 + V_base(5);
t120 = qJD(3) * t133 + V_base(4);
t125 = V_base(6) + qJD(1);
t182 = (t187 * t123 + t186 * t124) * t125 + (t190 * t123 + t188 * t124) * t120 + (t191 * t123 + t189 * t124) * t119;
t181 = (t197 * t123 + t196 * t124) * t125 + (t185 * t133 + t184 * t134) * t120 + (t184 * t133 - t185 * t134) * t119;
t130 = sin(pkin(6));
t177 = pkin(2) * t130;
t131 = cos(pkin(6));
t176 = t131 * pkin(2);
t175 = -t134 * rSges(5,2) + t183 * t133;
t174 = t133 * rSges(5,2) + t183 * t134;
t173 = t193 * t123 - t192 * t124;
t115 = t133 * pkin(1) - t134 * qJ(2);
t63 = -pkin(5) * t134 + t176 * t133;
t172 = -t115 - t63;
t171 = Icges(2,4) * t133;
t170 = Icges(3,4) * t130;
t169 = Icges(3,4) * t131;
t163 = qJD(4) * t123;
t162 = V_base(4) * t115 + V_base(3);
t161 = V_base(5) * pkin(4) + V_base(1);
t158 = qJD(2) * t133 + t161;
t157 = V_base(5) * t177 + t158;
t156 = rSges(3,1) * t131 - rSges(3,2) * t130;
t155 = rSges(4,1) * t124 - rSges(4,2) * t123;
t152 = Icges(3,1) * t131 - t170;
t149 = -Icges(3,2) * t130 + t169;
t146 = Icges(3,5) * t131 - Icges(3,6) * t130;
t117 = t134 * pkin(1) + t133 * qJ(2);
t143 = -qJD(2) * t134 + t125 * t117 + V_base(2);
t64 = pkin(5) * t133 + t176 * t134;
t140 = V_base(4) * t63 + (-t117 - t64) * V_base(5) + t162;
t139 = (Icges(3,5) * t130 + Icges(3,6) * t131) * t125 + (-Icges(3,3) * t134 + t146 * t133) * V_base(5) + (Icges(3,3) * t133 + t146 * t134) * V_base(4);
t138 = t125 * t64 + (-pkin(4) - t177) * V_base(4) + t143;
t106 = Icges(3,2) * t131 + t170;
t107 = Icges(3,1) * t130 + t169;
t83 = -Icges(3,6) * t134 + t149 * t133;
t84 = Icges(3,6) * t133 + t149 * t134;
t85 = -Icges(3,5) * t134 + t152 * t133;
t86 = Icges(3,5) * t133 + t152 * t134;
t135 = (-t130 * t84 + t131 * t86) * V_base(4) + (-t130 * t83 + t131 * t85) * V_base(5) + (-t106 * t130 + t107 * t131) * t125;
t127 = Icges(2,4) * t134;
t118 = t134 * rSges(2,1) - t133 * rSges(2,2);
t116 = t133 * rSges(2,1) + t134 * rSges(2,2);
t114 = Icges(2,1) * t134 - t171;
t113 = Icges(2,1) * t133 + t127;
t112 = -Icges(2,2) * t133 + t127;
t111 = Icges(2,2) * t134 + t171;
t110 = Icges(2,5) * t134 - Icges(2,6) * t133;
t109 = Icges(2,5) * t133 + Icges(2,6) * t134;
t108 = t130 * rSges(3,1) + t131 * rSges(3,2);
t104 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t103 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t102 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t100 = t123 * rSges(4,1) + t124 * rSges(4,2);
t88 = t133 * rSges(3,3) + t156 * t134;
t87 = -t134 * rSges(3,3) + t156 * t133;
t80 = t133 * rSges(4,3) + t155 * t134;
t78 = -t134 * rSges(4,3) + t155 * t133;
t62 = V_base(5) * rSges(2,3) - t125 * t116 + t161;
t61 = t125 * t118 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t60 = V_base(4) * t116 - V_base(5) * t118 + V_base(3);
t57 = V_base(5) * t108 + (-t115 - t87) * t125 + t158;
t56 = t125 * t88 + (-pkin(4) - t108) * V_base(4) + t143;
t55 = V_base(4) * t87 + (-t117 - t88) * V_base(5) + t162;
t54 = t119 * t100 + (-t78 + t172) * t125 + t157;
t53 = -t120 * t100 + t125 * t80 + t138;
t52 = -t119 * t80 + t120 * t78 + t140;
t51 = t134 * t163 + t173 * t119 + (t172 - t175) * t125 + t157;
t50 = -t173 * t120 + t174 * t125 + t133 * t163 + t138;
t49 = -qJD(4) * t124 - t174 * t119 + t175 * t120 + t140;
t1 = m(1) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(2) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(3) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + (t182 * t133 - t181 * t134) * t119 / 0.2e1 + (t181 * t133 + t182 * t134) * t120 / 0.2e1 + (t110 * t125 + t133 * t139 + t135 * t134 + (-t133 * t111 + t134 * t113 + Icges(1,4)) * V_base(5) + (-t133 * t112 + t134 * t114 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t109 * t125 + t135 * t133 - t134 * t139 + (t134 * t111 + t133 * t113 + Icges(1,2)) * V_base(5) + (t134 * t112 + t133 * t114 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t130 * t85 + t131 * t83 + t109) * V_base(5) + (t130 * t86 + t131 * t84 + t110) * V_base(4) + (t188 * t123 - t190 * t124) * t120 + (t189 * t123 - t191 * t124) * t119 + (t131 * t106 + t130 * t107 + t186 * t123 - t187 * t124 + Icges(2,3)) * t125) * t125 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
