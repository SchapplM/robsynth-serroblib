% Calculate kinetic energy for
% S4RRRP5
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:35
% EndTime: 2019-12-31 17:16:36
% DurationCPUTime: 1.52s
% Computational Cost: add. (686->177), mult. (845->249), div. (0->0), fcn. (677->6), ass. (0->96)
t203 = Icges(4,4) - Icges(5,5);
t202 = Icges(4,1) + Icges(5,1);
t201 = Icges(4,2) + Icges(5,3);
t131 = qJ(2) + qJ(3);
t128 = cos(t131);
t200 = t203 * t128;
t127 = sin(t131);
t199 = t203 * t127;
t198 = Icges(5,4) + Icges(4,5);
t197 = Icges(4,6) - Icges(5,6);
t196 = t201 * t127 - t200;
t195 = t202 * t128 - t199;
t194 = rSges(5,1) + pkin(3);
t193 = rSges(5,3) + qJ(4);
t133 = sin(qJ(1));
t135 = cos(qJ(1));
t192 = t196 * t133 + t197 * t135;
t191 = -t197 * t133 + t196 * t135;
t190 = t195 * t133 - t198 * t135;
t189 = t198 * t133 + t195 * t135;
t188 = -t201 * t128 - t199;
t187 = t202 * t127 + t200;
t186 = Icges(5,2) + Icges(4,3);
t185 = -t197 * t127 + t198 * t128;
t184 = t193 * t127 + t194 * t128;
t101 = V_base(5) + (-qJD(2) - qJD(3)) * t135;
t123 = qJD(2) * t133 + V_base(4);
t102 = qJD(3) * t133 + t123;
t125 = V_base(6) + qJD(1);
t183 = (t188 * t127 + t187 * t128) * t125 + (t191 * t127 + t189 * t128) * t102 + (t192 * t127 + t190 * t128) * t101;
t182 = (t198 * t127 + t197 * t128) * t125 + (t186 * t133 + t185 * t135) * t102 + (t185 * t133 - t186 * t135) * t101;
t132 = sin(qJ(2));
t178 = pkin(2) * t132;
t134 = cos(qJ(2));
t177 = pkin(2) * t134;
t175 = -rSges(5,2) * t135 + t184 * t133;
t174 = rSges(5,2) * t133 + t184 * t135;
t173 = t194 * t127 - t193 * t128;
t120 = t133 * pkin(1) - t135 * pkin(5);
t63 = -pkin(6) * t135 + t177 * t133;
t172 = -t120 - t63;
t171 = Icges(2,4) * t133;
t170 = Icges(3,4) * t132;
t169 = Icges(3,4) * t134;
t164 = qJD(4) * t127;
t163 = V_base(5) * pkin(4) + V_base(1);
t122 = -qJD(2) * t135 + V_base(5);
t160 = t122 * t178 + t163;
t159 = rSges(3,1) * t134 - rSges(3,2) * t132;
t158 = rSges(4,1) * t128 - rSges(4,2) * t127;
t155 = Icges(3,1) * t134 - t170;
t152 = -Icges(3,2) * t132 + t169;
t149 = Icges(3,5) * t134 - Icges(3,6) * t132;
t121 = t135 * pkin(1) + t133 * pkin(5);
t146 = -V_base(4) * pkin(4) + t125 * t121 + V_base(2);
t145 = V_base(4) * t120 - t121 * V_base(5) + V_base(3);
t142 = (Icges(3,5) * t132 + Icges(3,6) * t134) * t125 + (-Icges(3,3) * t135 + t149 * t133) * t122 + (Icges(3,3) * t133 + t149 * t135) * t123;
t64 = pkin(6) * t133 + t177 * t135;
t141 = -t122 * t64 + t123 * t63 + t145;
t140 = -t123 * t178 + t125 * t64 + t146;
t111 = Icges(3,2) * t134 + t170;
t114 = Icges(3,1) * t132 + t169;
t83 = -Icges(3,6) * t135 + t152 * t133;
t84 = Icges(3,6) * t133 + t152 * t135;
t85 = -Icges(3,5) * t135 + t155 * t133;
t86 = Icges(3,5) * t133 + t155 * t135;
t137 = (-t132 * t84 + t134 * t86) * t123 + (-t132 * t83 + t134 * t85) * t122 + (-t111 * t132 + t114 * t134) * t125;
t129 = Icges(2,4) * t135;
t119 = rSges(2,1) * t135 - rSges(2,2) * t133;
t118 = rSges(2,1) * t133 + rSges(2,2) * t135;
t117 = rSges(3,1) * t132 + rSges(3,2) * t134;
t116 = Icges(2,1) * t135 - t171;
t115 = Icges(2,1) * t133 + t129;
t113 = -Icges(2,2) * t133 + t129;
t112 = Icges(2,2) * t135 + t171;
t107 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t106 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t105 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t100 = rSges(4,1) * t127 + rSges(4,2) * t128;
t88 = rSges(3,3) * t133 + t159 * t135;
t87 = -rSges(3,3) * t135 + t159 * t133;
t80 = rSges(4,3) * t133 + t158 * t135;
t78 = -rSges(4,3) * t135 + t158 * t133;
t62 = V_base(5) * rSges(2,3) - t118 * t125 + t163;
t61 = t119 * t125 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t60 = t118 * V_base(4) - t119 * V_base(5) + V_base(3);
t57 = t117 * t122 + (-t120 - t87) * t125 + t163;
t56 = -t117 * t123 + t125 * t88 + t146;
t55 = -t122 * t88 + t123 * t87 + t145;
t54 = t100 * t101 + (-t78 + t172) * t125 + t160;
t53 = -t100 * t102 + t125 * t80 + t140;
t52 = -t101 * t80 + t102 * t78 + t141;
t51 = t135 * t164 + t173 * t101 + (t172 - t175) * t125 + t160;
t50 = -t173 * t102 + t174 * t125 + t133 * t164 + t140;
t49 = -qJD(4) * t128 - t174 * t101 + t175 * t102 + t141;
t1 = m(1) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(2) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(3) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + t123 * (t142 * t133 + t137 * t135) / 0.2e1 + t122 * (t137 * t133 - t142 * t135) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + (t183 * t133 - t182 * t135) * t101 / 0.2e1 + (t182 * t133 + t183 * t135) * t102 / 0.2e1 + ((-t112 * t133 + t115 * t135 + Icges(1,4)) * V_base(5) + (-t113 * t133 + t116 * t135 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t112 * t135 + t115 * t133 + Icges(1,2)) * V_base(5) + (t113 * t135 + t116 * t133 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t132 * t86 + t134 * t84) * t123 + (t132 * t85 + t134 * t83) * t122 + (t189 * t127 - t191 * t128) * t102 + (t190 * t127 - t192 * t128) * t101 + (t134 * t111 + t132 * t114 + t187 * t127 - t188 * t128 + Icges(2,3)) * t125) * t125 / 0.2e1 + V_base(4) * t125 * (Icges(2,5) * t135 - Icges(2,6) * t133) + V_base(5) * t125 * (Icges(2,5) * t133 + Icges(2,6) * t135) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
