% Calculate kinetic energy for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:07
% EndTime: 2019-12-31 17:15:08
% DurationCPUTime: 1.49s
% Computational Cost: add. (706->176), mult. (852->247), div. (0->0), fcn. (684->6), ass. (0->95)
t201 = Icges(4,4) + Icges(5,4);
t200 = Icges(4,1) + Icges(5,1);
t199 = Icges(4,2) + Icges(5,2);
t129 = qJ(2) + qJ(3);
t124 = cos(t129);
t198 = t201 * t124;
t123 = sin(t129);
t197 = t201 * t123;
t196 = Icges(4,5) + Icges(5,5);
t195 = Icges(4,6) + Icges(5,6);
t194 = -t199 * t123 + t198;
t193 = t200 * t124 - t197;
t192 = rSges(5,1) + pkin(3);
t131 = sin(qJ(1));
t133 = cos(qJ(1));
t191 = t194 * t131 - t195 * t133;
t190 = t195 * t131 + t194 * t133;
t189 = t193 * t131 - t196 * t133;
t188 = t196 * t131 + t193 * t133;
t187 = t199 * t124 + t197;
t186 = t200 * t123 + t198;
t185 = Icges(4,3) + Icges(5,3);
t184 = -t195 * t123 + t196 * t124;
t183 = rSges(5,3) + qJ(4);
t182 = -rSges(5,2) * t123 + t192 * t124;
t121 = V_base(6) + qJD(1);
t96 = V_base(5) + (-qJD(2) - qJD(3)) * t133;
t119 = qJD(2) * t131 + V_base(4);
t97 = qJD(3) * t131 + t119;
t181 = (-t190 * t123 + t188 * t124) * t97 + (-t191 * t123 + t189 * t124) * t96 + (-t187 * t123 + t186 * t124) * t121;
t180 = (t185 * t131 + t184 * t133) * t97 + (t184 * t131 - t185 * t133) * t96 + (t196 * t123 + t195 * t124) * t121;
t130 = sin(qJ(2));
t176 = pkin(2) * t130;
t132 = cos(qJ(2));
t175 = t132 * pkin(2);
t173 = t182 * t131 - t183 * t133;
t172 = t183 * t131 + t182 * t133;
t116 = t131 * pkin(1) - t133 * pkin(5);
t61 = -pkin(6) * t133 + t175 * t131;
t171 = -t116 - t61;
t170 = Icges(2,4) * t131;
t169 = Icges(3,4) * t130;
t168 = Icges(3,4) * t132;
t161 = V_base(5) * pkin(4) + V_base(1);
t158 = rSges(5,2) * t124 + t192 * t123;
t118 = -qJD(2) * t133 + V_base(5);
t157 = t118 * t176 + t161;
t156 = rSges(3,1) * t132 - rSges(3,2) * t130;
t155 = rSges(4,1) * t124 - rSges(4,2) * t123;
t153 = Icges(3,1) * t132 - t169;
t150 = -Icges(3,2) * t130 + t168;
t147 = Icges(3,5) * t132 - Icges(3,6) * t130;
t117 = t133 * pkin(1) + t131 * pkin(5);
t144 = -V_base(4) * pkin(4) + t121 * t117 + V_base(2);
t143 = V_base(4) * t116 - t117 * V_base(5) + V_base(3);
t140 = (Icges(3,5) * t130 + Icges(3,6) * t132) * t121 + (-Icges(3,3) * t133 + t147 * t131) * t118 + (Icges(3,3) * t131 + t147 * t133) * t119;
t62 = pkin(6) * t131 + t175 * t133;
t139 = -t118 * t62 + t119 * t61 + t143;
t138 = -t119 * t176 + t121 * t62 + t144;
t107 = Icges(3,2) * t132 + t169;
t110 = Icges(3,1) * t130 + t168;
t81 = -Icges(3,6) * t133 + t150 * t131;
t82 = Icges(3,6) * t131 + t150 * t133;
t83 = -Icges(3,5) * t133 + t153 * t131;
t84 = Icges(3,5) * t131 + t153 * t133;
t135 = (-t130 * t82 + t132 * t84) * t119 + (-t130 * t81 + t132 * t83) * t118 + (-t107 * t130 + t110 * t132) * t121;
t125 = Icges(2,4) * t133;
t115 = rSges(2,1) * t133 - rSges(2,2) * t131;
t114 = rSges(2,1) * t131 + rSges(2,2) * t133;
t113 = rSges(3,1) * t130 + rSges(3,2) * t132;
t112 = Icges(2,1) * t133 - t170;
t111 = Icges(2,1) * t131 + t125;
t109 = -Icges(2,2) * t131 + t125;
t108 = Icges(2,2) * t133 + t170;
t103 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t102 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t101 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t95 = rSges(4,1) * t123 + rSges(4,2) * t124;
t86 = rSges(3,3) * t131 + t156 * t133;
t85 = -rSges(3,3) * t133 + t156 * t131;
t78 = rSges(4,3) * t131 + t155 * t133;
t76 = -rSges(4,3) * t133 + t155 * t131;
t60 = V_base(5) * rSges(2,3) - t114 * t121 + t161;
t59 = t115 * t121 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t58 = t114 * V_base(4) - t115 * V_base(5) + V_base(3);
t53 = t113 * t118 + (-t116 - t85) * t121 + t161;
t52 = -t113 * t119 + t121 * t86 + t144;
t51 = -t118 * t86 + t119 * t85 + t143;
t50 = t95 * t96 + (-t76 + t171) * t121 + t157;
t49 = t121 * t78 - t95 * t97 + t138;
t48 = t76 * t97 - t78 * t96 + t139;
t47 = qJD(4) * t131 + t158 * t96 + (t171 - t173) * t121 + t157;
t46 = -qJD(4) * t133 + t172 * t121 - t158 * t97 + t138;
t45 = -t172 * t96 + t173 * t97 + t139;
t1 = m(1) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(2) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(3) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + t119 * (t140 * t131 + t135 * t133) / 0.2e1 + t118 * (t135 * t131 - t140 * t133) / 0.2e1 + m(4) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + m(5) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + (t181 * t131 - t180 * t133) * t96 / 0.2e1 + (t180 * t131 + t181 * t133) * t97 / 0.2e1 + ((-t108 * t131 + t111 * t133 + Icges(1,4)) * V_base(5) + (-t109 * t131 + t112 * t133 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t108 * t133 + t111 * t131 + Icges(1,2)) * V_base(5) + (t109 * t133 + t112 * t131 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t130 * t84 + t132 * t82) * t119 + (t130 * t83 + t132 * t81) * t118 + (t188 * t123 + t190 * t124) * t97 + (t189 * t123 + t191 * t124) * t96 + (t107 * t132 + t110 * t130 + t186 * t123 + t187 * t124 + Icges(2,3)) * t121) * t121 / 0.2e1 + t121 * V_base(4) * (Icges(2,5) * t133 - Icges(2,6) * t131) + V_base(5) * t121 * (Icges(2,5) * t131 + Icges(2,6) * t133) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
