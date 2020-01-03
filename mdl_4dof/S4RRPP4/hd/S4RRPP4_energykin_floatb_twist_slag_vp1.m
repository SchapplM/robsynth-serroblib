% Calculate kinetic energy for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:50
% EndTime: 2019-12-31 16:58:51
% DurationCPUTime: 1.27s
% Computational Cost: add. (462->148), mult. (814->197), div. (0->0), fcn. (646->4), ass. (0->78)
t193 = Icges(3,4) - Icges(5,4) - Icges(4,5);
t192 = Icges(3,1) + Icges(4,1) + Icges(5,1);
t191 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t127 = cos(qJ(2));
t190 = t193 * t127;
t125 = sin(qJ(2));
t189 = t193 * t125;
t188 = Icges(4,4) + Icges(3,5) - Icges(5,5);
t187 = Icges(3,6) - Icges(4,6) + Icges(5,6);
t186 = t191 * t125 - t190;
t185 = t192 * t127 - t189;
t184 = rSges(5,1) + pkin(3);
t126 = sin(qJ(1));
t128 = cos(qJ(1));
t183 = t186 * t126 + t187 * t128;
t182 = -t187 * t126 + t186 * t128;
t181 = t185 * t126 - t188 * t128;
t180 = t188 * t126 + t185 * t128;
t179 = -t191 * t127 - t189;
t178 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t177 = t192 * t125 + t190;
t176 = -t187 * t125 + t188 * t127;
t175 = rSges(5,3) + qJ(4);
t174 = rSges(5,2) * t125 + t184 * t127;
t118 = -qJD(2) * t128 + V_base(5);
t119 = qJD(2) * t126 + V_base(4);
t122 = V_base(6) + qJD(1);
t173 = (t179 * t125 + t177 * t127) * t122 + (t182 * t125 + t180 * t127) * t119 + (t183 * t125 + t181 * t127) * t118;
t172 = (t188 * t125 + t187 * t127) * t122 + (t178 * t126 + t176 * t128) * t119 + (t176 * t126 - t178 * t128) * t118;
t167 = t174 * t126 + t128 * t175;
t166 = -t126 * t175 + t174 * t128;
t116 = t126 * pkin(1) - pkin(5) * t128;
t148 = pkin(2) * t127 + qJ(3) * t125;
t86 = t148 * t126;
t165 = -t116 - t86;
t164 = Icges(2,4) * t126;
t157 = qJD(3) * t125;
t156 = V_base(5) * pkin(4) + V_base(1);
t153 = -rSges(5,2) * t127 + t184 * t125;
t110 = pkin(2) * t125 - qJ(3) * t127;
t152 = t118 * t110 + t128 * t157 + t156;
t151 = rSges(3,1) * t127 - rSges(3,2) * t125;
t150 = rSges(4,1) * t127 + rSges(4,3) * t125;
t117 = pkin(1) * t128 + t126 * pkin(5);
t138 = -V_base(4) * pkin(4) + t122 * t117 + V_base(2);
t137 = V_base(4) * t116 - t117 * V_base(5) + V_base(3);
t87 = t148 * t128;
t133 = t122 * t87 + t126 * t157 + t138;
t132 = -qJD(3) * t127 + t119 * t86 + t137;
t123 = Icges(2,4) * t128;
t115 = rSges(2,1) * t128 - t126 * rSges(2,2);
t114 = t126 * rSges(2,1) + rSges(2,2) * t128;
t113 = rSges(3,1) * t125 + rSges(3,2) * t127;
t112 = rSges(4,1) * t125 - rSges(4,3) * t127;
t109 = Icges(2,1) * t128 - t164;
t108 = Icges(2,1) * t126 + t123;
t104 = -Icges(2,2) * t126 + t123;
t103 = Icges(2,2) * t128 + t164;
t94 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t93 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t92 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t84 = t126 * rSges(3,3) + t151 * t128;
t83 = t126 * rSges(4,2) + t150 * t128;
t81 = -rSges(3,3) * t128 + t151 * t126;
t80 = -rSges(4,2) * t128 + t150 * t126;
t58 = V_base(5) * rSges(2,3) - t114 * t122 + t156;
t57 = t115 * t122 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t56 = t114 * V_base(4) - t115 * V_base(5) + V_base(3);
t55 = t113 * t118 + (-t116 - t81) * t122 + t156;
t54 = -t113 * t119 + t122 * t84 + t138;
t53 = -t118 * t84 + t119 * t81 + t137;
t52 = t112 * t118 + (-t80 + t165) * t122 + t152;
t51 = t122 * t83 + (-t110 - t112) * t119 + t133;
t50 = t119 * t80 + (-t83 - t87) * t118 + t132;
t49 = -qJD(4) * t126 + t153 * t118 + (t165 - t167) * t122 + t152;
t48 = qJD(4) * t128 + t166 * t122 + (-t110 - t153) * t119 + t133;
t47 = t167 * t119 + (-t87 - t166) * t118 + t132;
t1 = m(1) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(2) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(3) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + m(4) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(5) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + ((-t126 * t103 + t108 * t128 + Icges(1,4)) * V_base(5) + (-t126 * t104 + t109 * t128 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t103 * t128 + t126 * t108 + Icges(1,2)) * V_base(5) + (t104 * t128 + t126 * t109 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t173 * t126 - t172 * t128) * t118 / 0.2e1 + (t172 * t126 + t173 * t128) * t119 / 0.2e1 + ((t180 * t125 - t182 * t127) * t119 + (t181 * t125 - t183 * t127) * t118 + (t177 * t125 - t179 * t127 + Icges(2,3)) * t122) * t122 / 0.2e1 + V_base(4) * t122 * (Icges(2,5) * t128 - Icges(2,6) * t126) + V_base(5) * t122 * (Icges(2,5) * t126 + Icges(2,6) * t128) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
