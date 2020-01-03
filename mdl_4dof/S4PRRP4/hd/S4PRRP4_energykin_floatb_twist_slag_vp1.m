% Calculate kinetic energy for
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:43
% EndTime: 2019-12-31 16:27:44
% DurationCPUTime: 1.26s
% Computational Cost: add. (622->160), mult. (643->208), div. (0->0), fcn. (473->6), ass. (0->84)
t186 = Icges(4,4) - Icges(5,5);
t185 = Icges(4,1) + Icges(5,1);
t184 = Icges(4,2) + Icges(5,3);
t123 = cos(qJ(3));
t183 = t186 * t123;
t122 = sin(qJ(3));
t182 = t186 * t122;
t181 = Icges(5,4) + Icges(4,5);
t180 = Icges(4,6) - Icges(5,6);
t179 = t184 * t122 - t183;
t178 = t185 * t123 - t182;
t177 = rSges(5,1) + pkin(3);
t176 = rSges(5,3) + qJ(4);
t119 = pkin(6) + qJ(2);
t113 = sin(t119);
t114 = cos(t119);
t175 = t179 * t113 + t180 * t114;
t174 = -t180 * t113 + t179 * t114;
t173 = t178 * t113 - t181 * t114;
t172 = t181 * t113 + t178 * t114;
t171 = Icges(5,2) + Icges(4,3);
t170 = -t184 * t123 - t182;
t169 = t185 * t122 + t183;
t168 = -t180 * t122 + t181 * t123;
t167 = t176 * t122 + t177 * t123;
t116 = V_base(6) + qJD(2);
t97 = -qJD(3) * t114 + V_base(5);
t98 = qJD(3) * t113 + V_base(4);
t164 = (t174 * t122 + t172 * t123) * t98 + (t175 * t122 + t173 * t123) * t97 + (t170 * t122 + t169 * t123) * t116;
t163 = (t171 * t113 + t168 * t114) * t98 + (t168 * t113 - t171 * t114) * t97 + (t181 * t122 + t180 * t123) * t116;
t121 = cos(pkin(6));
t159 = pkin(1) * t121;
t158 = -pkin(4) - qJ(1);
t157 = -t114 * rSges(5,2) + t167 * t113;
t156 = t113 * rSges(5,2) + t167 * t114;
t120 = sin(pkin(6));
t155 = Icges(2,4) * t120;
t154 = Icges(3,4) * t113;
t149 = t177 * t122 - t176 * t123;
t148 = qJD(4) * t122;
t141 = pkin(1) * V_base(6);
t147 = t121 * t141 + V_base(2);
t146 = V_base(5) * qJ(1) + V_base(1);
t142 = qJD(1) + V_base(3);
t140 = V_base(4) * t120 * pkin(1) + t142;
t139 = rSges(4,1) * t123 - rSges(4,2) * t122;
t128 = V_base(5) * pkin(4) - t120 * t141 + t146;
t87 = pkin(2) * t114 + pkin(5) * t113;
t127 = t116 * t87 + t158 * V_base(4) + t147;
t86 = pkin(2) * t113 - pkin(5) * t114;
t126 = V_base(4) * t86 + (-t87 - t159) * V_base(5) + t140;
t115 = Icges(2,4) * t121;
t112 = Icges(3,4) * t114;
t109 = t122 * rSges(4,1) + rSges(4,2) * t123;
t100 = rSges(2,1) * t121 - rSges(2,2) * t120;
t99 = rSges(2,1) * t120 + rSges(2,2) * t121;
t96 = Icges(2,1) * t121 - t155;
t95 = Icges(2,1) * t120 + t115;
t94 = -Icges(2,2) * t120 + t115;
t93 = Icges(2,2) * t121 + t155;
t90 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t89 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t88 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t85 = rSges(3,1) * t114 - rSges(3,2) * t113;
t84 = rSges(3,1) * t113 + rSges(3,2) * t114;
t83 = Icges(3,1) * t114 - t154;
t82 = Icges(3,1) * t113 + t112;
t81 = -Icges(3,2) * t113 + t112;
t80 = Icges(3,2) * t114 + t154;
t73 = V_base(5) * rSges(2,3) - t99 * V_base(6) + t146;
t72 = t100 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t71 = t113 * rSges(4,3) + t139 * t114;
t69 = -t114 * rSges(4,3) + t139 * t113;
t55 = -t100 * V_base(5) + t99 * V_base(4) + t142;
t54 = V_base(5) * rSges(3,3) - t116 * t84 + t128;
t53 = t116 * t85 + (-rSges(3,3) + t158) * V_base(4) + t147;
t52 = t84 * V_base(4) + (-t85 - t159) * V_base(5) + t140;
t51 = t109 * t97 + (-t69 - t86) * t116 + t128;
t50 = -t109 * t98 + t116 * t71 + t127;
t49 = t69 * t98 - t71 * t97 + t126;
t48 = t114 * t148 + t149 * t97 + (-t86 - t157) * t116 + t128;
t47 = t113 * t148 + t156 * t116 - t149 * t98 + t127;
t46 = -qJD(4) * t123 - t156 * t97 + t157 * t98 + t126;
t1 = m(1) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(2) * (t55 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(3) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(4) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(5) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + (t113 * t164 - t163 * t114) * t97 / 0.2e1 + (t113 * t163 + t114 * t164) * t98 / 0.2e1 + ((t172 * t122 - t174 * t123) * t98 + (t173 * t122 - t175 * t123) * t97 + (t169 * t122 - t170 * t123 + Icges(3,3)) * t116) * t116 / 0.2e1 + ((-t113 * t80 + t114 * t82 - t120 * t93 + t121 * t95 + Icges(1,4)) * V_base(5) + (-t113 * t81 + t114 * t83 - t120 * t94 + t121 * t96 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t113 * t82 + t114 * t80 + t120 * t95 + t121 * t93 + Icges(1,2)) * V_base(5) + (t113 * t83 + t114 * t81 + t120 * t96 + t121 * t94 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t116 * (Icges(3,5) * t114 - Icges(3,6) * t113) + V_base(5) * t116 * (Icges(3,5) * t113 + Icges(3,6) * t114) + ((Icges(2,5) * t120 + Icges(2,6) * t121 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t121 - Icges(2,6) * t120 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
