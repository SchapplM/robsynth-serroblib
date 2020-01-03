% Calculate kinetic energy for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:35
% EndTime: 2019-12-31 16:43:37
% DurationCPUTime: 1.23s
% Computational Cost: add. (633->161), mult. (643->205), div. (0->0), fcn. (473->6), ass. (0->84)
t185 = Icges(4,4) - Icges(5,5);
t184 = Icges(4,1) + Icges(5,1);
t183 = Icges(4,2) + Icges(5,3);
t122 = cos(qJ(3));
t182 = t185 * t122;
t120 = sin(qJ(3));
t181 = t185 * t120;
t180 = Icges(5,4) + Icges(4,5);
t179 = Icges(4,6) - Icges(5,6);
t178 = t183 * t120 - t182;
t177 = t184 * t122 - t181;
t176 = rSges(5,1) + pkin(3);
t175 = rSges(5,3) + qJ(4);
t119 = qJ(1) + pkin(6);
t113 = sin(t119);
t114 = cos(t119);
t174 = t178 * t113 + t179 * t114;
t173 = -t179 * t113 + t178 * t114;
t172 = t177 * t113 - t180 * t114;
t171 = t180 * t113 + t177 * t114;
t170 = -t183 * t122 - t181;
t169 = Icges(5,2) + Icges(4,3);
t168 = t184 * t120 + t182;
t167 = -t179 * t120 + t180 * t122;
t166 = t175 * t120 + t176 * t122;
t115 = V_base(6) + qJD(1);
t91 = -qJD(3) * t114 + V_base(5);
t92 = qJD(3) * t113 + V_base(4);
t163 = (t173 * t120 + t171 * t122) * t92 + (t174 * t120 + t172 * t122) * t91 + (t170 * t120 + t168 * t122) * t115;
t162 = (t169 * t113 + t167 * t114) * t92 + (t167 * t113 - t169 * t114) * t91 + (t180 * t120 + t179 * t122) * t115;
t121 = sin(qJ(1));
t158 = pkin(1) * t121;
t123 = cos(qJ(1));
t157 = pkin(1) * t123;
t156 = -pkin(4) - qJ(2);
t155 = -rSges(5,2) * t114 + t113 * t166;
t154 = rSges(5,2) * t113 + t114 * t166;
t153 = Icges(2,4) * t121;
t152 = Icges(3,4) * t113;
t147 = t176 * t120 - t175 * t122;
t146 = qJD(4) * t120;
t145 = t115 * t157 + V_base(2);
t144 = V_base(5) * pkin(4) + V_base(1);
t86 = pkin(2) * t113 - pkin(5) * t114;
t141 = -t86 - t158;
t140 = V_base(5) * qJ(2) + t144;
t139 = V_base(4) * t158 + qJD(2) + V_base(3);
t138 = rSges(4,1) * t122 - rSges(4,2) * t120;
t87 = pkin(2) * t114 + pkin(5) * t113;
t127 = t115 * t87 + t156 * V_base(4) + t145;
t126 = V_base(4) * t86 + (-t87 - t157) * V_base(5) + t139;
t117 = Icges(2,4) * t123;
t112 = Icges(3,4) * t114;
t109 = rSges(2,1) * t123 - t121 * rSges(2,2);
t108 = t121 * rSges(2,1) + rSges(2,2) * t123;
t107 = rSges(4,1) * t120 + rSges(4,2) * t122;
t104 = Icges(2,1) * t123 - t153;
t103 = Icges(2,1) * t121 + t117;
t100 = -Icges(2,2) * t121 + t117;
t99 = Icges(2,2) * t123 + t153;
t90 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t89 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t88 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t85 = rSges(3,1) * t114 - rSges(3,2) * t113;
t84 = rSges(3,1) * t113 + rSges(3,2) * t114;
t83 = Icges(3,1) * t114 - t152;
t82 = Icges(3,1) * t113 + t112;
t81 = -Icges(3,2) * t113 + t112;
t80 = Icges(3,2) * t114 + t152;
t73 = rSges(4,3) * t113 + t114 * t138;
t71 = -rSges(4,3) * t114 + t113 * t138;
t69 = V_base(5) * rSges(2,3) - t108 * t115 + t144;
t68 = t109 * t115 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t55 = t108 * V_base(4) - t109 * V_base(5) + V_base(3);
t54 = V_base(5) * rSges(3,3) + (-t84 - t158) * t115 + t140;
t53 = t115 * t85 + (-rSges(3,3) + t156) * V_base(4) + t145;
t52 = V_base(4) * t84 + (-t85 - t157) * V_base(5) + t139;
t51 = t107 * t91 + (t141 - t71) * t115 + t140;
t50 = -t107 * t92 + t115 * t73 + t127;
t49 = t92 * t71 - t91 * t73 + t126;
t48 = t114 * t146 + t147 * t91 + (t141 - t155) * t115 + t140;
t47 = t113 * t146 + t115 * t154 - t147 * t92 + t127;
t46 = -qJD(4) * t122 - t154 * t91 + t155 * t92 + t126;
t1 = m(1) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(2) * (t55 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(3) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(4) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(5) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + (t163 * t113 - t162 * t114) * t91 / 0.2e1 + (t162 * t113 + t163 * t114) * t92 / 0.2e1 + ((t103 * t123 - t113 * t80 + t114 * t82 - t121 * t99 + Icges(1,4)) * V_base(5) + (-t121 * t100 + t123 * t104 - t113 * t81 + t114 * t83 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t121 * t103 + t113 * t82 + t114 * t80 + t123 * t99 + Icges(1,2)) * V_base(5) + (t100 * t123 + t121 * t104 + t113 * t83 + t114 * t81 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t171 * t120 - t173 * t122) * t92 + (t172 * t120 - t174 * t122) * t91 + (t168 * t120 - t170 * t122 + Icges(2,3) + Icges(3,3)) * t115) * t115 / 0.2e1 + t115 * V_base(5) * (Icges(2,5) * t121 + Icges(3,5) * t113 + Icges(2,6) * t123 + Icges(3,6) * t114) + t115 * V_base(4) * (Icges(2,5) * t123 + Icges(3,5) * t114 - Icges(2,6) * t121 - Icges(3,6) * t113) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
