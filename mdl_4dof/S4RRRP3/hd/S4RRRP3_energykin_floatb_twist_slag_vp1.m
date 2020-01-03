% Calculate kinetic energy for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:59
% EndTime: 2019-12-31 17:14:00
% DurationCPUTime: 1.27s
% Computational Cost: add. (654->160), mult. (643->211), div. (0->0), fcn. (473->6), ass. (0->84)
t183 = Icges(4,4) - Icges(5,5);
t182 = Icges(4,1) + Icges(5,1);
t181 = Icges(4,2) + Icges(5,3);
t123 = cos(qJ(3));
t180 = t183 * t123;
t121 = sin(qJ(3));
t179 = t183 * t121;
t178 = Icges(5,4) + Icges(4,5);
t177 = Icges(4,6) - Icges(5,6);
t176 = t181 * t121 - t180;
t175 = t182 * t123 - t179;
t174 = rSges(5,1) + pkin(3);
t173 = rSges(5,3) + qJ(4);
t120 = qJ(1) + qJ(2);
t115 = sin(t120);
t116 = cos(t120);
t172 = t176 * t115 + t177 * t116;
t171 = -t177 * t115 + t176 * t116;
t170 = t175 * t115 - t178 * t116;
t169 = t178 * t115 + t175 * t116;
t168 = -t181 * t123 - t179;
t167 = Icges(5,2) + Icges(4,3);
t166 = t182 * t121 + t180;
t165 = -t177 * t121 + t178 * t123;
t164 = t173 * t121 + t174 * t123;
t114 = V_base(6) + qJD(1);
t113 = qJD(2) + t114;
t91 = -qJD(3) * t116 + V_base(5);
t92 = qJD(3) * t115 + V_base(4);
t163 = (t171 * t121 + t169 * t123) * t92 + (t172 * t121 + t170 * t123) * t91 + (t168 * t121 + t166 * t123) * t113;
t162 = (t167 * t115 + t165 * t116) * t92 + (t165 * t115 - t167 * t116) * t91 + (t178 * t121 + t177 * t123) * t113;
t159 = -pkin(4) - pkin(5);
t122 = sin(qJ(1));
t157 = pkin(1) * t122;
t124 = cos(qJ(1));
t156 = pkin(1) * t124;
t155 = -rSges(5,2) * t116 + t164 * t115;
t154 = rSges(5,2) * t115 + t164 * t116;
t153 = Icges(2,4) * t122;
t152 = Icges(3,4) * t115;
t147 = t174 * t121 - t173 * t123;
t146 = qJD(4) * t121;
t145 = t114 * t156 + V_base(2);
t144 = V_base(4) * t157 + V_base(3);
t143 = V_base(5) * pkin(4) + V_base(1);
t140 = rSges(4,1) * t123 - rSges(4,2) * t121;
t131 = V_base(5) * pkin(5) - t114 * t157 + t143;
t87 = pkin(2) * t116 + pkin(6) * t115;
t128 = t113 * t87 + t159 * V_base(4) + t145;
t86 = pkin(2) * t115 - pkin(6) * t116;
t127 = V_base(4) * t86 + (-t87 - t156) * V_base(5) + t144;
t117 = Icges(2,4) * t124;
t112 = Icges(3,4) * t116;
t109 = rSges(2,1) * t124 - t122 * rSges(2,2);
t108 = t122 * rSges(2,1) + rSges(2,2) * t124;
t107 = rSges(4,1) * t121 + rSges(4,2) * t123;
t104 = Icges(2,1) * t124 - t153;
t103 = Icges(2,1) * t122 + t117;
t100 = -Icges(2,2) * t122 + t117;
t99 = Icges(2,2) * t124 + t153;
t90 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t89 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t88 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t85 = rSges(3,1) * t116 - rSges(3,2) * t115;
t84 = rSges(3,1) * t115 + rSges(3,2) * t116;
t83 = Icges(3,1) * t116 - t152;
t82 = Icges(3,1) * t115 + t112;
t81 = -Icges(3,2) * t115 + t112;
t80 = Icges(3,2) * t116 + t152;
t73 = rSges(4,3) * t115 + t116 * t140;
t71 = -rSges(4,3) * t116 + t115 * t140;
t57 = V_base(5) * rSges(2,3) - t108 * t114 + t143;
t56 = t109 * t114 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t55 = t108 * V_base(4) - t109 * V_base(5) + V_base(3);
t54 = V_base(5) * rSges(3,3) - t113 * t84 + t131;
t53 = t113 * t85 + (-rSges(3,3) + t159) * V_base(4) + t145;
t52 = V_base(4) * t84 + (-t85 - t156) * V_base(5) + t144;
t51 = t107 * t91 + (-t71 - t86) * t113 + t131;
t50 = -t107 * t92 + t113 * t73 + t128;
t49 = t92 * t71 - t91 * t73 + t127;
t48 = t116 * t146 + t147 * t91 + (-t86 - t155) * t113 + t131;
t47 = t113 * t154 + t115 * t146 - t147 * t92 + t128;
t46 = -qJD(4) * t123 - t154 * t91 + t155 * t92 + t127;
t1 = m(1) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(2) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(4) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(5) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + (t163 * t115 - t162 * t116) * t91 / 0.2e1 + (t162 * t115 + t163 * t116) * t92 / 0.2e1 + ((t169 * t121 - t171 * t123) * t92 + (t170 * t121 - t172 * t123) * t91 + (t166 * t121 - t168 * t123 + Icges(3,3)) * t113) * t113 / 0.2e1 + ((t103 * t124 - t115 * t80 + t116 * t82 - t122 * t99 + Icges(1,4)) * V_base(5) + (-t122 * t100 + t124 * t104 - t115 * t81 + t116 * t83 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t122 * t103 + t115 * t82 + t116 * t80 + t124 * t99 + Icges(1,2)) * V_base(5) + (t100 * t124 + t122 * t104 + t115 * t83 + t116 * t81 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t113 * (Icges(3,5) * t116 - Icges(3,6) * t115) + V_base(5) * t113 * (Icges(3,5) * t115 + Icges(3,6) * t116) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t122 + Icges(2,6) * t124) * V_base(5) + (Icges(2,5) * t124 - Icges(2,6) * t122) * V_base(4) + Icges(2,3) * t114 / 0.2e1) * t114;
T = t1;
