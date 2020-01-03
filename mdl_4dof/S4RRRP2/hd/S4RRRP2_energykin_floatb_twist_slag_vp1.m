% Calculate kinetic energy for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:51
% EndTime: 2019-12-31 17:12:52
% DurationCPUTime: 1.23s
% Computational Cost: add. (675->159), mult. (642->209), div. (0->0), fcn. (472->6), ass. (0->83)
t178 = Icges(4,4) + Icges(5,4);
t177 = Icges(4,1) + Icges(5,1);
t176 = Icges(4,2) + Icges(5,2);
t118 = cos(qJ(3));
t175 = t178 * t118;
t116 = sin(qJ(3));
t174 = t178 * t116;
t173 = Icges(4,5) + Icges(5,5);
t172 = Icges(4,6) + Icges(5,6);
t171 = -t176 * t116 + t175;
t170 = t177 * t118 - t174;
t169 = rSges(5,1) + pkin(3);
t114 = qJ(1) + qJ(2);
t109 = sin(t114);
t110 = cos(t114);
t168 = t171 * t109 - t172 * t110;
t167 = t172 * t109 + t171 * t110;
t166 = t170 * t109 - t173 * t110;
t165 = t173 * t109 + t170 * t110;
t164 = t176 * t118 + t174;
t163 = t177 * t116 + t175;
t162 = Icges(4,3) + Icges(5,3);
t161 = -t172 * t116 + t173 * t118;
t160 = rSges(5,3) + qJ(4);
t159 = -rSges(5,2) * t116 + t169 * t118;
t108 = V_base(6) + qJD(1);
t106 = qJD(2) + t108;
t85 = -qJD(3) * t110 + V_base(5);
t86 = qJD(3) * t109 + V_base(4);
t158 = (-t167 * t116 + t165 * t118) * t86 + (-t168 * t116 + t166 * t118) * t85 + (-t164 * t116 + t163 * t118) * t106;
t157 = (t162 * t109 + t161 * t110) * t86 + (t161 * t109 - t162 * t110) * t85 + (t173 * t116 + t172 * t118) * t106;
t154 = -pkin(4) - pkin(5);
t117 = sin(qJ(1));
t152 = pkin(1) * t117;
t119 = cos(qJ(1));
t151 = pkin(1) * t119;
t148 = t159 * t109 - t160 * t110;
t147 = t160 * t109 + t159 * t110;
t146 = Icges(2,4) * t117;
t145 = Icges(3,4) * t109;
t140 = t108 * t151 + V_base(2);
t139 = V_base(4) * t152 + V_base(3);
t138 = V_base(5) * pkin(4) + V_base(1);
t135 = rSges(5,2) * t118 + t169 * t116;
t134 = rSges(4,1) * t118 - rSges(4,2) * t116;
t126 = V_base(5) * pkin(5) - t108 * t152 + t138;
t81 = pkin(2) * t110 + pkin(6) * t109;
t123 = t106 * t81 + t154 * V_base(4) + t140;
t80 = pkin(2) * t109 - pkin(6) * t110;
t122 = V_base(4) * t80 + (-t81 - t151) * V_base(5) + t139;
t111 = Icges(2,4) * t119;
t105 = Icges(3,4) * t110;
t102 = rSges(2,1) * t119 - t117 * rSges(2,2);
t101 = t117 * rSges(2,1) + rSges(2,2) * t119;
t100 = rSges(4,1) * t116 + rSges(4,2) * t118;
t98 = Icges(2,1) * t119 - t146;
t97 = Icges(2,1) * t117 + t111;
t94 = -Icges(2,2) * t117 + t111;
t93 = Icges(2,2) * t119 + t146;
t84 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t83 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t82 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t79 = rSges(3,1) * t110 - rSges(3,2) * t109;
t78 = rSges(3,1) * t109 + rSges(3,2) * t110;
t77 = Icges(3,1) * t110 - t145;
t76 = Icges(3,1) * t109 + t105;
t75 = -Icges(3,2) * t109 + t105;
t74 = Icges(3,2) * t110 + t145;
t69 = rSges(4,3) * t109 + t110 * t134;
t67 = -rSges(4,3) * t110 + t109 * t134;
t53 = V_base(5) * rSges(2,3) - t101 * t108 + t138;
t52 = t102 * t108 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t51 = t101 * V_base(4) - t102 * V_base(5) + V_base(3);
t48 = V_base(5) * rSges(3,3) - t106 * t78 + t126;
t47 = t106 * t79 + (-rSges(3,3) + t154) * V_base(4) + t140;
t46 = V_base(4) * t78 + (-t79 - t151) * V_base(5) + t139;
t45 = t100 * t85 + (-t67 - t80) * t106 + t126;
t44 = -t100 * t86 + t106 * t69 + t123;
t43 = t86 * t67 - t85 * t69 + t122;
t42 = qJD(4) * t109 + t135 * t85 + (-t80 - t148) * t106 + t126;
t41 = -qJD(4) * t110 + t106 * t147 - t135 * t86 + t123;
t40 = -t147 * t85 + t148 * t86 + t122;
t1 = m(1) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(2) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(4) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + (t158 * t109 - t157 * t110) * t85 / 0.2e1 + (t157 * t109 + t158 * t110) * t86 / 0.2e1 + ((t116 * t165 + t118 * t167) * t86 + (t166 * t116 + t118 * t168) * t85 + (t163 * t116 + t164 * t118 + Icges(3,3)) * t106) * t106 / 0.2e1 + ((-t109 * t74 + t110 * t76 - t117 * t93 + t119 * t97 + Icges(1,4)) * V_base(5) + (-t109 * t75 + t110 * t77 - t117 * t94 + t119 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t109 * t76 + t110 * t74 + t117 * t97 + t119 * t93 + Icges(1,2)) * V_base(5) + (t109 * t77 + t110 * t75 + t117 * t98 + t119 * t94 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t106 * (Icges(3,5) * t110 - Icges(3,6) * t109) + V_base(5) * t106 * (Icges(3,5) * t109 + Icges(3,6) * t110) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t117 + Icges(2,6) * t119) * V_base(5) + (Icges(2,5) * t119 - Icges(2,6) * t117) * V_base(4) + Icges(2,3) * t108 / 0.2e1) * t108;
T = t1;
