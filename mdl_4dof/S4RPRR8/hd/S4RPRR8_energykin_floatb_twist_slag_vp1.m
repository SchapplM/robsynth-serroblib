% Calculate kinetic energy for
% S4RPRR8
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:03
% EndTime: 2019-12-31 16:55:04
% DurationCPUTime: 1.10s
% Computational Cost: add. (481->179), mult. (656->247), div. (0->0), fcn. (488->6), ass. (0->93)
t174 = Icges(2,4) + Icges(3,6);
t173 = Icges(2,1) + Icges(3,2);
t172 = -Icges(3,4) + Icges(2,5);
t171 = Icges(3,5) - Icges(2,6);
t170 = Icges(2,2) + Icges(3,3);
t119 = cos(qJ(1));
t169 = t174 * t119;
t117 = sin(qJ(1));
t168 = t174 * t117;
t167 = -t119 * t170 - t168;
t166 = t117 * t170 - t169;
t165 = t117 * t173 + t169;
t164 = t119 * t173 - t168;
t106 = V_base(6) + qJD(1);
t115 = qJ(3) + qJ(4);
t110 = sin(t115);
t111 = cos(t115);
t150 = Icges(5,4) * t110;
t128 = Icges(5,2) * t111 + t150;
t53 = Icges(5,6) * t119 + t117 * t128;
t54 = Icges(5,6) * t117 - t119 * t128;
t149 = Icges(5,4) * t111;
t130 = Icges(5,1) * t110 + t149;
t55 = Icges(5,5) * t119 + t117 * t130;
t56 = Icges(5,5) * t117 - t119 * t130;
t71 = -Icges(5,2) * t110 + t149;
t72 = Icges(5,1) * t111 - t150;
t103 = qJD(3) * t117 + V_base(5);
t74 = qJD(4) * t117 + t103;
t104 = qJD(3) * t119 + V_base(4);
t75 = qJD(4) * t119 + t104;
t161 = (t110 * t55 + t111 * t53) * t75 + (t110 * t56 + t111 * t54) * t74 + (t110 * t72 + t111 * t71) * t106;
t116 = sin(qJ(3));
t118 = cos(qJ(3));
t152 = Icges(4,4) * t116;
t129 = Icges(4,2) * t118 + t152;
t61 = Icges(4,6) * t119 + t117 * t129;
t62 = Icges(4,6) * t117 - t119 * t129;
t151 = Icges(4,4) * t118;
t131 = Icges(4,1) * t116 + t151;
t63 = Icges(4,5) * t119 + t117 * t131;
t64 = Icges(4,5) * t117 - t119 * t131;
t87 = -Icges(4,2) * t116 + t151;
t92 = Icges(4,1) * t118 - t152;
t160 = (t116 * t63 + t118 * t61) * t104 + (t116 * t64 + t118 * t62) * t103 + (t116 * t92 + t118 * t87) * t106;
t158 = pkin(3) * t116;
t157 = pkin(3) * t118;
t156 = t117 * pkin(5);
t155 = t119 * pkin(5);
t95 = pkin(1) * t117 - qJ(2) * t119;
t146 = V_base(4) * t95 + V_base(3);
t145 = V_base(5) * pkin(4) + V_base(1);
t142 = -t95 - t156;
t141 = qJD(2) * t117 + t145;
t140 = V_base(5) * pkin(2) + t141;
t139 = rSges(4,1) * t116 + rSges(4,2) * t118;
t138 = rSges(5,1) * t110 + rSges(5,2) * t111;
t127 = Icges(4,5) * t116 + Icges(4,6) * t118;
t126 = Icges(5,5) * t110 + Icges(5,6) * t111;
t99 = pkin(1) * t119 + qJ(2) * t117;
t125 = -qJD(2) * t119 + t106 * t99 + V_base(2);
t124 = (Icges(5,5) * t111 - Icges(5,6) * t110) * t106 + (Icges(5,3) * t119 + t117 * t126) * t75 + (Icges(5,3) * t117 - t119 * t126) * t74;
t123 = (Icges(4,3) * t117 - t119 * t127) * t103 + (Icges(4,3) * t119 + t117 * t127) * t104 + (Icges(4,5) * t118 - Icges(4,6) * t116) * t106;
t122 = V_base(4) * t156 + (-t99 - t155) * V_base(5) + t146;
t121 = t106 * t155 + (-pkin(2) - pkin(4)) * V_base(4) + t125;
t101 = rSges(2,1) * t119 - rSges(2,2) * t117;
t100 = -rSges(3,2) * t119 + rSges(3,3) * t117;
t98 = rSges(4,1) * t118 - rSges(4,2) * t116;
t97 = rSges(2,1) * t117 + rSges(2,2) * t119;
t96 = -rSges(3,2) * t117 - rSges(3,3) * t119;
t79 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t78 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t77 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t73 = rSges(5,1) * t111 - rSges(5,2) * t110;
t68 = pkin(6) * t119 + t117 * t158;
t67 = pkin(6) * t117 - t119 * t158;
t66 = rSges(4,3) * t117 - t119 * t139;
t65 = rSges(4,3) * t119 + t117 * t139;
t58 = rSges(5,3) * t117 - t119 * t138;
t57 = rSges(5,3) * t119 + t117 * t138;
t50 = V_base(5) * rSges(2,3) - t106 * t97 + t145;
t49 = t101 * t106 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t48 = -t101 * V_base(5) + t97 * V_base(4) + V_base(3);
t47 = V_base(5) * rSges(3,1) + (-t95 - t96) * t106 + t141;
t46 = t100 * t106 + (-rSges(3,1) - pkin(4)) * V_base(4) + t125;
t45 = t96 * V_base(4) + (-t100 - t99) * V_base(5) + t146;
t44 = t103 * t98 + (t142 - t66) * t106 + t140;
t43 = -t104 * t98 + t106 * t65 + t121;
t42 = -t103 * t65 + t104 * t66 + t122;
t41 = t103 * t157 + t73 * t74 + (t142 - t58 - t67) * t106 + t140;
t40 = -t104 * t157 - t73 * t75 + (t57 + t68) * t106 + t121;
t39 = -t103 * t68 + t104 * t67 - t57 * t74 + t58 * t75 + t122;
t1 = m(1) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(2) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + m(3) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + m(4) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + t104 * (t160 * t117 + t123 * t119) / 0.2e1 + t103 * (t123 * t117 - t160 * t119) / 0.2e1 + m(5) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) / 0.2e1 + t75 * (t161 * t117 + t124 * t119) / 0.2e1 + t74 * (t124 * t117 - t161 * t119) / 0.2e1 + ((t117 * t167 + t165 * t119 + Icges(1,4)) * V_base(5) + (t166 * t117 + t164 * t119 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t165 * t117 - t167 * t119 + Icges(1,2)) * V_base(5) + (t117 * t164 - t119 * t166 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t116 * t61 + t118 * t63) * t104 + (-t116 * t62 + t118 * t64) * t103 + (-t110 * t53 + t111 * t55) * t75 + (-t110 * t54 + t111 * t56) * t74 + (-t110 * t71 + t111 * t72 - t116 * t87 + t118 * t92 + Icges(3,1) + Icges(2,3)) * t106) * t106 / 0.2e1 + t106 * V_base(5) * (t117 * t172 - t119 * t171) + t106 * V_base(4) * (t117 * t171 + t119 * t172) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
