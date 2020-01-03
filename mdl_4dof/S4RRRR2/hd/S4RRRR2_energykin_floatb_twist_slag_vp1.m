% Calculate kinetic energy for
% S4RRRR2
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:06
% EndTime: 2019-12-31 17:23:07
% DurationCPUTime: 0.88s
% Computational Cost: add. (773->190), mult. (652->269), div. (0->0), fcn. (482->8), ass. (0->97)
t162 = -pkin(4) - pkin(5);
t126 = sin(qJ(1));
t160 = pkin(1) * t126;
t128 = cos(qJ(1));
t159 = pkin(1) * t128;
t125 = sin(qJ(3));
t158 = pkin(3) * t125;
t127 = cos(qJ(3));
t157 = pkin(3) * t127;
t155 = Icges(2,4) * t126;
t124 = qJ(1) + qJ(2);
t117 = sin(t124);
t154 = Icges(3,4) * t117;
t153 = Icges(4,4) * t125;
t152 = Icges(4,4) * t127;
t123 = qJ(3) + qJ(4);
t116 = sin(t123);
t151 = Icges(5,4) * t116;
t118 = cos(t123);
t150 = Icges(5,4) * t118;
t115 = V_base(6) + qJD(1);
t149 = t115 * t159 + V_base(2);
t148 = V_base(4) * t160 + V_base(3);
t147 = V_base(5) * pkin(4) + V_base(1);
t96 = qJD(3) * t117 + V_base(4);
t144 = rSges(4,1) * t127 - rSges(4,2) * t125;
t143 = rSges(5,1) * t118 - rSges(5,2) * t116;
t142 = Icges(4,1) * t127 - t153;
t141 = Icges(5,1) * t118 - t151;
t140 = -Icges(4,2) * t125 + t152;
t139 = -Icges(5,2) * t116 + t150;
t138 = Icges(4,5) * t127 - Icges(4,6) * t125;
t137 = Icges(5,5) * t118 - Icges(5,6) * t116;
t136 = V_base(5) * pkin(5) - t115 * t160 + t147;
t113 = qJD(2) + t115;
t119 = cos(t124);
t75 = V_base(5) + (-qJD(3) - qJD(4)) * t119;
t76 = qJD(4) * t117 + t96;
t135 = (Icges(5,5) * t116 + Icges(5,6) * t118) * t113 + (-Icges(5,3) * t119 + t117 * t137) * t75 + (Icges(5,3) * t117 + t119 * t137) * t76;
t95 = -qJD(3) * t119 + V_base(5);
t134 = (Icges(4,5) * t125 + Icges(4,6) * t127) * t113 + (-Icges(4,3) * t119 + t117 * t138) * t95 + (Icges(4,3) * t117 + t119 * t138) * t96;
t91 = t119 * pkin(2) + t117 * pkin(6);
t133 = t113 * t91 + t162 * V_base(4) + t149;
t90 = t117 * pkin(2) - t119 * pkin(6);
t132 = V_base(4) * t90 + (-t91 - t159) * V_base(5) + t148;
t57 = -Icges(5,6) * t119 + t117 * t139;
t58 = Icges(5,6) * t117 + t119 * t139;
t59 = -Icges(5,5) * t119 + t117 * t141;
t60 = Icges(5,5) * t117 + t119 * t141;
t81 = Icges(5,2) * t118 + t151;
t84 = Icges(5,1) * t116 + t150;
t131 = (-t116 * t58 + t118 * t60) * t76 + (-t116 * t57 + t118 * t59) * t75 + (-t116 * t81 + t118 * t84) * t113;
t100 = Icges(4,2) * t127 + t153;
t103 = Icges(4,1) * t125 + t152;
t68 = -Icges(4,6) * t119 + t117 * t140;
t69 = Icges(4,6) * t117 + t119 * t140;
t70 = -Icges(4,5) * t119 + t117 * t142;
t71 = Icges(4,5) * t117 + t119 * t142;
t130 = (-t125 * t69 + t127 * t71) * t96 + (-t125 * t68 + t127 * t70) * t95 + (-t100 * t125 + t103 * t127) * t113;
t120 = Icges(2,4) * t128;
t112 = Icges(3,4) * t119;
t108 = rSges(2,1) * t128 - rSges(2,2) * t126;
t107 = rSges(2,1) * t126 + rSges(2,2) * t128;
t106 = rSges(4,1) * t125 + rSges(4,2) * t127;
t105 = Icges(2,1) * t128 - t155;
t104 = Icges(2,1) * t126 + t120;
t102 = -Icges(2,2) * t126 + t120;
t101 = Icges(2,2) * t128 + t155;
t94 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t93 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t92 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t89 = rSges(3,1) * t119 - rSges(3,2) * t117;
t88 = rSges(3,1) * t117 + rSges(3,2) * t119;
t87 = rSges(5,1) * t116 + rSges(5,2) * t118;
t86 = Icges(3,1) * t119 - t154;
t85 = Icges(3,1) * t117 + t112;
t83 = -Icges(3,2) * t117 + t112;
t82 = Icges(3,2) * t119 + t154;
t73 = rSges(4,3) * t117 + t119 * t144;
t72 = -rSges(4,3) * t119 + t117 * t144;
t65 = V_base(5) * rSges(2,3) - t107 * t115 + t147;
t64 = t108 * t115 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t63 = t107 * V_base(4) - t108 * V_base(5) + V_base(3);
t62 = rSges(5,3) * t117 + t119 * t143;
t61 = -rSges(5,3) * t119 + t117 * t143;
t54 = pkin(7) * t117 + t119 * t157;
t53 = -pkin(7) * t119 + t117 * t157;
t52 = V_base(5) * rSges(3,3) - t113 * t88 + t136;
t51 = t113 * t89 + (-rSges(3,3) + t162) * V_base(4) + t149;
t50 = t88 * V_base(4) + (-t89 - t159) * V_base(5) + t148;
t49 = t106 * t95 + (-t72 - t90) * t113 + t136;
t48 = -t106 * t96 + t113 * t73 + t133;
t47 = t72 * t96 - t73 * t95 + t132;
t46 = t95 * t158 + t75 * t87 + (-t53 - t61 - t90) * t113 + t136;
t45 = -t96 * t158 - t76 * t87 + (t54 + t62) * t113 + t133;
t44 = t53 * t96 - t54 * t95 + t61 * t76 - t62 * t75 + t132;
t1 = m(1) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(2) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t96 * (t134 * t117 + t130 * t119) / 0.2e1 + t95 * (t130 * t117 - t134 * t119) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t76 * (t135 * t117 + t131 * t119) / 0.2e1 + t75 * (t131 * t117 - t135 * t119) / 0.2e1 + ((t125 * t71 + t127 * t69) * t96 + (t125 * t70 + t127 * t68) * t95 + (t116 * t60 + t118 * t58) * t76 + (t116 * t59 + t118 * t57) * t75 + (t127 * t100 + t125 * t103 + t116 * t84 + t118 * t81 + Icges(3,3)) * t113) * t113 / 0.2e1 + ((-t101 * t126 + t104 * t128 - t117 * t82 + t119 * t85 + Icges(1,4)) * V_base(5) + (-t126 * t102 + t128 * t105 - t117 * t83 + t119 * t86 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t128 * t101 + t126 * t104 + t117 * t85 + t119 * t82 + Icges(1,2)) * V_base(5) + (t102 * t128 + t105 * t126 + t117 * t86 + t119 * t83 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t113 * (Icges(3,5) * t119 - Icges(3,6) * t117) + V_base(5) * t113 * (Icges(3,5) * t117 + Icges(3,6) * t119) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t126 + Icges(2,6) * t128) * V_base(5) + (Icges(2,5) * t128 - Icges(2,6) * t126) * V_base(4) + Icges(2,3) * t115 / 0.2e1) * t115;
T = t1;
