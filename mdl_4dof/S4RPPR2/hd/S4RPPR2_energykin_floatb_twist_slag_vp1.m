% Calculate kinetic energy for
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:25
% EndTime: 2019-03-08 18:28:26
% DurationCPUTime: 0.83s
% Computational Cost: add. (441->154), mult. (546->184), div. (0->0), fcn. (458->6), ass. (0->80)
t151 = Icges(2,4) - Icges(3,5);
t150 = Icges(2,1) + Icges(3,1);
t149 = Icges(2,2) + Icges(3,3);
t112 = sin(qJ(1));
t148 = t151 * t112;
t113 = cos(qJ(1));
t147 = t151 * t113;
t146 = -t149 * t113 - t148;
t145 = t149 * t112 - t147;
t144 = t150 * t112 + t147;
t143 = t150 * t113 - t148;
t142 = V_base(4) / 0.2e1;
t141 = V_base(5) / 0.2e1;
t140 = Icges(3,4) + Icges(2,5);
t139 = Icges(2,6) - Icges(3,6);
t138 = rSges(5,3) + pkin(5);
t137 = pkin(2) * t113;
t136 = t112 * pkin(2);
t132 = cos(pkin(6));
t135 = t132 * pkin(3);
t131 = sin(pkin(6));
t121 = t112 * t131;
t74 = -t113 * t132 - t121;
t134 = Icges(4,4) * t74;
t128 = pkin(6) + qJ(4);
t118 = sin(t128);
t119 = cos(t128);
t70 = -t112 * t118 - t113 * t119;
t133 = Icges(5,4) * t70;
t92 = t112 * pkin(1) - qJ(2) * t113;
t127 = V_base(4) * t92 + V_base(3);
t126 = V_base(5) * pkin(4) + V_base(1);
t123 = -t92 - t136;
t95 = pkin(1) * t113 + t112 * qJ(2);
t122 = -t95 - t137;
t106 = V_base(6) + qJD(1);
t120 = t113 * t131;
t117 = qJD(2) * t112 + t126;
t116 = V_base(4) * t136 - qJD(3) + t127;
t115 = -qJD(2) * t113 + t106 * t95 + V_base(2);
t114 = V_base(4) * qJ(3) + t106 * t137 + t115;
t104 = -qJD(4) + t106;
t97 = rSges(2,1) * t113 - t112 * rSges(2,2);
t96 = rSges(3,1) * t113 + t112 * rSges(3,3);
t94 = t112 * rSges(2,1) + rSges(2,2) * t113;
t93 = t112 * rSges(3,1) - rSges(3,3) * t113;
t79 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t78 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t77 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t75 = t112 * t132 - t120;
t72 = Icges(4,4) * t75;
t71 = t112 * t119 - t113 * t118;
t69 = Icges(5,4) * t71;
t68 = pkin(3) * t121 + t135 * t113;
t67 = -pkin(3) * t120 + t135 * t112;
t66 = V_base(5) * rSges(2,3) - t106 * t94 + t126;
t65 = t106 * t97 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t64 = t94 * V_base(4) - t97 * V_base(5) + V_base(3);
t63 = -rSges(4,1) * t74 + rSges(4,2) * t75;
t62 = rSges(4,1) * t75 + rSges(4,2) * t74;
t61 = -Icges(4,1) * t74 + t72;
t60 = Icges(4,1) * t75 + t134;
t59 = Icges(4,2) * t75 - t134;
t58 = Icges(4,2) * t74 + t72;
t55 = -rSges(5,1) * t70 + rSges(5,2) * t71;
t54 = rSges(5,1) * t71 + rSges(5,2) * t70;
t53 = -Icges(5,1) * t70 + t69;
t52 = Icges(5,1) * t71 + t133;
t51 = Icges(5,2) * t71 - t133;
t50 = Icges(5,2) * t70 + t69;
t47 = V_base(5) * rSges(3,2) + (-t92 - t93) * t106 + t117;
t46 = t106 * t96 + (-rSges(3,2) - pkin(4)) * V_base(4) + t115;
t45 = V_base(4) * t93 + (-t95 - t96) * V_base(5) + t127;
t44 = (-qJ(3) - rSges(4,3)) * V_base(5) + (t123 - t62) * t106 + t117;
t43 = t106 * t63 + (rSges(4,3) - pkin(4)) * V_base(4) + t114;
t42 = V_base(4) * t62 + (t122 - t63) * V_base(5) + t116;
t41 = -t104 * t54 + (-qJ(3) - t138) * V_base(5) + (t123 - t67) * t106 + t117;
t40 = t104 * t55 + t106 * t68 + (-pkin(4) + t138) * V_base(4) + t114;
t39 = (t54 + t67) * V_base(4) + (t122 - t55 - t68) * V_base(5) + t116;
t1 = m(1) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(2) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(3) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + m(4) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + m(5) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + Icges(5,3) * t104 ^ 2 / 0.2e1 + (Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(4,3) / 0.2e1) * t106 ^ 2 + (Icges(1,6) * V_base(6) + (-Icges(5,5) * t71 - Icges(5,6) * t70) * t104 + (-Icges(4,5) * t75 - Icges(4,6) * t74 + t140 * t112 + t139 * t113) * t106 + Icges(1,2) * t141) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(5,5) * t70 - Icges(5,6) * t71) * t104 + (Icges(4,5) * t74 - Icges(4,6) * t75 - t139 * t112 + t140 * t113) * t106 + Icges(1,1) * t142) * V_base(4) + ((t146 * t112 + t144 * t113 + t50 * t71 - t52 * t70 + t58 * t75 - t60 * t74) * V_base(5) + (t145 * t112 + t143 * t113 + t71 * t51 - t70 * t53 + t75 * t59 - t74 * t61) * V_base(4)) * t142 + ((t144 * t112 - t146 * t113 + t70 * t50 + t71 * t52 + t74 * t58 + t75 * t60) * V_base(5) + (t143 * t112 - t145 * t113 + t51 * t70 + t53 * t71 + t59 * t74 + t61 * t75) * V_base(4)) * t141;
T  = t1;
