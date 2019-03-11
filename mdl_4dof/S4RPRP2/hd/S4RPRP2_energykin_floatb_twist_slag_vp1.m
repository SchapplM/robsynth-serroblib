% Calculate kinetic energy for
% S4RPRP2
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
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:41
% EndTime: 2019-03-08 18:30:41
% DurationCPUTime: 0.77s
% Computational Cost: add. (368->137), mult. (546->164), div. (0->0), fcn. (458->4), ass. (0->70)
t149 = Icges(2,4) - Icges(3,5);
t148 = Icges(4,4) + Icges(5,4);
t147 = Icges(2,1) + Icges(3,1);
t146 = Icges(4,1) + Icges(5,1);
t145 = Icges(2,2) + Icges(3,3);
t144 = Icges(4,2) + Icges(5,2);
t102 = sin(qJ(1));
t103 = cos(qJ(1));
t125 = sin(qJ(3));
t110 = t103 * t125;
t126 = cos(qJ(3));
t68 = t102 * t126 - t110;
t143 = t148 * t68;
t142 = t149 * t102;
t141 = t149 * t103;
t111 = t102 * t125;
t67 = -t103 * t126 - t111;
t140 = t148 * t67;
t139 = t144 * t67 + t143;
t138 = t144 * t68 - t140;
t137 = t146 * t68 + t140;
t136 = -t146 * t67 + t143;
t135 = -t145 * t103 - t142;
t134 = t145 * t102 - t141;
t133 = t147 * t102 + t141;
t132 = t147 * t103 - t142;
t131 = Icges(3,4) + Icges(2,5);
t130 = Icges(4,5) + Icges(5,5);
t129 = Icges(2,6) - Icges(3,6);
t128 = -Icges(4,6) - Icges(5,6);
t127 = t126 * pkin(3);
t124 = t102 * pkin(2);
t123 = t103 * pkin(2);
t122 = t68 * rSges(5,1) + t67 * rSges(5,2) - pkin(3) * t110 + t127 * t102;
t121 = -t67 * rSges(5,1) + t68 * rSges(5,2) + pkin(3) * t111 + t127 * t103;
t118 = rSges(5,3) + qJ(4);
t84 = t102 * pkin(1) - t103 * qJ(2);
t115 = V_base(4) * t84 + V_base(3);
t114 = V_base(5) * pkin(4) + V_base(1);
t87 = t103 * pkin(1) + t102 * qJ(2);
t109 = -t87 - t123;
t96 = V_base(6) + qJD(1);
t108 = V_base(4) * t124 + t115;
t107 = qJD(2) * t102 + t114;
t106 = -qJD(2) * t103 + t96 * t87 + V_base(2);
t105 = V_base(4) * pkin(5) + t96 * t123 + t106;
t104 = (-t84 - t124) * t96 + t107;
t92 = -qJD(3) + t96;
t89 = t103 * rSges(2,1) - t102 * rSges(2,2);
t88 = t103 * rSges(3,1) + t102 * rSges(3,3);
t86 = t102 * rSges(2,1) + t103 * rSges(2,2);
t85 = t102 * rSges(3,1) - t103 * rSges(3,3);
t71 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t70 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t69 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t60 = V_base(5) * rSges(2,3) - t96 * t86 + t114;
t59 = t96 * t89 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t58 = V_base(4) * t86 - V_base(5) * t89 + V_base(3);
t57 = -t67 * rSges(4,1) + t68 * rSges(4,2);
t55 = t68 * rSges(4,1) + t67 * rSges(4,2);
t41 = V_base(5) * rSges(3,2) + (-t84 - t85) * t96 + t107;
t40 = t96 * t88 + (-rSges(3,2) - pkin(4)) * V_base(4) + t106;
t39 = V_base(4) * t85 + (-t87 - t88) * V_base(5) + t115;
t38 = -t92 * t55 + (-pkin(5) - rSges(4,3)) * V_base(5) + t104;
t37 = t92 * t57 + (rSges(4,3) - pkin(4)) * V_base(4) + t105;
t36 = V_base(4) * t55 + (t109 - t57) * V_base(5) + t108;
t35 = -t122 * t92 + (-pkin(5) - t118) * V_base(5) + t104;
t34 = t121 * t92 + (-pkin(4) + t118) * V_base(4) + t105;
t33 = -qJD(4) + t122 * V_base(4) + (t109 - t121) * V_base(5) + t108;
t1 = m(1) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(3) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) / 0.2e1 + m(4) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + m(5) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * t96 + (t129 * V_base(5) + t131 * V_base(4)) * t103 + (-t129 * V_base(4) + t131 * V_base(5)) * t102) * t96 + ((Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t92 + (t128 * V_base(4) - t130 * V_base(5)) * t68 + (t128 * V_base(5) + t130 * V_base(4)) * t67) * t92 + ((t135 * t102 + t133 * t103 - t137 * t67 + t139 * t68) * V_base(5) + (t102 * t134 + t103 * t132 - t136 * t67 + t138 * t68) * V_base(4)) * V_base(4) / 0.2e1 + ((t133 * t102 - t135 * t103 + t137 * t68 + t139 * t67) * V_base(5) + (t102 * t132 - t103 * t134 + t136 * t68 + t138 * t67) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
