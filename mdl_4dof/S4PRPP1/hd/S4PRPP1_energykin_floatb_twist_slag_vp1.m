% Calculate kinetic energy for
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:53
% EndTime: 2019-03-08 18:17:54
% DurationCPUTime: 0.64s
% Computational Cost: add. (397->131), mult. (368->157), div. (0->0), fcn. (198->4), ass. (0->62)
t128 = Icges(3,4) + Icges(4,6) - Icges(5,6);
t127 = Icges(3,1) + Icges(4,2) + Icges(5,3);
t126 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t94 = pkin(5) + qJ(2);
t89 = cos(t94);
t125 = t128 * t89;
t88 = sin(t94);
t124 = t128 * t88;
t123 = t127 * t89 - t124;
t122 = t127 * t88 + t125;
t121 = -t126 * t89 - t124;
t120 = t126 * t88 - t125;
t119 = -Icges(4,4) + Icges(3,5) + Icges(5,5);
t118 = Icges(5,4) + Icges(4,5) - Icges(3,6);
t96 = cos(pkin(5));
t117 = pkin(1) * t96;
t116 = pkin(3) + rSges(5,1);
t115 = -pkin(4) - qJ(1);
t95 = sin(pkin(5));
t114 = Icges(2,4) * t95;
t104 = pkin(1) * V_base(6);
t109 = t96 * t104 + V_base(2);
t108 = V_base(5) * qJ(1) + V_base(1);
t105 = qJD(1) + V_base(3);
t68 = t89 * pkin(2) + t88 * qJ(3);
t103 = -t68 - t117;
t91 = V_base(6) + qJD(2);
t102 = t91 * t68 + t109;
t101 = -t89 * rSges(5,2) + (rSges(5,3) + qJ(4)) * t88;
t100 = V_base(4) * t95 * pkin(1) + t105;
t64 = t88 * pkin(2) - t89 * qJ(3);
t99 = V_base(4) * t64 + t100;
t98 = V_base(5) * pkin(4) - t95 * t104 + t108;
t97 = qJD(3) * t88 + t98;
t90 = Icges(2,4) * t96;
t82 = t96 * rSges(2,1) - t95 * rSges(2,2);
t81 = t95 * rSges(2,1) + t96 * rSges(2,2);
t80 = Icges(2,1) * t96 - t114;
t79 = Icges(2,1) * t95 + t90;
t78 = -Icges(2,2) * t95 + t90;
t77 = Icges(2,2) * t96 + t114;
t74 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t73 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t72 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t71 = t89 * rSges(3,1) - t88 * rSges(3,2);
t70 = -t89 * rSges(4,2) + t88 * rSges(4,3);
t67 = t88 * rSges(3,1) + t89 * rSges(3,2);
t66 = -t88 * rSges(4,2) - t89 * rSges(4,3);
t65 = t88 * rSges(5,2) + t89 * rSges(5,3);
t43 = V_base(5) * rSges(2,3) - V_base(6) * t81 + t108;
t42 = V_base(6) * t82 + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t41 = V_base(4) * t81 - V_base(5) * t82 + t105;
t40 = V_base(5) * rSges(3,3) - t91 * t67 + t98;
t39 = t91 * t71 + (-rSges(3,3) + t115) * V_base(4) + t109;
t38 = V_base(4) * t67 + (-t71 - t117) * V_base(5) + t100;
t37 = V_base(5) * rSges(4,1) + (-t64 - t66) * t91 + t97;
t36 = -qJD(3) * t89 + t91 * t70 + (-rSges(4,1) + t115) * V_base(4) + t102;
t35 = V_base(4) * t66 + (t103 - t70) * V_base(5) + t99;
t34 = qJD(4) * t89 + t116 * V_base(5) + (-t101 - t64) * t91 + t97;
t33 = qJD(4) * t88 + t91 * t65 + (qJ(4) * t91 - qJD(3)) * t89 + (t115 - t116) * V_base(4) + t102;
t32 = t101 * V_base(4) + (-qJ(4) * t89 + t103 - t65) * V_base(5) + t99;
t1 = m(1) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(3) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t35 ^ 2 + t36 ^ 2 + t37 ^ 2) / 0.2e1 + m(5) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + ((Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6) + (Icges(2,5) * t95 + Icges(2,6) * t96 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t96 - Icges(2,6) * t95 + Icges(1,5)) * V_base(4)) * V_base(6) + ((Icges(3,3) / 0.2e1 + Icges(4,1) / 0.2e1 + Icges(5,1) / 0.2e1) * t91 + (-t118 * V_base(5) + t119 * V_base(4)) * t89 + (t118 * V_base(4) + t119 * V_base(5)) * t88) * t91 + ((t121 * t88 + t122 * t89 - t95 * t77 + t96 * t79) * V_base(5) + (t120 * t88 + t123 * t89 - t95 * t78 + t96 * t80) * V_base(4)) * V_base(4) / 0.2e1 + ((-t121 * t89 + t122 * t88 + t96 * t77 + t95 * t79) * V_base(5) + (-t120 * t89 + t123 * t88 + t96 * t78 + t95 * t80) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
