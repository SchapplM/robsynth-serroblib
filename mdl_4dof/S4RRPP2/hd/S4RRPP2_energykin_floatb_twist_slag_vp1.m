% Calculate kinetic energy for
% S4RRPP2
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:58
% EndTime: 2019-03-08 18:33:59
% DurationCPUTime: 0.65s
% Computational Cost: add. (426->127), mult. (366->155), div. (0->0), fcn. (196->4), ass. (0->63)
t130 = Icges(3,4) - Icges(5,4) - Icges(4,5);
t129 = Icges(3,1) + Icges(4,1) + Icges(5,1);
t128 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t96 = qJ(1) + qJ(2);
t91 = sin(t96);
t127 = t130 * t91;
t92 = cos(t96);
t126 = t130 * t92;
t125 = rSges(5,1) + pkin(3);
t124 = -t128 * t92 - t127;
t123 = t128 * t91 - t126;
t122 = t129 * t91 + t126;
t121 = t129 * t92 - t127;
t120 = Icges(4,4) + Icges(3,5) - Icges(5,5);
t119 = Icges(3,6) - Icges(4,6) + Icges(5,6);
t118 = -pkin(4) - pkin(5);
t97 = sin(qJ(1));
t117 = pkin(1) * t97;
t98 = cos(qJ(1));
t116 = pkin(1) * t98;
t115 = Icges(2,4) * t97;
t111 = -qJ(4) - rSges(5,3);
t90 = V_base(6) + qJD(1);
t110 = t90 * t116 + V_base(2);
t109 = V_base(4) * t117 + V_base(3);
t108 = V_base(5) * pkin(4) + V_base(1);
t68 = t92 * pkin(2) + t91 * qJ(3);
t105 = -t68 - t116;
t104 = -t92 * rSges(5,2) + t125 * t91;
t103 = t91 * rSges(5,2) + t125 * t92;
t64 = t91 * pkin(2) - t92 * qJ(3);
t102 = V_base(4) * t64 + t109;
t89 = qJD(2) + t90;
t101 = -qJD(3) * t92 + t89 * t68 + t110;
t100 = V_base(5) * pkin(5) - t90 * t117 + t108;
t99 = qJD(3) * t91 + t100;
t93 = Icges(2,4) * t98;
t82 = t98 * rSges(2,1) - t97 * rSges(2,2);
t81 = t97 * rSges(2,1) + t98 * rSges(2,2);
t80 = Icges(2,1) * t98 - t115;
t79 = Icges(2,1) * t97 + t93;
t78 = -Icges(2,2) * t97 + t93;
t77 = Icges(2,2) * t98 + t115;
t74 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t73 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t72 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t71 = t92 * rSges(3,1) - t91 * rSges(3,2);
t70 = t92 * rSges(4,1) + t91 * rSges(4,3);
t67 = t91 * rSges(3,1) + t92 * rSges(3,2);
t66 = t91 * rSges(4,1) - t92 * rSges(4,3);
t43 = V_base(5) * rSges(2,3) - t90 * t81 + t108;
t42 = t90 * t82 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t41 = V_base(4) * t81 - V_base(5) * t82 + V_base(3);
t40 = V_base(5) * rSges(3,3) - t89 * t67 + t100;
t39 = t89 * t71 + (-rSges(3,3) + t118) * V_base(4) + t110;
t38 = V_base(4) * t67 + (-t71 - t116) * V_base(5) + t109;
t37 = V_base(5) * rSges(4,2) + (-t64 - t66) * t89 + t99;
t36 = t89 * t70 + (-rSges(4,2) + t118) * V_base(4) + t101;
t35 = V_base(4) * t66 + (t105 - t70) * V_base(5) + t102;
t34 = t111 * V_base(5) + (-t104 - t64) * t89 + t99;
t33 = t103 * t89 + (-t111 + t118) * V_base(4) + t101;
t32 = -qJD(4) + t104 * V_base(4) + (-t103 + t105) * V_base(5) + t102;
t1 = m(1) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(3) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t35 ^ 2 + t36 ^ 2 + t37 ^ 2) / 0.2e1 + m(5) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t97 + Icges(2,6) * t98) * V_base(5) + (Icges(2,5) * t98 - Icges(2,6) * t97) * V_base(4) + Icges(2,3) * t90 / 0.2e1) * t90 + ((Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1 + Icges(5,3) / 0.2e1) * t89 + (t119 * V_base(5) + t120 * V_base(4)) * t92 + (-t119 * V_base(4) + t120 * V_base(5)) * t91) * t89 + ((t122 * t92 + t124 * t91 - t97 * t77 + t98 * t79) * V_base(5) + (t121 * t92 + t123 * t91 - t97 * t78 + t98 * t80) * V_base(4)) * V_base(4) / 0.2e1 + ((t122 * t91 - t124 * t92 + t98 * t77 + t97 * t79) * V_base(5) + (t121 * t91 - t123 * t92 + t98 * t78 + t97 * t80) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
