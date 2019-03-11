% Calculate kinetic energy for
% S4PRRP2
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
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:58
% EndTime: 2019-03-08 18:23:59
% DurationCPUTime: 0.38s
% Computational Cost: add. (314->135), mult. (300->159), div. (0->0), fcn. (132->4), ass. (0->64)
t105 = rSges(5,1) + pkin(3);
t104 = Icges(4,5) + Icges(5,5);
t103 = Icges(4,6) + Icges(5,6);
t80 = cos(qJ(2));
t74 = Icges(3,4) * t80;
t79 = sin(qJ(2));
t102 = Icges(3,1) * t79 / 0.2e1 + t74 / 0.2e1;
t78 = qJ(2) + qJ(3);
t72 = cos(t78);
t99 = t72 / 0.2e1;
t98 = t80 / 0.2e1;
t97 = -pkin(4) - pkin(5);
t96 = pkin(2) * t79;
t95 = pkin(2) * t80;
t94 = Icges(3,4) * t79;
t71 = sin(t78);
t93 = Icges(4,4) * t71;
t92 = Icges(5,4) * t71;
t91 = qJ(4) + rSges(5,3);
t90 = V_base(4) * qJ(1) + V_base(3);
t89 = qJD(1) + V_base(2);
t88 = -pkin(1) - t95;
t87 = rSges(5,2) * t72 + t105 * t71;
t86 = -rSges(5,2) * t71 + t105 * t72;
t70 = V_base(6) + qJD(2);
t85 = V_base(4) * t96 + t90;
t84 = V_base(6) * pkin(1) + t89;
t83 = t70 * t95 + t84;
t82 = V_base(5) * pkin(4) - V_base(6) * qJ(1) + V_base(1);
t81 = V_base(5) * pkin(5) - t70 * t96 + t82;
t69 = qJD(3) + t70;
t68 = Icges(4,4) * t72;
t67 = Icges(5,4) * t72;
t64 = rSges(3,1) * t80 - t79 * rSges(3,2);
t63 = t79 * rSges(3,1) + rSges(3,2) * t80;
t62 = Icges(3,1) * t80 - t94;
t60 = -Icges(3,2) * t79 + t74;
t59 = Icges(3,2) * t80 + t94;
t56 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t55 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t54 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t53 = V_base(6) * rSges(2,1) + V_base(4) * rSges(2,2) + t89;
t52 = rSges(4,1) * t72 - rSges(4,2) * t71;
t50 = rSges(4,1) * t71 + rSges(4,2) * t72;
t48 = Icges(4,1) * t72 - t93;
t47 = Icges(4,1) * t71 + t68;
t46 = Icges(5,1) * t72 - t92;
t45 = Icges(5,1) * t71 + t67;
t44 = -Icges(4,2) * t71 + t68;
t43 = Icges(4,2) * t72 + t93;
t42 = -Icges(5,2) * t71 + t67;
t41 = Icges(5,2) * t72 + t92;
t36 = -V_base(5) * rSges(2,1) + V_base(4) * rSges(2,3) + t90;
t35 = -V_base(5) * rSges(2,2) + V_base(1) + (-rSges(2,3) - qJ(1)) * V_base(6);
t34 = V_base(5) * rSges(3,3) - t63 * t70 + t82;
t33 = t64 * t70 + (-rSges(3,3) - pkin(4)) * V_base(4) + t84;
t32 = t63 * V_base(4) + (-pkin(1) - t64) * V_base(5) + t90;
t31 = V_base(5) * rSges(4,3) - t50 * t69 + t81;
t30 = t52 * t69 + (-rSges(4,3) + t97) * V_base(4) + t83;
t29 = V_base(4) * t50 + (-t52 + t88) * V_base(5) + t85;
t28 = -t87 * t69 + t91 * V_base(5) + t81;
t27 = t86 * t69 + (-t91 + t97) * V_base(4) + t83;
t26 = qJD(4) + t87 * V_base(4) + (-t86 + t88) * V_base(5) + t85;
t1 = m(1) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + m(2) * (t35 ^ 2 + t36 ^ 2 + t53 ^ 2) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + Icges(3,3) * t70 ^ 2 / 0.2e1 + m(4) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,2) / 0.2e1) * V_base(6) ^ 2 + (Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t69 ^ 2 + ((Icges(3,5) * t79 + Icges(3,6) * t80) * t70 + (Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + t59 * t98 + t79 * t102 + (t43 + t41) * t99 + (t47 + t45) * t71 / 0.2e1) * V_base(5) + (Icges(1,6) - Icges(2,6)) * V_base(6) + (t103 * t72 + t104 * t71) * t69) * V_base(5) + ((Icges(3,5) * t80 - Icges(3,6) * t79) * t70 + (Icges(1,1) / 0.2e1 + Icges(2,1) / 0.2e1 - t79 * t60 / 0.2e1 + t62 * t98 + (t48 + t46) * t99 - (t44 + t42) * t71 / 0.2e1) * V_base(4) + (Icges(1,5) - Icges(2,4)) * V_base(6) + (Icges(1,4) + Icges(2,5) + (t102 + t60 / 0.2e1) * t80 + (-t59 / 0.2e1 + t62 / 0.2e1) * t79 + (t47 / 0.2e1 + t44 / 0.2e1 + t45 / 0.2e1 + t42 / 0.2e1) * t72 + (-t43 / 0.2e1 + t48 / 0.2e1 - t41 / 0.2e1 + t46 / 0.2e1) * t71) * V_base(5) + (-t103 * t71 + t104 * t72) * t69) * V_base(4);
T  = t1;
