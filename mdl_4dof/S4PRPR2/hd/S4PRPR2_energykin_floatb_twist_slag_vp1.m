% Calculate kinetic energy for
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:47
% EndTime: 2019-03-08 18:21:48
% DurationCPUTime: 0.32s
% Computational Cost: add. (347->142), mult. (300->175), div. (0->0), fcn. (132->6), ass. (0->67)
t85 = qJ(2) + pkin(6);
t78 = qJ(4) + t85;
t74 = cos(t78);
t70 = Icges(5,4) * t74;
t73 = sin(t78);
t110 = Icges(5,1) * t73 / 0.2e1 + t70 / 0.2e1;
t77 = cos(t85);
t72 = Icges(4,4) * t77;
t76 = sin(t85);
t109 = Icges(4,1) * t76 / 0.2e1 + t72 / 0.2e1;
t87 = cos(qJ(2));
t82 = Icges(3,4) * t87;
t86 = sin(qJ(2));
t108 = Icges(3,1) * t86 / 0.2e1 + t82 / 0.2e1;
t107 = t74 / 0.2e1;
t106 = t77 / 0.2e1;
t105 = t87 / 0.2e1;
t104 = pkin(2) * t86;
t103 = pkin(2) * t87;
t102 = pkin(3) * t76;
t101 = pkin(3) * t77;
t100 = pkin(5) + rSges(5,3);
t99 = -pkin(4) - qJ(3);
t98 = Icges(3,4) * t86;
t97 = Icges(4,4) * t76;
t96 = Icges(5,4) * t73;
t95 = V_base(4) * qJ(1) + V_base(3);
t94 = qJD(1) + V_base(2);
t93 = -pkin(1) - t103;
t79 = V_base(6) + qJD(2);
t92 = V_base(6) * pkin(1) + t94;
t91 = t79 * t103 + t92;
t90 = V_base(4) * t104 + qJD(3) + t95;
t89 = V_base(5) * pkin(4) - V_base(6) * qJ(1) + V_base(1);
t88 = V_base(5) * qJ(3) + t89;
t75 = qJD(4) + t79;
t68 = rSges(3,1) * t87 - t86 * rSges(3,2);
t67 = t86 * rSges(3,1) + rSges(3,2) * t87;
t66 = Icges(3,1) * t87 - t98;
t64 = -Icges(3,2) * t86 + t82;
t63 = Icges(3,2) * t87 + t98;
t60 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t59 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t58 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t57 = V_base(6) * rSges(2,1) + V_base(4) * rSges(2,2) + t94;
t56 = rSges(4,1) * t77 - rSges(4,2) * t76;
t55 = rSges(4,1) * t76 + rSges(4,2) * t77;
t54 = Icges(4,1) * t77 - t97;
t52 = -Icges(4,2) * t76 + t72;
t51 = Icges(4,2) * t77 + t97;
t48 = -V_base(5) * rSges(2,1) + V_base(4) * rSges(2,3) + t95;
t47 = -V_base(5) * rSges(2,2) + V_base(1) + (-rSges(2,3) - qJ(1)) * V_base(6);
t46 = rSges(5,1) * t74 - rSges(5,2) * t73;
t45 = rSges(5,1) * t73 + rSges(5,2) * t74;
t44 = Icges(5,1) * t74 - t96;
t42 = -Icges(5,2) * t73 + t70;
t41 = Icges(5,2) * t74 + t96;
t38 = V_base(5) * rSges(3,3) - t67 * t79 + t89;
t37 = t68 * t79 + (-rSges(3,3) - pkin(4)) * V_base(4) + t92;
t36 = t67 * V_base(4) + (-pkin(1) - t68) * V_base(5) + t95;
t35 = V_base(5) * rSges(4,3) + (-t55 - t104) * t79 + t88;
t34 = t56 * t79 + (-rSges(4,3) + t99) * V_base(4) + t91;
t33 = V_base(4) * t55 + (-t56 + t93) * V_base(5) + t90;
t32 = -t45 * t75 + t100 * V_base(5) + (-t102 - t104) * t79 + t88;
t31 = t79 * t101 + t46 * t75 + (t99 - t100) * V_base(4) + t91;
t30 = (t45 + t102) * V_base(4) + (-t46 + t93 - t101) * V_base(5) + t90;
t1 = m(1) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(2) * (t47 ^ 2 + t48 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + m(4) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) / 0.2e1 + m(5) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + Icges(5,3) * t75 ^ 2 / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,2) / 0.2e1) * V_base(6) ^ 2 + (Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1) * t79 ^ 2 + ((Icges(5,5) * t73 + Icges(5,6) * t74) * t75 + (Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + t63 * t105 + t86 * t108 + t51 * t106 + t76 * t109 + t41 * t107 + t73 * t110) * V_base(5) + (Icges(1,6) - Icges(2,6)) * V_base(6) + (Icges(3,5) * t86 + Icges(4,5) * t76 + Icges(3,6) * t87 + Icges(4,6) * t77) * t79) * V_base(5) + ((Icges(5,5) * t74 - Icges(5,6) * t73) * t75 + (Icges(1,1) / 0.2e1 + Icges(2,1) / 0.2e1 - t86 * t64 / 0.2e1 + t66 * t105 - t76 * t52 / 0.2e1 + t54 * t106 - t73 * t42 / 0.2e1 + t44 * t107) * V_base(4) + (Icges(1,5) - Icges(2,4)) * V_base(6) + (Icges(1,4) + Icges(2,5) + (t108 + t64 / 0.2e1) * t87 + (-t63 / 0.2e1 + t66 / 0.2e1) * t86 + (t109 + t52 / 0.2e1) * t77 + (-t51 / 0.2e1 + t54 / 0.2e1) * t76 + (t110 + t42 / 0.2e1) * t74 + (-t41 / 0.2e1 + t44 / 0.2e1) * t73) * V_base(5) + (Icges(3,5) * t87 + Icges(4,5) * t77 - Icges(3,6) * t86 - Icges(4,6) * t76) * t79) * V_base(4);
T  = t1;
