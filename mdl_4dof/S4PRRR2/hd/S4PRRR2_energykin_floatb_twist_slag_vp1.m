% Calculate kinetic energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:20
% EndTime: 2019-07-18 13:27:21
% DurationCPUTime: 0.40s
% Computational Cost: add. (341->134), mult. (282->175), div. (0->0), fcn. (132->6), ass. (0->65)
t73 = qJ(2) + qJ(3);
t72 = qJ(4) + t73;
t67 = cos(t72);
t98 = t67 / 0.2e1;
t70 = cos(t73);
t97 = t70 / 0.2e1;
t75 = cos(qJ(2));
t96 = t75 / 0.2e1;
t66 = sin(t72);
t81 = Icges(5,4) * t66;
t95 = Icges(5,2) * t98 + t81 / 0.2e1;
t69 = sin(t73);
t83 = Icges(4,4) * t69;
t94 = Icges(4,2) * t97 + t83 / 0.2e1;
t74 = sin(qJ(2));
t85 = Icges(3,4) * t74;
t93 = Icges(3,2) * t96 + t85 / 0.2e1;
t92 = -t66 / 0.2e1;
t91 = -t69 / 0.2e1;
t90 = -t74 / 0.2e1;
t89 = pkin(1) * t74;
t88 = pkin(1) * t75;
t87 = pkin(2) * t69;
t86 = pkin(2) * t70;
t84 = Icges(3,4) * t75;
t82 = Icges(4,4) * t70;
t80 = Icges(5,4) * t67;
t68 = V_base(5) - qJD(2);
t79 = t68 * t89 + V_base(3);
t78 = V_base(6) * qJ(1) + V_base(1);
t77 = -qJD(1) + V_base(2);
t76 = t68 * t88 + t78;
t65 = -qJD(3) + t68;
t64 = -qJD(4) + t65;
t61 = rSges(3,1) * t75 - t74 * rSges(3,2);
t60 = -t74 * rSges(3,1) - rSges(3,2) * t75;
t59 = Icges(3,1) * t75 - t85;
t58 = -Icges(3,1) * t74 - t84;
t57 = -Icges(3,2) * t74 + t84;
t53 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t52 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t51 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t50 = -V_base(4) * rSges(2,1) - V_base(6) * rSges(2,2) + t77;
t49 = rSges(4,1) * t70 - rSges(4,2) * t69;
t48 = -rSges(4,1) * t69 - rSges(4,2) * t70;
t47 = Icges(4,1) * t70 - t83;
t46 = -Icges(4,1) * t69 - t82;
t45 = -Icges(4,2) * t69 + t82;
t41 = rSges(5,1) * t67 - rSges(5,2) * t66;
t40 = -rSges(5,1) * t66 - rSges(5,2) * t67;
t39 = V_base(5) * rSges(2,1) + V_base(6) * rSges(2,3) + t78;
t38 = V_base(5) * rSges(2,2) + V_base(3) + (-qJ(1) - rSges(2,3)) * V_base(4);
t37 = Icges(5,1) * t67 - t81;
t36 = -Icges(5,1) * t66 - t80;
t35 = -Icges(5,2) * t66 + t80;
t31 = V_base(6) * rSges(3,3) + t61 * t68 + t78;
t30 = -t60 * t68 + V_base(3) + (-qJ(1) - rSges(3,3)) * V_base(4);
t29 = t60 * V_base(6) - t61 * V_base(4) + t77;
t28 = V_base(6) * rSges(4,3) + t49 * t65 + t76;
t27 = -t48 * t65 + (-qJ(1) - rSges(4,3)) * V_base(4) + t79;
t26 = V_base(6) * t48 - V_base(4) * t49 + (-V_base(6) * t74 - t75 * V_base(4)) * pkin(1) + t77;
t25 = V_base(6) * rSges(5,3) + t41 * t64 + t65 * t86 + t76;
t24 = t65 * t87 - t40 * t64 + (-qJ(1) - rSges(5,3)) * V_base(4) + t79;
t23 = (t40 - t87 - t89) * V_base(6) + (-t41 - t86 - t88) * V_base(4) + t77;
t1 = m(1) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + m(2) * (t38 ^ 2 + t39 ^ 2 + t50 ^ 2) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + Icges(3,3) * t68 ^ 2 / 0.2e1 + m(4) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + Icges(4,3) * t65 ^ 2 / 0.2e1 + m(5) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + Icges(5,3) * t64 ^ 2 / 0.2e1 + (Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(5) ^ 2 + ((-Icges(5,5) * t67 + Icges(5,6) * t66) * t64 + (-Icges(4,5) * t70 + Icges(4,6) * t69) * t65 + (-Icges(3,5) * t75 + Icges(3,6) * t74) * t68 + (Icges(1,3) / 0.2e1 + Icges(2,1) / 0.2e1 + t57 * t90 + t59 * t96 + t45 * t91 + t47 * t97 + t35 * t92 + t37 * t98) * V_base(6) + (-Icges(2,5) + Icges(1,6)) * V_base(5)) * V_base(6) + ((Icges(3,5) * t74 + Icges(3,6) * t75) * t68 + (Icges(4,5) * t69 + Icges(4,6) * t70) * t65 + (Icges(5,5) * t66 + Icges(5,6) * t67) * t64 + (Icges(1,1) / 0.2e1 + Icges(2,2) / 0.2e1 + t75 * t93 + t58 * t90 + t70 * t94 + t46 * t91 + t67 * t95 + t36 * t92) * V_base(4) + (Icges(1,5) - Icges(2,4) + (-t57 / 0.2e1 + t58 / 0.2e1) * t75 + (-t59 / 0.2e1 + t93) * t74 + (-t45 / 0.2e1 + t46 / 0.2e1) * t70 + (-t47 / 0.2e1 + t94) * t69 + (-t35 / 0.2e1 + t36 / 0.2e1) * t67 + (-t37 / 0.2e1 + t95) * t66) * V_base(6) + (Icges(1,4) + Icges(2,6)) * V_base(5)) * V_base(4);
T  = t1;
