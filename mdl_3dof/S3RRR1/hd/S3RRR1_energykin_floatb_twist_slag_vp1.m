% Calculate kinetic energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3RRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRR1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:47
% EndTime: 2019-03-08 18:07:47
% DurationCPUTime: 0.34s
% Computational Cost: add. (323->119), mult. (260->153), div. (0->0), fcn. (132->6), ass. (0->61)
t77 = qJ(1) + qJ(2);
t74 = qJ(3) + t77;
t69 = cos(t74);
t96 = t69 / 0.2e1;
t72 = cos(t77);
t95 = t72 / 0.2e1;
t79 = cos(qJ(1));
t94 = t79 / 0.2e1;
t93 = -pkin(3) - pkin(4);
t78 = sin(qJ(1));
t92 = pkin(1) * t78;
t91 = pkin(1) * t79;
t71 = sin(t77);
t90 = pkin(2) * t71;
t89 = pkin(2) * t72;
t88 = pkin(5) + rSges(4,3);
t87 = Icges(2,4) * t78;
t86 = Icges(3,4) * t71;
t68 = sin(t74);
t85 = Icges(4,4) * t68;
t70 = V_base(6) + qJD(1);
t84 = t70 * t91 + V_base(2);
t83 = V_base(4) * t92 + V_base(3);
t82 = V_base(5) * pkin(3) + V_base(1);
t67 = qJD(2) + t70;
t80 = V_base(5) * pkin(4) - t70 * t92 + t82;
t73 = Icges(2,4) * t79;
t66 = Icges(3,4) * t72;
t64 = qJD(3) + t67;
t63 = Icges(4,4) * t69;
t61 = t79 * rSges(2,1) - t78 * rSges(2,2);
t60 = t78 * rSges(2,1) + t79 * rSges(2,2);
t59 = Icges(2,1) * t79 - t87;
t58 = Icges(2,1) * t78 + t73;
t57 = -Icges(2,2) * t78 + t73;
t56 = Icges(2,2) * t79 + t87;
t53 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t52 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t51 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t50 = t72 * rSges(3,1) - t71 * rSges(3,2);
t49 = t71 * rSges(3,1) + t72 * rSges(3,2);
t48 = Icges(3,1) * t72 - t86;
t47 = Icges(3,1) * t71 + t66;
t46 = -Icges(3,2) * t71 + t66;
t45 = Icges(3,2) * t72 + t86;
t42 = t69 * rSges(4,1) - t68 * rSges(4,2);
t41 = t68 * rSges(4,1) + t69 * rSges(4,2);
t40 = Icges(4,1) * t69 - t85;
t39 = Icges(4,1) * t68 + t63;
t38 = -Icges(4,2) * t68 + t63;
t37 = Icges(4,2) * t69 + t85;
t34 = V_base(5) * rSges(2,3) - t70 * t60 + t82;
t33 = t70 * t61 + V_base(2) + (-pkin(3) - rSges(2,3)) * V_base(4);
t32 = V_base(4) * t60 - V_base(5) * t61 + V_base(3);
t31 = V_base(5) * rSges(3,3) - t67 * t49 + t80;
t30 = t67 * t50 + (-rSges(3,3) + t93) * V_base(4) + t84;
t29 = V_base(4) * t49 + (-t50 - t91) * V_base(5) + t83;
t28 = -t64 * t41 - t67 * t90 + t88 * V_base(5) + t80;
t27 = t67 * t89 + t64 * t42 + (-t88 + t93) * V_base(4) + t84;
t26 = (t41 + t90) * V_base(4) + (-t42 - t89 - t91) * V_base(5) + t83;
t1 = m(1) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + m(2) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + Icges(2,3) * t70 ^ 2 / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + Icges(3,3) * t67 ^ 2 / 0.2e1 + m(4) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + Icges(4,3) * t64 ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (Icges(4,5) * t68 + Icges(4,6) * t69) * t64 + (Icges(3,5) * t71 + Icges(3,6) * t72) * t67 + (Icges(2,5) * t78 + Icges(2,6) * t79) * t70 + (Icges(1,2) / 0.2e1 + t56 * t94 + t78 * t58 / 0.2e1 + t45 * t95 + t71 * t47 / 0.2e1 + t37 * t96 + t68 * t39 / 0.2e1) * V_base(5)) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t79 - Icges(2,6) * t78) * t70 + (Icges(3,5) * t72 - Icges(3,6) * t71) * t67 + (Icges(4,5) * t69 - Icges(4,6) * t68) * t64 + ((t57 + t58) * t79 + (-t56 + t59) * t78 + (t46 + t47) * t72 + (-t45 + t48) * t71 + (t38 + t39) * t69 + (-t37 + t40) * t68) * V_base(5) / 0.2e1 + (Icges(1,1) / 0.2e1 - t78 * t57 / 0.2e1 + t59 * t94 - t71 * t46 / 0.2e1 + t48 * t95 - t68 * t38 / 0.2e1 + t40 * t96) * V_base(4)) * V_base(4);
T  = t1;
