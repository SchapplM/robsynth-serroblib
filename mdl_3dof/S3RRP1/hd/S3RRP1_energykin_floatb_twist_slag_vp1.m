% Calculate kinetic energy for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:06
% EndTime: 2018-11-14 10:15:07
% DurationCPUTime: 0.52s
% Computational Cost: add. (287->105), mult. (264->133), div. (0->0), fcn. (138->4), ass. (0->52)
t102 = Icges(3,4) - Icges(4,5);
t101 = Icges(3,1) + Icges(4,1);
t100 = Icges(3,2) + Icges(4,3);
t73 = qJ(1) + qJ(2);
t68 = sin(t73);
t99 = t102 * t68;
t69 = cos(t73);
t98 = t102 * t69;
t97 = rSges(4,1) + pkin(2);
t96 = -t100 * t69 - t99;
t95 = t100 * t68 - t98;
t94 = t101 * t68 + t98;
t93 = t101 * t69 - t99;
t92 = rSges(4,3) + qJ(3);
t91 = Icges(4,4) + Icges(3,5);
t90 = Icges(3,6) - Icges(4,6);
t89 = -pkin(3) - pkin(4);
t74 = sin(qJ(1));
t88 = pkin(1) * t74;
t75 = cos(qJ(1));
t87 = pkin(1) * t75;
t86 = t97 * t68 - t92 * t69;
t85 = t92 * t68 + t97 * t69;
t84 = Icges(2,4) * t74;
t67 = V_base(6) + qJD(1);
t81 = t67 * t87 + V_base(2);
t80 = V_base(4) * t88 + V_base(3);
t79 = V_base(5) * pkin(3) + V_base(1);
t76 = V_base(5) * pkin(4) - t67 * t88 + t79;
t70 = Icges(2,4) * t75;
t66 = qJD(2) + t67;
t61 = rSges(2,1) * t75 - t74 * rSges(2,2);
t60 = t74 * rSges(2,1) + rSges(2,2) * t75;
t59 = Icges(2,1) * t75 - t84;
t58 = Icges(2,1) * t74 + t70;
t57 = -Icges(2,2) * t74 + t70;
t56 = Icges(2,2) * t75 + t84;
t53 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t52 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t51 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t50 = rSges(3,1) * t69 - rSges(3,2) * t68;
t47 = rSges(3,1) * t68 + rSges(3,2) * t69;
t32 = V_base(5) * rSges(2,3) - t60 * t67 + t79;
t31 = t61 * t67 + V_base(2) + (-rSges(2,3) - pkin(3)) * V_base(4);
t30 = t60 * V_base(4) - t61 * V_base(5) + V_base(3);
t29 = V_base(5) * rSges(3,3) - t47 * t66 + t76;
t28 = t50 * t66 + (-rSges(3,3) + t89) * V_base(4) + t81;
t27 = V_base(4) * t47 + (-t50 - t87) * V_base(5) + t80;
t26 = V_base(5) * rSges(4,2) + qJD(3) * t68 - t86 * t66 + t76;
t25 = -qJD(3) * t69 + t85 * t66 + (-rSges(4,2) + t89) * V_base(4) + t81;
t24 = t86 * V_base(4) + (-t85 - t87) * V_base(5) + t80;
t1 = m(1) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t27 ^ 2 + t28 ^ 2 + t29 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t74 + Icges(2,6) * t75) * V_base(5) + (Icges(2,5) * t75 - Icges(2,6) * t74) * V_base(4) + Icges(2,3) * t67 / 0.2e1) * t67 + ((Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t66 + (t90 * V_base(5) + t91 * V_base(4)) * t69 + (-t90 * V_base(4) + t91 * V_base(5)) * t68) * t66 + ((-t74 * t56 + t75 * t58 + t96 * t68 + t94 * t69) * V_base(5) + (-t74 * t57 + t75 * t59 + t95 * t68 + t93 * t69) * V_base(4)) * V_base(4) / 0.2e1 + ((t75 * t56 + t74 * t58 + t94 * t68 - t96 * t69) * V_base(5) + (t75 * t57 + t74 * t59 + t93 * t68 - t95 * t69) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
