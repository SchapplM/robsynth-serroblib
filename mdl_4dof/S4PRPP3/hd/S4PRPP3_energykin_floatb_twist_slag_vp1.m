% Calculate kinetic energy for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:01
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:08
% EndTime: 2018-11-14 14:01:09
% DurationCPUTime: 0.37s
% Computational Cost: add. (224->114), mult. (308->126), div. (0->0), fcn. (144->2), ass. (0->48)
t107 = Icges(3,4) / 0.2e1 - Icges(4,5) / 0.2e1 - Icges(5,4) / 0.2e1;
t103 = rSges(5,1) + pkin(3);
t102 = Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1 + Icges(3,2) / 0.2e1;
t101 = Icges(3,1) / 0.2e1 + Icges(5,1) / 0.2e1 + Icges(4,1) / 0.2e1;
t94 = Icges(4,4) + Icges(3,5) - Icges(5,5);
t93 = Icges(3,6) - Icges(4,6) + Icges(5,6);
t71 = sin(qJ(2));
t92 = t107 * t71;
t72 = cos(qJ(2));
t91 = t107 * t72;
t59 = pkin(2) * t72 + t71 * qJ(3);
t90 = -pkin(1) - t59;
t86 = -qJ(4) - rSges(5,3);
t85 = V_base(4) * qJ(1) + V_base(3);
t84 = qJD(1) + V_base(2);
t83 = -rSges(5,2) * t72 + t103 * t71;
t82 = t71 * rSges(5,2) + t103 * t72;
t55 = t71 * pkin(2) - qJ(3) * t72;
t81 = V_base(4) * t55 + t85;
t80 = V_base(6) * pkin(1) + t84;
t79 = t102 * t72 + t92;
t78 = t102 * t71 - t91;
t77 = t101 * t71 + t91;
t76 = t101 * t72 - t92;
t75 = V_base(5) * pkin(4) - V_base(6) * qJ(1) + V_base(1);
t74 = qJD(3) * t71 + t75;
t63 = V_base(6) + qJD(2);
t73 = -qJD(3) * t72 + t63 * t59 + t80;
t62 = rSges(3,1) * t72 - t71 * rSges(3,2);
t61 = rSges(4,1) * t72 + t71 * rSges(4,3);
t58 = t71 * rSges(3,1) + rSges(3,2) * t72;
t57 = t71 * rSges(4,1) - rSges(4,3) * t72;
t36 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t35 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t34 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t32 = V_base(6) * rSges(2,1) + V_base(4) * rSges(2,2) + t84;
t30 = -V_base(5) * rSges(2,1) + V_base(4) * rSges(2,3) + t85;
t29 = -V_base(5) * rSges(2,2) + V_base(1) + (-rSges(2,3) - qJ(1)) * V_base(6);
t28 = V_base(5) * rSges(3,3) - t58 * t63 + t75;
t27 = t62 * t63 + (-rSges(3,3) - pkin(4)) * V_base(4) + t80;
t26 = t58 * V_base(4) + (-pkin(1) - t62) * V_base(5) + t85;
t25 = V_base(5) * rSges(4,2) + (-t55 - t57) * t63 + t74;
t24 = t63 * t61 + (-rSges(4,2) - pkin(4)) * V_base(4) + t73;
t23 = t57 * V_base(4) + (-t61 + t90) * V_base(5) + t81;
t22 = t86 * V_base(5) + (-t55 - t83) * t63 + t74;
t21 = t82 * t63 + (-pkin(4) - t86) * V_base(4) + t73;
t20 = -qJD(4) + t83 * V_base(4) + (-t82 + t90) * V_base(5) + t81;
t1 = m(1) * (t34 ^ 2 + t35 ^ 2 + t36 ^ 2) / 0.2e1 + m(2) * (t29 ^ 2 + t30 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + m(5) * (t20 ^ 2 + t21 ^ 2 + t22 ^ 2) / 0.2e1 + (Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + t79 * t72 + t77 * t71) * V_base(5) ^ 2 + ((Icges(1,1) / 0.2e1 + Icges(2,1) / 0.2e1 + t76 * t72 + t78 * t71) * V_base(4) + (Icges(1,4) + Icges(2,5) + (t77 - t78) * t72 + (t76 - t79) * t71) * V_base(5)) * V_base(4) + ((Icges(1,3) / 0.2e1 + Icges(2,2) / 0.2e1) * V_base(6) + (Icges(1,6) - Icges(2,6)) * V_base(5) + (-Icges(2,4) + Icges(1,5)) * V_base(4)) * V_base(6) + ((Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1 + Icges(5,3) / 0.2e1) * t63 + (t93 * V_base(5) + t94 * V_base(4)) * t72 + (-t93 * V_base(4) + t94 * V_base(5)) * t71) * t63;
T  = t1;
