% Calculate kinetic energy for
% S4PPRP4
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
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:52
% EndTime: 2018-11-14 14:05:53
% DurationCPUTime: 0.29s
% Computational Cost: add. (185->112), mult. (248->112), div. (0->0), fcn. (84->2), ass. (0->43)
t84 = Icges(4,4) + Icges(5,4);
t81 = rSges(5,1) + pkin(3);
t80 = Icges(5,2) / 0.2e1 + Icges(4,2) / 0.2e1;
t79 = Icges(5,1) / 0.2e1 + Icges(4,1) / 0.2e1;
t74 = Icges(4,5) + Icges(5,5);
t73 = Icges(4,6) + Icges(5,6);
t52 = cos(qJ(3));
t72 = t84 * t52 / 0.2e1;
t51 = sin(qJ(3));
t71 = -t84 * t51 / 0.2e1;
t70 = -pkin(1) - pkin(4);
t69 = -pkin(2) - qJ(1);
t64 = rSges(5,3) + qJ(4);
t63 = V_base(6) * qJ(1) + V_base(2);
t62 = qJD(1) + V_base(1);
t61 = t52 * t80 - t71;
t60 = t51 * t80 - t72;
t59 = t51 * t79 + t72;
t58 = t52 * t79 + t71;
t57 = t52 * rSges(5,2) + t51 * t81;
t56 = -t51 * rSges(5,2) + t52 * t81;
t55 = V_base(6) * qJ(2) + t62;
t54 = V_base(4) * pkin(1) - qJD(2) + t63;
t53 = V_base(6) * pkin(2) + V_base(4) * pkin(4) + t54;
t45 = V_base(6) - qJD(3);
t44 = rSges(4,1) * t52 - rSges(4,2) * t51;
t42 = -rSges(4,1) * t51 - rSges(4,2) * t52;
t28 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t27 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t26 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t25 = -rSges(2,1) * V_base(5) - rSges(2,2) * V_base(6) + t62;
t24 = rSges(2,1) * V_base(4) + rSges(2,3) * V_base(6) + t63;
t23 = V_base(4) * rSges(2,2) + V_base(3) + (-qJ(1) - rSges(2,3)) * V_base(5);
t22 = V_base(3) + (-qJ(1) - rSges(3,1)) * V_base(5) + (-qJ(2) - rSges(3,3)) * V_base(4);
t21 = rSges(3,1) * V_base(6) - rSges(3,2) * V_base(4) + t54;
t20 = V_base(6) * rSges(3,3) + (-pkin(1) + rSges(3,2)) * V_base(5) + t55;
t19 = -t45 * t42 + (-rSges(4,3) + t70) * V_base(5) + t55;
t18 = V_base(4) * rSges(4,3) + t44 * t45 + t53;
t17 = V_base(3) + (-qJ(2) + t42) * V_base(4) + (-t44 + t69) * V_base(5);
t16 = t57 * t45 + (-t64 + t70) * V_base(5) + t55;
t15 = t45 * t56 + t64 * V_base(4) + t53;
t14 = -qJD(4) + V_base(3) + (-qJ(2) - t57) * V_base(4) + (-t56 + t69) * V_base(5);
t1 = m(1) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(2) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t22 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t15 ^ 2 + t16 ^ 2) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t45 ^ 2 + (Icges(1,3) / 0.2e1 + Icges(2,1) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) ^ 2 + ((Icges(1,2) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,3) / 0.2e1 + t61 * t52 + t59 * t51) * V_base(5) + (t51 * t74 + t52 * t73) * t45 + (-Icges(2,4) + Icges(1,6) - Icges(3,6)) * V_base(6)) * V_base(5) + ((Icges(1,1) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 + t58 * t52 + t60 * t51) * V_base(4) + (t51 * t73 - t52 * t74) * t45 + (Icges(1,5) - Icges(2,5) + Icges(3,4)) * V_base(6) + (Icges(1,4) + Icges(2,6) - Icges(3,5) + (-t59 + t60) * t52 + (-t58 + t61) * t51) * V_base(5)) * V_base(4);
T  = t1;
