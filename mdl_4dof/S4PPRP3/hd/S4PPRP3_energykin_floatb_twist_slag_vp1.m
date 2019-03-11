% Calculate kinetic energy for
% S4PPRP3
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:19
% EndTime: 2019-03-08 18:14:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (185->111), mult. (248->113), div. (0->0), fcn. (84->2), ass. (0->44)
t85 = Icges(4,4) + Icges(5,4);
t82 = rSges(5,1) + pkin(3);
t81 = Icges(5,2) / 0.2e1 + Icges(4,2) / 0.2e1;
t80 = Icges(5,1) / 0.2e1 + Icges(4,1) / 0.2e1;
t75 = -Icges(4,5) - Icges(5,5);
t74 = Icges(4,6) + Icges(5,6);
t54 = cos(qJ(3));
t73 = t85 * t54 / 0.2e1;
t53 = sin(qJ(3));
t72 = -t85 * t53 / 0.2e1;
t71 = -pkin(1) - pkin(4);
t68 = rSges(5,3) + qJ(4);
t67 = V_base(4) * qJ(1) + V_base(3);
t66 = qJD(1) + V_base(2);
t65 = qJD(2) + V_base(1);
t64 = t81 * t54 - t72;
t63 = t81 * t53 - t73;
t62 = t80 * t53 + t73;
t61 = t80 * t54 + t72;
t60 = t54 * rSges(5,2) + t82 * t53;
t59 = -t53 * rSges(5,2) + t82 * t54;
t58 = V_base(4) * pkin(2) + t67;
t57 = V_base(4) * pkin(1) + V_base(6) * qJ(2) + t66;
t56 = V_base(4) * pkin(4) + t57;
t55 = (-pkin(2) - qJ(1)) * V_base(6) + t65;
t45 = V_base(6) - qJD(3);
t44 = t54 * rSges(4,1) - t53 * rSges(4,2);
t42 = t53 * rSges(4,1) + t54 * rSges(4,2);
t28 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t27 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t26 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t25 = V_base(4) * rSges(2,1) - V_base(6) * rSges(2,2) + t66;
t24 = -V_base(5) * rSges(2,1) + V_base(1) + (-qJ(1) - rSges(2,3)) * V_base(6);
t23 = V_base(5) * rSges(2,2) + V_base(4) * rSges(2,3) + t67;
t22 = V_base(4) * rSges(3,1) + (-qJ(2) - rSges(3,3)) * V_base(5) + t67;
t21 = (-qJ(1) - rSges(3,1)) * V_base(6) + (-pkin(1) + rSges(3,2)) * V_base(5) + t65;
t20 = -V_base(4) * rSges(3,2) + V_base(6) * rSges(3,3) + t57;
t19 = V_base(4) * rSges(4,3) + t45 * t42 + t56;
t18 = -t45 * t44 + (-rSges(4,3) + t71) * V_base(5) + t55;
t17 = V_base(4) * t44 + (-qJ(2) - t42) * V_base(5) + t58;
t16 = t60 * t45 + t68 * V_base(4) + t56;
t15 = -t59 * t45 + (-t68 + t71) * V_base(5) + t55;
t14 = -qJD(4) + t59 * V_base(4) + (-qJ(2) - t60) * V_base(5) + t58;
t1 = m(1) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(2) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t22 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t15 ^ 2 + t16 ^ 2) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t45 ^ 2 + (Icges(1,3) / 0.2e1 + Icges(2,1) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) ^ 2 + ((Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 + t61 * t54 + t63 * t53) * V_base(5) + (t74 * t53 + t75 * t54) * t45 + (Icges(3,4) - Icges(2,5) + Icges(1,6)) * V_base(6)) * V_base(5) + ((Icges(1,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,3) / 0.2e1 + t64 * t54 + t62 * t53) * V_base(4) + (t75 * t53 - t74 * t54) * t45 + (Icges(1,5) + Icges(2,4) + Icges(3,6)) * V_base(6) + (Icges(1,4) - Icges(2,6) + Icges(3,5) + (t62 - t63) * t54 + (t61 - t64) * t53) * V_base(5)) * V_base(4);
T  = t1;
