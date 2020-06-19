% Calculate kinetic energy for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2RR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_energykin_floatb_twist_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_energykin_floatb_twist_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR3_energykin_floatb_twist_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:21
% EndTime: 2020-06-19 09:14:22
% DurationCPUTime: 0.46s
% Computational Cost: add. (163->81), mult. (172->110), div. (0->0), fcn. (84->4), ass. (0->38)
t44 = qJ(1) + qJ(2);
t41 = cos(t44);
t54 = t41 / 0.2e1;
t46 = cos(qJ(1));
t53 = t46 / 0.2e1;
t39 = V_base(6) + qJD(1);
t52 = pkin(1) * t39;
t51 = pkin(3) + rSges(3,3);
t45 = sin(qJ(1));
t50 = Icges(2,4) * t45;
t40 = sin(t44);
t49 = Icges(3,4) * t40;
t48 = V_base(5) * pkin(2) + V_base(1);
t42 = Icges(2,4) * t46;
t38 = qJD(2) + t39;
t37 = Icges(3,4) * t41;
t36 = t46 * rSges(2,1) - t45 * rSges(2,2);
t35 = t45 * rSges(2,1) + t46 * rSges(2,2);
t34 = Icges(2,1) * t46 - t50;
t33 = Icges(2,1) * t45 + t42;
t32 = -Icges(2,2) * t45 + t42;
t31 = Icges(2,2) * t46 + t50;
t28 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t27 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t26 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t25 = t41 * rSges(3,1) - t40 * rSges(3,2);
t24 = t40 * rSges(3,1) + t41 * rSges(3,2);
t23 = Icges(3,1) * t41 - t49;
t22 = Icges(3,1) * t40 + t37;
t21 = -Icges(3,2) * t40 + t37;
t20 = Icges(3,2) * t41 + t49;
t17 = V_base(5) * rSges(2,3) - t39 * t35 + t48;
t16 = t39 * t36 + V_base(2) + (-pkin(2) - rSges(2,3)) * V_base(4);
t15 = V_base(4) * t35 - V_base(5) * t36 + V_base(3);
t14 = -t38 * t24 - t45 * t52 + t51 * V_base(5) + t48;
t13 = t46 * t52 + t38 * t25 + V_base(2) + (-pkin(2) - t51) * V_base(4);
t12 = V_base(4) * t24 - V_base(5) * t25 + V_base(3) + (t45 * V_base(4) - t46 * V_base(5)) * pkin(1);
t1 = m(1) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + m(2) * (t15 ^ 2 + t16 ^ 2 + t17 ^ 2) / 0.2e1 + Icges(2,3) * t39 ^ 2 / 0.2e1 + m(3) * (t12 ^ 2 + t13 ^ 2 + t14 ^ 2) / 0.2e1 + Icges(3,3) * t38 ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (Icges(3,5) * t40 + Icges(3,6) * t41) * t38 + (Icges(2,5) * t45 + Icges(2,6) * t46) * t39 + (Icges(1,2) / 0.2e1 + t31 * t53 + t45 * t33 / 0.2e1 + t20 * t54 + t40 * t22 / 0.2e1) * V_base(5)) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t46 - Icges(2,6) * t45) * t39 + (Icges(3,5) * t41 - Icges(3,6) * t40) * t38 + ((t32 + t33) * t46 + (-t31 + t34) * t45 + (t21 + t22) * t41 + (-t20 + t23) * t40) * V_base(5) / 0.2e1 + (Icges(1,1) / 0.2e1 - t45 * t32 / 0.2e1 + t34 * t53 - t40 * t21 / 0.2e1 + t23 * t54) * V_base(4)) * V_base(4);
T = t1;
