% Calculate kinetic energy for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
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
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2RR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_energykin_floatb_twist_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_energykin_floatb_twist_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_energykin_floatb_twist_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 17:59:54
% EndTime: 2019-03-08 17:59:54
% DurationCPUTime: 0.35s
% Computational Cost: add. (168->92), mult. (294->140), div. (0->0), fcn. (212->4), ass. (0->42)
t46 = sin(qJ(1));
t60 = Icges(2,4) * t46;
t48 = cos(qJ(1));
t59 = Icges(2,4) * t48;
t45 = sin(qJ(2));
t58 = Icges(3,4) * t45;
t47 = cos(qJ(2));
t57 = Icges(3,4) * t47;
t54 = -rSges(3,1) * t47 + rSges(3,2) * t45;
t53 = -Icges(3,1) * t47 + t58;
t52 = Icges(3,2) * t45 - t57;
t51 = -Icges(3,5) * t47 + Icges(3,6) * t45;
t42 = qJD(2) * t46 + V_base(6);
t43 = -qJD(2) * t48 + V_base(4);
t44 = V_base(5) + qJD(1);
t50 = (-Icges(3,3) * t48 + t51 * t46) * t43 + (Icges(3,3) * t46 + t51 * t48) * t42 + (-Icges(3,5) * t45 - Icges(3,6) * t47) * t44;
t21 = -Icges(3,6) * t48 + t52 * t46;
t22 = Icges(3,6) * t46 + t52 * t48;
t23 = -Icges(3,5) * t48 + t53 * t46;
t24 = Icges(3,5) * t46 + t53 * t48;
t33 = -Icges(3,2) * t47 - t58;
t36 = -Icges(3,1) * t45 - t57;
t49 = (t21 * t45 - t23 * t47) * t43 + (t33 * t45 - t36 * t47) * t44 + (t22 * t45 - t24 * t47) * t42;
t41 = -rSges(2,1) * t48 + t46 * rSges(2,2);
t40 = -t46 * rSges(2,1) - rSges(2,2) * t48;
t39 = -rSges(3,1) * t45 - rSges(3,2) * t47;
t38 = -Icges(2,1) * t48 + t60;
t37 = -Icges(2,1) * t46 - t59;
t35 = Icges(2,2) * t46 - t59;
t34 = -Icges(2,2) * t48 - t60;
t29 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t28 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t27 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t26 = t46 * rSges(3,3) + t54 * t48;
t25 = -t48 * rSges(3,3) + t54 * t46;
t18 = -V_base(6) * rSges(2,3) + t41 * t44 + V_base(1);
t17 = V_base(4) * rSges(2,3) - t40 * t44 + V_base(3);
t16 = t40 * V_base(6) - t41 * V_base(4) + V_base(2);
t15 = -t39 * t42 + V_base(1) + (pkin(1) * t46 + t26) * t44;
t14 = t43 * t39 + V_base(3) + (pkin(1) * t48 - t25) * t44;
t13 = t42 * t25 - t43 * t26 + V_base(2) + (-t46 * V_base(4) - t48 * V_base(6)) * pkin(1);
t1 = m(1) * (t27 ^ 2 + t28 ^ 2 + t29 ^ 2) / 0.2e1 + m(2) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(3) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + t43 * (t49 * t46 - t50 * t48) / 0.2e1 + t42 * (t50 * t46 + t49 * t48) / 0.2e1 + ((-t21 * t47 - t23 * t45) * t43 + (-t22 * t47 - t24 * t45) * t42 + (-t47 * t33 - t45 * t36 + Icges(2,3)) * t44) * t44 / 0.2e1 + ((-t48 * t35 - t46 * t38 + Icges(1,5)) * V_base(6) + (-t48 * t34 - t46 * t37 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t46 * t35 - t48 * t38 + Icges(1,3)) * V_base(6) + (t46 * t34 - t48 * t37 + Icges(1,5)) * V_base(4)) * V_base(6) / 0.2e1 + t44 * V_base(6) * (-Icges(2,5) * t48 + Icges(2,6) * t46) + V_base(4) * t44 * (-Icges(2,5) * t46 - Icges(2,6) * t48) + (Icges(1,4) * V_base(4) + Icges(1,6) * V_base(6) + Icges(1,2) * V_base(5) / 0.2e1) * V_base(5);
T  = t1;
