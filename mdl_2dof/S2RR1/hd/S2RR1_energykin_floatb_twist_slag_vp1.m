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
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:19:05
% EndTime: 2020-01-03 11:19:06
% DurationCPUTime: 0.34s
% Computational Cost: add. (168->92), mult. (294->140), div. (0->0), fcn. (212->4), ass. (0->42)
t47 = sin(qJ(1));
t60 = Icges(2,4) * t47;
t46 = sin(qJ(2));
t59 = Icges(3,4) * t46;
t48 = cos(qJ(2));
t58 = Icges(3,4) * t48;
t55 = rSges(3,1) * t48 - rSges(3,2) * t46;
t54 = Icges(3,1) * t48 - t59;
t53 = -Icges(3,2) * t46 + t58;
t52 = Icges(3,5) * t48 - Icges(3,6) * t46;
t42 = -qJD(2) * t47 + V_base(6);
t49 = cos(qJ(1));
t43 = qJD(2) * t49 + V_base(4);
t44 = V_base(5) + qJD(1);
t51 = (Icges(3,3) * t49 + t52 * t47) * t43 + (-Icges(3,3) * t47 + t52 * t49) * t42 + (-Icges(3,5) * t46 - Icges(3,6) * t48) * t44;
t21 = Icges(3,6) * t49 + t53 * t47;
t22 = -Icges(3,6) * t47 + t53 * t49;
t23 = Icges(3,5) * t49 + t54 * t47;
t24 = -Icges(3,5) * t47 + t54 * t49;
t33 = -Icges(3,2) * t48 - t59;
t36 = -Icges(3,1) * t46 - t58;
t50 = (-t21 * t46 + t23 * t48) * t43 + (-t33 * t46 + t36 * t48) * t44 + (-t22 * t46 + t24 * t48) * t42;
t45 = Icges(2,4) * t49;
t41 = t49 * rSges(2,1) - t47 * rSges(2,2);
t40 = t47 * rSges(2,1) + t49 * rSges(2,2);
t39 = -t46 * rSges(3,1) - t48 * rSges(3,2);
t38 = Icges(2,1) * t49 - t60;
t37 = Icges(2,1) * t47 + t45;
t35 = -Icges(2,2) * t47 + t45;
t34 = Icges(2,2) * t49 + t60;
t29 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t28 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t27 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t26 = -t47 * rSges(3,3) + t55 * t49;
t25 = t49 * rSges(3,3) + t55 * t47;
t18 = -V_base(6) * rSges(2,3) + t44 * t41 + V_base(1);
t17 = V_base(4) * rSges(2,3) - t44 * t40 + V_base(3);
t16 = V_base(6) * t40 - V_base(4) * t41 + V_base(2);
t15 = -t42 * t39 + V_base(1) + (-pkin(1) * t47 + t26) * t44;
t14 = t43 * t39 + V_base(3) + (-pkin(1) * t49 - t25) * t44;
t13 = t42 * t25 - t43 * t26 + V_base(2) + (t47 * V_base(4) + t49 * V_base(6)) * pkin(1);
t1 = m(1) * (t27 ^ 2 + t28 ^ 2 + t29 ^ 2) / 0.2e1 + m(2) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(3) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + t43 * (t50 * t47 + t51 * t49) / 0.2e1 + t42 * (-t51 * t47 + t50 * t49) / 0.2e1 + ((-t48 * t21 - t46 * t23) * t43 + (-t48 * t22 - t46 * t24) * t42 + (-t48 * t33 - t46 * t36 + Icges(2,3)) * t44) * t44 / 0.2e1 + ((t49 * t35 + t47 * t38 + Icges(1,5)) * V_base(6) + (t49 * t34 + t47 * t37 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t47 * t35 + t49 * t38 + Icges(1,3)) * V_base(6) + (-t47 * t34 + t49 * t37 + Icges(1,5)) * V_base(4)) * V_base(6) / 0.2e1 + t44 * V_base(6) * (Icges(2,5) * t49 - Icges(2,6) * t47) + V_base(4) * t44 * (Icges(2,5) * t47 + Icges(2,6) * t49) + (Icges(1,4) * V_base(4) + Icges(1,6) * V_base(6) + Icges(1,2) * V_base(5) / 0.2e1) * V_base(5);
T = t1;
