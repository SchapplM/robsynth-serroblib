% Calculate kinetic energy for
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2019-03-08 18:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:11:05
% EndTime: 2019-03-08 18:11:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (218->127), mult. (248->142), div. (0->0), fcn. (84->4), ass. (0->51)
t59 = pkin(5) + qJ(4);
t51 = cos(t59);
t49 = Icges(5,4) * t51;
t50 = sin(t59);
t78 = Icges(5,1) * t50 / 0.2e1 + t49 / 0.2e1;
t61 = cos(pkin(5));
t52 = Icges(4,4) * t61;
t60 = sin(pkin(5));
t77 = Icges(4,1) * t60 / 0.2e1 + t52 / 0.2e1;
t76 = t51 / 0.2e1;
t75 = t61 / 0.2e1;
t74 = pkin(3) * t60;
t73 = pkin(3) * t61;
t72 = rSges(5,3) + pkin(4);
t71 = -pkin(1) - qJ(3);
t70 = -pkin(2) - qJ(1);
t69 = Icges(4,4) * t60;
t68 = Icges(5,4) * t50;
t67 = V_base(4) * qJ(1) + V_base(3);
t66 = qJD(1) + V_base(2);
t65 = qJD(2) + V_base(1);
t64 = V_base(4) * pkin(1) + V_base(6) * qJ(2) + t66;
t63 = V_base(4) * pkin(2) - qJD(3) + t67;
t62 = V_base(4) * qJ(3) + t64;
t53 = V_base(6) - qJD(4);
t48 = t61 * rSges(4,1) - t60 * rSges(4,2);
t47 = t60 * rSges(4,1) + t61 * rSges(4,2);
t46 = Icges(4,1) * t61 - t69;
t44 = -Icges(4,2) * t60 + t52;
t43 = Icges(4,2) * t61 + t69;
t40 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t39 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t38 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t37 = V_base(4) * rSges(2,1) - V_base(6) * rSges(2,2) + t66;
t36 = t51 * rSges(5,1) - t50 * rSges(5,2);
t35 = t50 * rSges(5,1) + t51 * rSges(5,2);
t34 = Icges(5,1) * t51 - t68;
t32 = -Icges(5,2) * t50 + t49;
t31 = Icges(5,2) * t51 + t68;
t28 = -V_base(5) * rSges(2,1) + V_base(1) + (-qJ(1) - rSges(2,3)) * V_base(6);
t27 = V_base(5) * rSges(2,2) + V_base(4) * rSges(2,3) + t67;
t26 = V_base(4) * rSges(3,1) + (-qJ(2) - rSges(3,3)) * V_base(5) + t67;
t25 = (-qJ(1) - rSges(3,1)) * V_base(6) + (-pkin(1) + rSges(3,2)) * V_base(5) + t65;
t24 = -V_base(4) * rSges(3,2) + V_base(6) * rSges(3,3) + t64;
t23 = V_base(4) * rSges(4,3) + V_base(6) * t47 + t62;
t22 = (-t48 + t70) * V_base(6) + (-rSges(4,3) + t71) * V_base(5) + t65;
t21 = V_base(4) * t48 + (-qJ(2) - t47) * V_base(5) + t63;
t20 = t53 * t35 + t72 * V_base(4) + V_base(6) * t74 + t62;
t19 = -t53 * t36 + (t70 - t73) * V_base(6) + (t71 - t72) * V_base(5) + t65;
t18 = (t36 + t73) * V_base(4) + (-qJ(2) - t35 - t74) * V_base(5) + t63;
t1 = m(1) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(2) * (t27 ^ 2 + t28 ^ 2 + t37 ^ 2) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t23 ^ 2) / 0.2e1 + m(5) * (t18 ^ 2 + t19 ^ 2 + t20 ^ 2) / 0.2e1 + Icges(5,3) * t53 ^ 2 / 0.2e1 + ((-Icges(5,5) * t51 + Icges(5,6) * t50) * t53 + (Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 - t60 * t44 / 0.2e1 + t46 * t75 - t50 * t32 / 0.2e1 + t34 * t76) * V_base(5)) * V_base(5) + ((-Icges(5,5) * t50 - Icges(5,6) * t51) * t53 + (Icges(1,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,3) / 0.2e1 + t43 * t75 + t60 * t77 + t31 * t76 + t50 * t78) * V_base(4) + (Icges(1,4) - Icges(2,6) + Icges(3,5) + (t44 / 0.2e1 + t77) * t61 + (t46 / 0.2e1 - t43 / 0.2e1) * t60 + (t32 / 0.2e1 + t78) * t51 + (t34 / 0.2e1 - t31 / 0.2e1) * t50) * V_base(5)) * V_base(4) + ((Icges(1,3) / 0.2e1 + Icges(2,1) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(4,3) / 0.2e1) * V_base(6) + (-Icges(4,5) * t61 + Icges(4,6) * t60 + Icges(3,4) - Icges(2,5) + Icges(1,6)) * V_base(5) + (-Icges(4,5) * t60 - Icges(4,6) * t61 + Icges(2,4) + Icges(1,5) + Icges(3,6)) * V_base(4)) * V_base(6);
T  = t1;
