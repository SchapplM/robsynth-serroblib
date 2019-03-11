% Calculate kinetic energy for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
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
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3PPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PPR1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:01:49
% EndTime: 2019-03-08 18:01:49
% DurationCPUTime: 0.14s
% Computational Cost: add. (120->88), mult. (162->98), div. (0->0), fcn. (40->2), ass. (0->30)
t34 = cos(qJ(3));
t31 = Icges(4,4) * t34;
t33 = sin(qJ(3));
t42 = Icges(4,1) * t33 / 0.2e1 + t31 / 0.2e1;
t41 = t34 / 0.2e1;
t40 = rSges(4,3) + pkin(3);
t39 = Icges(4,4) * t33;
t38 = V_base(4) * qJ(1) + V_base(3);
t37 = qJD(1) + V_base(2);
t36 = qJD(2) + V_base(1);
t35 = V_base(4) * pkin(1) + V_base(6) * qJ(2) + t37;
t28 = V_base(6) - qJD(3);
t27 = t34 * rSges(4,1) - t33 * rSges(4,2);
t26 = t33 * rSges(4,1) + t34 * rSges(4,2);
t25 = Icges(4,1) * t34 - t39;
t23 = -Icges(4,2) * t33 + t31;
t22 = Icges(4,2) * t34 + t39;
t19 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t18 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t17 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t16 = V_base(4) * rSges(2,1) - V_base(6) * rSges(2,2) + t37;
t15 = -V_base(5) * rSges(2,1) + V_base(1) + (-qJ(1) - rSges(2,3)) * V_base(6);
t14 = V_base(5) * rSges(2,2) + V_base(4) * rSges(2,3) + t38;
t13 = V_base(4) * rSges(3,1) + (-qJ(2) - rSges(3,3)) * V_base(5) + t38;
t12 = (-qJ(1) - rSges(3,1)) * V_base(6) + (-pkin(1) + rSges(3,2)) * V_base(5) + t36;
t11 = -V_base(4) * rSges(3,2) + V_base(6) * rSges(3,3) + t35;
t10 = t28 * t26 + t40 * V_base(4) + t35;
t9 = -t28 * t27 + (-pkin(2) - qJ(1)) * V_base(6) + (-pkin(1) - t40) * V_base(5) + t36;
t8 = (-qJ(2) - t26) * V_base(5) + (pkin(2) + t27) * V_base(4) + t38;
t1 = m(1) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + m(2) * (t14 ^ 2 + t15 ^ 2 + t16 ^ 2) / 0.2e1 + m(3) * (t11 ^ 2 + t12 ^ 2 + t13 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + Icges(4,3) * t28 ^ 2 / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,1) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) ^ 2 + ((-Icges(4,5) * t34 + Icges(4,6) * t33) * t28 + (Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 - t33 * t23 / 0.2e1 + t25 * t41) * V_base(5) + (Icges(3,4) - Icges(2,5) + Icges(1,6)) * V_base(6)) * V_base(5) + ((-Icges(4,5) * t33 - Icges(4,6) * t34) * t28 + (Icges(1,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,3) / 0.2e1 + t22 * t41 + t33 * t42) * V_base(4) + (Icges(1,5) + Icges(2,4) + Icges(3,6)) * V_base(6) + (Icges(1,4) - Icges(2,6) + Icges(3,5) + (t23 / 0.2e1 + t42) * t34 + (t25 / 0.2e1 - t22 / 0.2e1) * t33) * V_base(5)) * V_base(4);
T  = t1;
