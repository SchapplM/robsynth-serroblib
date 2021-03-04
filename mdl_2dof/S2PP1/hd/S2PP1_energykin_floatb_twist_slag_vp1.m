% Calculate kinetic energy for
% S2PP1
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
%   pkin=[a2]';
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
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2PP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2PP1_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2PP1_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S2PP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2PP1_energykin_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2PP1_energykin_floatb_twist_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2PP1_energykin_floatb_twist_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2PP1_energykin_floatb_twist_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2021-03-03 18:41:17
% EndTime: 2021-03-03 18:41:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (64->51), mult. (86->52), div. (0->0), fcn. (0->0), ass. (0->14)
t14 = pkin(1) - rSges(3,2);
t13 = qJ(2) + rSges(3,3);
t12 = V_base(6) * qJ(1) + V_base(2);
t11 = qJD(1) + V_base(1);
t9 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t8 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t7 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t6 = V_base(5) * rSges(2,1) + V_base(6) * rSges(2,2) + t11;
t5 = -V_base(4) * rSges(2,1) + V_base(6) * rSges(2,3) + t12;
t4 = -V_base(4) * rSges(2,2) + V_base(3) + (-qJ(1) - rSges(2,3)) * V_base(5);
t3 = V_base(3) + (-qJ(1) - rSges(3,1)) * V_base(5) + t13 * V_base(4);
t2 = V_base(6) * rSges(3,1) - t14 * V_base(4) + qJD(2) + t12;
t1 = -t13 * V_base(6) + t14 * V_base(5) + t11;
t10 = m(1) * (t7 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(2) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(3) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,1) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) ^ 2 + ((Icges(1,2) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,3) / 0.2e1) * V_base(5) + (-Icges(2,4) + Icges(1,6) - Icges(3,6)) * V_base(6)) * V_base(5) + ((Icges(1,1) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1) * V_base(4) + (-Icges(3,4) + Icges(1,5) + Icges(2,5)) * V_base(6) + (Icges(1,4) + Icges(3,5) - Icges(2,6)) * V_base(5)) * V_base(4);
T = t10;
