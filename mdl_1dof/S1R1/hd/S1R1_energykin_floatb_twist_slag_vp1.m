% Calculate kinetic energy for
% S1R1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
% m [2x1]
%   mass of all robot links (including the base)
% rSges [2x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [2x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S1R1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_energykin_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1R1_energykin_floatb_twist_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S1R1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_energykin_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1R1_energykin_floatb_twist_slag_vp1: m has to be [2x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [2,3]), ...
  'S1R1_energykin_floatb_twist_slag_vp1: rSges has to be [2x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [2 6]), ...
  'S1R1_energykin_floatb_twist_slag_vp1: Icges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:51
% EndTime: 2020-06-19 09:12:51
% DurationCPUTime: 0.16s
% Computational Cost: add. (63->44), mult. (94->64), div. (0->0), fcn. (40->2), ass. (0->20)
t22 = cos(qJ(1));
t26 = t22 / 0.2e1;
t25 = pkin(1) + rSges(2,3);
t21 = sin(qJ(1));
t24 = Icges(2,4) * t21;
t20 = Icges(2,4) * t22;
t19 = V_base(6) + qJD(1);
t18 = t22 * rSges(2,1) - t21 * rSges(2,2);
t17 = t21 * rSges(2,1) + t22 * rSges(2,2);
t16 = Icges(2,1) * t22 - t24;
t15 = Icges(2,1) * t21 + t20;
t14 = -Icges(2,2) * t21 + t20;
t13 = Icges(2,2) * t22 + t24;
t10 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t9 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t8 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t7 = -t19 * t17 + t25 * V_base(5) + V_base(1);
t6 = t19 * t18 - t25 * V_base(4) + V_base(2);
t5 = V_base(4) * t17 - V_base(5) * t18 + V_base(3);
t1 = m(1) * (t10 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + m(2) * (t5 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + Icges(2,3) * t19 ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (Icges(2,5) * t21 + Icges(2,6) * t22) * t19 + (Icges(1,2) / 0.2e1 + t13 * t26 + t21 * t15 / 0.2e1) * V_base(5)) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t22 - Icges(2,6) * t21) * t19 + (Icges(1,1) / 0.2e1 - t21 * t14 / 0.2e1 + t16 * t26) * V_base(4) + ((t14 + t15) * t22 + (-t13 + t16) * t21) * V_base(5) / 0.2e1) * V_base(4);
T = t1;
