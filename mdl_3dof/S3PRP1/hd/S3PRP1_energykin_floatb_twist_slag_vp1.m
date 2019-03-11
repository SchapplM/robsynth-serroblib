% Calculate kinetic energy for
% S3PRP1
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
%   pkin=[a2,a3,d2]';
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
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3PRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3PRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRP1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRP1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:02:55
% EndTime: 2019-03-08 18:02:55
% DurationCPUTime: 0.21s
% Computational Cost: add. (152->89), mult. (212->101), div. (0->0), fcn. (90->2), ass. (0->37)
t74 = Icges(4,5) / 0.2e1 - Icges(3,4) / 0.2e1;
t73 = rSges(4,1) + pkin(2);
t72 = Icges(4,3) / 0.2e1 + Icges(3,2) / 0.2e1;
t71 = Icges(4,1) / 0.2e1 + Icges(3,1) / 0.2e1;
t70 = rSges(4,3) + qJ(3);
t65 = Icges(4,4) + Icges(3,5);
t64 = Icges(3,6) - Icges(4,6);
t48 = sin(qJ(2));
t63 = t74 * t48;
t49 = cos(qJ(2));
t62 = t74 * t49;
t61 = t73 * t48 - t70 * t49;
t60 = t70 * t48 + t73 * t49;
t57 = V_base(4) * qJ(1) + V_base(3);
t56 = qJD(1) + V_base(2);
t55 = t72 * t49 - t63;
t54 = t72 * t48 + t62;
t53 = t71 * t48 - t62;
t52 = t71 * t49 + t63;
t51 = V_base(6) * pkin(1) + t56;
t50 = V_base(5) * pkin(3) - V_base(6) * qJ(1) + V_base(1);
t42 = V_base(6) + qJD(2);
t41 = t49 * rSges(3,1) - t48 * rSges(3,2);
t38 = t48 * rSges(3,1) + t49 * rSges(3,2);
t23 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t22 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t21 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t20 = V_base(6) * rSges(2,1) + V_base(4) * rSges(2,2) + t56;
t19 = -V_base(5) * rSges(2,1) + V_base(4) * rSges(2,3) + t57;
t18 = -V_base(5) * rSges(2,2) + V_base(1) + (-qJ(1) - rSges(2,3)) * V_base(6);
t17 = V_base(5) * rSges(3,3) - t42 * t38 + t50;
t16 = t42 * t41 + (-pkin(3) - rSges(3,3)) * V_base(4) + t51;
t15 = V_base(4) * t38 + (-pkin(1) - t41) * V_base(5) + t57;
t14 = V_base(5) * rSges(4,2) + qJD(3) * t48 - t61 * t42 + t50;
t13 = -qJD(3) * t49 + (-pkin(3) - rSges(4,2)) * V_base(4) + t60 * t42 + t51;
t12 = t61 * V_base(4) + (-pkin(1) - t60) * V_base(5) + t57;
t1 = m(1) * (t21 ^ 2 + t22 ^ 2 + t23 ^ 2) / 0.2e1 + m(2) * (t18 ^ 2 + t19 ^ 2 + t20 ^ 2) / 0.2e1 + m(3) * (t15 ^ 2 + t16 ^ 2 + t17 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t14 ^ 2) / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,2) / 0.2e1) * V_base(6) ^ 2 + (Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t42 ^ 2 + ((Icges(1,2) / 0.2e1 + Icges(2,3) / 0.2e1 + t55 * t49 + t53 * t48) * V_base(5) + (Icges(1,6) - Icges(2,6)) * V_base(6) + (t65 * t48 + t64 * t49) * t42) * V_base(5) + ((Icges(1,1) / 0.2e1 + Icges(2,1) / 0.2e1 + t52 * t49 + t54 * t48) * V_base(4) + (Icges(1,5) - Icges(2,4)) * V_base(6) + (Icges(1,4) + Icges(2,5) + (t53 - t54) * t49 + (t52 - t55) * t48) * V_base(5) + (-t64 * t48 + t65 * t49) * t42) * V_base(4);
T  = t1;
