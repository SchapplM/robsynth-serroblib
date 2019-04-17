% Calculate kinetic energy for
% S3PPP1
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
%   pkin=[a2,a3,theta1]';
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
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3PPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPP1_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3PPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPP1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PPP1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PPP1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-17 09:48:03
% EndTime: 2019-04-17 09:48:03
% DurationCPUTime: 0.41s
% Computational Cost: add. (168->93), mult. (270->109), div. (0->0), fcn. (146->2), ass. (0->42)
t88 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t87 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t86 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t61 = sin(pkin(3));
t85 = t88 * t61;
t62 = cos(pkin(3));
t84 = t88 * t62;
t83 = t87 * t62 - t85;
t82 = t87 * t61 + t84;
t81 = -t86 * t62 - t85;
t80 = t86 * t61 - t84;
t79 = V_base(4) / 0.2e1;
t78 = V_base(5) / 0.2e1;
t77 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t76 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t75 = pkin(2) + rSges(4,1);
t53 = pkin(1) * t62 + t61 * qJ(2);
t70 = V_base(6) * t53 + V_base(2);
t69 = V_base(5) * qJ(1) + V_base(1);
t66 = qJD(1) + V_base(3);
t65 = qJD(2) * t61 + t69;
t64 = -rSges(4,2) * t62 + (rSges(4,3) + qJ(3)) * t61;
t49 = t61 * pkin(1) - qJ(2) * t62;
t63 = V_base(4) * t49 + t66;
t56 = rSges(2,1) * t62 - t61 * rSges(2,2);
t55 = -rSges(3,2) * t62 + t61 * rSges(3,3);
t52 = t61 * rSges(2,1) + rSges(2,2) * t62;
t51 = -t61 * rSges(3,2) - rSges(3,3) * t62;
t50 = t61 * rSges(4,2) + rSges(4,3) * t62;
t30 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t29 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t28 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t25 = V_base(5) * rSges(2,3) - t52 * V_base(6) + t69;
t24 = t56 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t23 = t52 * V_base(4) - t56 * V_base(5) + t66;
t22 = V_base(5) * rSges(3,1) + (-t49 - t51) * V_base(6) + t65;
t21 = -qJD(2) * t62 + V_base(6) * t55 + (-rSges(3,1) - qJ(1)) * V_base(4) + t70;
t20 = t51 * V_base(4) + (-t53 - t55) * V_base(5) + t63;
t19 = qJD(3) * t62 + t75 * V_base(5) + (-t49 - t64) * V_base(6) + t65;
t18 = qJD(3) * t61 + V_base(6) * t50 + (qJ(3) * V_base(6) - qJD(2)) * t62 + (-qJ(1) - t75) * V_base(4) + t70;
t17 = t64 * V_base(4) + (-qJ(3) * t62 - t50 - t53) * V_base(5) + t63;
t1 = m(1) * (t28 ^ 2 + t29 ^ 2 + t30 ^ 2) / 0.2e1 + m(2) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t22 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(4,1) / 0.2e1) * V_base(6) ^ 2 + ((t77 * t61 - t76 * t62 + Icges(1,6)) * V_base(6) + Icges(1,2) * t78) * V_base(5) + (Icges(1,4) * V_base(5) + (t76 * t61 + t77 * t62 + Icges(1,5)) * V_base(6) + Icges(1,1) * t79) * V_base(4) + ((t81 * t61 + t82 * t62) * V_base(5) + (t80 * t61 + t83 * t62) * V_base(4)) * t79 + ((t82 * t61 - t81 * t62) * V_base(5) + (t83 * t61 - t80 * t62) * V_base(4)) * t78;
T  = t1;
