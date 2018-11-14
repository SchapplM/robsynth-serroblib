% Calculate kinetic energy for
% S3RPP1
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
%   pkin=[a2,a3,d1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3RPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_energykin_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPP1_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:32
% EndTime: 2018-11-14 10:13:33
% DurationCPUTime: 0.41s
% Computational Cost: add. (189->90), mult. (270->113), div. (0->0), fcn. (146->2), ass. (0->41)
t88 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t87 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t86 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t62 = sin(qJ(1));
t85 = t88 * t62;
t63 = cos(qJ(1));
t84 = t88 * t63;
t83 = rSges(4,3) + qJ(3);
t82 = t87 * t63 - t85;
t81 = t87 * t62 + t84;
t80 = -t86 * t63 - t85;
t79 = t86 * t62 - t84;
t78 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t77 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t76 = pkin(2) + rSges(4,1);
t49 = t62 * pkin(1) - qJ(2) * t63;
t71 = V_base(4) * t49 + V_base(3);
t70 = V_base(5) * pkin(3) + V_base(1);
t67 = qJD(2) * t62 + t70;
t66 = t62 * rSges(4,2) + t83 * t63;
t65 = -rSges(4,2) * t63 + t83 * t62;
t53 = pkin(1) * t63 + t62 * qJ(2);
t57 = V_base(6) + qJD(1);
t64 = -qJD(2) * t63 + t57 * t53 + V_base(2);
t56 = rSges(2,1) * t63 - t62 * rSges(2,2);
t55 = -rSges(3,2) * t63 + t62 * rSges(3,3);
t52 = t62 * rSges(2,1) + rSges(2,2) * t63;
t51 = -t62 * rSges(3,2) - rSges(3,3) * t63;
t30 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t29 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t28 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t25 = V_base(5) * rSges(2,3) - t52 * t57 + t70;
t24 = t56 * t57 + V_base(2) + (-rSges(2,3) - pkin(3)) * V_base(4);
t23 = t52 * V_base(4) - t56 * V_base(5) + V_base(3);
t22 = V_base(5) * rSges(3,1) + (-t49 - t51) * t57 + t67;
t21 = t57 * t55 + (-rSges(3,1) - pkin(3)) * V_base(4) + t64;
t20 = t51 * V_base(4) + (-t53 - t55) * V_base(5) + t71;
t19 = qJD(3) * t63 + t76 * V_base(5) + (-t49 - t65) * t57 + t67;
t18 = qJD(3) * t62 + t66 * t57 + (-pkin(3) - t76) * V_base(4) + t64;
t17 = t65 * V_base(4) + (-t53 - t66) * V_base(5) + t71;
t1 = m(1) * (t28 ^ 2 + t29 ^ 2 + t30 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t23 ^ 2 + t24 ^ 2 + t25 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t22 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(4,1) / 0.2e1) * t57 + (-t77 * V_base(5) + t78 * V_base(4)) * t63 + (t77 * V_base(4) + t78 * V_base(5)) * t62) * t57 + ((t80 * t62 + t81 * t63) * V_base(5) + (t79 * t62 + t82 * t63) * V_base(4)) * V_base(4) / 0.2e1 + ((t81 * t62 - t80 * t63) * V_base(5) + (t82 * t62 - t79 * t63) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
