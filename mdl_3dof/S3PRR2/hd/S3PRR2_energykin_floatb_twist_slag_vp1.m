% Calculate kinetic energy for
% S3PRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3PRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_energykin_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR2_energykin_floatb_twist_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3PRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR2_energykin_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR2_energykin_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRR2_energykin_floatb_twist_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:44
% EndTime: 2018-11-14 10:12:44
% DurationCPUTime: 0.19s
% Computational Cost: add. (194->103), mult. (208->131), div. (0->0), fcn. (84->4), ass. (0->46)
t51 = qJ(2) + qJ(3);
t46 = cos(t51);
t42 = Icges(4,4) * t46;
t45 = sin(t51);
t67 = Icges(4,1) * t45 / 0.2e1 + t42 / 0.2e1;
t53 = cos(qJ(2));
t48 = Icges(3,4) * t53;
t52 = sin(qJ(2));
t66 = Icges(3,1) * t52 / 0.2e1 + t48 / 0.2e1;
t65 = t46 / 0.2e1;
t64 = t53 / 0.2e1;
t63 = pkin(2) * t52;
t62 = pkin(2) * t53;
t61 = rSges(4,3) + pkin(4);
t60 = Icges(3,4) * t52;
t59 = Icges(4,4) * t45;
t58 = V_base(6) * qJ(1) + V_base(2);
t57 = V_base(4) * pkin(1) + V_base(3);
t56 = qJD(1) + V_base(1);
t44 = V_base(6) - qJD(2);
t55 = V_base(4) * pkin(3) + t58;
t54 = -V_base(6) * pkin(1) + t56;
t43 = -qJD(3) + t44;
t41 = t53 * rSges(3,1) - t52 * rSges(3,2);
t40 = t52 * rSges(3,1) + t53 * rSges(3,2);
t39 = Icges(3,1) * t53 - t60;
t37 = -Icges(3,2) * t52 + t48;
t36 = Icges(3,2) * t53 + t60;
t33 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t32 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t31 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t30 = -V_base(6) * rSges(2,1) + V_base(5) * rSges(2,2) + t56;
t29 = t46 * rSges(4,1) - t45 * rSges(4,2);
t28 = t45 * rSges(4,1) + t46 * rSges(4,2);
t27 = Icges(4,1) * t46 - t59;
t25 = -Icges(4,2) * t45 + t42;
t24 = Icges(4,2) * t46 + t59;
t21 = V_base(4) * rSges(2,1) + V_base(3) + (-qJ(1) - rSges(2,3)) * V_base(5);
t20 = -V_base(4) * rSges(2,2) + V_base(6) * rSges(2,3) + t58;
t19 = V_base(4) * rSges(3,3) + t44 * t40 + t55;
t18 = -t44 * t41 + (-pkin(3) - rSges(3,3)) * V_base(5) + t54;
t17 = V_base(4) * t41 + (-qJ(1) - t40) * V_base(5) + t57;
t16 = t43 * t28 + t44 * t63 + t61 * V_base(4) + t55;
t15 = -t44 * t62 - t43 * t29 + (-pkin(3) - t61) * V_base(5) + t54;
t14 = (t29 + t62) * V_base(4) + (-qJ(1) - t28 - t63) * V_base(5) + t57;
t1 = m(1) * (t31 ^ 2 + t32 ^ 2 + t33 ^ 2) / 0.2e1 + m(2) * (t20 ^ 2 + t21 ^ 2 + t30 ^ 2) / 0.2e1 + m(3) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + Icges(3,3) * t44 ^ 2 / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t16 ^ 2) / 0.2e1 + Icges(4,3) * t43 ^ 2 / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,2) / 0.2e1) * V_base(6) ^ 2 + ((-Icges(4,5) * t46 + Icges(4,6) * t45) * t43 + (-Icges(3,5) * t53 + Icges(3,6) * t52) * t44 + (Icges(1,2) / 0.2e1 + Icges(2,1) / 0.2e1 - t52 * t37 / 0.2e1 + t39 * t64 - t45 * t25 / 0.2e1 + t27 * t65) * V_base(5) + (Icges(2,4) + Icges(1,6)) * V_base(6)) * V_base(5) + ((-Icges(3,5) * t52 - Icges(3,6) * t53) * t44 + (-Icges(4,5) * t45 - Icges(4,6) * t46) * t43 + (Icges(1,1) / 0.2e1 + Icges(2,3) / 0.2e1 + t36 * t64 + t52 * t66 + t24 * t65 + t45 * t67) * V_base(4) + (Icges(1,5) + Icges(2,6)) * V_base(6) + (Icges(1,4) + Icges(2,5) + (t37 / 0.2e1 + t66) * t53 + (t39 / 0.2e1 - t36 / 0.2e1) * t52 + (t25 / 0.2e1 + t67) * t46 + (t27 / 0.2e1 - t24 / 0.2e1) * t45) * V_base(5)) * V_base(4);
T  = t1;
