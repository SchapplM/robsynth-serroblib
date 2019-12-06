% Calculate kinetic energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:19
% EndTime: 2019-12-05 16:02:20
% DurationCPUTime: 0.84s
% Computational Cost: add. (1741->136), mult. (2976->191), div. (0->0), fcn. (2376->10), ass. (0->53)
t65 = pkin(2) + pkin(7);
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t41 = t49 * V_base(5) + t51 * V_base(4);
t64 = pkin(6) * t41;
t40 = -t49 * V_base(4) + t51 * V_base(5);
t55 = sin(qJ(2));
t52 = cos(pkin(5));
t63 = cos(qJ(2));
t61 = t52 * t63;
t50 = sin(pkin(5));
t62 = t50 * t63;
t29 = -t40 * t61 + t41 * t55 - V_base(6) * t62;
t46 = V_base(5) * qJ(1) + V_base(1);
t47 = -V_base(4) * qJ(1) + V_base(2);
t37 = -t46 * t49 + t51 * t47;
t31 = V_base(6) * pkin(1) - t52 * t64 + t37;
t48 = V_base(3) + qJD(1);
t34 = -pkin(1) * t40 - t50 * t64 + t48;
t20 = -t31 * t50 + t52 * t34;
t59 = t40 * t52 + t50 * V_base(6);
t30 = t63 * t41 + t59 * t55;
t60 = -qJ(3) * t30 + t20;
t11 = t65 * t29 + t60;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t36 = -t40 * t50 + t52 * V_base(6) + qJD(2);
t38 = t51 * t46 + t49 * t47;
t27 = t59 * pkin(6) + t38;
t16 = -t55 * t27 + t31 * t61 + t34 * t62;
t58 = qJD(3) - t16;
t9 = t30 * pkin(3) - t65 * t36 + t58;
t6 = t57 * t11 + t54 * t9;
t17 = t63 * t27 + (t31 * t52 + t34 * t50) * t55;
t15 = -t36 * qJ(3) - t17;
t5 = -t54 * t11 + t57 * t9;
t22 = t29 * t57 - t54 * t36;
t12 = -pkin(3) * t29 - t15;
t56 = cos(qJ(5));
t53 = sin(qJ(5));
t28 = qJD(4) + t30;
t23 = t54 * t29 + t36 * t57;
t21 = qJD(5) - t22;
t19 = t23 * t56 + t28 * t53;
t18 = -t23 * t53 + t28 * t56;
t14 = -t36 * pkin(2) + t58;
t13 = pkin(2) * t29 + t60;
t7 = -pkin(4) * t22 - pkin(8) * t23 + t12;
t4 = pkin(8) * t28 + t6;
t3 = -t28 * pkin(4) - t5;
t2 = t4 * t56 + t53 * t7;
t1 = -t4 * t53 + t56 * t7;
t8 = m(2) * (t37 ^ 2 + t38 ^ 2 + t48 ^ 2) / 0.2e1 + m(3) * (t16 ^ 2 + t17 ^ 2 + t20 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t48 * mrSges(2,2) - t37 * mrSges(2,3) + Ifges(2,1) * t41 / 0.2e1) * t41 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t21 / 0.2e1) * t21 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t48 * mrSges(2,1) + t38 * mrSges(2,3) + Ifges(2,4) * t41 + Ifges(2,2) * t40 / 0.2e1) * t40 + (t12 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t23 / 0.2e1) * t23 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t21 + Ifges(6,1) * t19 / 0.2e1) * t19 + (-t12 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t23 + Ifges(5,6) * t28 + Ifges(5,2) * t22 / 0.2e1) * t22 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t19 + Ifges(6,6) * t21 + Ifges(6,2) * t18 / 0.2e1) * t18 + (t16 * mrSges(3,1) - t17 * mrSges(3,2) + t14 * mrSges(4,2) - t15 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t36) * t36 + (t14 * mrSges(4,1) + t20 * mrSges(3,2) - t16 * mrSges(3,3) - t13 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t30 + (-Ifges(4,4) + Ifges(3,5)) * t36) * t30 + (V_base(2) * mrSges(1,1) + t37 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t38 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t41 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t40 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t20 * mrSges(3,1) + t15 * mrSges(4,1) - t13 * mrSges(4,2) - t17 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t29 + (Ifges(4,5) - Ifges(3,6)) * t36 + (-Ifges(3,4) - Ifges(4,6)) * t30) * t29;
T = t8;
