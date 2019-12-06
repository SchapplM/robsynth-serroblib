% Calculate kinetic energy for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:42
% EndTime: 2019-12-05 15:44:43
% DurationCPUTime: 0.85s
% Computational Cost: add. (1795->132), mult. (2798->192), div. (0->0), fcn. (2196->10), ass. (0->49)
t52 = sin(pkin(8));
t54 = cos(pkin(8));
t42 = -t52 * V_base(4) + t54 * V_base(5);
t43 = t52 * V_base(5) + t54 * V_base(4);
t50 = V_base(3) + qJD(1);
t32 = -pkin(1) * t42 - pkin(5) * t43 + t50;
t47 = V_base(5) * qJ(1) + V_base(1);
t48 = -V_base(4) * qJ(1) + V_base(2);
t39 = t54 * t47 + t52 * t48;
t35 = V_base(6) * pkin(5) + t39;
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t25 = t60 * t32 - t35 * t57;
t37 = t43 * t60 + t57 * V_base(6);
t41 = qJD(2) - t42;
t21 = pkin(2) * t41 - qJ(3) * t37 + t25;
t26 = t57 * t32 + t60 * t35;
t36 = -t57 * t43 + t60 * V_base(6);
t24 = qJ(3) * t36 + t26;
t51 = sin(pkin(9));
t53 = cos(pkin(9));
t13 = t51 * t21 + t53 * t24;
t27 = t36 * t53 - t37 * t51;
t11 = pkin(6) * t27 + t13;
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t12 = t53 * t21 - t24 * t51;
t28 = t36 * t51 + t37 * t53;
t9 = pkin(3) * t41 - pkin(6) * t28 + t12;
t6 = t59 * t11 + t56 * t9;
t38 = -t52 * t47 + t48 * t54;
t5 = -t11 * t56 + t59 * t9;
t17 = t27 * t59 - t28 * t56;
t34 = -V_base(6) * pkin(1) - t38;
t29 = -pkin(2) * t36 + qJD(3) + t34;
t22 = -pkin(3) * t27 + t29;
t58 = cos(qJ(5));
t55 = sin(qJ(5));
t40 = qJD(4) + t41;
t18 = t27 * t56 + t28 * t59;
t16 = qJD(5) - t17;
t15 = t18 * t58 + t40 * t55;
t14 = -t18 * t55 + t40 * t58;
t7 = -pkin(4) * t17 - pkin(7) * t18 + t22;
t4 = pkin(7) * t40 + t6;
t3 = -pkin(4) * t40 - t5;
t2 = t4 * t58 + t55 * t7;
t1 = -t4 * t55 + t58 * t7;
t8 = m(2) * (t38 ^ 2 + t39 ^ 2 + t50 ^ 2) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t29 ^ 2) / 0.2e1 + m(5) * (t22 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t50 * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,1) * t43 / 0.2e1) * t43 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t40 / 0.2e1) * t40 + (t34 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t29 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,1) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t50 * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,4) * t43 + Ifges(2,2) * t42 / 0.2e1) * t42 + (-t34 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t29 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t28 + Ifges(4,2) * t27 / 0.2e1) * t27 + (t22 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t40 + Ifges(5,1) * t18 / 0.2e1) * t18 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t16 + Ifges(6,1) * t15 / 0.2e1) * t15 + (-t22 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t18 + Ifges(5,6) * t40 + Ifges(5,2) * t17 / 0.2e1) * t17 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t16 + Ifges(6,2) * t14 / 0.2e1) * t14 + (V_base(2) * mrSges(1,1) + t38 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t39 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t43 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t42 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t25 * mrSges(3,1) + t12 * mrSges(4,1) - t26 * mrSges(3,2) - t13 * mrSges(4,2) + Ifges(3,5) * t37 + Ifges(4,5) * t28 + Ifges(3,6) * t36 + Ifges(4,6) * t27 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t41) * t41;
T = t8;
