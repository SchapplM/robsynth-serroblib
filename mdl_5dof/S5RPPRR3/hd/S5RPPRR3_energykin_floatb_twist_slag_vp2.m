% Calculate kinetic energy for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:23
% EndTime: 2020-01-03 11:27:24
% DurationCPUTime: 0.96s
% Computational Cost: add. (2001->132), mult. (2966->191), div. (0->0), fcn. (2344->10), ass. (0->49)
t48 = V_base(6) * pkin(5) + V_base(2);
t49 = -V_base(5) * pkin(5) + V_base(3);
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t40 = t60 * t48 + t57 * t49;
t43 = t57 * V_base(5) - t60 * V_base(6);
t50 = V_base(4) + qJD(1);
t31 = pkin(1) * t50 - qJ(2) * t43 + t40;
t39 = t57 * t48 - t49 * t60;
t44 = t57 * V_base(6) + t60 * V_base(5);
t35 = qJ(2) * t44 + t39;
t52 = sin(pkin(8));
t54 = cos(pkin(8));
t27 = t52 * t31 + t54 * t35;
t22 = qJ(3) * t50 + t27;
t37 = t43 * t52 - t54 * t44;
t38 = t43 * t54 + t44 * t52;
t41 = -pkin(1) * t44 + qJD(2) + V_base(1);
t25 = pkin(2) * t37 - qJ(3) * t38 + t41;
t51 = sin(pkin(9));
t53 = cos(pkin(9));
t16 = t53 * t22 + t51 * t25;
t28 = -t38 * t51 + t50 * t53;
t12 = pkin(6) * t28 + t16;
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t15 = -t22 * t51 + t53 * t25;
t29 = t38 * t53 + t50 * t51;
t9 = pkin(3) * t37 - pkin(6) * t29 + t15;
t6 = t59 * t12 + t56 * t9;
t5 = -t12 * t56 + t59 * t9;
t26 = t31 * t54 - t52 * t35;
t36 = qJD(4) + t37;
t21 = -pkin(2) * t50 + qJD(3) - t26;
t17 = -pkin(3) * t28 + t21;
t61 = V_base(1) ^ 2;
t58 = cos(qJ(5));
t55 = sin(qJ(5));
t34 = qJD(5) + t36;
t19 = t28 * t56 + t29 * t59;
t18 = t28 * t59 - t29 * t56;
t14 = t18 * t55 + t19 * t58;
t13 = t18 * t58 - t19 * t55;
t10 = -pkin(4) * t18 + t17;
t4 = pkin(7) * t18 + t6;
t3 = pkin(4) * t36 - pkin(7) * t19 + t5;
t2 = t3 * t55 + t4 * t58;
t1 = t3 * t58 - t4 * t55;
t7 = m(2) * (t39 ^ 2 + t40 ^ 2 + t61) / 0.2e1 + m(1) * (V_base(2) ^ 2 + V_base(3) ^ 2 + t61) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t41 ^ 2) / 0.2e1 + m(5) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t21 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(1) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,2) * t44 / 0.2e1) * t44 + (t41 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t38 / 0.2e1) * t38 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t36 / 0.2e1) * t36 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (t21 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,1) * t29 / 0.2e1) * t29 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(1) * mrSges(2,2) - t40 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,1) * t43 / 0.2e1) * t43 + (-t21 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t29 + Ifges(4,2) * t28 / 0.2e1) * t28 + (t17 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t36 + Ifges(5,1) * t19 / 0.2e1) * t19 + (t10 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t17 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t19 + Ifges(5,6) * t36 + Ifges(5,2) * t18 / 0.2e1) * t18 + (-t10 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t34 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t40 * mrSges(2,1) + t26 * mrSges(3,1) - t39 * mrSges(2,2) - t27 * mrSges(3,2) + Ifges(2,5) * t43 + Ifges(3,5) * t38 + Ifges(2,6) * t44 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t50) * t50 + (t41 * mrSges(3,1) + t15 * mrSges(4,1) - t16 * mrSges(4,2) - t27 * mrSges(3,3) - Ifges(3,4) * t38 + Ifges(4,5) * t29 - Ifges(3,6) * t50 + Ifges(4,6) * t28 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t37) * t37;
T = t7;
