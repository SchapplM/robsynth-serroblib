% Calculate kinetic energy for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:00
% EndTime: 2020-01-03 11:20:01
% DurationCPUTime: 0.80s
% Computational Cost: add. (1833->132), mult. (2746->189), div. (0->0), fcn. (2148->10), ass. (0->48)
t48 = V_base(6) * pkin(5) + V_base(2);
t49 = -V_base(5) * pkin(5) + V_base(3);
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t39 = t58 * t48 + t56 * t49;
t42 = t56 * V_base(5) - t58 * V_base(6);
t50 = V_base(4) + qJD(1);
t32 = pkin(1) * t50 - qJ(2) * t42 + t39;
t38 = t56 * t48 - t49 * t58;
t43 = t56 * V_base(6) + t58 * V_base(5);
t35 = qJ(2) * t43 + t38;
t53 = sin(pkin(7));
t61 = cos(pkin(7));
t27 = t53 * t32 + t61 * t35;
t21 = qJ(3) * t50 + t27;
t36 = t42 * t53 - t61 * t43;
t37 = t61 * t42 + t53 * t43;
t40 = -pkin(1) * t43 + qJD(2) + V_base(1);
t23 = pkin(2) * t36 - qJ(3) * t37 + t40;
t52 = sin(pkin(8));
t60 = cos(pkin(8));
t15 = t60 * t21 + t52 * t23;
t10 = qJ(4) * t36 + t15;
t26 = t61 * t32 - t53 * t35;
t20 = -t50 * pkin(2) + qJD(3) - t26;
t29 = t37 * t52 - t60 * t50;
t30 = t60 * t37 + t52 * t50;
t13 = t29 * pkin(3) - t30 * qJ(4) + t20;
t51 = sin(pkin(9));
t54 = cos(pkin(9));
t6 = t54 * t10 + t51 * t13;
t5 = -t10 * t51 + t54 * t13;
t14 = -t52 * t21 + t60 * t23;
t9 = -t36 * pkin(3) + qJD(4) - t14;
t59 = V_base(1) ^ 2;
t57 = cos(qJ(5));
t55 = sin(qJ(5));
t28 = qJD(5) + t29;
t25 = t30 * t54 + t36 * t51;
t24 = -t30 * t51 + t36 * t54;
t17 = t24 * t55 + t25 * t57;
t16 = t24 * t57 - t25 * t55;
t7 = -t24 * pkin(4) + t9;
t4 = pkin(6) * t24 + t6;
t3 = pkin(4) * t29 - pkin(6) * t25 + t5;
t2 = t3 * t55 + t4 * t57;
t1 = t3 * t57 - t4 * t55;
t8 = m(2) * (t38 ^ 2 + t39 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(2) ^ 2 + V_base(3) ^ 2 + t59) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(1) * mrSges(2,1) + t38 * mrSges(2,3) + Ifges(2,2) * t43 / 0.2e1) * t43 + (t40 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t20 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t28 / 0.2e1) * t28 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(1) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,4) * t43 + Ifges(2,1) * t42 / 0.2e1) * t42 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,2) * t24 / 0.2e1) * t24 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t28 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t28 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t40 * mrSges(3,1) + t14 * mrSges(4,1) - t15 * mrSges(4,2) - t27 * mrSges(3,3) - Ifges(3,4) * t37 + Ifges(4,5) * t30 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t36) * t36 + (t39 * mrSges(2,1) + t26 * mrSges(3,1) - t38 * mrSges(2,2) - t27 * mrSges(3,2) + Ifges(2,5) * t42 + Ifges(3,5) * t37 + Ifges(2,6) * t43 - Ifges(3,6) * t36 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * t50) * t50 + (t20 * mrSges(4,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) - t15 * mrSges(4,3) - Ifges(4,4) * t30 + Ifges(5,5) * t25 - Ifges(4,6) * t36 + Ifges(5,6) * t24 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29) * t29;
T = t8;
