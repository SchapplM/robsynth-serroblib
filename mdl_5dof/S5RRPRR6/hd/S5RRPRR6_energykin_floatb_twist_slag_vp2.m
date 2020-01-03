% Calculate kinetic energy for
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:25
% EndTime: 2020-01-03 12:05:26
% DurationCPUTime: 0.89s
% Computational Cost: add. (1961->132), mult. (2746->193), div. (0->0), fcn. (2148->10), ass. (0->50)
t49 = V_base(6) * pkin(5) + V_base(2);
t50 = -V_base(5) * pkin(5) + V_base(3);
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t41 = t61 * t49 + t58 * t50;
t44 = t58 * V_base(5) - t61 * V_base(6);
t52 = V_base(4) + qJD(1);
t33 = pkin(1) * t52 - pkin(6) * t44 + t41;
t40 = t58 * t49 - t50 * t61;
t45 = t58 * V_base(6) + t61 * V_base(5);
t37 = pkin(6) * t45 + t40;
t57 = sin(qJ(2));
t63 = cos(qJ(2));
t27 = t57 * t33 + t63 * t37;
t51 = qJD(2) + t52;
t21 = qJ(3) * t51 + t27;
t38 = t44 * t57 - t63 * t45;
t39 = t63 * t44 + t57 * t45;
t42 = -pkin(1) * t45 + V_base(1);
t25 = pkin(2) * t38 - qJ(3) * t39 + t42;
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t15 = t54 * t21 + t53 * t25;
t10 = pkin(7) * t38 + t15;
t26 = t63 * t33 - t57 * t37;
t20 = -t51 * pkin(2) + qJD(3) - t26;
t30 = -t53 * t39 + t51 * t54;
t31 = t39 * t54 + t51 * t53;
t13 = -t30 * pkin(3) - t31 * pkin(7) + t20;
t56 = sin(qJ(4));
t60 = cos(qJ(4));
t6 = t60 * t10 + t56 * t13;
t5 = -t10 * t56 + t60 * t13;
t14 = -t53 * t21 + t25 * t54;
t29 = qJD(4) - t30;
t9 = -pkin(3) * t38 - t14;
t62 = V_base(1) ^ 2;
t59 = cos(qJ(5));
t55 = sin(qJ(5));
t28 = qJD(5) + t29;
t24 = t31 * t60 + t38 * t56;
t23 = -t31 * t56 + t38 * t60;
t17 = t23 * t55 + t24 * t59;
t16 = t23 * t59 - t24 * t55;
t7 = -pkin(4) * t23 + t9;
t4 = pkin(8) * t23 + t6;
t3 = pkin(4) * t29 - pkin(8) * t24 + t5;
t2 = t3 * t55 + t4 * t59;
t1 = t3 * t59 - t4 * t55;
t8 = m(2) * (t40 ^ 2 + t41 ^ 2 + t62) / 0.2e1 + m(1) * (V_base(2) ^ 2 + V_base(3) ^ 2 + t62) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t42 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t41 * mrSges(2,1) - t40 * mrSges(2,2) + Ifges(2,3) * t52 / 0.2e1) * t52 + (t26 * mrSges(3,1) - t27 * mrSges(3,2) + Ifges(3,3) * t51 / 0.2e1) * t51 + (t20 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,1) * t31 / 0.2e1) * t31 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t29 / 0.2e1) * t29 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t28 / 0.2e1) * t28 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(1) * mrSges(2,1) + t40 * mrSges(2,3) + Ifges(2,6) * t52 + Ifges(2,2) * t45 / 0.2e1) * t45 + (t42 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,5) * t51 + Ifges(3,1) * t39 / 0.2e1) * t39 + (-t20 * mrSges(4,1) + t15 * mrSges(4,3) + Ifges(4,4) * t31 + Ifges(4,2) * t30 / 0.2e1) * t30 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t29 + Ifges(5,1) * t24 / 0.2e1) * t24 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t28 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (V_base(1) * mrSges(2,2) - t41 * mrSges(2,3) + Ifges(2,4) * t45 + Ifges(2,5) * t52 + Ifges(2,1) * t44 / 0.2e1) * t44 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,6) * t29 + Ifges(5,2) * t23 / 0.2e1) * t23 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t28 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t42 * mrSges(3,1) + t14 * mrSges(4,1) - t15 * mrSges(4,2) - t27 * mrSges(3,3) - Ifges(3,4) * t39 + Ifges(4,5) * t31 - Ifges(3,6) * t51 + Ifges(4,6) * t30 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t38) * t38;
T = t8;
