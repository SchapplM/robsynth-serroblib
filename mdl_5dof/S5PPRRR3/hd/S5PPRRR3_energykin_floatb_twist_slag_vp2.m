% Calculate kinetic energy for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:17
% EndTime: 2019-12-05 15:16:18
% DurationCPUTime: 0.88s
% Computational Cost: add. (1607->132), mult. (2478->192), div. (0->0), fcn. (1912->10), ass. (0->49)
t54 = sin(pkin(8));
t62 = cos(pkin(8));
t44 = t54 * V_base(4) - t62 * V_base(5);
t45 = t54 * V_base(5) + t62 * V_base(4);
t52 = V_base(3) + qJD(1);
t33 = pkin(1) * t44 - qJ(2) * t45 + t52;
t49 = V_base(5) * qJ(1) + V_base(1);
t50 = -V_base(4) * qJ(1) + V_base(2);
t43 = t62 * t49 + t54 * t50;
t38 = V_base(6) * qJ(2) + t43;
t53 = sin(pkin(9));
t55 = cos(pkin(9));
t27 = t53 * t33 + t55 * t38;
t21 = pkin(5) * t44 + t27;
t42 = -t54 * t49 + t62 * t50;
t36 = -V_base(6) * pkin(1) + qJD(2) - t42;
t40 = -t45 * t53 + t55 * V_base(6);
t41 = t45 * t55 + t53 * V_base(6);
t23 = -t40 * pkin(2) - t41 * pkin(5) + t36;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t15 = t61 * t21 + t58 * t23;
t39 = qJD(3) - t40;
t10 = pkin(6) * t39 + t15;
t26 = t33 * t55 - t53 * t38;
t20 = -pkin(2) * t44 - t26;
t30 = -t58 * t41 + t44 * t61;
t31 = t41 * t61 + t58 * t44;
t13 = -pkin(3) * t30 - pkin(6) * t31 + t20;
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t6 = t60 * t10 + t57 * t13;
t5 = -t10 * t57 + t60 * t13;
t14 = -t58 * t21 + t23 * t61;
t29 = qJD(4) - t30;
t9 = -t39 * pkin(3) - t14;
t59 = cos(qJ(5));
t56 = sin(qJ(5));
t28 = qJD(5) + t29;
t25 = t31 * t60 + t39 * t57;
t24 = -t31 * t57 + t39 * t60;
t17 = t24 * t56 + t25 * t59;
t16 = t24 * t59 - t25 * t56;
t7 = -t24 * pkin(4) + t9;
t4 = pkin(7) * t24 + t6;
t3 = pkin(4) * t29 - pkin(7) * t25 + t5;
t2 = t3 * t56 + t4 * t59;
t1 = t3 * t59 - t4 * t56;
t8 = m(2) * (t42 ^ 2 + t43 ^ 2 + t52 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t52 * mrSges(2,2) - t42 * mrSges(2,3) + Ifges(2,1) * t45 / 0.2e1) * t45 + (t36 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t41 / 0.2e1) * t41 + (t14 * mrSges(4,1) - t15 * mrSges(4,2) + Ifges(4,3) * t39 / 0.2e1) * t39 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t29 / 0.2e1) * t29 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t36 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t41 + Ifges(3,2) * t40 / 0.2e1) * t40 + (t20 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,5) * t39 + Ifges(4,1) * t31 / 0.2e1) * t31 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t29 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t28 + Ifges(6,1) * t17 / 0.2e1) * t17 + (-t20 * mrSges(4,1) + t15 * mrSges(4,3) + Ifges(4,4) * t31 + Ifges(4,6) * t39 + Ifges(4,2) * t30 / 0.2e1) * t30 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t29 + Ifges(5,2) * t24 / 0.2e1) * t24 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t28 + Ifges(6,2) * t16 / 0.2e1) * t16 + (V_base(2) * mrSges(1,1) + t42 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t43 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t45 + Ifges(1,6) * V_base(5) + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t52 * mrSges(2,1) + t26 * mrSges(3,1) - t27 * mrSges(3,2) - t43 * mrSges(2,3) - Ifges(2,4) * t45 + Ifges(3,5) * t41 - Ifges(2,6) * V_base(6) + Ifges(3,6) * t40 + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t44) * t44;
T = t8;
