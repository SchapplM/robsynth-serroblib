% Calculate kinetic energy for
% S5RPPPR2
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
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:58:59
% EndTime: 2022-01-23 08:59:00
% DurationCPUTime: 0.80s
% Computational Cost: add. (1703->132), mult. (2358->189), div. (0->0), fcn. (1808->10), ass. (0->48)
t47 = V_base(5) * pkin(5) + V_base(1);
t48 = -V_base(4) * pkin(5) + V_base(2);
t56 = sin(qJ(1));
t61 = cos(qJ(1));
t38 = -t56 * t47 + t48 * t61;
t50 = V_base(6) + qJD(1);
t34 = -t50 * pkin(1) + qJD(2) - t38;
t42 = t56 * V_base(5) + t61 * V_base(4);
t53 = sin(pkin(7));
t60 = cos(pkin(7));
t36 = t42 * t53 - t50 * t60;
t37 = t42 * t60 + t53 * t50;
t20 = t36 * pkin(2) - t37 * qJ(3) + t34;
t41 = t56 * V_base(4) - t61 * V_base(5);
t31 = pkin(1) * t41 - qJ(2) * t42 + V_base(3);
t39 = t61 * t47 + t56 * t48;
t35 = qJ(2) * t50 + t39;
t27 = t53 * t31 + t60 * t35;
t22 = qJ(3) * t41 + t27;
t52 = sin(pkin(8));
t59 = cos(pkin(8));
t14 = t52 * t20 + t59 * t22;
t10 = qJ(4) * t36 + t14;
t26 = t31 * t60 - t53 * t35;
t21 = -t41 * pkin(2) + qJD(3) - t26;
t28 = t37 * t52 - t41 * t59;
t29 = t37 * t59 + t52 * t41;
t12 = t28 * pkin(3) - t29 * qJ(4) + t21;
t51 = sin(pkin(9));
t54 = cos(pkin(9));
t6 = t54 * t10 + t51 * t12;
t5 = -t10 * t51 + t12 * t54;
t24 = -t29 * t51 + t36 * t54;
t13 = t20 * t59 - t52 * t22;
t9 = -t36 * pkin(3) + qJD(4) - t13;
t58 = V_base(3) ^ 2;
t57 = cos(qJ(5));
t55 = sin(qJ(5));
t25 = t29 * t54 + t36 * t51;
t23 = qJD(5) - t24;
t16 = t25 * t57 + t28 * t55;
t15 = -t25 * t55 + t28 * t57;
t7 = -t24 * pkin(4) - t25 * pkin(6) + t9;
t4 = pkin(6) * t28 + t6;
t3 = -pkin(4) * t28 - t5;
t2 = t4 * t57 + t55 * t7;
t1 = -t4 * t55 + t57 * t7;
t8 = m(2) * (t38 ^ 2 + t39 ^ 2 + t58) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t58) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + Ifges(2,3) * t50 / 0.2e1) * t50 + (t34 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t21 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,1) * t29 / 0.2e1) * t29 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t25 / 0.2e1) * t25 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,5) * t50 + Ifges(2,1) * t42 / 0.2e1) * t42 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,2) * t24 / 0.2e1) * t24 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t23 + Ifges(6,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t23 + Ifges(6,2) * t15 / 0.2e1) * t15 + (t34 * mrSges(3,1) + t13 * mrSges(4,1) - t14 * mrSges(4,2) - t27 * mrSges(3,3) - Ifges(3,4) * t37 + Ifges(4,5) * t29 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t36) * t36 + (V_base(3) * mrSges(2,1) + t26 * mrSges(3,1) - t27 * mrSges(3,2) - t39 * mrSges(2,3) - Ifges(2,4) * t42 + Ifges(3,5) * t37 - Ifges(2,6) * t50 - Ifges(3,6) * t36 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t41) * t41 + (t21 * mrSges(4,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) - t14 * mrSges(4,3) - Ifges(4,4) * t29 + Ifges(5,5) * t25 - Ifges(4,6) * t36 + Ifges(5,6) * t24 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t28) * t28;
T = t8;
