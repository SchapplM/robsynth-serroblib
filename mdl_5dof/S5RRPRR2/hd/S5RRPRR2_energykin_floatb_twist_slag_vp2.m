% Calculate kinetic energy for
% S5RRPRR2
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:26:59
% EndTime: 2019-12-05 18:27:00
% DurationCPUTime: 0.97s
% Computational Cost: add. (2333->132), mult. (3134->193), div. (0->0), fcn. (2484->10), ass. (0->50)
t57 = sin(qJ(1));
t61 = cos(qJ(1));
t43 = -t57 * V_base(4) + t61 * V_base(5);
t44 = t57 * V_base(5) + t61 * V_base(4);
t32 = -pkin(1) * t43 - pkin(6) * t44 + V_base(3);
t48 = pkin(5) * V_base(5) + V_base(1);
t49 = -pkin(5) * V_base(4) + V_base(2);
t39 = t61 * t48 + t57 * t49;
t51 = V_base(6) + qJD(1);
t35 = pkin(6) * t51 + t39;
t56 = sin(qJ(2));
t60 = cos(qJ(2));
t25 = t60 * t32 - t35 * t56;
t37 = t44 * t60 + t51 * t56;
t42 = qJD(2) - t43;
t22 = pkin(2) * t42 - qJ(3) * t37 + t25;
t26 = t56 * t32 + t60 * t35;
t36 = -t44 * t56 + t51 * t60;
t24 = qJ(3) * t36 + t26;
t52 = sin(pkin(9));
t53 = cos(pkin(9));
t16 = t52 * t22 + t53 * t24;
t27 = t36 * t53 - t37 * t52;
t11 = pkin(7) * t27 + t16;
t55 = sin(qJ(4));
t59 = cos(qJ(4));
t15 = t53 * t22 - t24 * t52;
t28 = t36 * t52 + t37 * t53;
t9 = pkin(3) * t42 - pkin(7) * t28 + t15;
t6 = t59 * t11 + t55 * t9;
t5 = -t11 * t55 + t59 * t9;
t38 = -t57 * t48 + t49 * t61;
t34 = -pkin(1) * t51 - t38;
t41 = qJD(4) + t42;
t29 = -pkin(2) * t36 + qJD(3) + t34;
t19 = -pkin(3) * t27 + t29;
t62 = V_base(3) ^ 2;
t58 = cos(qJ(5));
t54 = sin(qJ(5));
t40 = qJD(5) + t41;
t18 = t27 * t55 + t28 * t59;
t17 = t27 * t59 - t28 * t55;
t14 = -pkin(4) * t17 + t19;
t13 = t17 * t54 + t18 * t58;
t12 = t17 * t58 - t18 * t54;
t4 = pkin(8) * t17 + t6;
t3 = pkin(4) * t41 - pkin(8) * t18 + t5;
t2 = t3 * t54 + t4 * t58;
t1 = t3 * t58 - t4 * t54;
t7 = m(2) * (t38 ^ 2 + t39 ^ 2 + t62) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t62) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t29 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + Ifges(2,3) * t51 / 0.2e1) * t51 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t41 / 0.2e1) * t41 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t40 / 0.2e1) * t40 + (t34 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t29 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,1) * t28 / 0.2e1) * t28 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,5) * t51 + Ifges(2,1) * t44 / 0.2e1) * t44 + (-t34 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t29 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t28 + Ifges(4,2) * t27 / 0.2e1) * t27 + (t19 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t41 + Ifges(5,1) * t18 / 0.2e1) * t18 + (t14 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t40 + Ifges(6,1) * t13 / 0.2e1) * t13 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,6) * t51 + Ifges(2,2) * t43 / 0.2e1) * t43 + (-t19 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t18 + Ifges(5,6) * t41 + Ifges(5,2) * t17 / 0.2e1) * t17 + (-t14 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t13 + Ifges(6,6) * t40 + Ifges(6,2) * t12 / 0.2e1) * t12 + (t25 * mrSges(3,1) + t15 * mrSges(4,1) - t26 * mrSges(3,2) - t16 * mrSges(4,2) + Ifges(3,5) * t37 + Ifges(4,5) * t28 + Ifges(3,6) * t36 + Ifges(4,6) * t27 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t42) * t42;
T = t7;
