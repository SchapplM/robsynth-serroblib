% Calculate kinetic energy for
% S5RRRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:18
% EndTime: 2019-12-31 22:21:19
% DurationCPUTime: 1.01s
% Computational Cost: add. (2137->132), mult. (2798->195), div. (0->0), fcn. (2196->10), ass. (0->51)
t56 = sin(qJ(1));
t61 = cos(qJ(1));
t43 = -t56 * V_base(4) + t61 * V_base(5);
t44 = t56 * V_base(5) + t61 * V_base(4);
t32 = -pkin(1) * t43 - pkin(6) * t44 + V_base(3);
t48 = V_base(5) * pkin(5) + V_base(1);
t49 = -V_base(4) * pkin(5) + V_base(2);
t39 = t61 * t48 + t56 * t49;
t51 = V_base(6) + qJD(1);
t35 = pkin(6) * t51 + t39;
t55 = sin(qJ(2));
t60 = cos(qJ(2));
t25 = t60 * t32 - t35 * t55;
t37 = t44 * t60 + t51 * t55;
t42 = qJD(2) - t43;
t22 = pkin(2) * t42 - pkin(7) * t37 + t25;
t26 = t55 * t32 + t60 * t35;
t36 = -t44 * t55 + t51 * t60;
t24 = pkin(7) * t36 + t26;
t54 = sin(qJ(3));
t59 = cos(qJ(3));
t13 = t54 * t22 + t59 * t24;
t27 = t36 * t59 - t37 * t54;
t11 = pkin(8) * t27 + t13;
t53 = sin(qJ(4));
t58 = cos(qJ(4));
t12 = t59 * t22 - t24 * t54;
t28 = t36 * t54 + t37 * t59;
t41 = qJD(3) + t42;
t9 = pkin(3) * t41 - pkin(8) * t28 + t12;
t6 = t58 * t11 + t53 * t9;
t38 = -t56 * t48 + t49 * t61;
t5 = -t11 * t53 + t58 * t9;
t17 = t27 * t58 - t28 * t53;
t34 = -pkin(1) * t51 - t38;
t29 = -pkin(2) * t36 + t34;
t19 = -pkin(3) * t27 + t29;
t62 = V_base(3) ^ 2;
t57 = cos(qJ(5));
t52 = sin(qJ(5));
t40 = qJD(4) + t41;
t18 = t27 * t53 + t28 * t58;
t16 = qJD(5) - t17;
t15 = t18 * t57 + t40 * t52;
t14 = -t18 * t52 + t40 * t57;
t7 = -pkin(4) * t17 - pkin(9) * t18 + t19;
t4 = pkin(9) * t40 + t6;
t3 = -pkin(4) * t40 - t5;
t2 = t4 * t57 + t52 * t7;
t1 = -t4 * t52 + t57 * t7;
t8 = m(2) * (t38 ^ 2 + t39 ^ 2 + t62) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t62) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + Ifges(2,3) * t51 / 0.2e1) * t51 + (t25 * mrSges(3,1) - t26 * mrSges(3,2) + Ifges(3,3) * t42 / 0.2e1) * t42 + (t12 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t40 / 0.2e1) * t40 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t16 / 0.2e1) * t16 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,5) * t51 + Ifges(2,1) * t44 / 0.2e1) * t44 + (t34 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,5) * t42 + Ifges(3,1) * t37 / 0.2e1) * t37 + (t29 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t28 / 0.2e1) * t28 + (t19 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t40 + Ifges(5,1) * t18 / 0.2e1) * t18 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t16 + Ifges(6,1) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,6) * t51 + Ifges(2,2) * t43 / 0.2e1) * t43 + (-t34 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,6) * t42 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t29 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t28 + Ifges(4,6) * t41 + Ifges(4,2) * t27 / 0.2e1) * t27 + (-t19 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t18 + Ifges(5,6) * t40 + Ifges(5,2) * t17 / 0.2e1) * t17 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t16 + Ifges(6,2) * t14 / 0.2e1) * t14;
T = t8;
