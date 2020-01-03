% Calculate kinetic energy for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:07
% EndTime: 2019-12-31 19:03:08
% DurationCPUTime: 0.90s
% Computational Cost: add. (1899->132), mult. (2746->193), div. (0->0), fcn. (2148->10), ass. (0->50)
t49 = pkin(5) * V_base(5) + V_base(1);
t50 = -pkin(5) * V_base(4) + V_base(2);
t57 = sin(qJ(1));
t61 = cos(qJ(1));
t41 = -t49 * t57 + t61 * t50;
t45 = t57 * V_base(5) + t61 * V_base(4);
t51 = V_base(6) + qJD(1);
t33 = pkin(1) * t51 - qJ(2) * t45 + t41;
t42 = t61 * t49 + t57 * t50;
t44 = -t57 * V_base(4) + t61 * V_base(5);
t37 = qJ(2) * t44 + t42;
t52 = sin(pkin(9));
t53 = cos(pkin(9));
t27 = t52 * t33 + t53 * t37;
t21 = pkin(6) * t51 + t27;
t39 = t44 * t53 - t45 * t52;
t40 = t44 * t52 + t45 * t53;
t43 = -pkin(1) * t44 + qJD(2) + V_base(3);
t23 = -pkin(2) * t39 - pkin(6) * t40 + t43;
t56 = sin(qJ(3));
t60 = cos(qJ(3));
t14 = t60 * t21 + t56 * t23;
t38 = qJD(3) - t39;
t10 = pkin(7) * t38 + t14;
t26 = t33 * t53 - t52 * t37;
t20 = -pkin(2) * t51 - t26;
t31 = -t56 * t40 + t51 * t60;
t32 = t40 * t60 + t51 * t56;
t15 = -pkin(3) * t31 - pkin(7) * t32 + t20;
t55 = sin(qJ(4));
t59 = cos(qJ(4));
t6 = t59 * t10 + t55 * t15;
t5 = -t10 * t55 + t59 * t15;
t13 = -t56 * t21 + t23 * t60;
t29 = qJD(4) - t31;
t9 = -pkin(3) * t38 - t13;
t62 = V_base(3) ^ 2;
t58 = cos(qJ(5));
t54 = sin(qJ(5));
t28 = qJD(5) + t29;
t25 = t32 * t59 + t38 * t55;
t24 = -t32 * t55 + t38 * t59;
t17 = t24 * t54 + t25 * t58;
t16 = t24 * t58 - t25 * t54;
t7 = -pkin(4) * t24 + t9;
t4 = pkin(8) * t24 + t6;
t3 = pkin(4) * t29 - pkin(8) * t25 + t5;
t2 = t3 * t54 + t4 * t58;
t1 = t3 * t58 - t4 * t54;
t8 = m(2) * (t41 ^ 2 + t42 ^ 2 + t62) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t62) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t43 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t41 * mrSges(2,3) + Ifges(2,1) * t45 / 0.2e1) * t45 + (t43 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t40 / 0.2e1) * t40 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t38 / 0.2e1) * t38 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t29 / 0.2e1) * t29 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t28 / 0.2e1) * t28 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t42 * mrSges(2,3) + Ifges(2,4) * t45 + Ifges(2,2) * t44 / 0.2e1) * t44 + (-t43 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t40 + Ifges(3,2) * t39 / 0.2e1) * t39 + (t20 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t38 + Ifges(4,1) * t32 / 0.2e1) * t32 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t29 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t28 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t20 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t38 + Ifges(4,2) * t31 / 0.2e1) * t31 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t29 + Ifges(5,2) * t24 / 0.2e1) * t24 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t28 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t41 * mrSges(2,1) + t26 * mrSges(3,1) - t42 * mrSges(2,2) - t27 * mrSges(3,2) + Ifges(2,5) * t45 + Ifges(3,5) * t40 + Ifges(2,6) * t44 + Ifges(3,6) * t39 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * t51) * t51;
T = t8;
