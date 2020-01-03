% Calculate kinetic energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:39
% EndTime: 2019-12-31 19:05:39
% DurationCPUTime: 0.76s
% Computational Cost: add. (1209->129), mult. (1566->178), div. (0->0), fcn. (1092->8), ass. (0->45)
t52 = sin(qJ(1));
t58 = cos(qJ(1));
t38 = t52 * V_base(4) - t58 * V_base(5);
t39 = t52 * V_base(5) + t58 * V_base(4);
t30 = t38 * pkin(1) - t39 * qJ(2) + V_base(3);
t22 = -pkin(2) * t38 - t30;
t51 = sin(qJ(3));
t55 = cos(qJ(3));
t28 = t38 * t55 - t51 * t39;
t29 = t38 * t51 + t39 * t55;
t10 = -pkin(3) * t28 - pkin(7) * t29 + t22;
t48 = V_base(6) + qJD(1);
t43 = V_base(5) * pkin(5) + V_base(1);
t44 = -V_base(4) * pkin(5) + V_base(2);
t33 = -t52 * t43 + t44 * t58;
t57 = qJD(2) - t33;
t19 = -t39 * pkin(6) + (-pkin(1) - pkin(2)) * t48 + t57;
t34 = t58 * t43 + t52 * t44;
t32 = t48 * qJ(2) + t34;
t25 = pkin(6) * t38 + t32;
t17 = t51 * t19 + t55 * t25;
t46 = qJD(3) - t48;
t13 = pkin(7) * t46 + t17;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t6 = t50 * t10 + t54 * t13;
t5 = t54 * t10 - t13 * t50;
t16 = t19 * t55 - t51 * t25;
t27 = qJD(4) - t28;
t12 = -pkin(3) * t46 - t16;
t56 = V_base(3) ^ 2;
t53 = cos(qJ(5));
t49 = sin(qJ(5));
t31 = -t48 * pkin(1) + t57;
t26 = qJD(5) + t27;
t21 = t29 * t54 + t46 * t50;
t20 = -t29 * t50 + t46 * t54;
t15 = t20 * t49 + t21 * t53;
t14 = t20 * t53 - t21 * t49;
t7 = -pkin(4) * t20 + t12;
t4 = pkin(8) * t20 + t6;
t3 = pkin(4) * t27 - pkin(8) * t21 + t5;
t2 = t3 * t49 + t4 * t53;
t1 = t3 * t53 - t4 * t49;
t8 = m(2) * (t33 ^ 2 + t34 ^ 2 + t56) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t56) / 0.2e1 + m(3) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t22 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t46 / 0.2e1) * t46 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t27 / 0.2e1) * t27 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t22 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t46 + Ifges(4,1) * t29 / 0.2e1) * t29 + (t12 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t27 + Ifges(5,1) * t21 / 0.2e1) * t21 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t26 + Ifges(6,1) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t22 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t29 + Ifges(4,6) * t46 + Ifges(4,2) * t28 / 0.2e1) * t28 + (-t12 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,6) * t27 + Ifges(5,2) * t20 / 0.2e1) * t20 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t26 + Ifges(6,2) * t14 / 0.2e1) * t14 + (t33 * mrSges(2,1) - t31 * mrSges(3,1) - t34 * mrSges(2,2) + t32 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * t48) * t48 + (V_base(3) * mrSges(2,2) + t31 * mrSges(3,2) - t33 * mrSges(2,3) - t30 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t39 + (Ifges(3,4) + Ifges(2,5)) * t48) * t39 + (V_base(3) * mrSges(2,1) + t30 * mrSges(3,1) - t32 * mrSges(3,2) - t34 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t38 + (-Ifges(2,6) + Ifges(3,6)) * t48 + (-Ifges(2,4) + Ifges(3,5)) * t39) * t38;
T = t8;
