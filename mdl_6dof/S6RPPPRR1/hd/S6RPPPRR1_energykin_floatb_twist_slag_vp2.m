% Calculate kinetic energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:47
% EndTime: 2019-03-09 01:29:48
% DurationCPUTime: 0.75s
% Computational Cost: add. (1475->151), mult. (2141->194), div. (0->0), fcn. (1576->8), ass. (0->49)
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t41 = -t55 * V_base(4) + t58 * V_base(5);
t42 = t55 * V_base(5) + t58 * V_base(4);
t52 = sin(pkin(9));
t62 = cos(pkin(9));
t35 = -t62 * t41 + t42 * t52;
t51 = V_base(6) + qJD(1);
t46 = V_base(5) * pkin(6) + V_base(1);
t47 = -V_base(4) * pkin(6) + V_base(2);
t37 = -t46 * t55 + t58 * t47;
t28 = pkin(1) * t51 - qJ(2) * t42 + t37;
t38 = t58 * t46 + t55 * t47;
t31 = qJ(2) * t41 + t38;
t22 = t52 * t28 + t62 * t31;
t17 = -t51 * qJ(3) - t22;
t61 = qJD(4) - t17;
t10 = -pkin(7) * t51 + (-pkin(3) - pkin(4)) * t35 + t61;
t32 = t35 * qJ(4);
t36 = t52 * t41 + t42 * t62;
t39 = -pkin(1) * t41 + qJD(2) + V_base(3);
t60 = pkin(2) * t35 + t39;
t12 = t32 + (pkin(7) - qJ(3)) * t36 + t60;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t6 = t54 * t10 + t57 * t12;
t21 = t62 * t28 - t52 * t31;
t16 = -t51 * pkin(2) + qJD(3) - t21;
t5 = t10 * t57 - t12 * t54;
t26 = t36 * t57 - t51 * t54;
t13 = t36 * pkin(3) - t51 * qJ(4) + t16;
t9 = -pkin(4) * t36 - t13;
t18 = -qJ(3) * t36 + t60;
t59 = V_base(3) ^ 2;
t56 = cos(qJ(6));
t53 = sin(qJ(6));
t34 = qJD(5) - t35;
t27 = t36 * t54 + t51 * t57;
t23 = qJD(6) - t26;
t20 = t27 * t56 + t34 * t53;
t19 = -t27 * t53 + t34 * t56;
t15 = t18 + t32;
t14 = -pkin(3) * t35 + t61;
t7 = -pkin(5) * t26 - pkin(8) * t27 + t9;
t4 = pkin(8) * t34 + t6;
t3 = -pkin(5) * t34 - t5;
t2 = t4 * t56 + t53 * t7;
t1 = -t4 * t53 + t56 * t7;
t8 = m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(2) * (t37 ^ 2 + t38 ^ 2 + t59) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t39 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t37 * mrSges(2,3) + Ifges(2,1) * t42 / 0.2e1) * t42 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t38 * mrSges(2,3) + Ifges(2,4) * t42 + Ifges(2,2) * t41 / 0.2e1) * t41 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t23 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t34 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t23 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t16 * mrSges(4,1) + t13 * mrSges(5,1) + t39 * mrSges(3,2) - t15 * mrSges(5,2) - t21 * mrSges(3,3) - t18 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t36) * t36 + (t39 * mrSges(3,1) + t17 * mrSges(4,1) - t14 * mrSges(5,1) - t18 * mrSges(4,2) - t22 * mrSges(3,3) + t15 * mrSges(5,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t35 + (-Ifges(3,4) - Ifges(4,6) + Ifges(5,6)) * t36) * t35 + (t37 * mrSges(2,1) + t21 * mrSges(3,1) - t38 * mrSges(2,2) - t22 * mrSges(3,2) + t16 * mrSges(4,2) + t14 * mrSges(5,2) - t17 * mrSges(4,3) - t13 * mrSges(5,3) + Ifges(2,5) * t42 + Ifges(2,6) * t41 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t51 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t36 + (Ifges(5,4) + Ifges(4,5) - Ifges(3,6)) * t35) * t51;
T  = t8;
