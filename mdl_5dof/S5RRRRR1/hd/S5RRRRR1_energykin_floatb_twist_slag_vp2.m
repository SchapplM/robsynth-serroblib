% Calculate kinetic energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:51
% EndTime: 2019-03-08 18:36:52
% DurationCPUTime: 0.92s
% Computational Cost: add. (1577->126), mult. (2130->189), div. (0->0), fcn. (1700->10), ass. (0->48)
t38 = V_base(5) * pkin(5) + V_base(1);
t39 = -V_base(4) * pkin(5) + V_base(2);
t47 = sin(qJ(1));
t52 = cos(qJ(1));
t30 = t38 * t52 + t39 * t47;
t35 = -t47 * V_base(4) + t52 * V_base(5);
t32 = -pkin(1) * t35 + V_base(3);
t46 = sin(qJ(2));
t51 = cos(qJ(2));
t24 = -t30 * t46 - t32 * t51;
t34 = qJD(2) + t35;
t20 = pkin(2) * t34 + t24;
t25 = t30 * t51 - t32 * t46;
t45 = sin(qJ(3));
t50 = cos(qJ(3));
t16 = t50 * t20 - t25 * t45;
t33 = qJD(3) + t34;
t11 = pkin(3) * t33 + t16;
t17 = t20 * t45 + t25 * t50;
t44 = sin(qJ(4));
t49 = cos(qJ(4));
t7 = t44 * t11 + t49 * t17;
t29 = -t38 * t47 + t52 * t39;
t42 = V_base(6) + qJD(1);
t26 = t42 * pkin(1) + t29;
t6 = t11 * t49 - t17 * t44;
t36 = t47 * V_base(5) + t52 * V_base(4);
t27 = -t36 * t46 - t42 * t51;
t28 = t36 * t51 - t42 * t46;
t21 = t27 * t50 - t28 * t45;
t22 = t27 * t45 + t28 * t50;
t13 = t21 * t49 - t22 * t44;
t23 = -pkin(2) * t27 + t26;
t18 = -pkin(3) * t21 + t23;
t53 = V_base(3) ^ 2;
t48 = cos(qJ(5));
t43 = sin(qJ(5));
t31 = qJD(4) + t33;
t14 = t21 * t44 + t22 * t49;
t12 = qJD(5) - t13;
t9 = t14 * t48 + t31 * t43;
t8 = -t14 * t43 + t31 * t48;
t5 = pkin(6) * t31 + t7;
t4 = -pkin(4) * t31 - t6;
t3 = -pkin(4) * t13 - pkin(6) * t14 + t18;
t2 = t3 * t43 + t48 * t5;
t1 = t3 * t48 - t43 * t5;
t10 = m(2) * (t29 ^ 2 + t30 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t18 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t23 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,1) * t9 / 0.2e1) * t9 + (t29 * mrSges(2,1) - t30 * mrSges(2,2) + Ifges(2,3) * t42 / 0.2e1) * t42 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t34 / 0.2e1) * t34 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t33 / 0.2e1) * t33 + (t6 * mrSges(5,1) - t7 * mrSges(5,2) + Ifges(5,3) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t9 + Ifges(6,2) * t8 / 0.2e1) * t8 + (V_base(3) * mrSges(2,2) - t29 * mrSges(2,3) + Ifges(2,5) * t42 + Ifges(2,1) * t36 / 0.2e1) * t36 + (t26 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,5) * t34 + Ifges(3,1) * t28 / 0.2e1) * t28 + (t23 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t33 + Ifges(4,1) * t22 / 0.2e1) * t22 + (t18 * mrSges(5,2) - t6 * mrSges(5,3) + Ifges(5,5) * t31 + Ifges(5,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t30 * mrSges(2,3) + Ifges(2,4) * t36 + Ifges(2,6) * t42 + Ifges(2,2) * t35 / 0.2e1) * t35 + (-t26 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t28 + Ifges(3,6) * t34 + Ifges(3,2) * t27 / 0.2e1) * t27 + (-t23 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t22 + Ifges(4,6) * t33 + Ifges(4,2) * t21 / 0.2e1) * t21 + (-t18 * mrSges(5,1) + t7 * mrSges(5,3) + Ifges(5,4) * t14 + Ifges(5,6) * t31 + Ifges(5,2) * t13 / 0.2e1) * t13 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t9 + Ifges(6,6) * t8 + Ifges(6,3) * t12 / 0.2e1) * t12;
T  = t10;
