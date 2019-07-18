% Calculate kinetic energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:04
% EndTime: 2019-07-18 17:17:05
% DurationCPUTime: 0.96s
% Computational Cost: add. (1489->124), mult. (1970->186), div. (0->0), fcn. (1588->10), ass. (0->47)
t39 = V_base(5) * pkin(4) + V_base(1);
t40 = -V_base(4) * pkin(4) + V_base(2);
t48 = sin(qJ(1));
t53 = cos(qJ(1));
t31 = t39 * t48 - t53 * t40;
t55 = t31 ^ 2;
t33 = t39 * t53 + t40 * t48;
t47 = sin(qJ(2));
t52 = cos(qJ(2));
t26 = -t33 * t47 + t52 * V_base(3);
t36 = -t48 * V_base(4) + t53 * V_base(5);
t35 = qJD(2) - t36;
t23 = pkin(1) * t35 + t26;
t27 = t33 * t52 + t47 * V_base(3);
t46 = sin(qJ(3));
t51 = cos(qJ(3));
t14 = t46 * t23 + t51 * t27;
t37 = t48 * V_base(5) + t53 * V_base(4);
t43 = V_base(6) + qJD(1);
t29 = -t37 * t47 + t43 * t52;
t30 = t37 * t52 + t43 * t47;
t19 = t29 * t51 - t46 * t30;
t20 = t29 * t46 + t30 * t51;
t22 = -pkin(1) * t29 + t31;
t10 = -pkin(2) * t19 - pkin(5) * t20 + t22;
t34 = qJD(3) + t35;
t12 = pkin(5) * t34 + t14;
t45 = sin(qJ(4));
t50 = cos(qJ(4));
t4 = t50 * t10 - t12 * t45;
t13 = t23 * t51 - t46 * t27;
t18 = qJD(4) - t19;
t11 = -pkin(2) * t34 - t13;
t54 = V_base(3) ^ 2;
t49 = cos(qJ(5));
t44 = sin(qJ(5));
t17 = qJD(5) + t18;
t16 = t20 * t50 + t34 * t45;
t15 = -t20 * t45 + t34 * t50;
t9 = t15 * t44 + t16 * t49;
t8 = t15 * t49 - t16 * t44;
t6 = -pkin(3) * t15 + t11;
t5 = t10 * t45 + t12 * t50;
t3 = pkin(3) * t18 + t4;
t2 = t3 * t44 + t49 * t5;
t1 = t3 * t49 - t44 * t5;
t7 = m(2) * (t33 ^ 2 + t54 + t55) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t54) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t55) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t22 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t6 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,1) * t9 / 0.2e1) * t9 + (-t31 * mrSges(2,1) - t33 * mrSges(2,2) + Ifges(2,3) * t43 / 0.2e1) * t43 + (t26 * mrSges(3,1) - t27 * mrSges(3,2) + Ifges(3,3) * t35 / 0.2e1) * t35 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t34 / 0.2e1) * t34 + (t4 * mrSges(5,1) - t5 * mrSges(5,2) + Ifges(5,3) * t18 / 0.2e1) * t18 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t6 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t9 + Ifges(6,2) * t8 / 0.2e1) * t8 + (V_base(3) * mrSges(2,2) + t31 * mrSges(2,3) + Ifges(2,5) * t43 + Ifges(2,1) * t37 / 0.2e1) * t37 + (t31 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,5) * t35 + Ifges(3,1) * t30 / 0.2e1) * t30 + (t22 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t34 + Ifges(4,1) * t20 / 0.2e1) * t20 + (t11 * mrSges(5,2) - t4 * mrSges(5,3) + Ifges(5,5) * t18 + Ifges(5,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t33 * mrSges(2,3) + Ifges(2,4) * t37 + Ifges(2,6) * t43 + Ifges(2,2) * t36 / 0.2e1) * t36 + (-t31 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t30 + Ifges(3,6) * t35 + Ifges(3,2) * t29 / 0.2e1) * t29 + (-t22 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t20 + Ifges(4,6) * t34 + Ifges(4,2) * t19 / 0.2e1) * t19 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t9 + Ifges(6,6) * t8 + Ifges(6,3) * t17 / 0.2e1) * t17 + (-t11 * mrSges(5,1) + t5 * mrSges(5,3) + Ifges(5,4) * t16 + Ifges(5,6) * t18 + Ifges(5,2) * t15 / 0.2e1) * t15;
T  = t7;
