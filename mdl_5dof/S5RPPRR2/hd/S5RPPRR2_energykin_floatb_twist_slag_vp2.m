% Calculate kinetic energy for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:30
% EndTime: 2019-12-05 17:39:30
% DurationCPUTime: 0.78s
% Computational Cost: add. (1289->129), mult. (1676->176), div. (0->0), fcn. (1176->8), ass. (0->46)
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t36 = t49 * V_base(5) + t52 * V_base(4);
t44 = V_base(6) + qJD(1);
t40 = V_base(5) * pkin(5) + V_base(1);
t41 = -V_base(4) * pkin(5) + V_base(2);
t31 = -t49 * t40 + t41 * t52;
t54 = qJD(2) - t31;
t56 = pkin(1) + qJ(3);
t22 = pkin(2) * t36 - t56 * t44 + t54;
t35 = t49 * V_base(4) - t52 * V_base(5);
t55 = -qJ(2) * t36 + V_base(3);
t24 = t56 * t35 + t55;
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t16 = t45 * t22 + t46 * t24;
t29 = t35 * t46 - t44 * t45;
t11 = pkin(6) * t29 + t16;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t15 = t46 * t22 - t24 * t45;
t30 = t35 * t45 + t44 * t46;
t9 = pkin(3) * t36 - pkin(6) * t30 + t15;
t6 = t51 * t11 + t48 * t9;
t32 = t52 * t40 + t49 * t41;
t28 = -t44 * qJ(2) - t32;
t5 = -t11 * t48 + t51 * t9;
t25 = -pkin(2) * t35 + qJD(3) - t28;
t17 = -pkin(3) * t29 + t25;
t34 = qJD(4) + t36;
t53 = V_base(3) ^ 2;
t50 = cos(qJ(5));
t47 = sin(qJ(5));
t33 = qJD(5) + t34;
t27 = -pkin(1) * t44 + t54;
t26 = pkin(1) * t35 + t55;
t19 = t29 * t48 + t30 * t51;
t18 = t29 * t51 - t30 * t48;
t14 = t18 * t47 + t19 * t50;
t13 = t18 * t50 - t19 * t47;
t12 = -pkin(4) * t18 + t17;
t4 = pkin(7) * t18 + t6;
t3 = pkin(4) * t34 - pkin(7) * t19 + t5;
t2 = t3 * t47 + t4 * t50;
t1 = t3 * t50 - t4 * t47;
t7 = m(2) * (t31 ^ 2 + t32 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t25 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t34 / 0.2e1) * t34 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t33 / 0.2e1) * t33 + (t25 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t25 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,2) * t29 / 0.2e1) * t29 + (t17 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t34 + Ifges(5,1) * t19 / 0.2e1) * t19 + (t12 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t33 + Ifges(6,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t17 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t19 + Ifges(5,6) * t34 + Ifges(5,2) * t18 / 0.2e1) * t18 + (-t12 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t33 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t31 * mrSges(2,1) - t32 * mrSges(2,2) + t27 * mrSges(3,2) - t28 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t44) * t44 + (V_base(3) * mrSges(2,1) + t28 * mrSges(3,1) - t26 * mrSges(3,2) - t32 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t35 + (Ifges(3,5) - Ifges(2,6)) * t44) * t35 + (t27 * mrSges(3,1) + t15 * mrSges(4,1) + V_base(3) * mrSges(2,2) - t16 * mrSges(4,2) - t31 * mrSges(2,3) - t26 * mrSges(3,3) + Ifges(4,5) * t30 + Ifges(4,6) * t29 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(2,1) / 0.2e1) * t36 + (-Ifges(3,4) + Ifges(2,5)) * t44 + (-Ifges(2,4) - Ifges(3,6)) * t35) * t36;
T = t7;
