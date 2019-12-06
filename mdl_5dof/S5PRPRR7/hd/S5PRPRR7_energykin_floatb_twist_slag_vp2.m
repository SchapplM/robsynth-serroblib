% Calculate kinetic energy for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:38
% EndTime: 2019-12-05 15:59:39
% DurationCPUTime: 0.73s
% Computational Cost: add. (1113->129), mult. (1604->177), div. (0->0), fcn. (1152->8), ass. (0->46)
t56 = pkin(2) + pkin(6);
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t38 = t46 * V_base(5) + t47 * V_base(4);
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t32 = t38 * t53 + t50 * V_base(6);
t37 = -t46 * V_base(4) + t47 * V_base(5);
t36 = qJD(2) - t37;
t45 = V_base(3) + qJD(1);
t24 = -pkin(1) * t37 - pkin(5) * t38 + t45;
t42 = V_base(5) * qJ(1) + V_base(1);
t43 = -V_base(4) * qJ(1) + V_base(2);
t34 = t47 * t42 + t46 * t43;
t28 = V_base(6) * pkin(5) + t34;
t19 = t24 * t53 - t50 * t28;
t55 = qJD(3) - t19;
t10 = t32 * pkin(3) - t56 * t36 + t55;
t31 = t38 * t50 - t53 * V_base(6);
t33 = -t46 * t42 + t43 * t47;
t27 = -V_base(6) * pkin(1) - t33;
t54 = -qJ(3) * t32 + t27;
t13 = t56 * t31 + t54;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t6 = t49 * t10 + t52 * t13;
t20 = t50 * t24 + t53 * t28;
t17 = -t36 * qJ(3) - t20;
t5 = t52 * t10 - t13 * t49;
t11 = -pkin(3) * t31 - t17;
t30 = qJD(4) + t32;
t51 = cos(qJ(5));
t48 = sin(qJ(5));
t29 = qJD(5) + t30;
t22 = t31 * t49 + t36 * t52;
t21 = t31 * t52 - t36 * t49;
t18 = pkin(2) * t31 + t54;
t16 = -t36 * pkin(2) + t55;
t15 = t21 * t48 + t22 * t51;
t14 = t21 * t51 - t22 * t48;
t7 = -pkin(4) * t21 + t11;
t4 = pkin(7) * t21 + t6;
t3 = pkin(4) * t30 - pkin(7) * t22 + t5;
t2 = t3 * t48 + t4 * t51;
t1 = t3 * t51 - t4 * t48;
t8 = m(2) * (t33 ^ 2 + t34 ^ 2 + t45 ^ 2) / 0.2e1 + m(3) * (t19 ^ 2 + t20 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t45 * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,1) * t38 / 0.2e1) * t38 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t30 / 0.2e1) * t30 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t45 * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t38 + Ifges(2,2) * t37 / 0.2e1) * t37 + (t11 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t30 + Ifges(5,1) * t22 / 0.2e1) * t22 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t15 / 0.2e1) * t15 + (-t11 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t22 + Ifges(5,6) * t30 + Ifges(5,2) * t21 / 0.2e1) * t21 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t29 + Ifges(6,2) * t14 / 0.2e1) * t14 + (t19 * mrSges(3,1) - t20 * mrSges(3,2) + t16 * mrSges(4,2) - t17 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t36) * t36 + (t16 * mrSges(4,1) + t27 * mrSges(3,2) - t19 * mrSges(3,3) - t18 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t32 + (-Ifges(4,4) + Ifges(3,5)) * t36) * t32 + (V_base(2) * mrSges(1,1) + t33 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t34 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t38 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t37 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t27 * mrSges(3,1) + t17 * mrSges(4,1) - t18 * mrSges(4,2) - t20 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t31 + (Ifges(4,5) - Ifges(3,6)) * t36 + (-Ifges(3,4) - Ifges(4,6)) * t32) * t31;
T = t8;
