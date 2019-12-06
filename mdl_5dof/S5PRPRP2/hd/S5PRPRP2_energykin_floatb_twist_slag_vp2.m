% Calculate kinetic energy for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:26
% EndTime: 2019-12-05 15:30:27
% DurationCPUTime: 0.73s
% Computational Cost: add. (1389->128), mult. (2108->175), div. (0->0), fcn. (1604->8), ass. (0->42)
t43 = V_base(5) * qJ(1) + V_base(1);
t44 = -V_base(4) * qJ(1) + V_base(2);
t48 = sin(pkin(7));
t50 = cos(pkin(7));
t34 = -t43 * t48 + t50 * t44;
t39 = t48 * V_base(5) + t50 * V_base(4);
t28 = V_base(6) * pkin(1) - pkin(5) * t39 + t34;
t35 = t50 * t43 + t48 * t44;
t38 = -t48 * V_base(4) + t50 * V_base(5);
t31 = pkin(5) * t38 + t35;
t52 = sin(qJ(2));
t54 = cos(qJ(2));
t22 = t54 * t28 - t52 * t31;
t45 = V_base(6) + qJD(2);
t17 = -t45 * pkin(2) + qJD(3) - t22;
t33 = t52 * t38 + t54 * t39;
t47 = sin(pkin(8));
t49 = cos(pkin(8));
t25 = -t33 * t47 + t45 * t49;
t26 = t33 * t49 + t45 * t47;
t13 = -t25 * pkin(3) - t26 * pkin(6) + t17;
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t23 = t52 * t28 + t54 * t31;
t18 = qJ(3) * t45 + t23;
t32 = -t54 * t38 + t39 * t52;
t46 = V_base(3) + qJD(1);
t36 = -pkin(1) * t38 + t46;
t19 = pkin(2) * t32 - qJ(3) * t33 + t36;
t12 = t49 * t18 + t47 * t19;
t8 = pkin(6) * t32 + t12;
t4 = t51 * t13 + t53 * t8;
t3 = t53 * t13 - t51 * t8;
t11 = -t47 * t18 + t19 * t49;
t7 = -pkin(3) * t32 - t11;
t24 = qJD(4) - t25;
t21 = t26 * t53 + t51 * t32;
t20 = -t51 * t26 + t32 * t53;
t5 = -pkin(4) * t20 + qJD(5) + t7;
t2 = qJ(5) * t20 + t4;
t1 = pkin(4) * t24 - qJ(5) * t21 + t3;
t6 = m(2) * (t34 ^ 2 + t35 ^ 2 + t46 ^ 2) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t36 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t17 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t22 * mrSges(3,1) - t23 * mrSges(3,2) + Ifges(3,3) * t45 / 0.2e1) * t45 + (t46 * mrSges(2,2) - t34 * mrSges(2,3) + Ifges(2,1) * t39 / 0.2e1) * t39 + (t17 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,1) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t46 * mrSges(2,1) + t35 * mrSges(2,3) + Ifges(2,4) * t39 + Ifges(2,2) * t38 / 0.2e1) * t38 + (t36 * mrSges(3,2) - t22 * mrSges(3,3) + Ifges(3,5) * t45 + Ifges(3,1) * t33 / 0.2e1) * t33 + (-t17 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t26 + Ifges(4,2) * t25 / 0.2e1) * t25 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t24) * t24 + (t7 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t21 + (Ifges(5,5) + Ifges(6,5)) * t24) * t21 + (V_base(2) * mrSges(1,1) + t34 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t35 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t39 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t38 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t36 * mrSges(3,1) + t11 * mrSges(4,1) - t12 * mrSges(4,2) - t23 * mrSges(3,3) - Ifges(3,4) * t33 + Ifges(4,5) * t26 - Ifges(3,6) * t45 + Ifges(4,6) * t25 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t32) * t32 + (-t7 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t20 + (Ifges(5,6) + Ifges(6,6)) * t24 + (Ifges(5,4) + Ifges(6,4)) * t21) * t20;
T = t6;
