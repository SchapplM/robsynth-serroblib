% Calculate kinetic energy for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:42
% EndTime: 2019-12-31 17:36:43
% DurationCPUTime: 0.80s
% Computational Cost: add. (1229->129), mult. (1858->175), div. (0->0), fcn. (1388->8), ass. (0->45)
t58 = -pkin(3) - pkin(4);
t57 = cos(qJ(2));
t43 = V_base(5) * qJ(1) + V_base(1);
t44 = -V_base(4) * qJ(1) + V_base(2);
t48 = sin(pkin(7));
t49 = cos(pkin(7));
t33 = -t43 * t48 + t49 * t44;
t38 = t48 * V_base(5) + t49 * V_base(4);
t25 = V_base(6) * pkin(1) - pkin(5) * t38 + t33;
t34 = t49 * t43 + t48 * t44;
t37 = -t48 * V_base(4) + t49 * V_base(5);
t28 = pkin(5) * t37 + t34;
t51 = sin(qJ(2));
t20 = t51 * t25 + t57 * t28;
t45 = V_base(6) + qJD(2);
t17 = qJ(3) * t45 + t20;
t31 = -t57 * t37 + t38 * t51;
t32 = t51 * t37 + t57 * t38;
t46 = V_base(3) + qJD(1);
t35 = -pkin(1) * t37 + t46;
t18 = pkin(2) * t31 - qJ(3) * t32 + t35;
t47 = sin(pkin(8));
t56 = cos(pkin(8));
t9 = t56 * t17 + t47 * t18;
t19 = t57 * t25 - t51 * t28;
t7 = t31 * qJ(4) + t9;
t8 = -t47 * t17 + t56 * t18;
t55 = pkin(2) * t45 - qJD(3) + t19;
t54 = qJD(4) - t8;
t22 = t56 * t32 + t47 * t45;
t53 = qJ(4) * t22 + t55;
t52 = cos(qJ(5));
t50 = sin(qJ(5));
t30 = qJD(5) - t31;
t21 = t32 * t47 - t56 * t45;
t12 = t50 * t21 + t22 * t52;
t11 = t21 * t52 - t50 * t22;
t10 = pkin(3) * t21 - t53;
t6 = -t31 * pkin(3) + t54;
t5 = t58 * t21 + t53;
t4 = pkin(6) * t21 + t7;
t3 = -t22 * pkin(6) + t58 * t31 + t54;
t2 = t50 * t3 + t4 * t52;
t1 = t3 * t52 - t50 * t4;
t13 = m(2) * (t33 ^ 2 + t34 ^ 2 + t46 ^ 2) / 0.2e1 + m(3) * (t19 ^ 2 + t20 ^ 2 + t35 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(4) * (t55 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t19 * mrSges(3,1) - t20 * mrSges(3,2) + Ifges(3,3) * t45 / 0.2e1) * t45 + (t46 * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,1) * t38 / 0.2e1) * t38 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t30 / 0.2e1) * t30 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t46 * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t38 + Ifges(2,2) * t37 / 0.2e1) * t37 + (t35 * mrSges(3,2) - t19 * mrSges(3,3) + Ifges(3,5) * t45 + Ifges(3,1) * t32 / 0.2e1) * t32 + (t5 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t30 + Ifges(6,1) * t12 / 0.2e1) * t12 + (-t5 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t12 + Ifges(6,6) * t30 + Ifges(6,2) * t11 / 0.2e1) * t11 + (-t55 * mrSges(4,2) + t6 * mrSges(5,2) - t8 * mrSges(4,3) - t10 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t22) * t22 + (-t55 * mrSges(4,1) + t10 * mrSges(5,1) - t7 * mrSges(5,2) - t9 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t21 + (-Ifges(4,4) + Ifges(5,5)) * t22) * t21 + (V_base(2) * mrSges(1,1) + t33 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t34 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t38 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t37 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t35 * mrSges(3,1) + t8 * mrSges(4,1) - t6 * mrSges(5,1) - t9 * mrSges(4,2) - t20 * mrSges(3,3) + t7 * mrSges(5,3) - Ifges(3,4) * t32 - Ifges(3,6) * t45 + (Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t31 + (Ifges(5,4) + Ifges(4,5)) * t22 + (-Ifges(4,6) + Ifges(5,6)) * t21) * t31;
T = t13;
