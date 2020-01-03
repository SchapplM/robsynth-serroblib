% Calculate kinetic energy for
% S5RPRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:50
% EndTime: 2019-12-31 18:18:51
% DurationCPUTime: 0.82s
% Computational Cost: add. (1861->132), mult. (2722->191), div. (0->0), fcn. (2132->10), ass. (0->49)
t46 = V_base(5) * pkin(5) + V_base(1);
t47 = -V_base(4) * pkin(5) + V_base(2);
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t38 = -t46 * t55 + t58 * t47;
t42 = t55 * V_base(5) + t58 * V_base(4);
t48 = V_base(6) + qJD(1);
t31 = pkin(1) * t48 - qJ(2) * t42 + t38;
t39 = t58 * t46 + t55 * t47;
t41 = -t55 * V_base(4) + t58 * V_base(5);
t34 = qJ(2) * t41 + t39;
t50 = sin(pkin(8));
t52 = cos(pkin(8));
t27 = t50 * t31 + t52 * t34;
t22 = pkin(6) * t48 + t27;
t36 = t41 * t52 - t42 * t50;
t37 = t41 * t50 + t42 * t52;
t40 = -pkin(1) * t41 + qJD(2) + V_base(3);
t25 = -pkin(2) * t36 - pkin(6) * t37 + t40;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t13 = t57 * t22 + t54 * t25;
t29 = -t37 * t54 + t48 * t57;
t11 = qJ(4) * t29 + t13;
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t12 = -t22 * t54 + t57 * t25;
t30 = t37 * t57 + t48 * t54;
t35 = qJD(3) - t36;
t9 = pkin(3) * t35 - qJ(4) * t30 + t12;
t6 = t51 * t11 + t49 * t9;
t26 = t31 * t52 - t50 * t34;
t5 = -t11 * t49 + t51 * t9;
t18 = t29 * t51 - t30 * t49;
t21 = -pkin(2) * t48 - t26;
t16 = -pkin(3) * t29 + qJD(4) + t21;
t59 = V_base(3) ^ 2;
t56 = cos(qJ(5));
t53 = sin(qJ(5));
t19 = t29 * t49 + t30 * t51;
t17 = qJD(5) - t18;
t15 = t19 * t56 + t35 * t53;
t14 = -t19 * t53 + t35 * t56;
t7 = -pkin(4) * t18 - pkin(7) * t19 + t16;
t4 = pkin(7) * t35 + t6;
t3 = -pkin(4) * t35 - t5;
t2 = t4 * t56 + t53 * t7;
t1 = -t4 * t53 + t56 * t7;
t8 = m(2) * (t38 ^ 2 + t39 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t16 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,1) * t42 / 0.2e1) * t42 + (t40 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t37 / 0.2e1) * t37 + (t21 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (t16 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t19 / 0.2e1) * t19 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t17 / 0.2e1) * t17 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,4) * t42 + Ifges(2,2) * t41 / 0.2e1) * t41 + (-t40 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t21 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,2) * t29 / 0.2e1) * t29 + (-t16 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t19 + Ifges(5,2) * t18 / 0.2e1) * t18 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t17 + Ifges(6,1) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t17 + Ifges(6,2) * t14 / 0.2e1) * t14 + (t38 * mrSges(2,1) + t26 * mrSges(3,1) - t39 * mrSges(2,2) - t27 * mrSges(3,2) + Ifges(2,5) * t42 + Ifges(3,5) * t37 + Ifges(2,6) * t41 + Ifges(3,6) * t36 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t48) * t48 + (t12 * mrSges(4,1) + t5 * mrSges(5,1) - t13 * mrSges(4,2) - t6 * mrSges(5,2) + Ifges(4,5) * t30 + Ifges(5,5) * t19 + Ifges(4,6) * t29 + Ifges(5,6) * t18 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t35) * t35;
T = t8;
