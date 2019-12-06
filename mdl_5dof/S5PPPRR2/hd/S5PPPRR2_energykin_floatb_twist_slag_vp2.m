% Calculate kinetic energy for
% S5PPPRR2
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPPRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:22
% EndTime: 2019-12-05 14:59:23
% DurationCPUTime: 0.80s
% Computational Cost: add. (1489->132), mult. (2358->190), div. (0->0), fcn. (1808->10), ass. (0->48)
t53 = sin(pkin(7));
t60 = cos(pkin(7));
t41 = t53 * V_base(4) - t60 * V_base(5);
t42 = t53 * V_base(5) + t60 * V_base(4);
t50 = V_base(3) + qJD(1);
t32 = pkin(1) * t41 - qJ(2) * t42 + t50;
t46 = V_base(5) * qJ(1) + V_base(1);
t47 = -V_base(4) * qJ(1) + V_base(2);
t40 = t60 * t46 + t53 * t47;
t36 = V_base(6) * qJ(2) + t40;
t52 = sin(pkin(8));
t59 = cos(pkin(8));
t27 = t52 * t32 + t59 * t36;
t20 = qJ(3) * t41 + t27;
t39 = -t53 * t46 + t60 * t47;
t35 = -V_base(6) * pkin(1) + qJD(2) - t39;
t37 = t42 * t52 - t59 * V_base(6);
t38 = t59 * t42 + t52 * V_base(6);
t22 = t37 * pkin(2) - t38 * qJ(3) + t35;
t51 = sin(pkin(9));
t54 = cos(pkin(9));
t14 = t54 * t20 + t51 * t22;
t10 = pkin(5) * t37 + t14;
t26 = t59 * t32 - t52 * t36;
t19 = -t41 * pkin(2) + qJD(3) - t26;
t29 = -t38 * t51 + t41 * t54;
t30 = t38 * t54 + t41 * t51;
t12 = -t29 * pkin(3) - t30 * pkin(5) + t19;
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t6 = t58 * t10 + t56 * t12;
t13 = -t51 * t20 + t22 * t54;
t5 = -t56 * t10 + t12 * t58;
t24 = -t56 * t30 + t37 * t58;
t9 = -pkin(3) * t37 - t13;
t57 = cos(qJ(5));
t55 = sin(qJ(5));
t28 = qJD(4) - t29;
t25 = t30 * t58 + t56 * t37;
t23 = qJD(5) - t24;
t16 = t25 * t57 + t28 * t55;
t15 = -t25 * t55 + t28 * t57;
t7 = -pkin(4) * t24 - pkin(6) * t25 + t9;
t4 = pkin(6) * t28 + t6;
t3 = -t28 * pkin(4) - t5;
t2 = t4 * t57 + t55 * t7;
t1 = -t4 * t55 + t57 * t7;
t8 = m(2) * (t39 ^ 2 + t40 ^ 2 + t50 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t35 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t50 * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,1) * t42 / 0.2e1) * t42 + (t35 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t38 / 0.2e1) * t38 + (t19 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t19 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,2) * t29 / 0.2e1) * t29 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t23 + Ifges(6,1) * t16 / 0.2e1) * t16 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t28 + Ifges(5,2) * t24 / 0.2e1) * t24 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t23 + Ifges(6,2) * t15 / 0.2e1) * t15 + (t50 * mrSges(2,1) + t26 * mrSges(3,1) - t27 * mrSges(3,2) - t40 * mrSges(2,3) - Ifges(2,4) * t42 + Ifges(3,5) * t38 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t41) * t41 + (V_base(2) * mrSges(1,1) + t39 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t40 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t42 + Ifges(1,6) * V_base(5) - Ifges(2,6) * t41 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t35 * mrSges(3,1) + t13 * mrSges(4,1) - t14 * mrSges(4,2) - t27 * mrSges(3,3) - Ifges(3,4) * t38 + Ifges(4,5) * t30 - Ifges(3,6) * t41 + Ifges(4,6) * t29 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t37) * t37;
T = t8;
