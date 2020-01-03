% Calculate kinetic energy for
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR13_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:03
% EndTime: 2019-12-31 20:32:04
% DurationCPUTime: 0.92s
% Computational Cost: add. (2119->132), mult. (2750->193), div. (0->0), fcn. (2144->10), ass. (0->50)
t57 = sin(qJ(1));
t61 = cos(qJ(1));
t43 = -t57 * V_base(4) + t61 * V_base(5);
t44 = t57 * V_base(5) + t61 * V_base(4);
t31 = -pkin(1) * t43 - pkin(6) * t44 + V_base(3);
t49 = V_base(5) * pkin(5) + V_base(1);
t50 = -V_base(4) * pkin(5) + V_base(2);
t41 = t61 * t49 + t57 * t50;
t51 = V_base(6) + qJD(1);
t36 = pkin(6) * t51 + t41;
t56 = sin(qJ(2));
t60 = cos(qJ(2));
t27 = t56 * t31 + t60 * t36;
t42 = qJD(2) - t43;
t24 = qJ(3) * t42 + t27;
t40 = -t57 * t49 + t50 * t61;
t35 = -pkin(1) * t51 - t40;
t38 = t44 * t56 - t60 * t51;
t39 = t44 * t60 + t51 * t56;
t25 = pkin(2) * t38 - qJ(3) * t39 + t35;
t52 = sin(pkin(9));
t53 = cos(pkin(9));
t16 = t53 * t24 + t52 * t25;
t28 = -t39 * t52 + t42 * t53;
t12 = pkin(7) * t28 + t16;
t55 = sin(qJ(4));
t59 = cos(qJ(4));
t15 = -t24 * t52 + t53 * t25;
t29 = t39 * t53 + t42 * t52;
t9 = pkin(3) * t38 - pkin(7) * t29 + t15;
t6 = t59 * t12 + t55 * t9;
t5 = -t12 * t55 + t59 * t9;
t26 = t31 * t60 - t56 * t36;
t37 = qJD(4) + t38;
t23 = -pkin(2) * t42 + qJD(3) - t26;
t17 = -pkin(3) * t28 + t23;
t62 = V_base(3) ^ 2;
t58 = cos(qJ(5));
t54 = sin(qJ(5));
t34 = qJD(5) + t37;
t19 = t28 * t55 + t29 * t59;
t18 = t28 * t59 - t29 * t55;
t14 = t18 * t54 + t19 * t58;
t13 = t18 * t58 - t19 * t54;
t10 = -pkin(4) * t18 + t17;
t4 = pkin(8) * t18 + t6;
t3 = pkin(4) * t37 - pkin(8) * t19 + t5;
t2 = t3 * t54 + t4 * t58;
t1 = t3 * t58 - t4 * t54;
t7 = m(2) * (t40 ^ 2 + t41 ^ 2 + t62) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t62) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t35 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t23 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) / 0.2e1 + m(5) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t40 * mrSges(2,1) - t41 * mrSges(2,2) + Ifges(2,3) * t51 / 0.2e1) * t51 + (t26 * mrSges(3,1) - t27 * mrSges(3,2) + Ifges(3,3) * t42 / 0.2e1) * t42 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t37 / 0.2e1) * t37 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (t23 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,1) * t29 / 0.2e1) * t29 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t40 * mrSges(2,3) + Ifges(2,5) * t51 + Ifges(2,1) * t44 / 0.2e1) * t44 + (t35 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,5) * t42 + Ifges(3,1) * t39 / 0.2e1) * t39 + (-t23 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t29 + Ifges(4,2) * t28 / 0.2e1) * t28 + (t17 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t37 + Ifges(5,1) * t19 / 0.2e1) * t19 + (t10 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t41 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,6) * t51 + Ifges(2,2) * t43 / 0.2e1) * t43 + (-t17 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t19 + Ifges(5,6) * t37 + Ifges(5,2) * t18 / 0.2e1) * t18 + (-t10 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t34 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t35 * mrSges(3,1) + t15 * mrSges(4,1) - t16 * mrSges(4,2) - t27 * mrSges(3,3) - Ifges(3,4) * t39 + Ifges(4,5) * t29 - Ifges(3,6) * t42 + Ifges(4,6) * t28 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t38) * t38;
T = t7;
