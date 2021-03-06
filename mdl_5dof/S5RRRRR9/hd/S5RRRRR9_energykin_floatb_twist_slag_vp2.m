% Calculate kinetic energy for
% S5RRRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:38
% EndTime: 2019-12-31 22:27:39
% DurationCPUTime: 1.02s
% Computational Cost: add. (2147->132), mult. (2750->195), div. (0->0), fcn. (2144->10), ass. (0->51)
t57 = sin(qJ(1));
t62 = cos(qJ(1));
t45 = -t57 * V_base(4) + t62 * V_base(5);
t46 = t57 * V_base(5) + t62 * V_base(4);
t31 = -pkin(1) * t45 - pkin(6) * t46 + V_base(3);
t50 = V_base(5) * pkin(5) + V_base(1);
t51 = -V_base(4) * pkin(5) + V_base(2);
t42 = t62 * t50 + t57 * t51;
t52 = V_base(6) + qJD(1);
t37 = pkin(6) * t52 + t42;
t56 = sin(qJ(2));
t61 = cos(qJ(2));
t27 = t56 * t31 + t61 * t37;
t44 = qJD(2) - t45;
t24 = pkin(7) * t44 + t27;
t41 = -t57 * t50 + t51 * t62;
t36 = -pkin(1) * t52 - t41;
t39 = -t56 * t46 + t52 * t61;
t40 = t46 * t61 + t52 * t56;
t25 = -pkin(2) * t39 - pkin(7) * t40 + t36;
t55 = sin(qJ(3));
t60 = cos(qJ(3));
t16 = t60 * t24 + t55 * t25;
t28 = -t40 * t55 + t44 * t60;
t12 = pkin(8) * t28 + t16;
t54 = sin(qJ(4));
t59 = cos(qJ(4));
t15 = -t24 * t55 + t60 * t25;
t29 = t40 * t60 + t44 * t55;
t38 = qJD(3) - t39;
t9 = pkin(3) * t38 - pkin(8) * t29 + t15;
t6 = t59 * t12 + t54 * t9;
t5 = -t12 * t54 + t59 * t9;
t26 = t31 * t61 - t56 * t37;
t23 = -pkin(2) * t44 - t26;
t35 = qJD(4) + t38;
t17 = -pkin(3) * t28 + t23;
t63 = V_base(3) ^ 2;
t58 = cos(qJ(5));
t53 = sin(qJ(5));
t34 = qJD(5) + t35;
t19 = t28 * t54 + t29 * t59;
t18 = t28 * t59 - t29 * t54;
t14 = t18 * t53 + t19 * t58;
t13 = t18 * t58 - t19 * t53;
t10 = -pkin(4) * t18 + t17;
t4 = pkin(9) * t18 + t6;
t3 = pkin(4) * t35 - pkin(9) * t19 + t5;
t2 = t3 * t53 + t4 * t58;
t1 = t3 * t58 - t4 * t53;
t7 = m(2) * (t41 ^ 2 + t42 ^ 2 + t63) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t63) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t36 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t23 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) / 0.2e1 + m(5) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t41 * mrSges(2,1) - t42 * mrSges(2,2) + Ifges(2,3) * t52 / 0.2e1) * t52 + (t26 * mrSges(3,1) - t27 * mrSges(3,2) + Ifges(3,3) * t44 / 0.2e1) * t44 + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t38 / 0.2e1) * t38 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t35 / 0.2e1) * t35 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t41 * mrSges(2,3) + Ifges(2,5) * t52 + Ifges(2,1) * t46 / 0.2e1) * t46 + (t36 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,5) * t44 + Ifges(3,1) * t40 / 0.2e1) * t40 + (t23 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,5) * t38 + Ifges(4,1) * t29 / 0.2e1) * t29 + (t17 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t35 + Ifges(5,1) * t19 / 0.2e1) * t19 + (t10 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t42 * mrSges(2,3) + Ifges(2,4) * t46 + Ifges(2,6) * t52 + Ifges(2,2) * t45 / 0.2e1) * t45 + (-t36 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t40 + Ifges(3,6) * t44 + Ifges(3,2) * t39 / 0.2e1) * t39 + (-t23 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t29 + Ifges(4,6) * t38 + Ifges(4,2) * t28 / 0.2e1) * t28 + (-t17 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t19 + Ifges(5,6) * t35 + Ifges(5,2) * t18 / 0.2e1) * t18 + (-t10 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t34 + Ifges(6,2) * t13 / 0.2e1) * t13;
T = t7;
