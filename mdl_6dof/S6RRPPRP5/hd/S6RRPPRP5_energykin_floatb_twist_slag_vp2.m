% Calculate kinetic energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:29
% EndTime: 2019-03-09 08:41:30
% DurationCPUTime: 0.79s
% Computational Cost: add. (1977->149), mult. (2497->194), div. (0->0), fcn. (1868->8), ass. (0->49)
t53 = sin(qJ(5));
t61 = cos(qJ(5));
t55 = sin(qJ(1));
t56 = cos(qJ(1));
t43 = t55 * V_base(5) + t56 * V_base(4);
t50 = V_base(6) + qJD(1);
t54 = sin(qJ(2));
t62 = cos(qJ(2));
t37 = t43 * t62 + t50 * t54;
t42 = -t55 * V_base(4) + t56 * V_base(5);
t41 = qJD(2) - t42;
t30 = -pkin(1) * t42 - pkin(7) * t43 + V_base(3);
t48 = pkin(6) * V_base(5) + V_base(1);
t49 = -pkin(6) * V_base(4) + V_base(2);
t39 = t48 * t56 + t49 * t55;
t34 = pkin(7) * t50 + t39;
t24 = t30 * t62 - t34 * t54;
t59 = qJD(3) - t24;
t60 = pkin(2) + qJ(4);
t15 = pkin(3) * t37 - t41 * t60 + t59;
t36 = t43 * t54 - t50 * t62;
t38 = -t48 * t55 + t49 * t56;
t33 = -pkin(1) * t50 - t38;
t58 = -qJ(3) * t37 + t33;
t18 = t36 * t60 + t58;
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t10 = t15 * t52 - t18 * t51;
t28 = t36 * t51 + t41 * t52;
t7 = pkin(4) * t37 - pkin(8) * t28 + t10;
t11 = t15 * t51 + t18 * t52;
t27 = t36 * t52 - t41 * t51;
t9 = pkin(8) * t27 + t11;
t4 = t53 * t7 + t61 * t9;
t25 = t30 * t54 + t34 * t62;
t22 = -qJ(3) * t41 - t25;
t3 = -t53 * t9 + t61 * t7;
t17 = -pkin(3) * t36 + qJD(4) - t22;
t12 = -pkin(4) * t27 + t17;
t57 = V_base(3) ^ 2;
t35 = qJD(5) + t37;
t23 = pkin(2) * t36 + t58;
t21 = -pkin(2) * t41 + t59;
t20 = t27 * t53 + t28 * t61;
t19 = -t27 * t61 + t28 * t53;
t5 = pkin(5) * t19 - qJ(6) * t20 + t12;
t2 = qJ(6) * t35 + t4;
t1 = -pkin(5) * t35 + qJD(6) - t3;
t6 = m(2) * (t38 ^ 2 + t39 ^ 2 + t57) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t33 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t23 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t17 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + Ifges(2,3) * t50 / 0.2e1) * t50 + (t17 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,1) * t28 / 0.2e1) * t28 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,5) * t50 + Ifges(2,1) * t43 / 0.2e1) * t43 + (-t17 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t28 + Ifges(5,2) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,4) * t43 + Ifges(2,6) * t50 + Ifges(2,2) * t42 / 0.2e1) * t42 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + t21 * mrSges(4,2) - t22 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t41) * t41 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35) * t35 + (t33 * mrSges(3,1) + t22 * mrSges(4,1) - t23 * mrSges(4,2) - t25 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t36 + (Ifges(4,5) - Ifges(3,6)) * t41) * t36 + (t12 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t20 + (Ifges(7,4) + Ifges(6,5)) * t35) * t20 + (t12 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t19 + (-Ifges(6,6) + Ifges(7,6)) * t35 + (-Ifges(6,4) + Ifges(7,5)) * t20) * t19 + (t21 * mrSges(4,1) + t10 * mrSges(5,1) + t33 * mrSges(3,2) - t11 * mrSges(5,2) - t24 * mrSges(3,3) - t23 * mrSges(4,3) + Ifges(5,5) * t28 + Ifges(5,6) * t27 + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t37 + (-Ifges(4,4) + Ifges(3,5)) * t41 + (-Ifges(3,4) - Ifges(4,6)) * t36) * t37;
T  = t6;
