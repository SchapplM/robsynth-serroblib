% Calculate kinetic energy for
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:16
% EndTime: 2019-03-09 04:46:17
% DurationCPUTime: 0.82s
% Computational Cost: add. (1899->149), mult. (2385->194), div. (0->0), fcn. (1728->8), ass. (0->49)
t62 = pkin(1) + pkin(7);
t51 = sin(pkin(9));
t60 = cos(pkin(9));
t54 = sin(qJ(1));
t61 = cos(qJ(1));
t42 = t54 * V_base(5) + t61 * V_base(4);
t50 = V_base(6) + qJD(1);
t46 = pkin(6) * V_base(5) + V_base(1);
t47 = -pkin(6) * V_base(4) + V_base(2);
t38 = -t46 * t54 + t47 * t61;
t58 = qJD(2) - t38;
t25 = pkin(2) * t42 - t50 * t62 + t58;
t41 = t54 * V_base(4) - t61 * V_base(5);
t59 = -qJ(2) * t42 + V_base(3);
t30 = t41 * t62 + t59;
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t17 = t25 * t53 + t30 * t56;
t40 = qJD(3) + t42;
t15 = pkin(8) * t40 + t17;
t39 = t46 * t61 + t47 * t54;
t34 = -qJ(2) * t50 - t39;
t31 = -pkin(2) * t41 - t34;
t36 = t41 * t56 - t50 * t53;
t37 = t41 * t53 + t50 * t56;
t22 = -pkin(3) * t36 - pkin(8) * t37 + t31;
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t10 = -t15 * t52 + t22 * t55;
t29 = t37 * t55 + t40 * t52;
t35 = qJD(4) - t36;
t7 = pkin(4) * t35 - qJ(5) * t29 + t10;
t11 = t15 * t55 + t22 * t52;
t28 = -t37 * t52 + t40 * t55;
t9 = qJ(5) * t28 + t11;
t4 = t51 * t7 + t60 * t9;
t16 = t25 * t56 - t30 * t53;
t14 = -pkin(3) * t40 - t16;
t3 = -t51 * t9 + t60 * t7;
t12 = -pkin(4) * t28 + qJD(5) + t14;
t57 = V_base(3) ^ 2;
t33 = -pkin(1) * t50 + t58;
t32 = pkin(1) * t41 + t59;
t19 = t28 * t51 + t29 * t60;
t18 = -t28 * t60 + t29 * t51;
t5 = pkin(5) * t18 - qJ(6) * t19 + t12;
t2 = qJ(6) * t35 + t4;
t1 = -pkin(5) * t35 + qJD(6) - t3;
t6 = m(2) * (t38 ^ 2 + t39 ^ 2 + t57) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t31 ^ 2) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t14 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t40 / 0.2e1) * t40 + (t14 * mrSges(5,2) - t10 * mrSges(5,3) + Ifges(5,1) * t29 / 0.2e1) * t29 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t31 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t40 + Ifges(4,1) * t37 / 0.2e1) * t37 + (-t14 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,2) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t31 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,6) * t40 + Ifges(4,2) * t36 / 0.2e1) * t36 + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + t33 * mrSges(3,2) - t34 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t50) * t50 + (t12 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t19) * t19 + (t33 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) - t32 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t42 + (-Ifges(3,4) + Ifges(2,5)) * t50) * t42 + (t12 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t18 + (-Ifges(6,4) + Ifges(7,5)) * t19) * t18 + (V_base(3) * mrSges(2,1) + t34 * mrSges(3,1) - t32 * mrSges(3,2) - t39 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t41 + (Ifges(3,5) - Ifges(2,6)) * t50 + (-Ifges(2,4) - Ifges(3,6)) * t42) * t41 + (t10 * mrSges(5,1) + t3 * mrSges(6,1) - t1 * mrSges(7,1) - t11 * mrSges(5,2) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + Ifges(5,5) * t29 + Ifges(5,6) * t28 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t35 + (Ifges(7,4) + Ifges(6,5)) * t19 + (-Ifges(6,6) + Ifges(7,6)) * t18) * t35;
T  = t6;
