% Calculate kinetic energy for
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:13
% EndTime: 2019-03-09 01:55:14
% DurationCPUTime: 0.83s
% Computational Cost: add. (1643->150), mult. (2113->194), div. (0->0), fcn. (1504->8), ass. (0->52)
t63 = pkin(4) + pkin(8);
t62 = cos(qJ(1));
t61 = cos(qJ(4));
t60 = pkin(1) + qJ(3);
t53 = sin(qJ(1));
t40 = t53 * V_base(5) + t62 * V_base(4);
t48 = V_base(6) + qJD(1);
t44 = pkin(6) * V_base(5) + V_base(1);
t45 = -pkin(6) * V_base(4) + V_base(2);
t36 = -t44 * t53 + t45 * t62;
t58 = qJD(2) - t36;
t26 = pkin(2) * t40 - t48 * t60 + t58;
t39 = t53 * V_base(4) - t62 * V_base(5);
t59 = -qJ(2) * t40 + V_base(3);
t28 = t39 * t60 + t59;
t49 = sin(pkin(9));
t50 = cos(pkin(9));
t16 = t26 * t50 - t28 * t49;
t35 = t39 * t49 + t48 * t50;
t12 = pkin(3) * t40 - pkin(7) * t35 + t16;
t17 = t26 * t49 + t28 * t50;
t34 = t39 * t50 - t48 * t49;
t15 = pkin(7) * t34 + t17;
t52 = sin(qJ(4));
t8 = t12 * t52 + t15 * t61;
t37 = t44 * t62 + t45 * t53;
t33 = -qJ(2) * t48 - t37;
t38 = qJD(4) + t40;
t6 = -qJ(5) * t38 - t8;
t7 = t12 * t61 - t15 * t52;
t57 = qJD(5) - t7;
t23 = t34 * t52 + t35 * t61;
t29 = -pkin(2) * t39 + qJD(3) - t33;
t20 = -pkin(3) * t34 + t29;
t56 = -qJ(5) * t23 + t20;
t55 = V_base(3) ^ 2;
t54 = cos(qJ(6));
t51 = sin(qJ(6));
t31 = -pkin(1) * t48 + t58;
t30 = pkin(1) * t39 + t59;
t22 = -t34 * t61 + t35 * t52;
t21 = qJD(6) + t23;
t19 = t22 * t51 + t38 * t54;
t18 = t22 * t54 - t38 * t51;
t10 = pkin(4) * t22 + t56;
t9 = t22 * t63 + t56;
t5 = -pkin(4) * t38 + t57;
t4 = -pkin(5) * t22 - t6;
t3 = t23 * pkin(5) - t38 * t63 + t57;
t2 = t3 * t51 + t54 * t9;
t1 = t3 * t54 - t51 * t9;
t11 = m(2) * (t36 ^ 2 + t37 ^ 2 + t55) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t55) / 0.2e1 + m(3) * (t30 ^ 2 + t31 ^ 2 + t33 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t29 ^ 2) / 0.2e1 + m(5) * (t20 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t29 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,1) * t35 / 0.2e1) * t35 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t21 / 0.2e1) * t21 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t29 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t35 + Ifges(4,2) * t34 / 0.2e1) * t34 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t21 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t21 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t36 * mrSges(2,1) - t37 * mrSges(2,2) + t31 * mrSges(3,2) - t33 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t48) * t48 + (t7 * mrSges(5,1) - t8 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t38) * t38 + (V_base(3) * mrSges(2,1) + t33 * mrSges(3,1) - t30 * mrSges(3,2) - t37 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t39 + (Ifges(3,5) - Ifges(2,6)) * t48) * t39 + (t5 * mrSges(6,1) + t20 * mrSges(5,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t23 + (-Ifges(6,4) + Ifges(5,5)) * t38) * t23 + (t20 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t22 + (Ifges(6,5) - Ifges(5,6)) * t38 + (-Ifges(5,4) - Ifges(6,6)) * t23) * t22 + (t31 * mrSges(3,1) + t16 * mrSges(4,1) + V_base(3) * mrSges(2,2) - t17 * mrSges(4,2) - t36 * mrSges(2,3) - t30 * mrSges(3,3) + Ifges(4,5) * t35 + Ifges(4,6) * t34 + (Ifges(2,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t40 + (-Ifges(3,4) + Ifges(2,5)) * t48 + (-Ifges(2,4) - Ifges(3,6)) * t39) * t40;
T  = t11;
