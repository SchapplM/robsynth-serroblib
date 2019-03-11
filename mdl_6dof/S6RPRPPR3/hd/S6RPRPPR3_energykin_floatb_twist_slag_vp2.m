% Calculate kinetic energy for
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:01
% EndTime: 2019-03-09 02:44:02
% DurationCPUTime: 0.79s
% Computational Cost: add. (1687->151), mult. (2421->194), div. (0->0), fcn. (1836->8), ass. (0->50)
t65 = -pkin(4) - pkin(8);
t49 = V_base(5) * pkin(6) + V_base(1);
t50 = -V_base(4) * pkin(6) + V_base(2);
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t39 = -t49 * t56 + t59 * t50;
t44 = t56 * V_base(5) + t59 * V_base(4);
t52 = V_base(6) + qJD(1);
t31 = pkin(1) * t52 - qJ(2) * t44 + t39;
t40 = t59 * t49 + t56 * t50;
t43 = -t56 * V_base(4) + t59 * V_base(5);
t35 = qJ(2) * t43 + t40;
t53 = sin(pkin(9));
t64 = cos(pkin(9));
t23 = t53 * t31 + t64 * t35;
t17 = pkin(7) * t52 + t23;
t37 = t64 * t43 - t44 * t53;
t38 = t53 * t43 + t64 * t44;
t41 = -pkin(1) * t43 + qJD(2) + V_base(3);
t19 = -pkin(2) * t37 - pkin(7) * t38 + t41;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t12 = t58 * t17 + t55 * t19;
t22 = t64 * t31 - t53 * t35;
t36 = qJD(3) - t37;
t10 = t36 * qJ(4) + t12;
t16 = -t52 * pkin(2) - t22;
t11 = -t55 * t17 + t19 * t58;
t29 = t38 * t55 - t58 * t52;
t30 = t38 * t58 + t52 * t55;
t13 = t29 * pkin(3) - t30 * qJ(4) + t16;
t63 = qJD(4) - t11;
t7 = -qJ(5) * t29 - t10;
t62 = qJD(5) - t13;
t61 = -qJ(5) * t30 + t63;
t60 = V_base(3) ^ 2;
t57 = cos(qJ(6));
t54 = sin(qJ(6));
t26 = qJD(6) + t30;
t21 = t29 * t57 - t36 * t54;
t20 = -t29 * t54 - t36 * t57;
t9 = -pkin(3) * t36 + t63;
t8 = -pkin(4) * t29 + t62;
t6 = pkin(5) * t36 - t7;
t5 = (-pkin(3) - pkin(4)) * t36 + t61;
t4 = pkin(5) * t30 + t65 * t29 + t62;
t3 = (-pkin(3) + t65) * t36 + t61;
t2 = t3 * t57 + t4 * t54;
t1 = -t3 * t54 + t4 * t57;
t14 = m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t60) / 0.2e1 + m(2) * (t39 ^ 2 + t40 ^ 2 + t60) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t41 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t16 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t13 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,1) * t44 / 0.2e1) * t44 + (t41 * mrSges(3,2) - t22 * mrSges(3,3) + Ifges(3,1) * t38 / 0.2e1) * t38 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t40 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,2) * t43 / 0.2e1) * t43 + (-t41 * mrSges(3,1) + t23 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,2) * t37 / 0.2e1) * t37 + (t6 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t26 + Ifges(7,1) * t21 / 0.2e1) * t21 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t6 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t26 + Ifges(7,2) * t20 / 0.2e1) * t20 + (t39 * mrSges(2,1) + t22 * mrSges(3,1) - t40 * mrSges(2,2) - t23 * mrSges(3,2) + Ifges(2,5) * t44 + Ifges(3,5) * t38 + Ifges(2,6) * t43 + Ifges(3,6) * t37 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t52) * t52 + (t11 * mrSges(4,1) - t9 * mrSges(5,1) - t7 * mrSges(6,1) - t12 * mrSges(4,2) + t5 * mrSges(6,2) + t10 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t36) * t36 + (t8 * mrSges(6,1) + t16 * mrSges(4,2) + t9 * mrSges(5,2) - t11 * mrSges(4,3) - t13 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t30 + (Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * t36) * t30 + (t16 * mrSges(4,1) + t13 * mrSges(5,1) - t10 * mrSges(5,2) + t8 * mrSges(6,2) - t12 * mrSges(4,3) - t7 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t29 + (-Ifges(6,5) - Ifges(4,6) + Ifges(5,6)) * t36 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t30) * t29;
T  = t14;
