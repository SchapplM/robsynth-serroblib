% Calculate kinetic energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:55
% EndTime: 2019-03-09 00:47:57
% DurationCPUTime: 1.48s
% Computational Cost: add. (5921->161), mult. (10173->239), div. (0->0), fcn. (8714->14), ass. (0->66)
t64 = V_base(5) * qJ(1) + V_base(1);
t65 = -V_base(4) * qJ(1) + V_base(2);
t67 = sin(pkin(12));
t69 = cos(pkin(12));
t57 = -t64 * t67 + t69 * t65;
t70 = cos(pkin(6));
t60 = t67 * V_base(5) + t69 * V_base(4);
t86 = pkin(7) * t60;
t52 = V_base(6) * pkin(1) - t70 * t86 + t57;
t59 = -t67 * V_base(4) + t69 * V_base(5);
t66 = V_base(3) + qJD(1);
t68 = sin(pkin(6));
t55 = -pkin(1) * t59 - t68 * t86 + t66;
t87 = t52 * t70 + t55 * t68;
t58 = t69 * t64 + t67 * t65;
t81 = t59 * t70 + t68 * V_base(6);
t48 = t81 * pkin(7) + t58;
t75 = sin(qJ(2));
t80 = cos(qJ(2));
t34 = -t75 * t48 + t87 * t80;
t50 = -t75 * t60 + t81 * t80;
t38 = -t52 * t68 + t70 * t55;
t51 = t60 * t80 + t81 * t75;
t29 = -pkin(2) * t50 - pkin(8) * t51 + t38;
t35 = t80 * t48 + t87 * t75;
t56 = -t59 * t68 + t70 * V_base(6) + qJD(2);
t33 = pkin(8) * t56 + t35;
t74 = sin(qJ(3));
t79 = cos(qJ(3));
t22 = t74 * t29 + t79 * t33;
t49 = qJD(3) - t50;
t20 = pkin(9) * t49 + t22;
t32 = -t56 * pkin(2) - t34;
t42 = -t74 * t51 + t56 * t79;
t43 = t51 * t79 + t56 * t74;
t25 = -t42 * pkin(3) - t43 * pkin(9) + t32;
t73 = sin(qJ(4));
t78 = cos(qJ(4));
t14 = t78 * t20 + t73 * t25;
t36 = -t43 * t73 + t49 * t78;
t11 = pkin(10) * t36 + t14;
t72 = sin(qJ(5));
t77 = cos(qJ(5));
t13 = -t20 * t73 + t78 * t25;
t37 = t43 * t78 + t49 * t73;
t41 = qJD(4) - t42;
t9 = pkin(4) * t41 - pkin(10) * t37 + t13;
t6 = t77 * t11 + t72 * t9;
t5 = -t11 * t72 + t77 * t9;
t21 = t29 * t79 - t74 * t33;
t19 = -pkin(3) * t49 - t21;
t40 = qJD(5) + t41;
t17 = -pkin(4) * t36 + t19;
t76 = cos(qJ(6));
t71 = sin(qJ(6));
t39 = qJD(6) + t40;
t27 = t36 * t72 + t37 * t77;
t26 = t36 * t77 - t37 * t72;
t16 = t26 * t71 + t27 * t76;
t15 = t26 * t76 - t27 * t71;
t12 = -pkin(5) * t26 + t17;
t4 = pkin(11) * t26 + t6;
t3 = pkin(5) * t40 - pkin(11) * t27 + t5;
t2 = t3 * t71 + t4 * t76;
t1 = t3 * t76 - t4 * t71;
t7 = (-t66 * mrSges(2,1) + t58 * mrSges(2,3) + Ifges(2,4) * t60 + Ifges(2,2) * t59 / 0.2e1) * t59 + (-t19 * mrSges(5,1) + t14 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,6) * t41 + Ifges(5,2) * t36 / 0.2e1) * t36 + (-t38 * mrSges(3,1) + t35 * mrSges(3,3) + Ifges(3,4) * t51 + Ifges(3,6) * t56 + Ifges(3,2) * t50 / 0.2e1) * t50 + (t66 * mrSges(2,2) - t57 * mrSges(2,3) + Ifges(2,1) * t60 / 0.2e1) * t60 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t13 * mrSges(5,1) - t14 * mrSges(5,2) + Ifges(5,3) * t41 / 0.2e1) * t41 + (-t32 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t49 + Ifges(4,2) * t42 / 0.2e1) * t42 + (-t17 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t40 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t40 / 0.2e1) * t40 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t32 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t49 + Ifges(4,1) * t43 / 0.2e1) * t43 + (V_base(2) * mrSges(1,1) + t57 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t58 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t60 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t59 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-t12 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t39 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t17 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t40 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t38 * mrSges(3,2) - t34 * mrSges(3,3) + Ifges(3,5) * t56 + Ifges(3,1) * t51 / 0.2e1) * t51 + (t12 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t39 + Ifges(7,1) * t16 / 0.2e1) * t16 + (t34 * mrSges(3,1) - t35 * mrSges(3,2) + Ifges(3,3) * t56 / 0.2e1) * t56 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t49 / 0.2e1) * t49 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t39 / 0.2e1) * t39 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t57 ^ 2 + t58 ^ 2 + t66 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) / 0.2e1 + (t19 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,5) * t41 + Ifges(5,1) * t37 / 0.2e1) * t37 + m(6) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t34 ^ 2 + t35 ^ 2 + t38 ^ 2) / 0.2e1;
T  = t7;
