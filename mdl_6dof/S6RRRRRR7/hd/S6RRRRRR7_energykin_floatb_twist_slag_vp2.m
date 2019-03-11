% Calculate kinetic energy for
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:28:55
% EndTime: 2019-03-10 04:28:56
% DurationCPUTime: 1.71s
% Computational Cost: add. (6799->161), mult. (10173->240), div. (0->0), fcn. (8714->14), ass. (0->67)
t64 = V_base(5) * pkin(7) + V_base(1);
t65 = -V_base(4) * pkin(7) + V_base(2);
t74 = sin(qJ(1));
t80 = cos(qJ(1));
t57 = -t64 * t74 + t80 * t65;
t66 = V_base(6) + qJD(1);
t68 = cos(pkin(6));
t60 = t74 * V_base(5) + t80 * V_base(4);
t87 = pkin(8) * t60;
t52 = pkin(1) * t66 - t68 * t87 + t57;
t59 = -t74 * V_base(4) + t80 * V_base(5);
t67 = sin(pkin(6));
t55 = -pkin(1) * t59 - t67 * t87 + V_base(3);
t88 = t52 * t68 + t55 * t67;
t58 = t80 * t64 + t74 * t65;
t82 = t59 * t68 + t66 * t67;
t49 = t82 * pkin(8) + t58;
t73 = sin(qJ(2));
t79 = cos(qJ(2));
t34 = -t73 * t49 + t88 * t79;
t50 = -t60 * t73 + t82 * t79;
t38 = -t52 * t67 + t68 * t55;
t51 = t60 * t79 + t82 * t73;
t29 = -pkin(2) * t50 - pkin(9) * t51 + t38;
t35 = t79 * t49 + t88 * t73;
t56 = -t59 * t67 + t66 * t68 + qJD(2);
t33 = pkin(9) * t56 + t35;
t72 = sin(qJ(3));
t78 = cos(qJ(3));
t22 = t72 * t29 + t78 * t33;
t48 = qJD(3) - t50;
t20 = pkin(10) * t48 + t22;
t32 = -pkin(2) * t56 - t34;
t42 = -t72 * t51 + t56 * t78;
t43 = t51 * t78 + t56 * t72;
t25 = -pkin(3) * t42 - pkin(10) * t43 + t32;
t71 = sin(qJ(4));
t77 = cos(qJ(4));
t14 = t77 * t20 + t71 * t25;
t36 = -t43 * t71 + t48 * t77;
t11 = pkin(11) * t36 + t14;
t70 = sin(qJ(5));
t76 = cos(qJ(5));
t13 = -t20 * t71 + t77 * t25;
t37 = t43 * t77 + t48 * t71;
t41 = qJD(4) - t42;
t9 = pkin(4) * t41 - pkin(11) * t37 + t13;
t6 = t76 * t11 + t70 * t9;
t5 = -t11 * t70 + t76 * t9;
t21 = t29 * t78 - t72 * t33;
t19 = -pkin(3) * t48 - t21;
t40 = qJD(5) + t41;
t17 = -pkin(4) * t36 + t19;
t81 = V_base(3) ^ 2;
t75 = cos(qJ(6));
t69 = sin(qJ(6));
t39 = qJD(6) + t40;
t27 = t36 * t70 + t37 * t76;
t26 = t36 * t76 - t37 * t70;
t16 = t26 * t69 + t27 * t75;
t15 = t26 * t75 - t27 * t69;
t12 = -pkin(5) * t26 + t17;
t4 = pkin(12) * t26 + t6;
t3 = pkin(5) * t40 - pkin(12) * t27 + t5;
t2 = t3 * t69 + t4 * t75;
t1 = t3 * t75 - t4 * t69;
t7 = (t38 * mrSges(3,2) - t34 * mrSges(3,3) + Ifges(3,5) * t56 + Ifges(3,1) * t51 / 0.2e1) * t51 + (-t17 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t40 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t19 * mrSges(5,1) + t14 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,6) * t41 + Ifges(5,2) * t36 / 0.2e1) * t36 + (-t12 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t39 + Ifges(7,2) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t39 / 0.2e1) * t39 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t40 / 0.2e1) * t40 + (-V_base(3) * mrSges(2,1) + t58 * mrSges(2,3) + Ifges(2,4) * t60 + Ifges(2,6) * t66 + Ifges(2,2) * t59 / 0.2e1) * t59 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (t12 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t39 + Ifges(7,1) * t16 / 0.2e1) * t16 + (t34 * mrSges(3,1) - t35 * mrSges(3,2) + Ifges(3,3) * t56 / 0.2e1) * t56 + (t32 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t43 / 0.2e1) * t43 + (t13 * mrSges(5,1) - t14 * mrSges(5,2) + Ifges(5,3) * t41 / 0.2e1) * t41 + (t17 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t40 + Ifges(6,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(2,2) - t57 * mrSges(2,3) + Ifges(2,5) * t66 + Ifges(2,1) * t60 / 0.2e1) * t60 + (t19 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,5) * t41 + Ifges(5,1) * t37 / 0.2e1) * t37 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t57 * mrSges(2,1) - t58 * mrSges(2,2) + Ifges(2,3) * t66 / 0.2e1) * t66 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t38 * mrSges(3,1) + t35 * mrSges(3,3) + Ifges(3,4) * t51 + Ifges(3,6) * t56 + Ifges(3,2) * t50 / 0.2e1) * t50 + (-t32 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t48 + Ifges(4,2) * t42 / 0.2e1) * t42 + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) / 0.2e1 + m(6) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t34 ^ 2 + t35 ^ 2 + t38 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t81) / 0.2e1 + m(2) * (t57 ^ 2 + t58 ^ 2 + t81) / 0.2e1;
T  = t7;
