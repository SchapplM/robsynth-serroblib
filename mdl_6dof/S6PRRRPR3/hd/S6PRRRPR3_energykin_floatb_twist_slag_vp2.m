% Calculate kinetic energy for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:10:51
% EndTime: 2019-03-08 23:10:52
% DurationCPUTime: 1.22s
% Computational Cost: add. (3947->158), mult. (6763->222), div. (0->0), fcn. (5696->12), ass. (0->62)
t57 = V_base(5) * qJ(1) + V_base(1);
t58 = -V_base(4) * qJ(1) + V_base(2);
t60 = sin(pkin(11));
t62 = cos(pkin(11));
t49 = -t57 * t60 + t62 * t58;
t63 = cos(pkin(6));
t53 = t60 * V_base(5) + t62 * V_base(4);
t78 = pkin(7) * t53;
t44 = V_base(6) * pkin(1) - t63 * t78 + t49;
t52 = -t60 * V_base(4) + t62 * V_base(5);
t59 = V_base(3) + qJD(1);
t61 = sin(pkin(6));
t47 = -pkin(1) * t52 - t61 * t78 + t59;
t80 = t44 * t63 + t47 * t61;
t50 = t62 * t57 + t60 * t58;
t72 = t52 * t63 + t61 * V_base(6);
t40 = pkin(7) * t72 + t50;
t67 = sin(qJ(2));
t70 = cos(qJ(2));
t30 = -t67 * t40 + t70 * t80;
t42 = -t67 * t53 + t70 * t72;
t79 = pkin(4) + pkin(10);
t77 = cos(qJ(4));
t32 = -t44 * t61 + t63 * t47;
t43 = t53 * t70 + t67 * t72;
t23 = -pkin(2) * t42 - pkin(8) * t43 + t32;
t31 = t70 * t40 + t80 * t67;
t48 = -t52 * t61 + t63 * V_base(6) + qJD(2);
t26 = pkin(8) * t48 + t31;
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t16 = t69 * t23 - t26 * t66;
t35 = t43 * t69 + t48 * t66;
t41 = qJD(3) - t42;
t12 = pkin(3) * t41 - pkin(9) * t35 + t16;
t17 = t66 * t23 + t69 * t26;
t34 = -t43 * t66 + t48 * t69;
t15 = pkin(9) * t34 + t17;
t65 = sin(qJ(4));
t8 = t65 * t12 + t77 * t15;
t39 = qJD(4) + t41;
t6 = -qJ(5) * t39 - t8;
t7 = t12 * t77 - t65 * t15;
t73 = qJD(5) - t7;
t29 = t65 * t34 + t35 * t77;
t25 = -t48 * pkin(2) - t30;
t18 = -t34 * pkin(3) + t25;
t71 = -t29 * qJ(5) + t18;
t68 = cos(qJ(6));
t64 = sin(qJ(6));
t28 = -t34 * t77 + t35 * t65;
t27 = qJD(6) + t29;
t20 = t28 * t64 + t39 * t68;
t19 = t28 * t68 - t39 * t64;
t10 = t28 * pkin(4) + t71;
t9 = t28 * t79 + t71;
t5 = -t39 * pkin(4) + t73;
t4 = -pkin(5) * t28 - t6;
t3 = t29 * pkin(5) - t39 * t79 + t73;
t2 = t3 * t64 + t68 * t9;
t1 = t3 * t68 - t64 * t9;
t11 = m(2) * (t49 ^ 2 + t50 ^ 2 + t59 ^ 2) / 0.2e1 + m(3) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t25 ^ 2) / 0.2e1 + m(5) * (t18 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t59 * mrSges(2,2) - t49 * mrSges(2,3) + Ifges(2,1) * t53 / 0.2e1) * t53 + (t30 * mrSges(3,1) - t31 * mrSges(3,2) + Ifges(3,3) * t48 / 0.2e1) * t48 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t59 * mrSges(2,1) + t50 * mrSges(2,3) + Ifges(2,4) * t53 + Ifges(2,2) * t52 / 0.2e1) * t52 + (t32 * mrSges(3,2) - t30 * mrSges(3,3) + Ifges(3,5) * t48 + Ifges(3,1) * t43 / 0.2e1) * t43 + (t25 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t35 / 0.2e1) * t35 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t27 + Ifges(7,1) * t20 / 0.2e1) * t20 + (-t32 * mrSges(3,1) + t31 * mrSges(3,3) + Ifges(3,4) * t43 + Ifges(3,6) * t48 + Ifges(3,2) * t42 / 0.2e1) * t42 + (-t25 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t35 + Ifges(4,6) * t41 + Ifges(4,2) * t34 / 0.2e1) * t34 + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t27 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t7 * mrSges(5,1) - t8 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t39) * t39 + (t5 * mrSges(6,1) + t18 * mrSges(5,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t29 + (-Ifges(6,4) + Ifges(5,5)) * t39) * t29 + (V_base(2) * mrSges(1,1) + t49 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t50 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t53 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t52 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t18 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t28 + (Ifges(6,5) - Ifges(5,6)) * t39 + (-Ifges(5,4) - Ifges(6,6)) * t29) * t28;
T  = t11;
