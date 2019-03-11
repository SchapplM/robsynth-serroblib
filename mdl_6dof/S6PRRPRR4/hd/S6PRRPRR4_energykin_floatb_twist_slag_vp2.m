% Calculate kinetic energy for
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:09
% EndTime: 2019-03-08 22:10:11
% DurationCPUTime: 1.16s
% Computational Cost: add. (3517->160), mult. (6001->225), div. (0->0), fcn. (5006->12), ass. (0->62)
t66 = sin(pkin(11));
t68 = cos(pkin(11));
t58 = t66 * V_base(5) + t68 * V_base(4);
t83 = pkin(7) * t58;
t63 = V_base(5) * qJ(1) + V_base(1);
t64 = -V_base(4) * qJ(1) + V_base(2);
t54 = -t63 * t66 + t68 * t64;
t69 = cos(pkin(6));
t46 = V_base(6) * pkin(1) - t69 * t83 + t54;
t57 = -t66 * V_base(4) + t68 * V_base(5);
t65 = V_base(3) + qJD(1);
t67 = sin(pkin(6));
t50 = -pkin(1) * t57 - t67 * t83 + t65;
t31 = -t46 * t67 + t69 * t50;
t73 = sin(qJ(2));
t76 = cos(qJ(2));
t79 = t67 * V_base(6);
t80 = t69 * t76;
t44 = t57 * t80 - t58 * t73 + t76 * t79;
t77 = t57 * t69 + t79;
t45 = t58 * t76 + t73 * t77;
t21 = -pkin(2) * t44 - pkin(8) * t45 + t31;
t55 = t68 * t63 + t66 * t64;
t42 = pkin(7) * t77 + t55;
t81 = t50 * t67;
t30 = t76 * t42 + (t46 * t69 + t81) * t73;
t53 = -t57 * t67 + t69 * V_base(6) + qJD(2);
t25 = pkin(8) * t53 + t30;
t72 = sin(qJ(3));
t82 = cos(qJ(3));
t16 = t72 * t21 + t82 * t25;
t43 = qJD(3) - t44;
t14 = t43 * qJ(4) + t16;
t34 = t45 * t72 - t53 * t82;
t11 = pkin(9) * t34 + t14;
t71 = sin(qJ(5));
t75 = cos(qJ(5));
t35 = t45 * t82 + t72 * t53;
t15 = t21 * t82 - t72 * t25;
t78 = qJD(4) - t15;
t9 = -t35 * pkin(9) + (-pkin(3) - pkin(4)) * t43 + t78;
t6 = t75 * t11 + t71 * t9;
t29 = -t73 * t42 + t46 * t80 + t76 * t81;
t24 = -t53 * pkin(2) - t29;
t5 = -t11 * t71 + t75 * t9;
t27 = t34 * t75 - t35 * t71;
t17 = t34 * pkin(3) - t35 * qJ(4) + t24;
t13 = -pkin(4) * t34 - t17;
t74 = cos(qJ(6));
t70 = sin(qJ(6));
t41 = qJD(5) - t43;
t28 = t34 * t71 + t35 * t75;
t26 = qJD(6) - t27;
t19 = t28 * t74 + t41 * t70;
t18 = -t28 * t70 + t41 * t74;
t12 = -t43 * pkin(3) + t78;
t7 = -pkin(5) * t27 - pkin(10) * t28 + t13;
t4 = pkin(10) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t4 * t74 + t7 * t70;
t1 = -t4 * t70 + t7 * t74;
t8 = m(2) * (t54 ^ 2 + t55 ^ 2 + t65 ^ 2) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t24 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t14 ^ 2 + t17 ^ 2) / 0.2e1 + m(6) * (t13 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t65 * mrSges(2,2) - t54 * mrSges(2,3) + Ifges(2,1) * t58 / 0.2e1) * t58 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t53 / 0.2e1) * t53 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t41 / 0.2e1) * t41 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t65 * mrSges(2,1) + t55 * mrSges(2,3) + Ifges(2,4) * t58 + Ifges(2,2) * t57 / 0.2e1) * t57 + (t31 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t53 + Ifges(3,1) * t45 / 0.2e1) * t45 + (t13 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t41 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t26 + Ifges(7,1) * t19 / 0.2e1) * t19 + (-t31 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t45 + Ifges(3,6) * t53 + Ifges(3,2) * t44 / 0.2e1) * t44 + (-t13 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t41 + Ifges(6,2) * t27 / 0.2e1) * t27 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t26 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t15 * mrSges(4,1) - t12 * mrSges(5,1) - t16 * mrSges(4,2) + t14 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t43) * t43 + (t24 * mrSges(4,2) + t12 * mrSges(5,2) - t15 * mrSges(4,3) - t17 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t35 + (Ifges(5,4) + Ifges(4,5)) * t43) * t35 + (V_base(2) * mrSges(1,1) + t54 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t55 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t58 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t57 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t24 * mrSges(4,1) + t17 * mrSges(5,1) - t14 * mrSges(5,2) - t16 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t34 + (-Ifges(4,6) + Ifges(5,6)) * t43 + (-Ifges(4,4) + Ifges(5,5)) * t35) * t34;
T  = t8;
