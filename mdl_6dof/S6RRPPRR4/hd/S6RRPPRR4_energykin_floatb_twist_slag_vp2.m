% Calculate kinetic energy for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:04
% EndTime: 2019-03-09 09:01:05
% DurationCPUTime: 1.22s
% Computational Cost: add. (4017->158), mult. (6325->221), div. (0->0), fcn. (5318->12), ass. (0->62)
t56 = V_base(5) * pkin(7) + V_base(1);
t57 = -V_base(4) * pkin(7) + V_base(2);
t65 = sin(qJ(1));
t69 = cos(qJ(1));
t49 = -t56 * t65 + t69 * t57;
t58 = V_base(6) + qJD(1);
t61 = cos(pkin(6));
t52 = t65 * V_base(5) + t69 * V_base(4);
t77 = pkin(8) * t52;
t42 = pkin(1) * t58 - t61 * t77 + t49;
t51 = -t65 * V_base(4) + t69 * V_base(5);
t60 = sin(pkin(6));
t46 = -pkin(1) * t51 - t60 * t77 + V_base(3);
t79 = t42 * t61 + t46 * t60;
t78 = pkin(3) + pkin(9);
t64 = sin(qJ(2));
t68 = cos(qJ(2));
t73 = t51 * t61 + t58 * t60;
t40 = -t52 * t64 + t68 * t73;
t41 = t52 * t68 + t64 * t73;
t59 = sin(pkin(11));
t74 = cos(pkin(11));
t32 = -t40 * t74 + t41 * t59;
t34 = -t42 * t60 + t61 * t46;
t28 = -pkin(2) * t40 + qJD(3) + t34;
t33 = t59 * t40 + t41 * t74;
t71 = -qJ(4) * t33 + t28;
t14 = t32 * t78 + t71;
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t48 = -t51 * t60 + t58 * t61 + qJD(2);
t50 = t69 * t56 + t65 * t57;
t39 = pkin(8) * t73 + t50;
t25 = -t39 * t64 + t79 * t68;
t21 = pkin(2) * t48 - qJ(3) * t41 + t25;
t26 = t68 * t39 + t79 * t64;
t24 = qJ(3) * t40 + t26;
t15 = t21 * t74 - t59 * t24;
t72 = qJD(4) - t15;
t9 = t33 * pkin(4) - t48 * t78 + t72;
t6 = t67 * t14 + t63 * t9;
t16 = t59 * t21 + t74 * t24;
t13 = -t48 * qJ(4) - t16;
t5 = -t14 * t63 + t67 * t9;
t29 = t32 * t67 - t48 * t63;
t10 = -pkin(4) * t32 - t13;
t70 = V_base(3) ^ 2;
t66 = cos(qJ(6));
t62 = sin(qJ(6));
t31 = qJD(5) + t33;
t30 = t32 * t63 + t48 * t67;
t27 = qJD(6) - t29;
t19 = t30 * t66 + t31 * t62;
t18 = -t30 * t62 + t31 * t66;
t17 = pkin(3) * t32 + t71;
t11 = -t48 * pkin(3) + t72;
t7 = -pkin(5) * t29 - pkin(10) * t30 + t10;
t4 = pkin(10) * t31 + t6;
t3 = -pkin(5) * t31 - t5;
t2 = t4 * t66 + t62 * t7;
t1 = -t4 * t62 + t66 * t7;
t8 = m(2) * (t49 ^ 2 + t50 ^ 2 + t70) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t70) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t13 ^ 2 + t17 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t49 * mrSges(2,1) - t50 * mrSges(2,2) + Ifges(2,3) * t58 / 0.2e1) * t58 + (t34 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,1) * t41 / 0.2e1) * t41 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t31 / 0.2e1) * t31 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t49 * mrSges(2,3) + Ifges(2,5) * t58 + Ifges(2,1) * t52 / 0.2e1) * t52 + (-t34 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t41 + Ifges(3,2) * t40 / 0.2e1) * t40 + (t10 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t31 + Ifges(6,1) * t30 / 0.2e1) * t30 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t27 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t50 * mrSges(2,3) + Ifges(2,4) * t52 + Ifges(2,6) * t58 + Ifges(2,2) * t51 / 0.2e1) * t51 + (-t10 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t30 + Ifges(6,6) * t31 + Ifges(6,2) * t29 / 0.2e1) * t29 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t27 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t11 * mrSges(5,1) + t28 * mrSges(4,2) - t15 * mrSges(4,3) - t17 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t33) * t33 + (t28 * mrSges(4,1) + t13 * mrSges(5,1) - t17 * mrSges(5,2) - t16 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t32 + (-Ifges(4,4) - Ifges(5,6)) * t33) * t32 + (t25 * mrSges(3,1) + t15 * mrSges(4,1) - t26 * mrSges(3,2) - t16 * mrSges(4,2) + t11 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(3,5) * t41 + Ifges(3,6) * t40 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t48 + (-Ifges(5,4) + Ifges(4,5)) * t33 + (Ifges(5,5) - Ifges(4,6)) * t32) * t48;
T  = t8;
