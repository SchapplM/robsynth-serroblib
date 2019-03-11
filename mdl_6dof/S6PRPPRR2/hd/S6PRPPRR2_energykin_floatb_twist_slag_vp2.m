% Calculate kinetic energy for
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:35
% EndTime: 2019-03-08 19:17:36
% DurationCPUTime: 1.13s
% Computational Cost: add. (3505->158), mult. (6325->220), div. (0->0), fcn. (5318->12), ass. (0->61)
t56 = V_base(5) * qJ(1) + V_base(1);
t57 = -V_base(4) * qJ(1) + V_base(2);
t60 = sin(pkin(10));
t62 = cos(pkin(10));
t49 = -t56 * t60 + t62 * t57;
t63 = cos(pkin(6));
t52 = t60 * V_base(5) + t62 * V_base(4);
t76 = pkin(7) * t52;
t42 = V_base(6) * pkin(1) - t63 * t76 + t49;
t51 = -t60 * V_base(4) + t62 * V_base(5);
t58 = V_base(3) + qJD(1);
t61 = sin(pkin(6));
t46 = -pkin(1) * t51 - t61 * t76 + t58;
t78 = t42 * t63 + t46 * t61;
t77 = pkin(3) + pkin(8);
t66 = sin(qJ(2));
t69 = cos(qJ(2));
t71 = t51 * t63 + t61 * V_base(6);
t40 = -t66 * t52 + t71 * t69;
t41 = t52 * t69 + t71 * t66;
t59 = sin(pkin(11));
t73 = cos(pkin(11));
t32 = -t73 * t40 + t41 * t59;
t34 = -t42 * t61 + t63 * t46;
t27 = -pkin(2) * t40 + qJD(3) + t34;
t33 = t59 * t40 + t73 * t41;
t70 = -qJ(4) * t33 + t27;
t14 = t77 * t32 + t70;
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t48 = -t51 * t61 + t63 * V_base(6) + qJD(2);
t50 = t62 * t56 + t60 * t57;
t39 = t71 * pkin(7) + t50;
t25 = -t39 * t66 + t78 * t69;
t21 = pkin(2) * t48 - qJ(3) * t41 + t25;
t26 = t69 * t39 + t78 * t66;
t24 = qJ(3) * t40 + t26;
t15 = t73 * t21 - t59 * t24;
t72 = qJD(4) - t15;
t9 = t33 * pkin(4) - t77 * t48 + t72;
t6 = t68 * t14 + t65 * t9;
t16 = t59 * t21 + t73 * t24;
t12 = -t48 * qJ(4) - t16;
t5 = -t14 * t65 + t68 * t9;
t29 = t32 * t68 - t48 * t65;
t10 = -pkin(4) * t32 - t12;
t67 = cos(qJ(6));
t64 = sin(qJ(6));
t31 = qJD(5) + t33;
t30 = t32 * t65 + t48 * t68;
t28 = qJD(6) - t29;
t19 = t30 * t67 + t31 * t64;
t18 = -t30 * t64 + t31 * t67;
t17 = pkin(3) * t32 + t70;
t11 = -t48 * pkin(3) + t72;
t7 = -pkin(5) * t29 - pkin(9) * t30 + t10;
t4 = pkin(9) * t31 + t6;
t3 = -pkin(5) * t31 - t5;
t2 = t4 * t67 + t64 * t7;
t1 = -t4 * t64 + t67 * t7;
t8 = m(2) * (t49 ^ 2 + t50 ^ 2 + t58 ^ 2) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t27 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t17 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t58 * mrSges(2,2) - t49 * mrSges(2,3) + Ifges(2,1) * t52 / 0.2e1) * t52 + (t34 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,1) * t41 / 0.2e1) * t41 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t31 / 0.2e1) * t31 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t58 * mrSges(2,1) + t50 * mrSges(2,3) + Ifges(2,4) * t52 + Ifges(2,2) * t51 / 0.2e1) * t51 + (-t34 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t41 + Ifges(3,2) * t40 / 0.2e1) * t40 + (t10 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t31 + Ifges(6,1) * t30 / 0.2e1) * t30 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t28 + Ifges(7,1) * t19 / 0.2e1) * t19 + (-t10 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t30 + Ifges(6,6) * t31 + Ifges(6,2) * t29 / 0.2e1) * t29 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t28 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t11 * mrSges(5,1) + t27 * mrSges(4,2) - t15 * mrSges(4,3) - t17 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t33) * t33 + (t27 * mrSges(4,1) + t12 * mrSges(5,1) - t17 * mrSges(5,2) - t16 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t32 + (-Ifges(4,4) - Ifges(5,6)) * t33) * t32 + (V_base(2) * mrSges(1,1) + t49 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t50 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t52 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t51 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t25 * mrSges(3,1) + t15 * mrSges(4,1) - t26 * mrSges(3,2) - t16 * mrSges(4,2) + t11 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(3,5) * t41 + Ifges(3,6) * t40 + (Ifges(5,1) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t48 + (-Ifges(5,4) + Ifges(4,5)) * t33 + (Ifges(5,5) - Ifges(4,6)) * t32) * t48;
T  = t8;
