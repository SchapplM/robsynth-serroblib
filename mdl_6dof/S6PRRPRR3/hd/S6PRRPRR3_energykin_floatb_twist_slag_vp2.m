% Calculate kinetic energy for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:16
% EndTime: 2019-03-08 22:02:18
% DurationCPUTime: 1.62s
% Computational Cost: add. (8943->166), mult. (16455->247), div. (0->0), fcn. (14382->16), ass. (0->71)
t66 = qJ(1) * V_base(5) + V_base(1);
t67 = -qJ(1) * V_base(4) + V_base(2);
t70 = sin(pkin(12));
t74 = cos(pkin(12));
t59 = -t66 * t70 + t74 * t67;
t76 = cos(pkin(6));
t62 = t70 * V_base(5) + t74 * V_base(4);
t92 = pkin(8) * t62;
t53 = V_base(6) * pkin(1) - t76 * t92 + t59;
t61 = -t70 * V_base(4) + t74 * V_base(5);
t68 = V_base(3) + qJD(1);
t72 = sin(pkin(6));
t57 = -pkin(1) * t61 - t72 * t92 + t68;
t94 = t53 * t76 + t57 * t72;
t60 = t74 * t66 + t70 * t67;
t85 = t61 * t76 + t72 * V_base(6);
t50 = t85 * pkin(8) + t60;
t80 = sin(qJ(2));
t84 = cos(qJ(2));
t41 = -t50 * t80 + t94 * t84;
t58 = -t61 * t72 + t76 * V_base(6) + qJD(2);
t75 = cos(pkin(7));
t52 = t62 * t84 + t85 * t80;
t91 = pkin(9) * t52;
t36 = pkin(2) * t58 - t75 * t91 + t41;
t45 = -t53 * t72 + t76 * t57;
t51 = -t80 * t62 + t85 * t84;
t71 = sin(pkin(7));
t40 = -pkin(2) * t51 - t71 * t91 + t45;
t93 = t36 * t75 + t40 * t71;
t42 = t84 * t50 + t94 * t80;
t86 = t51 * t75 + t58 * t71;
t35 = t86 * pkin(9) + t42;
t79 = sin(qJ(3));
t83 = cos(qJ(3));
t20 = -t35 * t79 + t93 * t83;
t44 = t52 * t83 + t86 * t79;
t46 = -t51 * t71 + t58 * t75 + qJD(3);
t16 = pkin(3) * t46 - qJ(4) * t44 + t20;
t21 = t83 * t35 + t93 * t79;
t43 = -t52 * t79 + t86 * t83;
t19 = qJ(4) * t43 + t21;
t69 = sin(pkin(13));
t73 = cos(pkin(13));
t12 = t69 * t16 + t73 * t19;
t10 = pkin(10) * t46 + t12;
t25 = -t36 * t71 + t75 * t40;
t24 = -pkin(3) * t43 + qJD(4) + t25;
t30 = t43 * t73 - t44 * t69;
t31 = t43 * t69 + t44 * t73;
t14 = -pkin(4) * t30 - pkin(10) * t31 + t24;
t78 = sin(qJ(5));
t82 = cos(qJ(5));
t6 = t82 * t10 + t78 * t14;
t11 = t16 * t73 - t69 * t19;
t5 = -t10 * t78 + t14 * t82;
t27 = -t31 * t78 + t46 * t82;
t9 = -pkin(4) * t46 - t11;
t81 = cos(qJ(6));
t77 = sin(qJ(6));
t29 = qJD(5) - t30;
t28 = t31 * t82 + t46 * t78;
t26 = qJD(6) - t27;
t23 = t28 * t81 + t29 * t77;
t22 = -t28 * t77 + t29 * t81;
t7 = -pkin(5) * t27 - pkin(11) * t28 + t9;
t4 = pkin(11) * t29 + t6;
t3 = -pkin(5) * t29 - t5;
t2 = t4 * t81 + t7 * t77;
t1 = -t4 * t77 + t7 * t81;
t8 = (t24 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,1) * t31 / 0.2e1) * t31 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (-t68 * mrSges(2,1) + t60 * mrSges(2,3) + Ifges(2,4) * t62 + Ifges(2,2) * t61 / 0.2e1) * t61 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t29 + Ifges(6,2) * t27 / 0.2e1) * t27 + (-t24 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t31 + Ifges(5,2) * t30 / 0.2e1) * t30 + m(2) * (t59 ^ 2 + t60 ^ 2 + t68 ^ 2) / 0.2e1 + m(3) * (t41 ^ 2 + t42 ^ 2 + t45 ^ 2) / 0.2e1 + m(4) * (t20 ^ 2 + t21 ^ 2 + t25 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t24 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t23 + Ifges(7,6) * t26 + Ifges(7,2) * t22 / 0.2e1) * t22 + (-t45 * mrSges(3,1) + t42 * mrSges(3,3) + Ifges(3,4) * t52 + Ifges(3,6) * t58 + Ifges(3,2) * t51 / 0.2e1) * t51 + (t68 * mrSges(2,2) - t59 * mrSges(2,3) + Ifges(2,1) * t62 / 0.2e1) * t62 + (t45 * mrSges(3,2) - t41 * mrSges(3,3) + Ifges(3,5) * t58 + Ifges(3,1) * t52 / 0.2e1) * t52 + (V_base(2) * mrSges(1,1) + t59 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t60 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t62 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t61 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t26 / 0.2e1) * t26 + (-t25 * mrSges(4,1) + t21 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,2) * t43 / 0.2e1) * t43 + (t20 * mrSges(4,1) + t11 * mrSges(5,1) - t21 * mrSges(4,2) - t12 * mrSges(5,2) + Ifges(4,5) * t44 + Ifges(5,5) * t31 + Ifges(4,6) * t43 + Ifges(5,6) * t30 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t46) * t46 + (t25 * mrSges(4,2) - t20 * mrSges(4,3) + Ifges(4,1) * t44 / 0.2e1) * t44 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t41 * mrSges(3,1) - t42 * mrSges(3,2) + Ifges(3,3) * t58 / 0.2e1) * t58 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t26 + Ifges(7,1) * t23 / 0.2e1) * t23;
T  = t8;
