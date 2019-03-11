% Calculate kinetic energy for
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:47:57
% EndTime: 2019-03-09 05:47:59
% DurationCPUTime: 1.48s
% Computational Cost: add. (6777->163), mult. (11045->233), div. (0->0), fcn. (9544->14), ass. (0->69)
t60 = V_base(5) * pkin(8) + V_base(1);
t61 = -V_base(4) * pkin(8) + V_base(2);
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t53 = -t60 * t72 + t75 * t61;
t62 = V_base(6) + qJD(1);
t68 = cos(pkin(6));
t56 = t72 * V_base(5) + t75 * V_base(4);
t82 = qJ(2) * t56;
t47 = pkin(1) * t62 - t68 * t82 + t53;
t55 = -t72 * V_base(4) + t75 * V_base(5);
t65 = sin(pkin(6));
t51 = -pkin(1) * t55 - t65 * t82 + V_base(3);
t92 = t47 * t68 + t51 * t65;
t54 = t75 * t60 + t72 * t61;
t79 = t55 * t68 + t62 * t65;
t44 = qJ(2) * t79 + t54;
t63 = sin(pkin(12));
t66 = cos(pkin(12));
t33 = -t44 * t63 + t66 * t92;
t52 = -t55 * t65 + t62 * t68;
t67 = cos(pkin(7));
t46 = t56 * t66 + t63 * t79;
t89 = pkin(9) * t46;
t29 = pkin(2) * t52 - t67 * t89 + t33;
t38 = -t47 * t65 + t68 * t51 + qJD(2);
t45 = -t56 * t63 + t66 * t79;
t64 = sin(pkin(7));
t32 = -pkin(2) * t45 - t64 * t89 + t38;
t91 = t29 * t67 + t32 * t64;
t34 = t66 * t44 + t63 * t92;
t80 = t45 * t67 + t52 * t64;
t26 = pkin(9) * t80 + t34;
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t17 = -t71 * t26 + t74 * t91;
t36 = -t46 * t71 + t74 * t80;
t90 = pkin(4) + pkin(11);
t88 = cos(qJ(4));
t18 = t74 * t26 + t71 * t91;
t40 = -t45 * t64 + t52 * t67 + qJD(3);
t14 = pkin(10) * t40 + t18;
t19 = -t29 * t64 + t67 * t32;
t37 = t46 * t74 + t71 * t80;
t16 = -pkin(3) * t36 - pkin(10) * t37 + t19;
t70 = sin(qJ(4));
t8 = t88 * t14 + t70 * t16;
t35 = qJD(4) - t36;
t6 = -qJ(5) * t35 - t8;
t7 = -t70 * t14 + t16 * t88;
t78 = qJD(5) - t7;
t28 = t37 * t88 + t70 * t40;
t13 = -pkin(3) * t40 - t17;
t77 = -qJ(5) * t28 + t13;
t76 = V_base(3) ^ 2;
t73 = cos(qJ(6));
t69 = sin(qJ(6));
t27 = t37 * t70 - t40 * t88;
t25 = qJD(6) + t28;
t21 = t27 * t69 + t35 * t73;
t20 = t27 * t73 - t35 * t69;
t10 = pkin(4) * t27 + t77;
t9 = t27 * t90 + t77;
t5 = -t35 * pkin(4) + t78;
t4 = -pkin(5) * t27 - t6;
t3 = t28 * pkin(5) - t35 * t90 + t78;
t2 = t3 * t69 + t73 * t9;
t1 = t3 * t73 - t69 * t9;
t11 = (t38 * mrSges(3,2) - t33 * mrSges(3,3) + Ifges(3,5) * t52 + Ifges(3,1) * t46 / 0.2e1) * t46 + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t25 + Ifges(7,2) * t20 / 0.2e1) * t20 + (-t19 * mrSges(4,1) + t18 * mrSges(4,3) + t37 * Ifges(4,4) + Ifges(4,6) * t40 + Ifges(4,2) * t36 / 0.2e1) * t36 + (t53 * mrSges(2,1) - t54 * mrSges(2,2) + Ifges(2,3) * t62 / 0.2e1) * t62 + (t19 * mrSges(4,2) - t17 * mrSges(4,3) + Ifges(4,5) * t40 + Ifges(4,1) * t37 / 0.2e1) * t37 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(2,1) + t54 * mrSges(2,3) + Ifges(2,4) * t56 + Ifges(2,6) * t62 + Ifges(2,2) * t55 / 0.2e1) * t55 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t21 / 0.2e1) * t21 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t13 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t27 + (Ifges(6,5) - Ifges(5,6)) * t35 + (-Ifges(5,4) - Ifges(6,6)) * t28) * t27 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(2,2) - t53 * mrSges(2,3) + Ifges(2,5) * t62 + Ifges(2,1) * t56 / 0.2e1) * t56 + (t33 * mrSges(3,1) - t34 * mrSges(3,2) + Ifges(3,3) * t52 / 0.2e1) * t52 + (t5 * mrSges(6,1) + t13 * mrSges(5,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t28 + (-Ifges(6,4) + Ifges(5,5)) * t35) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t38 * mrSges(3,1) + t34 * mrSges(3,3) + Ifges(3,4) * t46 + Ifges(3,6) * t52 + Ifges(3,2) * t45 / 0.2e1) * t45 + (t7 * mrSges(5,1) - t8 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t35) * t35 + (t17 * mrSges(4,1) - t18 * mrSges(4,2) + Ifges(4,3) * t40 / 0.2e1) * t40 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t38 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t76) / 0.2e1 + m(2) * (t53 ^ 2 + t54 ^ 2 + t76) / 0.2e1;
T  = t11;
