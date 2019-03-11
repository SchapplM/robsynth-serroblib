% Calculate kinetic energy for
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:37
% EndTime: 2019-03-08 18:59:39
% DurationCPUTime: 1.79s
% Computational Cost: add. (8669->166), mult. (16419->249), div. (0->0), fcn. (14354->16), ass. (0->72)
t67 = qJ(1) * V_base(5) + V_base(1);
t68 = -qJ(1) * V_base(4) + V_base(2);
t71 = sin(pkin(12));
t75 = cos(pkin(12));
t60 = -t67 * t71 + t75 * t68;
t77 = cos(pkin(6));
t63 = t71 * V_base(5) + t75 * V_base(4);
t89 = qJ(2) * t63;
t54 = V_base(6) * pkin(1) - t77 * t89 + t60;
t62 = -t71 * V_base(4) + t75 * V_base(5);
t69 = V_base(3) + qJD(1);
t73 = sin(pkin(6));
t58 = -pkin(1) * t62 - t73 * t89 + t69;
t96 = t54 * t77 + t58 * t73;
t61 = t75 * t67 + t71 * t68;
t86 = t62 * t77 + t73 * V_base(6);
t51 = qJ(2) * t86 + t61;
t70 = sin(pkin(13));
t74 = cos(pkin(13));
t39 = -t51 * t70 + t96 * t74;
t59 = -t62 * t73 + t77 * V_base(6);
t76 = cos(pkin(7));
t53 = t63 * t74 + t70 * t86;
t94 = pkin(8) * t53;
t33 = pkin(2) * t59 - t76 * t94 + t39;
t45 = -t54 * t73 + t77 * t58 + qJD(2);
t52 = -t63 * t70 + t74 * t86;
t72 = sin(pkin(7));
t38 = -pkin(2) * t52 - t72 * t94 + t45;
t95 = t33 * t76 + t38 * t72;
t40 = t74 * t51 + t96 * t70;
t87 = t52 * t76 + t59 * t72;
t32 = pkin(8) * t87 + t40;
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t21 = -t81 * t32 + t85 * t95;
t43 = -t81 * t53 + t85 * t87;
t22 = t85 * t32 + t95 * t81;
t46 = -t52 * t72 + t59 * t76 + qJD(3);
t17 = pkin(9) * t46 + t22;
t28 = -t33 * t72 + t76 * t38;
t44 = t53 * t85 + t81 * t87;
t20 = -pkin(3) * t43 - pkin(9) * t44 + t28;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t13 = t84 * t17 + t80 * t20;
t34 = -t44 * t80 + t46 * t84;
t11 = pkin(10) * t34 + t13;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t12 = -t17 * t80 + t84 * t20;
t35 = t44 * t84 + t46 * t80;
t42 = qJD(4) - t43;
t9 = pkin(4) * t42 - pkin(10) * t35 + t12;
t6 = t83 * t11 + t79 * t9;
t5 = -t11 * t79 + t83 * t9;
t26 = t34 * t83 - t35 * t79;
t16 = -t46 * pkin(3) - t21;
t14 = -t34 * pkin(4) + t16;
t82 = cos(qJ(6));
t78 = sin(qJ(6));
t41 = qJD(5) + t42;
t27 = t34 * t79 + t35 * t83;
t25 = qJD(6) - t26;
t24 = t27 * t82 + t41 * t78;
t23 = -t27 * t78 + t41 * t82;
t7 = -t26 * pkin(5) - t27 * pkin(11) + t14;
t4 = pkin(11) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t4 * t82 + t7 * t78;
t1 = -t4 * t78 + t7 * t82;
t8 = (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t28 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,6) * t46 + Ifges(4,2) * t43 / 0.2e1) * t43 + (t28 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t46 + Ifges(4,1) * t44 / 0.2e1) * t44 + (t45 * mrSges(3,2) - t39 * mrSges(3,3) + Ifges(3,5) * t59 + Ifges(3,1) * t53 / 0.2e1) * t53 + (-t69 * mrSges(2,1) + t61 * mrSges(2,3) + Ifges(2,4) * t63 + Ifges(2,2) * t62 / 0.2e1) * t62 + (-t16 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t35 + Ifges(5,6) * t42 + Ifges(5,2) * t34 / 0.2e1) * t34 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t41 / 0.2e1) * t41 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t41 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t69 * mrSges(2,2) - t60 * mrSges(2,3) + Ifges(2,1) * t63 / 0.2e1) * t63 + (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t42 / 0.2e1) * t42 + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t41 + Ifges(6,2) * t26 / 0.2e1) * t26 + m(2) * (t60 ^ 2 + t61 ^ 2 + t69 ^ 2) / 0.2e1 + m(3) * (t39 ^ 2 + t40 ^ 2 + t45 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t45 * mrSges(3,1) + t40 * mrSges(3,3) + Ifges(3,4) * t53 + Ifges(3,6) * t59 + Ifges(3,2) * t52 / 0.2e1) * t52 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t24 + Ifges(7,6) * t25 + Ifges(7,2) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t16 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t42 + Ifges(5,1) * t35 / 0.2e1) * t35 + (t39 * mrSges(3,1) - t40 * mrSges(3,2) + Ifges(3,3) * t59 / 0.2e1) * t59 + (V_base(2) * mrSges(1,1) + t60 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t61 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t63 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t62 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t24 / 0.2e1) * t24 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t46 / 0.2e1) * t46 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1;
T  = t8;
