% Calculate joint inertia matrix for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:00
% EndTime: 2019-03-09 01:44:01
% DurationCPUTime: 0.48s
% Computational Cost: add. (723->154), mult. (1253->222), div. (0->0), fcn. (1177->8), ass. (0->64)
t51 = cos(qJ(4));
t86 = t51 ^ 2;
t44 = sin(pkin(10));
t46 = cos(pkin(10));
t49 = sin(qJ(4));
t20 = t44 * t51 + t46 * t49;
t19 = t44 * t49 - t46 * t51;
t48 = sin(qJ(6));
t73 = t19 * t48;
t10 = -t20 * mrSges(7,2) + mrSges(7,3) * t73;
t50 = cos(qJ(6));
t72 = t19 * t50;
t11 = t20 * mrSges(7,1) + mrSges(7,3) * t72;
t55 = -t50 * t10 + t48 * t11;
t17 = t20 ^ 2;
t18 = t19 ^ 2;
t65 = t49 ^ 2 + t86;
t23 = m(5) * t65;
t85 = m(4) + t23 + m(6) * (t17 + t18);
t47 = cos(pkin(9));
t36 = -t47 * pkin(1) - pkin(2);
t29 = -pkin(7) + t36;
t64 = -qJ(5) + t29;
t14 = t64 * t49;
t60 = t64 * t51;
t6 = t44 * t14 - t46 * t60;
t84 = t6 ^ 2;
t45 = sin(pkin(9));
t31 = t45 * pkin(1) + qJ(3);
t83 = t31 ^ 2;
t82 = -2 * mrSges(6,3);
t80 = m(6) * pkin(4);
t79 = t48 / 0.2e1;
t78 = t19 * t6;
t77 = t6 * t20;
t75 = Ifges(7,4) * t48;
t74 = Ifges(7,4) * t50;
t25 = -t50 * mrSges(7,1) + t48 * mrSges(7,2);
t69 = -mrSges(6,1) + t25;
t68 = -Ifges(7,5) * t72 + Ifges(7,3) * t20;
t67 = Ifges(7,5) * t48 + Ifges(7,6) * t50;
t66 = t48 ^ 2 + t50 ^ 2;
t63 = t66 * t19;
t34 = t44 * pkin(4) + pkin(8);
t62 = t66 * t34;
t61 = t65 * mrSges(5,3);
t24 = t49 * pkin(4) + t31;
t5 = t20 * pkin(5) + t19 * pkin(8) + t24;
t8 = t46 * t14 + t44 * t60;
t1 = -t48 * t8 + t50 * t5;
t2 = t48 * t5 + t50 * t8;
t59 = t1 * t48 - t2 * t50;
t58 = t51 * mrSges(5,1) - t49 * mrSges(5,2);
t57 = t49 * mrSges(5,1) + t51 * mrSges(5,2);
t56 = -mrSges(7,1) * t48 - mrSges(7,2) * t50;
t54 = t19 * t46 - t20 * t44;
t35 = -t46 * pkin(4) - pkin(5);
t27 = Ifges(7,1) * t48 + t74;
t26 = Ifges(7,2) * t50 + t75;
t15 = t19 * mrSges(6,2);
t9 = t56 * t19;
t4 = Ifges(7,5) * t20 + (-Ifges(7,1) * t50 + t75) * t19;
t3 = Ifges(7,6) * t20 + (Ifges(7,2) * t48 - t74) * t19;
t7 = [Ifges(5,1) * t86 + 0.2e1 * t36 * mrSges(4,2) - 0.2e1 * t24 * t15 + 0.2e1 * t6 * t9 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t24 * mrSges(6,1) + Ifges(6,2) * t20 + t8 * t82 + t68) * t20 + (t6 * t82 + Ifges(6,1) * t19 + t48 * t3 - t50 * t4 + (Ifges(7,6) * t48 + (2 * Ifges(6,4))) * t20) * t19 + m(5) * (t65 * t29 ^ 2 + t83) + m(6) * (t24 ^ 2 + t8 ^ 2 + t84) + m(4) * (t36 ^ 2 + t83) + m(7) * (t1 ^ 2 + t2 ^ 2 + t84) + m(3) * (t45 ^ 2 + t47 ^ 2) * pkin(1) ^ 2 + (-0.2e1 * Ifges(5,4) * t51 + Ifges(5,2) * t49) * t49 + 0.2e1 * (t57 + mrSges(4,3)) * t31 + 0.2e1 * (t47 * mrSges(3,1) - t45 * mrSges(3,2)) * pkin(1) - 0.2e1 * t29 * t61; t20 * t9 + t55 * t19 + m(6) * (-t8 * t19 + t77) + m(7) * (t19 * t59 + t77); m(3) + m(7) * (t18 * t66 + t17) + t85; -t18 * mrSges(6,3) + t19 * t9 + mrSges(4,2) - t61 + (-mrSges(6,3) * t20 - t55) * t20 + m(7) * (-t20 * t59 + t78) + m(6) * (t20 * t8 + t78) + m(4) * t36 + t29 * t23; m(7) * (0.1e1 - t66) * t20 * t19; m(7) * (t17 * t66 + t18) + t85; Ifges(5,5) * t51 - Ifges(5,6) * t49 + t4 * t79 + t50 * t3 / 0.2e1 + t6 * t25 - t6 * mrSges(6,1) - t8 * mrSges(6,2) + t58 * t29 - t59 * mrSges(7,3) + (-Ifges(6,5) - t50 * t27 / 0.2e1 + t26 * t79) * t19 + (m(6) * (t44 * t8 - t46 * t6) + t54 * mrSges(6,3)) * pkin(4) + (m(7) * t6 + t9) * t35 + (-m(7) * t59 - t55) * t34 + (-Ifges(6,6) + t67 / 0.2e1) * t20; t15 + t69 * t20 - mrSges(7,3) * t63 + m(7) * (t35 * t20 - t34 * t63) + (-t19 * t44 - t20 * t46) * t80 - t57; t69 * t19 + (mrSges(7,3) * t66 - mrSges(6,2)) * t20 + m(7) * (t35 * t19 + t20 * t62) - t54 * t80 + t58; 0.2e1 * t35 * t25 + t50 * t26 + t48 * t27 + Ifges(5,3) + Ifges(6,3) + m(7) * (t34 ^ 2 * t66 + t35 ^ 2) + m(6) * (t44 ^ 2 + t46 ^ 2) * pkin(4) ^ 2 + 0.2e1 * (t46 * mrSges(6,1) - t44 * mrSges(6,2)) * pkin(4) + 0.2e1 * mrSges(7,3) * t62; t20 * mrSges(6,1) + t48 * t10 + t50 * t11 - t15 + m(7) * (t50 * t1 + t48 * t2) + m(6) * t24; 0; 0; 0; m(7) * t66 + m(6); t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,6) * t73 + t68; -t9; t56 * t20; t34 * t56 + t67; -t25; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
