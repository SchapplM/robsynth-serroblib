% Calculate joint inertia matrix for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:42
% EndTime: 2019-03-08 21:22:43
% DurationCPUTime: 0.80s
% Computational Cost: add. (1020->231), mult. (2152->322), div. (0->0), fcn. (2264->10), ass. (0->98)
t122 = Ifges(6,3) + Ifges(7,3);
t121 = 2 * mrSges(7,3);
t90 = m(6) / 0.2e1 + m(7) / 0.2e1;
t120 = 0.2e1 * t90;
t76 = sin(qJ(5));
t79 = cos(qJ(5));
t94 = t76 ^ 2 + t79 ^ 2;
t119 = 0.2e1 * t94;
t74 = cos(pkin(11));
t111 = pkin(3) * t74;
t60 = -pkin(4) - t111;
t46 = -pkin(5) * t79 + t60;
t62 = t76 * mrSges(7,2);
t47 = -t79 * mrSges(7,1) + t62;
t118 = m(7) * t46 + t47;
t73 = sin(pkin(6));
t78 = sin(qJ(2));
t102 = t73 * t78;
t75 = cos(pkin(6));
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t36 = -t77 * t102 + t75 * t80;
t37 = t80 * t102 + t75 * t77;
t72 = sin(pkin(11));
t14 = -t74 * t36 + t37 * t72;
t117 = t14 ^ 2;
t98 = -qJ(4) - pkin(8);
t50 = t98 * t80;
t89 = t98 * t77;
t28 = -t50 * t72 - t74 * t89;
t116 = t28 ^ 2;
t71 = t80 ^ 2;
t115 = 0.2e1 * t28;
t114 = m(5) * pkin(3);
t113 = m(7) * pkin(5);
t112 = pkin(3) * t72;
t110 = mrSges(6,2) * t76;
t109 = Ifges(6,4) * t76;
t108 = Ifges(6,4) * t79;
t107 = Ifges(7,4) * t76;
t106 = Ifges(7,4) * t79;
t105 = t14 * t28;
t44 = t72 * t80 + t74 * t77;
t104 = t44 * t76;
t103 = t44 * t79;
t81 = cos(qJ(2));
t101 = t73 * t81;
t100 = t76 * mrSges(7,3);
t99 = -Ifges(6,6) - Ifges(7,6);
t43 = t72 * t77 - t74 * t80;
t22 = -mrSges(7,2) * t43 - t44 * t100;
t23 = -mrSges(6,2) * t43 - mrSges(6,3) * t104;
t97 = t22 + t23;
t24 = mrSges(7,1) * t43 - mrSges(7,3) * t103;
t25 = mrSges(6,1) * t43 - mrSges(6,3) * t103;
t96 = t24 + t25;
t61 = -pkin(3) * t80 - pkin(2);
t21 = pkin(4) * t43 - pkin(9) * t44 + t61;
t30 = -t74 * t50 + t72 * t89;
t6 = t76 * t21 + t79 * t30;
t19 = mrSges(7,1) * t104 + mrSges(7,2) * t103;
t48 = -mrSges(6,1) * t79 + t110;
t95 = t48 - mrSges(5,1);
t93 = t77 ^ 2 + t71;
t92 = qJ(6) * t44;
t59 = pkin(9) + t112;
t91 = qJ(6) + t59;
t26 = t43 * mrSges(5,1) + t44 * mrSges(5,2);
t5 = t79 * t21 - t30 * t76;
t88 = t122 * t43 + (Ifges(6,5) + Ifges(7,5)) * t103;
t87 = mrSges(6,1) + mrSges(7,1) + t113;
t85 = mrSges(6,1) * t76 + mrSges(6,2) * t79;
t84 = -t36 * t77 + t37 * t80;
t67 = t73 ^ 2;
t66 = Ifges(6,5) * t76;
t65 = Ifges(7,5) * t76;
t64 = Ifges(6,6) * t79;
t63 = Ifges(7,6) * t79;
t57 = t67 * t81 ^ 2;
t54 = Ifges(6,1) * t76 + t108;
t53 = Ifges(7,1) * t76 + t106;
t52 = Ifges(6,2) * t79 + t109;
t51 = Ifges(7,2) * t79 + t107;
t49 = -mrSges(4,1) * t80 + mrSges(4,2) * t77;
t42 = t91 * t79;
t41 = t91 * t76;
t20 = t85 * t44;
t16 = t36 * t72 + t37 * t74;
t13 = pkin(5) * t104 + t28;
t12 = Ifges(6,5) * t43 + (Ifges(6,1) * t79 - t109) * t44;
t11 = Ifges(7,5) * t43 + (Ifges(7,1) * t79 - t107) * t44;
t10 = Ifges(6,6) * t43 + (-Ifges(6,2) * t76 + t108) * t44;
t9 = Ifges(7,6) * t43 + (-Ifges(7,2) * t76 + t106) * t44;
t8 = -t76 * t101 + t16 * t79;
t7 = -t79 * t101 - t16 * t76;
t4 = -t76 * t92 + t6;
t3 = pkin(5) * t43 - t79 * t92 + t5;
t1 = [m(2) + m(5) * (t16 ^ 2 + t117 + t57) + m(4) * (t36 ^ 2 + t37 ^ 2 + t57) + m(3) * (t67 * t78 ^ 2 + t75 ^ 2 + t57) + (t7 ^ 2 + t8 ^ 2 + t117) * t120; -t16 * t43 * mrSges(5,3) + t97 * t8 + t96 * t7 + t84 * mrSges(4,3) + (t44 * mrSges(5,3) + t19 + t20) * t14 + (-t78 * mrSges(3,2) + (mrSges(3,1) - t26 - t49) * t81) * t73 + m(7) * (t13 * t14 + t3 * t7 + t4 * t8) + m(6) * (t5 * t7 + t6 * t8 + t105) + m(5) * (-t61 * t101 + t16 * t30 + t105) + m(4) * (pkin(2) * t101 + t84 * pkin(8)); Ifges(4,2) * t71 - 0.2e1 * pkin(2) * t49 + 0.2e1 * t13 * t19 + t20 * t115 + 0.2e1 * t4 * t22 + 0.2e1 * t6 * t23 + 0.2e1 * t3 * t24 + 0.2e1 * t5 * t25 + 0.2e1 * t61 * t26 + Ifges(3,3) + (Ifges(4,1) * t77 + 0.2e1 * Ifges(4,4) * t80) * t77 + 0.2e1 * t93 * pkin(8) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t30 + Ifges(5,2) * t43 + t88) * t43 + m(4) * (t93 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t30 ^ 2 + t61 ^ 2 + t116) + m(7) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2 + t116) + (mrSges(5,3) * t115 + Ifges(5,1) * t44 - 0.2e1 * Ifges(5,4) * t43 + (t11 + t12) * t79 + (t99 * t43 - t10 - t9) * t76) * t44; t36 * mrSges(4,1) - t37 * mrSges(4,2) + m(7) * (-t41 * t7 + t42 * t8) + (t72 * t114 - mrSges(5,2)) * t16 + (m(6) * t60 - t74 * t114 + t118 + t95) * t14 + (m(6) * t59 + mrSges(6,3) + mrSges(7,3)) * (-t7 * t76 + t79 * t8); -t30 * mrSges(5,2) + Ifges(4,5) * t77 + Ifges(4,6) * t80 + t13 * t47 + t46 * t19 + t60 * t20 + t42 * t22 - t41 * t24 + (t65 / 0.2e1 + t63 / 0.2e1 + t66 / 0.2e1 + t64 / 0.2e1 - Ifges(5,6) - mrSges(5,3) * t112) * t43 + t95 * t28 + (-mrSges(4,1) * t77 - mrSges(4,2) * t80) * pkin(8) + (t59 * t23 + t4 * mrSges(7,3) + t6 * mrSges(6,3) + t9 / 0.2e1 + t10 / 0.2e1) * t79 + (-t59 * t25 - t3 * mrSges(7,3) - t5 * mrSges(6,3) + t11 / 0.2e1 + t12 / 0.2e1) * t76 + m(6) * (t28 * t60 + (-t5 * t76 + t6 * t79) * t59) + m(7) * (t13 * t46 - t3 * t41 + t4 * t42) + (-t28 * t74 + t30 * t72) * t114 + (-mrSges(5,3) * t111 + Ifges(5,5) + (t53 / 0.2e1 + t54 / 0.2e1) * t79 + (-t51 / 0.2e1 - t52 / 0.2e1) * t76) * t44; 0.2e1 * t46 * t47 + 0.2e1 * t60 * t48 + Ifges(4,3) + Ifges(5,3) + (t42 * t121 + t51 + t52) * t79 + (t41 * t121 + t53 + t54) * t76 + m(7) * (t41 ^ 2 + t42 ^ 2 + t46 ^ 2) + m(6) * (t94 * t59 ^ 2 + t60 ^ 2) + m(5) * (t72 ^ 2 + t74 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t74 - mrSges(5,2) * t72) * pkin(3) + t59 * mrSges(6,3) * t119; -m(5) * t101 + (t7 * t79 + t76 * t8) * t120; t96 * t79 + t97 * t76 + m(7) * (t3 * t79 + t4 * t76) + m(6) * (t5 * t79 + t6 * t76) + m(5) * t61 + t26; m(7) * (-t41 * t79 + t42 * t76); t90 * t119 + m(5); (-mrSges(6,2) - mrSges(7,2)) * t8 + t87 * t7; mrSges(6,1) * t5 + mrSges(7,1) * t3 - mrSges(6,2) * t6 - mrSges(7,2) * t4 + t99 * t104 + (m(7) * t3 + t24) * pkin(5) + t88; -mrSges(7,1) * t41 - mrSges(7,2) * t42 + t63 + t64 + t65 + t66 - t85 * t59 + (-m(7) * t41 - t100) * pkin(5); t87 * t79 - t110 - t62; (0.2e1 * mrSges(7,1) + t113) * pkin(5) + t122; m(7) * t14; m(7) * t13 + t19; t118; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
