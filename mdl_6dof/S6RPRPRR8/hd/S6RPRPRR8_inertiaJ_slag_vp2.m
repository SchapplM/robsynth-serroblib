% Calculate joint inertia matrix for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:47
% EndTime: 2019-03-09 03:57:49
% DurationCPUTime: 0.82s
% Computational Cost: add. (1511->238), mult. (2690->349), div. (0->0), fcn. (2865->8), ass. (0->86)
t82 = cos(qJ(3));
t119 = t82 ^ 2;
t75 = sin(pkin(10));
t76 = cos(pkin(10));
t79 = sin(qJ(3));
t52 = t75 * t79 - t76 * t82;
t81 = cos(qJ(5));
t106 = t52 * t81;
t53 = t75 * t82 + t76 * t79;
t118 = -Ifges(6,5) * t106 + Ifges(6,3) * t53;
t83 = -pkin(1) - pkin(7);
t100 = -qJ(4) + t83;
t57 = t100 * t79;
t94 = t100 * t82;
t33 = t57 * t75 - t76 * t94;
t117 = t33 ^ 2;
t116 = t52 ^ 2;
t115 = -2 * mrSges(5,3);
t113 = m(7) * pkin(5);
t64 = pkin(3) * t75 + pkin(8);
t111 = pkin(9) + t64;
t78 = sin(qJ(5));
t110 = Ifges(6,4) * t78;
t109 = Ifges(6,4) * t81;
t108 = t33 * t52;
t107 = t52 * t78;
t66 = t79 * pkin(3) + qJ(2);
t29 = pkin(4) * t53 + pkin(8) * t52 + t66;
t35 = t76 * t57 + t75 * t94;
t10 = t78 * t29 + t81 * t35;
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t56 = t77 * t81 + t78 * t80;
t90 = t77 * t78 - t80 * t81;
t105 = Ifges(7,5) * t56 - Ifges(7,6) * t90;
t59 = -t81 * mrSges(6,1) + t78 * mrSges(6,2);
t104 = t59 - mrSges(5,1);
t103 = Ifges(6,5) * t78 + Ifges(6,6) * t81;
t102 = t78 ^ 2 + t81 ^ 2;
t101 = t79 ^ 2 + t119;
t20 = t56 * t52;
t22 = t90 * t52;
t99 = Ifges(7,5) * t22 + Ifges(7,6) * t20 + Ifges(7,3) * t53;
t65 = -pkin(3) * t76 - pkin(4);
t98 = m(4) * t101;
t97 = t102 * t64;
t96 = t101 * mrSges(4,3);
t19 = t56 * t53;
t21 = t90 * t53;
t95 = -t19 * mrSges(7,1) + t21 * mrSges(7,2);
t36 = mrSges(7,1) * t90 + t56 * mrSges(7,2);
t9 = t81 * t29 - t35 * t78;
t93 = t10 * t81 - t78 * t9;
t92 = -mrSges(6,1) * t78 - mrSges(6,2) * t81;
t91 = t76 * t52 - t75 * t53;
t44 = t111 * t78;
t45 = t111 * t81;
t27 = -t44 * t80 - t45 * t77;
t28 = -t44 * t77 + t45 * t80;
t89 = t27 * mrSges(7,1) - t28 * mrSges(7,2) + t105;
t4 = pkin(5) * t53 + pkin(9) * t106 + t9;
t7 = pkin(9) * t107 + t10;
t2 = t4 * t80 - t7 * t77;
t3 = t4 * t77 + t7 * t80;
t88 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t99;
t87 = (mrSges(7,1) * t80 - mrSges(7,2) * t77) * pkin(5);
t84 = qJ(2) ^ 2;
t61 = Ifges(6,1) * t78 + t109;
t60 = Ifges(6,2) * t81 + t110;
t58 = -pkin(5) * t81 + t65;
t50 = t53 ^ 2;
t40 = t52 * mrSges(5,2);
t38 = Ifges(7,1) * t56 - Ifges(7,4) * t90;
t37 = Ifges(7,4) * t56 - Ifges(7,2) * t90;
t31 = mrSges(6,1) * t53 + mrSges(6,3) * t106;
t30 = -mrSges(6,2) * t53 + mrSges(6,3) * t107;
t26 = t92 * t52;
t15 = -pkin(5) * t107 + t33;
t14 = Ifges(6,5) * t53 + (-Ifges(6,1) * t81 + t110) * t52;
t13 = Ifges(6,6) * t53 + (Ifges(6,2) * t78 - t109) * t52;
t12 = mrSges(7,1) * t53 - mrSges(7,3) * t22;
t11 = -mrSges(7,2) * t53 + mrSges(7,3) * t20;
t8 = -mrSges(7,1) * t20 + mrSges(7,2) * t22;
t6 = Ifges(7,1) * t22 + Ifges(7,4) * t20 + Ifges(7,5) * t53;
t5 = Ifges(7,4) * t22 + Ifges(7,2) * t20 + Ifges(7,6) * t53;
t1 = [Ifges(4,1) * t119 - 0.2e1 * t66 * t40 + 0.2e1 * t33 * t26 + t22 * t6 + 0.2e1 * t10 * t30 + 0.2e1 * t9 * t31 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t15 * t8 + t20 * t5 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t83 * t96 + (0.2e1 * mrSges(5,1) * t66 + Ifges(5,2) * t53 + t35 * t115 + t118 + t99) * t53 + (t33 * t115 + Ifges(5,1) * t52 + t78 * t13 - t81 * t14 + (Ifges(6,6) * t78 + (2 * Ifges(5,4))) * t53) * t52 + m(4) * (t101 * t83 ^ 2 + t84) + m(3) * ((pkin(1) ^ 2) + t84) + m(5) * (t35 ^ 2 + t66 ^ 2 + t117) + m(6) * (t10 ^ 2 + t9 ^ 2 + t117) + m(7) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) + (-0.2e1 * Ifges(4,4) * t82 + Ifges(4,2) * t79) * t79 + 0.2e1 * (mrSges(4,1) * t79 + mrSges(4,2) * t82 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t116 * mrSges(5,3) - t21 * t11 - t19 * t12 + mrSges(3,2) + (t26 + t8) * t52 - t96 + (-mrSges(5,3) * t53 + t81 * t30 - t78 * t31) * t53 + m(7) * (t15 * t52 - t19 * t2 - t21 * t3) + m(6) * (t53 * t93 + t108) + m(5) * (t35 * t53 + t108) + t83 * t98; m(3) + m(7) * (t19 ^ 2 + t21 ^ 2 + t116) + m(6) * (t102 * t50 + t116) + m(5) * (t50 + t116) + t98; t65 * t26 - Ifges(5,6) * t53 - t90 * t5 / 0.2e1 + t56 * t6 / 0.2e1 + t58 * t8 - t35 * mrSges(5,2) + t15 * t36 + t20 * t37 / 0.2e1 + t22 * t38 / 0.2e1 + t27 * t12 + t28 * t11 + (t83 * mrSges(4,1) + Ifges(4,5)) * t82 + (-t83 * mrSges(4,2) - Ifges(4,6)) * t79 + t104 * t33 + (-t2 * t56 - t3 * t90) * mrSges(7,3) + (t64 * t30 + t10 * mrSges(6,3) + t13 / 0.2e1) * t81 + (-t64 * t31 - t9 * mrSges(6,3) + t14 / 0.2e1) * t78 + m(6) * (t33 * t65 + t64 * t93) + m(7) * (t15 * t58 + t2 * t27 + t28 * t3) + (-t81 * t61 / 0.2e1 + t78 * t60 / 0.2e1 - Ifges(5,5)) * t52 + (m(5) * (-t33 * t76 + t35 * t75) + t91 * mrSges(5,3)) * pkin(3) + (t103 + t105) * t53 / 0.2e1; t82 * mrSges(4,1) - t79 * mrSges(4,2) + (t19 * t56 + t21 * t90) * mrSges(7,3) + (mrSges(6,3) * t102 - mrSges(5,2)) * t53 + (t36 + t104) * t52 + m(7) * (-t19 * t27 - t21 * t28 + t52 * t58) + m(6) * (t52 * t65 + t53 * t97) - m(5) * t91 * pkin(3); 0.2e1 * t58 * t36 - t90 * t37 + t56 * t38 + 0.2e1 * t65 * t59 + t81 * t60 + t78 * t61 + Ifges(4,3) + Ifges(5,3) + m(7) * (t27 ^ 2 + t28 ^ 2 + t58 ^ 2) + m(6) * (t102 * t64 ^ 2 + t65 ^ 2) + m(5) * (t75 ^ 2 + t76 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t76 - mrSges(5,2) * t75) * pkin(3) + 0.2e1 * (-t27 * t56 - t28 * t90) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t97; t53 * mrSges(5,1) + t56 * t11 - t90 * t12 + t78 * t30 + t81 * t31 - t40 + m(7) * (-t2 * t90 + t3 * t56) + m(6) * (t10 * t78 + t81 * t9) + m(5) * t66; m(7) * (t19 * t90 - t21 * t56); m(7) * (-t27 * t90 + t28 * t56); m(5) + m(6) * t102 + m(7) * (t56 ^ 2 + t90 ^ 2); Ifges(6,6) * t107 + t9 * mrSges(6,1) - t10 * mrSges(6,2) + (t77 * t11 + t80 * t12 + m(7) * (t2 * t80 + t3 * t77)) * pkin(5) + t88 + t118; t92 * t53 + (-t19 * t80 - t21 * t77) * t113 + t95; t92 * t64 + (m(7) * (t27 * t80 + t28 * t77) + (-t56 * t80 - t77 * t90) * mrSges(7,3)) * pkin(5) + t89 + t103; (t56 * t77 - t80 * t90) * t113 - t59 - t36; Ifges(6,3) + Ifges(7,3) + m(7) * (t77 ^ 2 + t80 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t87; t88; t95; t89; -t36; Ifges(7,3) + t87; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
