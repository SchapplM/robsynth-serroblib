% Calculate joint inertia matrix for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14V3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:28
% EndTime: 2019-04-12 15:03:30
% DurationCPUTime: 0.89s
% Computational Cost: add. (630->267), mult. (1772->379), div. (0->0), fcn. (1704->8), ass. (0->122)
t87 = cos(qJ(4));
t145 = 0.2e1 * t87;
t88 = cos(qJ(2));
t144 = 0.2e1 * t88;
t84 = sin(qJ(2));
t109 = t84 * t87;
t82 = sin(qJ(5));
t86 = cos(qJ(5));
t41 = t82 * t109 + t88 * t86;
t108 = t86 * t87;
t42 = t84 * t108 - t88 * t82;
t143 = (t88 * mrSges(5,1) + t41 * mrSges(6,1) + t42 * mrSges(6,2) + mrSges(5,3) * t109) * qJ(3);
t75 = t82 ^ 2;
t79 = t86 ^ 2;
t102 = t75 + t79;
t81 = sin(qJ(6));
t85 = cos(qJ(6));
t44 = (mrSges(7,1) * t81 + mrSges(7,2) * t85) * t82;
t142 = t102 * mrSges(6,3) + t44 * t82 - mrSges(5,2);
t141 = 2 * mrSges(4,3);
t83 = sin(qJ(4));
t112 = t83 * t84;
t12 = Ifges(6,1) * t42 - Ifges(6,4) * t41 + Ifges(6,5) * t112;
t140 = t12 / 0.2e1;
t119 = Ifges(7,4) * t85;
t27 = -Ifges(7,6) * t86 + (-Ifges(7,2) * t81 + t119) * t82;
t139 = t27 / 0.2e1;
t120 = Ifges(7,4) * t81;
t30 = -Ifges(7,5) * t86 + (Ifges(7,1) * t85 - t120) * t82;
t138 = t30 / 0.2e1;
t122 = Ifges(6,4) * t82;
t31 = -Ifges(6,5) * t87 + (Ifges(6,1) * t86 - t122) * t83;
t137 = t31 / 0.2e1;
t110 = t83 * t86;
t90 = t81 * t110 + t85 * t87;
t136 = -t90 / 0.2e1;
t43 = t85 * t110 - t81 * t87;
t135 = t43 / 0.2e1;
t55 = Ifges(7,5) * t81 + Ifges(7,6) * t85;
t134 = t55 / 0.2e1;
t57 = Ifges(7,2) * t85 + t120;
t133 = t57 / 0.2e1;
t60 = Ifges(7,1) * t81 + t119;
t132 = t60 / 0.2e1;
t121 = Ifges(6,4) * t86;
t61 = Ifges(6,1) * t82 + t121;
t131 = t61 / 0.2e1;
t130 = -t81 / 0.2e1;
t129 = t81 / 0.2e1;
t128 = t85 / 0.2e1;
t111 = t83 * t85;
t17 = t84 * t111 - t42 * t81;
t18 = t81 * t112 + t42 * t85;
t127 = mrSges(6,1) * t112 + t17 * mrSges(7,1) - t18 * mrSges(7,2) - t42 * mrSges(6,3);
t126 = mrSges(6,2) * t86;
t125 = mrSges(7,3) * t82;
t124 = Ifges(5,4) * t83;
t123 = Ifges(5,4) * t87;
t118 = Ifges(5,6) * t88;
t20 = -mrSges(6,2) * t112 - t41 * mrSges(6,3);
t117 = t20 * t86;
t113 = t82 * t83;
t47 = t87 * mrSges(6,2) - mrSges(6,3) * t113;
t116 = t47 * t86;
t77 = t84 ^ 2;
t89 = qJ(3) ^ 2;
t115 = t77 * t89;
t80 = t87 ^ 2;
t114 = t80 * t89;
t107 = -t86 * mrSges(6,1) + t82 * mrSges(6,2) - mrSges(5,1);
t106 = t85 * mrSges(7,1) - t81 * mrSges(7,2) + mrSges(6,1);
t105 = t87 * mrSges(6,1) + mrSges(7,1) * t90 + t43 * mrSges(7,2) + mrSges(6,3) * t110;
t56 = Ifges(6,5) * t82 + Ifges(6,6) * t86;
t103 = t81 ^ 2 + t85 ^ 2;
t101 = qJ(3) * t84;
t100 = qJ(3) * t87;
t99 = 0.2e1 * t82;
t76 = t83 ^ 2;
t98 = t76 * t115;
t1 = Ifges(7,5) * t18 + Ifges(7,6) * t17 + Ifges(7,3) * t41;
t7 = Ifges(7,5) * t43 - Ifges(7,6) * t90 + Ifges(7,3) * t113;
t8 = Ifges(6,5) * t42 - Ifges(6,6) * t41 + Ifges(6,3) * t112;
t97 = t83 * t101;
t10 = Ifges(6,4) * t42 - Ifges(6,2) * t41 + Ifges(6,6) * t112;
t96 = -t10 / 0.2e1 + t1 / 0.2e1;
t28 = -Ifges(6,6) * t87 + (-Ifges(6,2) * t82 + t121) * t83;
t95 = -t28 / 0.2e1 + t7 / 0.2e1;
t25 = -Ifges(7,3) * t86 + (Ifges(7,5) * t85 - Ifges(7,6) * t81) * t82;
t58 = Ifges(6,2) * t86 + t122;
t94 = t25 / 0.2e1 - t58 / 0.2e1;
t23 = t90 * t101;
t24 = t43 * t101;
t92 = -t23 * t81 - t24 * t85;
t33 = (-t81 * t108 + t111) * qJ(3);
t34 = (t85 * t108 + t81 * t83) * qJ(3);
t91 = -t33 * t81 + t34 * t85;
t73 = Ifges(5,5) * t83;
t71 = Ifges(5,6) * t87;
t69 = t76 * t89;
t68 = Ifges(5,5) * t109;
t67 = Ifges(6,5) * t110;
t64 = t77 * t114;
t63 = t75 * t114;
t62 = Ifges(5,1) * t83 + t123;
t59 = Ifges(5,2) * t87 + t124;
t52 = t75 * t98;
t49 = -t86 * mrSges(7,1) - t85 * t125;
t48 = t88 * mrSges(5,2) - mrSges(5,3) * t112;
t46 = t86 * mrSges(7,2) - t81 * t125;
t45 = (mrSges(6,1) * t82 + t126) * t83;
t32 = -Ifges(5,5) * t88 + (Ifges(5,1) * t87 - t124) * t84;
t29 = -t118 + (-Ifges(5,2) * t83 + t123) * t84;
t26 = -Ifges(6,6) * t113 - Ifges(6,3) * t87 + t67;
t22 = mrSges(7,1) * t113 - t43 * mrSges(7,3);
t19 = -mrSges(7,2) * t113 - mrSges(7,3) * t90;
t11 = Ifges(7,1) * t43 - Ifges(7,4) * t90 + Ifges(7,5) * t113;
t9 = Ifges(7,4) * t43 - Ifges(7,2) * t90 + Ifges(7,6) * t113;
t6 = t41 * mrSges(7,1) - t18 * mrSges(7,3);
t5 = -t41 * mrSges(7,2) + t17 * mrSges(7,3);
t3 = Ifges(7,1) * t18 + Ifges(7,4) * t17 + Ifges(7,5) * t41;
t2 = Ifges(7,4) * t18 + Ifges(7,2) * t17 + Ifges(7,6) * t41;
t4 = [m(4) * t115 + t42 * t12 + t17 * t2 + t18 * t3 + 0.2e1 * t23 * t6 - 0.2e1 * t24 * t5 + Ifges(2,3) + (-t10 + t1) * t41 + m(5) * (t64 + t98) + m(6) * (t79 * t98 + t52 + t64) + m(7) * (t23 ^ 2 + t24 ^ 2 + t52) + (-t68 + (Ifges(5,3) + Ifges(4,3) + Ifges(3,2)) * t88) * t88 + (t87 * t32 + (Ifges(3,1) + Ifges(4,1)) * t84 + (-t29 + t8 + t118) * t83 + (mrSges(4,1) * t144 + t84 * t141 + (t127 * t99 - 0.2e1 * t117 - 0.2e1 * t48) * t83) * qJ(3) + (Ifges(3,4) - Ifges(4,5)) * t144 + t143 * t145) * t84; m(7) * (t33 * t23 - t34 * t24) + t42 * t137 + t3 * t135 + t33 * t6 + t34 * t5 + t2 * t136 + t23 * t22 - t24 * t19 + t17 * t9 / 0.2e1 + t18 * t11 / 0.2e1 + (Ifges(4,4) + Ifges(3,5)) * t84 + t95 * t41 + (qJ(3) * mrSges(4,2) - t73 / 0.2e1 - t71 / 0.2e1 - Ifges(4,6) + Ifges(3,6)) * t88 + (t84 * t62 / 0.2e1 + t29 / 0.2e1 - t8 / 0.2e1 + (-t127 * t82 + t84 * t45 + t117 + t48) * qJ(3)) * t87 + (t86 * t140 + t32 / 0.2e1 + t96 * t82 + t143 + (-t59 / 0.2e1 + t26 / 0.2e1 + (-t105 * t82 - t116) * qJ(3) + (-m(7) * t75 / 0.2e1 + m(6) * (0.1e1 - t102) / 0.2e1) * t89 * t145) * t84) * t83; m(4) * t89 + t43 * t11 + 0.2e1 * t34 * t19 + 0.2e1 * t33 * t22 - t90 * t9 + Ifges(4,2) + Ifges(3,3) + (-t26 + t59) * t87 + m(7) * (t33 ^ 2 + t34 ^ 2 + t63) + m(6) * (t79 * t114 + t63 + t69) + m(5) * (t69 + t114) + (t86 * t31 + t62 + (-t28 + t7) * t82) * t83 + (0.2e1 * t83 * t45 + t141 + 0.2e1 * (t76 + t80) * mrSges(5,3) + (t105 * t99 + 0.2e1 * t116) * t87) * qJ(3); t127 * t86 + (t83 * mrSges(5,1) + t87 * mrSges(5,2) + mrSges(4,2)) * t84 + (m(7) * (t86 * t97 + t92) + t85 * t5 - t81 * t6 + t20) * t82; -t87 * mrSges(5,1) + t83 * mrSges(5,2) - mrSges(4,1) - t105 * t86 + (m(7) * (-t86 * t100 + t91) + t85 * t19 - t81 * t22 + t47) * t82; m(4) + m(5) + m(6) * t102 + m(7) * (t103 * t75 + t79); t68 - Ifges(5,3) * t88 + t17 * t139 - t24 * t46 + t23 * t49 + t18 * t138 + t42 * t131 - t96 * t86 + t94 * t41 + (t3 * t128 + t2 * t130 + t140) * t82 + ((t56 / 0.2e1 - Ifges(5,6)) * t83 + (t107 * t87 - t142 * t83) * qJ(3)) * t84; t34 * t46 + t33 * t49 + t27 * t136 + t30 * t135 - t87 * t56 / 0.2e1 + t73 + t71 + (t83 * t131 - t95) * t86 + (t11 * t128 + t9 * t130 + t94 * t83 + t137) * t82 + (t107 * t83 + t142 * t87) * qJ(3); -t86 * t44 + (t46 * t85 - t49 * t81) * t82; Ifges(5,3) + (-t25 + t58) * t86 + (-t27 * t81 + t30 * t85 + t61) * t82; t3 * t129 + t2 * t128 + t17 * t133 + t41 * t134 + t18 * t132 + t92 * mrSges(7,3) + (t106 * t82 + t126) * t97 + t8; -t90 * t133 + t43 * t132 + t9 * t128 + t11 * t129 + t67 + (-qJ(3) * t126 - Ifges(6,3)) * t87 + t91 * mrSges(7,3) + ((-Ifges(6,6) + t134) * t83 - t106 * t100) * t82; t106 * t86 + (t103 * mrSges(7,3) - mrSges(6,2)) * t82; -t86 * t55 / 0.2e1 + (t82 * t132 + t139) * t85 + (t138 - t82 * t57 / 0.2e1) * t81 + t56; t85 * t57 + t81 * t60 + Ifges(6,3); t23 * mrSges(7,1) + t24 * mrSges(7,2) + t1; t33 * mrSges(7,1) - t34 * mrSges(7,2) + t7; -t44; t25; t55; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t4(1) t4(2) t4(4) t4(7) t4(11) t4(16); t4(2) t4(3) t4(5) t4(8) t4(12) t4(17); t4(4) t4(5) t4(6) t4(9) t4(13) t4(18); t4(7) t4(8) t4(9) t4(10) t4(14) t4(19); t4(11) t4(12) t4(13) t4(14) t4(15) t4(20); t4(16) t4(17) t4(18) t4(19) t4(20) t4(21);];
Mq  = res;
