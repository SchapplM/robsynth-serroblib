% Calculate joint inertia matrix for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:47
% EndTime: 2019-03-08 18:59:49
% DurationCPUTime: 0.79s
% Computational Cost: add. (1264->223), mult. (3154->339), div. (0->0), fcn. (3618->14), ass. (0->97)
t84 = sin(qJ(6));
t88 = cos(qJ(6));
t58 = -mrSges(7,1) * t88 + mrSges(7,2) * t84;
t139 = t58 - mrSges(6,1);
t104 = mrSges(7,1) * t84 + mrSges(7,2) * t88;
t85 = sin(qJ(5));
t86 = sin(qJ(4));
t89 = cos(qJ(5));
t90 = cos(qJ(4));
t55 = t85 * t90 + t86 * t89;
t28 = t104 * t55;
t131 = -pkin(10) - pkin(9);
t107 = t131 * t86;
t63 = t131 * t90;
t40 = -t89 * t107 - t85 * t63;
t138 = m(7) * t40 + t28;
t54 = t85 * t86 - t89 * t90;
t69 = -pkin(4) * t90 - pkin(3);
t29 = pkin(5) * t54 - pkin(11) * t55 + t69;
t42 = t85 * t107 - t89 * t63;
t10 = t29 * t88 - t42 * t84;
t11 = t29 * t84 + t42 * t88;
t125 = t11 * t88;
t117 = t84 * mrSges(7,3);
t31 = -mrSges(7,2) * t54 - t55 * t117;
t121 = t55 * t88;
t32 = mrSges(7,1) * t54 - mrSges(7,3) * t121;
t137 = m(7) * (-t10 * t84 + t125) + t88 * t31 - t84 * t32;
t81 = cos(pkin(13));
t82 = cos(pkin(7));
t118 = t81 * t82;
t79 = sin(pkin(7));
t87 = sin(qJ(3));
t120 = t79 * t87;
t78 = sin(pkin(13));
t80 = sin(pkin(6));
t83 = cos(pkin(6));
t91 = cos(qJ(3));
t35 = t83 * t120 + (t87 * t118 + t78 * t91) * t80;
t45 = -t79 * t80 * t81 + t82 * t83;
t14 = -t35 * t86 + t45 * t90;
t15 = t35 * t90 + t45 * t86;
t6 = -t89 * t14 + t15 * t85;
t136 = t6 ^ 2;
t47 = -t86 * t120 + t82 * t90;
t48 = t90 * t120 + t82 * t86;
t24 = -t89 * t47 + t48 * t85;
t135 = t24 ^ 2;
t119 = t79 * t91;
t33 = -t83 * t119 + (-t118 * t91 + t78 * t87) * t80;
t30 = t33 ^ 2;
t134 = t40 ^ 2;
t77 = t90 ^ 2;
t133 = 0.2e1 * t40;
t132 = m(6) * pkin(4);
t130 = t24 * t6;
t8 = t14 * t85 + t15 * t89;
t3 = t33 * t84 + t8 * t88;
t129 = t3 * t88;
t128 = t40 * t6;
t127 = Ifges(7,4) * t84;
t126 = Ifges(7,4) * t88;
t26 = t47 * t85 + t48 * t89;
t19 = -t84 * t119 + t26 * t88;
t124 = t19 * t88;
t123 = t24 * t40;
t122 = t55 * t84;
t114 = Ifges(7,5) * t121 + Ifges(7,3) * t54;
t113 = Ifges(7,5) * t84 + Ifges(7,6) * t88;
t112 = t84 ^ 2 + t88 ^ 2;
t111 = t86 ^ 2 + t77;
t110 = t33 * t119;
t36 = mrSges(6,1) * t54 + mrSges(6,2) * t55;
t59 = -mrSges(5,1) * t90 + mrSges(5,2) * t86;
t109 = mrSges(4,1) - t36 - t59;
t60 = Ifges(7,2) * t88 + t127;
t61 = Ifges(7,1) * t84 + t126;
t108 = t88 * t60 + t84 * t61 + Ifges(6,3);
t67 = pkin(4) * t85 + pkin(11);
t106 = t112 * t67;
t2 = t33 * t88 - t8 * t84;
t105 = -t2 * t84 + t129;
t102 = -t14 * t86 + t15 * t90;
t18 = -t88 * t119 - t26 * t84;
t101 = -t18 * t84 + t124;
t100 = -t47 * t86 + t48 * t90;
t99 = 0.2e1 * t112 * mrSges(7,3);
t98 = (mrSges(6,1) * t89 - mrSges(6,2) * t85) * pkin(4);
t97 = -t8 * mrSges(6,2) + mrSges(7,3) * t129 - t2 * t117 + t139 * t6;
t96 = -t26 * mrSges(6,2) + mrSges(7,3) * t124 - t18 * t117 + t139 * t24;
t21 = Ifges(7,6) * t54 + (-Ifges(7,2) * t84 + t126) * t55;
t22 = Ifges(7,5) * t54 + (Ifges(7,1) * t88 - t127) * t55;
t95 = -t42 * mrSges(6,2) + mrSges(7,3) * t125 - t10 * t117 + t84 * t22 / 0.2e1 + t88 * t21 / 0.2e1 - t60 * t122 / 0.2e1 + t61 * t121 / 0.2e1 + Ifges(6,5) * t55 + (t113 / 0.2e1 - Ifges(6,6)) * t54 + t139 * t40;
t72 = t79 ^ 2;
t68 = -pkin(4) * t89 - pkin(5);
t65 = t72 * t91 ^ 2;
t1 = [m(2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t136) + m(6) * (t8 ^ 2 + t136 + t30) + m(5) * (t14 ^ 2 + t15 ^ 2 + t30) + m(4) * (t35 ^ 2 + t45 ^ 2 + t30) + m(3) * (t83 ^ 2 + (t78 ^ 2 + t81 ^ 2) * t80 ^ 2); m(3) * t83 + m(7) * (t18 * t2 + t19 * t3 + t130) + m(6) * (t26 * t8 - t110 + t130) + m(5) * (t14 * t47 + t15 * t48 - t110) + m(4) * (t82 * t45 + (-t33 * t91 + t35 * t87) * t79); m(3) + m(7) * (t18 ^ 2 + t19 ^ 2 + t135) + m(6) * (t26 ^ 2 + t135 + t65) + m(5) * (t47 ^ 2 + t48 ^ 2 + t65) + m(4) * (t72 * t87 ^ 2 + t82 ^ 2 + t65); -t35 * mrSges(4,2) + t2 * t32 + t6 * t28 + t3 * t31 + (-t54 * t8 + t55 * t6) * mrSges(6,3) + t102 * mrSges(5,3) - t109 * t33 + m(7) * (t10 * t2 + t11 * t3 + t128) + m(6) * (t33 * t69 + t42 * t8 + t128) + m(5) * (-pkin(3) * t33 + t102 * pkin(9)); t18 * t32 + t19 * t31 + t24 * t28 + (t24 * t55 - t26 * t54) * mrSges(6,3) + t100 * mrSges(5,3) + (-t87 * mrSges(4,2) + t109 * t91) * t79 + m(7) * (t10 * t18 + t11 * t19 + t123) + m(6) * (-t69 * t119 + t26 * t42 + t123) + m(5) * (pkin(3) * t119 + t100 * pkin(9)); Ifges(5,2) * t77 - 0.2e1 * pkin(3) * t59 + 0.2e1 * t10 * t32 + 0.2e1 * t11 * t31 + t28 * t133 + 0.2e1 * t69 * t36 + Ifges(4,3) + (Ifges(5,1) * t86 + 0.2e1 * Ifges(5,4) * t90) * t86 + 0.2e1 * t111 * pkin(9) * mrSges(5,3) + (-0.2e1 * mrSges(6,3) * t42 + Ifges(6,2) * t54 + t114) * t54 + m(7) * (t10 ^ 2 + t11 ^ 2 + t134) + m(6) * (t42 ^ 2 + t69 ^ 2 + t134) + m(5) * (t111 * pkin(9) ^ 2 + pkin(3) ^ 2) + (mrSges(6,3) * t133 + Ifges(6,1) * t55 - t84 * t21 + t88 * t22 + (-Ifges(7,6) * t84 - (2 * Ifges(6,4))) * t54) * t55; t14 * mrSges(5,1) - t15 * mrSges(5,2) + m(7) * (t105 * t67 + t68 * t6) + (-t6 * t89 + t8 * t85) * t132 + t97; t47 * mrSges(5,1) - t48 * mrSges(5,2) + m(7) * (t101 * t67 + t68 * t24) + (-t24 * t89 + t26 * t85) * t132 + t96; (m(6) * (-t40 * t89 + t42 * t85) + (-t54 * t85 - t55 * t89) * mrSges(6,3)) * pkin(4) + (-mrSges(5,1) * t86 - mrSges(5,2) * t90) * pkin(9) + Ifges(5,6) * t90 + Ifges(5,5) * t86 + t95 + t138 * t68 + t137 * t67; 0.2e1 * t68 * t58 + Ifges(5,3) + 0.2e1 * t98 + t67 * t99 + m(7) * (t112 * t67 ^ 2 + t68 ^ 2) + m(6) * (t85 ^ 2 + t89 ^ 2) * pkin(4) ^ 2 + t108; m(7) * (-pkin(5) * t6 + t105 * pkin(11)) + t97; m(7) * (-pkin(5) * t24 + t101 * pkin(11)) + t96; -t138 * pkin(5) + t137 * pkin(11) + t95; m(7) * (-pkin(5) * t68 + pkin(11) * t106) + (-pkin(5) + t68) * t58 + t98 + (t112 * pkin(11) + t106) * mrSges(7,3) + t108; -0.2e1 * pkin(5) * t58 + m(7) * (t112 * pkin(11) ^ 2 + pkin(5) ^ 2) + pkin(11) * t99 + t108; mrSges(7,1) * t2 - mrSges(7,2) * t3; mrSges(7,1) * t18 - mrSges(7,2) * t19; mrSges(7,1) * t10 - mrSges(7,2) * t11 - Ifges(7,6) * t122 + t114; -t104 * t67 + t113; -t104 * pkin(11) + t113; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
