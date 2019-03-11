% Calculate joint inertia matrix for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:50
% EndTime: 2019-03-09 07:04:51
% DurationCPUTime: 0.80s
% Computational Cost: add. (3099->209), mult. (5829->308), div. (0->0), fcn. (6993->10), ass. (0->105)
t85 = sin(qJ(6));
t121 = t85 * mrSges(7,3);
t83 = sin(pkin(11));
t84 = cos(pkin(11));
t88 = sin(qJ(3));
t92 = cos(qJ(3));
t55 = -t83 * t88 + t84 * t92;
t56 = t83 * t92 + t84 * t88;
t87 = sin(qJ(4));
t91 = cos(qJ(4));
t39 = t55 * t91 - t56 * t87;
t40 = t55 * t87 + t56 * t91;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t30 = -t90 * t39 + t40 * t86;
t31 = t39 * t86 + t40 * t90;
t15 = -mrSges(7,2) * t30 - t31 * t121;
t89 = cos(qJ(6));
t123 = t31 * t89;
t16 = mrSges(7,1) * t30 - mrSges(7,3) * t123;
t142 = t89 * t15 - t85 * t16;
t63 = -mrSges(7,1) * t89 + mrSges(7,2) * t85;
t130 = pkin(5) * t63;
t81 = t85 ^ 2;
t128 = mrSges(7,3) * t81;
t73 = pkin(10) * t128;
t82 = t89 ^ 2;
t127 = mrSges(7,3) * t82;
t74 = pkin(10) * t127;
t141 = t73 + t74 - t130;
t104 = mrSges(7,1) * t85 + mrSges(7,2) * t89;
t14 = t104 * t31;
t117 = pkin(7) + qJ(2);
t59 = t117 * t83;
t60 = t117 * t84;
t41 = -t59 * t92 - t60 * t88;
t101 = -pkin(8) * t56 + t41;
t42 = -t59 * t88 + t60 * t92;
t102 = pkin(8) * t55 + t42;
t21 = t87 * t101 + t91 * t102;
t18 = pkin(9) * t39 + t21;
t20 = t91 * t101 - t87 * t102;
t97 = -t40 * pkin(9) + t20;
t6 = t18 * t86 - t90 * t97;
t140 = m(7) * t6 + t14;
t131 = pkin(4) * t90;
t71 = -pkin(5) - t131;
t53 = t71 * t63;
t70 = pkin(4) * t86 + pkin(10);
t61 = t70 * t128;
t62 = t70 * t127;
t75 = mrSges(6,1) * t131;
t139 = t53 + t61 + t62 + t75;
t69 = -pkin(2) * t84 - pkin(1);
t46 = -pkin(3) * t55 + t69;
t32 = -pkin(4) * t39 + t46;
t13 = pkin(5) * t30 - pkin(10) * t31 + t32;
t8 = t90 * t18 + t86 * t97;
t3 = t13 * t85 + t8 * t89;
t129 = t3 * t89;
t2 = t13 * t89 - t8 * t85;
t105 = -t2 * t85 + t129;
t138 = m(7) * t105 + t142;
t137 = t6 ^ 2;
t136 = 0.2e1 * t6;
t80 = t84 ^ 2;
t135 = 0.2e1 * t32;
t134 = 0.2e1 * t39;
t133 = 0.2e1 * t55;
t132 = pkin(3) * t87;
t126 = Ifges(7,4) * t85;
t125 = Ifges(7,4) * t89;
t124 = t31 * t85;
t72 = pkin(3) * t91 + pkin(4);
t52 = t90 * t132 + t86 * t72;
t122 = t52 * mrSges(6,2);
t119 = t86 * mrSges(6,2);
t116 = Ifges(7,5) * t123 + Ifges(7,3) * t30;
t115 = Ifges(7,5) * t85 + Ifges(7,6) * t89;
t114 = t83 ^ 2 + t80;
t113 = t81 + t82;
t112 = pkin(4) * t119;
t64 = Ifges(7,2) * t89 + t126;
t65 = Ifges(7,1) * t85 + t125;
t111 = t89 * t64 + t85 * t65 + Ifges(6,3);
t110 = t113 * t70;
t109 = -t84 * mrSges(3,1) + t83 * mrSges(3,2);
t108 = -t55 * mrSges(4,1) + t56 * mrSges(4,2);
t107 = -t39 * mrSges(5,1) + t40 * mrSges(5,2);
t106 = Ifges(5,3) + t111;
t51 = -t86 * t132 + t72 * t90;
t103 = (mrSges(5,1) * t91 - mrSges(5,2) * t87) * pkin(3);
t49 = -pkin(5) - t51;
t43 = t49 * t63;
t50 = pkin(10) + t52;
t44 = t50 * t128;
t45 = t50 * t127;
t47 = t51 * mrSges(6,1);
t100 = t111 + t43 + t44 + t45 + t47 - t122;
t11 = Ifges(7,6) * t30 + (-Ifges(7,2) * t85 + t125) * t31;
t12 = Ifges(7,5) * t30 + (Ifges(7,1) * t89 - t126) * t31;
t99 = -t8 * mrSges(6,2) + mrSges(7,3) * t129 - t2 * t121 + t89 * t11 / 0.2e1 - t64 * t124 / 0.2e1 + t65 * t123 / 0.2e1 + Ifges(6,5) * t31 + t85 * t12 / 0.2e1 + (-mrSges(6,1) + t63) * t6 + (t115 / 0.2e1 - Ifges(6,6)) * t30;
t98 = t20 * mrSges(5,1) - t21 * mrSges(5,2) + Ifges(5,5) * t40 + Ifges(5,6) * t39 + t99;
t26 = t31 * mrSges(6,2);
t1 = [Ifges(3,2) * t80 - 0.2e1 * pkin(1) * t109 + 0.2e1 * t69 * t108 + Ifges(4,2) * t55 ^ 2 + Ifges(5,2) * t39 ^ 2 + 0.2e1 * t46 * t107 + t42 * mrSges(4,3) * t133 + t21 * mrSges(5,3) * t134 + t26 * t135 + 0.2e1 * t2 * t16 + t14 * t136 + 0.2e1 * t3 * t15 + Ifges(2,3) + (Ifges(3,1) * t83 + 0.2e1 * Ifges(3,4) * t84) * t83 + (mrSges(6,1) * t135 - 0.2e1 * t8 * mrSges(6,3) + Ifges(6,2) * t30 + t116) * t30 + 0.2e1 * t114 * qJ(2) * mrSges(3,3) + (-0.2e1 * t41 * mrSges(4,3) + Ifges(4,1) * t56 + Ifges(4,4) * t133) * t56 + (-0.2e1 * t20 * mrSges(5,3) + Ifges(5,1) * t40 + Ifges(5,4) * t134) * t40 + (mrSges(6,3) * t136 + Ifges(6,1) * t31 - t85 * t11 + t89 * t12 + (-Ifges(7,6) * t85 - (2 * Ifges(6,4))) * t30) * t31 + m(3) * (t114 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t41 ^ 2 + t42 ^ 2 + t69 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2 + t46 ^ 2) + m(6) * (t32 ^ 2 + t8 ^ 2 + t137) + m(7) * (t2 ^ 2 + t3 ^ 2 + t137); -m(3) * pkin(1) + t30 * mrSges(6,1) + t85 * t15 + t89 * t16 + t26 + m(7) * (t2 * t89 + t3 * t85) + m(6) * t32 + m(5) * t46 + m(4) * t69 + t107 + t108 + t109; m(7) * t113 + m(3) + m(4) + m(5) + m(6); Ifges(4,6) * t55 + Ifges(4,5) * t56 + t41 * mrSges(4,1) - t42 * mrSges(4,2) + t49 * t14 + (m(5) * (t20 * t91 + t21 * t87) + (t87 * t39 - t91 * t40) * mrSges(5,3)) * pkin(3) + t98 + m(7) * (t105 * t50 + t49 * t6) + (-t52 * t30 - t51 * t31) * mrSges(6,3) + m(6) * (-t51 * t6 + t52 * t8) + t142 * t50; 0; -0.2e1 * t122 + Ifges(4,3) + 0.2e1 * t43 + 0.2e1 * t44 + 0.2e1 * t45 + 0.2e1 * t47 + 0.2e1 * t103 + m(7) * (t113 * t50 ^ 2 + t49 ^ 2) + m(6) * (t51 ^ 2 + t52 ^ 2) + m(5) * (t87 ^ 2 + t91 ^ 2) * pkin(3) ^ 2 + t106; (m(6) * (-t6 * t90 + t8 * t86) + (-t86 * t30 - t90 * t31) * mrSges(6,3)) * pkin(4) + t98 + t140 * t71 + t138 * t70; 0; t103 + t100 + m(7) * (t50 * t110 + t49 * t71) + (m(6) * (t51 * t90 + t52 * t86) - t119) * pkin(4) + Ifges(5,3) + t139; -0.2e1 * t112 + 0.2e1 * t53 + 0.2e1 * t61 + 0.2e1 * t62 + 0.2e1 * t75 + m(7) * (t113 * t70 ^ 2 + t71 ^ 2) + m(6) * (t86 ^ 2 + t90 ^ 2) * pkin(4) ^ 2 + t106; -t140 * pkin(5) + t138 * pkin(10) + t99; 0; m(7) * (t113 * t50 * pkin(10) - pkin(5) * t49) + t100 + t141; m(7) * (-pkin(5) * t71 + pkin(10) * t110) - t112 + t111 + t139 + t141; -0.2e1 * t130 + m(7) * (t113 * pkin(10) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t74 + 0.2e1 * t73 + t111; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t124 + t116; -t63; -t104 * t50 + t115; -t104 * t70 + t115; -t104 * pkin(10) + t115; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
