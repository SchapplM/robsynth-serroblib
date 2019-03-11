% Calculate joint inertia matrix for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:16
% EndTime: 2019-03-09 15:21:18
% DurationCPUTime: 1.05s
% Computational Cost: add. (2970->271), mult. (5561->391), div. (0->0), fcn. (6442->10), ass. (0->102)
t111 = sin(pkin(10));
t113 = cos(pkin(10));
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t116 = sin(qJ(2));
t150 = -pkin(8) - pkin(7);
t132 = t150 * t116;
t119 = cos(qJ(2));
t133 = t150 * t119;
t67 = t115 * t133 + t118 * t132;
t89 = t115 * t119 + t116 * t118;
t124 = -t89 * qJ(4) + t67;
t68 = t115 * t132 - t118 * t133;
t88 = -t115 * t116 + t118 * t119;
t50 = qJ(4) * t88 + t68;
t30 = t111 * t50 - t113 * t124;
t156 = t30 ^ 2;
t155 = 0.2e1 * t30;
t110 = sin(pkin(11));
t112 = cos(pkin(11));
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t86 = -t110 * t114 + t112 * t117;
t87 = t110 * t117 + t112 * t114;
t61 = -t86 * mrSges(7,1) + mrSges(7,2) * t87;
t154 = 0.2e1 * t61;
t103 = -pkin(2) * t119 - pkin(1);
t70 = -pkin(3) * t88 + t103;
t153 = 0.2e1 * t70;
t93 = -t112 * mrSges(6,1) + mrSges(6,2) * t110;
t152 = 0.2e1 * t93;
t149 = pkin(2) * t115;
t148 = pkin(5) * t112;
t102 = pkin(2) * t118 + pkin(3);
t74 = t102 * t113 - t111 * t149;
t147 = t74 * mrSges(5,1);
t75 = t102 * t111 + t113 * t149;
t146 = t75 * mrSges(5,2);
t59 = t111 * t89 - t113 * t88;
t60 = t111 * t88 + t113 * t89;
t28 = pkin(4) * t59 - qJ(5) * t60 + t70;
t32 = t111 * t124 + t113 * t50;
t12 = t110 * t28 + t112 * t32;
t139 = t112 * t60;
t141 = t110 * t60;
t38 = mrSges(6,1) * t141 + mrSges(6,2) * t139;
t145 = Ifges(7,5) * t87 + Ifges(7,6) * t86;
t144 = Ifges(6,4) * t110;
t143 = Ifges(6,4) * t112;
t11 = -t110 * t32 + t112 * t28;
t142 = t11 * t110;
t140 = t112 * t12;
t72 = qJ(5) + t75;
t138 = t112 * t72;
t137 = t110 ^ 2 + t112 ^ 2;
t136 = t116 ^ 2 + t119 ^ 2;
t135 = 2 * mrSges(7,3);
t36 = t87 * t60;
t37 = t86 * t60;
t134 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t59;
t101 = -pkin(3) * t113 - pkin(4);
t13 = t36 * mrSges(7,1) + mrSges(7,2) * t37;
t100 = pkin(3) * t111 + qJ(5);
t131 = t137 * t100;
t73 = -pkin(4) - t74;
t130 = t113 * mrSges(5,1) - t111 * mrSges(5,2);
t129 = t140 - t142;
t128 = 0.2e1 * t137 * mrSges(6,3);
t62 = Ifges(7,4) * t87 + Ifges(7,2) * t86;
t63 = Ifges(7,1) * t87 + Ifges(7,4) * t86;
t94 = Ifges(6,2) * t112 + t144;
t95 = Ifges(6,1) * t110 + t143;
t127 = t110 * t95 + t112 * t94 + t62 * t86 + t63 * t87 + Ifges(4,3) + Ifges(5,3);
t126 = (mrSges(4,1) * t118 - mrSges(4,2) * t115) * pkin(2);
t125 = t61 + t93;
t15 = pkin(5) * t141 + t30;
t4 = pkin(5) * t59 - pkin(9) * t139 + t11;
t5 = -pkin(9) * t141 + t12;
t2 = -t114 * t5 + t117 * t4;
t22 = t59 * Ifges(6,6) + (-Ifges(6,2) * t110 + t143) * t60;
t23 = t59 * Ifges(6,5) + (Ifges(6,1) * t112 - t144) * t60;
t3 = t114 * t4 + t117 * t5;
t8 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t59;
t9 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t59;
t123 = -t68 * mrSges(4,2) - t32 * mrSges(5,2) + mrSges(6,3) * t140 + t15 * t61 + t110 * t23 / 0.2e1 + t112 * t22 / 0.2e1 - t36 * t62 / 0.2e1 + t37 * t63 / 0.2e1 - t94 * t141 / 0.2e1 + t95 * t139 / 0.2e1 - Ifges(5,6) * t59 + Ifges(5,5) * t60 + t86 * t8 / 0.2e1 + t67 * mrSges(4,1) + t87 * t9 / 0.2e1 + Ifges(4,6) * t88 + Ifges(4,5) * t89 + (t93 - mrSges(5,1)) * t30 + (Ifges(6,5) * t110 + Ifges(6,6) * t112 + t145) * t59 / 0.2e1 + (-t2 * t87 + t3 * t86) * mrSges(7,3);
t105 = t112 * pkin(9);
t92 = t101 - t148;
t78 = t100 * t112 + t105;
t77 = (-pkin(9) - t100) * t110;
t69 = t73 - t148;
t65 = t105 + t138;
t64 = (-pkin(9) - t72) * t110;
t55 = t60 * mrSges(5,2);
t52 = t114 * t77 + t117 * t78;
t51 = -t114 * t78 + t117 * t77;
t43 = t114 * t64 + t117 * t65;
t42 = -t114 * t65 + t117 * t64;
t40 = mrSges(6,1) * t59 - mrSges(6,3) * t139;
t39 = -mrSges(6,2) * t59 - mrSges(6,3) * t141;
t17 = mrSges(7,1) * t59 - mrSges(7,3) * t37;
t16 = -mrSges(7,2) * t59 - mrSges(7,3) * t36;
t1 = [0.2e1 * t136 * pkin(7) * mrSges(3,3) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t119 + mrSges(3,2) * t116) + t116 * (Ifges(3,1) * t116 + Ifges(3,4) * t119) + t119 * (Ifges(3,4) * t116 + Ifges(3,2) * t119) + 0.2e1 * (-t67 * t89 + t68 * t88) * mrSges(4,3) + 0.2e1 * t103 * (-mrSges(4,1) * t88 + mrSges(4,2) * t89) + t88 * (Ifges(4,4) * t89 + Ifges(4,2) * t88) + t89 * (Ifges(4,1) * t89 + Ifges(4,4) * t88) - t36 * t8 + t37 * t9 + 0.2e1 * t12 * t39 + 0.2e1 * t11 * t40 + 0.2e1 * t15 * t13 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + m(4) * (t103 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(7) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) + t55 * t153 + t38 * t155 + Ifges(2,3) + m(3) * (pkin(7) ^ 2 * t136 + pkin(1) ^ 2) + (mrSges(5,1) * t153 - 0.2e1 * t32 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t59 + t134) * t59 + (mrSges(5,3) * t155 + Ifges(5,1) * t60 - t110 * t22 + t112 * t23 + (Ifges(6,5) * t112 - Ifges(6,6) * t110 - (2 * Ifges(5,4))) * t59) * t60 + m(5) * (t32 ^ 2 + t70 ^ 2 + t156) + m(6) * (t11 ^ 2 + t12 ^ 2 + t156); m(5) * (-t30 * t74 + t32 * t75) + m(7) * (t15 * t69 + t2 * t42 + t3 * t43) + Ifges(3,6) * t119 + (-mrSges(3,1) * t116 - mrSges(3,2) * t119) * pkin(7) + (-t59 * t75 - t60 * t74) * mrSges(5,3) + t39 * t138 + t123 + Ifges(3,5) * t116 + t69 * t13 + t73 * t38 + t42 * t17 + t43 * t16 + m(6) * (t129 * t72 + t30 * t73) + (-t11 * mrSges(6,3) - t72 * t40) * t110 + (m(4) * (t115 * t68 + t118 * t67) + (t115 * t88 - t118 * t89) * mrSges(4,3)) * pkin(2); m(7) * (t42 ^ 2 + t43 ^ 2 + t69 ^ 2) + m(5) * (t74 ^ 2 + t75 ^ 2) + t72 * t128 + t127 + (-t42 * t87 + t43 * t86) * t135 + m(4) * (t115 ^ 2 + t118 ^ 2) * pkin(2) ^ 2 + t73 * t152 + t69 * t154 + 0.2e1 * t147 - 0.2e1 * t146 + m(6) * (t137 * t72 ^ 2 + t73 ^ 2) + 0.2e1 * t126 + Ifges(3,3); -mrSges(6,3) * t142 + t123 + (-t110 * t40 + t112 * t39) * t100 + m(7) * (t15 * t92 + t2 * t51 + t3 * t52) + t92 * t13 + t101 * t38 + t51 * t17 + t52 * t16 + (m(5) * (t111 * t32 - t113 * t30) + (-t111 * t59 - t113 * t60) * mrSges(5,3)) * pkin(3) + m(6) * (t100 * t129 + t101 * t30); t147 - t146 + (t73 + t101) * t93 + (t92 + t69) * t61 + t126 + m(6) * (t101 * t73 + t131 * t72) + m(7) * (t42 * t51 + t43 * t52 + t69 * t92) + (m(5) * (t111 * t75 + t113 * t74) + t130) * pkin(3) + ((-t42 - t51) * t87 + (t43 + t52) * t86) * mrSges(7,3) + (t137 * t72 + t131) * mrSges(6,3) + t127; t101 * t152 + t92 * t154 + (-t51 * t87 + t52 * t86) * t135 + t100 * t128 + m(7) * (t51 ^ 2 + t52 ^ 2 + t92 ^ 2) + m(6) * (t100 ^ 2 * t137 + t101 ^ 2) + t127 + (0.2e1 * t130 + m(5) * (t111 ^ 2 + t113 ^ 2) * pkin(3)) * pkin(3); t59 * mrSges(5,1) + t110 * t39 + t112 * t40 + t87 * t16 + t86 * t17 + t55 + m(7) * (t2 * t86 + t3 * t87) + m(6) * (t11 * t112 + t110 * t12) + m(5) * t70; m(7) * (t42 * t86 + t43 * t87); m(7) * (t51 * t86 + t52 * t87); m(5) + m(6) * t137 + m(7) * (t86 ^ 2 + t87 ^ 2); m(6) * t30 + m(7) * t15 + t13 + t38; m(6) * t73 + m(7) * t69 + t125; m(6) * t101 + m(7) * t92 + t125; 0; m(6) + m(7); mrSges(7,1) * t2 - mrSges(7,2) * t3 + t134; mrSges(7,1) * t42 - mrSges(7,2) * t43 + t145; mrSges(7,1) * t51 - mrSges(7,2) * t52 + t145; -t61; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
