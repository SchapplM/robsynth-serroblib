% Calculate joint inertia matrix for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:29
% EndTime: 2019-03-09 13:12:31
% DurationCPUTime: 0.90s
% Computational Cost: add. (3296->224), mult. (6146->330), div. (0->0), fcn. (7293->10), ass. (0->107)
t88 = sin(qJ(6));
t122 = t88 * mrSges(7,3);
t86 = sin(pkin(11));
t87 = cos(pkin(11));
t91 = sin(qJ(2));
t95 = cos(qJ(2));
t60 = -t86 * t91 + t87 * t95;
t61 = t86 * t95 + t87 * t91;
t90 = sin(qJ(4));
t94 = cos(qJ(4));
t48 = t60 * t94 - t61 * t90;
t49 = t60 * t90 + t61 * t94;
t89 = sin(qJ(5));
t93 = cos(qJ(5));
t30 = -t93 * t48 + t49 * t89;
t31 = t48 * t89 + t49 * t93;
t15 = -mrSges(7,2) * t30 - t31 * t122;
t92 = cos(qJ(6));
t126 = t31 * t92;
t16 = mrSges(7,1) * t30 - mrSges(7,3) * t126;
t145 = t92 * t15 - t88 * t16;
t66 = -mrSges(7,1) * t92 + mrSges(7,2) * t88;
t133 = pkin(5) * t66;
t82 = t88 ^ 2;
t131 = mrSges(7,3) * t82;
t77 = pkin(10) * t131;
t84 = t92 ^ 2;
t130 = mrSges(7,3) * t84;
t78 = pkin(10) * t130;
t144 = t77 + t78 - t133;
t106 = mrSges(7,1) * t88 + mrSges(7,2) * t92;
t14 = t106 * t31;
t118 = -qJ(3) - pkin(7);
t67 = t118 * t91;
t68 = t118 * t95;
t50 = t67 * t87 + t68 * t86;
t104 = -pkin(8) * t61 + t50;
t51 = t67 * t86 - t68 * t87;
t105 = pkin(8) * t60 + t51;
t20 = t104 * t94 - t105 * t90;
t100 = -t49 * pkin(9) + t20;
t21 = t90 * t104 + t94 * t105;
t18 = pkin(9) * t48 + t21;
t6 = -t100 * t93 + t18 * t89;
t143 = m(7) * t6 + t14;
t134 = pkin(4) * t93;
t75 = -pkin(5) - t134;
t58 = t75 * t66;
t74 = pkin(4) * t89 + pkin(10);
t64 = t74 * t131;
t65 = t74 * t130;
t79 = mrSges(6,1) * t134;
t142 = t58 + t64 + t65 + t79;
t76 = -pkin(2) * t95 - pkin(1);
t53 = -pkin(3) * t60 + t76;
t32 = -pkin(4) * t48 + t53;
t13 = pkin(5) * t30 - pkin(10) * t31 + t32;
t8 = t100 * t89 + t93 * t18;
t3 = t13 * t88 + t8 * t92;
t132 = t3 * t92;
t2 = t13 * t92 - t8 * t88;
t107 = -t2 * t88 + t132;
t141 = m(7) * t107 + t145;
t140 = t6 ^ 2;
t139 = 0.2e1 * t6;
t138 = 0.2e1 * t32;
t137 = 0.2e1 * t48;
t136 = 0.2e1 * t60;
t135 = pkin(2) * t86;
t129 = Ifges(7,4) * t88;
t128 = Ifges(7,4) * t92;
t127 = t31 * t88;
t73 = pkin(2) * t87 + pkin(3);
t56 = -t90 * t135 + t94 * t73;
t55 = pkin(4) + t56;
t57 = t94 * t135 + t73 * t90;
t43 = t89 * t55 + t93 * t57;
t125 = t43 * mrSges(6,2);
t124 = t56 * mrSges(5,1);
t123 = t57 * mrSges(5,2);
t120 = t89 * mrSges(6,2);
t117 = Ifges(7,5) * t126 + Ifges(7,3) * t30;
t116 = Ifges(7,5) * t88 + Ifges(7,6) * t92;
t115 = t82 + t84;
t114 = t91 ^ 2 + t95 ^ 2;
t113 = pkin(4) * t120;
t69 = Ifges(7,2) * t92 + t129;
t70 = Ifges(7,1) * t88 + t128;
t112 = t92 * t69 + t88 * t70 + Ifges(6,3);
t111 = t115 * t74;
t110 = -t60 * mrSges(4,1) + t61 * mrSges(4,2);
t109 = -t48 * mrSges(5,1) + t49 * mrSges(5,2);
t108 = Ifges(5,3) + t112;
t42 = t55 * t93 - t57 * t89;
t40 = -pkin(5) - t42;
t33 = t40 * t66;
t41 = pkin(10) + t43;
t34 = t41 * t131;
t35 = t41 * t130;
t38 = t42 * mrSges(6,1);
t103 = t112 + t33 + t34 + t35 + t38 - t125;
t11 = Ifges(7,6) * t30 + (-Ifges(7,2) * t88 + t128) * t31;
t12 = Ifges(7,5) * t30 + (Ifges(7,1) * t92 - t129) * t31;
t102 = -t8 * mrSges(6,2) + mrSges(7,3) * t132 - t2 * t122 + t92 * t11 / 0.2e1 - t69 * t127 / 0.2e1 + t70 * t126 / 0.2e1 + Ifges(6,5) * t31 + t88 * t12 / 0.2e1 + (-mrSges(6,1) + t66) * t6 + (t116 / 0.2e1 - Ifges(6,6)) * t30;
t101 = t20 * mrSges(5,1) - t21 * mrSges(5,2) + Ifges(5,5) * t49 + Ifges(5,6) * t48 + t102;
t26 = t31 * mrSges(6,2);
t1 = [t91 * (Ifges(3,1) * t91 + Ifges(3,4) * t95) + t95 * (Ifges(3,4) * t91 + Ifges(3,2) * t95) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t95 + mrSges(3,2) * t91) + 0.2e1 * t76 * t110 + t51 * mrSges(4,3) * t136 + t21 * mrSges(5,3) * t137 + Ifges(4,2) * t60 ^ 2 + Ifges(5,2) * t48 ^ 2 + 0.2e1 * t53 * t109 + t26 * t138 + t14 * t139 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + (mrSges(6,1) * t138 - 0.2e1 * t8 * mrSges(6,3) + Ifges(6,2) * t30 + t117) * t30 + 0.2e1 * t114 * pkin(7) * mrSges(3,3) + (-0.2e1 * t50 * mrSges(4,3) + Ifges(4,1) * t61 + Ifges(4,4) * t136) * t61 + (-0.2e1 * t20 * mrSges(5,3) + Ifges(5,1) * t49 + Ifges(5,4) * t137) * t49 + (mrSges(6,3) * t139 + Ifges(6,1) * t31 - t88 * t11 + t92 * t12 + (-Ifges(7,6) * t88 - (2 * Ifges(6,4))) * t30) * t31 + m(3) * (t114 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t50 ^ 2 + t51 ^ 2 + t76 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2 + t53 ^ 2) + m(6) * (t32 ^ 2 + t8 ^ 2 + t140) + m(7) * (t2 ^ 2 + t3 ^ 2 + t140); m(7) * (t107 * t41 + t40 * t6) + Ifges(3,5) * t91 + Ifges(3,6) * t95 + t145 * t41 + m(6) * (-t42 * t6 + t43 * t8) + m(5) * (t20 * t56 + t21 * t57) + t101 + Ifges(4,6) * t60 + Ifges(4,5) * t61 + t50 * mrSges(4,1) - t51 * mrSges(4,2) + t40 * t14 + (-t91 * mrSges(3,1) - t95 * mrSges(3,2)) * pkin(7) + (-t43 * t30 - t42 * t31) * mrSges(6,3) + (t57 * t48 - t56 * t49) * mrSges(5,3) + (m(4) * (t50 * t87 + t51 * t86) + (t60 * t86 - t61 * t87) * mrSges(4,3)) * pkin(2); 0.2e1 * t124 - 0.2e1 * t123 - 0.2e1 * t125 + Ifges(3,3) + Ifges(4,3) + 0.2e1 * t33 + 0.2e1 * t34 + 0.2e1 * t35 + 0.2e1 * t38 + m(7) * (t115 * t41 ^ 2 + t40 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t56 ^ 2 + t57 ^ 2) + t108 + (0.2e1 * t87 * mrSges(4,1) - 0.2e1 * t86 * mrSges(4,2) + m(4) * (t86 ^ 2 + t87 ^ 2) * pkin(2)) * pkin(2); t30 * mrSges(6,1) + t88 * t15 + t92 * t16 + t26 + m(7) * (t2 * t92 + t3 * t88) + m(6) * t32 + m(5) * t53 + m(4) * t76 + t109 + t110; 0; m(7) * t115 + m(4) + m(5) + m(6); t101 + (m(6) * (-t6 * t93 + t8 * t89) + (-t89 * t30 - t93 * t31) * mrSges(6,3)) * pkin(4) + t143 * t75 + t141 * t74; t103 + t124 - t123 + m(7) * (t41 * t111 + t40 * t75) + (m(6) * (t42 * t93 + t43 * t89) - t120) * pkin(4) + Ifges(5,3) + t142; 0; -0.2e1 * t113 + 0.2e1 * t58 + 0.2e1 * t64 + 0.2e1 * t65 + 0.2e1 * t79 + m(7) * (t115 * t74 ^ 2 + t75 ^ 2) + m(6) * (t89 ^ 2 + t93 ^ 2) * pkin(4) ^ 2 + t108; -pkin(5) * t143 + pkin(10) * t141 + t102; m(7) * (t115 * t41 * pkin(10) - pkin(5) * t40) + t103 + t144; 0; m(7) * (-pkin(5) * t75 + pkin(10) * t111) - t113 + t112 + t142 + t144; -0.2e1 * t133 + m(7) * (t115 * pkin(10) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t78 + 0.2e1 * t77 + t112; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t127 + t117; -t106 * t41 + t116; -t66; -t106 * t74 + t116; -t106 * pkin(10) + t116; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
