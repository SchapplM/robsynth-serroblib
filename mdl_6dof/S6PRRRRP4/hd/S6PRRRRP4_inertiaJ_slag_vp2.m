% Calculate joint inertia matrix for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:10
% EndTime: 2019-03-09 00:13:13
% DurationCPUTime: 1.17s
% Computational Cost: add. (1307->299), mult. (2833->411), div. (0->0), fcn. (2839->10), ass. (0->111)
t144 = -mrSges(6,1) - mrSges(7,1);
t143 = -mrSges(6,2) + mrSges(7,3);
t142 = 2 * pkin(8);
t141 = mrSges(7,2) + mrSges(6,3);
t102 = sin(qJ(3));
t106 = cos(qJ(3));
t103 = sin(qJ(2));
t98 = sin(pkin(6));
t128 = t103 * t98;
t99 = cos(pkin(6));
t57 = t102 * t128 - t99 * t106;
t56 = t57 ^ 2;
t140 = 2 * mrSges(7,3);
t139 = -pkin(10) - pkin(9);
t138 = pkin(8) * t106;
t92 = t102 * pkin(8);
t136 = Ifges(7,2) + Ifges(6,3);
t100 = sin(qJ(5));
t104 = cos(qJ(5));
t101 = sin(qJ(4));
t105 = cos(qJ(4));
t124 = t102 * t105;
t73 = -pkin(3) * t106 - pkin(9) * t102 - pkin(2);
t66 = t105 * t73;
t24 = -pkin(10) * t124 + t66 + (-pkin(8) * t101 - pkin(4)) * t106;
t125 = t101 * t102;
t46 = t101 * t73 + t105 * t138;
t34 = -pkin(10) * t125 + t46;
t7 = t100 * t24 + t104 * t34;
t68 = t100 * t105 + t101 * t104;
t54 = t68 * t102;
t41 = -mrSges(7,2) * t54 - mrSges(7,3) * t106;
t42 = mrSges(6,2) * t106 - mrSges(6,3) * t54;
t135 = t41 + t42;
t67 = t100 * t101 - t104 * t105;
t55 = t67 * t102;
t43 = -mrSges(6,1) * t106 + mrSges(6,3) * t55;
t44 = t106 * mrSges(7,1) - t55 * mrSges(7,2);
t134 = -t43 + t44;
t74 = -mrSges(5,1) * t105 + mrSges(5,2) * t101;
t133 = t74 - mrSges(4,1);
t72 = pkin(4) * t125 + t92;
t132 = t101 ^ 2 + t105 ^ 2;
t131 = Ifges(5,4) * t101;
t130 = Ifges(5,4) * t105;
t129 = t102 * t57;
t59 = t102 * t99 + t106 * t128;
t127 = t106 * t59;
t107 = cos(qJ(2));
t126 = t107 * t98;
t121 = t139 * t101;
t78 = t139 * t105;
t38 = -t100 * t78 - t104 * t121;
t40 = t100 * t121 - t104 * t78;
t123 = t38 ^ 2 + t40 ^ 2;
t122 = Ifges(5,3) + t136;
t87 = -pkin(4) * t105 - pkin(3);
t32 = -t101 * t59 - t105 * t126;
t33 = -t101 * t126 + t105 * t59;
t11 = t100 * t33 - t104 * t32;
t13 = t100 * t32 + t104 * t33;
t120 = t11 * t38 + t40 * t13;
t118 = (Ifges(7,4) + Ifges(6,5)) * t55 + (Ifges(6,6) - Ifges(7,6)) * t54;
t116 = mrSges(5,1) * t101 + mrSges(5,2) * t105;
t6 = -t100 * t34 + t104 * t24;
t115 = -t101 * t32 + t105 * t33;
t114 = t144 * t11 + t143 * t13;
t113 = (mrSges(6,1) * t104 - mrSges(6,2) * t100) * pkin(4);
t3 = -qJ(6) * t106 + t7;
t4 = pkin(5) * t106 - t6;
t112 = t6 * mrSges(6,1) - t4 * mrSges(7,1) - t7 * mrSges(6,2) + t3 * mrSges(7,3) - t118;
t61 = Ifges(7,6) * t67;
t62 = Ifges(6,6) * t67;
t63 = Ifges(6,5) * t68;
t64 = Ifges(7,4) * t68;
t111 = t143 * t40 + t144 * t38 + t61 - t62 + t63 + t64;
t109 = pkin(8) ^ 2;
t97 = t106 ^ 2;
t95 = t102 ^ 2;
t93 = t98 ^ 2;
t91 = t95 * t109;
t90 = Ifges(5,5) * t101;
t89 = Ifges(5,6) * t105;
t86 = -pkin(4) * t104 - pkin(5);
t84 = pkin(4) * t100 + qJ(6);
t82 = t93 * t107 ^ 2;
t79 = Ifges(5,5) * t124;
t77 = Ifges(5,1) * t101 + t130;
t76 = Ifges(5,2) * t105 + t131;
t75 = -mrSges(4,1) * t106 + mrSges(4,2) * t102;
t71 = -mrSges(5,1) * t106 - mrSges(5,3) * t124;
t70 = mrSges(5,2) * t106 - mrSges(5,3) * t125;
t60 = t116 * t102;
t53 = -Ifges(5,5) * t106 + (Ifges(5,1) * t105 - t131) * t102;
t52 = -Ifges(5,6) * t106 + (-Ifges(5,2) * t101 + t130) * t102;
t45 = -t101 * t138 + t66;
t31 = Ifges(6,1) * t68 - Ifges(6,4) * t67;
t30 = Ifges(7,1) * t68 + Ifges(7,5) * t67;
t29 = Ifges(6,4) * t68 - Ifges(6,2) * t67;
t28 = Ifges(7,5) * t68 + Ifges(7,3) * t67;
t27 = mrSges(6,1) * t67 + mrSges(6,2) * t68;
t26 = mrSges(7,1) * t67 - mrSges(7,3) * t68;
t22 = pkin(5) * t67 - qJ(6) * t68 + t87;
t20 = mrSges(6,1) * t54 - mrSges(6,2) * t55;
t19 = mrSges(7,1) * t54 + mrSges(7,3) * t55;
t18 = -Ifges(6,1) * t55 - Ifges(6,4) * t54 - Ifges(6,5) * t106;
t17 = -Ifges(7,1) * t55 - Ifges(7,4) * t106 + Ifges(7,5) * t54;
t16 = -Ifges(6,4) * t55 - Ifges(6,2) * t54 - Ifges(6,6) * t106;
t15 = -Ifges(7,5) * t55 - Ifges(7,6) * t106 + Ifges(7,3) * t54;
t14 = pkin(5) * t54 + qJ(6) * t55 + t72;
t1 = [m(2) + m(5) * (t32 ^ 2 + t33 ^ 2 + t56) + m(4) * (t59 ^ 2 + t56 + t82) + m(3) * (t103 ^ 2 * t93 + t99 ^ 2 + t82) + (m(7) + m(6)) * (t11 ^ 2 + t13 ^ 2 + t56); mrSges(4,3) * t127 + t32 * t71 + t33 * t70 + t135 * t13 + t134 * t11 + (-t103 * mrSges(3,2) + (mrSges(3,1) - t75) * t107) * t98 + (t102 * mrSges(4,3) + t19 + t20 + t60) * t57 + m(7) * (t11 * t4 + t13 * t3 + t14 * t57) + m(6) * (-t11 * t6 + t13 * t7 + t57 * t72) + m(5) * (pkin(8) * t129 + t32 * t45 + t33 * t46) + m(4) * (pkin(2) * t126 + (t127 + t129) * pkin(8)); -0.2e1 * pkin(2) * t75 + 0.2e1 * t14 * t19 + 0.2e1 * t72 * t20 + 0.2e1 * t3 * t41 + 0.2e1 * t4 * t44 + 0.2e1 * t7 * t42 + 0.2e1 * t6 * t43 + 0.2e1 * t45 * t71 + 0.2e1 * t46 * t70 + Ifges(3,3) - (t17 + t18) * t55 + (t15 - t16) * t54 + (t95 + t97) * mrSges(4,3) * t142 + (Ifges(4,1) * t102 - t101 * t52 + t105 * t53 + t142 * t60) * t102 + m(4) * (pkin(2) ^ 2 + t109 * t97 + t91) + m(5) * (t45 ^ 2 + t46 ^ 2 + t91) + m(6) * (t6 ^ 2 + t7 ^ 2 + t72 ^ 2) + m(7) * (t14 ^ 2 + t3 ^ 2 + t4 ^ 2) + (-t79 + (Ifges(4,2) + t122) * t106 + (Ifges(5,6) * t101 + (2 * Ifges(4,4))) * t102 + t118) * t106; -t59 * mrSges(4,2) + t115 * mrSges(5,3) + (t26 + t27 + t133) * t57 + m(7) * (t22 * t57 + t120) + m(6) * (t57 * t87 + t120) + m(5) * (-pkin(3) * t57 + pkin(9) * t115) + t141 * (t11 * t68 - t13 * t67); -pkin(3) * t60 + t14 * t26 + t22 * t19 + t87 * t20 + t72 * t27 - (t30 / 0.2e1 + t31 / 0.2e1) * t55 + (t28 / 0.2e1 - t29 / 0.2e1) * t54 + t135 * t40 + t134 * t38 + (-t90 / 0.2e1 - t89 / 0.2e1 - t63 / 0.2e1 + t62 / 0.2e1 - t64 / 0.2e1 - t61 / 0.2e1 + Ifges(4,6) - pkin(8) * mrSges(4,2)) * t106 + (t52 / 0.2e1 + pkin(9) * t70 + t46 * mrSges(5,3)) * t105 + (t53 / 0.2e1 - pkin(9) * t71 - t45 * mrSges(5,3)) * t101 + (Ifges(4,5) - t101 * t76 / 0.2e1 + t105 * t77 / 0.2e1 + t133 * pkin(8)) * t102 + m(5) * (-pkin(3) * t92 + (-t101 * t45 + t105 * t46) * pkin(9)) + m(6) * (-t38 * t6 + t40 * t7 + t72 * t87) + m(7) * (t14 * t22 + t3 * t40 + t38 * t4) + (t17 / 0.2e1 + t18 / 0.2e1 + t4 * mrSges(7,2) - t6 * mrSges(6,3)) * t68 + (t15 / 0.2e1 - t16 / 0.2e1 - t3 * mrSges(7,2) - t7 * mrSges(6,3)) * t67; -0.2e1 * pkin(3) * t74 + t101 * t77 + t105 * t76 + 0.2e1 * t22 * t26 + 0.2e1 * t87 * t27 + Ifges(4,3) + m(7) * (t22 ^ 2 + t123) + m(6) * (t87 ^ 2 + t123) + m(5) * (pkin(9) ^ 2 * t132 + pkin(3) ^ 2) + (t30 + t31) * t68 + (t28 - t29) * t67 + 0.2e1 * t132 * pkin(9) * mrSges(5,3) + 0.2e1 * (t38 * t68 - t40 * t67) * t141; t32 * mrSges(5,1) - t33 * mrSges(5,2) + m(7) * (t11 * t86 + t13 * t84) + m(6) * (t100 * t13 - t104 * t11) * pkin(4) + t114; t112 + (m(6) * (t100 * t7 + t104 * t6) + t104 * t43 + t100 * t42) * pkin(4) + m(7) * (t3 * t84 + t4 * t86) - t122 * t106 + t79 + t84 * t41 + t86 * t44 + t45 * mrSges(5,1) - t46 * mrSges(5,2) - Ifges(5,6) * t125; m(7) * (t38 * t86 + t40 * t84) + t90 + t89 - t116 * pkin(9) + (-t67 * t84 + t68 * t86) * mrSges(7,2) + (m(6) * (t100 * t40 - t104 * t38) + (-t100 * t67 - t104 * t68) * mrSges(6,3)) * pkin(4) + t111; -0.2e1 * t86 * mrSges(7,1) + t84 * t140 + 0.2e1 * t113 + m(7) * (t84 ^ 2 + t86 ^ 2) + m(6) * (t100 ^ 2 + t104 ^ 2) * pkin(4) ^ 2 + t122; m(7) * (-pkin(5) * t11 + qJ(6) * t13) + t114; m(7) * (-pkin(5) * t4 + qJ(6) * t3) + qJ(6) * t41 - pkin(5) * t44 - t136 * t106 + t112; m(7) * (-pkin(5) * t38 + qJ(6) * t40) + (-pkin(5) * t68 - qJ(6) * t67) * mrSges(7,2) + t111; m(7) * (-pkin(5) * t86 + qJ(6) * t84) + t113 + (qJ(6) + t84) * mrSges(7,3) + (pkin(5) - t86) * mrSges(7,1) + t136; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t140 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t136; m(7) * t11; m(7) * t4 + t44; m(7) * t38 + t68 * mrSges(7,2); m(7) * t86 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
