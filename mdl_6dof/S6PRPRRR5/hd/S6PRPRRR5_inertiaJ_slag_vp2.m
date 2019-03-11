% Calculate joint inertia matrix for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:41:07
% EndTime: 2019-03-08 20:41:09
% DurationCPUTime: 0.72s
% Computational Cost: add. (899->196), mult. (1805->284), div. (0->0), fcn. (1803->10), ass. (0->88)
t72 = sin(qJ(6));
t76 = cos(qJ(6));
t102 = t72 ^ 2 + t76 ^ 2;
t129 = mrSges(7,3) * t102;
t50 = -t76 * mrSges(7,1) + t72 * mrSges(7,2);
t128 = t50 - mrSges(6,1);
t108 = t72 * mrSges(7,3);
t73 = sin(qJ(5));
t74 = sin(qJ(4));
t77 = cos(qJ(5));
t78 = cos(qJ(4));
t46 = t73 * t74 - t77 * t78;
t47 = t73 * t78 + t77 * t74;
t19 = -t47 * mrSges(7,2) + t46 * t108;
t111 = t46 * t76;
t20 = t47 * mrSges(7,1) + mrSges(7,3) * t111;
t127 = t76 * t19 - t72 * t20;
t92 = -mrSges(7,1) * t72 - mrSges(7,2) * t76;
t17 = t92 * t46;
t80 = -pkin(2) - pkin(8);
t120 = -pkin(9) + t80;
t49 = t120 * t74;
t98 = t120 * t78;
t23 = t73 * t49 - t77 * t98;
t126 = m(7) * t23 + t17;
t57 = t74 * pkin(4) + qJ(3);
t18 = t47 * pkin(5) + t46 * pkin(10) + t57;
t25 = t77 * t49 + t73 * t98;
t2 = t76 * t18 - t72 * t25;
t3 = t72 * t18 + t76 * t25;
t94 = -t2 * t72 + t3 * t76;
t125 = m(7) * t94 + t127;
t70 = sin(pkin(6));
t79 = cos(qJ(2));
t109 = t70 * t79;
t71 = cos(pkin(6));
t35 = -t78 * t109 - t71 * t74;
t36 = -t74 * t109 + t71 * t78;
t13 = -t77 * t35 + t73 * t36;
t124 = t13 ^ 2;
t123 = t23 ^ 2;
t43 = t46 ^ 2;
t68 = t78 ^ 2;
t122 = -2 * mrSges(6,3);
t121 = m(6) * pkin(4);
t118 = mrSges(7,3) * t76;
t117 = Ifges(7,4) * t72;
t116 = Ifges(7,4) * t76;
t115 = t23 * t13;
t114 = t46 * t13;
t113 = t46 * t23;
t112 = t46 * t72;
t75 = sin(qJ(2));
t110 = t70 * t75;
t105 = -Ifges(7,5) * t111 + Ifges(7,3) * t47;
t104 = t74 * mrSges(5,1) + t78 * mrSges(5,2) + mrSges(4,3);
t103 = Ifges(7,5) * t72 + Ifges(7,6) * t76;
t101 = t74 ^ 2 + t68;
t52 = Ifges(7,2) * t76 + t117;
t53 = Ifges(7,1) * t72 + t116;
t100 = t76 * t52 + t72 * t53 + Ifges(6,3);
t99 = m(5) * t101;
t97 = t102 * pkin(10);
t59 = t73 * pkin(4) + pkin(10);
t96 = t102 * t59;
t95 = t101 * mrSges(5,3);
t15 = t73 * t35 + t77 * t36;
t7 = t76 * t110 - t72 * t15;
t8 = t72 * t110 + t76 * t15;
t93 = -t7 * t72 + t76 * t8;
t91 = t47 * t15 + t114;
t90 = t78 * t35 + t74 * t36;
t89 = t46 * t77 - t47 * t73;
t88 = 0.2e1 * t129;
t87 = t128 * t46 + (-mrSges(6,2) + t129) * t47;
t86 = (t77 * mrSges(6,1) - t73 * mrSges(6,2)) * pkin(4);
t85 = -t15 * mrSges(6,2) - t7 * t108 + t8 * t118 + t128 * t13;
t10 = Ifges(7,6) * t47 + (Ifges(7,2) * t72 - t116) * t46;
t11 = Ifges(7,5) * t47 + (-Ifges(7,1) * t76 + t117) * t46;
t84 = -t25 * mrSges(6,2) - t2 * t108 + t3 * t118 + t52 * t112 / 0.2e1 - t53 * t111 / 0.2e1 - Ifges(6,5) * t46 + t72 * t11 / 0.2e1 + t76 * t10 / 0.2e1 + (t103 / 0.2e1 - Ifges(6,6)) * t47 + t128 * t23;
t81 = qJ(3) ^ 2;
t64 = t70 ^ 2;
t60 = -t77 * pkin(4) - pkin(5);
t56 = t64 * t75 ^ 2;
t54 = qJ(3) * t110;
t42 = t47 ^ 2;
t22 = t47 * mrSges(6,1) - t46 * mrSges(6,2);
t1 = [m(2) + m(7) * (t7 ^ 2 + t8 ^ 2 + t124) + m(6) * (t15 ^ 2 + t124 + t56) + m(5) * (t35 ^ 2 + t36 ^ 2 + t56) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * (t64 * t79 ^ 2 + t71 ^ 2 + t56); t13 * t17 + t8 * t19 + t7 * t20 - t91 * mrSges(6,3) - t90 * mrSges(5,3) + ((mrSges(3,1) - mrSges(4,2)) * t79 + (-mrSges(3,2) + t22 + t104) * t75) * t70 + m(7) * (t2 * t7 + t3 * t8 + t115) + m(6) * (t57 * t110 + t25 * t15 + t115) + m(5) * (t90 * t80 + t54) + m(4) * (pkin(2) * t109 + t54); Ifges(5,1) * t68 - 0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t23 * t17 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t57 * t22 + Ifges(4,1) + Ifges(3,3) + (-0.2e1 * Ifges(5,4) * t78 + Ifges(5,2) * t74) * t74 + (Ifges(6,2) * t47 + t25 * t122 + t105) * t47 + (t23 * t122 + Ifges(6,1) * t46 + t72 * t10 - t76 * t11 + (Ifges(7,6) * t72 + (2 * Ifges(6,4))) * t47) * t46 + m(7) * (t2 ^ 2 + t3 ^ 2 + t123) + m(6) * (t25 ^ 2 + t57 ^ 2 + t123) + m(5) * (t101 * t80 ^ 2 + t81) + m(4) * (pkin(2) ^ 2 + t81) + 0.2e1 * t104 * qJ(3) - 0.2e1 * t80 * t95; -m(4) * t109 + m(7) * (t93 * t47 + t114) + m(6) * t91 + m(5) * t90; -m(4) * pkin(2) - t43 * mrSges(6,3) + t46 * t17 + mrSges(4,2) - t95 + (-mrSges(6,3) * t47 + t127) * t47 + m(7) * (t94 * t47 + t113) + m(6) * (t47 * t25 + t113) + t80 * t99; m(4) + t99 + m(6) * (t42 + t43) + m(7) * (t102 * t42 + t43); t35 * mrSges(5,1) - t36 * mrSges(5,2) + m(7) * (t60 * t13 + t93 * t59) + (-t13 * t77 + t15 * t73) * t121 + t85; (m(6) * (-t23 * t77 + t25 * t73) + t89 * mrSges(6,3)) * pkin(4) + (t80 * mrSges(5,1) + Ifges(5,5)) * t78 + (-t80 * mrSges(5,2) - Ifges(5,6)) * t74 + t84 + t126 * t60 + t125 * t59; t78 * mrSges(5,1) - t74 * mrSges(5,2) + m(7) * (t60 * t46 + t47 * t96) - t89 * t121 + t87; 0.2e1 * t60 * t50 + Ifges(5,3) + 0.2e1 * t86 + t59 * t88 + m(7) * (t102 * t59 ^ 2 + t60 ^ 2) + m(6) * (t73 ^ 2 + t77 ^ 2) * pkin(4) ^ 2 + t100; m(7) * (-pkin(5) * t13 + t93 * pkin(10)) + t85; -t126 * pkin(5) + t125 * pkin(10) + t84; m(7) * (-pkin(5) * t46 + t47 * t97) + t87; m(7) * (-pkin(5) * t60 + pkin(10) * t96) + (-pkin(5) + t60) * t50 + t86 + (t96 + t97) * mrSges(7,3) + t100; -0.2e1 * pkin(5) * t50 + m(7) * (t102 * pkin(10) ^ 2 + pkin(5) ^ 2) + pkin(10) * t88 + t100; t7 * mrSges(7,1) - t8 * mrSges(7,2); t2 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,6) * t112 + t105; t92 * t47; t92 * t59 + t103; t92 * pkin(10) + t103; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
