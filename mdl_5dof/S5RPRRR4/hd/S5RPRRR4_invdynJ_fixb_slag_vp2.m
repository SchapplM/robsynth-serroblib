% Calculate vector of inverse dynamics joint torques for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:19
% EndTime: 2019-12-05 18:14:24
% DurationCPUTime: 1.59s
% Computational Cost: add. (2691->229), mult. (4880->301), div. (0->0), fcn. (2672->16), ass. (0->117)
t101 = cos(qJ(5));
t97 = sin(qJ(5));
t70 = -mrSges(6,1) * t101 + mrSges(6,2) * t97;
t159 = t97 / 0.2e1;
t161 = m(6) * pkin(8);
t166 = mrSges(6,3) + t161;
t139 = qJD(5) * t97;
t92 = qJDD(1) + qJDD(3);
t87 = qJDD(4) + t92;
t93 = qJD(1) + qJD(3);
t90 = qJD(4) + t93;
t55 = t101 * t87 - t139 * t90;
t138 = qJD(5) * t101;
t56 = t138 * t90 + t87 * t97;
t33 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t102 = cos(qJ(4));
t103 = cos(qJ(3));
t95 = sin(pkin(9));
t157 = pkin(1) * t95;
t134 = qJD(1) * t157;
t96 = cos(pkin(9));
t82 = pkin(1) * t96 + pkin(2);
t69 = t82 * qJD(1);
t99 = sin(qJ(3));
t47 = t103 * t134 + t69 * t99;
t141 = t102 * t47;
t46 = t103 * t69 - t134 * t99;
t39 = pkin(3) * t93 + t46;
t98 = sin(qJ(4));
t18 = t39 * t98 + t141;
t136 = qJD(1) * qJD(3);
t140 = qJD(3) * t69;
t68 = t82 * qJDD(1);
t28 = -t99 * t140 + t103 * t68 + (-qJDD(1) * t99 - t103 * t136) * t157;
t20 = pkin(3) * t92 + t28;
t27 = t103 * t140 + t68 * t99 + (qJDD(1) * t103 - t136 * t99) * t157;
t9 = -qJD(4) * t18 + t102 * t20 - t27 * t98;
t6 = -pkin(4) * t87 - t9;
t165 = m(6) * t6 + t33;
t164 = m(6) * pkin(4) + mrSges(5,1) - t70;
t149 = t47 * t98;
t17 = t102 * t39 - t149;
t128 = mrSges(6,1) * t97 + mrSges(6,2) * t101;
t15 = -pkin(4) * t90 - t17;
t163 = t15 * t128 + qJD(5) * (Ifges(6,5) * t101 - Ifges(6,6) * t97) / 0.2e1;
t162 = t28 * mrSges(4,1) + Ifges(4,3) * t92;
t58 = t103 * t82 - t99 * t157;
t57 = pkin(3) + t58;
t59 = t103 * t157 + t82 * t99;
t32 = t102 * t59 + t98 * t57;
t16 = pkin(8) * t90 + t18;
t13 = qJD(2) * t101 - t16 * t97;
t14 = qJD(2) * t97 + t101 * t16;
t125 = t101 * t14 - t13 * t97;
t158 = m(5) + m(6);
t94 = qJ(1) + pkin(9);
t91 = qJ(3) + t94;
t83 = qJ(4) + t91;
t77 = sin(t83);
t156 = g(2) * t77;
t8 = t17 * qJD(4) + t102 * t27 + t98 * t20;
t5 = pkin(8) * t87 + t8;
t3 = -qJD(5) * t14 + qJDD(2) * t101 - t5 * t97;
t155 = t3 * t97;
t152 = Ifges(6,4) * t97;
t2 = qJD(5) * t13 + qJDD(2) * t97 + t101 * t5;
t151 = t101 * t2;
t148 = t90 * t97;
t146 = Ifges(6,4) * t101;
t145 = Ifges(6,2) * t101;
t143 = t101 * t90;
t137 = m(4) + t158;
t133 = m(3) + t137;
t54 = t70 * t90;
t132 = t90 * mrSges(5,1) - t54;
t130 = pkin(3) * t158 + mrSges(4,1);
t80 = sin(t91);
t81 = cos(t91);
t129 = g(2) * t81 + g(3) * t80;
t127 = t145 + t152;
t124 = t101 * t13 + t14 * t97;
t41 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t55;
t42 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t56;
t123 = t101 * t41 - t97 * t42;
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t148;
t63 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t143;
t122 = t101 * t63 - t97 * t62;
t31 = t102 * t57 - t59 * t98;
t120 = pkin(2) * t137 + mrSges(3,1);
t118 = pkin(1) * t133 + mrSges(2,1);
t116 = t97 * (Ifges(6,1) * t101 - t152);
t115 = t90 * mrSges(5,2) - t122;
t78 = cos(t83);
t114 = -t77 * mrSges(5,2) + t164 * t78;
t113 = -qJD(5) * t124 - t155;
t112 = -t80 * mrSges(4,2) + t114;
t111 = -t13 * t138 - t139 * t14 + t151 - t155;
t44 = Ifges(6,6) * qJD(5) + t127 * t90;
t72 = Ifges(6,4) * t143;
t45 = Ifges(6,1) * t148 + Ifges(6,5) * qJD(5) + t72;
t110 = t9 * mrSges(5,1) + mrSges(6,3) * t151 + t6 * t70 + (Ifges(6,1) * t56 + Ifges(6,4) * t55) * t159 + t101 * (Ifges(6,4) * t56 + Ifges(6,2) * t55) / 0.2e1 + t55 * t127 / 0.2e1 + t56 * (Ifges(6,1) * t97 + t146) / 0.2e1 - t44 * t139 / 0.2e1 + Ifges(5,3) * t87 + (t45 + t90 * (-Ifges(6,2) * t97 + t146)) * t138 / 0.2e1 + (0.2e1 * Ifges(6,5) * t159 + Ifges(6,6) * t101) * qJDD(5) + (t116 * t90 / 0.2e1 + t163) * qJD(5);
t109 = t164 * t77 + (mrSges(5,2) - t166) * t78;
t108 = t81 * mrSges(4,2) + t109;
t107 = m(6) * t111 - t138 * t62 - t139 * t63 + t123;
t106 = (t113 + t156) * mrSges(6,3) + t110 - t8 * mrSges(5,2);
t104 = cos(qJ(1));
t100 = sin(qJ(1));
t89 = cos(t94);
t88 = sin(t94);
t53 = t59 * qJD(3);
t52 = t58 * qJD(3);
t29 = -pkin(4) - t31;
t22 = t102 * t46 - t149;
t21 = t46 * t98 + t141;
t11 = t32 * qJD(4) + t102 * t53 + t52 * t98;
t10 = qJD(4) * t31 + t102 * t52 - t53 * t98;
t1 = [(-t10 * t90 - t32 * t87 - t8) * mrSges(5,2) + (-t52 * t93 - t59 * t92 - t27) * mrSges(4,2) + (-t11 * t90 + t31 * t87) * mrSges(5,1) + (-t53 * t93 + t58 * t92) * mrSges(4,1) + t113 * mrSges(6,3) + t122 * t10 + t107 * (pkin(8) + t32) + (mrSges(2,2) * t104 + mrSges(3,2) * t89 + t100 * t118 + t120 * t88 + t130 * t80 + t108) * g(3) + t110 + t11 * t54 + t29 * t33 + m(6) * (t125 * t10 + t11 * t15 + t29 * t6) + (-mrSges(2,2) * t100 - mrSges(3,2) * t88 + t118 * t104 + t120 * t89 + t130 * t81 + t166 * t77 + t112) * g(2) + m(4) * (t27 * t59 + t28 * t58 - t46 * t53 + t47 * t52) + m(5) * (t10 * t18 - t11 * t17 + t31 * t9 + t32 * t8) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t96 - 0.2e1 * mrSges(3,2) * t95 + m(3) * (t95 ^ 2 + t96 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t162; t63 * t138 - t62 * t139 + t97 * t41 + t101 * t42 + m(6) * (qJD(5) * t125 + t101 * t3 + t2 * t97) + (m(3) + m(4) + m(5)) * qJDD(2) - t133 * g(1); -m(6) * (t125 * t22 + t15 * t21) + (t81 * mrSges(4,1) + t161 * t77 + t112) * g(2) + (t80 * mrSges(4,1) + t108) * g(3) + (t46 * t93 - t27) * mrSges(4,2) - m(5) * (-t17 * t21 + t18 * t22) + t47 * t93 * mrSges(4,1) + t106 + t115 * t22 + t132 * t21 + ((mrSges(5,1) * t102 - mrSges(5,2) * t98) * t87 + t129 * m(6) + (t102 * t9 + t8 * t98 + t129) * m(5) + ((-m(5) * t17 + m(6) * t15 - t132) * t98 + (m(5) * t18 + m(6) * t125 - t115) * t102) * qJD(4)) * pkin(3) + t107 * (pkin(3) * t98 + pkin(8)) + t162 + t165 * (-pkin(3) * t102 - pkin(4)); -m(6) * (t125 * t17 + t15 * t18) + t114 * g(2) + t109 * g(3) + t132 * t18 + t115 * t17 + t106 + ((-t101 * t62 - t97 * t63) * qJD(5) + (t111 + t156) * m(6) + t123) * pkin(8) - t165 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t56 + Ifges(6,6) * t55 + Ifges(6,3) * qJDD(5) + g(1) * t70 - t13 * t63 + t14 * t62 + (t44 * t159 + (-t116 / 0.2e1 + t145 * t159) * t90 + t124 * mrSges(6,3) - (t45 + t72) * t101 / 0.2e1 - t163) * t90 + (g(3) * t78 - t156) * t128;];
tau = t1;
