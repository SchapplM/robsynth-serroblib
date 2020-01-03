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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:52:03
% EndTime: 2020-01-03 11:52:08
% DurationCPUTime: 1.94s
% Computational Cost: add. (2691->247), mult. (4880->323), div. (0->0), fcn. (2672->16), ass. (0->121)
t109 = sin(qJ(5));
t168 = t109 / 0.2e1;
t105 = qJD(1) + qJD(3);
t100 = qJD(4) + t105;
t113 = cos(qJ(5));
t146 = qJD(5) * t109;
t104 = qJDD(1) + qJDD(3);
t97 = qJDD(4) + t104;
t55 = -t100 * t146 + t113 * t97;
t145 = qJD(5) * t113;
t56 = t100 * t145 + t109 * t97;
t33 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t110 = sin(qJ(4));
t114 = cos(qJ(4));
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t107 = sin(pkin(9));
t166 = pkin(1) * t107;
t140 = qJD(1) * t166;
t108 = cos(pkin(9));
t92 = pkin(1) * t108 + pkin(2);
t69 = t92 * qJD(1);
t47 = t111 * t69 + t115 * t140;
t153 = t114 * t47;
t46 = -t111 * t140 + t115 * t69;
t39 = pkin(3) * t105 + t46;
t18 = t110 * t39 + t153;
t143 = qJD(1) * qJD(3);
t150 = qJD(3) * t69;
t68 = t92 * qJDD(1);
t28 = -t111 * t150 + t115 * t68 + (-qJDD(1) * t111 - t115 * t143) * t166;
t20 = pkin(3) * t104 + t28;
t27 = t115 * t150 + t111 * t68 + (qJDD(1) * t115 - t111 * t143) * t166;
t9 = -t18 * qJD(4) - t110 * t27 + t114 * t20;
t6 = -pkin(4) * t97 - t9;
t177 = m(6) * t6 + t33;
t156 = t110 * t47;
t17 = t114 * t39 - t156;
t70 = -mrSges(6,1) * t113 + t109 * mrSges(6,2);
t176 = -mrSges(5,1) + t70;
t134 = mrSges(6,1) * t109 + mrSges(6,2) * t113;
t15 = -pkin(4) * t100 - t17;
t175 = t15 * t134 + qJD(5) * (Ifges(6,5) * t113 - Ifges(6,6) * t109) / 0.2e1;
t174 = m(3) * pkin(1);
t173 = t28 * mrSges(4,1) + Ifges(4,3) * t104;
t16 = pkin(8) * t100 + t18;
t14 = qJD(2) * t109 + t113 * t16;
t147 = t14 * qJD(5);
t8 = t17 * qJD(4) + t110 * t20 + t114 * t27;
t5 = pkin(8) * t97 + t8;
t3 = qJDD(2) * t113 - t109 * t5 - t147;
t172 = -t3 - t147;
t58 = -t111 * t166 + t115 * t92;
t13 = qJD(2) * t113 - t109 * t16;
t2 = t13 * qJD(5) + qJDD(2) * t109 + t113 * t5;
t165 = t113 * t2;
t171 = -t13 * t145 + t165;
t155 = t113 * t14;
t130 = -t109 * t13 + t155;
t106 = qJ(1) + pkin(9);
t101 = qJ(3) + t106;
t93 = qJ(4) + t101;
t86 = cos(t93);
t169 = g(3) * t86;
t164 = t3 * t109;
t57 = pkin(3) + t58;
t59 = t111 * t92 + t115 * t166;
t32 = t110 * t57 + t114 * t59;
t85 = sin(t93);
t163 = t86 * pkin(4) + t85 * pkin(8);
t161 = Ifges(6,4) * t109;
t160 = Ifges(6,4) * t113;
t159 = Ifges(6,2) * t113;
t112 = sin(qJ(1));
t98 = sin(t106);
t152 = t112 * pkin(1) + pkin(2) * t98;
t116 = cos(qJ(1));
t99 = cos(t106);
t151 = t116 * pkin(1) + pkin(2) * t99;
t149 = t100 * t109;
t148 = t100 * t113;
t144 = m(3) + m(4) + m(5);
t91 = cos(t101);
t84 = pkin(3) * t91;
t141 = t84 + t151;
t139 = -mrSges(2,1) - t174;
t54 = t70 * t100;
t136 = t100 * mrSges(5,1) - t54;
t78 = t85 * pkin(4);
t90 = sin(t101);
t83 = pkin(3) * t90;
t135 = -pkin(8) * t86 + t78 + t83;
t133 = t159 + t161;
t131 = t109 * t14 + t113 * t13;
t41 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t55;
t42 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t56;
t129 = -t109 * t42 + t113 * t41;
t31 = -t110 * t59 + t114 * t57;
t127 = -t86 * mrSges(5,2) + t176 * t85;
t125 = t109 * (Ifges(6,1) * t113 - t161);
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t149;
t63 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t148;
t124 = t100 * mrSges(5,2) + t109 * t62 - t113 * t63;
t123 = -t90 * mrSges(4,1) - t91 * mrSges(4,2) + t127;
t122 = (mrSges(5,2) - mrSges(6,3)) * t85 + t176 * t86;
t121 = -t91 * mrSges(4,1) + t90 * mrSges(4,2) + t122;
t120 = -t14 * t146 - t164 + t171;
t44 = Ifges(6,6) * qJD(5) + t133 * t100;
t72 = Ifges(6,4) * t148;
t45 = Ifges(6,1) * t149 + Ifges(6,5) * qJD(5) + t72;
t119 = t9 * mrSges(5,1) + mrSges(6,3) * t165 + t6 * t70 + (Ifges(6,1) * t56 + Ifges(6,4) * t55) * t168 + t113 * (Ifges(6,4) * t56 + Ifges(6,2) * t55) / 0.2e1 + t55 * t133 / 0.2e1 + t56 * (Ifges(6,1) * t109 + t160) / 0.2e1 - t44 * t146 / 0.2e1 + Ifges(5,3) * t97 + (t45 + t100 * (-Ifges(6,2) * t109 + t160)) * t145 / 0.2e1 + (0.2e1 * Ifges(6,5) * t168 + Ifges(6,6) * t113) * qJDD(5) + (t125 * t100 / 0.2e1 + t175) * qJD(5);
t118 = (-t131 * qJD(5) - t164 + t169) * mrSges(6,3) - t8 * mrSges(5,2) + t119;
t53 = t59 * qJD(3);
t52 = t58 * qJD(3);
t30 = pkin(8) + t32;
t29 = -pkin(4) - t31;
t22 = t114 * t46 - t156;
t21 = t110 * t46 + t153;
t11 = t32 * qJD(4) + t110 * t52 + t114 * t53;
t10 = t31 * qJD(4) - t110 * t53 + t114 * t52;
t1 = [(t172 * mrSges(6,3) + (-m(6) * t13 - t62) * t10 + (m(6) * t172 - qJD(5) * t63 - t42) * t30) * t109 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t108 * mrSges(3,1) - 0.2e1 * t107 * mrSges(3,2) + (t107 ^ 2 + t108 ^ 2) * t174) * pkin(1)) * qJDD(1) + t11 * t54 + t29 * t33 + (-t11 * t100 + t31 * t97) * mrSges(5,1) + (t58 * t104 - t53 * t105) * mrSges(4,1) + m(6) * (t10 * t155 + t11 * t15 + t171 * t30 + t29 * t6) + (mrSges(2,2) * t112 - mrSges(3,1) * t99 + mrSges(3,2) * t98 - m(6) * (t141 + t163) - m(5) * t141 - m(4) * t151 + t139 * t116 + t121) * g(2) + (-mrSges(2,2) * t116 - m(6) * (t135 + t152) + t86 * mrSges(6,3) - mrSges(3,1) * t98 - mrSges(3,2) * t99 - m(4) * t152 - m(5) * (t83 + t152) + t139 * t112 + t123) * g(3) + (t10 * t63 + t30 * t41 + (-t13 * mrSges(6,3) - t30 * t62) * qJD(5)) * t113 + m(4) * (t27 * t59 + t28 * t58 - t46 * t53 + t47 * t52) + m(5) * (t10 * t18 - t11 * t17 + t31 * t9 + t32 * t8) + t119 + (-t59 * t104 - t52 * t105 - t27) * mrSges(4,2) + (-t10 * t100 - t32 * t97 - t8) * mrSges(5,2) + t173; t63 * t145 - t62 * t146 + t109 * t41 + t113 * t42 + m(6) * (t130 * qJD(5) + t109 * t2 + t113 * t3) + t144 * qJDD(2) + (-m(6) - t144) * g(1); t124 * t22 + t136 * t21 + (m(6) * t120 - t62 * t145 - t63 * t146 + t129) * (pkin(3) * t110 + pkin(8)) + ((mrSges(5,1) * t114 - mrSges(5,2) * t110) * t97 + (-g(2) * t91 - g(3) * t90 + t110 * t8 + t114 * t9) * m(5) + ((-m(5) * t17 + m(6) * t15 - t136) * t110 + (m(5) * t18 + m(6) * t130 - t124) * t114) * qJD(4)) * pkin(3) - m(6) * (t130 * t22 + t15 * t21) + t47 * t105 * mrSges(4,1) + (t105 * t46 - t27) * mrSges(4,2) + (-m(6) * t135 + t123) * g(3) + t118 - m(5) * (-t17 * t21 + t18 * t22) + (-m(6) * (t84 + t163) + t121) * g(2) + t173 + t177 * (-pkin(3) * t114 - pkin(4)); t136 * t18 + t124 * t17 + (-m(6) * t78 + t127) * g(3) + ((-t109 * t63 - t113 * t62) * qJD(5) + (t120 + t169) * m(6) + t129) * pkin(8) + (-m(6) * t163 + t122) * g(2) - m(6) * (t130 * t17 + t15 * t18) + t118 - t177 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t56 + Ifges(6,6) * t55 + Ifges(6,3) * qJDD(5) + g(1) * t70 - t13 * t63 + t14 * t62 + (t44 * t168 + (-t125 / 0.2e1 + t159 * t168) * t100 + t131 * mrSges(6,3) - (t45 + t72) * t113 / 0.2e1 - t175) * t100 + (g(2) * t85 - t169) * t134;];
tau = t1;
