% Calculate vector of inverse dynamics joint torques for
% S5RRPPR4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:44
% DurationCPUTime: 1.79s
% Computational Cost: add. (1907->257), mult. (2472->336), div. (0->0), fcn. (1184->10), ass. (0->121)
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t109 = qJD(1) + qJD(2);
t111 = sin(pkin(8));
t112 = cos(pkin(8));
t117 = cos(qJ(2));
t160 = pkin(1) * qJD(1);
t146 = t117 * t160;
t136 = qJD(3) - t146;
t173 = pkin(2) + pkin(3);
t55 = -t173 * t109 + t136;
t114 = sin(qJ(2));
t147 = t114 * t160;
t73 = qJ(3) * t109 + t147;
t20 = t111 * t55 + t112 * t73;
t18 = -pkin(7) * t109 + t20;
t13 = qJD(4) * t116 - t113 * t18;
t14 = qJD(4) * t113 + t116 * t18;
t130 = -t113 * t13 + t116 * t14;
t167 = mrSges(3,1) + mrSges(4,1);
t181 = mrSges(3,2) - mrSges(4,3);
t180 = mrSges(5,2) - mrSges(6,3);
t79 = mrSges(6,1) * t116 - mrSges(6,2) * t113;
t179 = mrSges(5,1) + t79;
t178 = m(5) + m(6);
t56 = t111 * t146 - t112 * t147;
t177 = qJD(3) * t111 - t56;
t176 = t167 * t114;
t175 = t130 * t112;
t174 = m(3) * pkin(1);
t108 = qJDD(1) + qJDD(2);
t159 = pkin(1) * qJD(2);
t142 = qJD(1) * t159;
t156 = pkin(1) * qJDD(1);
t66 = -t114 * t142 + t117 * t156;
t127 = qJDD(3) - t66;
t27 = -t173 * t108 + t127;
t67 = t114 * t156 + t117 * t142;
t30 = qJ(3) * t108 + qJD(3) * t109 + t67;
t9 = t111 * t27 + t112 * t30;
t6 = -pkin(7) * t108 + t9;
t2 = t13 * qJD(5) + qJDD(4) * t113 + t116 * t6;
t172 = t2 * mrSges(6,3);
t171 = -t113 / 0.2e1;
t170 = -t116 / 0.2e1;
t3 = -t14 * qJD(5) + qJDD(4) * t116 - t113 * t6;
t169 = t113 * t3;
t168 = -qJD(5) / 0.2e1;
t100 = -pkin(1) * t117 - pkin(2);
t89 = -pkin(3) + t100;
t96 = pkin(1) * t114 + qJ(3);
t39 = t111 * t89 + t112 * t96;
t110 = qJ(1) + qJ(2);
t105 = sin(t110);
t106 = cos(t110);
t166 = t106 * pkin(2) + t105 * qJ(3);
t163 = Ifges(6,4) * t113;
t162 = Ifges(6,4) * t116;
t161 = Ifges(6,2) * t116;
t76 = t112 * qJ(3) - t111 * t173;
t155 = t109 * t113;
t154 = t109 * t116;
t152 = qJD(5) * t113;
t151 = qJD(5) * t116;
t150 = m(4) + t178;
t98 = t106 * pkin(3);
t149 = t98 + t166;
t118 = cos(qJ(1));
t107 = t118 * pkin(1);
t148 = t107 + t166;
t145 = t114 * t159;
t141 = -t151 / 0.2e1;
t140 = -mrSges(5,1) * t111 - mrSges(4,3);
t61 = -t108 * t116 + t109 * t152;
t31 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t61;
t72 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t155;
t139 = -qJD(5) * t72 + t31;
t62 = -t108 * t113 - t109 * t151;
t32 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t62;
t74 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t154;
t138 = -qJD(5) * t74 - t32;
t57 = (t111 * t114 + t112 * t117) * t160;
t137 = qJD(3) * t112 - t57;
t135 = mrSges(6,1) * t113 + mrSges(6,2) * t116;
t134 = -t161 - t163;
t133 = -Ifges(6,5) * t116 + Ifges(6,6) * t113;
t19 = -t111 * t73 + t112 * t55;
t132 = -t111 * t19 + t112 * t20;
t8 = -t111 * t30 + t112 * t27;
t38 = -t111 * t96 + t112 * t89;
t131 = -t113 * t14 - t116 * t13;
t129 = -t113 * t72 + t116 * t74;
t75 = -t111 * qJ(3) - t112 * t173;
t58 = -t105 * t111 - t106 * t112;
t59 = -t105 * t112 + t106 * t111;
t128 = -t58 * pkin(4) + pkin(7) * t59 + t149;
t17 = pkin(4) * t109 - t19;
t126 = t17 * t135;
t125 = t113 * (-Ifges(6,1) * t116 + t163);
t124 = t131 * qJD(5) + t116 * t2 - t169;
t123 = t181 * t105 - t167 * t106 + t179 * t58 + t180 * t59;
t122 = m(6) * t124;
t42 = -pkin(2) * t108 + t127;
t45 = Ifges(6,6) * qJD(5) + t134 * t109;
t85 = Ifges(6,4) * t154;
t46 = -Ifges(6,1) * t155 + Ifges(6,5) * qJD(5) - t85;
t5 = pkin(4) * t108 - t8;
t121 = -t42 * mrSges(4,1) - t8 * mrSges(5,1) - t67 * mrSges(3,2) + t9 * mrSges(5,2) + t5 * t79 + (Ifges(6,1) * t62 + Ifges(6,4) * t61 + Ifges(6,5) * qJDD(5)) * t171 + (Ifges(6,4) * t62 + Ifges(6,2) * t61 + Ifges(6,6) * qJDD(5)) * t170 + t61 * t134 / 0.2e1 + t62 * (-Ifges(6,1) * t113 - t162) / 0.2e1 + t30 * mrSges(4,3) + t45 * t152 / 0.2e1 + t46 * t141 + t66 * mrSges(3,1) + qJDD(5) * (-Ifges(6,5) * t113 - Ifges(6,6) * t116) / 0.2e1 + (t125 * t168 + (Ifges(6,2) * t113 - t162) * t141) * t109 + (Ifges(3,3) + Ifges(4,2)) * t108 + (-t126 + qJD(5) * t133 / 0.2e1) * qJD(5) + (t13 * t151 + t14 * t152 + t169) * mrSges(6,3);
t120 = (m(4) * pkin(2) + t178 * t173 + t167) * t105 + (-m(6) * pkin(7) + t180) * t58 + (-m(6) * pkin(4) - t179) * t59 + (-t150 * qJ(3) + t181) * t106;
t115 = sin(qJ(1));
t86 = t117 * t159 + qJD(3);
t71 = -pkin(7) + t76;
t70 = pkin(4) - t75;
t69 = -pkin(2) * t109 + t136;
t60 = t79 * t109;
t41 = t111 * t145 + t112 * t86;
t40 = t111 * t86 - t112 * t145;
t34 = -pkin(7) + t39;
t33 = pkin(4) - t38;
t21 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t1 = [t121 + t34 * t122 + (t118 * mrSges(2,2) + (mrSges(2,1) + (m(3) + t150) * pkin(1)) * t115 + t120) * g(1) + t40 * t60 + t33 * t21 + (-mrSges(4,1) * t100 - mrSges(5,1) * t38 + mrSges(5,2) * t39 + mrSges(4,3) * t96 + Ifges(5,3) + (mrSges(3,1) * t117 - mrSges(3,2) * t114) * pkin(1)) * t108 + (t40 * mrSges(5,1) + t41 * mrSges(5,2) + t86 * mrSges(4,3) + (-mrSges(3,2) * t117 - t176) * t159) * t109 + m(6) * (t130 * t41 + t17 * t40 + t33 * t5) + (t139 * t34 + t41 * t74 - t172) * t116 + (t138 * t34 - t41 * t72) * t113 + (t114 * t67 + t117 * t66) * t174 + Ifges(2,3) * qJDD(1) + m(4) * (t100 * t42 + t69 * t145 + t30 * t96 + t73 * t86) + m(5) * (-t19 * t40 + t20 * t41 + t38 * t8 + t39 * t9) + (-m(6) * (t107 + t128) + t115 * mrSges(2,2) - m(4) * t148 - m(5) * (t98 + t148) + (-mrSges(2,1) - t174) * t118 + t123) * g(2); t121 + (-mrSges(5,1) * t56 - mrSges(5,2) * t57 + (mrSges(5,2) * t112 - t140) * qJD(3) + (t181 * t117 + t176) * t160) * t109 + t120 * g(1) + t71 * t122 + t177 * t60 + (t137 * t74 + t139 * t71 - t172) * t116 + (-t137 * t72 + t138 * t71) * t113 + (mrSges(4,1) * pkin(2) - mrSges(5,1) * t75 + mrSges(5,2) * t76 + mrSges(4,3) * qJ(3) + Ifges(5,3)) * t108 + t70 * t21 + t123 * g(2) + (-t128 * g(2) + t175 * qJD(3) - t130 * t57 + t177 * t17 + t5 * t70) * m(6) + (-t149 * g(2) + qJD(3) * t132 + t19 * t56 - t20 * t57 + t75 * t8 + t76 * t9) * m(5) + (-pkin(2) * t42 + qJ(3) * t30 + qJD(3) * t73 - (t114 * t69 + t117 * t73) * t160 - t166 * g(2)) * m(4); -t112 * t21 + (-mrSges(5,1) * t112 - mrSges(4,1)) * t108 + (t108 * mrSges(5,2) - t113 * t32 + t116 * t31 + (-t113 * t74 - t116 * t72) * qJD(5)) * t111 + m(6) * (t111 * t124 - t112 * t5) + m(5) * (t111 * t9 + t112 * t8) + m(4) * t42 + (-t111 * t60 - m(4) * t73 + t140 * t109 + (-mrSges(5,2) * t109 - t129) * t112 - m(6) * (t111 * t17 + t175) - m(5) * t132) * t109 + (-t105 * g(1) + t106 * g(2)) * t150; t113 * t31 + t116 * t32 + t129 * qJD(5) + (qJD(5) * t130 + t113 * t2 + t116 * t3 + g(3)) * m(6) + (qJDD(4) + g(3)) * m(5); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t62 + Ifges(6,6) * t61 + Ifges(6,3) * qJDD(5) + g(3) * t79 - t13 * t74 + t14 * t72 + (t126 + t116 * t46 / 0.2e1 + t45 * t171 + t85 * t170 + t133 * t168 + (t125 / 0.2e1 + t113 * t161 / 0.2e1) * t109 + t131 * mrSges(6,3)) * t109 + (-g(1) * t58 - g(2) * t59) * t135;];
tau = t1;
