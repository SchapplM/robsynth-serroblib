% Calculate vector of inverse dynamics joint torques for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:29
% DurationCPUTime: 1.19s
% Computational Cost: add. (1447->203), mult. (2322->284), div. (0->0), fcn. (1359->12), ass. (0->106)
t109 = sin(pkin(7));
t110 = cos(pkin(7));
t124 = -t110 * mrSges(4,1) + t109 * mrSges(4,2);
t163 = -mrSges(3,1) + t124;
t139 = t109 ^ 2 + t110 ^ 2;
t132 = t139 * mrSges(4,3);
t106 = pkin(7) + qJ(4);
t96 = sin(t106);
t97 = cos(t106);
t162 = t97 * mrSges(5,1) - t96 * mrSges(5,2);
t161 = mrSges(3,2) - mrSges(5,3) - mrSges(4,3);
t160 = -t162 + t163;
t107 = qJD(1) + qJD(2);
t116 = cos(qJ(2));
t159 = t107 * t116;
t112 = sin(qJ(4));
t115 = cos(qJ(4));
t67 = -t109 * t112 + t110 * t115;
t54 = t67 * t107;
t68 = t109 * t115 + t110 * t112;
t55 = t68 * t107;
t158 = -mrSges(5,1) * t54 + mrSges(5,2) * t55 + t163 * t107;
t157 = m(3) * pkin(1);
t155 = t55 / 0.2e1;
t153 = Ifges(5,4) * t55;
t152 = pkin(1) * t116;
t111 = -pkin(6) - qJ(3);
t113 = sin(qJ(2));
t146 = pkin(1) * qJD(2);
t136 = qJD(1) * t146;
t142 = pkin(1) * qJDD(1);
t66 = t113 * t142 + t116 * t136;
t108 = qJ(1) + qJ(2);
t98 = sin(t108);
t99 = cos(t108);
t148 = t99 * pkin(2) + t98 * qJ(3);
t147 = pkin(1) * qJD(1);
t103 = qJDD(1) + qJDD(2);
t141 = t103 * t109;
t140 = t103 * t110;
t138 = t116 * t147;
t137 = t113 * t146;
t91 = pkin(3) * t110 + pkin(2);
t60 = t67 * qJD(4);
t28 = t68 * t103 + t107 * t60;
t61 = t68 * qJD(4);
t29 = t67 * t103 - t107 * t61;
t6 = -t29 * mrSges(5,1) + t28 * mrSges(5,2);
t44 = qJ(3) * t103 + qJD(3) * t107 + t66;
t135 = pkin(6) * t103 + t44;
t70 = qJ(3) * t107 + t113 * t147;
t134 = pkin(6) * t107 + t70;
t133 = -t111 * t98 + t99 * t91;
t131 = t139 * t103;
t130 = t139 * t107;
t129 = qJD(3) * t139;
t58 = -mrSges(4,1) * t140 + mrSges(4,2) * t141;
t128 = -g(1) * t98 + g(2) * t99;
t65 = -t113 * t136 + t116 * t142;
t126 = qJD(3) - t138;
t45 = t134 * t109;
t46 = t134 * t110;
t18 = -t112 * t46 - t115 * t45;
t19 = -t112 * t45 + t115 * t46;
t90 = pkin(1) * t113 + qJ(3);
t63 = (-pkin(6) - t90) * t109;
t100 = t110 * pkin(6);
t64 = t110 * t90 + t100;
t33 = -t112 * t64 + t115 * t63;
t34 = t112 * t63 + t115 * t64;
t75 = t111 * t109;
t76 = qJ(3) * t110 + t100;
t40 = -t112 * t76 + t115 * t75;
t41 = t112 * t75 + t115 * t76;
t123 = qJDD(3) - t65;
t120 = t160 * t99 + t161 * t98;
t119 = (m(4) * pkin(2) + m(5) * t91 - t160) * t98 + (-m(4) * qJ(3) + m(5) * t111 + t161) * t99;
t31 = t135 * t109;
t32 = t135 * t110;
t2 = t18 * qJD(4) - t112 * t31 + t115 * t32;
t22 = Ifges(5,2) * t54 + Ifges(5,6) * qJD(4) + t153;
t51 = Ifges(5,4) * t54;
t23 = Ifges(5,1) * t55 + Ifges(5,5) * qJD(4) + t51;
t3 = -t19 * qJD(4) - t112 * t32 - t115 * t31;
t39 = -t91 * t103 + t123;
t52 = -t91 * t107 + t126;
t53 = -pkin(2) * t103 + t123;
t118 = (Ifges(4,4) * t109 + Ifges(4,2) * t110) * t140 + (Ifges(4,1) * t109 + Ifges(4,4) * t110) * t141 + Ifges(3,3) * t103 + t65 * mrSges(3,1) + t53 * t124 - t66 * mrSges(3,2) + t60 * t23 / 0.2e1 + qJD(4) * (Ifges(5,5) * t60 - Ifges(5,6) * t61) / 0.2e1 + t52 * (mrSges(5,1) * t61 + mrSges(5,2) * t60) + t54 * (Ifges(5,4) * t60 - Ifges(5,2) * t61) / 0.2e1 - t61 * t22 / 0.2e1 + (Ifges(5,1) * t60 - Ifges(5,4) * t61) * t155 + t44 * t132 + (-t18 * t60 - t19 * t61) * mrSges(5,3) + (t39 * mrSges(5,2) - t3 * mrSges(5,3) + Ifges(5,1) * t28 + Ifges(5,4) * t29 + Ifges(5,5) * qJDD(4)) * t68 + (-t39 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t28 + Ifges(5,2) * t29 + Ifges(5,6) * qJDD(4)) * t67;
t117 = cos(qJ(1));
t114 = sin(qJ(1));
t101 = t117 * pkin(1);
t93 = -pkin(2) - t152;
t80 = t116 * t146 + qJD(3);
t74 = -t91 - t152;
t69 = -pkin(2) * t107 + t126;
t48 = t67 * t138;
t47 = t68 * t138;
t43 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t55;
t42 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t54;
t25 = -t68 * qJD(3) - t41 * qJD(4);
t24 = t67 * qJD(3) + t40 * qJD(4);
t21 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t29;
t20 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t28;
t11 = -t34 * qJD(4) - t68 * t80;
t10 = t33 * qJD(4) + t67 * t80;
t1 = [(t113 * t66 + t116 * t65) * t157 + t118 + m(5) * (t10 * t19 + t11 * t18 + t52 * t137 + t2 * t34 + t3 * t33 + t39 * t74) + (t80 * t130 + t90 * t131) * mrSges(4,3) + (t117 * mrSges(2,2) + (mrSges(2,1) + (m(3) + m(4) + m(5)) * pkin(1)) * t114 + t119) * g(1) + t93 * t58 + t74 * t6 + t10 * t42 + t11 * t43 + (t114 * mrSges(2,2) - m(4) * (t101 + t148) - m(5) * (t101 + t133) + (-mrSges(2,1) - t157) * t117 + t120) * g(2) + t33 * t20 + t34 * t21 + m(4) * (t69 * t137 + t53 * t93 + t139 * (t44 * t90 + t70 * t80)) + ((mrSges(3,1) * t116 - mrSges(3,2) * t113) * t103 + (-mrSges(3,2) * t159 + t113 * t158) * qJD(2)) * pkin(1) + Ifges(2,3) * qJDD(1); t118 + (qJ(3) * t131 + t107 * t129) * mrSges(4,3) + t119 * g(1) + (-m(4) * t148 - m(5) * t133 + t120) * g(2) - m(5) * (-t18 * t47 + t19 * t48) + m(5) * (t18 * t25 + t19 * t24 + t2 * t41 + t3 * t40 - t39 * t91) - t91 * t6 - pkin(2) * t58 + t40 * t20 + t41 * t21 + (t25 + t47) * t43 + (t24 - t48) * t42 + m(4) * (t139 * t44 * qJ(3) - pkin(2) * t53 + t70 * t129) + ((mrSges(3,2) - t132) * t159 - m(4) * t139 * t116 * t70 + (-m(4) * t69 - m(5) * t52 - t158) * t113) * t147; -t107 ^ 2 * t132 - t54 * t42 + t55 * t43 + t58 + t6 + (t18 * t55 - t19 * t54 + t128 + t39) * m(5) + (-t70 * t130 + t128 + t53) * m(4); Ifges(5,5) * t28 + Ifges(5,6) * t29 + Ifges(5,3) * qJDD(4) - t2 * mrSges(5,2) + t3 * mrSges(5,1) - t52 * (mrSges(5,1) * t55 + mrSges(5,2) * t54) - t55 * (Ifges(5,1) * t54 - t153) / 0.2e1 + t22 * t155 - qJD(4) * (Ifges(5,5) * t54 - Ifges(5,6) * t55) / 0.2e1 - t18 * t42 + t19 * t43 - g(3) * t162 + (t18 * t54 + t19 * t55) * mrSges(5,3) - (-Ifges(5,2) * t55 + t23 + t51) * t54 / 0.2e1 + (g(1) * t99 + g(2) * t98) * (mrSges(5,1) * t96 + mrSges(5,2) * t97);];
tau = t1;
