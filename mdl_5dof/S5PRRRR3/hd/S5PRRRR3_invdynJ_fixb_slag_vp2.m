% Calculate vector of inverse dynamics joint torques for
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:09
% DurationCPUTime: 1.28s
% Computational Cost: add. (1768->203), mult. (2617->277), div. (0->0), fcn. (1292->12), ass. (0->103)
t87 = sin(qJ(5));
t144 = t87 / 0.2e1;
t153 = mrSges(5,2) - mrSges(6,3);
t120 = qJD(5) * t87;
t84 = qJDD(2) + qJDD(3);
t79 = qJDD(4) + t84;
t86 = qJD(2) + qJD(3);
t82 = qJD(4) + t86;
t90 = cos(qJ(5));
t34 = -t120 * t82 + t79 * t90;
t119 = qJD(5) * t90;
t35 = t119 * t82 + t79 * t87;
t13 = -mrSges(6,1) * t34 + mrSges(6,2) * t35;
t124 = pkin(2) * qJD(2);
t89 = sin(qJ(3));
t117 = t89 * t124;
t92 = cos(qJ(3));
t116 = t92 * t124;
t54 = t86 * pkin(3) + t116;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t27 = t117 * t91 + t54 * t88;
t141 = pkin(2) * t92;
t47 = -qJD(3) * t117 + qJDD(2) * t141;
t36 = pkin(3) * t84 + t47;
t123 = qJD(3) * t92;
t48 = (qJD(2) * t123 + qJDD(2) * t89) * pkin(2);
t9 = -qJD(4) * t27 + t36 * t91 - t48 * t88;
t6 = -pkin(4) * t79 - t9;
t152 = m(6) * t6 + t13;
t130 = t82 * t87;
t49 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t130;
t129 = t82 * t90;
t50 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t129;
t100 = t82 * mrSges(5,2) + t87 * t49 - t90 * t50;
t21 = pkin(8) * t82 + t27;
t14 = qJD(1) * t90 - t21 * t87;
t15 = qJD(1) * t87 + t21 * t90;
t107 = -t14 * t87 + t15 * t90;
t151 = m(5) * t27 + m(6) * t107 - t100;
t55 = -t90 * mrSges(6,1) + mrSges(6,2) * t87;
t112 = t88 * t117;
t121 = qJD(4) * t91;
t8 = -qJD(4) * t112 + t54 * t121 + t88 * t36 + t91 * t48;
t5 = pkin(8) * t79 + t8;
t2 = qJD(5) * t14 + qJDD(1) * t87 + t5 * t90;
t3 = -qJD(5) * t15 + qJDD(1) * t90 - t5 * t87;
t150 = t2 * t90 - t3 * t87;
t149 = -mrSges(5,1) + t55;
t111 = mrSges(6,1) * t87 + mrSges(6,2) * t90;
t26 = t54 * t91 - t112;
t20 = -pkin(4) * t82 - t26;
t148 = t20 * t111 + qJD(5) * (Ifges(6,5) * t90 - Ifges(6,6) * t87) / 0.2e1 - (t14 * t90 + t15 * t87) * mrSges(6,3);
t147 = t47 * mrSges(4,1) + Ifges(4,3) * t84;
t142 = m(5) + m(6);
t85 = pkin(9) + qJ(2);
t83 = qJ(3) + t85;
t71 = sin(t83);
t140 = g(1) * t71;
t80 = sin(t85);
t139 = g(1) * t80;
t135 = Ifges(6,4) * t87;
t134 = Ifges(6,4) * t90;
t133 = Ifges(6,2) * t90;
t128 = t88 * t89;
t127 = t89 * t91;
t73 = qJ(4) + t83;
t66 = sin(t73);
t67 = cos(t73);
t125 = t67 * pkin(4) + t66 * pkin(8);
t76 = pkin(3) + t141;
t43 = pkin(2) * t127 + t88 * t76;
t122 = qJD(4) * t88;
t72 = cos(t83);
t65 = pkin(3) * t72;
t118 = t65 + t125;
t115 = m(2) + m(3) + m(4) + m(5);
t33 = t55 * t82;
t113 = t82 * mrSges(5,1) - t33;
t110 = t133 + t135;
t104 = t88 * t92 + t127;
t103 = t91 * t92 - t128;
t42 = -pkin(2) * t128 + t76 * t91;
t101 = t87 * (Ifges(6,1) * t90 - t135);
t99 = t149 * t67 + t153 * t66;
t98 = -t72 * mrSges(4,1) + t71 * mrSges(4,2) + t99;
t97 = (m(6) * pkin(4) - t149) * t66 + (-m(6) * pkin(8) + t153) * t67;
t96 = t72 * mrSges(4,2) + t97;
t28 = Ifges(6,6) * qJD(5) + t110 * t82;
t56 = Ifges(6,4) * t129;
t29 = Ifges(6,1) * t130 + Ifges(6,5) * qJD(5) + t56;
t95 = (Ifges(6,1) * t35 + Ifges(6,4) * t34) * t144 + t90 * (Ifges(6,4) * t35 + Ifges(6,2) * t34) / 0.2e1 + t34 * t110 / 0.2e1 + t35 * (Ifges(6,1) * t87 + t134) / 0.2e1 - t28 * t120 / 0.2e1 + t6 * t55 + Ifges(5,3) * t79 + t9 * mrSges(5,1) + (t29 + t82 * (-Ifges(6,2) * t87 + t134)) * t119 / 0.2e1 + t150 * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t144 + Ifges(6,6) * t90) * qJDD(5) + (t101 * t82 / 0.2e1 + t148) * qJD(5);
t24 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t34;
t25 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t35;
t94 = -t49 * t119 - t50 * t120 + m(6) * (-t119 * t14 - t120 * t15 + t150) + t90 * t24 - t87 * t25;
t93 = -t8 * mrSges(5,2) + t95;
t81 = cos(t85);
t70 = pkin(2) * t81;
t40 = t103 * t124;
t39 = t104 * t124;
t37 = -pkin(4) - t42;
t17 = t76 * t122 + (qJD(3) * t104 + t121 * t89) * pkin(2);
t1 = [t50 * t119 - t49 * t120 + t87 * t24 + t90 * t25 + m(6) * (qJD(5) * t107 + t2 * t87 + t3 * t90) + t115 * qJDD(1) + (-m(6) - t115) * g(3); m(5) * (-t17 * t26 + t42 * t9 + t43 * t8) + (t142 * t139 + (-t123 * t86 - t84 * t89) * mrSges(4,2) + (-qJD(3) * t86 * t89 + t84 * t92) * mrSges(4,1) + (-g(2) * t81 + t47 * t92 + t48 * t89 + t139) * m(4)) * pkin(2) + (-m(6) * (t70 + t118) - t81 * mrSges(3,1) + t80 * mrSges(3,2) - m(5) * (t65 + t70) + t98) * g(2) + t94 * (pkin(8) + t43) + t95 + (t80 * mrSges(3,1) + t81 * mrSges(3,2) + (pkin(3) * t142 + mrSges(4,1)) * t71 + t96) * g(1) + m(6) * (t17 * t20 + t37 * t6) + (-t43 * t79 - t8) * mrSges(5,2) + t37 * t13 - t48 * mrSges(4,2) + t17 * t33 + (-t17 * t82 + t42 * t79) * mrSges(5,1) + Ifges(3,3) * qJDD(2) + t151 * (t76 * t121 + (qJD(3) * t103 - t122 * t89) * pkin(2)) + t147; (-m(6) * t118 + t98) * g(2) - m(5) * (-t26 * t39 + t27 * t40) + t100 * t40 + t113 * t39 + t94 * (pkin(3) * t88 + pkin(8)) + (m(6) * t140 + (t91 * mrSges(5,1) - t88 * mrSges(5,2)) * t79 + (-g(2) * t72 + t8 * t88 + t9 * t91 + t140) * m(5) + ((-m(5) * t26 + m(6) * t20 - t113) * t88 + t151 * t91) * qJD(4)) * pkin(3) + t93 - m(6) * (t107 * t40 + t20 * t39) + (t116 * t86 - t48) * mrSges(4,2) + (t71 * mrSges(4,1) + t96) * g(1) + t86 * mrSges(4,1) * t117 + t147 + t152 * (-pkin(3) * t91 - pkin(4)); t97 * g(1) + t94 * pkin(8) + t113 * t27 + t100 * t26 + (-m(6) * t125 + t99) * g(2) + t93 - m(6) * (t107 * t26 + t20 * t27) - t152 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t35 + Ifges(6,6) * t34 + Ifges(6,3) * qJDD(5) + g(3) * t55 - t14 * t50 + t15 * t49 + (t28 * t144 + (-t101 / 0.2e1 + t133 * t144) * t82 - (t29 + t56) * t90 / 0.2e1 - t148) * t82 + (g(1) * t67 + g(2) * t66) * t111;];
tau = t1;
