% Calculate vector of inverse dynamics joint torques for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:30
% DurationCPUTime: 0.99s
% Computational Cost: add. (954->192), mult. (1635->269), div. (0->0), fcn. (824->12), ass. (0->96)
t84 = sin(qJ(4));
t131 = t84 / 0.2e1;
t134 = mrSges(4,2) - mrSges(5,3);
t87 = cos(qJ(4));
t112 = t87 * mrSges(5,1);
t133 = -mrSges(4,1) - t112;
t103 = mrSges(5,1) * t84 + mrSges(5,2) * t87;
t111 = pkin(1) * qJD(1);
t85 = sin(qJ(2));
t105 = t85 * t111;
t80 = qJD(1) + qJD(2);
t88 = cos(qJ(2));
t52 = pkin(2) * t80 + t88 * t111;
t82 = sin(pkin(7));
t83 = cos(pkin(7));
t22 = -t82 * t105 + t52 * t83;
t16 = -pkin(3) * t80 - t22;
t132 = t16 * t103 + qJD(4) * (Ifges(5,5) * t87 - Ifges(5,6) * t84) / 0.2e1;
t23 = t83 * t105 + t82 * t52;
t17 = pkin(6) * t80 + t23;
t13 = qJD(3) * t87 - t17 * t84;
t14 = qJD(3) * t84 + t17 * t87;
t99 = -t13 * t84 + t14 * t87;
t129 = m(4) + m(5);
t128 = pkin(1) * t88;
t81 = qJ(1) + qJ(2);
t77 = cos(t81);
t70 = pkin(2) * t77;
t127 = pkin(2) * t82;
t126 = pkin(2) * t83;
t86 = sin(qJ(1));
t125 = g(1) * t86;
t47 = -qJD(2) * t105 + qJDD(1) * t128;
t79 = qJDD(1) + qJDD(2);
t31 = pkin(2) * t79 + t47;
t110 = qJD(2) * t88;
t48 = (qJD(1) * t110 + qJDD(1) * t85) * pkin(1);
t12 = t82 * t31 + t83 * t48;
t6 = pkin(6) * t79 + t12;
t2 = qJD(4) * t13 + qJDD(3) * t84 + t6 * t87;
t124 = t2 * t87;
t109 = qJD(4) * t14;
t3 = qJDD(3) * t87 - t6 * t84 - t109;
t123 = t3 * t84;
t122 = Ifges(5,4) * t84;
t121 = Ifges(5,4) * t87;
t120 = Ifges(5,2) * t87;
t117 = t80 * t84;
t116 = t80 * t87;
t115 = t82 * t85;
t114 = t83 * t85;
t113 = t84 * mrSges(5,2);
t71 = pkin(2) + t128;
t39 = pkin(1) * t114 + t82 * t71;
t108 = qJD(4) * t84;
t107 = qJD(4) * t87;
t75 = pkin(7) + t81;
t66 = sin(t75);
t67 = cos(t75);
t106 = t67 * pkin(3) + t66 * pkin(6) + t70;
t55 = -t112 + t113;
t102 = t120 + t122;
t100 = t13 * t87 + t14 * t84;
t11 = t31 * t83 - t48 * t82;
t50 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t117;
t51 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t116;
t98 = -t84 * t50 + t87 * t51;
t38 = -pkin(1) * t115 + t71 * t83;
t97 = pkin(1) * (t82 * t88 + t114);
t96 = pkin(1) * (t83 * t88 - t115);
t94 = t84 * (Ifges(5,1) * t87 - t122);
t76 = sin(t81);
t93 = -t77 * mrSges(3,1) + t76 * mrSges(3,2) + t133 * t67 + t134 * t66;
t92 = -t100 * qJD(4) - t123;
t91 = t77 * mrSges(3,2) + (t129 * pkin(2) + mrSges(3,1)) * t76 - t66 * t113 + t134 * t67;
t29 = Ifges(5,6) * qJD(4) + t102 * t80;
t57 = Ifges(5,4) * t116;
t30 = Ifges(5,1) * t117 + Ifges(5,5) * qJD(4) + t57;
t42 = -t80 * t108 + t79 * t87;
t43 = t80 * t107 + t79 * t84;
t5 = -pkin(3) * t79 - t11;
t90 = -t48 * mrSges(3,2) + mrSges(5,3) * t124 + t5 * t55 + t11 * mrSges(4,1) + t42 * t102 / 0.2e1 + t43 * (Ifges(5,1) * t84 + t121) / 0.2e1 - t29 * t108 / 0.2e1 + t47 * mrSges(3,1) + (Ifges(5,1) * t43 + Ifges(5,4) * t42) * t131 + t87 * (Ifges(5,4) * t43 + Ifges(5,2) * t42) / 0.2e1 + (Ifges(4,3) + Ifges(3,3)) * t79 + (t30 + t80 * (-Ifges(5,2) * t84 + t121)) * t107 / 0.2e1 + (0.2e1 * Ifges(5,5) * t131 + Ifges(5,6) * t87) * qJDD(4) + (t94 * t80 / 0.2e1 + t132) * qJD(4);
t89 = cos(qJ(1));
t78 = t89 * pkin(1);
t69 = -pkin(3) - t126;
t68 = pkin(6) + t127;
t62 = t67 * pkin(6);
t41 = t55 * t80;
t36 = qJD(1) * t96;
t35 = qJD(2) * t97;
t34 = qJD(1) * t97;
t32 = -pkin(3) - t38;
t26 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t43;
t25 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t42;
t15 = -mrSges(5,1) * t42 + mrSges(5,2) * t43;
t1 = [(-m(5) * (t78 + t106) + t67 * t113 - t89 * mrSges(2,1) + t86 * mrSges(2,2) - m(4) * (t70 + t78) + t93) * g(2) + (-t50 * t107 - t51 * t108 - t84 * t26 + t87 * t25 + m(5) * (-t13 * t107 - t14 * t108 - t123 + t124)) * (pkin(6) + t39) + (t129 * t125 + (-t80 * t110 - t85 * t79) * mrSges(3,2) + (-qJD(2) * t85 * t80 + t88 * t79) * mrSges(3,1) + (-g(2) * t89 + t47 * t88 + t48 * t85 + t125) * m(3)) * pkin(1) + m(5) * (t16 * t35 + t32 * t5) + (-m(5) * t62 + t86 * mrSges(2,1) + t89 * mrSges(2,2) + (m(5) * pkin(3) - t133) * t66 + t91) * g(1) + t92 * mrSges(5,3) + (-t39 * t79 - t12) * mrSges(4,2) + (-t35 * t80 + t38 * t79) * mrSges(4,1) + t90 + t32 * t15 + t35 * t41 + m(4) * (t11 * t38 + t12 * t39 - t22 * t35) + Ifges(2,3) * qJDD(1) + (m(4) * t23 + t99 * m(5) - t80 * mrSges(4,2) + t98) * qJD(2) * t96; (g(2) * mrSges(5,2) * t67 + t36 * t50 + (-qJD(4) * t51 - t26) * t68 + (-t3 - t109) * mrSges(5,3)) * t84 + (t34 * mrSges(4,1) + t36 * mrSges(4,2) + (mrSges(3,1) * t85 + mrSges(3,2) * t88) * t111) * t80 + (g(1) * mrSges(5,1) * t66 + t68 * t25 - t36 * t51 + (-t13 * mrSges(5,3) - t50 * t68) * qJD(4)) * t87 + (t66 * mrSges(4,1) + t91) * g(1) + (-t79 * t127 - t12) * mrSges(4,2) + t93 * g(2) + t90 + t69 * t15 - t34 * t41 + t79 * mrSges(4,1) * t126 + ((t11 * t83 + t12 * t82) * pkin(2) - t70 * g(2) + t22 * t34 - t23 * t36) * m(4) + ((t92 + t124) * t68 + t5 * t69 + (pkin(3) * t66 - t62) * g(1) - t16 * t34 - t99 * t36 - t106 * g(2)) * m(5); t84 * t25 + t87 * t26 + t98 * qJD(4) + (t99 * qJD(4) + t2 * t84 + t3 * t87 - g(3)) * m(5) + (qJDD(3) - g(3)) * m(4); t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t43 + Ifges(5,6) * t42 + Ifges(5,3) * qJDD(4) + g(3) * t55 - t13 * t51 + t14 * t50 + (t29 * t131 + (-t94 / 0.2e1 + t120 * t131) * t80 + t100 * mrSges(5,3) - (t30 + t57) * t87 / 0.2e1 - t132) * t80 + (g(1) * t67 + g(2) * t66) * t103;];
tau = t1;
