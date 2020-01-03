% Calculate vector of inverse dynamics joint torques for
% S4RPRR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:33
% EndTime: 2019-12-31 16:51:35
% DurationCPUTime: 1.21s
% Computational Cost: add. (950->171), mult. (1480->229), div. (0->0), fcn. (660->6), ass. (0->82)
t60 = sin(qJ(4));
t62 = cos(qJ(4));
t38 = -t62 * mrSges(5,1) + t60 * mrSges(5,2);
t129 = mrSges(4,1) - t38;
t128 = m(5) * pkin(3) + t129;
t65 = -pkin(1) - pkin(2);
t43 = qJD(1) * t65 + qJD(2);
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t90 = qJ(2) * qJD(1);
t18 = t43 * t63 - t61 * t90;
t127 = qJD(3) * t18;
t115 = t60 / 0.2e1;
t57 = -qJD(1) + qJD(3);
t10 = -pkin(3) * t57 - t18;
t75 = mrSges(5,1) * t60 + mrSges(5,2) * t62;
t126 = t10 * t75;
t124 = mrSges(3,1) + mrSges(2,1);
t123 = -mrSges(3,3) + mrSges(2,2);
t36 = t63 * qJ(2) + t61 * t65;
t42 = qJDD(1) * t65 + qJDD(2);
t88 = qJD(1) * qJD(2);
t45 = qJDD(1) * qJ(2) + t88;
t19 = t43 * t61 + t63 * t90;
t93 = qJD(3) * t19;
t6 = t42 * t63 - t45 * t61 - t93;
t80 = t129 * t57;
t56 = -qJDD(1) + qJDD(3);
t92 = qJD(4) * t60;
t21 = t56 * t62 - t57 * t92;
t91 = qJD(4) * t62;
t22 = t56 * t60 + t57 * t91;
t122 = t62 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t21) - t60 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t22);
t101 = t60 * mrSges(5,3);
t33 = qJD(4) * mrSges(5,1) - t101 * t57;
t98 = t62 * mrSges(5,3);
t34 = -qJD(4) * mrSges(5,2) + t57 * t98;
t69 = t57 * mrSges(4,2) + t60 * t33 - t62 * t34;
t4 = -pkin(3) * t56 - t6;
t9 = -mrSges(5,1) * t21 + mrSges(5,2) * t22;
t121 = m(5) * t4 + t9;
t11 = pkin(6) * t57 + t19;
t81 = t11 * (t60 ^ 2 + t62 ^ 2);
t113 = sin(qJ(1));
t64 = cos(qJ(1));
t24 = -t113 * t61 - t64 * t63;
t25 = -t113 * t63 + t61 * t64;
t120 = t25 * mrSges(4,2) + t128 * t24;
t119 = t24 * mrSges(4,2) - t128 * t25;
t118 = -m(5) * t10 + t80;
t117 = -m(5) * t81 + t69;
t5 = t61 * t42 + t63 * t45 + t127;
t3 = pkin(6) * t56 + t5;
t1 = -t11 * t92 + t3 * t62;
t2 = -t11 * t91 - t3 * t60;
t78 = t1 * t62 - t2 * t60;
t116 = m(5) * t78 - t33 * t91 - t34 * t92 + t122;
t112 = Ifges(5,4) * t60;
t111 = Ifges(5,4) * t62;
t110 = Ifges(5,2) * t62;
t107 = t56 * mrSges(4,1);
t106 = t56 * mrSges(4,2);
t103 = t57 * t60;
t102 = t57 * t62;
t95 = t64 * pkin(1) + t113 * qJ(2);
t89 = qJDD(1) * mrSges(3,1);
t86 = t64 * pkin(2) + t95;
t83 = -m(5) * pkin(6) - mrSges(5,3);
t79 = -pkin(1) * t113 + t64 * qJ(2);
t74 = Ifges(5,1) * t62 - t112;
t40 = t110 + t112;
t73 = t62 * t33 + t60 * t34;
t35 = -qJ(2) * t61 + t63 * t65;
t27 = qJD(4) * (Ifges(5,5) * t62 - Ifges(5,6) * t60);
t68 = -pkin(2) * t113 + t79;
t14 = Ifges(5,6) * qJD(4) + t40 * t57;
t44 = Ifges(5,4) * t102;
t15 = Ifges(5,1) * t103 + Ifges(5,5) * qJD(4) + t44;
t67 = Ifges(4,3) * t56 + qJD(4) * t126 + t4 * t38 + t6 * mrSges(4,1) + t1 * t98 + t21 * t40 / 0.2e1 + t22 * (Ifges(5,1) * t60 + t111) / 0.2e1 - t5 * mrSges(4,2) + (Ifges(5,1) * t22 + Ifges(5,4) * t21) * t115 + t62 * (Ifges(5,4) * t22 + Ifges(5,2) * t21) / 0.2e1 - t14 * t92 / 0.2e1 + t15 * t91 / 0.2e1 - t2 * t101 + (0.2e1 * Ifges(5,5) * t115 + Ifges(5,6) * t62) * qJDD(4) + (t27 + t74 * t103 + (-Ifges(5,2) * t60 + t111) * t102) * qJD(4) / 0.2e1;
t66 = qJD(1) ^ 2;
t50 = -qJDD(1) * pkin(1) + qJDD(2);
t7 = [t35 * t107 + pkin(1) * t89 - t67 + m(4) * (t35 * t6 + t36 * t5) - t50 * mrSges(3,1) - t36 * t106 + m(3) * (-pkin(1) * t50 + (t45 + t88) * qJ(2)) + t121 * (pkin(3) - t35) + (-m(4) * t18 - t118) * (qJD(2) * t61 + t36 * qJD(3)) + (Ifges(3,2) + Ifges(2,3)) * qJDD(1) + (m(4) * t19 - t117) * (qJD(2) * t63 + qJD(3) * t35) + 0.2e1 * t45 * mrSges(3,3) + t116 * (-pkin(6) + t36) + (-m(3) * t95 - m(5) * (pkin(6) * t25 + t86) - t25 * mrSges(5,3) - m(4) * t86 - t124 * t64 + t123 * t113 + t120) * g(2) + (-m(4) * t68 - m(5) * (t24 * pkin(6) + t68) - t24 * mrSges(5,3) - m(3) * t79 + t123 * t64 + t124 * t113 + t119) * g(1); -t89 - t66 * mrSges(3,3) + (-t66 * qJ(2) + t50) * m(3) + (t107 - t9 - t69 * qJD(3) + m(4) * (t6 + t93) + m(5) * (qJD(3) * t81 - t4)) * t63 + (-t106 - t73 * qJD(4) - t80 * qJD(3) + m(4) * (t5 - t127) + m(5) * (qJD(3) * t10 + t78) + t122) * t61 + (t80 * t61 + t69 * t63 - m(4) * (-t18 * t61 + t19 * t63) - m(5) * (t10 * t61 + t63 * t81)) * qJD(1) + (m(5) + m(4) + m(3)) * (-g(1) * t113 + g(2) * t64); t67 + (-t25 * t83 - t120) * g(2) + (-t24 * t83 - t119) * g(1) - t121 * pkin(3) + t118 * t19 + t117 * t18 + t116 * pkin(6); t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t22 + Ifges(5,6) * t21 + Ifges(5,3) * qJDD(4) - g(3) * t38 + t73 * t11 + (-t126 + t14 * t115 - t27 / 0.2e1 + (-t60 * t74 / 0.2e1 + t110 * t115) * t57 - (t15 + t44) * t62 / 0.2e1) * t57 + (-g(1) * t24 - g(2) * t25) * t75;];
tau = t7;
