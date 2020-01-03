% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR16_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:48
% EndTime: 2019-12-31 20:44:52
% DurationCPUTime: 1.57s
% Computational Cost: add. (15313->175), mult. (34668->231), div. (0->0), fcn. (24225->10), ass. (0->93)
t129 = -2 * qJD(3);
t82 = sin(pkin(5));
t114 = qJD(1) * t82;
t86 = sin(qJ(2));
t106 = t86 * t114;
t83 = cos(pkin(5));
t77 = qJD(1) * t83 + qJD(2);
t128 = (pkin(2) * t77 + t129) * t106;
t122 = t82 * t86;
t87 = sin(qJ(1));
t91 = cos(qJ(1));
t105 = t87 * g(1) - g(2) * t91;
t92 = qJD(1) ^ 2;
t60 = pkin(7) * t82 * t92 + qJDD(1) * pkin(1) + t105;
t124 = t60 * t83;
t104 = -g(1) * t91 - g(2) * t87;
t111 = qJDD(1) * t82;
t61 = -pkin(1) * t92 + pkin(7) * t111 + t104;
t90 = cos(qJ(2));
t101 = -g(3) * t122 + t86 * t124 + t90 * t61;
t113 = qJD(1) * t90;
t107 = t82 * t113;
t62 = (-pkin(2) * t90 - qJ(3) * t86) * t114;
t75 = t77 ^ 2;
t76 = qJDD(1) * t83 + qJDD(2);
t127 = t75 * pkin(2) - t76 * qJ(3) - t62 * t107 + t129 * t77 - t101;
t126 = -pkin(2) - pkin(8);
t125 = t83 * g(3);
t123 = t82 ^ 2 * t92;
t121 = t82 * t90;
t120 = mrSges(3,1) - mrSges(4,2);
t119 = mrSges(3,3) + mrSges(4,1);
t109 = t90 ^ 2 * t123;
t65 = pkin(3) * t106 - pkin(8) * t77;
t66 = (qJD(2) * t113 + qJDD(1) * t86) * t82;
t67 = -qJD(2) * t106 + t111 * t90;
t22 = -pkin(3) * t109 - t125 - t66 * qJ(3) + t126 * t67 + (-t60 + (-qJ(3) * t77 * t90 - t65 * t86) * qJD(1)) * t82 + t128;
t117 = g(3) * t121 + t86 * t61;
t99 = -t75 * qJ(3) + t62 * t106 + qJDD(3) + t117;
t24 = t66 * pkin(3) + t126 * t76 + (-pkin(3) * t114 * t77 - pkin(8) * t123 * t86 - t124) * t90 + t99;
t85 = sin(qJ(4));
t89 = cos(qJ(4));
t118 = t22 * t89 + t24 * t85;
t59 = mrSges(4,1) * t106 + mrSges(4,2) * t77;
t116 = mrSges(3,1) * t77 - mrSges(3,3) * t106 - t59;
t63 = (mrSges(4,2) * t90 - mrSges(4,3) * t86) * t114;
t115 = t63 + (-mrSges(3,1) * t90 + mrSges(3,2) * t86) * t114;
t50 = -t107 * t89 - t77 * t85;
t51 = -t107 * t85 + t77 * t89;
t39 = -pkin(4) * t50 - pkin(9) * t51;
t55 = qJDD(4) + t66;
t71 = qJD(4) + t106;
t69 = t71 ^ 2;
t17 = -pkin(4) * t69 + pkin(9) * t55 + t39 * t50 + t118;
t36 = -qJD(4) * t51 - t67 * t89 - t76 * t85;
t37 = qJD(4) * t50 - t67 * t85 + t76 * t89;
t93 = t67 * pkin(3) - pkin(8) * t109 + t77 * t65 - t127;
t18 = (-t50 * t71 - t37) * pkin(9) + (t51 * t71 - t36) * pkin(4) + t93;
t84 = sin(qJ(5));
t88 = cos(qJ(5));
t40 = -t51 * t84 + t71 * t88;
t26 = qJD(5) * t40 + t37 * t88 + t55 * t84;
t41 = t51 * t88 + t71 * t84;
t30 = -mrSges(6,1) * t40 + mrSges(6,2) * t41;
t49 = qJD(5) - t50;
t31 = -mrSges(6,2) * t49 + mrSges(6,3) * t40;
t34 = qJDD(5) - t36;
t14 = m(6) * (-t17 * t84 + t18 * t88) - t26 * mrSges(6,3) + t34 * mrSges(6,1) - t41 * t30 + t49 * t31;
t25 = -qJD(5) * t41 - t37 * t84 + t55 * t88;
t32 = mrSges(6,1) * t49 - mrSges(6,3) * t41;
t15 = m(6) * (t17 * t88 + t18 * t84) + t25 * mrSges(6,3) - t34 * mrSges(6,2) + t40 * t30 - t49 * t32;
t38 = -mrSges(5,1) * t50 + mrSges(5,2) * t51;
t43 = mrSges(5,1) * t71 - mrSges(5,3) * t51;
t10 = m(5) * t118 - mrSges(5,2) * t55 + mrSges(5,3) * t36 - t14 * t84 + t15 * t88 + t38 * t50 - t43 * t71;
t103 = -t82 * t60 - t125;
t102 = -t22 * t85 + t24 * t89;
t42 = -mrSges(5,2) * t71 + mrSges(5,3) * t50;
t94 = m(6) * (-pkin(4) * t55 - pkin(9) * t69 + t39 * t51 - t102) - t25 * mrSges(6,1) + t26 * mrSges(6,2) - t40 * t31 + t41 * t32;
t11 = m(5) * t102 + mrSges(5,1) * t55 - mrSges(5,3) * t37 - t38 * t51 + t42 * t71 - t94;
t58 = -mrSges(4,1) * t107 - mrSges(4,3) * t77;
t100 = t89 * t10 - t85 * t11 + m(4) * (-t67 * pkin(2) + (-t107 * t77 - t66) * qJ(3) + t103 + t128) + t58 * t107 - t66 * mrSges(4,3);
t57 = -mrSges(3,2) * t77 + mrSges(3,3) * t107;
t5 = m(3) * t103 + t66 * mrSges(3,2) - t120 * t67 + (t116 * t86 - t57 * t90) * t114 + t100;
t108 = t90 * t124;
t98 = -m(4) * (-t76 * pkin(2) - t108 + t99) - t85 * t10 - t89 * t11;
t6 = m(3) * (t108 - t117) + (t57 - t58) * t77 + t120 * t76 - t119 * t66 - t115 * t106 + t98;
t96 = m(5) * t93 - t36 * mrSges(5,1) + t37 * mrSges(5,2) + t14 * t88 + t15 * t84 - t50 * t42 + t51 * t43;
t95 = -m(4) * t127 + t96;
t8 = m(3) * t101 - t116 * t77 + (-mrSges(3,2) + mrSges(4,3)) * t76 + t119 * t67 + t115 * t107 + t95;
t110 = t121 * t6 + t122 * t8 + t5 * t83;
t2 = m(2) * t104 - mrSges(2,1) * t92 - qJDD(1) * mrSges(2,2) - t6 * t86 + t8 * t90;
t1 = m(2) * t105 + qJDD(1) * mrSges(2,1) - t92 * mrSges(2,2) - t82 * t5 + (t6 * t90 + t8 * t86) * t83;
t3 = [-m(1) * g(1) - t1 * t87 + t2 * t91, t2, t8, t67 * mrSges(4,2) - t106 * t59 + t100, t10, t15; -m(1) * g(2) + t1 * t91 + t2 * t87, t1, t6, -t67 * mrSges(4,1) - t76 * mrSges(4,3) - t107 * t63 - t77 * t59 - t95, t11, t14; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t5, t66 * mrSges(4,1) + t76 * mrSges(4,2) + t106 * t63 + t77 * t58 - t98, t96, t94;];
f_new = t3;
