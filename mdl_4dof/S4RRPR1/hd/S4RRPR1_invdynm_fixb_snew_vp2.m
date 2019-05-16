% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR1
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
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:21:31
% EndTime: 2019-05-04 19:21:32
% DurationCPUTime: 0.91s
% Computational Cost: add. (15839->120), mult. (21234->147), div. (0->0), fcn. (11538->8), ass. (0->58)
t120 = sin(qJ(2));
t123 = cos(qJ(2));
t114 = qJD(1) + qJD(2);
t112 = t114 ^ 2;
t113 = qJDD(1) + qJDD(2);
t117 = sin(pkin(7));
t118 = cos(pkin(7));
t119 = sin(qJ(4));
t122 = cos(qJ(4));
t108 = qJD(4) + t114;
t106 = t108 ^ 2;
t107 = qJDD(4) + t113;
t121 = sin(qJ(1));
t124 = cos(qJ(1));
t101 = t121 * g(1) - t124 * g(2);
t98 = qJDD(1) * pkin(1) + t101;
t102 = -t124 * g(1) - t121 * g(2);
t125 = qJD(1) ^ 2;
t99 = -t125 * pkin(1) + t102;
t93 = -t120 * t99 + t123 * t98;
t90 = t113 * pkin(2) + t93;
t94 = t120 * t98 + t123 * t99;
t91 = -t112 * pkin(2) + t94;
t85 = -t117 * t91 + t118 * t90;
t82 = t113 * pkin(3) + t85;
t86 = t117 * t90 + t118 * t91;
t83 = -t112 * pkin(3) + t86;
t80 = -t119 * t83 + t122 * t82;
t77 = m(5) * t80 + t107 * mrSges(5,1) - t106 * mrSges(5,2);
t81 = t119 * t82 + t122 * t83;
t78 = m(5) * t81 - t106 * mrSges(5,1) - t107 * mrSges(5,2);
t71 = t119 * t78 + t122 * t77;
t68 = m(4) * t85 + t113 * mrSges(4,1) - t112 * mrSges(4,2) + t71;
t131 = -t119 * t77 + t122 * t78;
t69 = m(4) * t86 - t112 * mrSges(4,1) - t113 * mrSges(4,2) + t131;
t62 = t117 * t69 + t118 * t68;
t59 = m(3) * t93 + t113 * mrSges(3,1) - t112 * mrSges(3,2) + t62;
t132 = -t117 * t68 + t118 * t69;
t60 = m(3) * t94 - t112 * mrSges(3,1) - t113 * mrSges(3,2) + t132;
t55 = t120 * t60 + t123 * t59;
t116 = -g(3) + qJDD(3);
t133 = (m(4) + m(5)) * t116;
t130 = -t120 * t59 + t123 * t60;
t129 = mrSges(5,1) * t80 - mrSges(5,2) * t81 + Ifges(5,3) * t107;
t128 = mrSges(4,1) * t85 - mrSges(4,2) * t86 + Ifges(4,3) * t113 + pkin(3) * t71 + t129;
t127 = mrSges(3,1) * t93 - mrSges(3,2) * t94 + Ifges(3,3) * t113 + pkin(2) * t62 + t128;
t126 = mrSges(2,1) * t101 - mrSges(2,2) * t102 + Ifges(2,3) * qJDD(1) + pkin(1) * t55 + t127;
t73 = mrSges(5,2) * t116 - mrSges(5,3) * t80 + Ifges(5,5) * t107 - t106 * Ifges(5,6);
t72 = -mrSges(5,1) * t116 + mrSges(5,3) * t81 + t106 * Ifges(5,5) + Ifges(5,6) * t107;
t64 = mrSges(4,2) * t116 - mrSges(4,3) * t85 + Ifges(4,5) * t113 - t112 * Ifges(4,6) - pkin(6) * t71 - t119 * t72 + t122 * t73;
t63 = Ifges(4,6) * t113 + t112 * Ifges(4,5) + mrSges(4,3) * t86 + t119 * t73 + t122 * t72 + pkin(6) * t131 + (-pkin(3) * m(5) - mrSges(4,1)) * t116;
t53 = m(2) * t102 - t125 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t130;
t52 = m(2) * t101 + qJDD(1) * mrSges(2,1) - t125 * mrSges(2,2) + t55;
t51 = -mrSges(3,2) * g(3) - mrSges(3,3) * t93 + Ifges(3,5) * t113 - t112 * Ifges(3,6) - qJ(3) * t62 - t117 * t63 + t118 * t64;
t50 = mrSges(3,1) * g(3) + mrSges(3,3) * t94 + t112 * Ifges(3,5) + Ifges(3,6) * t113 - pkin(2) * t133 + qJ(3) * t132 + t117 * t64 + t118 * t63;
t49 = -mrSges(2,2) * g(3) - mrSges(2,3) * t101 + Ifges(2,5) * qJDD(1) - t125 * Ifges(2,6) - pkin(5) * t55 - t120 * t50 + t123 * t51;
t48 = Ifges(2,6) * qJDD(1) + t125 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t102 + t120 * t51 + t123 * t50 - pkin(1) * (-m(3) * g(3) + t133) + pkin(5) * t130;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t124 * t49 - t121 * t48 - pkin(4) * (t121 * t53 + t124 * t52), t49, t51, t64, t73; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t121 * t49 + t124 * t48 + pkin(4) * (-t121 * t52 + t124 * t53), t48, t50, t63, t72; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t126, t126, t127, t128, t129;];
m_new  = t1;
