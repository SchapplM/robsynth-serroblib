% Calculate vector of cutting torques with Newton-Euler for
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2019-05-04 19:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:59:47
% EndTime: 2019-05-04 18:59:48
% DurationCPUTime: 0.42s
% Computational Cost: add. (5903->104), mult. (7537->123), div. (0->0), fcn. (3888->6), ass. (0->49)
t121 = -mrSges(1,2) - mrSges(2,3);
t102 = sin(pkin(6));
t103 = cos(pkin(6));
t110 = qJD(2) ^ 2;
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t101 = -g(2) + qJDD(1);
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t86 = t105 * g(1) + t107 * t101;
t82 = qJDD(2) * pkin(2) + t86;
t87 = -t107 * g(1) + t105 * t101;
t83 = -t110 * pkin(2) + t87;
t77 = -t102 * t83 + t103 * t82;
t74 = qJDD(2) * pkin(3) + t77;
t78 = t102 * t82 + t103 * t83;
t75 = -t110 * pkin(3) + t78;
t72 = -t104 * t75 + t106 * t74;
t97 = qJD(2) + qJD(4);
t95 = t97 ^ 2;
t96 = qJDD(2) + qJDD(4);
t69 = m(5) * t72 + t96 * mrSges(5,1) - t95 * mrSges(5,2);
t73 = t104 * t74 + t106 * t75;
t70 = m(5) * t73 - t95 * mrSges(5,1) - t96 * mrSges(5,2);
t63 = t104 * t70 + t106 * t69;
t60 = m(4) * t77 + qJDD(2) * mrSges(4,1) - t110 * mrSges(4,2) + t63;
t118 = -t104 * t69 + t106 * t70;
t61 = m(4) * t78 - t110 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t118;
t54 = t102 * t61 + t103 * t60;
t100 = -g(3) + qJDD(3);
t120 = (m(4) + m(5)) * t100;
t119 = -t102 * t60 + t103 * t61;
t51 = m(3) * t86 + qJDD(2) * mrSges(3,1) - t110 * mrSges(3,2) + t54;
t52 = m(3) * t87 - t110 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t119;
t117 = -t105 * t51 + t107 * t52;
t116 = mrSges(5,1) * t72 - mrSges(5,2) * t73 + Ifges(5,3) * t96;
t64 = -mrSges(5,1) * t100 + mrSges(5,3) * t73 + t95 * Ifges(5,5) + Ifges(5,6) * t96;
t65 = mrSges(5,2) * t100 - mrSges(5,3) * t72 + Ifges(5,5) * t96 - t95 * Ifges(5,6);
t55 = Ifges(4,6) * qJDD(2) + t110 * Ifges(4,5) + mrSges(4,3) * t78 + t104 * t65 + t106 * t64 + pkin(5) * t118 + (-pkin(3) * m(5) - mrSges(4,1)) * t100;
t56 = mrSges(4,2) * t100 - mrSges(4,3) * t77 + Ifges(4,5) * qJDD(2) - t110 * Ifges(4,6) - pkin(5) * t63 - t104 * t64 + t106 * t65;
t44 = mrSges(3,1) * g(3) + mrSges(3,3) * t87 + t110 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t120 + qJ(3) * t119 + t102 * t56 + t103 * t55;
t47 = -mrSges(3,2) * g(3) - mrSges(3,3) * t86 + Ifges(3,5) * qJDD(2) - t110 * Ifges(3,6) - qJ(3) * t54 - t102 * t55 + t103 * t56;
t115 = pkin(4) * t117 + t105 * t47 + t107 * t44 + mrSges(2,2) * g(1) + mrSges(2,1) * g(3) + pkin(1) * (m(3) * g(3) - t120);
t49 = t105 * t52 + t107 * t51;
t114 = mrSges(2,2) * t101 - pkin(4) * t49 - t105 * t44 + t107 * t47;
t113 = mrSges(4,1) * t77 - mrSges(4,2) * t78 + Ifges(4,3) * qJDD(2) + pkin(3) * t63 + t116;
t112 = mrSges(3,1) * t86 - mrSges(3,2) * t87 + Ifges(3,3) * qJDD(2) + pkin(2) * t54 + t113;
t111 = mrSges(2,1) * t101 + pkin(1) * t49 + t112;
t1 = [mrSges(1,3) * g(2) + qJ(1) * t120 + (qJ(1) * (-m(2) - m(3)) + t121) * g(3) + t114, -mrSges(2,3) * g(3) + t114, t47, t56, t65; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t115, -mrSges(2,3) * g(1) - t111, t44, t55, t64; t111 - qJ(1) * t117 + (qJ(1) * m(2) - t121) * g(1) - mrSges(1,1) * g(2), t115, t112, t113, t116;];
m_new  = t1;
