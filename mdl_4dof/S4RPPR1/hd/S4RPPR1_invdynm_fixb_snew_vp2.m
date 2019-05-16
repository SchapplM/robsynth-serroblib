% Calculate vector of cutting torques with Newton-Euler for
% S4RPPR1
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-05-04 19:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:08:27
% EndTime: 2019-05-04 19:08:27
% DurationCPUTime: 0.40s
% Computational Cost: add. (4616->123), mult. (6902->139), div. (0->0), fcn. (2808->6), ass. (0->54)
t127 = -pkin(2) - pkin(3);
t108 = sin(pkin(6));
t109 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t110 = sin(qJ(4));
t112 = cos(qJ(4));
t111 = sin(qJ(1));
t113 = cos(qJ(1));
t90 = t111 * g(1) - t113 * g(2);
t86 = qJDD(1) * pkin(1) + t90;
t91 = -t113 * g(1) - t111 * g(2);
t87 = -t114 * pkin(1) + t91;
t82 = t108 * t86 + t109 * t87;
t124 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t82;
t73 = t127 * t114 + t124;
t81 = -t108 * t87 + t109 * t86;
t121 = -t114 * qJ(3) + qJDD(3) - t81;
t76 = t127 * qJDD(1) + t121;
t71 = -t110 * t73 + t112 * t76;
t97 = -qJD(1) + qJD(4);
t95 = t97 ^ 2;
t96 = -qJDD(1) + qJDD(4);
t68 = m(5) * t71 + t96 * mrSges(5,1) - (t95 * mrSges(5,2));
t72 = t110 * t76 + t112 * t73;
t69 = m(5) * t72 - t95 * mrSges(5,1) - t96 * mrSges(5,2);
t63 = -t110 * t68 + t112 * t69;
t77 = -t114 * pkin(2) + t124;
t123 = m(4) * t77 + qJDD(1) * mrSges(4,3) + t63;
t58 = m(3) * t82 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t114 + t123;
t62 = t110 * t69 + t112 * t68;
t79 = -qJDD(1) * pkin(2) + t121;
t119 = -m(4) * t79 + qJDD(1) * mrSges(4,1) + t114 * mrSges(4,3) - t62;
t59 = m(3) * t81 + qJDD(1) * mrSges(3,1) - t114 * mrSges(3,2) + t119;
t52 = t108 * t58 + t109 * t59;
t126 = -m(5) * pkin(3) - mrSges(4,1);
t125 = -t108 * t59 + t109 * t58;
t122 = mrSges(5,1) * t71 - mrSges(5,2) * t72 + Ifges(5,3) * t96;
t105 = g(3) - qJDD(2);
t65 = -mrSges(5,1) * t105 + mrSges(5,3) * t72 + (t95 * Ifges(5,5)) + Ifges(5,6) * t96;
t66 = mrSges(5,2) * t105 - mrSges(5,3) * t71 + Ifges(5,5) * t96 - t95 * Ifges(5,6);
t120 = mrSges(4,2) * t79 + Ifges(4,4) * qJDD(1) + t114 * Ifges(4,6) - pkin(5) * t62 - t110 * t65 + t112 * t66;
t118 = -mrSges(4,2) * t77 + pkin(5) * t63 + t110 * t66 + t112 * t65;
t117 = -mrSges(4,1) * t79 + mrSges(4,3) * t77 + Ifges(4,2) * qJDD(1) - pkin(3) * t62 - t122;
t116 = -mrSges(3,2) * t82 + Ifges(3,3) * qJDD(1) + t117 + qJ(3) * (-t114 * mrSges(4,1) + t123) + pkin(2) * t119 + mrSges(3,1) * t81;
t115 = mrSges(2,1) * t90 - mrSges(2,2) * t91 + Ifges(2,3) * qJDD(1) + pkin(1) * t52 + t116;
t92 = m(4) * t105;
t89 = -m(5) * t105 - t92;
t54 = -mrSges(3,3) * t81 + Ifges(3,5) * qJDD(1) - t114 * Ifges(3,6) - qJ(3) * t89 + (-mrSges(3,2) + mrSges(4,3)) * t105 + t120;
t53 = mrSges(3,3) * t82 - pkin(2) * t89 + (Ifges(4,4) + Ifges(3,5)) * t114 + (Ifges(3,6) - Ifges(4,6)) * qJDD(1) + (mrSges(3,1) - t126) * t105 - t118;
t50 = m(2) * t91 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t125;
t49 = m(2) * t90 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) + t52;
t48 = -mrSges(2,2) * g(3) - mrSges(2,3) * t90 + Ifges(2,5) * qJDD(1) - t114 * Ifges(2,6) - qJ(2) * t52 - t108 * t53 + t109 * t54;
t47 = Ifges(2,6) * qJDD(1) + t114 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t91 + t108 * t54 + t109 * t53 - pkin(1) * (-t92 + (-m(3) - m(5)) * t105) + qJ(2) * t125;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t113 * t48 - t111 * t47 - pkin(4) * (t111 * t50 + t113 * t49), t48, t54, mrSges(4,3) * t105 + t120, t66; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t111 * t48 + t113 * t47 + pkin(4) * (-t111 * t49 + t113 * t50), t47, t53, t117, t65; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t115, t115, t116, -t114 * Ifges(4,4) + Ifges(4,6) * qJDD(1) + t126 * t105 + t118, t122;];
m_new  = t1;
