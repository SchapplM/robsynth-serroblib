% Calculate vector of cutting torques with Newton-Euler for
% S4RPRP1
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-05-04 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:13:50
% EndTime: 2019-05-04 19:13:51
% DurationCPUTime: 0.50s
% Computational Cost: add. (5856->120), mult. (8288->138), div. (0->0), fcn. (3974->6), ass. (0->50)
t104 = (qJD(1) + qJD(3));
t102 = t104 ^ 2;
t103 = qJDD(1) + qJDD(3);
t110 = sin(qJ(3));
t112 = cos(qJ(3));
t108 = sin(pkin(6));
t109 = cos(pkin(6));
t111 = sin(qJ(1));
t113 = cos(qJ(1));
t90 = t111 * g(1) - t113 * g(2);
t87 = qJDD(1) * pkin(1) + t90;
t114 = qJD(1) ^ 2;
t91 = -t113 * g(1) - t111 * g(2);
t88 = -t114 * pkin(1) + t91;
t82 = -t108 * t88 + t109 * t87;
t79 = qJDD(1) * pkin(2) + t82;
t83 = t108 * t87 + t109 * t88;
t80 = -t114 * pkin(2) + t83;
t76 = t110 * t79 + t112 * t80;
t71 = -t102 * pkin(3) + t103 * qJ(4) + (2 * qJD(4) * t104) + t76;
t126 = mrSges(5,2) * t71;
t125 = -mrSges(4,1) - mrSges(5,1);
t124 = m(5) * t71 + t103 * mrSges(5,3);
t64 = m(4) * t76 - t103 * mrSges(4,2) + t125 * t102 + t124;
t75 = -t110 * t80 + t112 * t79;
t73 = -t103 * pkin(3) - t102 * qJ(4) + qJDD(4) - t75;
t119 = -m(5) * t73 + t103 * mrSges(5,1) + t102 * mrSges(5,3);
t66 = m(4) * t75 + t103 * mrSges(4,1) - t102 * mrSges(4,2) + t119;
t59 = t110 * t64 + t112 * t66;
t56 = m(3) * t82 + qJDD(1) * mrSges(3,1) - t114 * mrSges(3,2) + t59;
t120 = -t110 * t66 + t112 * t64;
t57 = m(3) * t83 - t114 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t120;
t51 = t108 * t57 + t109 * t56;
t107 = -g(3) + qJDD(2);
t123 = (m(4) + m(5)) * t107;
t122 = mrSges(5,2) * t73 + Ifges(5,4) * t103 + t102 * Ifges(5,6);
t121 = -t108 * t56 + t109 * t57;
t118 = -mrSges(5,1) * t73 + mrSges(5,3) * t71 + Ifges(5,2) * t103;
t117 = -mrSges(4,2) * t76 + t118 + qJ(4) * (-t102 * mrSges(5,1) + t124) + pkin(3) * t119 + mrSges(4,1) * t75 + Ifges(4,3) * t103;
t116 = mrSges(3,1) * t82 - mrSges(3,2) * t83 + Ifges(3,3) * qJDD(1) + pkin(2) * t59 + t117;
t115 = mrSges(2,1) * t90 - mrSges(2,2) * t91 + Ifges(2,3) * qJDD(1) + pkin(1) * t51 + t116;
t61 = -mrSges(4,3) * t75 + Ifges(4,5) * t103 - t102 * Ifges(4,6) + (-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t107 + t122;
t60 = t126 + mrSges(4,3) * t76 + (Ifges(4,6) - Ifges(5,6)) * t103 + (Ifges(5,4) + Ifges(4,5)) * t102 + (-m(5) * pkin(3) + t125) * t107;
t52 = mrSges(3,2) * t107 - mrSges(3,3) * t82 + Ifges(3,5) * qJDD(1) - t114 * Ifges(3,6) - pkin(5) * t59 - t110 * t60 + t112 * t61;
t49 = -mrSges(3,1) * t107 + mrSges(3,3) * t83 + t114 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t123 + pkin(5) * t120 + t110 * t61 + t112 * t60;
t48 = m(2) * t91 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t121;
t47 = m(2) * t90 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) + t51;
t46 = -mrSges(2,2) * g(3) - mrSges(2,3) * t90 + Ifges(2,5) * qJDD(1) - t114 * Ifges(2,6) - qJ(2) * t51 - t108 * t49 + t109 * t52;
t45 = Ifges(2,6) * qJDD(1) + t114 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t91 + t108 * t52 + t109 * t49 - pkin(1) * (m(3) * t107 + t123) + qJ(2) * t121;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t113 * t46 - t111 * t45 - pkin(4) * (t111 * t48 + t113 * t47), t46, t52, t61, -mrSges(5,3) * t107 + t122; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t111 * t46 + t113 * t45 + pkin(4) * (-t111 * t47 + t113 * t48), t45, t49, t60, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t115, t115, t116, t117, mrSges(5,1) * t107 - t102 * Ifges(5,4) + Ifges(5,6) * t103 - t126;];
m_new  = t1;
