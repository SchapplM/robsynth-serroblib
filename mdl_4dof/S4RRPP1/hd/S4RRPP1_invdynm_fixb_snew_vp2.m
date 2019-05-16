% Calculate vector of cutting torques with Newton-Euler for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-05-04 19:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:19:11
% EndTime: 2019-05-04 19:19:11
% DurationCPUTime: 0.46s
% Computational Cost: add. (6422->120), mult. (8288->138), div. (0->0), fcn. (3974->6), ass. (0->50)
t104 = (qJD(1) + qJD(2));
t102 = t104 ^ 2;
t103 = qJDD(1) + qJDD(2);
t107 = sin(pkin(6));
t108 = cos(pkin(6));
t109 = sin(qJ(2));
t111 = cos(qJ(2));
t110 = sin(qJ(1));
t112 = cos(qJ(1));
t89 = t110 * g(1) - t112 * g(2);
t86 = qJDD(1) * pkin(1) + t89;
t113 = qJD(1) ^ 2;
t90 = -t112 * g(1) - t110 * g(2);
t87 = -t113 * pkin(1) + t90;
t81 = -t109 * t87 + t111 * t86;
t78 = t103 * pkin(2) + t81;
t82 = t109 * t86 + t111 * t87;
t79 = -t102 * pkin(2) + t82;
t75 = t107 * t78 + t108 * t79;
t70 = -t102 * pkin(3) + t103 * qJ(4) + (2 * qJD(4) * t104) + t75;
t125 = mrSges(5,2) * t70;
t124 = -mrSges(4,1) - mrSges(5,1);
t123 = m(5) * t70 + t103 * mrSges(5,3);
t63 = m(4) * t75 - t103 * mrSges(4,2) + t124 * t102 + t123;
t74 = -t107 * t79 + t108 * t78;
t72 = -t103 * pkin(3) - t102 * qJ(4) + qJDD(4) - t74;
t118 = -m(5) * t72 + t103 * mrSges(5,1) + t102 * mrSges(5,3);
t65 = m(4) * t74 + t103 * mrSges(4,1) - t102 * mrSges(4,2) + t118;
t58 = t107 * t63 + t108 * t65;
t55 = m(3) * t81 + t103 * mrSges(3,1) - t102 * mrSges(3,2) + t58;
t120 = -t107 * t65 + t108 * t63;
t56 = m(3) * t82 - t102 * mrSges(3,1) - t103 * mrSges(3,2) + t120;
t49 = t109 * t56 + t111 * t55;
t106 = -g(3) + qJDD(3);
t122 = (m(4) + m(5)) * t106;
t121 = mrSges(5,2) * t72 + Ifges(5,4) * t103 + t102 * Ifges(5,6);
t119 = -t109 * t55 + t111 * t56;
t117 = -mrSges(5,1) * t72 + mrSges(5,3) * t70 + Ifges(5,2) * t103;
t116 = -mrSges(4,2) * t75 + t117 + qJ(4) * (-t102 * mrSges(5,1) + t123) + pkin(3) * t118 + mrSges(4,1) * t74 + Ifges(4,3) * t103;
t115 = mrSges(3,1) * t81 - mrSges(3,2) * t82 + Ifges(3,3) * t103 + pkin(2) * t58 + t116;
t114 = mrSges(2,1) * t89 - mrSges(2,2) * t90 + Ifges(2,3) * qJDD(1) + pkin(1) * t49 + t115;
t60 = -mrSges(4,3) * t74 + Ifges(4,5) * t103 - t102 * Ifges(4,6) + (-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t106 + t121;
t59 = t125 + mrSges(4,3) * t75 + (Ifges(4,6) - Ifges(5,6)) * t103 + (Ifges(5,4) + Ifges(4,5)) * t102 + (-m(5) * pkin(3) + t124) * t106;
t51 = -mrSges(3,2) * g(3) - mrSges(3,3) * t81 + Ifges(3,5) * t103 - t102 * Ifges(3,6) - qJ(3) * t58 - t107 * t59 + t108 * t60;
t50 = mrSges(3,1) * g(3) + mrSges(3,3) * t82 + t102 * Ifges(3,5) + Ifges(3,6) * t103 - pkin(2) * t122 + qJ(3) * t120 + t107 * t60 + t108 * t59;
t47 = m(2) * t90 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t119;
t46 = m(2) * t89 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) + t49;
t45 = -mrSges(2,2) * g(3) - mrSges(2,3) * t89 + Ifges(2,5) * qJDD(1) - t113 * Ifges(2,6) - pkin(5) * t49 - t109 * t50 + t111 * t51;
t44 = Ifges(2,6) * qJDD(1) + t113 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t90 + t109 * t51 + t111 * t50 - pkin(1) * (-m(3) * g(3) + t122) + pkin(5) * t119;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t112 * t45 - t110 * t44 - pkin(4) * (t110 * t47 + t112 * t46), t45, t51, t60, -mrSges(5,3) * t106 + t121; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t110 * t45 + t112 * t44 + pkin(4) * (-t110 * t46 + t112 * t47), t44, t50, t59, t117; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t114, t114, t115, t116, mrSges(5,1) * t106 - t102 * Ifges(5,4) + Ifges(5,6) * t103 - t125;];
m_new  = t1;
