% Calculate vector of cutting torques with Newton-Euler for
% S4RPRP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-05-04 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:15:07
% EndTime: 2019-05-04 19:15:07
% DurationCPUTime: 0.26s
% Computational Cost: add. (2902->119), mult. (3814->128), div. (0->0), fcn. (1058->4), ass. (0->50)
t126 = -pkin(1) - pkin(2);
t101 = g(3) + qJDD(4);
t125 = m(5) * t101;
t124 = mrSges(2,1) + mrSges(3,1);
t123 = (-mrSges(4,2) - mrSges(5,2));
t104 = sin(qJ(3));
t106 = cos(qJ(3));
t109 = qJD(1) ^ 2;
t105 = sin(qJ(1));
t107 = cos(qJ(1));
t82 = -t107 * g(1) - t105 * g(2);
t117 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t82;
t72 = t126 * t109 + t117;
t81 = t105 * g(1) - t107 * g(2);
t116 = -t109 * qJ(2) + qJDD(2) - t81;
t75 = t126 * qJDD(1) + t116;
t67 = -t104 * t72 + t106 * t75;
t93 = -qJDD(1) + qJDD(3);
t64 = t93 * pkin(3) + t67;
t122 = m(5) * t64 + t93 * mrSges(5,1);
t68 = t104 * t75 + t106 * t72;
t94 = -qJD(1) + qJD(3);
t92 = t94 ^ 2;
t65 = -t92 * pkin(3) + t68;
t121 = mrSges(5,3) * t65 + (t92 * Ifges(5,5)) + Ifges(5,6) * t93;
t120 = mrSges(5,1) * t64 - mrSges(5,2) * t65 + Ifges(5,3) * t93;
t56 = m(4) * t67 + t93 * mrSges(4,1) + (t123 * t92) + t122;
t63 = m(5) * t65;
t57 = m(4) * t68 + t63 + t123 * t93 + (-mrSges(4,1) - mrSges(5,1)) * t92;
t51 = -t104 * t56 + t106 * t57;
t119 = mrSges(5,2) * t101 - mrSges(5,3) * t64 + Ifges(5,5) * t93;
t50 = t104 * t57 + t106 * t56;
t76 = -t109 * pkin(1) + t117;
t118 = m(3) * t76 + qJDD(1) * mrSges(3,3) + t51;
t78 = -qJDD(1) * pkin(1) + t116;
t115 = -m(3) * t78 + qJDD(1) * mrSges(3,1) + t109 * mrSges(3,3) - t50;
t52 = Ifges(4,6) * t93 + (t92 * Ifges(4,5)) - mrSges(4,1) * g(3) + mrSges(4,3) * t68 + qJ(4) * (-t92 * mrSges(5,1) - t93 * mrSges(5,2) + t63) + (-m(5) * pkin(3) - mrSges(5,1)) * t101 + t121;
t58 = -t92 * mrSges(5,2) + t122;
t54 = mrSges(4,2) * g(3) - mrSges(4,3) * t67 + Ifges(4,5) * t93 - qJ(4) * t58 + (-Ifges(4,6) - Ifges(5,6)) * t92 + t119;
t114 = mrSges(3,2) * t78 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t109 * Ifges(3,6) - pkin(5) * t50 - t104 * t52 + t106 * t54;
t113 = mrSges(4,1) * t67 - mrSges(4,2) * t68 + Ifges(4,3) * t93 + pkin(3) * t58 + t120;
t112 = mrSges(3,2) * t76 - pkin(2) * (-m(4) * g(3) - t125) - pkin(5) * t51 - t104 * t54 - t106 * t52;
t111 = -mrSges(3,1) * t78 + mrSges(3,3) * t76 + Ifges(3,2) * qJDD(1) - pkin(2) * t50 - t113;
t110 = -mrSges(2,2) * t82 + t111 + qJ(2) * (-t109 * mrSges(3,1) + t118) + pkin(1) * t115 + mrSges(2,1) * t81 + Ifges(2,3) * qJDD(1);
t80 = -t125 + (-m(3) - m(4)) * g(3);
t47 = m(2) * t81 + qJDD(1) * mrSges(2,1) - t109 * mrSges(2,2) + t115;
t46 = m(2) * t82 - qJDD(1) * mrSges(2,2) - t124 * t109 + t118;
t45 = -mrSges(2,2) * g(3) - mrSges(2,3) * t81 + Ifges(2,5) * qJDD(1) - t109 * Ifges(2,6) - qJ(2) * t80 + t114;
t44 = mrSges(2,3) * t82 - pkin(1) * t80 + (Ifges(3,4) + Ifges(2,5)) * t109 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t124 * g(3) + t112;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t107 * t45 - t105 * t44 - pkin(4) * (t105 * t46 + t107 * t47), t45, t114, t54, -t92 * Ifges(5,6) + t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t105 * t45 + t107 * t44 + pkin(4) * (-t105 * t47 + t107 * t46), t44, t111, t52, -mrSges(5,1) * t101 + t121; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t110, t110, -mrSges(3,1) * g(3) - t109 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t112, t113, t120;];
m_new  = t1;
