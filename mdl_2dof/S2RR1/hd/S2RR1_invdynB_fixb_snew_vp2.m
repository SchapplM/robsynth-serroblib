% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S2RR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynB_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynB_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynB_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynB_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_invdynB_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_invdynB_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_invdynB_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:07
% EndTime: 2020-01-03 11:19:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (345->84), mult. (664->119), div. (0->0), fcn. (312->4), ass. (0->33)
t109 = sin(qJ(2));
t118 = qJD(1) * t109;
t111 = cos(qJ(2));
t117 = qJD(1) * t111;
t116 = qJD(1) * qJD(2);
t110 = sin(qJ(1));
t112 = cos(qJ(1));
t107 = -t112 * g(1) + t110 * g(3);
t100 = (mrSges(3,1) * t111 - mrSges(3,2) * t109) * qJD(1);
t102 = -t109 * qJDD(1) - t111 * t116;
t105 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t117;
t106 = -t110 * g(1) - t112 * g(3);
t99 = -qJDD(1) * pkin(1) + t106;
t94 = t111 * g(2) - t109 * t99;
t92 = m(3) * t94 + qJDD(2) * mrSges(3,1) - t102 * mrSges(3,3) + qJD(2) * t105 + t100 * t118;
t103 = -t111 * qJDD(1) + t109 * t116;
t104 = qJD(2) * mrSges(3,1) + mrSges(3,3) * t118;
t95 = t109 * g(2) + t111 * t99;
t93 = m(3) * t95 - qJDD(2) * mrSges(3,2) + t103 * mrSges(3,3) - qJD(2) * t104 - t100 * t117;
t115 = -t109 * t92 + t111 * t93;
t114 = -t109 * t93 - t111 * t92;
t113 = qJD(1) ^ 2;
t101 = -t113 * pkin(1) + t107;
t98 = Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t109 - Ifges(3,4) * t111) * qJD(1);
t97 = Ifges(3,6) * qJD(2) + (-Ifges(3,4) * t109 - Ifges(3,2) * t111) * qJD(1);
t96 = Ifges(3,3) * qJD(2) + (-Ifges(3,5) * t109 - Ifges(3,6) * t111) * qJD(1);
t90 = m(2) * t107 + m(3) * t101 + qJDD(1) * mrSges(2,1) - t103 * mrSges(3,1) - t113 * mrSges(2,2) + t102 * mrSges(3,2) + (-t104 * t109 + t105 * t111) * qJD(1);
t89 = mrSges(3,2) * t101 - mrSges(3,3) * t94 + Ifges(3,1) * t102 + Ifges(3,4) * t103 + Ifges(3,5) * qJDD(2) - qJD(2) * t97 - t96 * t117;
t88 = -mrSges(3,1) * t101 + mrSges(3,3) * t95 + Ifges(3,4) * t102 + Ifges(3,2) * t103 + Ifges(3,6) * qJDD(2) + qJD(2) * t98 + t96 * t118;
t87 = mrSges(2,1) * g(2) + mrSges(3,1) * t94 - mrSges(3,2) * t95 + mrSges(2,3) * t106 + t113 * Ifges(2,5) + Ifges(3,5) * t102 + Ifges(2,6) * qJDD(1) + Ifges(3,6) * t103 + Ifges(3,3) * qJDD(2) + (-t109 * t97 + t111 * t98) * qJD(1);
t86 = m(2) * t106 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t115;
t85 = -mrSges(2,2) * g(2) - mrSges(2,3) * t107 + Ifges(2,5) * qJDD(1) - t113 * Ifges(2,6) + pkin(1) * t114 - t109 * t88 + t111 * t89;
t1 = [-m(1) * g(1) + t110 * t86 + t112 * t90; (-m(1) - m(2)) * g(2) + t114; -m(1) * g(3) - t110 * t90 + t112 * t86; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t110 * t85 + t112 * t87; mrSges(1,1) * g(3) + mrSges(2,1) * t107 - mrSges(2,2) * t106 - mrSges(1,3) * g(1) + Ifges(2,3) * qJDD(1) - pkin(1) * t115 - t109 * t89 - t111 * t88; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t110 * t87 + t112 * t85;];
tauB = t1;
