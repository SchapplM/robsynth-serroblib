% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S2RR3
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
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
% tauJB [(6+2)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S2RR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR3_invdynJB_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:25
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (477->61), mult. (700->77), div. (0->0), fcn. (332->4), ass. (0->29)
t135 = sin(qJ(1));
t137 = cos(qJ(1));
t126 = t135 * g(1) - t137 * g(2);
t123 = qJDD(1) * pkin(1) + t126;
t127 = -t137 * g(1) - t135 * g(2);
t138 = qJD(1) ^ 2;
t124 = -t138 * pkin(1) + t127;
t134 = sin(qJ(2));
t136 = cos(qJ(2));
t121 = t136 * t123 - t134 * t124;
t132 = qJD(1) + qJD(2);
t130 = t132 ^ 2;
t131 = qJDD(1) + qJDD(2);
t118 = m(3) * t121 + t131 * mrSges(3,1) - t130 * mrSges(3,2);
t122 = t134 * t123 + t136 * t124;
t119 = m(3) * t122 - t130 * mrSges(3,1) - t131 * mrSges(3,2);
t112 = t136 * t118 + t134 * t119;
t109 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t138 * mrSges(2,2) + t112;
t141 = -t134 * t118 + t136 * t119;
t110 = m(2) * t127 - t138 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t141;
t143 = t137 * t109 + t135 * t110;
t142 = -t135 * t109 + t137 * t110;
t140 = mrSges(3,1) * t121 - mrSges(3,2) * t122 + Ifges(3,3) * t131;
t139 = mrSges(2,1) * t126 - mrSges(2,2) * t127 + Ifges(2,3) * qJDD(1) + pkin(1) * t112 + t140;
t117 = -mrSges(3,2) * g(3) - mrSges(3,3) * t121 + Ifges(3,5) * t131 - t130 * Ifges(3,6);
t116 = mrSges(3,1) * g(3) + mrSges(3,3) * t122 + t130 * Ifges(3,5) + Ifges(3,6) * t131;
t105 = -mrSges(2,2) * g(3) - mrSges(2,3) * t126 + Ifges(2,5) * qJDD(1) - t138 * Ifges(2,6) - pkin(3) * t112 - t134 * t116 + t136 * t117;
t104 = Ifges(2,6) * qJDD(1) + t138 * Ifges(2,5) + mrSges(2,3) * t127 + t134 * t117 + t136 * t116 + pkin(3) * t141 + (pkin(1) * m(3) + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t142; -m(1) * g(2) + t143; (-m(1) - m(2) - m(3)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(2) * t143 - t135 * t104 + t137 * t105; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(2) * t142 + t137 * t104 + t135 * t105; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t139; t139; t140;];
tauJB = t1;
