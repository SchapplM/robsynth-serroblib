% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
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
% tauJB [(6+2)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S2RR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_invdynJB_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:08
% EndTime: 2020-01-03 11:19:08
% DurationCPUTime: 0.17s
% Computational Cost: add. (425->84), mult. (819->118), div. (0->0), fcn. (388->4), ass. (0->35)
t151 = sin(qJ(2));
t161 = qJD(1) * t151;
t153 = cos(qJ(2));
t160 = qJD(1) * t153;
t159 = qJD(1) * qJD(2);
t152 = sin(qJ(1));
t154 = cos(qJ(1));
t147 = -t154 * g(1) + t152 * g(3);
t146 = -t152 * g(1) - t154 * g(3);
t138 = -qJDD(1) * pkin(1) + t146;
t130 = t153 * g(2) - t151 * t138;
t139 = (mrSges(3,1) * t153 - mrSges(3,2) * t151) * qJD(1);
t141 = -t151 * qJDD(1) - t153 * t159;
t145 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t160;
t127 = m(3) * t130 + qJDD(2) * mrSges(3,1) - t141 * mrSges(3,3) + qJD(2) * t145 + t139 * t161;
t131 = t151 * g(2) + t153 * t138;
t142 = -t153 * qJDD(1) + t151 * t159;
t144 = qJD(2) * mrSges(3,1) + mrSges(3,3) * t161;
t128 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t142 * mrSges(3,3) - qJD(2) * t144 - t139 * t160;
t121 = -t151 * t127 + t153 * t128;
t158 = -t153 * t127 - t151 * t128;
t134 = Ifges(3,6) * qJD(2) + (-Ifges(3,4) * t151 - Ifges(3,2) * t153) * qJD(1);
t135 = Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t151 - Ifges(3,4) * t153) * qJD(1);
t157 = mrSges(3,1) * t130 - mrSges(3,2) * t131 + Ifges(3,5) * t141 + Ifges(3,6) * t142 + Ifges(3,3) * qJDD(2) - t134 * t161 + t135 * t160;
t133 = Ifges(3,3) * qJD(2) + (-Ifges(3,5) * t151 - Ifges(3,6) * t153) * qJD(1);
t155 = qJD(1) ^ 2;
t140 = -t155 * pkin(1) + t147;
t123 = -mrSges(3,1) * t140 + mrSges(3,3) * t131 + Ifges(3,4) * t141 + Ifges(3,2) * t142 + Ifges(3,6) * qJDD(2) + qJD(2) * t135 + t133 * t161;
t124 = mrSges(3,2) * t140 - mrSges(3,3) * t130 + Ifges(3,1) * t141 + Ifges(3,4) * t142 + Ifges(3,5) * qJDD(2) - qJD(2) * t134 - t133 * t160;
t156 = mrSges(2,1) * t147 - mrSges(2,2) * t146 + Ifges(2,3) * qJDD(1) - pkin(1) * t121 - t153 * t123 - t151 * t124;
t125 = m(2) * t147 + m(3) * t140 + qJDD(1) * mrSges(2,1) - t142 * mrSges(3,1) - t155 * mrSges(2,2) + t141 * mrSges(3,2) + (-t144 * t151 + t145 * t153) * qJD(1);
t122 = mrSges(2,1) * g(2) + mrSges(2,3) * t146 + t155 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + t157;
t120 = m(2) * t146 - t155 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t121;
t119 = -mrSges(2,2) * g(2) - mrSges(2,3) * t147 + Ifges(2,5) * qJDD(1) - t155 * Ifges(2,6) + pkin(1) * t158 - t151 * t123 + t153 * t124;
t1 = [-m(1) * g(1) + t152 * t120 + t154 * t125; (-m(1) - m(2)) * g(2) + t158; -m(1) * g(3) + t154 * t120 - t152 * t125; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t152 * t119 + t154 * t122; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t156; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t154 * t119 - t152 * t122; t156; t157;];
tauJB = t1;
