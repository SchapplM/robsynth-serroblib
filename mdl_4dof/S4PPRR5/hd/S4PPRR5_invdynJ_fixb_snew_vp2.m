% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:43
% EndTime: 2019-12-31 16:19:43
% DurationCPUTime: 0.13s
% Computational Cost: add. (266->67), mult. (433->95), div. (0->0), fcn. (226->6), ass. (0->31)
t149 = sin(pkin(6));
t150 = cos(pkin(6));
t143 = -g(1) * t149 + g(2) * t150 + qJDD(2);
t148 = -g(3) + qJDD(1);
t152 = sin(qJ(3));
t154 = cos(qJ(3));
t133 = t152 * t143 + t154 * t148;
t151 = sin(qJ(4));
t160 = qJD(3) * t151;
t153 = cos(qJ(4));
t159 = qJD(3) * t153;
t158 = qJD(3) * qJD(4);
t155 = qJD(3) ^ 2;
t131 = -pkin(3) * t155 + qJDD(3) * pkin(5) + t133;
t144 = -g(1) * t150 - g(2) * t149;
t128 = -t131 * t151 + t144 * t153;
t129 = t131 * t153 + t144 * t151;
t140 = (-t153 * mrSges(5,1) + t151 * mrSges(5,2)) * qJD(3);
t141 = qJDD(3) * t151 + t153 * t158;
t142 = qJDD(3) * t153 - t151 * t158;
t145 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t160;
t146 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t159;
t157 = -(m(5) * t128 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t141 + qJD(4) * t146 - t140 * t160) * t151 + t153 * (m(5) * t129 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t142 - qJD(4) * t145 + t140 * t159);
t132 = t143 * t154 - t148 * t152;
t130 = -qJDD(3) * pkin(3) - pkin(5) * t155 - t132;
t156 = -m(5) * t130 + t142 * mrSges(5,1) - mrSges(5,2) * t141 - t145 * t160 + t146 * t159;
t136 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t151 + Ifges(5,4) * t153) * qJD(3);
t135 = Ifges(5,6) * qJD(4) + (t151 * Ifges(5,4) + t153 * Ifges(5,2)) * qJD(3);
t125 = m(4) * t132 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t155 + t156;
t124 = m(4) * t133 - mrSges(4,1) * t155 - qJDD(3) * mrSges(4,2) + t157;
t1 = [t124 * t154 - t125 * t152 + (m(2) + m(3)) * t148; m(3) * t143 + t124 * t152 + t125 * t154; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t132 - mrSges(4,2) * t133 + t151 * (mrSges(5,2) * t130 - mrSges(5,3) * t128 + Ifges(5,1) * t141 + Ifges(5,4) * t142 + Ifges(5,5) * qJDD(4) - qJD(4) * t135) + t153 * (-mrSges(5,1) * t130 + mrSges(5,3) * t129 + Ifges(5,4) * t141 + Ifges(5,2) * t142 + Ifges(5,6) * qJDD(4) + qJD(4) * t136) + pkin(3) * t156 + pkin(5) * t157; mrSges(5,1) * t128 - mrSges(5,2) * t129 + Ifges(5,5) * t141 + Ifges(5,6) * t142 + Ifges(5,3) * qJDD(4) + (t151 * t135 - t153 * t136) * qJD(3);];
tauJ = t1;
