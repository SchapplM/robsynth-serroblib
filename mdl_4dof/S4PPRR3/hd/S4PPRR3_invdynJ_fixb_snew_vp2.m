% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:24
% EndTime: 2019-12-31 16:17:24
% DurationCPUTime: 0.10s
% Computational Cost: add. (239->68), mult. (406->95), div. (0->0), fcn. (220->6), ass. (0->31)
t139 = sin(pkin(6));
t140 = cos(pkin(6));
t134 = -t139 * g(1) + t140 * g(2) + qJDD(2);
t135 = -t140 * g(1) - t139 * g(2);
t142 = sin(qJ(3));
t144 = cos(qJ(3));
t123 = t142 * t134 + t144 * t135;
t141 = sin(qJ(4));
t150 = qJD(3) * t141;
t143 = cos(qJ(4));
t149 = qJD(3) * t143;
t148 = qJD(3) * qJD(4);
t145 = qJD(3) ^ 2;
t121 = -t145 * pkin(3) + qJDD(3) * pkin(5) + t123;
t138 = g(3) - qJDD(1);
t118 = -t141 * t121 + t143 * t138;
t131 = (-t143 * mrSges(5,1) + t141 * mrSges(5,2)) * qJD(3);
t132 = t141 * qJDD(3) + t143 * t148;
t137 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t149;
t116 = m(5) * t118 + qJDD(4) * mrSges(5,1) - t132 * mrSges(5,3) + qJD(4) * t137 - t131 * t150;
t119 = t143 * t121 + t141 * t138;
t133 = t143 * qJDD(3) - t141 * t148;
t136 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t150;
t117 = m(5) * t119 - qJDD(4) * mrSges(5,2) + t133 * mrSges(5,3) - qJD(4) * t136 + t131 * t149;
t147 = -t141 * t116 + t143 * t117;
t122 = t144 * t134 - t142 * t135;
t120 = -qJDD(3) * pkin(3) - t145 * pkin(5) - t122;
t146 = -m(5) * t120 + t133 * mrSges(5,1) - t132 * mrSges(5,2) - t136 * t150 + t137 * t149;
t126 = Ifges(5,5) * qJD(4) + (t141 * Ifges(5,1) + t143 * Ifges(5,4)) * qJD(3);
t125 = Ifges(5,6) * qJD(4) + (t141 * Ifges(5,4) + t143 * Ifges(5,2)) * qJD(3);
t1 = [-t143 * t116 - t141 * t117 + (-m(2) - m(3) - m(4)) * t138; m(3) * t134 + t142 * (m(4) * t123 - t145 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t147) + t144 * (m(4) * t122 + qJDD(3) * mrSges(4,1) - t145 * mrSges(4,2) + t146); Ifges(4,3) * qJDD(3) + mrSges(4,1) * t122 - mrSges(4,2) * t123 + t141 * (mrSges(5,2) * t120 - mrSges(5,3) * t118 + Ifges(5,1) * t132 + Ifges(5,4) * t133 + Ifges(5,5) * qJDD(4) - qJD(4) * t125) + t143 * (-mrSges(5,1) * t120 + mrSges(5,3) * t119 + Ifges(5,4) * t132 + Ifges(5,2) * t133 + Ifges(5,6) * qJDD(4) + qJD(4) * t126) + pkin(3) * t146 + pkin(5) * t147; mrSges(5,1) * t118 - mrSges(5,2) * t119 + Ifges(5,5) * t132 + Ifges(5,6) * t133 + Ifges(5,3) * qJDD(4) + (t141 * t125 - t143 * t126) * qJD(3);];
tauJ = t1;
