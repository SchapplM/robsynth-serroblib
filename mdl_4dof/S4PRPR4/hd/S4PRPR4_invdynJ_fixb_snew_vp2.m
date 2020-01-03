% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPR4
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:58
% EndTime: 2019-12-31 16:21:58
% DurationCPUTime: 0.14s
% Computational Cost: add. (294->74), mult. (498->100), div. (0->0), fcn. (250->6), ass. (0->36)
t182 = -pkin(2) - pkin(5);
t167 = sin(pkin(6));
t168 = cos(pkin(6));
t162 = t167 * g(1) - t168 * g(2);
t163 = -t168 * g(1) - t167 * g(2);
t170 = sin(qJ(2));
t172 = cos(qJ(2));
t181 = t170 * t162 + t172 * t163;
t169 = sin(qJ(4));
t180 = qJD(2) * t169;
t171 = cos(qJ(4));
t179 = qJD(2) * t171;
t178 = qJD(2) * qJD(4);
t177 = t172 * t162 - t170 * t163;
t173 = qJD(2) ^ 2;
t175 = -t173 * qJ(3) + qJDD(3) - t177;
t150 = t182 * qJDD(2) + t175;
t166 = -g(3) + qJDD(1);
t147 = t171 * t150 - t169 * t166;
t159 = (t169 * mrSges(5,1) + t171 * mrSges(5,2)) * qJD(2);
t161 = t171 * qJDD(2) - t169 * t178;
t164 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t180;
t145 = m(5) * t147 + qJDD(4) * mrSges(5,1) - t161 * mrSges(5,3) + qJD(4) * t164 - t159 * t179;
t148 = t169 * t150 + t171 * t166;
t160 = -t169 * qJDD(2) - t171 * t178;
t165 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t179;
t146 = m(5) * t148 - qJDD(4) * mrSges(5,2) + t160 * mrSges(5,3) - qJD(4) * t165 - t159 * t180;
t176 = t171 * t145 + t169 * t146;
t174 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t181;
t155 = Ifges(5,5) * qJD(4) + (t171 * Ifges(5,1) - t169 * Ifges(5,4)) * qJD(2);
t154 = Ifges(5,6) * qJD(4) + (t171 * Ifges(5,4) - t169 * Ifges(5,2)) * qJD(2);
t152 = -qJDD(2) * pkin(2) + t175;
t151 = t173 * pkin(2) - t174;
t149 = t182 * t173 + t174;
t144 = m(4) * t152 + qJDD(2) * mrSges(4,2) - t173 * mrSges(4,3) + t176;
t1 = [-t169 * t145 + t171 * t146 + (m(2) + m(3) + m(4)) * t166; mrSges(3,1) * t177 - mrSges(3,2) * t181 + mrSges(4,2) * t152 - mrSges(4,3) * t151 + t171 * (mrSges(5,2) * t149 - mrSges(5,3) * t147 + Ifges(5,1) * t161 + Ifges(5,4) * t160 + Ifges(5,5) * qJDD(4) - qJD(4) * t154) - t169 * (-mrSges(5,1) * t149 + mrSges(5,3) * t148 + Ifges(5,4) * t161 + Ifges(5,2) * t160 + Ifges(5,6) * qJDD(4) + qJD(4) * t155) - pkin(5) * t176 - pkin(2) * t144 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (-m(4) * t151 + m(5) * t149 - t160 * mrSges(5,1) + t173 * mrSges(4,2) + t161 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + (t164 * t169 + t165 * t171) * qJD(2)) * qJ(3); t144; mrSges(5,1) * t147 - mrSges(5,2) * t148 + Ifges(5,5) * t161 + Ifges(5,6) * t160 + Ifges(5,3) * qJDD(4) + (t171 * t154 + t169 * t155) * qJD(2);];
tauJ = t1;
