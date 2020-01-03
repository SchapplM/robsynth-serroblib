% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPR7
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:39
% EndTime: 2019-12-31 16:25:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (336->79), mult. (548->105), div. (0->0), fcn. (260->6), ass. (0->36)
t192 = -pkin(2) - pkin(5);
t177 = sin(pkin(6));
t178 = cos(pkin(6));
t169 = -t178 * g(1) - t177 * g(2);
t174 = -g(3) + qJDD(1);
t180 = sin(qJ(2));
t182 = cos(qJ(2));
t156 = t182 * t169 + t180 * t174;
t179 = sin(qJ(4));
t191 = qJD(2) * t179;
t181 = cos(qJ(4));
t190 = qJD(2) * t181;
t189 = qJD(2) * qJD(4);
t155 = -t180 * t169 + t182 * t174;
t183 = qJD(2) ^ 2;
t187 = -t183 * qJ(3) + qJDD(3) - t155;
t152 = t192 * qJDD(2) + t187;
t168 = -t177 * g(1) + t178 * g(2);
t148 = t181 * t152 - t179 * t168;
t149 = t179 * t152 + t181 * t168;
t165 = (t179 * mrSges(5,1) + t181 * mrSges(5,2)) * qJD(2);
t166 = -t179 * qJDD(2) - t181 * t189;
t167 = t181 * qJDD(2) - t179 * t189;
t170 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t191;
t171 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t190;
t188 = t181 * (m(5) * t148 + qJDD(4) * mrSges(5,1) - t167 * mrSges(5,3) + qJD(4) * t170 - t165 * t190) + t179 * (m(5) * t149 - qJDD(4) * mrSges(5,2) + t166 * mrSges(5,3) - qJD(4) * t171 - t165 * t191);
t186 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t156;
t154 = -qJDD(2) * pkin(2) + t187;
t185 = -m(4) * t154 + t183 * mrSges(4,3) - t188;
t151 = t192 * t183 + t186;
t153 = t183 * pkin(2) - t186;
t184 = -m(4) * t153 + m(5) * t151 - t166 * mrSges(5,1) + t183 * mrSges(4,2) + t167 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t170 * t191 + t171 * t190;
t159 = Ifges(5,5) * qJD(4) + (t181 * Ifges(5,1) - t179 * Ifges(5,4)) * qJD(2);
t158 = Ifges(5,6) * qJD(4) + (t181 * Ifges(5,4) - t179 * Ifges(5,2)) * qJD(2);
t145 = qJDD(2) * mrSges(4,2) - t185;
t1 = [m(2) * t174 + t180 * (m(3) * t156 - t183 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t184) + t182 * (m(3) * t155 - t183 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t185); mrSges(3,1) * t155 - mrSges(3,2) * t156 + mrSges(4,2) * t154 - mrSges(4,3) * t153 + t181 * (mrSges(5,2) * t151 - mrSges(5,3) * t148 + Ifges(5,1) * t167 + Ifges(5,4) * t166 + Ifges(5,5) * qJDD(4) - qJD(4) * t158) - t179 * (-mrSges(5,1) * t151 + mrSges(5,3) * t149 + Ifges(5,4) * t167 + Ifges(5,2) * t166 + Ifges(5,6) * qJDD(4) + qJD(4) * t159) - pkin(5) * t188 - pkin(2) * t145 + qJ(3) * t184 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t145; mrSges(5,1) * t148 - mrSges(5,2) * t149 + Ifges(5,5) * t167 + Ifges(5,6) * t166 + Ifges(5,3) * qJDD(4) + (t181 * t158 + t179 * t159) * qJD(2);];
tauJ = t1;
