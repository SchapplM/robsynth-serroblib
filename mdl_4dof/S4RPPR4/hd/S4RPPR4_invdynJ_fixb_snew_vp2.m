% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:50
% EndTime: 2019-12-31 16:38:51
% DurationCPUTime: 0.19s
% Computational Cost: add. (435->87), mult. (732->113), div. (0->0), fcn. (324->6), ass. (0->40)
t199 = -pkin(2) - pkin(5);
t185 = sin(qJ(1));
t187 = cos(qJ(1));
t195 = t185 * g(1) - t187 * g(2);
t170 = qJDD(1) * pkin(1) + t195;
t188 = qJD(1) ^ 2;
t194 = -t187 * g(1) - t185 * g(2);
t172 = -t188 * pkin(1) + t194;
t182 = sin(pkin(6));
t183 = cos(pkin(6));
t160 = t182 * t170 + t183 * t172;
t184 = sin(qJ(4));
t198 = qJD(1) * t184;
t186 = cos(qJ(4));
t197 = qJD(1) * t186;
t196 = qJD(1) * qJD(4);
t159 = t183 * t170 - t182 * t172;
t192 = -t188 * qJ(3) + qJDD(3) - t159;
t156 = t199 * qJDD(1) + t192;
t179 = -g(3) + qJDD(2);
t152 = t186 * t156 - t184 * t179;
t171 = (t184 * mrSges(5,1) + t186 * mrSges(5,2)) * qJD(1);
t174 = t186 * qJDD(1) - t184 * t196;
t175 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t198;
t150 = m(5) * t152 + qJDD(4) * mrSges(5,1) - t174 * mrSges(5,3) + qJD(4) * t175 - t171 * t197;
t153 = t184 * t156 + t186 * t179;
t173 = -t184 * qJDD(1) - t186 * t196;
t176 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t197;
t151 = m(5) * t153 - qJDD(4) * mrSges(5,2) + t173 * mrSges(5,3) - qJD(4) * t176 - t171 * t198;
t193 = t186 * t150 + t184 * t151;
t191 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t160;
t158 = -qJDD(1) * pkin(2) + t192;
t190 = m(4) * t158 - t188 * mrSges(4,3) + t193;
t155 = t199 * t188 + t191;
t157 = t188 * pkin(2) - t191;
t189 = -m(4) * t157 + m(5) * t155 - t173 * mrSges(5,1) + t188 * mrSges(4,2) + t174 * mrSges(5,2) + qJDD(1) * mrSges(4,3) + t175 * t198 + t176 * t197;
t166 = Ifges(5,5) * qJD(4) + (t186 * Ifges(5,1) - t184 * Ifges(5,4)) * qJD(1);
t165 = Ifges(5,6) * qJD(4) + (t186 * Ifges(5,4) - t184 * Ifges(5,2)) * qJD(1);
t149 = qJDD(1) * mrSges(4,2) + t190;
t1 = [mrSges(2,1) * t195 - mrSges(2,2) * t194 + pkin(1) * (t182 * (m(3) * t160 - t188 * mrSges(3,1) + t189) + t183 * (m(3) * t159 - t188 * mrSges(3,2) - t190)) - pkin(5) * t193 + mrSges(3,1) * t159 - mrSges(3,2) * t160 - pkin(2) * t149 + qJ(3) * t189 + mrSges(4,2) * t158 - mrSges(4,3) * t157 + t186 * (mrSges(5,2) * t155 - mrSges(5,3) * t152 + Ifges(5,1) * t174 + Ifges(5,4) * t173 + Ifges(5,5) * qJDD(4) - qJD(4) * t165) - t184 * (-mrSges(5,1) * t155 + mrSges(5,3) * t153 + Ifges(5,4) * t174 + Ifges(5,2) * t173 + Ifges(5,6) * qJDD(4) + qJD(4) * t166) + (pkin(1) * (-t182 * mrSges(3,2) + t183 * (mrSges(3,1) - mrSges(4,2))) + Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1); -t184 * t150 + t186 * t151 + (m(3) + m(4)) * t179; t149; mrSges(5,1) * t152 - mrSges(5,2) * t153 + Ifges(5,5) * t174 + Ifges(5,6) * t173 + Ifges(5,3) * qJDD(4) + (t186 * t165 + t184 * t166) * qJD(1);];
tauJ = t1;
