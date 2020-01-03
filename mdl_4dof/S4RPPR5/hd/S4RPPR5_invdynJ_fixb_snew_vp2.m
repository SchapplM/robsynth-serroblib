% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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

function tauJ = S4RPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:45
% EndTime: 2019-12-31 16:39:45
% DurationCPUTime: 0.21s
% Computational Cost: add. (683->89), mult. (1136->117), div. (0->0), fcn. (406->6), ass. (0->42)
t190 = -pkin(1) - pkin(2);
t179 = qJD(1) ^ 2;
t176 = sin(qJ(1));
t178 = cos(qJ(1));
t184 = -g(1) * t178 - g(2) * t176;
t182 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t184;
t155 = t190 * t179 + t182;
t186 = g(1) * t176 - t178 * g(2);
t181 = -qJ(2) * t179 + qJDD(2) - t186;
t156 = t190 * qJDD(1) + t181;
t173 = sin(pkin(6));
t174 = cos(pkin(6));
t152 = t174 * t155 + t173 * t156;
t175 = sin(qJ(4));
t189 = qJD(1) * t175;
t177 = cos(qJ(4));
t188 = qJD(1) * t177;
t187 = qJD(1) * qJD(4);
t150 = -pkin(3) * t179 - qJDD(1) * pkin(5) + t152;
t172 = g(3) + qJDD(3);
t147 = -t150 * t175 + t172 * t177;
t164 = (t177 * mrSges(5,1) - t175 * mrSges(5,2)) * qJD(1);
t165 = -qJDD(1) * t175 - t177 * t187;
t168 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t188;
t145 = m(5) * t147 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t165 + qJD(4) * t168 + t164 * t189;
t148 = t150 * t177 + t172 * t175;
t166 = -qJDD(1) * t177 + t175 * t187;
t167 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t189;
t146 = m(5) * t148 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t166 - qJD(4) * t167 - t164 * t188;
t185 = -t175 * t145 + t177 * t146;
t142 = m(4) * t152 - mrSges(4,1) * t179 + qJDD(1) * mrSges(4,2) + t185;
t151 = -t155 * t173 + t156 * t174;
t149 = qJDD(1) * pkin(3) - pkin(5) * t179 - t151;
t180 = -m(5) * t149 + t166 * mrSges(5,1) - mrSges(5,2) * t165 + t167 * t189 - t168 * t188;
t143 = m(4) * t151 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t179 + t180;
t183 = t142 * t173 + t143 * t174;
t161 = (Ifges(5,5) * qJD(4)) + (-t175 * Ifges(5,1) - t177 * Ifges(5,4)) * qJD(1);
t160 = (Ifges(5,6) * qJD(4)) + (-t175 * Ifges(5,4) - t177 * Ifges(5,2)) * qJD(1);
t158 = -qJDD(1) * pkin(1) + t181;
t157 = -pkin(1) * t179 + t182;
t141 = m(3) * t158 - qJDD(1) * mrSges(3,1) - mrSges(3,3) * t179 + t183;
t1 = [-pkin(1) * t141 + qJ(2) * (m(3) * t157 - mrSges(3,1) * t179 + t142 * t174 - t143 * t173) + mrSges(2,1) * t186 - mrSges(2,2) * t184 - pkin(2) * t183 - mrSges(3,1) * t158 + mrSges(3,3) * t157 - mrSges(4,1) * t151 + mrSges(4,2) * t152 - t175 * (mrSges(5,2) * t149 - mrSges(5,3) * t147 + Ifges(5,1) * t165 + Ifges(5,4) * t166 + Ifges(5,5) * qJDD(4) - qJD(4) * t160) - t177 * (-mrSges(5,1) * t149 + mrSges(5,3) * t148 + Ifges(5,4) * t165 + Ifges(5,2) * t166 + Ifges(5,6) * qJDD(4) + qJD(4) * t161) - pkin(3) * t180 - pkin(5) * t185 + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); t141; m(4) * t172 + t145 * t177 + t146 * t175; mrSges(5,1) * t147 - mrSges(5,2) * t148 + Ifges(5,5) * t165 + Ifges(5,6) * t166 + Ifges(5,3) * qJDD(4) + (-t175 * t160 + t177 * t161) * qJD(1);];
tauJ = t1;
