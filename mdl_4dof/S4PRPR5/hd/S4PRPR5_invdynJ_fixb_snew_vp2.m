% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:01
% EndTime: 2019-12-31 16:23:01
% DurationCPUTime: 0.15s
% Computational Cost: add. (532->82), mult. (870->114), div. (0->0), fcn. (494->8), ass. (0->40)
t183 = sin(pkin(6));
t185 = cos(pkin(6));
t176 = -t185 * g(1) - t183 * g(2);
t181 = -g(3) + qJDD(1);
t187 = sin(qJ(2));
t189 = cos(qJ(2));
t164 = -t187 * t176 + t189 * t181;
t162 = qJDD(2) * pkin(2) + t164;
t165 = t189 * t176 + t187 * t181;
t190 = qJD(2) ^ 2;
t163 = -t190 * pkin(2) + t165;
t182 = sin(pkin(7));
t184 = cos(pkin(7));
t159 = t182 * t162 + t184 * t163;
t157 = -t190 * pkin(3) + qJDD(2) * pkin(5) + t159;
t175 = -t183 * g(1) + t185 * g(2) + qJDD(3);
t186 = sin(qJ(4));
t188 = cos(qJ(4));
t154 = -t186 * t157 + t188 * t175;
t172 = (-t188 * mrSges(5,1) + t186 * mrSges(5,2)) * qJD(2);
t193 = qJD(2) * qJD(4);
t173 = t186 * qJDD(2) + t188 * t193;
t194 = qJD(2) * t188;
t178 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t194;
t195 = qJD(2) * t186;
t152 = m(5) * t154 + qJDD(4) * mrSges(5,1) - t173 * mrSges(5,3) + qJD(4) * t178 - t172 * t195;
t155 = t188 * t157 + t186 * t175;
t174 = t188 * qJDD(2) - t186 * t193;
t177 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t195;
t153 = m(5) * t155 - qJDD(4) * mrSges(5,2) + t174 * mrSges(5,3) - qJD(4) * t177 + t172 * t194;
t192 = -t186 * t152 + t188 * t153;
t148 = m(4) * t159 - t190 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t192;
t158 = t184 * t162 - t182 * t163;
t156 = -qJDD(2) * pkin(3) - t190 * pkin(5) - t158;
t191 = -m(5) * t156 + t174 * mrSges(5,1) - t173 * mrSges(5,2) - t177 * t195 + t178 * t194;
t150 = m(4) * t158 + qJDD(2) * mrSges(4,1) - t190 * mrSges(4,2) + t191;
t196 = t182 * t148 + t184 * t150;
t168 = Ifges(5,5) * qJD(4) + (t186 * Ifges(5,1) + t188 * Ifges(5,4)) * qJD(2);
t167 = Ifges(5,6) * qJD(4) + (t186 * Ifges(5,4) + t188 * Ifges(5,2)) * qJD(2);
t1 = [m(2) * t181 + t187 * (m(3) * t165 - t190 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t184 * t148 - t182 * t150) + t189 * (m(3) * t164 + qJDD(2) * mrSges(3,1) - t190 * mrSges(3,2) + t196); mrSges(3,1) * t164 - mrSges(3,2) * t165 + mrSges(4,1) * t158 - mrSges(4,2) * t159 + t186 * (mrSges(5,2) * t156 - mrSges(5,3) * t154 + Ifges(5,1) * t173 + Ifges(5,4) * t174 + Ifges(5,5) * qJDD(4) - qJD(4) * t167) + t188 * (-mrSges(5,1) * t156 + mrSges(5,3) * t155 + Ifges(5,4) * t173 + Ifges(5,2) * t174 + Ifges(5,6) * qJDD(4) + qJD(4) * t168) + pkin(3) * t191 + pkin(5) * t192 + pkin(2) * t196 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); m(4) * t175 + t188 * t152 + t186 * t153; mrSges(5,1) * t154 - mrSges(5,2) * t155 + Ifges(5,5) * t173 + Ifges(5,6) * t174 + Ifges(5,3) * qJDD(4) + (t186 * t167 - t188 * t168) * qJD(2);];
tauJ = t1;
