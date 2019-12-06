% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:31
% EndTime: 2019-12-05 14:59:32
% DurationCPUTime: 0.20s
% Computational Cost: add. (844->81), mult. (1282->116), div. (0->0), fcn. (918->10), ass. (0->45)
t179 = sin(pkin(7));
t182 = cos(pkin(7));
t173 = -t182 * g(1) - t179 * g(2);
t176 = -g(3) + qJDD(1);
t178 = sin(pkin(8));
t181 = cos(pkin(8));
t163 = t181 * t173 + t178 * t176;
t172 = -t179 * g(1) + t182 * g(2) + qJDD(2);
t177 = sin(pkin(9));
t180 = cos(pkin(9));
t160 = t180 * t163 + t177 * t172;
t190 = -t178 * t173 + t181 * t176;
t162 = qJDD(3) - t190;
t184 = sin(qJ(4));
t186 = cos(qJ(4));
t157 = t186 * t160 + t184 * t162;
t183 = sin(qJ(5));
t194 = qJD(4) * t183;
t185 = cos(qJ(5));
t193 = qJD(4) * t185;
t192 = qJD(4) * qJD(5);
t187 = qJD(4) ^ 2;
t155 = -t187 * pkin(4) + qJDD(4) * pkin(6) + t157;
t159 = t177 * t163 - t180 * t172;
t152 = -t183 * t155 + t185 * t159;
t169 = (-t185 * mrSges(6,1) + t183 * mrSges(6,2)) * qJD(4);
t170 = t183 * qJDD(4) + t185 * t192;
t175 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t193;
t150 = m(6) * t152 + qJDD(5) * mrSges(6,1) - t170 * mrSges(6,3) + qJD(5) * t175 - t169 * t194;
t153 = t185 * t155 + t183 * t159;
t171 = t185 * qJDD(4) - t183 * t192;
t174 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t194;
t151 = m(6) * t153 - qJDD(5) * mrSges(6,2) + t171 * mrSges(6,3) - qJD(5) * t174 + t169 * t193;
t191 = -t183 * t150 + t185 * t151;
t156 = -t184 * t160 + t186 * t162;
t147 = m(5) * t157 - t187 * mrSges(5,1) - qJDD(4) * mrSges(5,2) + t191;
t154 = -qJDD(4) * pkin(4) - t187 * pkin(6) - t156;
t188 = -m(6) * t154 + t171 * mrSges(6,1) - t170 * mrSges(6,2) - t174 * t194 + t175 * t193;
t148 = m(5) * t156 + qJDD(4) * mrSges(5,1) - t187 * mrSges(5,2) + t188;
t189 = m(4) * t162 + t184 * t147 + t186 * t148;
t166 = Ifges(6,5) * qJD(5) + (t183 * Ifges(6,1) + t185 * Ifges(6,4)) * qJD(4);
t165 = Ifges(6,6) * qJD(5) + (t183 * Ifges(6,4) + t185 * Ifges(6,2)) * qJD(4);
t146 = -t185 * t150 - t183 * t151 + (-m(4) - m(5)) * t159;
t145 = m(4) * t160 + t186 * t147 - t184 * t148;
t1 = [m(2) * t176 + t178 * (m(3) * t163 + t180 * t145 - t177 * t146) + t181 * (m(3) * t190 - t189); m(3) * t172 + t177 * t145 + t180 * t146; t189; Ifges(5,3) * qJDD(4) + mrSges(5,1) * t156 - mrSges(5,2) * t157 + t183 * (mrSges(6,2) * t154 - mrSges(6,3) * t152 + Ifges(6,1) * t170 + Ifges(6,4) * t171 + Ifges(6,5) * qJDD(5) - qJD(5) * t165) + t185 * (-mrSges(6,1) * t154 + mrSges(6,3) * t153 + Ifges(6,4) * t170 + Ifges(6,2) * t171 + Ifges(6,6) * qJDD(5) + qJD(5) * t166) + pkin(4) * t188 + pkin(6) * t191; mrSges(6,1) * t152 - mrSges(6,2) * t153 + Ifges(6,5) * t170 + Ifges(6,6) * t171 + Ifges(6,3) * qJDD(5) + (t183 * t165 - t185 * t166) * qJD(4);];
tauJ = t1;
