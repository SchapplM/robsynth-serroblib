% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:35
% EndTime: 2019-12-31 16:51:35
% DurationCPUTime: 0.23s
% Computational Cost: add. (1032->92), mult. (1347->122), div. (0->0), fcn. (480->6), ass. (0->44)
t198 = -pkin(1) - pkin(2);
t180 = -qJD(1) + qJD(3);
t182 = sin(qJ(4));
t197 = t180 * t182;
t185 = cos(qJ(4));
t196 = t180 * t185;
t188 = qJD(1) ^ 2;
t184 = sin(qJ(1));
t187 = cos(qJ(1));
t193 = -t187 * g(1) - t184 * g(2);
t191 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t193;
t162 = t198 * t188 + t191;
t194 = t184 * g(1) - t187 * g(2);
t190 = -t188 * qJ(2) + qJDD(2) - t194;
t163 = t198 * qJDD(1) + t190;
t183 = sin(qJ(3));
t186 = cos(qJ(3));
t159 = t186 * t162 + t183 * t163;
t195 = qJD(4) * t180;
t178 = t180 ^ 2;
t179 = -qJDD(1) + qJDD(3);
t157 = -t178 * pkin(3) + t179 * pkin(6) + t159;
t154 = t185 * g(3) - t182 * t157;
t155 = t182 * g(3) + t185 * t157;
t171 = (-mrSges(5,1) * t185 + mrSges(5,2) * t182) * t180;
t172 = t182 * t179 + t185 * t195;
t173 = t185 * t179 - t182 * t195;
t174 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t197;
t175 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t196;
t147 = -t182 * (m(5) * t154 + qJDD(4) * mrSges(5,1) - t172 * mrSges(5,3) + qJD(4) * t175 - t171 * t197) + t185 * (m(5) * t155 - qJDD(4) * mrSges(5,2) + t173 * mrSges(5,3) - qJD(4) * t174 + t171 * t196);
t146 = m(4) * t159 - t178 * mrSges(4,1) - t179 * mrSges(4,2) + t147;
t158 = -t183 * t162 + t186 * t163;
t156 = -t179 * pkin(3) - t178 * pkin(6) - t158;
t151 = -m(5) * t156 + t173 * mrSges(5,1) - t172 * mrSges(5,2) - t174 * t197 + t175 * t196;
t150 = m(4) * t158 + t179 * mrSges(4,1) - t178 * mrSges(4,2) + t151;
t192 = t183 * t146 + t186 * t150;
t164 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t182 + Ifges(5,6) * t185) * t180;
t165 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t182 + Ifges(5,2) * t185) * t180;
t166 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t182 + Ifges(5,4) * t185) * t180;
t189 = mrSges(4,1) * t158 - mrSges(4,2) * t159 + Ifges(4,3) * t179 + pkin(3) * t151 + pkin(6) * t147 + t185 * (-mrSges(5,1) * t156 + mrSges(5,3) * t155 + Ifges(5,4) * t172 + Ifges(5,2) * t173 + Ifges(5,6) * qJDD(4) + qJD(4) * t166 - t164 * t197) + t182 * (mrSges(5,2) * t156 - mrSges(5,3) * t154 + Ifges(5,1) * t172 + Ifges(5,4) * t173 + Ifges(5,5) * qJDD(4) - qJD(4) * t165 + t164 * t196);
t170 = -qJDD(1) * pkin(1) + t190;
t167 = -t188 * pkin(1) + t191;
t145 = m(3) * t170 - qJDD(1) * mrSges(3,1) - t188 * mrSges(3,3) + t192;
t1 = [qJ(2) * (m(3) * t167 - t188 * mrSges(3,1) + t186 * t146 - t183 * t150) - pkin(1) * t145 + mrSges(2,1) * t194 - mrSges(2,2) * t193 - pkin(2) * t192 - mrSges(3,1) * t170 + mrSges(3,3) * t167 + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) - t189; t145; t189; mrSges(5,1) * t154 - mrSges(5,2) * t155 + Ifges(5,5) * t172 + Ifges(5,6) * t173 + Ifges(5,3) * qJDD(4) + (t165 * t182 - t166 * t185) * t180;];
tauJ = t1;
