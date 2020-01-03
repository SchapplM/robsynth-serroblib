% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:48
% EndTime: 2019-12-31 16:27:49
% DurationCPUTime: 0.41s
% Computational Cost: add. (474->110), mult. (902->139), div. (0->0), fcn. (440->6), ass. (0->48)
t223 = sin(qJ(3));
t225 = cos(qJ(3));
t241 = Ifges(4,4) - Ifges(5,5);
t249 = t223 * (Ifges(4,1) + Ifges(5,1)) + t225 * t241;
t248 = t223 * t241 + t225 * (Ifges(4,2) + Ifges(5,3));
t240 = Ifges(5,4) + Ifges(4,5);
t246 = Ifges(5,6) - Ifges(4,6);
t244 = t223 * (-t248 * qJD(2) + t246 * qJD(3)) + t225 * (t249 * qJD(2) + t240 * qJD(3));
t228 = qJD(2) ^ 2;
t221 = sin(pkin(6));
t222 = cos(pkin(6));
t210 = t221 * g(1) - t222 * g(2);
t211 = -t222 * g(1) - t221 * g(2);
t224 = sin(qJ(2));
t226 = cos(qJ(2));
t231 = t226 * t210 - t224 * t211;
t191 = -qJDD(2) * pkin(2) - t228 * pkin(5) - t231;
t232 = qJD(2) * qJD(3);
t206 = t223 * qJDD(2) + t225 * t232;
t207 = t225 * qJDD(2) - t223 * t232;
t184 = -t207 * pkin(3) - t206 * qJ(4) + (-0.2e1 * qJD(4) * t223 + (pkin(3) * t223 - qJ(4) * t225) * qJD(3)) * qJD(2) + t191;
t243 = m(5) * t184;
t242 = mrSges(4,3) + mrSges(5,2);
t220 = -g(3) + qJDD(1);
t238 = t225 * t220;
t235 = t224 * t210 + t226 * t211;
t192 = -t228 * pkin(2) + qJDD(2) * pkin(5) + t235;
t189 = t225 * t192 + t223 * t220;
t234 = qJD(2) * t223;
t233 = qJD(2) * t225;
t203 = (-t225 * pkin(3) - t223 * qJ(4)) * qJD(2);
t227 = qJD(3) ^ 2;
t186 = -t227 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t203 * t233 + t189;
t204 = (-t225 * mrSges(5,1) - t223 * mrSges(5,3)) * qJD(2);
t213 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t234;
t230 = m(5) * t186 + qJDD(3) * mrSges(5,3) + qJD(3) * t213 + t204 * t233;
t187 = -qJDD(3) * pkin(3) - t227 * qJ(4) - t238 + qJDD(4) + (qJD(2) * t203 + t192) * t223;
t215 = mrSges(5,2) * t233 + qJD(3) * mrSges(5,3);
t229 = -m(5) * t187 + qJDD(3) * mrSges(5,1) + qJD(3) * t215;
t214 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t233;
t212 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t234;
t205 = (-t225 * mrSges(4,1) + t223 * mrSges(4,2)) * qJD(2);
t188 = -t223 * t192 + t238;
t183 = t206 * mrSges(5,2) + t204 * t234 - t229;
t182 = t243 - t207 * mrSges(5,1) - t206 * mrSges(5,3) + (-t213 * t223 - t215 * t225) * qJD(2);
t181 = m(4) * t188 + qJDD(3) * mrSges(4,1) + qJD(3) * t214 - t242 * t206 + (-t204 - t205) * t234 + t229;
t180 = m(4) * t189 - qJDD(3) * mrSges(4,2) - qJD(3) * t212 + t205 * t233 + t242 * t207 + t230;
t1 = [t223 * t180 + t225 * t181 + (m(2) + m(3)) * t220; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t231 - mrSges(3,2) * t235 + t223 * (mrSges(4,2) * t191 + mrSges(5,2) * t187 - mrSges(4,3) * t188 - mrSges(5,3) * t184 - qJ(4) * t182) + t225 * (-mrSges(4,1) * t191 - mrSges(5,1) * t184 + mrSges(5,2) * t186 + mrSges(4,3) * t189 - pkin(3) * t182) + pkin(5) * (t225 * t180 - t223 * t181) + (t223 * t240 - t225 * t246) * qJDD(3) + t244 * qJD(3) + t248 * t207 + t249 * t206 + (-m(4) * t191 - t243 + (mrSges(4,1) + mrSges(5,1)) * t207 + (-mrSges(4,2) + mrSges(5,3)) * t206 + ((t214 + t215) * t225 + (-t212 + t213) * t223) * qJD(2)) * pkin(2); mrSges(4,1) * t188 - mrSges(4,2) * t189 - mrSges(5,1) * t187 + mrSges(5,3) * t186 - pkin(3) * t183 + qJ(4) * t230 + (qJ(4) * mrSges(5,2) - t246) * t207 + t240 * t206 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) - t244 * qJD(2); t183;];
tauJ = t1;
