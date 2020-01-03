% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRP6
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:01
% EndTime: 2019-12-31 16:46:02
% DurationCPUTime: 0.54s
% Computational Cost: add. (553->119), mult. (1044->142), div. (0->0), fcn. (405->4), ass. (0->49)
t246 = sin(qJ(3));
t248 = cos(qJ(3));
t271 = Ifges(4,4) + Ifges(5,4);
t277 = t248 * (Ifges(4,1) + Ifges(5,1)) - t246 * t271;
t276 = t248 * t271 - t246 * (Ifges(4,2) + Ifges(5,2));
t270 = Ifges(4,5) + Ifges(5,5);
t269 = Ifges(4,6) + Ifges(5,6);
t273 = (t276 * qJD(1) + t269 * qJD(3)) * t248;
t247 = sin(qJ(1));
t249 = cos(qJ(1));
t255 = -t249 * g(1) - t247 * g(2);
t252 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t255;
t272 = -pkin(1) - pkin(5);
t250 = qJD(1) ^ 2;
t258 = t247 * g(1) - t249 * g(2);
t251 = -t250 * qJ(2) + qJDD(2) - t258;
t217 = t272 * qJDD(1) + t251;
t213 = -t248 * g(3) + t246 * t217;
t263 = qJD(1) * qJD(3);
t232 = -t246 * qJDD(1) - t248 * t263;
t264 = t248 * qJD(1);
t237 = qJD(3) * pkin(3) - qJ(4) * t264;
t245 = t246 ^ 2;
t261 = -0.2e1 * qJD(1) * qJD(4);
t209 = -t245 * t250 * pkin(3) + t232 * qJ(4) - qJD(3) * t237 + t246 * t261 + t213;
t268 = m(5) * t209 + t232 * mrSges(5,3);
t266 = t277 * qJD(1) + t270 * qJD(3);
t212 = t246 * g(3) + t248 * t217;
t265 = qJD(1) * t246;
t259 = t246 * t263;
t233 = t248 * qJDD(1) - t259;
t208 = t248 * t261 + (-t233 - t259) * qJ(4) + (-t246 * t248 * t250 + qJDD(3)) * pkin(3) + t212;
t235 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t265;
t260 = m(5) * t208 + qJDD(3) * mrSges(5,1) + qJD(3) * t235;
t230 = (t246 * mrSges(5,1) + t248 * mrSges(5,2)) * qJD(1);
t257 = qJD(1) * (-t230 - (t246 * mrSges(4,1) + t248 * mrSges(4,2)) * qJD(1));
t211 = t237 * t264 - t232 * pkin(3) + qJDD(4) + (-qJ(4) * t245 + t272) * t250 + t252;
t238 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t264;
t256 = m(5) * t211 + t233 * mrSges(5,2) + t235 * t265 + t238 * t264;
t236 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t265;
t239 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t264;
t254 = t248 * (m(4) * t212 + qJDD(3) * mrSges(4,1) + qJD(3) * t236 + (-mrSges(4,3) - mrSges(5,3)) * t233 + t248 * t257 + t260) + t246 * (m(4) * t213 + t232 * mrSges(4,3) + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t238 - t239) * qJD(3) + t246 * t257 + t268);
t219 = -qJDD(1) * pkin(1) + t251;
t218 = t250 * pkin(1) - t252;
t216 = t272 * t250 + t252;
t205 = -t233 * mrSges(5,3) - t230 * t264 + t260;
t204 = -t232 * mrSges(5,1) + t256;
t201 = m(3) * t219 + qJDD(1) * mrSges(3,2) - t250 * mrSges(3,3) + t254;
t1 = [mrSges(2,1) * t258 - mrSges(2,2) * t255 + mrSges(3,2) * t219 - mrSges(3,3) * t218 + t248 * (mrSges(4,2) * t216 + mrSges(5,2) * t211 - mrSges(4,3) * t212 - mrSges(5,3) * t208 - qJ(4) * t205) - t246 * (-mrSges(4,1) * t216 + mrSges(4,3) * t213 - mrSges(5,1) * t211 + mrSges(5,3) * t209 - pkin(3) * t204 + qJ(4) * (-t230 * t265 + t268)) - pkin(5) * t254 - pkin(1) * t201 + qJ(2) * (-m(3) * t218 + m(4) * t216 + t250 * mrSges(3,2) + t236 * t265 + t239 * t264 + t256) + (qJ(2) * mrSges(4,2) + t277) * t233 + (t248 * t270 - t246 * (-qJ(4) * mrSges(5,2) + t269)) * qJDD(3) + (qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t273 - t246 * (-qJ(4) * t238 + t266)) * qJD(3) + (qJ(2) * (-mrSges(4,1) - mrSges(5,1)) + t276) * t232; t201; mrSges(4,1) * t212 + mrSges(5,1) * t208 - mrSges(4,2) * t213 - mrSges(5,2) * t209 + pkin(3) * t205 + t270 * t233 + t269 * t232 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t266 * t246 + t273) * qJD(1); t204;];
tauJ = t1;
