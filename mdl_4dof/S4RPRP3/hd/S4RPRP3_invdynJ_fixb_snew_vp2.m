% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:35
% EndTime: 2019-12-31 16:42:36
% DurationCPUTime: 0.58s
% Computational Cost: add. (675->126), mult. (1287->154), div. (0->0), fcn. (578->6), ass. (0->59)
t250 = cos(qJ(3));
t268 = Ifges(4,4) + Ifges(5,4);
t277 = t250 * t268;
t276 = Ifges(4,1) + Ifges(5,1);
t267 = Ifges(4,5) + Ifges(5,5);
t275 = Ifges(4,2) + Ifges(5,2);
t274 = Ifges(4,6) + Ifges(5,6);
t247 = cos(pkin(6));
t273 = pkin(1) * t247 + pkin(2);
t248 = sin(qJ(3));
t272 = -t267 * qJD(3) + (-t276 * t248 - t277) * qJD(1);
t252 = qJD(1) ^ 2;
t270 = pkin(3) * t252;
t269 = -mrSges(4,2) - mrSges(5,2);
t249 = sin(qJ(1));
t251 = cos(qJ(1));
t257 = t249 * g(1) - t251 * g(2);
t228 = qJDD(1) * pkin(1) + t257;
t255 = -t251 * g(1) - t249 * g(2);
t231 = -t252 * pkin(1) + t255;
t246 = sin(pkin(6));
t212 = t246 * t228 + t247 * t231;
t210 = -t252 * pkin(2) + qJDD(1) * pkin(5) + t212;
t245 = -g(3) + qJDD(2);
t207 = t250 * t210 + t248 * t245;
t266 = t274 * qJD(3) + (t268 * t248 + t250 * t275) * qJD(1);
t263 = t248 * qJD(1);
t236 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t263;
t265 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t263 - t236;
t264 = qJD(1) * t250;
t262 = qJD(1) * qJD(3);
t261 = qJD(1) * qJD(4);
t258 = t250 * t262;
t232 = t248 * qJDD(1) + t258;
t241 = t250 * t245;
t203 = qJDD(3) * pkin(3) + t241 + (-t232 + t258) * qJ(4) + (t250 * t270 - t210 - 0.2e1 * t261) * t248;
t238 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t264;
t260 = m(5) * t203 + qJDD(3) * mrSges(5,1) + qJD(3) * t238;
t233 = t250 * qJDD(1) - t248 * t262;
t235 = qJD(3) * pkin(3) - qJ(4) * t263;
t244 = t250 ^ 2;
t204 = t233 * qJ(4) - qJD(3) * t235 - t244 * t270 + 0.2e1 * t250 * t261 + t207;
t229 = (-t250 * mrSges(5,1) + t248 * mrSges(5,2)) * qJD(1);
t259 = m(5) * t204 + t233 * mrSges(5,3) + t229 * t264;
t211 = t247 * t228 - t246 * t231;
t254 = -qJDD(1) * pkin(2) - t211;
t205 = t235 * t263 - t233 * pkin(3) + qJDD(4) + (-qJ(4) * t244 - pkin(5)) * t252 + t254;
t256 = m(5) * t205 - t233 * mrSges(5,1) - t238 * t264;
t209 = -t252 * pkin(5) + t254;
t239 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t264;
t253 = -m(4) * t209 + t233 * mrSges(4,1) + t239 * t264 - t256;
t230 = (-t250 * mrSges(4,1) + t248 * mrSges(4,2)) * qJD(1);
t206 = -t248 * t210 + t241;
t200 = t232 * mrSges(5,2) + t236 * t263 + t256;
t199 = -t232 * mrSges(5,3) - t229 * t263 + t260;
t198 = m(4) * t207 + t233 * mrSges(4,3) + t265 * qJD(3) + t269 * qJDD(3) + t230 * t264 + t259;
t197 = m(4) * t206 + qJDD(3) * mrSges(4,1) + qJD(3) * t239 + (-mrSges(4,3) - mrSges(5,3)) * t232 + (-t229 - t230) * t263 + t260;
t196 = t250 * t198;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t257 - mrSges(2,2) * t255 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t211 - mrSges(3,2) * t212 + t250 * (-mrSges(4,1) * t209 - mrSges(5,1) * t205 + mrSges(4,3) * t207 + mrSges(5,3) * t204 - pkin(3) * t200 + qJ(4) * t259 + t275 * t233 + (-qJ(4) * mrSges(5,2) + t274) * qJDD(3) + (-qJ(4) * t236 - t272) * qJD(3)) + pkin(2) * t253 + pkin(5) * t196 + pkin(1) * (t246 * (m(3) * t212 - t252 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t196) + t247 * (m(3) * t211 + qJDD(1) * mrSges(3,1) - t252 * mrSges(3,2) + t253)) + (t273 * t269 + t277) * t232 + (mrSges(4,2) * t209 + mrSges(5,2) * t205 - mrSges(4,3) * t206 - mrSges(5,3) * t203 - qJ(4) * t199 + t268 * t233 + t276 * t232 + (-pkin(1) * t246 - pkin(5)) * t197 + t267 * qJDD(3) - t266 * qJD(3) + t273 * qJD(1) * t265) * t248; m(3) * t245 + t250 * t197 + t248 * t198; mrSges(4,1) * t206 + mrSges(5,1) * t203 - mrSges(4,2) * t207 - mrSges(5,2) * t204 + pkin(3) * t199 + t274 * t233 + t267 * t232 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t266 * t248 + t272 * t250) * qJD(1); t200;];
tauJ = t1;
