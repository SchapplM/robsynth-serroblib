% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRR2
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:13
% EndTime: 2019-12-31 17:23:13
% DurationCPUTime: 0.40s
% Computational Cost: add. (3516->149), mult. (4505->199), div. (0->0), fcn. (2547->8), ass. (0->66)
t254 = qJD(1) + qJD(2);
t250 = t254 ^ 2;
t277 = pkin(3) * t250;
t259 = sin(qJ(1));
t263 = cos(qJ(1));
t272 = t259 * g(1) - g(2) * t263;
t242 = qJDD(1) * pkin(1) + t272;
t269 = -g(1) * t263 - t259 * g(2);
t243 = -qJD(1) ^ 2 * pkin(1) + t269;
t258 = sin(qJ(2));
t262 = cos(qJ(2));
t226 = t258 * t242 + t262 * t243;
t252 = qJDD(1) + qJDD(2);
t223 = -pkin(2) * t250 + pkin(6) * t252 + t226;
t257 = sin(qJ(3));
t276 = t223 * t257;
t275 = t254 * t257;
t261 = cos(qJ(3));
t274 = t254 * t261;
t273 = qJD(3) * t254;
t237 = t252 * t257 + t261 * t273;
t205 = qJDD(3) * pkin(3) - pkin(7) * t237 - t276 + (pkin(7) * t273 + t257 * t277 - g(3)) * t261;
t215 = -g(3) * t257 + t261 * t223;
t238 = t252 * t261 - t257 * t273;
t246 = qJD(3) * pkin(3) - pkin(7) * t275;
t255 = t261 ^ 2;
t206 = pkin(7) * t238 - qJD(3) * t246 - t255 * t277 + t215;
t256 = sin(qJ(4));
t260 = cos(qJ(4));
t203 = t205 * t260 - t206 * t256;
t232 = (-t256 * t257 + t260 * t261) * t254;
t213 = qJD(4) * t232 + t237 * t260 + t238 * t256;
t233 = (t256 * t261 + t257 * t260) * t254;
t221 = -mrSges(5,1) * t232 + mrSges(5,2) * t233;
t253 = qJD(3) + qJD(4);
t227 = -mrSges(5,2) * t253 + mrSges(5,3) * t232;
t251 = qJDD(3) + qJDD(4);
t200 = m(5) * t203 + mrSges(5,1) * t251 - mrSges(5,3) * t213 - t221 * t233 + t227 * t253;
t204 = t205 * t256 + t206 * t260;
t212 = -qJD(4) * t233 - t237 * t256 + t238 * t260;
t228 = mrSges(5,1) * t253 - mrSges(5,3) * t233;
t201 = m(5) * t204 - mrSges(5,2) * t251 + mrSges(5,3) * t212 + t221 * t232 - t228 * t253;
t193 = t260 * t200 + t256 * t201;
t214 = -g(3) * t261 - t276;
t236 = (-mrSges(4,1) * t261 + mrSges(4,2) * t257) * t254;
t244 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t275;
t245 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t274;
t270 = -t200 * t256 + t260 * t201;
t271 = -(m(4) * t214 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t237 + qJD(3) * t245 - t236 * t275 + t193) * t257 + t261 * (m(4) * t215 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t238 - qJD(3) * t244 + t236 * t274 + t270);
t225 = t242 * t262 - t258 * t243;
t268 = -pkin(2) * t252 - t225;
t207 = t246 * t275 - pkin(3) * t238 + (-pkin(7) * t255 - pkin(6)) * t250 + t268;
t216 = Ifges(5,5) * t233 + Ifges(5,6) * t232 + Ifges(5,3) * t253;
t218 = Ifges(5,1) * t233 + Ifges(5,4) * t232 + Ifges(5,5) * t253;
t194 = -mrSges(5,1) * t207 + mrSges(5,3) * t204 + Ifges(5,4) * t213 + Ifges(5,2) * t212 + Ifges(5,6) * t251 - t216 * t233 + t218 * t253;
t217 = Ifges(5,4) * t233 + Ifges(5,2) * t232 + Ifges(5,6) * t253;
t195 = mrSges(5,2) * t207 - mrSges(5,3) * t203 + Ifges(5,1) * t213 + Ifges(5,4) * t212 + Ifges(5,5) * t251 + t216 * t232 - t217 * t253;
t222 = -pkin(6) * t250 + t268;
t229 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t257 + Ifges(4,6) * t261) * t254;
t230 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t257 + Ifges(4,2) * t261) * t254;
t231 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t257 + Ifges(4,4) * t261) * t254;
t266 = m(5) * t207 - t212 * mrSges(5,1) + mrSges(5,2) * t213 - t232 * t227 + t228 * t233;
t264 = -m(4) * t222 + t238 * mrSges(4,1) - mrSges(4,2) * t237 - t244 * t275 + t245 * t274 - t266;
t267 = -mrSges(3,2) * t226 + t261 * (-mrSges(4,1) * t222 + mrSges(4,3) * t215 + Ifges(4,4) * t237 + Ifges(4,2) * t238 + Ifges(4,6) * qJDD(3) - pkin(3) * t266 + pkin(7) * t270 + qJD(3) * t231 + t260 * t194 + t256 * t195 - t229 * t275) + t257 * (mrSges(4,2) * t222 - mrSges(4,3) * t214 + Ifges(4,1) * t237 + Ifges(4,4) * t238 + Ifges(4,5) * qJDD(3) - pkin(7) * t193 - qJD(3) * t230 - t194 * t256 + t195 * t260 + t229 * t274) + pkin(6) * t271 + pkin(2) * t264 + mrSges(3,1) * t225 + Ifges(3,3) * t252;
t265 = mrSges(5,1) * t203 - mrSges(5,2) * t204 + Ifges(5,5) * t213 + Ifges(5,6) * t212 + Ifges(5,3) * t251 + t233 * t217 - t218 * t232;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t272 - mrSges(2,2) * t269 + pkin(1) * (t258 * (m(3) * t226 - mrSges(3,1) * t250 - mrSges(3,2) * t252 + t271) + t262 * (m(3) * t225 + mrSges(3,1) * t252 - mrSges(3,2) * t250 + t264)) + t267; t267; mrSges(4,1) * t214 - mrSges(4,2) * t215 + Ifges(4,5) * t237 + Ifges(4,6) * t238 + Ifges(4,3) * qJDD(3) + pkin(3) * t193 + (t230 * t257 - t231 * t261) * t254 + t265; t265;];
tauJ = t1;
