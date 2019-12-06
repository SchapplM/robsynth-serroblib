% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:06
% EndTime: 2019-12-05 18:53:08
% DurationCPUTime: 0.74s
% Computational Cost: add. (5612->193), mult. (7623->255), div. (0->0), fcn. (5228->10), ass. (0->78)
t264 = qJD(1) + qJD(2);
t266 = sin(qJ(4));
t267 = sin(qJ(3));
t271 = cos(qJ(4));
t272 = cos(qJ(3));
t247 = (t266 * t267 - t271 * t272) * t264;
t287 = t264 * t267;
t286 = t264 * t272;
t269 = sin(qJ(1));
t274 = cos(qJ(1));
t282 = t269 * g(1) - t274 * g(2);
t253 = qJDD(1) * pkin(1) + t282;
t281 = -t274 * g(1) - t269 * g(2);
t254 = -qJD(1) ^ 2 * pkin(1) + t281;
t268 = sin(qJ(2));
t273 = cos(qJ(2));
t238 = t268 * t253 + t273 * t254;
t234 = -t267 * g(3) + t272 * t238;
t260 = t264 ^ 2;
t228 = (-t260 * t272 ^ 2 - qJD(3) ^ 2) * pkin(2) + t234;
t233 = -t272 * g(3) - t267 * t238;
t278 = (t260 * t267 * t272 + qJDD(3)) * pkin(2) + t233;
t211 = t271 * t228 + t266 * t278;
t237 = -t273 * t253 + t268 * t254;
t262 = qJDD(1) + qJDD(2);
t284 = qJD(3) * t264;
t283 = t267 * t284;
t251 = t272 * t262 - t283;
t224 = (-t251 + t283) * pkin(2) + t237;
t265 = sin(qJ(5));
t270 = cos(qJ(5));
t207 = -t265 * t211 + t270 * t224;
t250 = t267 * t262 + t272 * t284;
t223 = -t247 * qJD(4) + t271 * t250 + t266 * t251;
t248 = (t266 * t272 + t267 * t271) * t264;
t263 = qJD(3) + qJD(4);
t239 = -t265 * t248 + t270 * t263;
t261 = qJDD(3) + qJDD(4);
t213 = t239 * qJD(5) + t270 * t223 + t265 * t261;
t240 = t270 * t248 + t265 * t263;
t218 = -t239 * mrSges(6,1) + t240 * mrSges(6,2);
t222 = -t248 * qJD(4) - t266 * t250 + t271 * t251;
t221 = qJDD(5) - t222;
t243 = qJD(5) + t247;
t226 = -t243 * mrSges(6,2) + t239 * mrSges(6,3);
t205 = m(6) * t207 + t221 * mrSges(6,1) - t213 * mrSges(6,3) - t240 * t218 + t243 * t226;
t208 = t270 * t211 + t265 * t224;
t212 = -t240 * qJD(5) - t265 * t223 + t270 * t261;
t227 = t243 * mrSges(6,1) - t240 * mrSges(6,3);
t206 = m(6) * t208 - t221 * mrSges(6,2) + t212 * mrSges(6,3) + t239 * t218 - t243 * t227;
t232 = t247 * mrSges(5,1) + t248 * mrSges(5,2);
t242 = t263 * mrSges(5,1) - t248 * mrSges(5,3);
t197 = m(5) * t211 - t261 * mrSges(5,2) + t222 * mrSges(5,3) - t265 * t205 + t270 * t206 - t247 * t232 - t263 * t242;
t210 = t266 * t228 - t271 * t278;
t241 = -t263 * mrSges(5,2) - t247 * mrSges(5,3);
t204 = t261 * mrSges(5,1) + t212 * mrSges(6,1) - t213 * mrSges(6,2) - t223 * mrSges(5,3) + t239 * t226 - t240 * t227 - t248 * t232 + t263 * t241 + (-m(5) - m(6)) * t210;
t285 = t266 * t197 + t271 * t204;
t214 = Ifges(6,5) * t240 + Ifges(6,6) * t239 + Ifges(6,3) * t243;
t216 = Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t243;
t201 = -mrSges(6,1) * t210 + mrSges(6,3) * t208 + Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t221 - t240 * t214 + t243 * t216;
t215 = Ifges(6,4) * t240 + Ifges(6,2) * t239 + Ifges(6,6) * t243;
t202 = mrSges(6,2) * t210 - mrSges(6,3) * t207 + Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t221 + t239 * t214 - t243 * t215;
t229 = Ifges(5,5) * t248 - Ifges(5,6) * t247 + Ifges(5,3) * t263;
t230 = Ifges(5,4) * t248 - Ifges(5,2) * t247 + Ifges(5,6) * t263;
t195 = mrSges(5,2) * t224 + mrSges(5,3) * t210 + Ifges(5,1) * t223 + Ifges(5,4) * t222 + Ifges(5,5) * t261 - t265 * t201 + t270 * t202 - t247 * t229 - t263 * t230;
t231 = Ifges(5,1) * t248 - Ifges(5,4) * t247 + Ifges(5,5) * t263;
t276 = mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t213 + Ifges(6,6) * t212 + Ifges(6,3) * t221 + t240 * t215 - t239 * t216;
t198 = -mrSges(5,1) * t224 + mrSges(5,3) * t211 + Ifges(5,4) * t223 + Ifges(5,2) * t222 + Ifges(5,6) * t261 - t248 * t229 + t263 * t231 - t276;
t244 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t267 + Ifges(4,6) * t272) * t264;
t245 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t267 + Ifges(4,2) * t272) * t264;
t246 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t267 + Ifges(4,4) * t272) * t264;
t275 = m(5) * t224 - t222 * mrSges(5,1) + t223 * mrSges(5,2) + t270 * t205 + t265 * t206 + t247 * t241 + t248 * t242;
t279 = -mrSges(3,2) * t238 + t272 * (-mrSges(4,1) * t237 + mrSges(4,3) * t234 + Ifges(4,4) * t250 + Ifges(4,2) * t251 + Ifges(4,6) * qJDD(3) - pkin(2) * t275 + qJD(3) * t246 + t266 * t195 + t271 * t198 - t244 * t287) + t267 * (mrSges(4,2) * t237 - mrSges(4,3) * t233 + Ifges(4,1) * t250 + Ifges(4,4) * t251 + Ifges(4,5) * qJDD(3) - qJD(3) * t245 + t271 * t195 - t266 * t198 + t244 * t286) - mrSges(3,1) * t237 + Ifges(3,3) * t262;
t277 = -mrSges(5,1) * t210 - mrSges(5,2) * t211 + Ifges(5,5) * t223 + Ifges(5,6) * t222 + Ifges(5,3) * t261 + t270 * t201 + t265 * t202 + t248 * t230 + t247 * t231;
t256 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t286;
t255 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t287;
t249 = (-mrSges(4,1) * t272 + mrSges(4,2) * t267) * t264;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t282 - mrSges(2,2) * t281 + pkin(1) * (t268 * (m(3) * t238 - t262 * mrSges(3,2) - t260 * mrSges(3,1) + t272 * (m(4) * t234 - qJDD(3) * mrSges(4,2) + t251 * mrSges(4,3) - qJD(3) * t255 + t271 * t197 - t266 * t204 + t249 * t286) - t267 * (m(4) * t233 + qJDD(3) * mrSges(4,1) - t250 * mrSges(4,3) + qJD(3) * t256 - t249 * t287 + t285)) + t273 * (t262 * mrSges(3,1) + t251 * mrSges(4,1) - t260 * mrSges(3,2) - t250 * mrSges(4,2) + (-t255 * t267 + t256 * t272) * t264 + (-m(3) - m(4)) * t237 - t275)) + t279; t279; t277 + pkin(2) * t285 + Ifges(4,3) * qJDD(3) + (t267 * t245 - t272 * t246) * t264 + mrSges(4,1) * t233 - mrSges(4,2) * t234 + Ifges(4,5) * t250 + Ifges(4,6) * t251; t277; t276;];
tauJ = t1;
