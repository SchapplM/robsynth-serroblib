% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:34
% EndTime: 2019-12-05 15:16:34
% DurationCPUTime: 0.42s
% Computational Cost: add. (2221->148), mult. (4045->199), div. (0->0), fcn. (2650->10), ass. (0->68)
t265 = sin(pkin(8));
t267 = cos(pkin(8));
t255 = -t267 * g(1) - t265 * g(2);
t263 = -g(3) + qJDD(1);
t264 = sin(pkin(9));
t266 = cos(pkin(9));
t242 = t266 * t255 + t264 * t263;
t254 = -t265 * g(1) + t267 * g(2) + qJDD(2);
t270 = sin(qJ(3));
t273 = cos(qJ(3));
t232 = t273 * t242 + t270 * t254;
t274 = qJD(3) ^ 2;
t227 = -t274 * pkin(3) + qJDD(3) * pkin(6) + t232;
t241 = t264 * t255 - t266 * t263;
t269 = sin(qJ(4));
t272 = cos(qJ(4));
t217 = -t269 * t227 + t272 * t241;
t282 = qJD(3) * qJD(4);
t281 = t272 * t282;
t252 = t269 * qJDD(3) + t281;
t214 = (-t252 + t281) * pkin(7) + (t269 * t272 * t274 + qJDD(4)) * pkin(4) + t217;
t218 = t272 * t227 + t269 * t241;
t253 = t272 * qJDD(3) - t269 * t282;
t283 = t269 * qJD(3);
t258 = qJD(4) * pkin(4) - pkin(7) * t283;
t262 = t272 ^ 2;
t215 = -t262 * t274 * pkin(4) + t253 * pkin(7) - qJD(4) * t258 + t218;
t268 = sin(qJ(5));
t271 = cos(qJ(5));
t212 = t271 * t214 - t268 * t215;
t246 = (-t269 * t268 + t272 * t271) * qJD(3);
t225 = t246 * qJD(5) + t271 * t252 + t268 * t253;
t247 = (t272 * t268 + t269 * t271) * qJD(3);
t234 = -t246 * mrSges(6,1) + t247 * mrSges(6,2);
t261 = qJD(4) + qJD(5);
t239 = -t261 * mrSges(6,2) + t246 * mrSges(6,3);
t260 = qJDD(4) + qJDD(5);
t209 = m(6) * t212 + t260 * mrSges(6,1) - t225 * mrSges(6,3) - t247 * t234 + t261 * t239;
t213 = t268 * t214 + t271 * t215;
t224 = -t247 * qJD(5) - t268 * t252 + t271 * t253;
t240 = t261 * mrSges(6,1) - t247 * mrSges(6,3);
t210 = m(6) * t213 - t260 * mrSges(6,2) + t224 * mrSges(6,3) + t246 * t234 - t261 * t240;
t202 = t271 * t209 + t268 * t210;
t284 = qJD(3) * t272;
t251 = (-t272 * mrSges(5,1) + t269 * mrSges(5,2)) * qJD(3);
t257 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t284;
t200 = m(5) * t217 + qJDD(4) * mrSges(5,1) - t252 * mrSges(5,3) + qJD(4) * t257 - t251 * t283 + t202;
t256 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t283;
t279 = -t268 * t209 + t271 * t210;
t201 = m(5) * t218 - qJDD(4) * mrSges(5,2) + t253 * mrSges(5,3) - qJD(4) * t256 + t251 * t284 + t279;
t280 = -t269 * t200 + t272 * t201;
t231 = -t270 * t242 + t273 * t254;
t278 = -qJDD(3) * pkin(3) - t231;
t216 = t258 * t283 - t253 * pkin(4) + (-pkin(7) * t262 - pkin(6)) * t274 + t278;
t277 = m(6) * t216 - t224 * mrSges(6,1) + t225 * mrSges(6,2) - t246 * t239 + t247 * t240;
t229 = Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t261;
t230 = Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t261;
t276 = mrSges(6,1) * t212 - mrSges(6,2) * t213 + Ifges(6,5) * t225 + Ifges(6,6) * t224 + Ifges(6,3) * t260 + t247 * t229 - t246 * t230;
t226 = -t274 * pkin(6) + t278;
t275 = -m(5) * t226 + t253 * mrSges(5,1) - t252 * mrSges(5,2) - t256 * t283 + t257 * t284 - t277;
t245 = Ifges(5,5) * qJD(4) + (t269 * Ifges(5,1) + t272 * Ifges(5,4)) * qJD(3);
t244 = Ifges(5,6) * qJD(4) + (t269 * Ifges(5,4) + t272 * Ifges(5,2)) * qJD(3);
t228 = Ifges(6,5) * t247 + Ifges(6,6) * t246 + Ifges(6,3) * t261;
t205 = m(4) * t231 + qJDD(3) * mrSges(4,1) - t274 * mrSges(4,2) + t275;
t204 = mrSges(6,2) * t216 - mrSges(6,3) * t212 + Ifges(6,1) * t225 + Ifges(6,4) * t224 + Ifges(6,5) * t260 + t246 * t228 - t261 * t229;
t203 = -mrSges(6,1) * t216 + mrSges(6,3) * t213 + Ifges(6,4) * t225 + Ifges(6,2) * t224 + Ifges(6,6) * t260 - t247 * t228 + t261 * t230;
t198 = m(4) * t232 - t274 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t280;
t1 = [m(2) * t263 + t264 * (m(3) * t242 + t273 * t198 - t270 * t205) + t266 * (-t272 * t200 - t269 * t201 + (-m(3) - m(4)) * t241); m(3) * t254 + t270 * t198 + t273 * t205; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t231 - mrSges(4,2) * t232 + t269 * (mrSges(5,2) * t226 - mrSges(5,3) * t217 + Ifges(5,1) * t252 + Ifges(5,4) * t253 + Ifges(5,5) * qJDD(4) - pkin(7) * t202 - qJD(4) * t244 - t268 * t203 + t271 * t204) + t272 * (-mrSges(5,1) * t226 + mrSges(5,3) * t218 + Ifges(5,4) * t252 + Ifges(5,2) * t253 + Ifges(5,6) * qJDD(4) - pkin(4) * t277 + pkin(7) * t279 + qJD(4) * t245 + t271 * t203 + t268 * t204) + pkin(3) * t275 + pkin(6) * t280; mrSges(5,1) * t217 - mrSges(5,2) * t218 + Ifges(5,5) * t252 + Ifges(5,6) * t253 + Ifges(5,3) * qJDD(4) + pkin(4) * t202 + (t269 * t244 - t272 * t245) * qJD(3) + t276; t276;];
tauJ = t1;
