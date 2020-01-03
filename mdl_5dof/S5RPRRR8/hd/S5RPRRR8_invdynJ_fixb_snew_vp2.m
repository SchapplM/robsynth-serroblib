% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:51
% EndTime: 2019-12-31 19:05:52
% DurationCPUTime: 0.57s
% Computational Cost: add. (5767->167), mult. (7442->216), div. (0->0), fcn. (3520->8), ass. (0->75)
t318 = -pkin(1) - pkin(2);
t292 = -qJD(1) + qJD(3);
t296 = sin(qJ(4));
t317 = t292 * t296;
t300 = cos(qJ(4));
t316 = t292 * t300;
t303 = qJD(1) ^ 2;
t298 = sin(qJ(1));
t302 = cos(qJ(1));
t311 = -t302 * g(1) - t298 * g(2);
t308 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t311;
t266 = t318 * t303 + t308;
t313 = t298 * g(1) - t302 * g(2);
t307 = -t303 * qJ(2) + qJDD(2) - t313;
t267 = t318 * qJDD(1) + t307;
t297 = sin(qJ(3));
t301 = cos(qJ(3));
t255 = t301 * t266 + t297 * t267;
t288 = t292 ^ 2;
t290 = -qJDD(1) + qJDD(3);
t253 = -t288 * pkin(3) + t290 * pkin(7) + t255;
t244 = t300 * g(3) - t296 * t253;
t315 = qJD(4) * t292;
t314 = t300 * t315;
t278 = t296 * t290 + t314;
t240 = (-t278 + t314) * pkin(8) + (t288 * t296 * t300 + qJDD(4)) * pkin(4) + t244;
t245 = t296 * g(3) + t300 * t253;
t279 = t300 * t290 - t296 * t315;
t282 = qJD(4) * pkin(4) - pkin(8) * t317;
t294 = t300 ^ 2;
t241 = -t294 * t288 * pkin(4) + t279 * pkin(8) - qJD(4) * t282 + t245;
t295 = sin(qJ(5));
t299 = cos(qJ(5));
t238 = t299 * t240 - t295 * t241;
t271 = (-t295 * t296 + t299 * t300) * t292;
t251 = t271 * qJD(5) + t299 * t278 + t295 * t279;
t272 = (t295 * t300 + t296 * t299) * t292;
t260 = -t271 * mrSges(6,1) + t272 * mrSges(6,2);
t291 = qJD(4) + qJD(5);
t261 = -t291 * mrSges(6,2) + t271 * mrSges(6,3);
t289 = qJDD(4) + qJDD(5);
t235 = m(6) * t238 + t289 * mrSges(6,1) - t251 * mrSges(6,3) - t272 * t260 + t291 * t261;
t239 = t295 * t240 + t299 * t241;
t250 = -t272 * qJD(5) - t295 * t278 + t299 * t279;
t262 = t291 * mrSges(6,1) - t272 * mrSges(6,3);
t236 = m(6) * t239 - t289 * mrSges(6,2) + t250 * mrSges(6,3) + t271 * t260 - t291 * t262;
t227 = t299 * t235 + t295 * t236;
t277 = (-mrSges(5,1) * t300 + mrSges(5,2) * t296) * t292;
t280 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t317;
t281 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t316;
t312 = -t295 * t235 + t299 * t236;
t224 = -t296 * (m(5) * t244 + qJDD(4) * mrSges(5,1) - t278 * mrSges(5,3) + qJD(4) * t281 - t277 * t317 + t227) + t300 * (m(5) * t245 - qJDD(4) * mrSges(5,2) + t279 * mrSges(5,3) - qJD(4) * t280 + t277 * t316 + t312);
t254 = -t297 * t266 + t301 * t267;
t223 = m(4) * t255 - t288 * mrSges(4,1) - t290 * mrSges(4,2) + t224;
t309 = -t290 * pkin(3) - t254;
t252 = -t288 * pkin(7) + t309;
t242 = t282 * t317 - t279 * pkin(4) + (-pkin(8) * t294 - pkin(7)) * t288 + t309;
t306 = m(6) * t242 - t250 * mrSges(6,1) + t251 * mrSges(6,2) - t271 * t261 + t272 * t262;
t231 = -m(5) * t252 + t279 * mrSges(5,1) - t278 * mrSges(5,2) - t280 * t317 + t281 * t316 - t306;
t230 = m(4) * t254 + t290 * mrSges(4,1) - t288 * mrSges(4,2) + t231;
t310 = t297 * t223 + t301 * t230;
t257 = Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t291;
t258 = Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t291;
t305 = mrSges(6,1) * t238 - mrSges(6,2) * t239 + Ifges(6,5) * t251 + Ifges(6,6) * t250 + Ifges(6,3) * t289 + t272 * t257 - t271 * t258;
t256 = Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t291;
t228 = -mrSges(6,1) * t242 + mrSges(6,3) * t239 + Ifges(6,4) * t251 + Ifges(6,2) * t250 + Ifges(6,6) * t289 - t272 * t256 + t291 * t258;
t229 = mrSges(6,2) * t242 - mrSges(6,3) * t238 + Ifges(6,1) * t251 + Ifges(6,4) * t250 + Ifges(6,5) * t289 + t271 * t256 - t291 * t257;
t268 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t296 + Ifges(5,6) * t300) * t292;
t269 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t296 + Ifges(5,2) * t300) * t292;
t270 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t296 + Ifges(5,4) * t300) * t292;
t304 = mrSges(4,1) * t254 - mrSges(4,2) * t255 + Ifges(4,3) * t290 + pkin(3) * t231 + pkin(7) * t224 + t300 * (-mrSges(5,1) * t252 + mrSges(5,3) * t245 + Ifges(5,4) * t278 + Ifges(5,2) * t279 + Ifges(5,6) * qJDD(4) - pkin(4) * t306 + pkin(8) * t312 + qJD(4) * t270 + t299 * t228 + t295 * t229 - t268 * t317) + t296 * (mrSges(5,2) * t252 - mrSges(5,3) * t244 + Ifges(5,1) * t278 + Ifges(5,4) * t279 + Ifges(5,5) * qJDD(4) - pkin(8) * t227 - qJD(4) * t269 - t295 * t228 + t299 * t229 + t268 * t316);
t276 = -qJDD(1) * pkin(1) + t307;
t273 = -t303 * pkin(1) + t308;
t221 = m(3) * t276 - qJDD(1) * mrSges(3,1) - t303 * mrSges(3,3) + t310;
t1 = [-pkin(1) * t221 + qJ(2) * (m(3) * t273 - t303 * mrSges(3,1) + t301 * t223 - t297 * t230) + mrSges(2,1) * t313 - mrSges(2,2) * t311 - pkin(2) * t310 - mrSges(3,1) * t276 + mrSges(3,3) * t273 + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) - t304; t221; t304; mrSges(5,1) * t244 - mrSges(5,2) * t245 + Ifges(5,5) * t278 + Ifges(5,6) * t279 + Ifges(5,3) * qJDD(4) + pkin(4) * t227 + (t269 * t296 - t270 * t300) * t292 + t305; t305;];
tauJ = t1;
