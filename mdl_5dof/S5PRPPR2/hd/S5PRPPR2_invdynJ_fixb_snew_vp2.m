% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:10
% EndTime: 2019-12-05 15:24:11
% DurationCPUTime: 0.47s
% Computational Cost: add. (2204->129), mult. (4176->172), div. (0->0), fcn. (2658->10), ass. (0->68)
t293 = cos(pkin(9));
t288 = t293 ^ 2;
t318 = 0.2e1 * t293;
t300 = qJD(2) ^ 2;
t317 = pkin(4) * t293;
t290 = sin(pkin(9));
t316 = t290 * mrSges(5,2);
t315 = t288 * t300;
t292 = sin(pkin(7));
t295 = cos(pkin(7));
t279 = -t295 * g(1) - t292 * g(2);
t289 = -g(3) + qJDD(1);
t297 = sin(qJ(2));
t299 = cos(qJ(2));
t269 = -t297 * t279 + t299 * t289;
t265 = qJDD(2) * pkin(2) + t269;
t270 = t299 * t279 + t297 * t289;
t266 = -t300 * pkin(2) + t270;
t291 = sin(pkin(8));
t294 = cos(pkin(8));
t253 = t291 * t265 + t294 * t266;
t251 = -t300 * pkin(3) + qJDD(2) * qJ(4) + t253;
t278 = -t292 * g(1) + t295 * g(2) + qJDD(3);
t310 = qJD(2) * qJD(4);
t312 = t293 * t278 - 0.2e1 * t290 * t310;
t246 = -t290 * t251 + t312;
t303 = mrSges(5,3) * qJDD(2) + t300 * (-t293 * mrSges(5,1) + t316);
t244 = (-pkin(6) * qJDD(2) + t300 * t317 - t251) * t290 + t312;
t247 = t293 * t251 + t290 * t278 + t310 * t318;
t309 = t293 * qJDD(2);
t245 = -pkin(4) * t315 + pkin(6) * t309 + t247;
t296 = sin(qJ(5));
t298 = cos(qJ(5));
t242 = t298 * t244 - t296 * t245;
t304 = -t290 * t296 + t293 * t298;
t271 = t304 * qJD(2);
t305 = t290 * t298 + t293 * t296;
t272 = t305 * qJD(2);
t258 = -t271 * mrSges(6,1) + t272 * mrSges(6,2);
t261 = t271 * qJD(5) + t305 * qJDD(2);
t267 = -qJD(5) * mrSges(6,2) + t271 * mrSges(6,3);
t239 = m(6) * t242 + qJDD(5) * mrSges(6,1) - t261 * mrSges(6,3) + qJD(5) * t267 - t272 * t258;
t243 = t296 * t244 + t298 * t245;
t260 = -t272 * qJD(5) + t304 * qJDD(2);
t268 = qJD(5) * mrSges(6,1) - t272 * mrSges(6,3);
t240 = m(6) * t243 - qJDD(5) * mrSges(6,2) + t260 * mrSges(6,3) - qJD(5) * t268 + t271 * t258;
t313 = t298 * t239 + t296 * t240;
t230 = m(5) * t246 - t303 * t290 + t313;
t307 = -t296 * t239 + t298 * t240;
t231 = m(5) * t247 + t303 * t293 + t307;
t308 = -t290 * t230 + t293 * t231;
t228 = m(4) * t253 - t300 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t308;
t252 = t294 * t265 - t291 * t266;
t306 = qJDD(4) - t252;
t250 = -qJDD(2) * pkin(3) - t300 * qJ(4) + t306;
t287 = t290 ^ 2;
t248 = (-pkin(3) - t317) * qJDD(2) + (-qJ(4) + (-t287 - t288) * pkin(6)) * t300 + t306;
t302 = m(6) * t248 - t260 * mrSges(6,1) + t261 * mrSges(6,2) - t271 * t267 + t272 * t268;
t301 = -m(5) * t250 + mrSges(5,1) * t309 - t302 + (t287 * t300 + t315) * mrSges(5,3);
t235 = m(4) * t252 - t300 * mrSges(4,2) + (mrSges(4,1) - t316) * qJDD(2) + t301;
t314 = t291 * t228 + t294 * t235;
t256 = Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * qJD(5);
t255 = Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * qJD(5);
t254 = Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * qJD(5);
t241 = qJDD(2) * t316 - t301;
t233 = mrSges(6,2) * t248 - mrSges(6,3) * t242 + Ifges(6,1) * t261 + Ifges(6,4) * t260 + Ifges(6,5) * qJDD(5) - qJD(5) * t255 + t271 * t254;
t232 = -mrSges(6,1) * t248 + mrSges(6,3) * t243 + Ifges(6,4) * t261 + Ifges(6,2) * t260 + Ifges(6,6) * qJDD(5) + qJD(5) * t256 - t272 * t254;
t1 = [m(2) * t289 + t297 * (m(3) * t270 - t300 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t294 * t228 - t291 * t235) + t299 * (m(3) * t269 + qJDD(2) * mrSges(3,1) - t300 * mrSges(3,2) + t314); mrSges(3,1) * t269 - mrSges(3,2) * t270 + mrSges(4,1) * t252 - mrSges(4,2) * t253 + t290 * (mrSges(5,2) * t250 - mrSges(5,3) * t246 - pkin(6) * t313 - t296 * t232 + t298 * t233) + t293 * (-mrSges(5,1) * t250 + mrSges(5,3) * t247 - pkin(4) * t302 + pkin(6) * t307 + t298 * t232 + t296 * t233) - pkin(3) * t241 + qJ(4) * t308 + pkin(2) * t314 + (Ifges(5,2) * t288 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t290 + Ifges(5,4) * t318) * t290) * qJDD(2); m(4) * t278 + t293 * t230 + t290 * t231; t241; mrSges(6,1) * t242 - mrSges(6,2) * t243 + Ifges(6,5) * t261 + Ifges(6,6) * t260 + Ifges(6,3) * qJDD(5) + t272 * t255 - t271 * t256;];
tauJ = t1;
