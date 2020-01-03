% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:58
% EndTime: 2020-01-03 12:01:59
% DurationCPUTime: 0.68s
% Computational Cost: add. (6933->167), mult. (8931->220), div. (0->0), fcn. (5058->10), ass. (0->77)
t312 = qJD(1) + qJD(2);
t318 = sin(qJ(4));
t338 = t312 * t318;
t322 = cos(qJ(4));
t337 = t312 * t322;
t320 = sin(qJ(1));
t324 = cos(qJ(1));
t330 = -t324 * g(2) - t320 * g(3);
t297 = qJDD(1) * pkin(1) + t330;
t333 = -t320 * g(2) + t324 * g(3);
t298 = -qJD(1) ^ 2 * pkin(1) + t333;
t319 = sin(qJ(2));
t323 = cos(qJ(2));
t280 = t323 * t297 - t319 * t298;
t310 = qJDD(1) + qJDD(2);
t277 = t310 * pkin(2) + t280;
t281 = t319 * t297 + t323 * t298;
t308 = t312 ^ 2;
t278 = -t308 * pkin(2) + t281;
t315 = sin(pkin(9));
t316 = cos(pkin(9));
t262 = t315 * t277 + t316 * t278;
t259 = -t308 * pkin(3) + t310 * pkin(7) + t262;
t314 = -g(1) + qJDD(3);
t255 = -t318 * t259 + t322 * t314;
t335 = qJD(4) * t312;
t334 = t322 * t335;
t292 = t318 * t310 + t334;
t252 = (-t292 + t334) * pkin(8) + (t308 * t318 * t322 + qJDD(4)) * pkin(4) + t255;
t256 = t322 * t259 + t318 * t314;
t293 = t322 * t310 - t318 * t335;
t301 = qJD(4) * pkin(4) - pkin(8) * t338;
t313 = t322 ^ 2;
t253 = -t313 * t308 * pkin(4) + t293 * pkin(8) - qJD(4) * t301 + t256;
t317 = sin(qJ(5));
t321 = cos(qJ(5));
t250 = t321 * t252 - t317 * t253;
t287 = (-t317 * t318 + t321 * t322) * t312;
t268 = t287 * qJD(5) + t321 * t292 + t317 * t293;
t288 = (t317 * t322 + t318 * t321) * t312;
t273 = -t287 * mrSges(6,1) + t288 * mrSges(6,2);
t311 = qJD(4) + qJD(5);
t282 = -t311 * mrSges(6,2) + t287 * mrSges(6,3);
t309 = qJDD(4) + qJDD(5);
t247 = m(6) * t250 + t309 * mrSges(6,1) - t268 * mrSges(6,3) - t288 * t273 + t311 * t282;
t251 = t317 * t252 + t321 * t253;
t267 = -t288 * qJD(5) - t317 * t292 + t321 * t293;
t283 = t311 * mrSges(6,1) - t288 * mrSges(6,3);
t248 = m(6) * t251 - t309 * mrSges(6,2) + t267 * mrSges(6,3) + t287 * t273 - t311 * t283;
t238 = t321 * t247 + t317 * t248;
t291 = (-mrSges(5,1) * t322 + mrSges(5,2) * t318) * t312;
t300 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t337;
t236 = m(5) * t255 + qJDD(4) * mrSges(5,1) - t292 * mrSges(5,3) + qJD(4) * t300 - t291 * t338 + t238;
t299 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t338;
t331 = -t317 * t247 + t321 * t248;
t237 = m(5) * t256 - qJDD(4) * mrSges(5,2) + t293 * mrSges(5,3) - qJD(4) * t299 + t291 * t337 + t331;
t332 = -t318 * t236 + t322 * t237;
t233 = m(4) * t262 - t308 * mrSges(4,1) - t310 * mrSges(4,2) + t332;
t261 = t316 * t277 - t315 * t278;
t329 = -t310 * pkin(3) - t261;
t258 = -t308 * pkin(7) + t329;
t254 = t301 * t338 - t293 * pkin(4) + (-pkin(8) * t313 - pkin(7)) * t308 + t329;
t328 = m(6) * t254 - t267 * mrSges(6,1) + t268 * mrSges(6,2) - t287 * t282 + t288 * t283;
t325 = -m(5) * t258 + t293 * mrSges(5,1) - t292 * mrSges(5,2) - t299 * t338 + t300 * t337 - t328;
t242 = m(4) * t261 + t310 * mrSges(4,1) - t308 * mrSges(4,2) + t325;
t336 = t315 * t233 + t316 * t242;
t270 = Ifges(6,4) * t288 + Ifges(6,2) * t287 + Ifges(6,6) * t311;
t271 = Ifges(6,1) * t288 + Ifges(6,4) * t287 + Ifges(6,5) * t311;
t327 = mrSges(6,1) * t250 - mrSges(6,2) * t251 + Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * t309 + t288 * t270 - t287 * t271;
t269 = Ifges(6,5) * t288 + Ifges(6,6) * t287 + Ifges(6,3) * t311;
t239 = -mrSges(6,1) * t254 + mrSges(6,3) * t251 + Ifges(6,4) * t268 + Ifges(6,2) * t267 + Ifges(6,6) * t309 - t288 * t269 + t311 * t271;
t240 = mrSges(6,2) * t254 - mrSges(6,3) * t250 + Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * t309 + t287 * t269 - t311 * t270;
t284 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t318 + Ifges(5,6) * t322) * t312;
t285 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t318 + Ifges(5,2) * t322) * t312;
t286 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t318 + Ifges(5,4) * t322) * t312;
t326 = -mrSges(3,2) * t281 - mrSges(4,2) * t262 + t322 * (-mrSges(5,1) * t258 + mrSges(5,3) * t256 + Ifges(5,4) * t292 + Ifges(5,2) * t293 + Ifges(5,6) * qJDD(4) - pkin(4) * t328 + pkin(8) * t331 + qJD(4) * t286 + t321 * t239 + t317 * t240 - t284 * t338) + pkin(2) * t336 + t318 * (mrSges(5,2) * t258 - mrSges(5,3) * t255 + Ifges(5,1) * t292 + Ifges(5,4) * t293 + Ifges(5,5) * qJDD(4) - pkin(8) * t238 - qJD(4) * t285 - t317 * t239 + t321 * t240 + t284 * t337) + pkin(7) * t332 + pkin(3) * t325 + mrSges(4,1) * t261 + mrSges(3,1) * t280 + (Ifges(4,3) + Ifges(3,3)) * t310;
t1 = [Ifges(2,3) * qJDD(1) - mrSges(2,2) * t333 + pkin(1) * (t319 * (m(3) * t281 - t308 * mrSges(3,1) - t310 * mrSges(3,2) + t316 * t233 - t315 * t242) + t323 * (m(3) * t280 + t310 * mrSges(3,1) - t308 * mrSges(3,2) + t336)) + mrSges(2,1) * t330 + t326; t326; m(4) * t314 + t322 * t236 + t318 * t237; mrSges(5,1) * t255 - mrSges(5,2) * t256 + Ifges(5,5) * t292 + Ifges(5,6) * t293 + Ifges(5,3) * qJDD(4) + pkin(4) * t238 + (t285 * t318 - t286 * t322) * t312 + t327; t327;];
tauJ = t1;
