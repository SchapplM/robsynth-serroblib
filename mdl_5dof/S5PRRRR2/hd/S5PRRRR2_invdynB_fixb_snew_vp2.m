% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:16
% EndTime: 2019-07-18 13:30:16
% DurationCPUTime: 0.72s
% Computational Cost: add. (10642->161), mult. (11950->202), div. (0->0), fcn. (6252->8), ass. (0->71)
t338 = -m(1) - m(2);
t314 = qJD(2) + qJD(3);
t307 = qJD(4) + t314;
t316 = sin(qJ(5));
t337 = t307 * t316;
t320 = cos(qJ(5));
t336 = t307 * t320;
t319 = sin(qJ(2));
t323 = cos(qJ(2));
t303 = t319 * g(1) - t323 * g(2);
t301 = qJDD(2) * pkin(2) + t303;
t304 = -t323 * g(1) - t319 * g(2);
t324 = qJD(2) ^ 2;
t302 = -t324 * pkin(2) + t304;
t318 = sin(qJ(3));
t322 = cos(qJ(3));
t288 = t322 * t301 - t318 * t302;
t313 = qJDD(2) + qJDD(3);
t286 = t313 * pkin(3) + t288;
t289 = t318 * t301 + t322 * t302;
t312 = t314 ^ 2;
t287 = -t312 * pkin(3) + t289;
t317 = sin(qJ(4));
t321 = cos(qJ(4));
t283 = t317 * t286 + t321 * t287;
t305 = t307 ^ 2;
t306 = qJDD(4) + t313;
t280 = t306 * pkin(6) + t283;
t315 = -g(3) + qJDD(1);
t278 = -t316 * t280 + t320 * t315;
t293 = (-mrSges(6,1) * t320 + mrSges(6,2) * t316) * t307;
t333 = qJD(5) * t307;
t294 = t316 * t306 + t320 * t333;
t300 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t336;
t276 = m(6) * t278 + qJDD(5) * mrSges(6,1) - t294 * mrSges(6,3) + qJD(5) * t300 - t293 * t337;
t279 = t320 * t280 + t316 * t315;
t295 = t320 * t306 - t316 * t333;
t299 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t337;
t277 = m(6) * t279 - qJDD(5) * mrSges(6,2) + t295 * mrSges(6,3) - qJD(5) * t299 + t293 * t336;
t328 = -t316 * t276 + t320 * t277;
t267 = m(5) * t283 - t305 * mrSges(5,1) - t306 * mrSges(5,2) + t328;
t282 = t321 * t286 - t317 * t287;
t281 = -t305 * pkin(6) - t282;
t272 = m(5) * t282 - m(6) * t281 + t306 * mrSges(5,1) + t295 * mrSges(6,1) - t305 * mrSges(5,2) - t294 * mrSges(6,2) + (-t299 * t316 + t300 * t320) * t307;
t264 = t317 * t267 + t321 * t272;
t262 = m(4) * t288 + t313 * mrSges(4,1) - t312 * mrSges(4,2) + t264;
t329 = t321 * t267 - t317 * t272;
t263 = m(4) * t289 - t312 * mrSges(4,1) - t313 * mrSges(4,2) + t329;
t257 = t322 * t262 + t318 * t263;
t255 = m(3) * t303 + qJDD(2) * mrSges(3,1) - t324 * mrSges(3,2) + t257;
t330 = -t318 * t262 + t322 * t263;
t256 = m(3) * t304 - t324 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t330;
t335 = t323 * t255 + t319 * t256;
t334 = t320 * t276 + t316 * t277;
t332 = m(5) * t315 + t334;
t331 = -t319 * t255 + t323 * t256;
t327 = m(4) * t315 + t332;
t326 = qJ(1) * m(2) + mrSges(1,3) + mrSges(2,3);
t325 = m(3) * t315 + t327;
t292 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t316 + Ifges(6,4) * t320) * t307;
t291 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t316 + Ifges(6,2) * t320) * t307;
t290 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t316 + Ifges(6,6) * t320) * t307;
t270 = mrSges(6,2) * t281 - mrSges(6,3) * t278 + Ifges(6,1) * t294 + Ifges(6,4) * t295 + Ifges(6,5) * qJDD(5) - qJD(5) * t291 + t290 * t336;
t269 = -mrSges(6,1) * t281 + mrSges(6,3) * t279 + Ifges(6,4) * t294 + Ifges(6,2) * t295 + Ifges(6,6) * qJDD(5) + qJD(5) * t292 - t290 * t337;
t268 = -mrSges(5,1) * t315 - mrSges(6,1) * t278 + mrSges(6,2) * t279 + mrSges(5,3) * t283 + t305 * Ifges(5,5) - Ifges(6,5) * t294 + Ifges(5,6) * t306 - Ifges(6,6) * t295 - Ifges(6,3) * qJDD(5) + (-t291 * t316 + t292 * t320) * t307;
t258 = mrSges(5,2) * t315 - mrSges(5,3) * t282 + Ifges(5,5) * t306 - t305 * Ifges(5,6) - pkin(6) * t334 - t316 * t269 + t320 * t270;
t251 = mrSges(4,2) * t315 - mrSges(4,3) * t288 + Ifges(4,5) * t313 - t312 * Ifges(4,6) - pkin(5) * t264 + t321 * t258 - t317 * t268;
t250 = -mrSges(4,1) * t315 + mrSges(4,3) * t289 + t312 * Ifges(4,5) + Ifges(4,6) * t313 - pkin(3) * t332 + pkin(5) * t329 + t317 * t258 + t321 * t268;
t249 = mrSges(3,2) * t315 - mrSges(3,3) * t303 + Ifges(3,5) * qJDD(2) - t324 * Ifges(3,6) - pkin(4) * t257 - t318 * t250 + t322 * t251;
t248 = -mrSges(3,1) * t315 + mrSges(3,3) * t304 + t324 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t327 + pkin(4) * t330 + t322 * t250 + t318 * t251;
t1 = [t338 * g(1) + t331; t338 * g(2) + t335; -m(1) * g(3) + m(2) * t315 + t325; -mrSges(1,2) * g(3) + mrSges(2,2) * t315 + t326 * g(2) - qJ(1) * t335 - t319 * t248 + t323 * t249; mrSges(1,1) * g(3) - mrSges(2,1) * t315 - pkin(1) * t325 - t326 * g(1) + qJ(1) * t331 + t323 * t248 + t319 * t249; pkin(1) * t335 + pkin(2) * t257 + mrSges(3,1) * t303 - mrSges(3,2) * t304 + pkin(3) * t264 + mrSges(4,1) * t288 - mrSges(4,2) * t289 + t320 * t269 + pkin(6) * t328 + mrSges(5,1) * t282 - mrSges(5,2) * t283 + t316 * t270 + Ifges(5,3) * t306 + Ifges(4,3) * t313 + Ifges(3,3) * qJDD(2) + (-mrSges(1,1) - mrSges(2,1)) * g(2) + (mrSges(2,2) + mrSges(1,2)) * g(1);];
tauB  = t1;
