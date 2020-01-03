% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:19
% EndTime: 2019-12-31 18:42:20
% DurationCPUTime: 0.99s
% Computational Cost: add. (3855->196), mult. (7294->237), div. (0->0), fcn. (4190->8), ass. (0->84)
t372 = Ifges(5,1) + Ifges(6,1);
t364 = Ifges(5,4) + Ifges(6,4);
t363 = Ifges(5,5) + Ifges(6,5);
t371 = Ifges(5,2) + Ifges(6,2);
t362 = Ifges(5,6) + Ifges(6,6);
t335 = sin(qJ(4));
t338 = cos(qJ(4));
t336 = sin(qJ(3));
t357 = t336 * qJD(1);
t319 = t338 * qJD(3) - t335 * t357;
t339 = cos(qJ(3));
t355 = qJD(1) * qJD(3);
t351 = t339 * t355;
t325 = t336 * qJDD(1) + t351;
t296 = t319 * qJD(4) + t335 * qJDD(3) + t338 * t325;
t320 = t335 * qJD(3) + t338 * t357;
t300 = -t319 * mrSges(6,1) + t320 * mrSges(6,2);
t337 = sin(qJ(1));
t340 = cos(qJ(1));
t350 = t337 * g(1) - t340 * g(2);
t321 = qJDD(1) * pkin(1) + t350;
t342 = qJD(1) ^ 2;
t346 = -t340 * g(1) - t337 * g(2);
t323 = -t342 * pkin(1) + t346;
t333 = sin(pkin(8));
t334 = cos(pkin(8));
t297 = t334 * t321 - t333 * t323;
t282 = -qJDD(1) * pkin(2) - t342 * pkin(6) - t297;
t352 = t336 * t355;
t326 = t339 * qJDD(1) - t352;
t273 = (-t325 - t351) * pkin(7) + (-t326 + t352) * pkin(3) + t282;
t298 = t333 * t321 + t334 * t323;
t283 = -t342 * pkin(2) + qJDD(1) * pkin(6) + t298;
t332 = -g(3) + qJDD(2);
t278 = t339 * t283 + t336 * t332;
t324 = (-t339 * pkin(3) - t336 * pkin(7)) * qJD(1);
t341 = qJD(3) ^ 2;
t356 = t339 * qJD(1);
t276 = -t341 * pkin(3) + qJDD(3) * pkin(7) + t324 * t356 + t278;
t269 = t338 * t273 - t335 * t276;
t318 = qJDD(4) - t326;
t329 = qJD(4) - t356;
t265 = -0.2e1 * qJD(5) * t320 + (t319 * t329 - t296) * qJ(5) + (t319 * t320 + t318) * pkin(4) + t269;
t302 = -t329 * mrSges(6,2) + t319 * mrSges(6,3);
t354 = m(6) * t265 + t318 * mrSges(6,1) + t329 * t302;
t262 = -t296 * mrSges(6,3) - t320 * t300 + t354;
t270 = t335 * t273 + t338 * t276;
t295 = -t320 * qJD(4) + t338 * qJDD(3) - t335 * t325;
t304 = t329 * pkin(4) - t320 * qJ(5);
t317 = t319 ^ 2;
t267 = -t317 * pkin(4) + t295 * qJ(5) + 0.2e1 * qJD(5) * t319 - t329 * t304 + t270;
t360 = -t371 * t319 - t364 * t320 - t362 * t329;
t366 = t364 * t319 + t372 * t320 + t363 * t329;
t368 = Ifges(5,3) + Ifges(6,3);
t370 = mrSges(5,1) * t269 + mrSges(6,1) * t265 - mrSges(5,2) * t270 - mrSges(6,2) * t267 + pkin(4) * t262 + t362 * t295 + t363 * t296 + t368 * t318 - t366 * t319 - t360 * t320;
t365 = -mrSges(5,2) - mrSges(6,2);
t361 = -t362 * t319 - t363 * t320 - t368 * t329;
t305 = t329 * mrSges(6,1) - t320 * mrSges(6,3);
t358 = -t329 * mrSges(5,1) + t320 * mrSges(5,3) - t305;
t353 = m(6) * t267 + t295 * mrSges(6,3) + t319 * t300;
t322 = (-t339 * mrSges(4,1) + t336 * mrSges(4,2)) * qJD(1);
t327 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t357;
t301 = -t319 * mrSges(5,1) + t320 * mrSges(5,2);
t303 = -t329 * mrSges(5,2) + t319 * mrSges(5,3);
t258 = m(5) * t269 + t318 * mrSges(5,1) + t329 * t303 + (-t300 - t301) * t320 + (-mrSges(5,3) - mrSges(6,3)) * t296 + t354;
t260 = m(5) * t270 + t295 * mrSges(5,3) + t319 * t301 + t365 * t318 + t358 * t329 + t353;
t348 = -t335 * t258 + t338 * t260;
t255 = m(4) * t278 - qJDD(3) * mrSges(4,2) + t326 * mrSges(4,3) - qJD(3) * t327 + t322 * t356 + t348;
t277 = -t336 * t283 + t339 * t332;
t328 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t356;
t275 = -qJDD(3) * pkin(3) - t341 * pkin(7) + t324 * t357 - t277;
t268 = -t295 * pkin(4) - t317 * qJ(5) + t320 * t304 + qJDD(5) + t275;
t347 = -m(6) * t268 + t295 * mrSges(6,1) + t319 * t302;
t343 = -m(5) * t275 + t295 * mrSges(5,1) + t365 * t296 + t319 * t303 + t358 * t320 + t347;
t261 = m(4) * t277 + qJDD(3) * mrSges(4,1) - t325 * mrSges(4,3) + qJD(3) * t328 - t322 * t357 + t343;
t349 = t339 * t255 - t336 * t261;
t257 = t338 * t258 + t335 * t260;
t344 = -m(4) * t282 + t326 * mrSges(4,1) - t325 * mrSges(4,2) - t327 * t357 + t328 * t356 - t257;
t313 = Ifges(4,5) * qJD(3) + (t336 * Ifges(4,1) + t339 * Ifges(4,4)) * qJD(1);
t312 = Ifges(4,6) * qJD(3) + (t336 * Ifges(4,4) + t339 * Ifges(4,2)) * qJD(1);
t263 = t296 * mrSges(6,2) + t320 * t305 - t347;
t256 = mrSges(5,2) * t275 + mrSges(6,2) * t268 - mrSges(5,3) * t269 - mrSges(6,3) * t265 - qJ(5) * t262 + t364 * t295 + t372 * t296 + t363 * t318 - t361 * t319 + t360 * t329;
t253 = -mrSges(5,1) * t275 + mrSges(5,3) * t270 - mrSges(6,1) * t268 + mrSges(6,3) * t267 - pkin(4) * t263 + qJ(5) * t353 + (-qJ(5) * t305 + t366) * t329 + t361 * t320 + (-qJ(5) * mrSges(6,2) + t362) * t318 + t364 * t296 + t371 * t295;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t350 - mrSges(2,2) * t346 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t297 - mrSges(3,2) * t298 + t336 * (mrSges(4,2) * t282 - mrSges(4,3) * t277 + Ifges(4,1) * t325 + Ifges(4,4) * t326 + Ifges(4,5) * qJDD(3) - pkin(7) * t257 - qJD(3) * t312 - t335 * t253 + t338 * t256) + t339 * (-mrSges(4,1) * t282 + mrSges(4,3) * t278 + Ifges(4,4) * t325 + Ifges(4,2) * t326 + Ifges(4,6) * qJDD(3) - pkin(3) * t257 + qJD(3) * t313 - t370) + pkin(2) * t344 + pkin(6) * t349 + pkin(1) * (t333 * (m(3) * t298 - t342 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t349) + t334 * (m(3) * t297 + qJDD(1) * mrSges(3,1) - t342 * mrSges(3,2) + t344)); m(3) * t332 + t336 * t255 + t339 * t261; Ifges(4,5) * t325 + Ifges(4,6) * t326 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t277 - mrSges(4,2) * t278 + t335 * t256 + t338 * t253 + pkin(3) * t343 + pkin(7) * t348 + (t336 * t312 - t339 * t313) * qJD(1); t370; t263;];
tauJ = t1;
