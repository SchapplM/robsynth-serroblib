% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:49
% EndTime: 2019-12-31 16:27:50
% DurationCPUTime: 0.77s
% Computational Cost: add. (3762->161), mult. (7139->203), div. (0->0), fcn. (3528->6), ass. (0->67)
t349 = Ifges(4,1) + Ifges(5,1);
t344 = Ifges(4,4) - Ifges(5,5);
t343 = Ifges(5,4) + Ifges(4,5);
t348 = Ifges(4,2) + Ifges(5,3);
t347 = Ifges(5,6) - Ifges(4,6);
t346 = Ifges(4,3) + Ifges(5,2);
t345 = mrSges(4,3) + mrSges(5,2);
t318 = -g(3) + qJDD(1);
t323 = cos(qJ(3));
t341 = t323 * t318;
t319 = sin(pkin(6));
t320 = cos(pkin(6));
t307 = t319 * g(1) - t320 * g(2);
t308 = -t320 * g(1) - t319 * g(2);
t322 = sin(qJ(2));
t324 = cos(qJ(2));
t284 = t322 * t307 + t324 * t308;
t326 = qJD(2) ^ 2;
t282 = -t326 * pkin(2) + qJDD(2) * pkin(5) + t284;
t321 = sin(qJ(3));
t279 = t323 * t282 + t321 * t318;
t302 = (-mrSges(4,1) * t323 + mrSges(4,2) * t321) * qJD(2);
t334 = qJD(2) * qJD(3);
t304 = t323 * qJDD(2) - t321 * t334;
t336 = qJD(2) * t321;
t309 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t336;
t300 = (-pkin(3) * t323 - qJ(4) * t321) * qJD(2);
t325 = qJD(3) ^ 2;
t335 = qJD(2) * t323;
t276 = -t325 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t300 * t335 + t279;
t301 = (-mrSges(5,1) * t323 - mrSges(5,3) * t321) * qJD(2);
t310 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t336;
t329 = m(5) * t276 + qJDD(3) * mrSges(5,3) + qJD(3) * t310 + t301 * t335;
t271 = m(4) * t279 - qJDD(3) * mrSges(4,2) - qJD(3) * t309 + t302 * t335 + t345 * t304 + t329;
t278 = -t321 * t282 + t341;
t303 = t321 * qJDD(2) + t323 * t334;
t311 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t335;
t277 = -qJDD(3) * pkin(3) - t325 * qJ(4) - t341 + qJDD(4) + (qJD(2) * t300 + t282) * t321;
t312 = mrSges(5,2) * t335 + qJD(3) * mrSges(5,3);
t328 = -m(5) * t277 + qJDD(3) * mrSges(5,1) + qJD(3) * t312;
t272 = m(4) * t278 + qJDD(3) * mrSges(4,1) + qJD(3) * t311 - t345 * t303 + (-t301 - t302) * t336 + t328;
t330 = t323 * t271 - t321 * t272;
t264 = m(3) * t284 - t326 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t330;
t283 = t324 * t307 - t322 * t308;
t281 = -qJDD(2) * pkin(2) - t326 * pkin(5) - t283;
t274 = -t304 * pkin(3) - t303 * qJ(4) + (-0.2e1 * qJD(4) * t321 + (pkin(3) * t321 - qJ(4) * t323) * qJD(3)) * qJD(2) + t281;
t273 = m(5) * t274 - t304 * mrSges(5,1) - t303 * mrSges(5,3) - t310 * t336 - t312 * t335;
t327 = -m(4) * t281 + t304 * mrSges(4,1) - t303 * mrSges(4,2) - t309 * t336 + t311 * t335 - t273;
t267 = m(3) * t283 + qJDD(2) * mrSges(3,1) - t326 * mrSges(3,2) + t327;
t259 = t322 * t264 + t324 * t267;
t257 = m(2) * t307 + t259;
t331 = t324 * t264 - t322 * t267;
t258 = m(2) * t308 + t331;
t340 = t320 * t257 + t319 * t258;
t265 = t321 * t271 + t323 * t272;
t339 = t347 * qJD(3) + (-t344 * t321 - t348 * t323) * qJD(2);
t338 = t346 * qJD(3) + (t343 * t321 - t347 * t323) * qJD(2);
t337 = t343 * qJD(3) + (t349 * t321 + t344 * t323) * qJD(2);
t333 = m(3) * t318 + t265;
t332 = -t319 * t257 + t320 * t258;
t261 = mrSges(4,2) * t281 + mrSges(5,2) * t277 - mrSges(4,3) * t278 - mrSges(5,3) * t274 - qJ(4) * t273 + t339 * qJD(3) + t343 * qJDD(3) + t349 * t303 + t344 * t304 + t338 * t335;
t260 = -mrSges(4,1) * t281 - mrSges(5,1) * t274 + mrSges(5,2) * t276 + mrSges(4,3) * t279 - pkin(3) * t273 + t337 * qJD(3) - qJDD(3) * t347 + t344 * t303 + t348 * t304 - t338 * t336;
t253 = Ifges(3,6) * qJDD(2) + t326 * Ifges(3,5) - mrSges(3,1) * t318 + mrSges(3,3) * t284 - mrSges(4,1) * t278 + mrSges(4,2) * t279 + mrSges(5,1) * t277 - mrSges(5,3) * t276 - pkin(3) * t328 - qJ(4) * t329 - pkin(2) * t265 + (-qJ(4) * mrSges(5,2) + t347) * t304 + (pkin(3) * mrSges(5,2) - t343) * t303 - t346 * qJDD(3) + (t337 * t323 + (pkin(3) * t301 + t339) * t321) * qJD(2);
t252 = mrSges(3,2) * t318 - mrSges(3,3) * t283 + Ifges(3,5) * qJDD(2) - t326 * Ifges(3,6) - pkin(5) * t265 - t321 * t260 + t323 * t261;
t251 = mrSges(2,2) * t318 - mrSges(2,3) * t307 - pkin(4) * t259 + t324 * t252 - t322 * t253;
t250 = -mrSges(2,1) * t318 + mrSges(2,3) * t308 - pkin(1) * t333 + pkin(4) * t331 + t322 * t252 + t324 * t253;
t1 = [-m(1) * g(1) + t332; -m(1) * g(2) + t340; -m(1) * g(3) + m(2) * t318 + t333; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t340 - t319 * t250 + t320 * t251; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t332 + t320 * t250 + t319 * t251; -mrSges(1,1) * g(2) + mrSges(2,1) * t307 + mrSges(3,1) * t283 + mrSges(1,2) * g(1) - mrSges(2,2) * t308 - mrSges(3,2) * t284 + Ifges(3,3) * qJDD(2) + pkin(1) * t259 + pkin(2) * t327 + pkin(5) * t330 + t323 * t260 + t321 * t261;];
tauB = t1;
