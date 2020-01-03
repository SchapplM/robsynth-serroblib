% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:16
% EndTime: 2019-12-31 16:30:17
% DurationCPUTime: 0.77s
% Computational Cost: add. (3409->162), mult. (6323->203), div. (0->0), fcn. (3048->6), ass. (0->66)
t348 = Ifges(4,1) + Ifges(5,1);
t343 = Ifges(4,4) - Ifges(5,5);
t342 = Ifges(5,4) + Ifges(4,5);
t347 = Ifges(4,2) + Ifges(5,3);
t346 = Ifges(5,6) - Ifges(4,6);
t345 = Ifges(4,3) + Ifges(5,2);
t344 = mrSges(4,3) + mrSges(5,2);
t340 = cos(pkin(6));
t319 = sin(pkin(6));
t307 = t319 * g(1) - t340 * g(2);
t322 = cos(qJ(3));
t339 = t322 * t307;
t308 = -t340 * g(1) - t319 * g(2);
t318 = -g(3) + qJDD(1);
t321 = sin(qJ(2));
t323 = cos(qJ(2));
t284 = t323 * t308 + t321 * t318;
t325 = qJD(2) ^ 2;
t282 = -t325 * pkin(2) + qJDD(2) * pkin(5) + t284;
t320 = sin(qJ(3));
t279 = t322 * t282 - t320 * t307;
t302 = (-mrSges(4,1) * t322 + mrSges(4,2) * t320) * qJD(2);
t332 = qJD(2) * qJD(3);
t304 = t322 * qJDD(2) - t320 * t332;
t334 = qJD(2) * t320;
t309 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t334;
t300 = (-pkin(3) * t322 - qJ(4) * t320) * qJD(2);
t324 = qJD(3) ^ 2;
t333 = qJD(2) * t322;
t276 = -t324 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t300 * t333 + t279;
t301 = (-mrSges(5,1) * t322 - mrSges(5,3) * t320) * qJD(2);
t310 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t334;
t328 = m(5) * t276 + qJDD(3) * mrSges(5,3) + qJD(3) * t310 + t301 * t333;
t271 = m(4) * t279 - qJDD(3) * mrSges(4,2) - qJD(3) * t309 + t302 * t333 + t344 * t304 + t328;
t278 = -t320 * t282 - t339;
t303 = t320 * qJDD(2) + t322 * t332;
t311 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t333;
t277 = -qJDD(3) * pkin(3) - t324 * qJ(4) + t339 + qJDD(4) + (qJD(2) * t300 + t282) * t320;
t312 = mrSges(5,2) * t333 + qJD(3) * mrSges(5,3);
t327 = -m(5) * t277 + qJDD(3) * mrSges(5,1) + qJD(3) * t312;
t272 = m(4) * t278 + qJDD(3) * mrSges(4,1) + qJD(3) * t311 - t344 * t303 + (-t301 - t302) * t334 + t327;
t329 = t322 * t271 - t320 * t272;
t264 = m(3) * t284 - t325 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t329;
t283 = -t321 * t308 + t323 * t318;
t281 = -qJDD(2) * pkin(2) - t325 * pkin(5) - t283;
t274 = -t304 * pkin(3) - t303 * qJ(4) + (-0.2e1 * qJD(4) * t320 + (pkin(3) * t320 - qJ(4) * t322) * qJD(3)) * qJD(2) + t281;
t273 = m(5) * t274 - t304 * mrSges(5,1) - t303 * mrSges(5,3) - t310 * t334 - t312 * t333;
t326 = -m(4) * t281 + t304 * mrSges(4,1) - t303 * mrSges(4,2) - t309 * t334 + t311 * t333 - t273;
t269 = m(3) * t283 + qJDD(2) * mrSges(3,1) - t325 * mrSges(3,2) + t326;
t330 = t323 * t264 - t321 * t269;
t258 = m(2) * t308 + t330;
t267 = t320 * t271 + t322 * t272;
t266 = (m(2) + m(3)) * t307 - t267;
t338 = t319 * t258 + t340 * t266;
t259 = t321 * t264 + t323 * t269;
t337 = t346 * qJD(3) + (-t343 * t320 - t347 * t322) * qJD(2);
t336 = t345 * qJD(3) + (t342 * t320 - t346 * t322) * qJD(2);
t335 = t342 * qJD(3) + (t348 * t320 + t343 * t322) * qJD(2);
t331 = t340 * t258 - t319 * t266;
t261 = mrSges(4,2) * t281 + mrSges(5,2) * t277 - mrSges(4,3) * t278 - mrSges(5,3) * t274 - qJ(4) * t273 + t337 * qJD(3) + t342 * qJDD(3) + t348 * t303 + t343 * t304 + t336 * t333;
t260 = -mrSges(4,1) * t281 - mrSges(5,1) * t274 + mrSges(5,2) * t276 + mrSges(4,3) * t279 - pkin(3) * t273 + t335 * qJD(3) - qJDD(3) * t346 + t343 * t303 + t347 * t304 - t336 * t334;
t255 = Ifges(3,6) * qJDD(2) + t325 * Ifges(3,5) + mrSges(3,1) * t307 + mrSges(3,3) * t284 - mrSges(4,1) * t278 + mrSges(4,2) * t279 + mrSges(5,1) * t277 - mrSges(5,3) * t276 - pkin(3) * t327 - qJ(4) * t328 - pkin(2) * t267 + (-qJ(4) * mrSges(5,2) + t346) * t304 + (pkin(3) * mrSges(5,2) - t342) * t303 - t345 * qJDD(3) + (t335 * t322 + (pkin(3) * t301 + t337) * t320) * qJD(2);
t254 = -mrSges(3,2) * t307 - mrSges(3,3) * t283 + Ifges(3,5) * qJDD(2) - t325 * Ifges(3,6) - pkin(5) * t267 - t320 * t260 + t322 * t261;
t253 = -mrSges(2,1) * t318 - mrSges(3,1) * t283 + mrSges(3,2) * t284 + mrSges(2,3) * t308 - Ifges(3,3) * qJDD(2) - pkin(1) * t259 - pkin(2) * t326 - pkin(5) * t329 - t322 * t260 - t320 * t261;
t252 = mrSges(2,2) * t318 - mrSges(2,3) * t307 - pkin(4) * t259 + t323 * t254 - t321 * t255;
t1 = [-m(1) * g(1) + t331; -m(1) * g(2) + t338; -m(1) * g(3) + m(2) * t318 + t259; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t338 + t340 * t252 - t319 * t253; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t331 + t319 * t252 + t340 * t253; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t307 - mrSges(2,2) * t308 + t321 * t254 + t323 * t255 + pkin(1) * (m(3) * t307 - t267) + pkin(4) * t330;];
tauB = t1;
