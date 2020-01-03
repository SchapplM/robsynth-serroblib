% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:07
% EndTime: 2019-12-31 17:38:08
% DurationCPUTime: 0.74s
% Computational Cost: add. (7445->170), mult. (11575->206), div. (0->0), fcn. (5210->8), ass. (0->73)
t345 = m(3) + m(4);
t344 = -pkin(2) - pkin(3);
t343 = mrSges(3,1) + mrSges(4,1);
t342 = Ifges(4,4) + Ifges(3,5);
t341 = Ifges(3,6) - Ifges(4,6);
t340 = cos(pkin(7));
t319 = sin(pkin(7));
t307 = -t340 * g(1) - t319 * g(2);
t316 = -g(3) + qJDD(1);
t322 = sin(qJ(2));
t324 = cos(qJ(2));
t293 = t324 * t307 + t322 * t316;
t325 = qJD(2) ^ 2;
t330 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t293;
t290 = -t325 * pkin(2) + t330;
t287 = t344 * t325 + t330;
t292 = -t322 * t307 + t324 * t316;
t328 = -t325 * qJ(3) + qJDD(3) - t292;
t289 = t344 * qJDD(2) + t328;
t318 = sin(pkin(8));
t320 = cos(pkin(8));
t284 = t320 * t287 + t318 * t289;
t282 = -t325 * pkin(4) - qJDD(2) * pkin(6) + t284;
t306 = t319 * g(1) - t340 * g(2);
t305 = qJDD(4) + t306;
t321 = sin(qJ(5));
t323 = cos(qJ(5));
t279 = -t321 * t282 + t323 * t305;
t302 = (mrSges(6,1) * t323 - mrSges(6,2) * t321) * qJD(2);
t336 = qJD(2) * qJD(5);
t303 = -t321 * qJDD(2) - t323 * t336;
t337 = qJD(2) * t323;
t309 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t337;
t338 = qJD(2) * t321;
t277 = m(6) * t279 + qJDD(5) * mrSges(6,1) - t303 * mrSges(6,3) + qJD(5) * t309 + t302 * t338;
t280 = t323 * t282 + t321 * t305;
t304 = -t323 * qJDD(2) + t321 * t336;
t308 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t338;
t278 = m(6) * t280 - qJDD(5) * mrSges(6,2) + t304 * mrSges(6,3) - qJD(5) * t308 - t302 * t337;
t331 = -t321 * t277 + t323 * t278;
t266 = m(5) * t284 - t325 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + t331;
t283 = -t318 * t287 + t320 * t289;
t281 = qJDD(2) * pkin(4) - t325 * pkin(6) - t283;
t326 = -m(6) * t281 + t304 * mrSges(6,1) - t303 * mrSges(6,2) + t308 * t338 - t309 * t337;
t273 = m(5) * t283 - qJDD(2) * mrSges(5,1) - t325 * mrSges(5,2) + t326;
t332 = t320 * t266 - t318 * t273;
t329 = m(4) * t290 + qJDD(2) * mrSges(4,3) + t332;
t261 = m(3) * t293 - qJDD(2) * mrSges(3,2) - t343 * t325 + t329;
t264 = t318 * t266 + t320 * t273;
t291 = -qJDD(2) * pkin(2) + t328;
t327 = -m(4) * t291 + qJDD(2) * mrSges(4,1) + t325 * mrSges(4,3) - t264;
t262 = m(3) * t292 + qJDD(2) * mrSges(3,1) - t325 * mrSges(3,2) + t327;
t333 = t324 * t261 - t322 * t262;
t255 = m(2) * t307 + t333;
t270 = t323 * t277 + t321 * t278;
t335 = -m(5) * t305 - t270;
t268 = (m(2) + t345) * t306 - t335;
t339 = t319 * t255 + t340 * t268;
t256 = t322 * t261 + t324 * t262;
t334 = t340 * t255 - t319 * t268;
t296 = (Ifges(6,5) * qJD(5)) + (-Ifges(6,1) * t321 - Ifges(6,4) * t323) * qJD(2);
t295 = (Ifges(6,6) * qJD(5)) + (-Ifges(6,4) * t321 - Ifges(6,2) * t323) * qJD(2);
t294 = (Ifges(6,3) * qJD(5)) + (-Ifges(6,5) * t321 - Ifges(6,6) * t323) * qJD(2);
t272 = mrSges(6,2) * t281 - mrSges(6,3) * t279 + Ifges(6,1) * t303 + Ifges(6,4) * t304 + Ifges(6,5) * qJDD(5) - qJD(5) * t295 - t294 * t337;
t271 = -mrSges(6,1) * t281 + mrSges(6,3) * t280 + Ifges(6,4) * t303 + Ifges(6,2) * t304 + Ifges(6,6) * qJDD(5) + qJD(5) * t296 + t294 * t338;
t269 = -m(4) * t306 + t335;
t263 = -mrSges(5,1) * t305 - mrSges(6,1) * t279 + mrSges(6,2) * t280 + mrSges(5,3) * t284 + t325 * Ifges(5,5) - Ifges(6,5) * t303 - Ifges(5,6) * qJDD(2) - Ifges(6,6) * t304 - Ifges(6,3) * qJDD(5) - pkin(4) * t270 + (t295 * t321 - t296 * t323) * qJD(2);
t257 = mrSges(5,2) * t305 - mrSges(5,3) * t283 - Ifges(5,5) * qJDD(2) - t325 * Ifges(5,6) - pkin(6) * t270 - t321 * t271 + t323 * t272;
t252 = mrSges(4,2) * t291 - mrSges(3,3) * t292 - qJ(3) * t269 - qJ(4) * t264 + t320 * t257 - t318 * t263 - t341 * t325 + (-mrSges(3,2) + mrSges(4,3)) * t306 + t342 * qJDD(2);
t251 = mrSges(4,2) * t290 + mrSges(3,3) * t293 - pkin(2) * t269 - pkin(3) * t335 - qJ(4) * t332 + t341 * qJDD(2) - t318 * t257 - t320 * t263 + t343 * t306 + t342 * t325;
t250 = -pkin(1) * t256 - pkin(2) * t327 - qJ(3) * (-t325 * mrSges(4,1) + t329) - mrSges(3,1) * t292 + mrSges(3,2) * t293 + pkin(3) * t264 + mrSges(4,1) * t291 - mrSges(4,3) * t290 + t321 * t272 + t323 * t271 + pkin(4) * t326 + pkin(6) * t331 + mrSges(5,1) * t283 - mrSges(5,2) * t284 + mrSges(2,3) * t307 - mrSges(2,1) * t316 + (-Ifges(3,3) - Ifges(4,2) - Ifges(5,3)) * qJDD(2);
t249 = mrSges(2,2) * t316 - mrSges(2,3) * t306 - pkin(5) * t256 - t322 * t251 + t324 * t252;
t1 = [-m(1) * g(1) + t334; -m(1) * g(2) + t339; -m(1) * g(3) + m(2) * t316 + t256; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t339 + t340 * t249 - t319 * t250; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t334 + t319 * t249 + t340 * t250; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t307 + t322 * t252 + t324 * t251 - pkin(1) * t335 + pkin(5) * t333 + (pkin(1) * t345 + mrSges(2,1)) * t306;];
tauB = t1;
