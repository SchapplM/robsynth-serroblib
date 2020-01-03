% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:51
% EndTime: 2019-12-31 17:36:51
% DurationCPUTime: 0.48s
% Computational Cost: add. (1318->134), mult. (2875->172), div. (0->0), fcn. (1740->8), ass. (0->69)
t319 = cos(pkin(8));
t352 = 0.2e1 * t319;
t325 = qJD(2) ^ 2;
t318 = sin(pkin(7));
t320 = cos(pkin(7));
t309 = t318 * g(1) - t320 * g(2);
t310 = -t320 * g(1) - t318 * g(2);
t322 = sin(qJ(2));
t324 = cos(qJ(2));
t343 = t322 * t309 + t324 * t310;
t293 = -t325 * pkin(2) + qJDD(2) * qJ(3) + t343;
t316 = -g(3) + qJDD(1);
t317 = sin(pkin(8));
t339 = (qJD(2) * qJD(3));
t280 = t319 * t316 + (-t293 - (2 * t339)) * t317;
t332 = -t319 * mrSges(5,1) - t317 * mrSges(5,3);
t305 = t332 * qJD(2);
t351 = (t305 + (-t319 * mrSges(4,1) + t317 * mrSges(4,2)) * qJD(2)) * qJD(2) + (mrSges(5,2) + mrSges(4,3)) * qJDD(2);
t315 = t319 ^ 2;
t350 = pkin(4) * t325;
t348 = pkin(6) * qJDD(2);
t347 = t317 * qJ(4);
t346 = t325 * qJ(3);
t331 = -t319 * pkin(3) - t347;
t304 = t331 * qJD(2);
t340 = t317 * qJD(2);
t277 = t304 * t340 + qJDD(4) - t280;
t274 = (-t319 * t350 - t348) * t317 + t277;
t281 = t319 * t293 + t317 * t316 + t339 * t352;
t278 = t319 * qJD(2) * t304 + t281;
t275 = -t315 * t350 - t319 * t348 + t278;
t321 = sin(qJ(5));
t323 = cos(qJ(5));
t272 = t323 * t274 - t321 * t275;
t329 = -t317 * t321 - t319 * t323;
t298 = t329 * qJD(2);
t330 = t317 * t323 - t319 * t321;
t299 = t330 * qJD(2);
t288 = -t298 * mrSges(6,1) + t299 * mrSges(6,2);
t292 = t298 * qJD(5) + t330 * qJDD(2);
t294 = -qJD(5) * mrSges(6,2) + t298 * mrSges(6,3);
t270 = m(6) * t272 + qJDD(5) * mrSges(6,1) - t292 * mrSges(6,3) + qJD(5) * t294 - t299 * t288;
t273 = t321 * t274 + t323 * t275;
t291 = -t299 * qJD(5) + t329 * qJDD(2);
t295 = qJD(5) * mrSges(6,1) - t299 * mrSges(6,3);
t271 = m(6) * t273 - qJDD(5) * mrSges(6,2) + t291 * mrSges(6,3) - qJD(5) * t295 + t298 * t288;
t345 = t323 * t270 + t321 * t271;
t344 = t324 * t309 - t322 * t310;
t341 = -t317 ^ 2 - t315;
t337 = -qJDD(3) + t344;
t336 = t341 * mrSges(5,2);
t335 = -t321 * t270 + t323 * t271;
t333 = m(5) * t277 + t345;
t328 = -0.2e1 * qJD(4) * t340 - t337;
t276 = (t341 * pkin(6) + qJ(3)) * t325 + (t347 + pkin(2) + (pkin(3) + pkin(4)) * t319) * qJDD(2) - t328;
t327 = -m(6) * t276 + t291 * mrSges(6,1) - t292 * mrSges(6,2) + t298 * t294 - t299 * t295;
t279 = -t346 + (-pkin(2) + t331) * qJDD(2) + t328;
t326 = m(5) * t279 + t327;
t290 = -qJDD(2) * pkin(2) - t337 - t346;
t284 = Ifges(6,1) * t299 + Ifges(6,4) * t298 + Ifges(6,5) * qJD(5);
t283 = Ifges(6,4) * t299 + Ifges(6,2) * t298 + Ifges(6,6) * qJD(5);
t282 = Ifges(6,5) * t299 + Ifges(6,6) * t298 + Ifges(6,3) * qJD(5);
t266 = t332 * qJDD(2) + t325 * t336 + t326;
t265 = m(4) * t290 + (t341 * mrSges(4,3) + t336) * t325 + ((-mrSges(4,1) - mrSges(5,1)) * t319 + (mrSges(4,2) - mrSges(5,3)) * t317) * qJDD(2) + t326;
t264 = mrSges(6,2) * t276 - mrSges(6,3) * t272 + Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * qJDD(5) - qJD(5) * t283 + t298 * t282;
t263 = -mrSges(6,1) * t276 + mrSges(6,3) * t273 + Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * qJDD(5) + qJD(5) * t284 - t299 * t282;
t262 = m(4) * t281 + m(5) * t278 + t319 * t351 + t335;
t261 = m(4) * t280 - t317 * t351 - t333;
t1 = [t319 * t261 + t317 * t262 + (m(2) + m(3)) * t316; mrSges(3,1) * t344 - mrSges(3,2) * t343 + t317 * (mrSges(4,2) * t290 + mrSges(5,2) * t277 - mrSges(4,3) * t280 - mrSges(5,3) * t279 - pkin(6) * t345 - qJ(4) * t266 - t321 * t263 + t323 * t264) + t319 * (-mrSges(4,1) * t290 - mrSges(5,1) * t279 + mrSges(5,2) * t278 + mrSges(4,3) * t281 - pkin(3) * t266 - pkin(4) * t327 - pkin(6) * t335 - t323 * t263 - t321 * t264) - pkin(2) * t265 + qJ(3) * (-t317 * t261 + t319 * t262) + (Ifges(3,3) + (Ifges(4,2) + Ifges(5,3)) * t315 + ((Ifges(4,1) + Ifges(5,1)) * t317 + (Ifges(4,4) - Ifges(5,5)) * t352) * t317) * qJDD(2); t265; (qJDD(2) * mrSges(5,2) + qJD(2) * t305) * t317 + t333; mrSges(6,1) * t272 - mrSges(6,2) * t273 + Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * qJDD(5) + t299 * t283 - t298 * t284;];
tauJ = t1;
