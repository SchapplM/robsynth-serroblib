% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:39
% EndTime: 2019-12-31 16:25:40
% DurationCPUTime: 0.44s
% Computational Cost: add. (2584->137), mult. (4258->168), div. (0->0), fcn. (2046->6), ass. (0->63)
t332 = sin(pkin(6));
t352 = cos(pkin(6));
t321 = -t352 * g(1) - t332 * g(2);
t329 = -g(3) + qJDD(1);
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t303 = -t334 * t321 + t336 * t329;
t337 = qJD(2) ^ 2;
t342 = -t337 * qJ(3) + qJDD(3) - t303;
t357 = -pkin(2) - pkin(5);
t300 = t357 * qJDD(2) + t342;
t320 = t332 * g(1) - t352 * g(2);
t333 = sin(qJ(4));
t335 = cos(qJ(4));
t296 = t335 * t300 + t333 * t320;
t317 = (mrSges(5,1) * t333 + mrSges(5,2) * t335) * qJD(2);
t347 = qJD(2) * qJD(4);
t319 = t335 * qJDD(2) - t333 * t347;
t349 = qJD(2) * t333;
t322 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t349;
t348 = qJD(2) * t335;
t293 = m(5) * t296 + qJDD(4) * mrSges(5,1) - t319 * mrSges(5,3) + qJD(4) * t322 - t317 * t348;
t297 = t333 * t300 - t335 * t320;
t318 = -t333 * qJDD(2) - t335 * t347;
t323 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t348;
t294 = m(5) * t297 - qJDD(4) * mrSges(5,2) + t318 * mrSges(5,3) - qJD(4) * t323 - t317 * t349;
t284 = t335 * t293 + t333 * t294;
t302 = -qJDD(2) * pkin(2) + t342;
t340 = -m(4) * t302 + t337 * mrSges(4,3) - t284;
t281 = qJDD(2) * mrSges(4,2) - t340;
t304 = t336 * t321 + t334 * t329;
t341 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t304;
t299 = t357 * t337 + t341;
t307 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t335 - Ifges(5,6) * t333) * qJD(2);
t309 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t335 - Ifges(5,4) * t333) * qJD(2);
t288 = -mrSges(5,1) * t299 + mrSges(5,3) * t297 + Ifges(5,4) * t319 + Ifges(5,2) * t318 + Ifges(5,6) * qJDD(4) + qJD(4) * t309 - t307 * t348;
t308 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t335 - Ifges(5,2) * t333) * qJD(2);
t289 = mrSges(5,2) * t299 - mrSges(5,3) * t296 + Ifges(5,1) * t319 + Ifges(5,4) * t318 + Ifges(5,5) * qJDD(4) - qJD(4) * t308 - t307 * t349;
t301 = t337 * pkin(2) - t341;
t343 = -m(5) * t299 + t318 * mrSges(5,1) - t319 * mrSges(5,2) - t322 * t349 - t323 * t348;
t290 = -m(4) * t301 + t337 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t343;
t359 = mrSges(3,1) * t303 - mrSges(3,2) * t304 + mrSges(4,2) * t302 - mrSges(4,3) * t301 - pkin(2) * t281 - pkin(5) * t284 + qJ(3) * t290 - t333 * t288 + t335 * t289 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t358 = m(3) + m(4);
t356 = mrSges(3,1) - mrSges(4,2);
t355 = -Ifges(4,4) + Ifges(3,5);
t354 = Ifges(4,5) - Ifges(3,6);
t279 = m(3) * t303 - t337 * mrSges(3,2) + t356 * qJDD(2) + t340;
t287 = m(3) * t304 - t337 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t290;
t344 = -t334 * t279 + t336 * t287;
t276 = m(2) * t321 + t344;
t350 = -t333 * t293 + t335 * t294;
t282 = (m(2) + t358) * t320 - t350;
t351 = t332 * t276 + t352 * t282;
t277 = t336 * t279 + t334 * t287;
t346 = m(2) * t329 + t277;
t345 = t352 * t276 - t332 * t282;
t339 = mrSges(5,1) * t296 - mrSges(5,2) * t297 + Ifges(5,5) * t319 + Ifges(5,6) * t318 + Ifges(5,3) * qJDD(4) + t308 * t348 + t309 * t349;
t283 = -m(4) * t320 + t350;
t273 = mrSges(4,1) * t302 - mrSges(3,3) * t303 + pkin(3) * t284 - qJ(3) * t283 + t354 * t337 + (-mrSges(3,2) + mrSges(4,3)) * t320 + t355 * qJDD(2) + t339;
t272 = -mrSges(4,1) * t301 + mrSges(3,3) * t304 - pkin(2) * t283 - pkin(3) * t343 - pkin(5) * t350 - t354 * qJDD(2) - t335 * t288 - t333 * t289 + t356 * t320 + t355 * t337;
t271 = -mrSges(2,1) * t329 + mrSges(2,3) * t321 - pkin(1) * t277 - t359;
t270 = mrSges(2,2) * t329 - mrSges(2,3) * t320 - pkin(4) * t277 - t334 * t272 + t336 * t273;
t1 = [-m(1) * g(1) + t345; -m(1) * g(2) + t351; -m(1) * g(3) + t346; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t351 + t352 * t270 - t332 * t271; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t345 + t332 * t270 + t352 * t271; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t321 + t334 * t273 + t336 * t272 - pkin(1) * t350 + pkin(4) * t344 + (pkin(1) * t358 + mrSges(2,1)) * t320; t346; t359; t281; t339;];
tauJB = t1;
