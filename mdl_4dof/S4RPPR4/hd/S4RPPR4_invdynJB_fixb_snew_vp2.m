% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:51
% EndTime: 2019-12-31 16:38:51
% DurationCPUTime: 0.48s
% Computational Cost: add. (3342->146), mult. (5623->177), div. (0->0), fcn. (2458->6), ass. (0->64)
t378 = -pkin(2) - pkin(5);
t377 = mrSges(3,1) - mrSges(4,2);
t376 = -Ifges(4,4) + Ifges(3,5);
t375 = Ifges(4,5) - Ifges(3,6);
t357 = sin(qJ(1));
t359 = cos(qJ(1));
t341 = t357 * g(1) - g(2) * t359;
t333 = qJDD(1) * pkin(1) + t341;
t342 = -g(1) * t359 - g(2) * t357;
t360 = qJD(1) ^ 2;
t335 = -pkin(1) * t360 + t342;
t354 = sin(pkin(6));
t355 = cos(pkin(6));
t318 = t355 * t333 - t354 * t335;
t366 = -t360 * qJ(3) + qJDD(3) - t318;
t313 = t378 * qJDD(1) + t366;
t351 = -g(3) + qJDD(2);
t356 = sin(qJ(4));
t358 = cos(qJ(4));
t309 = t313 * t358 - t351 * t356;
t334 = (mrSges(5,1) * t356 + mrSges(5,2) * t358) * qJD(1);
t371 = qJD(1) * qJD(4);
t337 = qJDD(1) * t358 - t356 * t371;
t373 = qJD(1) * t356;
t339 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t373;
t372 = qJD(1) * t358;
t306 = m(5) * t309 + qJDD(4) * mrSges(5,1) - t337 * mrSges(5,3) + qJD(4) * t339 - t334 * t372;
t310 = t313 * t356 + t351 * t358;
t336 = -qJDD(1) * t356 - t358 * t371;
t340 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t372;
t307 = m(5) * t310 - qJDD(4) * mrSges(5,2) + t336 * mrSges(5,3) - qJD(4) * t340 - t334 * t373;
t297 = t358 * t306 + t356 * t307;
t316 = -qJDD(1) * pkin(2) + t366;
t364 = -m(4) * t316 + t360 * mrSges(4,3) - t297;
t292 = m(3) * t318 - t360 * mrSges(3,2) + t377 * qJDD(1) + t364;
t319 = t354 * t333 + t355 * t335;
t365 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t319;
t314 = pkin(2) * t360 - t365;
t312 = t378 * t360 + t365;
t367 = -m(5) * t312 + mrSges(5,1) * t336 - t337 * mrSges(5,2) - t339 * t373 - t340 * t372;
t362 = -m(4) * t314 + t360 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t367;
t300 = m(3) * t319 - mrSges(3,1) * t360 - qJDD(1) * mrSges(3,2) + t362;
t290 = t355 * t292 + t354 * t300;
t287 = m(2) * t341 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t360 + t290;
t369 = -t292 * t354 + t355 * t300;
t288 = m(2) * t342 - mrSges(2,1) * t360 - qJDD(1) * mrSges(2,2) + t369;
t374 = t359 * t287 + t357 * t288;
t370 = -t287 * t357 + t359 * t288;
t368 = -t356 * t306 + t358 * t307;
t296 = m(4) * t351 + t368;
t295 = m(3) * t351 + t296;
t326 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t358 - Ifges(5,2) * t356) * qJD(1);
t327 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t358 - Ifges(5,4) * t356) * qJD(1);
t363 = mrSges(5,1) * t309 - mrSges(5,2) * t310 + Ifges(5,5) * t337 + Ifges(5,6) * t336 + Ifges(5,3) * qJDD(4) + t326 * t372 + t327 * t373;
t294 = qJDD(1) * mrSges(4,2) - t364;
t325 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t358 - Ifges(5,6) * t356) * qJD(1);
t302 = -mrSges(5,1) * t312 + mrSges(5,3) * t310 + Ifges(5,4) * t337 + Ifges(5,2) * t336 + Ifges(5,6) * qJDD(4) + qJD(4) * t327 - t325 * t372;
t303 = mrSges(5,2) * t312 - mrSges(5,3) * t309 + Ifges(5,1) * t337 + Ifges(5,4) * t336 + Ifges(5,5) * qJDD(4) - qJD(4) * t326 - t325 * t373;
t361 = mrSges(2,1) * t341 + mrSges(3,1) * t318 - mrSges(2,2) * t342 - mrSges(3,2) * t319 + mrSges(4,2) * t316 - mrSges(4,3) * t314 + pkin(1) * t290 - pkin(2) * t294 - pkin(5) * t297 + qJ(3) * t362 - t302 * t356 + t358 * t303 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1);
t283 = mrSges(4,1) * t316 - mrSges(3,3) * t318 + pkin(3) * t297 - qJ(3) * t296 + t375 * t360 + (mrSges(3,2) - mrSges(4,3)) * t351 + t376 * qJDD(1) + t363;
t282 = -mrSges(4,1) * t314 + mrSges(3,3) * t319 - pkin(2) * t296 - pkin(3) * t367 - pkin(5) * t368 - t375 * qJDD(1) - t358 * t302 - t356 * t303 - t377 * t351 + t376 * t360;
t281 = -mrSges(2,2) * g(3) - mrSges(2,3) * t341 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t360 - qJ(2) * t290 - t282 * t354 + t283 * t355;
t280 = mrSges(2,1) * g(3) + mrSges(2,3) * t342 + t360 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t295 + qJ(2) * t369 + t355 * t282 + t354 * t283;
t1 = [-m(1) * g(1) + t370; -m(1) * g(2) + t374; (-m(1) - m(2)) * g(3) + t295; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t374 - t357 * t280 + t359 * t281; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t370 + t359 * t280 + t357 * t281; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361; t361; t295; t294; t363;];
tauJB = t1;
