% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRP6
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:02
% EndTime: 2019-12-31 16:46:03
% DurationCPUTime: 0.57s
% Computational Cost: add. (2141->170), mult. (4023->196), div. (0->0), fcn. (1555->4), ass. (0->67)
t380 = Ifges(4,1) + Ifges(5,1);
t373 = Ifges(4,4) + Ifges(5,4);
t371 = Ifges(4,5) + Ifges(5,5);
t379 = -Ifges(4,2) - Ifges(5,2);
t369 = Ifges(4,6) + Ifges(5,6);
t378 = (Ifges(4,3) + Ifges(5,3));
t343 = sin(qJ(1));
t345 = cos(qJ(1));
t332 = -t345 * g(1) - t343 * g(2);
t351 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t332;
t346 = qJD(1) ^ 2;
t376 = (-pkin(1) - pkin(5));
t304 = (t376 * t346) + t351;
t342 = sin(qJ(3));
t344 = cos(qJ(3));
t361 = qJD(1) * qJD(3);
t323 = -t342 * qJDD(1) - t344 * t361;
t357 = t342 * t361;
t324 = t344 * qJDD(1) - t357;
t363 = qJD(1) * t342;
t327 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t363;
t362 = qJD(1) * t344;
t330 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t362;
t328 = (qJD(3) * pkin(3)) - qJ(4) * t362;
t339 = t342 ^ 2;
t298 = t328 * t362 - t323 * pkin(3) + qJDD(4) + (-qJ(4) * t339 + t376) * t346 + t351;
t326 = -(qJD(3) * mrSges(5,2)) - mrSges(5,3) * t363;
t329 = (qJD(3) * mrSges(5,1)) - mrSges(5,3) * t362;
t353 = m(5) * t298 + t324 * mrSges(5,2) + t326 * t363 + t329 * t362;
t377 = -m(4) * t304 - t324 * mrSges(4,2) + (mrSges(4,1) + mrSges(5,1)) * t323 - t327 * t363 - t330 * t362 - t353;
t375 = mrSges(2,1) - mrSges(3,2);
t372 = Ifges(2,5) - Ifges(3,4);
t370 = (-Ifges(2,6) + Ifges(3,5));
t331 = t343 * g(1) - t345 * g(2);
t350 = -t346 * qJ(2) + qJDD(2) - t331;
t305 = t376 * qJDD(1) + t350;
t299 = t342 * g(3) + t344 * t305;
t321 = (mrSges(5,1) * t342 + mrSges(5,2) * t344) * qJD(1);
t354 = qJD(1) * (-t321 - (mrSges(4,1) * t342 + mrSges(4,2) * t344) * qJD(1));
t359 = -2 * qJD(1) * qJD(4);
t295 = t344 * t359 + (-t324 - t357) * qJ(4) + (-t342 * t344 * t346 + qJDD(3)) * pkin(3) + t299;
t358 = m(5) * t295 + qJDD(3) * mrSges(5,1) + qJD(3) * t326;
t290 = m(4) * t299 + qJDD(3) * mrSges(4,1) + qJD(3) * t327 + (-mrSges(4,3) - mrSges(5,3)) * t324 + t344 * t354 + t358;
t300 = -t344 * g(3) + t342 * t305;
t296 = -t339 * t346 * pkin(3) + t323 * qJ(4) - qJD(3) * t328 + t342 * t359 + t300;
t367 = m(5) * t296 + t323 * mrSges(5,3);
t291 = m(4) * t300 + t323 * mrSges(4,3) + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t329 - t330) * qJD(3) + t342 * t354 + t367;
t285 = t344 * t290 + t342 * t291;
t307 = -qJDD(1) * pkin(1) + t350;
t348 = -m(3) * t307 + (t346 * mrSges(3,3)) - t285;
t283 = m(2) * t331 - (t346 * mrSges(2,2)) + t375 * qJDD(1) + t348;
t306 = t346 * pkin(1) - t351;
t347 = -m(3) * t306 + (t346 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t377;
t288 = m(2) * t332 - (t346 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t347;
t368 = t345 * t283 + t343 * t288;
t366 = -(t378 * qJD(3)) + (t369 * t342 - t371 * t344) * qJD(1);
t365 = t369 * qJD(3) + (t379 * t342 + t373 * t344) * qJD(1);
t364 = t371 * qJD(3) + (-t373 * t342 + t380 * t344) * qJD(1);
t356 = -t343 * t283 + t345 * t288;
t355 = -t342 * t290 + t344 * t291;
t292 = -t324 * mrSges(5,3) - t321 * t362 + t358;
t284 = -m(3) * g(3) + t355;
t281 = mrSges(4,2) * t304 + mrSges(5,2) * t298 - mrSges(4,3) * t299 - mrSges(5,3) * t295 - qJ(4) * t292 - t365 * qJD(3) + t371 * qJDD(3) + t373 * t323 + t380 * t324 + t366 * t363;
t280 = -mrSges(4,1) * t304 + mrSges(4,3) * t300 - mrSges(5,1) * t298 + mrSges(5,3) * t296 - pkin(3) * t353 + qJ(4) * t367 + t373 * t324 + ((pkin(3) * mrSges(5,1)) - t379) * t323 + (-qJ(4) * mrSges(5,2) + t369) * qJDD(3) + (-qJ(4) * t329 + t364) * qJD(3) + (-qJ(4) * t321 * t342 + t366 * t344) * qJD(1);
t279 = mrSges(3,1) * t307 + mrSges(4,1) * t299 + mrSges(5,1) * t295 - mrSges(4,2) * t300 - mrSges(5,2) * t296 - mrSges(2,3) * t331 + pkin(2) * t285 + pkin(3) * t292 - qJ(2) * t284 + (t370 * t346) + t371 * t324 + t369 * t323 + t378 * qJDD(3) + t372 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t364 * t342 + t365 * t344) * qJD(1);
t278 = -mrSges(3,1) * t306 + mrSges(2,3) * t332 - pkin(1) * t284 - pkin(2) * t377 - pkin(5) * t355 + t375 * g(3) - t370 * qJDD(1) - t344 * t280 - t342 * t281 + t372 * t346;
t1 = [-m(1) * g(1) + t356; -m(1) * g(2) + t368; (-m(1) - m(2) - m(3)) * g(3) + t355; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t368 - t343 * t278 + t345 * t279; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t356 + t345 * t278 + t343 * t279; -t342 * t280 - pkin(5) * t285 + mrSges(2,1) * t331 - mrSges(2,2) * t332 + pkin(1) * t348 + qJ(2) * t347 + mrSges(3,2) * t307 - mrSges(3,3) * t306 + t344 * t281 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
