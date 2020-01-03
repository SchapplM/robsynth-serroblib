% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:39
% EndTime: 2019-12-31 16:33:39
% DurationCPUTime: 0.66s
% Computational Cost: add. (7031->144), mult. (8905->187), div. (0->0), fcn. (5046->8), ass. (0->67)
t364 = sin(pkin(7));
t386 = cos(pkin(7));
t354 = -t386 * g(1) - t364 * g(2);
t363 = -g(3) + qJDD(1);
t367 = sin(qJ(2));
t370 = cos(qJ(2));
t340 = -t367 * t354 + t370 * t363;
t338 = qJDD(2) * pkin(2) + t340;
t341 = t370 * t354 + t367 * t363;
t371 = qJD(2) ^ 2;
t339 = -t371 * pkin(2) + t341;
t366 = sin(qJ(3));
t369 = cos(qJ(3));
t335 = t366 * t338 + t369 * t339;
t362 = qJD(2) + qJD(3);
t360 = t362 ^ 2;
t361 = qJDD(2) + qJDD(3);
t332 = -t360 * pkin(3) + t361 * pkin(6) + t335;
t353 = t364 * g(1) - t386 * g(2);
t365 = sin(qJ(4));
t368 = cos(qJ(4));
t329 = -t365 * t332 - t368 * t353;
t330 = t368 * t332 - t365 * t353;
t343 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t365 + Ifges(5,2) * t368) * t362;
t344 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t365 + Ifges(5,4) * t368) * t362;
t382 = qJD(4) * t362;
t348 = t365 * t361 + t368 * t382;
t349 = t368 * t361 - t365 * t382;
t388 = mrSges(5,1) * t329 - mrSges(5,2) * t330 + Ifges(5,5) * t348 + Ifges(5,6) * t349 + Ifges(5,3) * qJDD(4) + (t343 * t365 - t344 * t368) * t362;
t387 = m(3) + m(4);
t385 = t362 * t365;
t384 = t362 * t368;
t347 = (-mrSges(5,1) * t368 + mrSges(5,2) * t365) * t362;
t352 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t384;
t325 = m(5) * t329 + qJDD(4) * mrSges(5,1) - t348 * mrSges(5,3) + qJD(4) * t352 - t347 * t385;
t351 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t385;
t326 = m(5) * t330 - qJDD(4) * mrSges(5,2) + t349 * mrSges(5,3) - qJD(4) * t351 + t347 * t384;
t377 = -t365 * t325 + t368 * t326;
t310 = m(4) * t335 - t360 * mrSges(4,1) - t361 * mrSges(4,2) + t377;
t334 = t369 * t338 - t366 * t339;
t331 = -t361 * pkin(3) - t360 * pkin(6) - t334;
t374 = -m(5) * t331 + t349 * mrSges(5,1) - t348 * mrSges(5,2) - t351 * t385 + t352 * t384;
t321 = m(4) * t334 + t361 * mrSges(4,1) - t360 * mrSges(4,2) + t374;
t307 = t366 * t310 + t369 * t321;
t305 = m(3) * t340 + qJDD(2) * mrSges(3,1) - t371 * mrSges(3,2) + t307;
t378 = t369 * t310 - t366 * t321;
t306 = m(3) * t341 - t371 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t378;
t379 = -t367 * t305 + t370 * t306;
t298 = m(2) * t354 + t379;
t314 = t368 * t325 + t365 * t326;
t312 = (m(2) + t387) * t353 - t314;
t383 = t364 * t298 + t386 * t312;
t299 = t370 * t305 + t367 * t306;
t381 = m(2) * t363 + t299;
t380 = t386 * t298 - t364 * t312;
t342 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t365 + Ifges(5,6) * t368) * t362;
t318 = -mrSges(5,1) * t331 + mrSges(5,3) * t330 + Ifges(5,4) * t348 + Ifges(5,2) * t349 + Ifges(5,6) * qJDD(4) + qJD(4) * t344 - t342 * t385;
t319 = mrSges(5,2) * t331 - mrSges(5,3) * t329 + Ifges(5,1) * t348 + Ifges(5,4) * t349 + Ifges(5,5) * qJDD(4) - qJD(4) * t343 + t342 * t384;
t375 = mrSges(4,1) * t334 - mrSges(4,2) * t335 + Ifges(4,3) * t361 + pkin(3) * t374 + pkin(6) * t377 + t368 * t318 + t365 * t319;
t372 = mrSges(3,1) * t340 - mrSges(3,2) * t341 + Ifges(3,3) * qJDD(2) + pkin(2) * t307 + t375;
t301 = mrSges(4,1) * t353 + mrSges(4,3) * t335 + t360 * Ifges(4,5) + Ifges(4,6) * t361 - pkin(3) * t314 - t388;
t300 = -mrSges(4,2) * t353 - mrSges(4,3) * t334 + Ifges(4,5) * t361 - t360 * Ifges(4,6) - pkin(6) * t314 - t365 * t318 + t368 * t319;
t295 = -mrSges(3,2) * t353 - mrSges(3,3) * t340 + Ifges(3,5) * qJDD(2) - t371 * Ifges(3,6) - pkin(5) * t307 + t369 * t300 - t366 * t301;
t294 = Ifges(3,6) * qJDD(2) + t371 * Ifges(3,5) + mrSges(3,1) * t353 + mrSges(3,3) * t341 + t366 * t300 + t369 * t301 - pkin(2) * (-m(4) * t353 + t314) + pkin(5) * t378;
t293 = -mrSges(2,1) * t363 + mrSges(2,3) * t354 - pkin(1) * t299 - t372;
t292 = mrSges(2,2) * t363 - mrSges(2,3) * t353 - pkin(4) * t299 - t367 * t294 + t370 * t295;
t1 = [-m(1) * g(1) + t380; -m(1) * g(2) + t383; -m(1) * g(3) + t381; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t383 + t386 * t292 - t364 * t293; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t380 + t364 * t292 + t386 * t293; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t354 + t367 * t295 + t370 * t294 - pkin(1) * t314 + pkin(4) * t379 + (pkin(1) * t387 + mrSges(2,1)) * t353; t381; t372; t375; t388;];
tauJB = t1;
