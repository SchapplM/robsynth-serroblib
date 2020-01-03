% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR5
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:31
% EndTime: 2019-12-31 17:03:31
% DurationCPUTime: 0.52s
% Computational Cost: add. (4942->149), mult. (5919->179), div. (0->0), fcn. (2586->6), ass. (0->67)
t391 = -m(3) - m(4);
t390 = -pkin(2) - pkin(6);
t389 = mrSges(3,1) - mrSges(4,2);
t388 = -Ifges(4,4) + Ifges(3,5);
t387 = Ifges(4,5) - Ifges(3,6);
t362 = qJD(1) + qJD(2);
t365 = sin(qJ(4));
t386 = t362 * t365;
t368 = cos(qJ(4));
t385 = t362 * t368;
t367 = sin(qJ(1));
t370 = cos(qJ(1));
t352 = t367 * g(1) - t370 * g(2);
t347 = qJDD(1) * pkin(1) + t352;
t353 = -t370 * g(1) - t367 * g(2);
t371 = qJD(1) ^ 2;
t348 = -t371 * pkin(1) + t353;
t366 = sin(qJ(2));
t369 = cos(qJ(2));
t329 = t369 * t347 - t366 * t348;
t360 = t362 ^ 2;
t361 = qJDD(1) + qJDD(2);
t377 = -t360 * qJ(3) + qJDD(3) - t329;
t324 = t390 * t361 + t377;
t320 = t365 * g(3) + t368 * t324;
t341 = (mrSges(5,1) * t365 + mrSges(5,2) * t368) * t362;
t383 = qJD(4) * t362;
t343 = t368 * t361 - t365 * t383;
t349 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t386;
t317 = m(5) * t320 + qJDD(4) * mrSges(5,1) - t343 * mrSges(5,3) + qJD(4) * t349 - t341 * t385;
t321 = -t368 * g(3) + t365 * t324;
t342 = -t365 * t361 - t368 * t383;
t350 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t385;
t318 = m(5) * t321 - qJDD(4) * mrSges(5,2) + t342 * mrSges(5,3) - qJD(4) * t350 - t341 * t386;
t308 = t368 * t317 + t365 * t318;
t327 = -t361 * pkin(2) + t377;
t376 = -m(4) * t327 + t360 * mrSges(4,3) - t308;
t304 = m(3) * t329 - t360 * mrSges(3,2) + t389 * t361 + t376;
t330 = t366 * t347 + t369 * t348;
t378 = t361 * qJ(3) + 0.2e1 * qJD(3) * t362 + t330;
t325 = t360 * pkin(2) - t378;
t323 = t390 * t360 + t378;
t379 = -m(5) * t323 + t342 * mrSges(5,1) - t343 * mrSges(5,2) - t349 * t386 - t350 * t385;
t374 = -m(4) * t325 + t360 * mrSges(4,2) + t361 * mrSges(4,3) - t379;
t311 = m(3) * t330 - t360 * mrSges(3,1) - t361 * mrSges(3,2) + t374;
t302 = t369 * t304 + t366 * t311;
t299 = m(2) * t352 + qJDD(1) * mrSges(2,1) - t371 * mrSges(2,2) + t302;
t381 = -t366 * t304 + t369 * t311;
t300 = m(2) * t353 - t371 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t381;
t384 = t370 * t299 + t367 * t300;
t382 = -t367 * t299 + t370 * t300;
t380 = t365 * t317 - t368 * t318;
t334 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t368 - Ifges(5,2) * t365) * t362;
t335 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t368 - Ifges(5,4) * t365) * t362;
t375 = mrSges(5,1) * t320 - mrSges(5,2) * t321 + Ifges(5,5) * t343 + Ifges(5,6) * t342 + Ifges(5,3) * qJDD(4) + t334 * t385 + t335 * t386;
t306 = t361 * mrSges(4,2) - t376;
t333 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t368 - Ifges(5,6) * t365) * t362;
t313 = -mrSges(5,1) * t323 + mrSges(5,3) * t321 + Ifges(5,4) * t343 + Ifges(5,2) * t342 + Ifges(5,6) * qJDD(4) + qJD(4) * t335 - t333 * t385;
t314 = mrSges(5,2) * t323 - mrSges(5,3) * t320 + Ifges(5,1) * t343 + Ifges(5,4) * t342 + Ifges(5,5) * qJDD(4) - qJD(4) * t334 - t333 * t386;
t373 = mrSges(3,1) * t329 - mrSges(3,2) * t330 + mrSges(4,2) * t327 - mrSges(4,3) * t325 - pkin(2) * t306 - pkin(6) * t308 + qJ(3) * t374 - t365 * t313 + t368 * t314 + (Ifges(3,3) + Ifges(4,1)) * t361;
t372 = mrSges(2,1) * t352 - mrSges(2,2) * t353 + Ifges(2,3) * qJDD(1) + pkin(1) * t302 + t373;
t307 = -m(4) * g(3) - t380;
t295 = mrSges(4,1) * t327 - mrSges(3,3) * t329 + pkin(3) * t308 - qJ(3) * t307 + t388 * t361 + t387 * t360 + (-mrSges(3,2) + mrSges(4,3)) * g(3) + t375;
t294 = -mrSges(4,1) * t325 + mrSges(3,3) * t330 - pkin(2) * t307 - pkin(3) * t379 + pkin(6) * t380 + t389 * g(3) - t368 * t313 - t365 * t314 + t388 * t360 - t387 * t361;
t293 = -mrSges(2,2) * g(3) - mrSges(2,3) * t352 + Ifges(2,5) * qJDD(1) - t371 * Ifges(2,6) - pkin(5) * t302 - t366 * t294 + t369 * t295;
t292 = Ifges(2,6) * qJDD(1) + t371 * Ifges(2,5) + mrSges(2,3) * t353 + t366 * t295 + t369 * t294 + pkin(1) * t380 + pkin(5) * t381 + (-pkin(1) * t391 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t382; -m(1) * g(2) + t384; (-m(1) - m(2) + t391) * g(3) - t380; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t384 - t367 * t292 + t370 * t293; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t382 + t370 * t292 + t367 * t293; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t372; t372; t373; t306; t375;];
tauJB = t1;
