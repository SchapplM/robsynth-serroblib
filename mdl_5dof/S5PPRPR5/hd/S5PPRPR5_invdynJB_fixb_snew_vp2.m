% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:27
% EndTime: 2019-12-31 17:33:28
% DurationCPUTime: 0.50s
% Computational Cost: add. (3498->150), mult. (5515->177), div. (0->0), fcn. (2748->6), ass. (0->68)
t358 = sin(pkin(7));
t359 = cos(pkin(7));
t347 = t358 * g(1) - t359 * g(2);
t345 = qJDD(2) - t347;
t348 = -t359 * g(1) - t358 * g(2);
t361 = sin(qJ(3));
t363 = cos(qJ(3));
t327 = t363 * t345 - t361 * t348;
t364 = qJD(3) ^ 2;
t370 = -t364 * qJ(4) + qJDD(4) - t327;
t385 = -pkin(3) - pkin(6);
t324 = t385 * qJDD(3) + t370;
t355 = g(3) - qJDD(1);
t360 = sin(qJ(5));
t362 = cos(qJ(5));
t320 = t362 * t324 - t360 * t355;
t342 = (mrSges(6,1) * t360 + mrSges(6,2) * t362) * qJD(3);
t375 = qJD(3) * qJD(5);
t344 = t362 * qJDD(3) - t360 * t375;
t377 = qJD(3) * t360;
t349 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t377;
t376 = qJD(3) * t362;
t317 = m(6) * t320 + qJDD(5) * mrSges(6,1) - t344 * mrSges(6,3) + qJD(5) * t349 - t342 * t376;
t321 = t360 * t324 + t362 * t355;
t343 = -t360 * qJDD(3) - t362 * t375;
t350 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t376;
t318 = m(6) * t321 - qJDD(5) * mrSges(6,2) + t343 * mrSges(6,3) - qJD(5) * t350 - t342 * t377;
t309 = t362 * t317 + t360 * t318;
t326 = -qJDD(3) * pkin(3) + t370;
t367 = -m(5) * t326 + t364 * mrSges(5,3) - t309;
t306 = qJDD(3) * mrSges(5,2) - t367;
t328 = t361 * t345 + t363 * t348;
t368 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t328;
t323 = t385 * t364 + t368;
t331 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t362 - Ifges(6,6) * t360) * qJD(3);
t333 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t362 - Ifges(6,4) * t360) * qJD(3);
t312 = -mrSges(6,1) * t323 + mrSges(6,3) * t321 + Ifges(6,4) * t344 + Ifges(6,2) * t343 + Ifges(6,6) * qJDD(5) + qJD(5) * t333 - t331 * t376;
t332 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t362 - Ifges(6,2) * t360) * qJD(3);
t313 = mrSges(6,2) * t323 - mrSges(6,3) * t320 + Ifges(6,1) * t344 + Ifges(6,4) * t343 + Ifges(6,5) * qJDD(5) - qJD(5) * t332 - t331 * t377;
t325 = t364 * pkin(3) - t368;
t371 = -m(6) * t323 + t343 * mrSges(6,1) - t344 * mrSges(6,2) - t349 * t377 - t350 * t376;
t314 = -m(5) * t325 + t364 * mrSges(5,2) + qJDD(3) * mrSges(5,3) - t371;
t387 = mrSges(4,1) * t327 - mrSges(4,2) * t328 + mrSges(5,2) * t326 - mrSges(5,3) * t325 - pkin(3) * t306 - pkin(6) * t309 + qJ(4) * t314 - t360 * t312 + t362 * t313 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3);
t386 = -m(4) - m(5);
t384 = mrSges(4,1) - mrSges(5,2);
t383 = -mrSges(2,2) + mrSges(3,3);
t382 = -Ifges(5,4) + Ifges(4,5);
t381 = Ifges(5,5) - Ifges(4,6);
t305 = m(4) * t327 - t364 * mrSges(4,2) + t384 * qJDD(3) + t367;
t311 = m(4) * t328 - t364 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t314;
t304 = t363 * t305 + t361 * t311;
t303 = m(3) * t345 + t304;
t301 = m(2) * t347 - t303;
t373 = -t361 * t305 + t363 * t311;
t372 = m(3) * t348 + t373;
t302 = m(2) * t348 + t372;
t379 = t359 * t301 + t358 * t302;
t378 = -t360 * t317 + t362 * t318;
t374 = -t358 * t301 + t359 * t302;
t307 = -t378 + (-m(3) + t386) * t355;
t369 = -m(2) * t355 + t307;
t366 = mrSges(6,1) * t320 - mrSges(6,2) * t321 + Ifges(6,5) * t344 + Ifges(6,6) * t343 + Ifges(6,3) * qJDD(5) + t332 * t376 + t333 * t377;
t308 = m(5) * t355 + t378;
t297 = mrSges(5,1) * t326 - mrSges(4,3) * t327 + pkin(4) * t309 - qJ(4) * t308 + t381 * t364 + (mrSges(4,2) - mrSges(5,3)) * t355 + t382 * qJDD(3) + t366;
t296 = -mrSges(5,1) * t325 + mrSges(4,3) * t328 - pkin(3) * t308 - pkin(4) * t371 - pkin(6) * t378 - t381 * qJDD(3) - t362 * t312 - t360 * t313 - t384 * t355 + t382 * t364;
t295 = mrSges(3,2) * t345 - mrSges(2,3) * t347 - pkin(5) * t304 - qJ(2) * t307 - t361 * t296 + t363 * t297 + t383 * t355;
t294 = -t361 * t297 - t363 * t296 + pkin(2) * t378 - pkin(5) * t373 - pkin(1) * t307 + (mrSges(3,2) + mrSges(2,3)) * t348 + (-pkin(2) * t386 + mrSges(2,1) + mrSges(3,1)) * t355;
t1 = [-m(1) * g(1) + t374; -m(1) * g(2) + t379; -m(1) * g(3) + t369; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t379 - t358 * t294 + t359 * t295; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t374 + t359 * t294 + t358 * t295; -mrSges(1,1) * g(2) + mrSges(2,1) * t347 - mrSges(3,1) * t345 + mrSges(1,2) * g(1) - pkin(1) * t303 - pkin(2) * t304 + qJ(2) * t372 + t383 * t348 - t387; t369; t303; t387; t306; t366;];
tauJB = t1;
