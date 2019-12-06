% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:46
% EndTime: 2019-12-05 15:12:47
% DurationCPUTime: 1.61s
% Computational Cost: add. (20170->163), mult. (26280->212), div. (0->0), fcn. (18190->10), ass. (0->79)
t427 = sin(pkin(8));
t429 = cos(pkin(8));
t418 = t427 * g(1) - t429 * g(2);
t417 = qJDD(2) - t418;
t419 = -t429 * g(1) - t427 * g(2);
t425 = -g(3) + qJDD(1);
t426 = sin(pkin(9));
t428 = cos(pkin(9));
t404 = -t426 * t419 + t428 * t425;
t405 = t428 * t419 + t426 * t425;
t432 = sin(qJ(3));
t435 = cos(qJ(3));
t399 = t435 * t404 - t432 * t405;
t397 = qJDD(3) * pkin(3) + t399;
t400 = t432 * t404 + t435 * t405;
t436 = qJD(3) ^ 2;
t398 = -t436 * pkin(3) + t400;
t431 = sin(qJ(4));
t434 = cos(qJ(4));
t394 = t431 * t397 + t434 * t398;
t424 = qJD(3) + qJD(4);
t422 = t424 ^ 2;
t423 = qJDD(3) + qJDD(4);
t391 = -t422 * pkin(4) + t423 * pkin(7) + t394;
t430 = sin(qJ(5));
t433 = cos(qJ(5));
t388 = -t430 * t391 + t433 * t417;
t411 = (-mrSges(6,1) * t433 + mrSges(6,2) * t430) * t424;
t449 = qJD(5) * t424;
t412 = t430 * t423 + t433 * t449;
t451 = t424 * t433;
t416 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t451;
t452 = t424 * t430;
t384 = m(6) * t388 + qJDD(5) * mrSges(6,1) - t412 * mrSges(6,3) + qJD(5) * t416 - t411 * t452;
t389 = t433 * t391 + t430 * t417;
t413 = t433 * t423 - t430 * t449;
t415 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t452;
t385 = m(6) * t389 - qJDD(5) * mrSges(6,2) + t413 * mrSges(6,3) - qJD(5) * t415 + t411 * t451;
t373 = t433 * t384 + t430 * t385;
t448 = m(5) * t417 + t373;
t371 = (m(3) + m(4)) * t417 + t448;
t407 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t430 + Ifges(6,2) * t433) * t424;
t408 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t430 + Ifges(6,4) * t433) * t424;
t454 = mrSges(6,1) * t388 - mrSges(6,2) * t389 + Ifges(6,5) * t412 + Ifges(6,6) * t413 + Ifges(6,3) * qJDD(5) + (t407 * t430 - t408 * t433) * t424;
t442 = -t430 * t384 + t433 * t385;
t368 = m(5) * t394 - t422 * mrSges(5,1) - t423 * mrSges(5,2) + t442;
t393 = t434 * t397 - t431 * t398;
t390 = -t423 * pkin(4) - t422 * pkin(7) - t393;
t439 = -m(6) * t390 + t413 * mrSges(6,1) - t412 * mrSges(6,2) - t415 * t452 + t416 * t451;
t380 = m(5) * t393 + t423 * mrSges(5,1) - t422 * mrSges(5,2) + t439;
t365 = t431 * t368 + t434 * t380;
t363 = m(4) * t399 + qJDD(3) * mrSges(4,1) - t436 * mrSges(4,2) + t365;
t443 = t434 * t368 - t431 * t380;
t364 = m(4) * t400 - t436 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t443;
t357 = t435 * t363 + t432 * t364;
t355 = m(3) * t404 + t357;
t444 = -t432 * t363 + t435 * t364;
t356 = m(3) * t405 + t444;
t445 = -t426 * t355 + t428 * t356;
t348 = m(2) * t419 + t445;
t370 = m(2) * t418 - t371;
t450 = t427 * t348 + t429 * t370;
t349 = t428 * t355 + t426 * t356;
t447 = m(2) * t425 + t349;
t446 = t429 * t348 - t427 * t370;
t406 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t430 + Ifges(6,6) * t433) * t424;
t377 = -mrSges(6,1) * t390 + mrSges(6,3) * t389 + Ifges(6,4) * t412 + Ifges(6,2) * t413 + Ifges(6,6) * qJDD(5) + qJD(5) * t408 - t406 * t452;
t378 = mrSges(6,2) * t390 - mrSges(6,3) * t388 + Ifges(6,1) * t412 + Ifges(6,4) * t413 + Ifges(6,5) * qJDD(5) - qJD(5) * t407 + t406 * t451;
t440 = mrSges(5,1) * t393 - mrSges(5,2) * t394 + Ifges(5,3) * t423 + pkin(4) * t439 + pkin(7) * t442 + t433 * t377 + t430 * t378;
t437 = mrSges(4,1) * t399 - mrSges(4,2) * t400 + Ifges(4,3) * qJDD(3) + pkin(3) * t365 + t440;
t359 = -mrSges(5,1) * t417 + mrSges(5,3) * t394 + t422 * Ifges(5,5) + Ifges(5,6) * t423 - pkin(4) * t373 - t454;
t358 = mrSges(5,2) * t417 - mrSges(5,3) * t393 + Ifges(5,5) * t423 - t422 * Ifges(5,6) - pkin(7) * t373 - t430 * t377 + t433 * t378;
t351 = mrSges(4,2) * t417 - mrSges(4,3) * t399 + Ifges(4,5) * qJDD(3) - t436 * Ifges(4,6) - pkin(6) * t365 + t434 * t358 - t431 * t359;
t350 = -mrSges(4,1) * t417 + mrSges(4,3) * t400 + t436 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t448 + pkin(6) * t443 + t431 * t358 + t434 * t359;
t345 = -mrSges(2,1) * t425 - mrSges(3,1) * t404 + mrSges(3,2) * t405 + mrSges(2,3) * t419 - pkin(1) * t349 - pkin(2) * t357 - t437;
t344 = mrSges(3,2) * t417 - mrSges(3,3) * t404 - pkin(5) * t357 - t432 * t350 + t435 * t351;
t343 = -mrSges(3,1) * t417 + mrSges(3,3) * t405 + t432 * t351 + t435 * t350 - pkin(2) * (m(4) * t417 + t448) + pkin(5) * t444;
t342 = mrSges(2,2) * t425 - mrSges(2,3) * t418 - qJ(2) * t349 - t426 * t343 + t428 * t344;
t1 = [-m(1) * g(1) + t446; -m(1) * g(2) + t450; -m(1) * g(3) + t447; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t450 + t429 * t342 - t427 * t345; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t446 + t427 * t342 + t429 * t345; -mrSges(1,1) * g(2) + mrSges(2,1) * t418 + mrSges(1,2) * g(1) - mrSges(2,2) * t419 - pkin(1) * t371 + qJ(2) * t445 + t428 * t343 + t426 * t344; t447; t371; t437; t440; t454;];
tauJB = t1;
