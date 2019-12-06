% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:15
% EndTime: 2019-12-05 18:20:18
% DurationCPUTime: 2.46s
% Computational Cost: add. (34334->207), mult. (45201->270), div. (0->0), fcn. (25113->10), ass. (0->98)
t460 = 2 * qJD(4);
t428 = sin(qJ(1));
t431 = cos(qJ(1));
t409 = t431 * g(2) + t428 * g(3);
t403 = qJDD(1) * pkin(1) + t409;
t408 = t428 * g(2) - t431 * g(3);
t432 = qJD(1) ^ 2;
t404 = -t432 * pkin(1) + t408;
t427 = sin(qJ(2));
t430 = cos(qJ(2));
t389 = t430 * t403 - t427 * t404;
t419 = qJDD(1) + qJDD(2);
t384 = t419 * pkin(2) + t389;
t390 = t427 * t403 + t430 * t404;
t420 = (qJD(1) + qJD(2));
t418 = t420 ^ 2;
t385 = -t418 * pkin(2) + t390;
t423 = sin(pkin(8));
t425 = cos(pkin(8));
t380 = t423 * t384 + t425 * t385;
t378 = -t418 * pkin(3) + t419 * qJ(4) + t380;
t459 = (t420 * t460) + t378;
t422 = sin(pkin(9));
t458 = mrSges(5,2) * t422;
t456 = mrSges(5,3) * t419;
t455 = t422 * t420;
t426 = sin(qJ(5));
t454 = t422 * t426;
t429 = cos(qJ(5));
t453 = t422 * t429;
t424 = cos(pkin(9));
t452 = t424 * t419;
t451 = t424 * t420;
t421 = -g(1) + qJDD(3);
t450 = t424 * t421;
t374 = t422 * t421 + t459 * t424;
t397 = (-mrSges(5,1) * t424 + t458) * t420;
t439 = -pkin(4) * t424 - pkin(7) * t422;
t399 = t439 * t420;
t372 = t399 * t451 + t374;
t379 = t425 * t384 - t423 * t385;
t434 = -t418 * qJ(4) + qJDD(4) - t379;
t375 = (-pkin(3) + t439) * t419 + t434;
t369 = -t426 * t372 + t429 * t375;
t406 = qJD(5) - t451;
t447 = t420 * t454;
t392 = -t406 * mrSges(6,2) - mrSges(6,3) * t447;
t394 = (mrSges(6,1) * t426 + mrSges(6,2) * t429) * t455;
t448 = qJD(5) * t420;
t396 = (t419 * t429 - t426 * t448) * t422;
t405 = qJDD(5) - t452;
t446 = t420 * t453;
t367 = m(6) * t369 + t405 * mrSges(6,1) - t396 * mrSges(6,3) + t406 * t392 - t394 * t446;
t370 = t429 * t372 + t426 * t375;
t393 = t406 * mrSges(6,1) - mrSges(6,3) * t446;
t395 = (-t419 * t426 - t429 * t448) * t422;
t368 = m(6) * t370 - t405 * mrSges(6,2) + t395 * mrSges(6,3) - t406 * t393 - t394 * t447;
t440 = -t426 * t367 + t429 * t368;
t360 = m(5) * t374 + (t397 * t420 + t456) * t424 + t440;
t373 = -t459 * t422 + t450;
t371 = -t450 + (t378 + (t460 + t399) * t420) * t422;
t435 = -m(6) * t371 + t395 * mrSges(6,1) - t396 * mrSges(6,2);
t365 = m(5) * t373 + (-t456 + (-t392 * t426 - t393 * t429 - t397) * t420) * t422 + t435;
t441 = t424 * t360 - t422 * t365;
t354 = m(4) * t380 - t418 * mrSges(4,1) - t419 * mrSges(4,2) + t441;
t361 = t429 * t367 + t426 * t368;
t377 = -t419 * pkin(3) + t434;
t433 = -m(5) * t377 + mrSges(5,1) * t452 - t361 + (t422 ^ 2 + t424 ^ 2) * mrSges(5,3) * t418;
t357 = m(4) * t379 - t418 * mrSges(4,2) + (mrSges(4,1) - t458) * t419 + t433;
t349 = t423 * t354 + t425 * t357;
t347 = m(3) * t389 + t419 * mrSges(3,1) - t418 * mrSges(3,2) + t349;
t442 = t425 * t354 - t423 * t357;
t348 = m(3) * t390 - t418 * mrSges(3,1) - t419 * mrSges(3,2) + t442;
t342 = t430 * t347 + t427 * t348;
t355 = t422 * t360 + t424 * t365;
t445 = m(4) * t421 + t355;
t443 = -t427 * t347 + t430 * t348;
t340 = m(2) * t408 - t432 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t443;
t341 = m(2) * t409 + qJDD(1) * mrSges(2,1) - t432 * mrSges(2,2) + t342;
t444 = t431 * t340 - t428 * t341;
t438 = Ifges(5,1) * t422 + Ifges(5,4) * t424;
t437 = Ifges(5,5) * t422 + Ifges(5,6) * t424;
t436 = -t428 * t340 - t431 * t341;
t398 = t437 * t420;
t388 = Ifges(6,5) * t406 + (Ifges(6,1) * t429 - Ifges(6,4) * t426) * t455;
t387 = Ifges(6,6) * t406 + (Ifges(6,4) * t429 - Ifges(6,2) * t426) * t455;
t386 = Ifges(6,3) * t406 + (Ifges(6,5) * t429 - Ifges(6,6) * t426) * t455;
t363 = mrSges(6,2) * t371 - mrSges(6,3) * t369 + Ifges(6,1) * t396 + Ifges(6,4) * t395 + Ifges(6,5) * t405 - t386 * t447 - t406 * t387;
t362 = -mrSges(6,1) * t371 + mrSges(6,3) * t370 + Ifges(6,4) * t396 + Ifges(6,2) * t395 + Ifges(6,6) * t405 - t386 * t446 + t406 * t388;
t351 = Ifges(5,2) * t452 - mrSges(5,1) * t377 - mrSges(6,1) * t369 + mrSges(6,2) * t370 + mrSges(5,3) * t374 - Ifges(6,5) * t396 - Ifges(6,6) * t395 - Ifges(6,3) * t405 - pkin(4) * t361 + (Ifges(5,4) * t419 + (-t387 * t429 - t388 * t426 - t398) * t420) * t422;
t350 = mrSges(5,2) * t377 - mrSges(5,3) * t373 - pkin(7) * t361 - t426 * t362 + t429 * t363 + t398 * t451 + t438 * t419;
t343 = t418 * Ifges(4,5) - mrSges(4,1) * t421 + mrSges(4,3) * t380 - mrSges(5,1) * t373 + mrSges(5,2) * t374 - t426 * t363 - t429 * t362 - pkin(4) * t435 - pkin(7) * t440 - pkin(3) * t355 + (Ifges(4,6) - t437) * t419 + (-pkin(4) * (-t392 * t454 - t393 * t453) + (-t422 * (Ifges(5,4) * t422 + Ifges(5,2) * t424) + t424 * t438) * t420) * t420;
t338 = mrSges(4,2) * t421 - mrSges(4,3) * t379 + Ifges(4,5) * t419 - t418 * Ifges(4,6) - qJ(4) * t355 + t424 * t350 - t422 * t351;
t337 = -mrSges(3,2) * g(1) - mrSges(3,3) * t389 + Ifges(3,5) * t419 - t418 * Ifges(3,6) - qJ(3) * t349 + t425 * t338 - t423 * t343;
t336 = mrSges(3,1) * g(1) + mrSges(3,3) * t390 + t418 * Ifges(3,5) + Ifges(3,6) * t419 - pkin(2) * t445 + qJ(3) * t442 + t423 * t338 + t425 * t343;
t335 = -mrSges(2,2) * g(1) - mrSges(2,3) * t409 + Ifges(2,5) * qJDD(1) - t432 * Ifges(2,6) - pkin(6) * t342 - t427 * t336 + t430 * t337;
t334 = Ifges(2,6) * qJDD(1) + t432 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t408 + t427 * t337 + t430 * t336 - pkin(1) * (-m(3) * g(1) + t445) + pkin(6) * t443;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t445; -m(1) * g(2) + t436; -m(1) * g(3) + t444; pkin(1) * t342 + pkin(2) * t349 + mrSges(3,1) * t389 - mrSges(3,2) * t390 + qJ(4) * t441 + t422 * t350 + t424 * t351 + pkin(3) * t433 + mrSges(4,1) * t379 - mrSges(4,2) * t380 + mrSges(2,1) * t409 - mrSges(2,2) * t408 + Ifges(2,3) * qJDD(1) - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(3) * t458 + Ifges(3,3) + Ifges(4,3)) * t419; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t444 - t431 * t334 - t428 * t335; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t436 - t428 * t334 + t431 * t335;];
tauB = t1;
