% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:33
% EndTime: 2020-01-03 12:07:36
% DurationCPUTime: 2.48s
% Computational Cost: add. (46436->187), mult. (53945->233), div. (0->0), fcn. (28876->10), ass. (0->85)
t492 = sin(qJ(1));
t496 = cos(qJ(1));
t473 = -t496 * g(2) - t492 * g(3);
t469 = qJDD(1) * pkin(1) + t473;
t472 = -t492 * g(2) + t496 * g(3);
t497 = qJD(1) ^ 2;
t470 = -t497 * pkin(1) + t472;
t491 = sin(qJ(2));
t495 = cos(qJ(2));
t454 = t495 * t469 - t491 * t470;
t483 = qJDD(1) + qJDD(2);
t451 = t483 * pkin(2) + t454;
t455 = t491 * t469 + t495 * t470;
t484 = qJD(1) + qJD(2);
t482 = t484 ^ 2;
t452 = -t482 * pkin(2) + t455;
t490 = sin(qJ(3));
t494 = cos(qJ(3));
t446 = t494 * t451 - t490 * t452;
t478 = qJDD(3) + t483;
t443 = t478 * pkin(3) + t446;
t447 = t490 * t451 + t494 * t452;
t479 = qJD(3) + t484;
t477 = t479 ^ 2;
t444 = -t477 * pkin(3) + t447;
t487 = sin(pkin(9));
t488 = cos(pkin(9));
t440 = t487 * t443 + t488 * t444;
t437 = -t477 * pkin(4) + t478 * pkin(8) + t440;
t486 = -g(1) + qJDD(4);
t489 = sin(qJ(5));
t493 = cos(qJ(5));
t434 = -t489 * t437 + t493 * t486;
t435 = t493 * t437 + t489 * t486;
t457 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t489 + Ifges(6,2) * t493) * t479;
t458 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t489 + Ifges(6,4) * t493) * t479;
t509 = qJD(5) * t479;
t462 = t489 * t478 + t493 * t509;
t463 = t493 * t478 - t489 * t509;
t514 = mrSges(6,1) * t434 - mrSges(6,2) * t435 + Ifges(6,5) * t462 + Ifges(6,6) * t463 + Ifges(6,3) * qJDD(5) + (t457 * t489 - t458 * t493) * t479;
t513 = -m(3) - m(4);
t512 = t479 * t489;
t511 = t479 * t493;
t461 = (-mrSges(6,1) * t493 + mrSges(6,2) * t489) * t479;
t468 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t511;
t432 = m(6) * t434 + qJDD(5) * mrSges(6,1) - t462 * mrSges(6,3) + qJD(5) * t468 - t461 * t512;
t467 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t512;
t433 = m(6) * t435 - qJDD(5) * mrSges(6,2) + t463 * mrSges(6,3) - qJD(5) * t467 + t461 * t511;
t504 = -t489 * t432 + t493 * t433;
t418 = m(5) * t440 - t477 * mrSges(5,1) - t478 * mrSges(5,2) + t504;
t439 = t488 * t443 - t487 * t444;
t436 = -t478 * pkin(4) - t477 * pkin(8) - t439;
t502 = -m(6) * t436 + t463 * mrSges(6,1) - t462 * mrSges(6,2) - t467 * t512 + t468 * t511;
t427 = m(5) * t439 + t478 * mrSges(5,1) - t477 * mrSges(5,2) + t502;
t415 = t487 * t418 + t488 * t427;
t411 = m(4) * t446 + t478 * mrSges(4,1) - t477 * mrSges(4,2) + t415;
t505 = t488 * t418 - t487 * t427;
t412 = m(4) * t447 - t477 * mrSges(4,1) - t478 * mrSges(4,2) + t505;
t406 = t494 * t411 + t490 * t412;
t403 = m(3) * t454 + t483 * mrSges(3,1) - t482 * mrSges(3,2) + t406;
t506 = -t490 * t411 + t494 * t412;
t404 = m(3) * t455 - t482 * mrSges(3,1) - t483 * mrSges(3,2) + t506;
t507 = -t491 * t403 + t495 * t404;
t394 = m(2) * t472 - t497 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t507;
t397 = t495 * t403 + t491 * t404;
t395 = m(2) * t473 + qJDD(1) * mrSges(2,1) - t497 * mrSges(2,2) + t397;
t510 = t492 * t394 + t496 * t395;
t421 = t493 * t432 + t489 * t433;
t419 = m(5) * t486 + t421;
t508 = -t496 * t394 + t492 * t395;
t456 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t489 + Ifges(6,6) * t493) * t479;
t424 = -mrSges(6,1) * t436 + mrSges(6,3) * t435 + Ifges(6,4) * t462 + Ifges(6,2) * t463 + Ifges(6,6) * qJDD(5) + qJD(5) * t458 - t456 * t512;
t425 = mrSges(6,2) * t436 - mrSges(6,3) * t434 + Ifges(6,1) * t462 + Ifges(6,4) * t463 + Ifges(6,5) * qJDD(5) - qJD(5) * t457 + t456 * t511;
t500 = mrSges(4,1) * t446 + mrSges(5,1) * t439 - mrSges(4,2) * t447 - mrSges(5,2) * t440 + pkin(3) * t415 + pkin(4) * t502 + pkin(8) * t504 + t493 * t424 + t489 * t425 + (Ifges(5,3) + Ifges(4,3)) * t478;
t499 = mrSges(3,1) * t454 - mrSges(3,2) * t455 + Ifges(3,3) * t483 + pkin(2) * t406 + t500;
t498 = mrSges(2,1) * t473 - mrSges(2,2) * t472 + Ifges(2,3) * qJDD(1) + pkin(1) * t397 + t499;
t413 = -mrSges(5,1) * t486 + mrSges(5,3) * t440 + t477 * Ifges(5,5) + Ifges(5,6) * t478 - pkin(4) * t421 - t514;
t407 = mrSges(5,2) * t486 - mrSges(5,3) * t439 + Ifges(5,5) * t478 - t477 * Ifges(5,6) - pkin(8) * t421 - t489 * t424 + t493 * t425;
t399 = -mrSges(4,2) * g(1) - mrSges(4,3) * t446 + Ifges(4,5) * t478 - t477 * Ifges(4,6) - qJ(4) * t415 + t488 * t407 - t487 * t413;
t398 = mrSges(4,1) * g(1) + mrSges(4,3) * t447 + t477 * Ifges(4,5) + Ifges(4,6) * t478 - pkin(3) * t419 + qJ(4) * t505 + t487 * t407 + t488 * t413;
t390 = -mrSges(3,2) * g(1) - mrSges(3,3) * t454 + Ifges(3,5) * t483 - t482 * Ifges(3,6) - pkin(7) * t406 - t490 * t398 + t494 * t399;
t389 = Ifges(3,6) * t483 + t482 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t455 + t490 * t399 + t494 * t398 - pkin(2) * (-m(4) * g(1) + t419) + pkin(7) * t506;
t388 = -mrSges(2,2) * g(1) - mrSges(2,3) * t473 + Ifges(2,5) * qJDD(1) - t497 * Ifges(2,6) - pkin(6) * t397 - t491 * t389 + t495 * t390;
t387 = Ifges(2,6) * qJDD(1) + t497 * Ifges(2,5) + mrSges(2,3) * t472 + t491 * t390 + t495 * t389 - pkin(1) * t419 + pkin(6) * t507 + (-pkin(1) * t513 + mrSges(2,1)) * g(1);
t1 = [(-m(1) - m(2) + t513) * g(1) + t419; -m(1) * g(2) + t510; -m(1) * g(3) + t508; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t498; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t508 + t496 * t387 + t492 * t388; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t510 + t492 * t387 - t496 * t388; t498; t499; t500; t419; t514;];
tauJB = t1;
