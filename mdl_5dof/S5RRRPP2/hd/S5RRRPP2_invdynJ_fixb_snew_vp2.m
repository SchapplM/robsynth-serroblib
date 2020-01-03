% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:35
% EndTime: 2019-12-31 20:51:36
% DurationCPUTime: 0.67s
% Computational Cost: add. (2889->173), mult. (3608->212), div. (0->0), fcn. (1553->6), ass. (0->70)
t408 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t399 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t398 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t407 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t397 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t406 = Ifges(4,3) + Ifges(5,2) + Ifges(6,3);
t405 = 2 * qJD(4);
t368 = qJD(1) + qJD(2);
t366 = t368 ^ 2;
t404 = t366 * pkin(7);
t403 = mrSges(4,3) + mrSges(5,2);
t375 = sin(qJ(3));
t402 = t368 * t375;
t378 = cos(qJ(3));
t401 = t368 * t378;
t377 = sin(qJ(1));
t380 = cos(qJ(1));
t391 = t377 * g(1) - t380 * g(2);
t352 = qJDD(1) * pkin(1) + t391;
t387 = -t380 * g(1) - t377 * g(2);
t353 = -qJD(1) ^ 2 * pkin(1) + t387;
t376 = sin(qJ(2));
t379 = cos(qJ(2));
t313 = t376 * t352 + t379 * t353;
t367 = qJDD(1) + qJDD(2);
t310 = -t366 * pkin(2) + t367 * pkin(7) + t313;
t305 = -t378 * g(3) - t375 * t310;
t312 = t379 * t352 - t376 * t353;
t400 = qJD(3) * t378;
t396 = -0.2e1 * qJD(5) * t368;
t395 = (-t375 * t398 - t378 * t397) * t368 - t406 * qJD(3);
t394 = (-t375 * t399 - t407 * t378) * t368 - t397 * qJD(3);
t393 = (t408 * t375 + t378 * t399) * t368 + t398 * qJD(3);
t392 = t368 * t400;
t306 = -t375 * g(3) + t378 * t310;
t339 = (-mrSges(5,1) * t378 - mrSges(5,3) * t375) * t368;
t340 = (mrSges(6,1) * t378 + mrSges(6,2) * t375) * t368;
t341 = (-mrSges(4,1) * t378 + mrSges(4,2) * t375) * t368;
t342 = t375 * t367 + t392;
t343 = -qJD(3) * t402 + t378 * t367;
t356 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t402;
t359 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t401;
t338 = (-pkin(3) * t378 - qJ(4) * t375) * t368;
t381 = qJD(3) ^ 2;
t304 = -qJDD(3) * pkin(3) - t381 * qJ(4) + t338 * t402 + qJDD(4) - t305;
t300 = t375 * t396 + (-t342 + t392) * qJ(5) + (-t366 * t375 * t378 - qJDD(3)) * pkin(4) + t304;
t358 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t401;
t295 = m(6) * t300 - qJDD(3) * mrSges(6,1) - t342 * mrSges(6,3) - qJD(3) * t358 - t340 * t402;
t360 = mrSges(5,2) * t401 + qJD(3) * mrSges(5,3);
t383 = -m(5) * t304 + qJDD(3) * mrSges(5,1) + qJD(3) * t360 - t295;
t303 = -t381 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t405 + t338 * t401 + t306;
t357 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t402;
t354 = -qJD(3) * pkin(4) - qJ(5) * t402;
t374 = t378 ^ 2;
t299 = -t374 * t366 * pkin(4) - t343 * qJ(5) + qJD(3) * t354 + t378 * t396 + t303;
t355 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t402;
t389 = m(6) * t299 + qJDD(3) * mrSges(6,2) - t343 * mrSges(6,3) + qJD(3) * t355;
t385 = m(5) * t303 + qJDD(3) * mrSges(5,3) + qJD(3) * t357 + t339 * t401 + t389;
t390 = -t375 * (m(4) * t305 + qJDD(3) * mrSges(4,1) + qJD(3) * t359 + (-t339 - t341) * t402 - t403 * t342 + t383) + t378 * (m(4) * t306 - qJDD(3) * mrSges(4,2) - qJD(3) * t356 + (-t340 + t341) * t401 + t403 * t343 + t385);
t388 = t367 * pkin(2) + t312;
t386 = -t342 * qJ(4) - t388;
t297 = qJDD(5) + (-qJ(5) * t374 + pkin(7)) * t366 + (pkin(3) + pkin(4)) * t343 + (qJ(4) * t400 + (-pkin(3) * qJD(3) + t354 + t405) * t375) * t368 - t386;
t294 = m(6) * t297 + t343 * mrSges(6,1) + t342 * mrSges(6,2) + t355 * t402 + t358 * t401;
t301 = -t343 * pkin(3) - t404 + (-0.2e1 * qJD(4) * t375 + (pkin(3) * t375 - qJ(4) * t378) * qJD(3)) * t368 + t386;
t292 = m(5) * t301 - t343 * mrSges(5,1) - t342 * mrSges(5,3) - t357 * t402 - t360 * t401 - t294;
t309 = -t388 - t404;
t382 = -m(4) * t309 + t343 * mrSges(4,1) - t342 * mrSges(4,2) - t356 * t402 + t359 * t401 - t292;
t384 = -mrSges(3,2) * t313 + t378 * (-mrSges(4,1) * t309 + mrSges(4,3) * t306 - mrSges(5,1) * t301 + mrSges(5,2) * t303 + mrSges(6,1) * t297 - mrSges(6,3) * t299 + pkin(4) * t294 - qJ(5) * t389 - pkin(3) * t292 + (qJ(5) * t340 * t378 + t395 * t375) * t368 + t407 * t343 + t399 * t342 + t397 * qJDD(3) + t393 * qJD(3)) + t375 * (mrSges(4,2) * t309 + mrSges(5,2) * t304 + mrSges(6,2) * t297 - mrSges(4,3) * t305 - mrSges(5,3) * t301 - mrSges(6,3) * t300 - qJ(4) * t292 - qJ(5) * t295 + t394 * qJD(3) + t398 * qJDD(3) + t408 * t342 + t399 * t343 - t395 * t401) + pkin(7) * t390 + pkin(2) * t382 + mrSges(3,1) * t312 + Ifges(3,3) * t367;
t293 = t342 * mrSges(5,2) + t339 * t402 - t383;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t391 - mrSges(2,2) * t387 + pkin(1) * (t376 * (m(3) * t313 - t366 * mrSges(3,1) - t367 * mrSges(3,2) + t390) + t379 * (m(3) * t312 + t367 * mrSges(3,1) - t366 * mrSges(3,2) + t382)) + t384; t384; mrSges(4,1) * t305 - mrSges(4,2) * t306 - mrSges(5,1) * t304 + mrSges(5,3) * t303 - mrSges(6,1) * t300 + mrSges(6,2) * t299 - pkin(4) * t295 - pkin(3) * t293 + qJ(4) * t385 + (qJ(4) * mrSges(5,2) + t397) * t343 + t398 * t342 + t406 * qJDD(3) + (-t394 * t375 + (-qJ(4) * t340 - t393) * t378) * t368; t293; t294;];
tauJ = t1;
