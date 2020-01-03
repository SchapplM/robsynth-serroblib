% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:50
% EndTime: 2019-12-31 21:50:52
% DurationCPUTime: 1.16s
% Computational Cost: add. (8850->199), mult. (11141->246), div. (0->0), fcn. (6549->8), ass. (0->85)
t400 = Ifges(5,1) + Ifges(6,1);
t390 = Ifges(5,4) - Ifges(6,5);
t398 = Ifges(6,4) + Ifges(5,5);
t399 = Ifges(5,2) + Ifges(6,3);
t397 = Ifges(5,6) - Ifges(6,6);
t396 = Ifges(5,3) + Ifges(6,2);
t366 = sin(qJ(4));
t364 = qJD(1) + qJD(2);
t370 = cos(qJ(3));
t388 = t364 * t370;
t367 = sin(qJ(3));
t389 = t364 * t367;
t393 = cos(qJ(4));
t335 = t366 * t389 - t393 * t388;
t336 = (t366 * t370 + t393 * t367) * t364;
t363 = qJD(3) + qJD(4);
t395 = t399 * t335 - t390 * t336 - t397 * t363;
t394 = -t390 * t335 + t400 * t336 + t398 * t363;
t360 = t364 ^ 2;
t392 = pkin(3) * t360;
t391 = -mrSges(5,3) - mrSges(6,2);
t369 = sin(qJ(1));
t372 = cos(qJ(1));
t382 = t369 * g(1) - t372 * g(2);
t347 = qJDD(1) * pkin(1) + t382;
t378 = -t372 * g(1) - t369 * g(2);
t348 = -qJD(1) ^ 2 * pkin(1) + t378;
t368 = sin(qJ(2));
t371 = cos(qJ(2));
t325 = t368 * t347 + t371 * t348;
t362 = qJDD(1) + qJDD(2);
t322 = -t360 * pkin(2) + t362 * pkin(7) + t325;
t387 = t367 * t322;
t384 = qJD(3) * t364;
t342 = t367 * t362 + t370 * t384;
t292 = qJDD(3) * pkin(3) - t342 * pkin(8) - t387 + (pkin(8) * t384 + t367 * t392 - g(3)) * t370;
t307 = -t367 * g(3) + t370 * t322;
t343 = t370 * t362 - t367 * t384;
t351 = qJD(3) * pkin(3) - pkin(8) * t389;
t365 = t370 ^ 2;
t293 = t343 * pkin(8) - qJD(3) * t351 - t365 * t392 + t307;
t289 = t366 * t292 + t393 * t293;
t304 = t336 * qJD(4) + t366 * t342 - t393 * t343;
t329 = t363 * mrSges(5,1) - t336 * mrSges(5,3);
t361 = qJDD(3) + qJDD(4);
t318 = t335 * pkin(4) - t336 * qJ(5);
t359 = t363 ^ 2;
t285 = -t359 * pkin(4) + t361 * qJ(5) + 0.2e1 * qJD(5) * t363 - t335 * t318 + t289;
t330 = -t363 * mrSges(6,1) + t336 * mrSges(6,2);
t383 = m(6) * t285 + t361 * mrSges(6,3) + t363 * t330;
t319 = t335 * mrSges(6,1) - t336 * mrSges(6,3);
t385 = -t335 * mrSges(5,1) - t336 * mrSges(5,2) - t319;
t274 = m(5) * t289 - t361 * mrSges(5,2) + t391 * t304 - t363 * t329 + t385 * t335 + t383;
t288 = t393 * t292 - t366 * t293;
t305 = -t335 * qJD(4) + t393 * t342 + t366 * t343;
t328 = -t363 * mrSges(5,2) - t335 * mrSges(5,3);
t286 = -t361 * pkin(4) - t359 * qJ(5) + t336 * t318 + qJDD(5) - t288;
t331 = -t335 * mrSges(6,2) + t363 * mrSges(6,3);
t379 = -m(6) * t286 + t361 * mrSges(6,1) + t363 * t331;
t276 = m(5) * t288 + t361 * mrSges(5,1) + t391 * t305 + t363 * t328 + t385 * t336 + t379;
t270 = t366 * t274 + t393 * t276;
t386 = t397 * t335 - t398 * t336 - t396 * t363;
t306 = -t370 * g(3) - t387;
t341 = (-mrSges(4,1) * t370 + mrSges(4,2) * t367) * t364;
t349 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t389;
t350 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t388;
t380 = t393 * t274 - t366 * t276;
t381 = -t367 * (m(4) * t306 + qJDD(3) * mrSges(4,1) - t342 * mrSges(4,3) + qJD(3) * t350 - t341 * t389 + t270) + t370 * (m(4) * t307 - qJDD(3) * mrSges(4,2) + t343 * mrSges(4,3) - qJD(3) * t349 + t341 * t388 + t380);
t324 = t371 * t347 - t368 * t348;
t377 = -t362 * pkin(2) - t324;
t294 = -t343 * pkin(3) + t351 * t389 + (-pkin(8) * t365 - pkin(7)) * t360 + t377;
t282 = -0.2e1 * qJD(5) * t336 + (t335 * t363 - t305) * qJ(5) + (t336 * t363 + t304) * pkin(4) + t294;
t277 = m(6) * t282 + t304 * mrSges(6,1) - t305 * mrSges(6,3) - t336 * t330 + t335 * t331;
t266 = -mrSges(5,1) * t294 - mrSges(6,1) * t282 + mrSges(6,2) * t285 + mrSges(5,3) * t289 - pkin(4) * t277 - t399 * t304 + t390 * t305 + t386 * t336 + t397 * t361 + t394 * t363;
t267 = mrSges(5,2) * t294 + mrSges(6,2) * t286 - mrSges(5,3) * t288 - mrSges(6,3) * t282 - qJ(5) * t277 - t390 * t304 + t400 * t305 + t386 * t335 + t398 * t361 + t395 * t363;
t321 = -t360 * pkin(7) + t377;
t332 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t367 + Ifges(4,6) * t370) * t364;
t333 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t367 + Ifges(4,2) * t370) * t364;
t334 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t367 + Ifges(4,4) * t370) * t364;
t375 = m(5) * t294 + t304 * mrSges(5,1) + t305 * mrSges(5,2) + t335 * t328 + t336 * t329 + t277;
t373 = -m(4) * t321 + t343 * mrSges(4,1) - t342 * mrSges(4,2) - t349 * t389 + t350 * t388 - t375;
t376 = -mrSges(3,2) * t325 + t370 * (-mrSges(4,1) * t321 + mrSges(4,3) * t307 + Ifges(4,4) * t342 + Ifges(4,2) * t343 + Ifges(4,6) * qJDD(3) - pkin(3) * t375 + pkin(8) * t380 + qJD(3) * t334 + t393 * t266 + t366 * t267 - t332 * t389) + t367 * (mrSges(4,2) * t321 - mrSges(4,3) * t306 + Ifges(4,1) * t342 + Ifges(4,4) * t343 + Ifges(4,5) * qJDD(3) - pkin(8) * t270 - qJD(3) * t333 - t366 * t266 + t393 * t267 + t332 * t388) + pkin(7) * t381 + pkin(2) * t373 + mrSges(3,1) * t324 + Ifges(3,3) * t362;
t280 = t305 * mrSges(6,2) + t336 * t319 - t379;
t374 = mrSges(5,1) * t288 - mrSges(6,1) * t286 - mrSges(5,2) * t289 + mrSges(6,3) * t285 - pkin(4) * t280 + qJ(5) * t383 + t396 * t361 - t395 * t336 + (-qJ(5) * t319 + t394) * t335 + t398 * t305 + (-qJ(5) * mrSges(6,2) - t397) * t304;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t382 - mrSges(2,2) * t378 + pkin(1) * (t368 * (m(3) * t325 - t360 * mrSges(3,1) - t362 * mrSges(3,2) + t381) + t371 * (m(3) * t324 + t362 * mrSges(3,1) - t360 * mrSges(3,2) + t373)) + t376; t376; t374 + Ifges(4,3) * qJDD(3) + (t367 * t333 - t370 * t334) * t364 + Ifges(4,5) * t342 + Ifges(4,6) * t343 + mrSges(4,1) * t306 - mrSges(4,2) * t307 + pkin(3) * t270; t374; t280;];
tauJ = t1;
