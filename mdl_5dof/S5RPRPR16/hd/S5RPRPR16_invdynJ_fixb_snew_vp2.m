% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR16_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:51
% EndTime: 2019-12-31 18:38:52
% DurationCPUTime: 1.04s
% Computational Cost: add. (2405->198), mult. (4698->236), div. (0->0), fcn. (2201->6), ass. (0->88)
t357 = sin(qJ(3));
t360 = cos(qJ(3));
t392 = Ifges(4,4) + Ifges(5,6);
t403 = t360 * (Ifges(4,1) + Ifges(5,2)) - t357 * t392;
t402 = t360 * t392 + t357 * (-Ifges(4,2) - Ifges(5,3));
t391 = Ifges(4,5) - Ifges(5,4);
t390 = Ifges(4,6) - Ifges(5,5);
t399 = (t402 * qJD(1) + t390 * qJD(3)) * t360;
t358 = sin(qJ(1));
t361 = cos(qJ(1));
t376 = -g(1) * t361 - g(2) * t358;
t372 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t376;
t398 = -2 * qJD(4);
t397 = -pkin(1) - pkin(6);
t396 = pkin(3) + pkin(7);
t363 = qJD(1) ^ 2;
t395 = pkin(7) * t363;
t394 = g(3) * t357;
t393 = mrSges(4,1) - mrSges(5,2);
t379 = g(1) * t358 - t361 * g(2);
t371 = -qJ(2) * t363 + qJDD(2) - t379;
t319 = t397 * qJDD(1) + t371;
t389 = t319 * t360;
t387 = t403 * qJD(1) + t391 * qJD(3);
t385 = qJD(1) * t357;
t343 = mrSges(5,1) * t385 - qJD(3) * mrSges(5,3);
t386 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t385 - t343;
t384 = qJD(1) * t360;
t383 = qJD(1) * qJD(3);
t381 = t397 * t363;
t380 = t360 * t383;
t350 = t357 * t383;
t314 = -g(3) * t360 + t357 * t319;
t338 = qJDD(1) * t357 + t380;
t345 = pkin(4) * t384 - qJD(3) * pkin(7);
t355 = t357 ^ 2;
t339 = qJDD(1) * t360 - t350;
t365 = pkin(3) * t380 + t384 * t398 + t372 + (-t339 + t350) * qJ(4);
t294 = -t345 * t384 + t396 * t338 + (-pkin(4) * t355 + t397) * t363 + t365;
t335 = (pkin(3) * t357 - qJ(4) * t360) * qJD(1);
t362 = qJD(3) ^ 2;
t370 = -qJ(4) * t362 + t335 * t384 + qJDD(4) - t389;
t297 = pkin(4) * t339 - t396 * qJDD(3) + (pkin(4) * t383 + t360 * t395 - g(3)) * t357 + t370;
t356 = sin(qJ(5));
t359 = cos(qJ(5));
t292 = -t294 * t356 + t297 * t359;
t333 = -qJD(3) * t356 + t359 * t385;
t311 = qJD(5) * t333 + qJDD(3) * t359 + t338 * t356;
t334 = qJD(3) * t359 + t356 * t385;
t312 = -mrSges(6,1) * t333 + mrSges(6,2) * t334;
t348 = qJD(5) + t384;
t315 = -mrSges(6,2) * t348 + mrSges(6,3) * t333;
t332 = qJDD(5) + t339;
t289 = m(6) * t292 + mrSges(6,1) * t332 - mrSges(6,3) * t311 - t312 * t334 + t315 * t348;
t293 = t294 * t359 + t297 * t356;
t310 = -qJD(5) * t334 - qJDD(3) * t356 + t338 * t359;
t316 = mrSges(6,1) * t348 - mrSges(6,3) * t334;
t290 = m(6) * t293 - mrSges(6,2) * t332 + mrSges(6,3) * t310 + t312 * t333 - t316 * t348;
t378 = -t289 * t356 + t359 * t290;
t336 = (-mrSges(5,2) * t357 - mrSges(5,3) * t360) * qJD(1);
t377 = qJD(1) * (-t336 - (mrSges(4,1) * t357 + mrSges(4,2) * t360) * qJD(1));
t313 = t389 + t394;
t342 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t384;
t366 = -pkin(3) * t362 + qJDD(3) * qJ(4) - t335 * t385 + t314;
t300 = qJD(3) * t398 - t366;
t344 = mrSges(5,1) * t384 + qJD(3) * mrSges(5,2);
t296 = -t355 * t395 - pkin(4) * t338 + ((2 * qJD(4)) + t345) * qJD(3) + t366;
t369 = -m(6) * t296 + mrSges(6,1) * t310 - t311 * mrSges(6,2) + t315 * t333 - t334 * t316;
t364 = -m(5) * t300 + qJDD(3) * mrSges(5,3) + qJD(3) * t344 - t369;
t284 = t289 * t359 + t290 * t356;
t301 = -qJDD(3) * pkin(3) + t370 - t394;
t368 = -m(5) * t301 - t339 * mrSges(5,1) - t284;
t375 = (m(4) * t313 - mrSges(4,3) * t339 + t386 * qJD(3) + t393 * qJDD(3) + t360 * t377 + t368) * t360 + (m(4) * t314 - qJDD(3) * mrSges(4,2) - qJD(3) * t342 + (-mrSges(4,3) - mrSges(5,1)) * t338 + t357 * t377 + t364) * t357;
t299 = pkin(3) * t338 + t365 + t381;
t373 = m(5) * t299 - t339 * mrSges(5,3) + t378;
t304 = Ifges(6,4) * t334 + Ifges(6,2) * t333 + Ifges(6,6) * t348;
t305 = Ifges(6,1) * t334 + Ifges(6,4) * t333 + Ifges(6,5) * t348;
t367 = mrSges(6,1) * t292 - mrSges(6,2) * t293 + Ifges(6,5) * t311 + Ifges(6,6) * t310 + Ifges(6,3) * t332 + t334 * t304 - t333 * t305;
t321 = -qJDD(1) * pkin(1) + t371;
t320 = pkin(1) * t363 - t372;
t318 = t381 + t372;
t303 = Ifges(6,5) * t334 + Ifges(6,6) * t333 + Ifges(6,3) * t348;
t286 = mrSges(6,2) * t296 - mrSges(6,3) * t292 + Ifges(6,1) * t311 + Ifges(6,4) * t310 + Ifges(6,5) * t332 + t303 * t333 - t304 * t348;
t285 = -mrSges(6,1) * t296 + mrSges(6,3) * t293 + Ifges(6,4) * t311 + Ifges(6,2) * t310 + Ifges(6,6) * t332 - t303 * t334 + t305 * t348;
t283 = qJDD(3) * mrSges(5,2) + qJD(3) * t343 + t336 * t384 - t368;
t282 = -mrSges(5,2) * t338 + (-t343 * t357 - t344 * t360) * qJD(1) + t373;
t280 = m(3) * t321 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t363 + t375;
t1 = [mrSges(2,1) * t379 - mrSges(2,2) * t376 + mrSges(3,2) * t321 - mrSges(3,3) * t320 + t360 * (mrSges(5,1) * t301 + mrSges(4,2) * t318 - mrSges(4,3) * t313 - mrSges(5,3) * t299 + pkin(4) * t284 - qJ(4) * t282 + t367) - t357 * (-mrSges(4,1) * t318 - mrSges(5,1) * t300 + mrSges(5,2) * t299 + mrSges(4,3) * t314 - pkin(3) * t282 - pkin(4) * t369 - pkin(7) * t378 - t359 * t285 - t356 * t286) - pkin(6) * t375 - pkin(1) * t280 + t403 * t339 + (-t357 * t390 + t360 * t391) * qJDD(3) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t357 * t387 - t399) * qJD(3) - t402 * t338 + (-m(3) * t320 + m(4) * t318 + mrSges(3,2) * t363 + t373 + mrSges(4,2) * t339 + qJDD(1) * mrSges(3,3) + t393 * t338 + ((t342 - t344) * t360 + t386 * t357) * qJD(1)) * qJ(2); t280; mrSges(4,1) * t313 - mrSges(4,2) * t314 + mrSges(5,2) * t301 - mrSges(5,3) * t300 + t359 * t286 - t356 * t285 - pkin(7) * t284 - pkin(3) * t283 + qJ(4) * t364 + t391 * t339 + (-qJ(4) * mrSges(5,1) - t390) * t338 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + (t399 + (-qJ(4) * t336 + t387) * t357) * qJD(1); t283; t367;];
tauJ = t1;
