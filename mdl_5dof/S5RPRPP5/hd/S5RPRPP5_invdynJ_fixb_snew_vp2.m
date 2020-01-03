% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:12
% EndTime: 2019-12-31 18:16:14
% DurationCPUTime: 0.81s
% Computational Cost: add. (1122->173), mult. (2161->200), div. (0->0), fcn. (821->4), ass. (0->68)
t357 = sin(qJ(3));
t359 = cos(qJ(3));
t379 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t400 = t359 * (Ifges(4,1) + Ifges(5,1) + Ifges(6,1)) - t357 * t379;
t399 = (Ifges(4,2) + Ifges(5,3) + Ifges(6,2)) * t357 - t359 * t379;
t398 = -2 * qJD(1);
t378 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t377 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t393 = t359 * (t399 * qJD(1) - t377 * qJD(3));
t362 = qJD(1) ^ 2;
t358 = sin(qJ(1));
t360 = cos(qJ(1));
t384 = -t360 * g(1) - t358 * g(2);
t306 = (t362 * pkin(1)) - qJDD(1) * qJ(2) + (qJD(2) * t398) - t384;
t392 = 2 * qJD(4);
t391 = -pkin(3) - pkin(4);
t390 = pkin(4) * t362;
t389 = -mrSges(6,2) - mrSges(5,3);
t388 = -mrSges(4,3) - mrSges(5,2);
t387 = qJ(4) * t357;
t381 = qJD(1) * qJD(3);
t331 = t359 * qJDD(1) - t357 * t381;
t386 = t331 * qJ(4);
t382 = t359 * qJD(1);
t339 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t382;
t341 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t382;
t385 = -t339 - t341;
t383 = t357 * qJD(1);
t375 = t400 * qJD(1) + t378 * qJD(3);
t374 = t358 * g(1) - t360 * g(2);
t367 = -t362 * qJ(2) + qJDD(2) - t374;
t305 = (-pkin(1) - pkin(6)) * qJDD(1) + t367;
t302 = -t359 * g(3) + t357 * t305;
t327 = (t357 * mrSges(5,1) - t359 * mrSges(5,3)) * qJD(1);
t373 = qJD(1) * (-t327 - (t357 * mrSges(4,1) + t359 * mrSges(4,2)) * qJD(1));
t330 = t357 * qJDD(1) + t359 * t381;
t338 = -qJD(3) * pkin(4) - qJ(5) * t382;
t355 = t357 ^ 2;
t292 = t386 + qJDD(5) + (-qJ(5) * t355 + pkin(6)) * t362 + t391 * t330 + (-qJD(3) * t387 + (-pkin(3) * qJD(3) + t338 + t392) * t359) * qJD(1) + t306;
t336 = qJD(3) * mrSges(6,2) + mrSges(6,3) * t383;
t371 = m(6) * t292 - t330 * mrSges(6,1) - t336 * t383;
t301 = t357 * g(3) + t359 * t305;
t326 = (t357 * pkin(3) - t359 * qJ(4)) * qJD(1);
t361 = qJD(3) ^ 2;
t370 = -t361 * qJ(4) + t326 * t382 + qJDD(4);
t365 = -t361 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t392 + t302;
t294 = -t355 * t390 + t330 * qJ(5) + qJD(3) * t338 + ((2 * qJD(5)) - t326) * t383 + t365;
t328 = (-t357 * mrSges(6,1) + t359 * mrSges(6,2)) * qJD(1);
t369 = m(6) * t294 + qJDD(3) * mrSges(6,2) + t330 * mrSges(6,3) + qJD(3) * t339 + t328 * t383;
t337 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t383;
t340 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t382;
t295 = -t331 * qJ(5) + ((qJD(5) * t398) - t305) * t359 + t391 * qJDD(3) + (-qJ(5) * t381 + t359 * t390 - g(3)) * t357 + t370;
t291 = m(6) * t295 - qJDD(3) * mrSges(6,1) - t331 * mrSges(6,3) - qJD(3) * t336 - t328 * t382;
t300 = -qJDD(3) * pkin(3) - t301 + t370;
t342 = -mrSges(5,2) * t383 + qJD(3) * mrSges(5,3);
t363 = -m(5) * t300 + qJDD(3) * mrSges(5,1) + qJD(3) * t342 - t291;
t299 = -t326 * t383 + t365;
t364 = m(5) * t299 + qJDD(3) * mrSges(5,3) + qJD(3) * t341 + t369;
t368 = t357 * (m(4) * t302 - qJDD(3) * mrSges(4,2) - qJD(3) * t340 + t388 * t330 + t357 * t373 + t364) + t359 * (m(4) * t301 + qJDD(3) * mrSges(4,1) + qJD(3) * t337 + t388 * t331 + t359 * t373 + t363);
t304 = -(t362 * pkin(6)) - t306;
t297 = t330 * pkin(3) - t386 + (-0.2e1 * qJD(4) * t359 + (pkin(3) * t359 + t387) * qJD(3)) * qJD(1) + t304;
t366 = m(5) * t297 + t330 * mrSges(5,1) + t342 * t383 - t371;
t307 = -qJDD(1) * pkin(1) + t367;
t290 = t331 * mrSges(6,2) + t339 * t382 + t371;
t289 = t331 * mrSges(5,2) + t327 * t382 - t363;
t288 = t389 * t331 + t385 * t382 + t366;
t285 = m(3) * t307 + qJDD(1) * mrSges(3,2) - (t362 * mrSges(3,3)) + t368;
t1 = [mrSges(2,1) * t374 - mrSges(2,2) * t384 + mrSges(3,2) * t307 - mrSges(3,3) * t306 + t359 * (mrSges(4,2) * t304 + mrSges(5,2) * t300 + mrSges(6,2) * t292 - mrSges(4,3) * t301 - mrSges(5,3) * t297 - mrSges(6,3) * t295 - qJ(4) * t288 - qJ(5) * t291) - t357 * (-mrSges(4,1) * t304 - mrSges(5,1) * t297 + mrSges(6,1) * t292 + mrSges(5,2) * t299 + mrSges(4,3) * t302 - mrSges(6,3) * t294 - pkin(3) * t288 + pkin(4) * t290 - qJ(5) * t369) - pkin(6) * t368 - pkin(1) * t285 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t399 * t330 + (-t357 * t377 + t359 * t378) * qJDD(3) + (-t357 * t375 + t393) * qJD(3) + t400 * t331 + (-m(3) * t306 + m(4) * t304 + t362 * mrSges(3,2) + t366 + qJDD(1) * mrSges(3,3) + t330 * mrSges(4,1) + (mrSges(4,2) + t389) * t331 + (t337 * t357 + (t340 + t385) * t359) * qJD(1)) * qJ(2); t285; mrSges(4,1) * t301 - mrSges(4,2) * t302 - mrSges(5,1) * t300 + mrSges(5,3) * t299 - mrSges(6,1) * t295 + mrSges(6,2) * t294 - pkin(4) * t291 - pkin(3) * t289 + qJ(4) * t364 + t378 * t331 + (-qJ(4) * mrSges(5,2) - t377) * t330 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * qJDD(3) + (-t393 + (-qJ(4) * t327 + t375) * t357) * qJD(1); t289; t290;];
tauJ = t1;
