% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:40
% EndTime: 2019-12-31 19:12:41
% DurationCPUTime: 1.08s
% Computational Cost: add. (7483->216), mult. (14654->277), div. (0->0), fcn. (9143->8), ass. (0->88)
t361 = sin(qJ(1));
t365 = cos(qJ(1));
t376 = -t365 * g(1) - t361 * g(2);
t372 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t376;
t359 = sin(qJ(4));
t360 = sin(qJ(3));
t363 = cos(qJ(4));
t364 = cos(qJ(3));
t340 = (t364 * t359 + t360 * t363) * qJD(1);
t385 = -pkin(1) - pkin(6);
t366 = qJD(1) ^ 2;
t379 = t361 * g(1) - t365 * g(2);
t371 = -t366 * qJ(2) + qJDD(2) - t379;
t333 = t385 * qJDD(1) + t371;
t324 = t360 * g(3) + t364 * t333;
t382 = qJD(1) * qJD(3);
t380 = t360 * t382;
t345 = t364 * qJDD(1) - t380;
t304 = (-t345 - t380) * pkin(7) + (-t360 * t364 * t366 + qJDD(3)) * pkin(3) + t324;
t325 = -t364 * g(3) + t360 * t333;
t344 = -t360 * qJDD(1) - t364 * t382;
t383 = t364 * qJD(1);
t348 = qJD(3) * pkin(3) - pkin(7) * t383;
t357 = t360 ^ 2;
t305 = -t357 * t366 * pkin(3) + t344 * pkin(7) - qJD(3) * t348 + t325;
t294 = t359 * t304 + t363 * t305;
t341 = (-t360 * t359 + t364 * t363) * qJD(1);
t314 = -t341 * qJD(4) + t363 * t344 - t359 * t345;
t322 = t340 * mrSges(5,1) + t341 * mrSges(5,2);
t355 = qJD(3) + qJD(4);
t331 = t355 * mrSges(5,1) - t341 * mrSges(5,3);
t354 = qJDD(3) + qJDD(4);
t309 = -t344 * pkin(3) + t348 * t383 + (-pkin(7) * t357 + t385) * t366 + t372;
t315 = -t340 * qJD(4) + t359 * t344 + t363 * t345;
t289 = (t340 * t355 - t315) * pkin(8) + (t341 * t355 - t314) * pkin(4) + t309;
t323 = t340 * pkin(4) - t341 * pkin(8);
t353 = t355 ^ 2;
t291 = -t353 * pkin(4) + t354 * pkin(8) - t340 * t323 + t294;
t358 = sin(qJ(5));
t362 = cos(qJ(5));
t287 = t362 * t289 - t358 * t291;
t326 = -t358 * t341 + t362 * t355;
t297 = t326 * qJD(5) + t362 * t315 + t358 * t354;
t327 = t362 * t341 + t358 * t355;
t306 = -t326 * mrSges(6,1) + t327 * mrSges(6,2);
t313 = qJDD(5) - t314;
t336 = qJD(5) + t340;
t316 = -t336 * mrSges(6,2) + t326 * mrSges(6,3);
t284 = m(6) * t287 + t313 * mrSges(6,1) - t297 * mrSges(6,3) - t327 * t306 + t336 * t316;
t288 = t358 * t289 + t362 * t291;
t296 = -t327 * qJD(5) - t358 * t315 + t362 * t354;
t317 = t336 * mrSges(6,1) - t327 * mrSges(6,3);
t285 = m(6) * t288 - t313 * mrSges(6,2) + t296 * mrSges(6,3) + t326 * t306 - t336 * t317;
t377 = -t358 * t284 + t362 * t285;
t272 = m(5) * t294 - t354 * mrSges(5,2) + t314 * mrSges(5,3) - t340 * t322 - t355 * t331 + t377;
t293 = t363 * t304 - t359 * t305;
t330 = -t355 * mrSges(5,2) - t340 * mrSges(5,3);
t290 = -t354 * pkin(4) - t353 * pkin(8) + t341 * t323 - t293;
t370 = -m(6) * t290 + t296 * mrSges(6,1) - t297 * mrSges(6,2) + t326 * t316 - t327 * t317;
t280 = m(5) * t293 + t354 * mrSges(5,1) - t315 * mrSges(5,3) - t341 * t322 + t355 * t330 + t370;
t269 = t359 * t272 + t363 * t280;
t274 = t362 * t284 + t358 * t285;
t384 = qJD(1) * t360;
t378 = t363 * t272 - t359 * t280;
t343 = (t360 * mrSges(4,1) + t364 * mrSges(4,2)) * qJD(1);
t346 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t384;
t347 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t383;
t375 = t364 * (m(4) * t324 + qJDD(3) * mrSges(4,1) - t345 * mrSges(4,3) + qJD(3) * t346 - t343 * t383 + t269) + t360 * (m(4) * t325 - qJDD(3) * mrSges(4,2) + t344 * mrSges(4,3) - qJD(3) * t347 - t343 * t384 + t378);
t369 = m(5) * t309 - t314 * mrSges(5,1) + t315 * mrSges(5,2) + t340 * t330 + t341 * t331 + t274;
t298 = Ifges(6,5) * t327 + Ifges(6,6) * t326 + Ifges(6,3) * t336;
t300 = Ifges(6,1) * t327 + Ifges(6,4) * t326 + Ifges(6,5) * t336;
t277 = -mrSges(6,1) * t290 + mrSges(6,3) * t288 + Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t313 - t327 * t298 + t336 * t300;
t299 = Ifges(6,4) * t327 + Ifges(6,2) * t326 + Ifges(6,6) * t336;
t278 = mrSges(6,2) * t290 - mrSges(6,3) * t287 + Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t313 + t326 * t298 - t336 * t299;
t319 = Ifges(5,4) * t341 - Ifges(5,2) * t340 + Ifges(5,6) * t355;
t320 = Ifges(5,1) * t341 - Ifges(5,4) * t340 + Ifges(5,5) * t355;
t368 = mrSges(5,1) * t293 - mrSges(5,2) * t294 + Ifges(5,5) * t315 + Ifges(5,6) * t314 + Ifges(5,3) * t354 + pkin(4) * t370 + pkin(8) * t377 + t362 * t277 + t358 * t278 + t341 * t319 + t340 * t320;
t367 = mrSges(6,1) * t287 - mrSges(6,2) * t288 + Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t313 + t327 * t299 - t326 * t300;
t339 = (Ifges(4,5) * qJD(3)) + (t364 * Ifges(4,1) - t360 * Ifges(4,4)) * qJD(1);
t338 = (Ifges(4,6) * qJD(3)) + (t364 * Ifges(4,4) - t360 * Ifges(4,2)) * qJD(1);
t335 = -qJDD(1) * pkin(1) + t371;
t334 = t366 * pkin(1) - t372;
t332 = t385 * t366 + t372;
t318 = Ifges(5,5) * t341 - Ifges(5,6) * t340 + Ifges(5,3) * t355;
t266 = -mrSges(5,1) * t309 + mrSges(5,3) * t294 + Ifges(5,4) * t315 + Ifges(5,2) * t314 + Ifges(5,6) * t354 - pkin(4) * t274 - t341 * t318 + t355 * t320 - t367;
t265 = mrSges(5,2) * t309 - mrSges(5,3) * t293 + Ifges(5,1) * t315 + Ifges(5,4) * t314 + Ifges(5,5) * t354 - pkin(8) * t274 - t358 * t277 + t362 * t278 - t340 * t318 - t355 * t319;
t264 = m(3) * t335 + qJDD(1) * mrSges(3,2) - (t366 * mrSges(3,3)) + t375;
t1 = [mrSges(2,1) * t379 - mrSges(2,2) * t376 + mrSges(3,2) * t335 - mrSges(3,3) * t334 + t364 * (mrSges(4,2) * t332 - mrSges(4,3) * t324 + Ifges(4,1) * t345 + Ifges(4,4) * t344 + Ifges(4,5) * qJDD(3) - pkin(7) * t269 - qJD(3) * t338 + t363 * t265 - t359 * t266) - t360 * (-mrSges(4,1) * t332 + mrSges(4,3) * t325 + Ifges(4,4) * t345 + Ifges(4,2) * t344 + Ifges(4,6) * qJDD(3) - pkin(3) * t369 + pkin(7) * t378 + qJD(3) * t339 + t359 * t265 + t363 * t266) - pkin(6) * t375 - pkin(1) * t264 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t334 + m(4) * t332 - t344 * mrSges(4,1) + t366 * mrSges(3,2) + t345 * mrSges(4,2) + t369 + qJDD(1) * mrSges(3,3) + (t346 * t360 + t347 * t364) * qJD(1)) * qJ(2); t264; (t364 * t338 + t360 * t339) * qJD(1) + t368 + Ifges(4,3) * qJDD(3) + Ifges(4,5) * t345 + Ifges(4,6) * t344 + mrSges(4,1) * t324 - mrSges(4,2) * t325 + pkin(3) * t269; t368; t367;];
tauJ = t1;
