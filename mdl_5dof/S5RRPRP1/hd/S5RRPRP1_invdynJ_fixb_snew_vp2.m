% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:43
% EndTime: 2022-01-20 10:19:44
% DurationCPUTime: 0.58s
% Computational Cost: add. (3017->146), mult. (3889->177), div. (0->0), fcn. (1891->8), ass. (0->71)
t345 = Ifges(5,1) + Ifges(6,1);
t340 = Ifges(5,4) + Ifges(6,4);
t339 = Ifges(5,5) + Ifges(6,5);
t344 = Ifges(5,2) + Ifges(6,2);
t338 = Ifges(5,6) + Ifges(6,6);
t343 = Ifges(5,3) + Ifges(6,3);
t307 = qJD(1) + qJD(2);
t305 = t307 ^ 2;
t342 = pkin(4) * t305;
t341 = -mrSges(5,2) - mrSges(6,2);
t313 = sin(qJ(4));
t337 = t307 * t313;
t316 = cos(qJ(4));
t336 = t307 * t316;
t315 = sin(qJ(1));
t318 = cos(qJ(1));
t325 = t315 * g(1) - t318 * g(2);
t293 = qJDD(1) * pkin(1) + t325;
t322 = -t318 * g(1) - t315 * g(2);
t294 = -qJD(1) ^ 2 * pkin(1) + t322;
t314 = sin(qJ(2));
t317 = cos(qJ(2));
t271 = t317 * t293 - t314 * t294;
t306 = qJDD(1) + qJDD(2);
t268 = t306 * pkin(2) + t271;
t272 = t314 * t293 + t317 * t294;
t269 = -t305 * pkin(2) + t272;
t311 = sin(pkin(8));
t312 = cos(pkin(8));
t264 = t311 * t268 + t312 * t269;
t261 = -t305 * pkin(3) + t306 * pkin(7) + t264;
t310 = -g(3) + qJDD(3);
t301 = t316 * t310;
t257 = -t313 * t261 + t301;
t285 = (-mrSges(6,1) * t316 + mrSges(6,2) * t313) * t307;
t286 = (-mrSges(5,1) * t316 + mrSges(5,2) * t313) * t307;
t330 = qJD(4) * t307;
t326 = t316 * t330;
t287 = t313 * t306 + t326;
t299 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t336;
t329 = qJD(5) * t307;
t254 = qJDD(4) * pkin(4) + t301 + (-t287 + t326) * qJ(5) + (t316 * t342 - t261 - 0.2e1 * t329) * t313;
t298 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t336;
t328 = m(6) * t254 + qJDD(4) * mrSges(6,1) + qJD(4) * t298;
t248 = m(5) * t257 + qJDD(4) * mrSges(5,1) + qJD(4) * t299 + (-t285 - t286) * t337 + (-mrSges(5,3) - mrSges(6,3)) * t287 + t328;
t258 = t316 * t261 + t313 * t310;
t288 = t316 * t306 - t313 * t330;
t295 = qJD(4) * pkin(4) - qJ(5) * t337;
t309 = t316 ^ 2;
t255 = t288 * qJ(5) - qJD(4) * t295 - t309 * t342 + 0.2e1 * t316 * t329 + t258;
t327 = m(6) * t255 + t288 * mrSges(6,3) + t285 * t336;
t296 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t337;
t331 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t337 - t296;
t249 = m(5) * t258 + t288 * mrSges(5,3) + t331 * qJD(4) + t341 * qJDD(4) + t286 * t336 + t327;
t324 = -t313 * t248 + t316 * t249;
t242 = m(4) * t264 - t305 * mrSges(4,1) - t306 * mrSges(4,2) + t324;
t263 = t312 * t268 - t311 * t269;
t321 = -t306 * pkin(3) - t263;
t260 = -t305 * pkin(7) + t321;
t256 = t295 * t337 - t288 * pkin(4) + qJDD(5) + (-qJ(5) * t309 - pkin(7)) * t305 + t321;
t323 = -m(6) * t256 + t288 * mrSges(6,1) + t298 * t336;
t319 = -m(5) * t260 + t288 * mrSges(5,1) + t341 * t287 + t299 * t336 + t331 * t337 + t323;
t245 = m(4) * t263 + t306 * mrSges(4,1) - t305 * mrSges(4,2) + t319;
t335 = t311 * t242 + t312 * t245;
t334 = (t313 * t339 + t338 * t316) * t307 + t343 * qJD(4);
t333 = (t313 * t340 + t344 * t316) * t307 + t338 * qJD(4);
t332 = (-t345 * t313 - t316 * t340) * t307 - t339 * qJD(4);
t250 = t287 * mrSges(6,2) + t296 * t337 - t323;
t251 = -t287 * mrSges(6,3) - t285 * t337 + t328;
t320 = -mrSges(3,2) * t272 - mrSges(4,2) * t264 + pkin(2) * t335 + t316 * (-mrSges(5,1) * t260 + mrSges(5,3) * t258 - mrSges(6,1) * t256 + mrSges(6,3) * t255 - pkin(4) * t250 + qJ(5) * t327 - t334 * t337 + t344 * t288 + t340 * t287 + (-qJ(5) * mrSges(6,2) + t338) * qJDD(4) + (-qJ(5) * t296 - t332) * qJD(4)) + t313 * (mrSges(5,2) * t260 + mrSges(6,2) * t256 - mrSges(5,3) * t257 - mrSges(6,3) * t254 - qJ(5) * t251 - t333 * qJD(4) + t339 * qJDD(4) + t345 * t287 + t340 * t288 + t334 * t336) + pkin(7) * t324 + pkin(3) * t319 + mrSges(4,1) * t263 + mrSges(3,1) * t271 + (Ifges(4,3) + Ifges(3,3)) * t306;
t1 = [Ifges(2,3) * qJDD(1) + pkin(1) * (t314 * (m(3) * t272 - t305 * mrSges(3,1) - t306 * mrSges(3,2) + t312 * t242 - t311 * t245) + t317 * (m(3) * t271 + t306 * mrSges(3,1) - t305 * mrSges(3,2) + t335)) + mrSges(2,1) * t325 - mrSges(2,2) * t322 + t320; t320; m(4) * t310 + t316 * t248 + t313 * t249; mrSges(5,1) * t257 + mrSges(6,1) * t254 - mrSges(5,2) * t258 - mrSges(6,2) * t255 + pkin(4) * t251 + t338 * t288 + t339 * t287 + t343 * qJDD(4) + (t333 * t313 + t332 * t316) * t307; t250;];
tauJ = t1;
