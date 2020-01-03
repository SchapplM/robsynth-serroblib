% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:15
% EndTime: 2020-01-03 11:45:16
% DurationCPUTime: 0.55s
% Computational Cost: add. (2571->147), mult. (3501->178), div. (0->0), fcn. (1699->8), ass. (0->72)
t344 = Ifges(5,1) + Ifges(6,1);
t339 = Ifges(5,4) + Ifges(6,4);
t338 = Ifges(5,5) + Ifges(6,5);
t343 = Ifges(5,2) + Ifges(6,2);
t337 = Ifges(5,6) + Ifges(6,6);
t342 = Ifges(5,3) + Ifges(6,3);
t305 = qJD(1) + qJD(3);
t303 = t305 ^ 2;
t341 = pkin(4) * t303;
t340 = -mrSges(5,2) - mrSges(6,2);
t311 = sin(qJ(4));
t336 = t305 * t311;
t314 = cos(qJ(4));
t335 = t305 * t314;
t313 = sin(qJ(1));
t316 = cos(qJ(1));
t321 = -t316 * g(2) - t313 * g(3);
t292 = qJDD(1) * pkin(1) + t321;
t317 = qJD(1) ^ 2;
t324 = -t313 * g(2) + t316 * g(3);
t293 = -t317 * pkin(1) + t324;
t309 = sin(pkin(8));
t310 = cos(pkin(8));
t270 = t310 * t292 - t309 * t293;
t268 = qJDD(1) * pkin(2) + t270;
t271 = t309 * t292 + t310 * t293;
t269 = -t317 * pkin(2) + t271;
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t264 = t312 * t268 + t315 * t269;
t304 = qJDD(1) + qJDD(3);
t261 = -t303 * pkin(3) + t304 * pkin(7) + t264;
t308 = -g(1) + qJDD(2);
t300 = t314 * t308;
t257 = -t311 * t261 + t300;
t284 = (-mrSges(6,1) * t314 + mrSges(6,2) * t311) * t305;
t285 = (-mrSges(5,1) * t314 + mrSges(5,2) * t311) * t305;
t329 = qJD(4) * t305;
t325 = t314 * t329;
t286 = t311 * t304 + t325;
t298 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t335;
t328 = qJD(5) * t305;
t254 = qJDD(4) * pkin(4) + t300 + (-t286 + t325) * qJ(5) + (t314 * t341 - t261 - 0.2e1 * t328) * t311;
t297 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t335;
t327 = m(6) * t254 + qJDD(4) * mrSges(6,1) + qJD(4) * t297;
t248 = m(5) * t257 + qJDD(4) * mrSges(5,1) + qJD(4) * t298 + (-t284 - t285) * t336 + (-mrSges(5,3) - mrSges(6,3)) * t286 + t327;
t258 = t314 * t261 + t311 * t308;
t287 = t314 * t304 - t311 * t329;
t294 = qJD(4) * pkin(4) - qJ(5) * t336;
t307 = t314 ^ 2;
t255 = t287 * qJ(5) - qJD(4) * t294 - t307 * t341 + 0.2e1 * t314 * t328 + t258;
t326 = m(6) * t255 + t287 * mrSges(6,3) + t284 * t335;
t295 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t336;
t330 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t336 - t295;
t249 = m(5) * t258 + t287 * mrSges(5,3) + t330 * qJD(4) + t340 * qJDD(4) + t285 * t335 + t326;
t323 = -t311 * t248 + t314 * t249;
t242 = m(4) * t264 - t303 * mrSges(4,1) - t304 * mrSges(4,2) + t323;
t263 = t315 * t268 - t312 * t269;
t320 = -t304 * pkin(3) - t263;
t260 = -t303 * pkin(7) + t320;
t256 = t294 * t336 - t287 * pkin(4) + qJDD(5) + (-qJ(5) * t307 - pkin(7)) * t303 + t320;
t322 = -m(6) * t256 + t287 * mrSges(6,1) + t297 * t335;
t318 = -m(5) * t260 + t287 * mrSges(5,1) + t340 * t286 + t298 * t335 + t330 * t336 + t322;
t245 = m(4) * t263 + t304 * mrSges(4,1) - t303 * mrSges(4,2) + t318;
t334 = t312 * t242 + t315 * t245;
t333 = (t311 * t338 + t337 * t314) * t305 + t342 * qJD(4);
t332 = (t311 * t339 + t343 * t314) * t305 + t337 * qJD(4);
t331 = (-t344 * t311 - t314 * t339) * t305 - t338 * qJD(4);
t250 = t286 * mrSges(6,2) + t295 * t336 - t322;
t251 = -t286 * mrSges(6,3) - t284 * t336 + t327;
t319 = -mrSges(4,2) * t264 + t314 * (-mrSges(5,1) * t260 + mrSges(5,3) * t258 - mrSges(6,1) * t256 + mrSges(6,3) * t255 - pkin(4) * t250 + qJ(5) * t326 - t333 * t336 + t343 * t287 + t339 * t286 + (-qJ(5) * mrSges(6,2) + t337) * qJDD(4) + (-qJ(5) * t295 - t331) * qJD(4)) + t311 * (mrSges(5,2) * t260 + mrSges(6,2) * t256 - mrSges(5,3) * t257 - mrSges(6,3) * t254 - qJ(5) * t251 - t332 * qJD(4) + t338 * qJDD(4) + t344 * t286 + t339 * t287 + t333 * t335) + pkin(7) * t323 + pkin(3) * t318 + mrSges(4,1) * t263 + Ifges(4,3) * t304;
t1 = [pkin(1) * (t309 * (m(3) * t271 - t317 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t315 * t242 - t312 * t245) + t310 * (m(3) * t270 + qJDD(1) * mrSges(3,1) - t317 * mrSges(3,2) + t334)) + mrSges(2,1) * t321 - mrSges(2,2) * t324 + pkin(2) * t334 + mrSges(3,1) * t270 - mrSges(3,2) * t271 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t319; t314 * t248 + t311 * t249 + (m(3) + m(4)) * t308; t319; mrSges(5,1) * t257 + mrSges(6,1) * t254 - mrSges(5,2) * t258 - mrSges(6,2) * t255 + pkin(4) * t251 + t337 * t287 + t338 * t286 + t342 * qJDD(4) + (t311 * t332 + t314 * t331) * t305; t250;];
tauJ = t1;
