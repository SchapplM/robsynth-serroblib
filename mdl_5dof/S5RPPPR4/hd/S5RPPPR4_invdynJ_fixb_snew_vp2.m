% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:11
% EndTime: 2019-12-31 17:45:12
% DurationCPUTime: 0.64s
% Computational Cost: add. (1857->131), mult. (3622->170), div. (0->0), fcn. (2022->8), ass. (0->67)
t299 = sin(pkin(8));
t292 = t299 ^ 2;
t301 = cos(pkin(8));
t293 = t301 ^ 2;
t326 = -t292 - t293;
t335 = t326 * mrSges(5,3);
t304 = sin(qJ(1));
t306 = cos(qJ(1));
t319 = t304 * g(1) - g(2) * t306;
t283 = qJDD(1) * pkin(1) + t319;
t307 = qJD(1) ^ 2;
t317 = -g(1) * t306 - g(2) * t304;
t284 = -pkin(1) * t307 + t317;
t300 = sin(pkin(7));
t302 = cos(pkin(7));
t272 = t283 * t302 - t300 * t284;
t311 = -qJ(3) * t307 + qJDD(3) - t272;
t329 = -pkin(2) - qJ(4);
t334 = -(2 * qJD(1) * qJD(4)) + t329 * qJDD(1) + t311;
t273 = t300 * t283 + t302 * t284;
t333 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t273;
t262 = pkin(2) * t307 - t333;
t332 = -m(4) * t262 + t307 * mrSges(4,2) + qJDD(1) * mrSges(4,3);
t330 = pkin(4) * t307;
t296 = -g(3) + qJDD(2);
t322 = t301 * qJDD(1);
t327 = t334 * t301;
t251 = -pkin(6) * t322 + (-t301 * t330 - t296) * t299 + t327;
t256 = t301 * t296 + t334 * t299;
t323 = t299 * qJDD(1);
t252 = -pkin(6) * t323 - t292 * t330 + t256;
t303 = sin(qJ(5));
t305 = cos(qJ(5));
t249 = t251 * t305 - t252 * t303;
t315 = -t299 * t305 - t301 * t303;
t276 = t315 * qJD(1);
t314 = -t299 * t303 + t301 * t305;
t277 = t314 * qJD(1);
t268 = -mrSges(6,1) * t276 + mrSges(6,2) * t277;
t271 = qJD(5) * t276 + t314 * qJDD(1);
t274 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t276;
t247 = m(6) * t249 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t271 + qJD(5) * t274 - t268 * t277;
t250 = t251 * t303 + t252 * t305;
t270 = -qJD(5) * t277 + t315 * qJDD(1);
t275 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t277;
t248 = m(6) * t250 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t270 - qJD(5) * t275 + t268 * t276;
t328 = t305 * t247 + t303 * t248;
t318 = -t247 * t303 + t305 * t248;
t255 = -t296 * t299 + t327;
t312 = -qJDD(1) * mrSges(5,3) - t307 * (t299 * mrSges(5,1) + t301 * mrSges(5,2));
t240 = m(5) * t255 + t312 * t301 + t328;
t241 = m(5) * t256 + t312 * t299 + t318;
t316 = t240 * t301 + t241 * t299;
t313 = qJDD(4) + t333;
t254 = pkin(4) * t323 + (t326 * pkin(6) + t329) * t307 + t313;
t310 = m(6) * t254 - mrSges(6,1) * t270 + t271 * mrSges(6,2) - t274 * t276 + t277 * t275;
t260 = t329 * t307 + t313;
t309 = m(5) * t260 + mrSges(5,1) * t323 + mrSges(5,2) * t322 + t310;
t263 = -qJDD(1) * pkin(2) + t311;
t239 = m(4) * t263 + qJDD(1) * mrSges(4,2) - t307 * mrSges(4,3) + t316;
t308 = t307 * t335 + t309;
t266 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * qJD(5);
t265 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * qJD(5);
t264 = Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * qJD(5);
t243 = mrSges(6,2) * t254 - mrSges(6,3) * t249 + Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * qJDD(5) - qJD(5) * t265 + t264 * t276;
t242 = -mrSges(6,1) * t254 + mrSges(6,3) * t250 + Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * qJDD(5) + qJD(5) * t266 - t264 * t277;
t1 = [pkin(1) * (t300 * (m(3) * t273 + t309 + t332) + t302 * (m(3) * t272 - t239) + (t300 * (-mrSges(3,1) + t335) - t302 * mrSges(3,2)) * t307) + mrSges(2,1) * t319 - mrSges(2,2) * t317 - pkin(2) * t239 + qJ(3) * (t308 + t332) - t299 * (-mrSges(5,1) * t260 + mrSges(5,3) * t256 - pkin(4) * t310 + pkin(6) * t318 + t305 * t242 + t303 * t243) - qJ(4) * t316 + mrSges(3,1) * t272 - mrSges(3,2) * t273 + mrSges(4,2) * t263 - mrSges(4,3) * t262 + t301 * (mrSges(5,2) * t260 - mrSges(5,3) * t255 - pkin(6) * t328 - t303 * t242 + t305 * t243) + (pkin(1) * (t302 * mrSges(3,1) - t300 * mrSges(3,2)) + Ifges(5,1) * t293 + Ifges(2,3) + Ifges(3,3) + Ifges(4,1) + (-0.2e1 * Ifges(5,4) * t301 + Ifges(5,2) * t299) * t299) * qJDD(1); -t240 * t299 + t241 * t301 + (m(3) + m(4)) * t296; t239; t308; mrSges(6,1) * t249 - mrSges(6,2) * t250 + Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * qJDD(5) + t265 * t277 - t266 * t276;];
tauJ = t1;
