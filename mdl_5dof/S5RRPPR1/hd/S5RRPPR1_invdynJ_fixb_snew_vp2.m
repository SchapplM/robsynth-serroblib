% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR1
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:31
% EndTime: 2022-01-20 09:51:31
% DurationCPUTime: 0.72s
% Computational Cost: add. (5699->135), mult. (7766->178), div. (0->0), fcn. (4499->10), ass. (0->72)
t312 = cos(pkin(9));
t307 = t312 ^ 2;
t308 = qJD(1) + qJD(2);
t304 = t308 ^ 2;
t339 = pkin(4) * t312;
t310 = sin(pkin(9));
t338 = mrSges(5,2) * t310;
t336 = t304 * t307;
t305 = qJDD(1) + qJDD(2);
t335 = t305 * t312;
t316 = sin(qJ(1));
t319 = cos(qJ(1));
t330 = t316 * g(1) - t319 * g(2);
t292 = qJDD(1) * pkin(1) + t330;
t327 = -t319 * g(1) - t316 * g(2);
t293 = -qJD(1) ^ 2 * pkin(1) + t327;
t315 = sin(qJ(2));
t318 = cos(qJ(2));
t281 = t318 * t292 - t315 * t293;
t278 = t305 * pkin(2) + t281;
t282 = t315 * t292 + t318 * t293;
t279 = -t304 * pkin(2) + t282;
t311 = sin(pkin(8));
t313 = cos(pkin(8));
t266 = t311 * t278 + t313 * t279;
t263 = -t304 * pkin(3) + t305 * qJ(4) + t266;
t309 = -g(3) + qJDD(3);
t331 = qJD(4) * t308;
t332 = t312 * t309 - 0.2e1 * t310 * t331;
t259 = -t310 * t263 + t332;
t325 = mrSges(5,3) * t305 + (-mrSges(5,1) * t312 + t338) * t304;
t256 = (-pkin(7) * t305 + t304 * t339 - t263) * t310 + t332;
t260 = t310 * t309 + (t263 + 0.2e1 * t331) * t312;
t257 = -pkin(4) * t336 + pkin(7) * t335 + t260;
t314 = sin(qJ(5));
t317 = cos(qJ(5));
t254 = t317 * t256 - t314 * t257;
t323 = -t310 * t314 + t312 * t317;
t285 = t323 * t308;
t324 = t310 * t317 + t312 * t314;
t286 = t324 * t308;
t272 = -t285 * mrSges(6,1) + t286 * mrSges(6,2);
t274 = t285 * qJD(5) + t324 * t305;
t283 = -qJD(5) * mrSges(6,2) + t285 * mrSges(6,3);
t252 = m(6) * t254 + qJDD(5) * mrSges(6,1) - t274 * mrSges(6,3) + qJD(5) * t283 - t286 * t272;
t255 = t314 * t256 + t317 * t257;
t273 = -t286 * qJD(5) + t323 * t305;
t284 = qJD(5) * mrSges(6,1) - t286 * mrSges(6,3);
t253 = m(6) * t255 - qJDD(5) * mrSges(6,2) + t273 * mrSges(6,3) - qJD(5) * t284 + t285 * t272;
t333 = t317 * t252 + t314 * t253;
t241 = m(5) * t259 - t325 * t310 + t333;
t328 = -t314 * t252 + t317 * t253;
t242 = m(5) * t260 + t325 * t312 + t328;
t329 = -t310 * t241 + t312 * t242;
t238 = m(4) * t266 - t304 * mrSges(4,1) - t305 * mrSges(4,2) + t329;
t265 = t313 * t278 - t311 * t279;
t326 = qJDD(4) - t265;
t262 = -t305 * pkin(3) - t304 * qJ(4) + t326;
t306 = t310 ^ 2;
t258 = (-pkin(3) - t339) * t305 + (-qJ(4) + (-t306 - t307) * pkin(7)) * t304 + t326;
t322 = m(6) * t258 - t273 * mrSges(6,1) + t274 * mrSges(6,2) - t285 * t283 + t286 * t284;
t320 = -m(5) * t262 + mrSges(5,1) * t335 - t322 + (t304 * t306 + t336) * mrSges(5,3);
t246 = m(4) * t265 - t304 * mrSges(4,2) + (mrSges(4,1) - t338) * t305 + t320;
t334 = t311 * t238 + t313 * t246;
t267 = Ifges(6,5) * t286 + Ifges(6,6) * t285 + Ifges(6,3) * qJD(5);
t269 = Ifges(6,1) * t286 + Ifges(6,4) * t285 + Ifges(6,5) * qJD(5);
t243 = -mrSges(6,1) * t258 + mrSges(6,3) * t255 + Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * qJDD(5) + qJD(5) * t269 - t286 * t267;
t268 = Ifges(6,4) * t286 + Ifges(6,2) * t285 + Ifges(6,6) * qJD(5);
t244 = mrSges(6,2) * t258 - mrSges(6,3) * t254 + Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * qJDD(5) - qJD(5) * t268 + t285 * t267;
t248 = t305 * t338 - t320;
t321 = -mrSges(3,2) * t282 - mrSges(4,2) * t266 + t312 * (-mrSges(5,1) * t262 + mrSges(5,3) * t260 - pkin(4) * t322 + pkin(7) * t328 + t317 * t243 + t314 * t244) + pkin(2) * t334 + t310 * (mrSges(5,2) * t262 - mrSges(5,3) * t259 - pkin(7) * t333 - t314 * t243 + t317 * t244) + qJ(4) * t329 - pkin(3) * t248 + mrSges(4,1) * t265 + mrSges(3,1) * t281 + (Ifges(5,2) * t307 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t310 + 0.2e1 * Ifges(5,4) * t312) * t310) * t305;
t1 = [pkin(1) * (t315 * (m(3) * t282 - t304 * mrSges(3,1) - t305 * mrSges(3,2) + t313 * t238 - t311 * t246) + t318 * (m(3) * t281 + t305 * mrSges(3,1) - t304 * mrSges(3,2) + t334)) - mrSges(2,2) * t327 + mrSges(2,1) * t330 + t321 + Ifges(2,3) * qJDD(1); t321; m(4) * t309 + t312 * t241 + t310 * t242; t248; mrSges(6,1) * t254 - mrSges(6,2) * t255 + Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * qJDD(5) + t286 * t268 - t285 * t269;];
tauJ = t1;
