% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:32:00
% EndTime: 2019-12-05 18:32:01
% DurationCPUTime: 0.71s
% Computational Cost: add. (6933->167), mult. (8931->220), div. (0->0), fcn. (5058->10), ass. (0->77)
t318 = qJD(1) + qJD(2);
t324 = sin(qJ(4));
t344 = t318 * t324;
t328 = cos(qJ(4));
t343 = t318 * t328;
t326 = sin(qJ(1));
t330 = cos(qJ(1));
t341 = t330 * g(2) + t326 * g(3);
t301 = qJDD(1) * pkin(1) + t341;
t338 = t326 * g(2) - t330 * g(3);
t302 = -qJD(1) ^ 2 * pkin(1) + t338;
t325 = sin(qJ(2));
t329 = cos(qJ(2));
t284 = t329 * t301 - t325 * t302;
t316 = qJDD(1) + qJDD(2);
t281 = t316 * pkin(2) + t284;
t285 = t325 * t301 + t329 * t302;
t314 = t318 ^ 2;
t282 = -t314 * pkin(2) + t285;
t321 = sin(pkin(9));
t322 = cos(pkin(9));
t266 = t321 * t281 + t322 * t282;
t263 = -t314 * pkin(3) + t316 * pkin(7) + t266;
t320 = -g(1) + qJDD(3);
t259 = -t324 * t263 + t328 * t320;
t340 = qJD(4) * t318;
t339 = t328 * t340;
t296 = t324 * t316 + t339;
t256 = (-t296 + t339) * pkin(8) + (t314 * t324 * t328 + qJDD(4)) * pkin(4) + t259;
t260 = t328 * t263 + t324 * t320;
t297 = t328 * t316 - t324 * t340;
t305 = qJD(4) * pkin(4) - pkin(8) * t344;
t319 = t328 ^ 2;
t257 = -t319 * t314 * pkin(4) + t297 * pkin(8) - qJD(4) * t305 + t260;
t323 = sin(qJ(5));
t327 = cos(qJ(5));
t254 = t327 * t256 - t323 * t257;
t291 = (-t323 * t324 + t327 * t328) * t318;
t272 = t291 * qJD(5) + t327 * t296 + t323 * t297;
t292 = (t323 * t328 + t324 * t327) * t318;
t277 = -t291 * mrSges(6,1) + t292 * mrSges(6,2);
t317 = qJD(4) + qJD(5);
t286 = -t317 * mrSges(6,2) + t291 * mrSges(6,3);
t315 = qJDD(4) + qJDD(5);
t251 = m(6) * t254 + t315 * mrSges(6,1) - t272 * mrSges(6,3) - t292 * t277 + t317 * t286;
t255 = t323 * t256 + t327 * t257;
t271 = -t292 * qJD(5) - t323 * t296 + t327 * t297;
t287 = t317 * mrSges(6,1) - t292 * mrSges(6,3);
t252 = m(6) * t255 - t315 * mrSges(6,2) + t271 * mrSges(6,3) + t291 * t277 - t317 * t287;
t242 = t327 * t251 + t323 * t252;
t295 = (-mrSges(5,1) * t328 + mrSges(5,2) * t324) * t318;
t304 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t343;
t240 = m(5) * t259 + qJDD(4) * mrSges(5,1) - t296 * mrSges(5,3) + qJD(4) * t304 - t295 * t344 + t242;
t303 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t344;
t336 = -t323 * t251 + t327 * t252;
t241 = m(5) * t260 - qJDD(4) * mrSges(5,2) + t297 * mrSges(5,3) - qJD(4) * t303 + t295 * t343 + t336;
t337 = -t324 * t240 + t328 * t241;
t237 = m(4) * t266 - t314 * mrSges(4,1) - t316 * mrSges(4,2) + t337;
t265 = t322 * t281 - t321 * t282;
t335 = -t316 * pkin(3) - t265;
t262 = -t314 * pkin(7) + t335;
t258 = t305 * t344 - t297 * pkin(4) + (-pkin(8) * t319 - pkin(7)) * t314 + t335;
t334 = m(6) * t258 - t271 * mrSges(6,1) + t272 * mrSges(6,2) - t291 * t286 + t292 * t287;
t331 = -m(5) * t262 + t297 * mrSges(5,1) - t296 * mrSges(5,2) - t303 * t344 + t304 * t343 - t334;
t246 = m(4) * t265 + t316 * mrSges(4,1) - t314 * mrSges(4,2) + t331;
t342 = t321 * t237 + t322 * t246;
t274 = Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t317;
t275 = Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t317;
t333 = mrSges(6,1) * t254 - mrSges(6,2) * t255 + Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t315 + t292 * t274 - t291 * t275;
t273 = Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * t317;
t243 = -mrSges(6,1) * t258 + mrSges(6,3) * t255 + Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t315 - t292 * t273 + t317 * t275;
t244 = mrSges(6,2) * t258 - mrSges(6,3) * t254 + Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t315 + t291 * t273 - t317 * t274;
t288 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t324 + Ifges(5,6) * t328) * t318;
t289 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t324 + Ifges(5,2) * t328) * t318;
t290 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t324 + Ifges(5,4) * t328) * t318;
t332 = -mrSges(3,2) * t285 - mrSges(4,2) * t266 + t328 * (-mrSges(5,1) * t262 + mrSges(5,3) * t260 + Ifges(5,4) * t296 + Ifges(5,2) * t297 + Ifges(5,6) * qJDD(4) - pkin(4) * t334 + pkin(8) * t336 + qJD(4) * t290 + t327 * t243 + t323 * t244 - t288 * t344) + pkin(2) * t342 + t324 * (mrSges(5,2) * t262 - mrSges(5,3) * t259 + Ifges(5,1) * t296 + Ifges(5,4) * t297 + Ifges(5,5) * qJDD(4) - pkin(8) * t242 - qJD(4) * t289 - t323 * t243 + t327 * t244 + t288 * t343) + pkin(7) * t337 + pkin(3) * t331 + mrSges(4,1) * t265 + mrSges(3,1) * t284 + (Ifges(4,3) + Ifges(3,3)) * t316;
t1 = [t332 + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t341 + pkin(1) * (t325 * (m(3) * t285 - t314 * mrSges(3,1) - t316 * mrSges(3,2) + t322 * t237 - t321 * t246) + t329 * (m(3) * t284 + t316 * mrSges(3,1) - t314 * mrSges(3,2) + t342)) - mrSges(2,2) * t338; t332; m(4) * t320 + t328 * t240 + t324 * t241; mrSges(5,1) * t259 - mrSges(5,2) * t260 + Ifges(5,5) * t296 + Ifges(5,6) * t297 + Ifges(5,3) * qJDD(4) + pkin(4) * t242 + (t289 * t324 - t290 * t328) * t318 + t333; t333;];
tauJ = t1;
