% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 14:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:43:11
% EndTime: 2019-05-06 14:43:38
% DurationCPUTime: 14.27s
% Computational Cost: add. (235902->381), mult. (523798->473), div. (0->0), fcn. (347006->10), ass. (0->145)
t362 = sin(qJ(1));
t366 = cos(qJ(1));
t335 = -t366 * g(1) - t362 * g(2);
t368 = qJD(1) ^ 2;
t308 = -t368 * pkin(1) + qJDD(1) * pkin(7) + t335;
t361 = sin(qJ(2));
t365 = cos(qJ(2));
t286 = -t365 * g(3) - t361 * t308;
t287 = -t361 * g(3) + t365 * t308;
t299 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t361 - Ifges(4,3) * t365) * qJD(1);
t302 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t361 + Ifges(3,2) * t365) * qJD(1);
t322 = (-mrSges(4,1) * t365 - mrSges(4,3) * t361) * qJD(1);
t389 = qJD(1) * qJD(2);
t387 = t365 * t389;
t324 = t361 * qJDD(1) + t387;
t388 = t361 * t389;
t325 = t365 * qJDD(1) - t388;
t321 = (-pkin(2) * t365 - qJ(3) * t361) * qJD(1);
t367 = qJD(2) ^ 2;
t390 = qJD(1) * t365;
t397 = 2 * qJD(3);
t262 = -t367 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t397 + t321 * t390 + t287;
t391 = qJD(1) * t361;
t333 = -qJD(2) * pkin(3) - pkin(8) * t391;
t394 = t365 ^ 2 * t368;
t255 = -pkin(3) * t394 - t325 * pkin(8) + qJD(2) * t333 + t262;
t270 = -qJDD(2) * pkin(2) - t367 * qJ(3) + t321 * t391 + qJDD(3) - t286;
t256 = (-t324 + t387) * pkin(8) + (-t361 * t365 * t368 - qJDD(2)) * pkin(3) + t270;
t360 = sin(qJ(4));
t364 = cos(qJ(4));
t223 = -t360 * t255 + t364 * t256;
t305 = (-t360 * t361 - t364 * t365) * qJD(1);
t272 = t305 * qJD(4) + t364 * t324 - t360 * t325;
t306 = (-t360 * t365 + t361 * t364) * qJD(1);
t346 = -qJDD(2) + qJDD(4);
t347 = -qJD(2) + qJD(4);
t215 = (t305 * t347 - t272) * qJ(5) + (t305 * t306 + t346) * pkin(4) + t223;
t224 = t364 * t255 + t360 * t256;
t271 = -t306 * qJD(4) - t360 * t324 - t364 * t325;
t289 = t347 * pkin(4) - t306 * qJ(5);
t298 = t305 ^ 2;
t217 = -t298 * pkin(4) + t271 * qJ(5) - t347 * t289 + t224;
t356 = sin(pkin(10));
t357 = cos(pkin(10));
t283 = t357 * t305 - t356 * t306;
t396 = 2 * qJD(5);
t212 = t356 * t215 + t357 * t217 + t283 * t396;
t238 = t357 * t271 - t356 * t272;
t284 = t356 * t305 + t357 * t306;
t253 = -t283 * mrSges(6,1) + t284 * mrSges(6,2);
t274 = t347 * mrSges(6,1) - t284 * mrSges(6,3);
t254 = -t283 * pkin(5) - t284 * pkin(9);
t345 = t347 ^ 2;
t208 = -t345 * pkin(5) + t346 * pkin(9) + t283 * t254 + t212;
t334 = t362 * g(1) - t366 * g(2);
t307 = -qJDD(1) * pkin(1) - t368 * pkin(7) - t334;
t381 = -t325 * pkin(2) + t307 + (-t324 - t387) * qJ(3);
t242 = -pkin(2) * t388 + t325 * pkin(3) - pkin(8) * t394 - t381 + (t333 + t397) * t391;
t219 = -t271 * pkin(4) - t298 * qJ(5) + t306 * t289 + qJDD(5) + t242;
t239 = t356 * t271 + t357 * t272;
t213 = t219 + (-t283 * t347 - t239) * pkin(9) + (t284 * t347 - t238) * pkin(5);
t359 = sin(qJ(6));
t363 = cos(qJ(6));
t205 = -t359 * t208 + t363 * t213;
t268 = -t359 * t284 + t363 * t347;
t226 = t268 * qJD(6) + t363 * t239 + t359 * t346;
t237 = qJDD(6) - t238;
t269 = t363 * t284 + t359 * t347;
t240 = -t268 * mrSges(7,1) + t269 * mrSges(7,2);
t276 = qJD(6) - t283;
t243 = -t276 * mrSges(7,2) + t268 * mrSges(7,3);
t201 = m(7) * t205 + t237 * mrSges(7,1) - t226 * mrSges(7,3) - t269 * t240 + t276 * t243;
t206 = t363 * t208 + t359 * t213;
t225 = -t269 * qJD(6) - t359 * t239 + t363 * t346;
t244 = t276 * mrSges(7,1) - t269 * mrSges(7,3);
t202 = m(7) * t206 - t237 * mrSges(7,2) + t225 * mrSges(7,3) + t268 * t240 - t276 * t244;
t384 = -t359 * t201 + t363 * t202;
t187 = m(6) * t212 - t346 * mrSges(6,2) + t238 * mrSges(6,3) + t283 * t253 - t347 * t274 + t384;
t383 = -t357 * t215 + t356 * t217;
t211 = -0.2e1 * qJD(5) * t284 - t383;
t273 = -t347 * mrSges(6,2) + t283 * mrSges(6,3);
t207 = -t346 * pkin(5) - t345 * pkin(9) + (t396 + t254) * t284 + t383;
t378 = -m(7) * t207 + t225 * mrSges(7,1) - t226 * mrSges(7,2) + t268 * t243 - t269 * t244;
t197 = m(6) * t211 + t346 * mrSges(6,1) - t239 * mrSges(6,3) - t284 * t253 + t347 * t273 + t378;
t181 = t356 * t187 + t357 * t197;
t285 = -t305 * mrSges(5,1) + t306 * mrSges(5,2);
t288 = -t347 * mrSges(5,2) + t305 * mrSges(5,3);
t178 = m(5) * t223 + t346 * mrSges(5,1) - t272 * mrSges(5,3) - t306 * t285 + t347 * t288 + t181;
t290 = t347 * mrSges(5,1) - t306 * mrSges(5,3);
t385 = t357 * t187 - t356 * t197;
t179 = m(5) * t224 - t346 * mrSges(5,2) + t271 * mrSges(5,3) + t305 * t285 - t347 * t290 + t385;
t173 = t364 * t178 + t360 * t179;
t278 = Ifges(5,4) * t306 + Ifges(5,2) * t305 + Ifges(5,6) * t347;
t279 = Ifges(5,1) * t306 + Ifges(5,4) * t305 + Ifges(5,5) * t347;
t227 = Ifges(7,5) * t269 + Ifges(7,6) * t268 + Ifges(7,3) * t276;
t229 = Ifges(7,1) * t269 + Ifges(7,4) * t268 + Ifges(7,5) * t276;
t194 = -mrSges(7,1) * t207 + mrSges(7,3) * t206 + Ifges(7,4) * t226 + Ifges(7,2) * t225 + Ifges(7,6) * t237 - t269 * t227 + t276 * t229;
t228 = Ifges(7,4) * t269 + Ifges(7,2) * t268 + Ifges(7,6) * t276;
t195 = mrSges(7,2) * t207 - mrSges(7,3) * t205 + Ifges(7,1) * t226 + Ifges(7,4) * t225 + Ifges(7,5) * t237 + t268 * t227 - t276 * t228;
t246 = Ifges(6,4) * t284 + Ifges(6,2) * t283 + Ifges(6,6) * t347;
t247 = Ifges(6,1) * t284 + Ifges(6,4) * t283 + Ifges(6,5) * t347;
t376 = mrSges(6,1) * t211 - mrSges(6,2) * t212 + Ifges(6,5) * t239 + Ifges(6,6) * t238 + Ifges(6,3) * t346 + pkin(5) * t378 + pkin(9) * t384 + t363 * t194 + t359 * t195 + t284 * t246 - t283 * t247;
t372 = mrSges(5,1) * t223 - mrSges(5,2) * t224 + Ifges(5,5) * t272 + Ifges(5,6) * t271 + Ifges(5,3) * t346 + pkin(4) * t181 + t306 * t278 - t305 * t279 + t376;
t371 = -mrSges(4,1) * t270 + mrSges(4,3) * t262 + Ifges(4,4) * t324 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t325 - pkin(3) * t173 - t372;
t332 = mrSges(4,2) * t390 + qJD(2) * mrSges(4,3);
t377 = -m(4) * t270 + qJDD(2) * mrSges(4,1) + qJD(2) * t332 - t173;
t174 = -t360 * t178 + t364 * t179;
t330 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t391;
t380 = m(4) * t262 + qJDD(2) * mrSges(4,3) + qJD(2) * t330 + t322 * t390 + t174;
t303 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t361 - Ifges(4,5) * t365) * qJD(1);
t392 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t361 + Ifges(3,4) * t365) * qJD(1) + t303;
t399 = -qJD(1) * ((t299 - t302) * t361 + t392 * t365) + mrSges(3,1) * t286 - mrSges(3,2) * t287 + Ifges(3,5) * t324 + Ifges(3,6) * t325 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t324 * mrSges(4,2) - t322 * t391 + t377) + qJ(3) * (t325 * mrSges(4,2) + t380) + t371;
t395 = mrSges(3,3) + mrSges(4,2);
t190 = t363 * t201 + t359 * t202;
t323 = (-mrSges(3,1) * t365 + mrSges(3,2) * t361) * qJD(1);
t329 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t391;
t169 = m(3) * t287 - qJDD(2) * mrSges(3,2) - qJD(2) * t329 + t323 * t390 + t395 * t325 + t380;
t331 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t390;
t170 = m(3) * t286 + qJDD(2) * mrSges(3,1) + qJD(2) * t331 - t395 * t324 + (-t322 - t323) * t391 + t377;
t386 = t365 * t169 - t361 * t170;
t382 = m(6) * t219 - t238 * mrSges(6,1) + t239 * mrSges(6,2) - t283 * t273 + t284 * t274 + t190;
t188 = -m(5) * t242 + t271 * mrSges(5,1) - t272 * mrSges(5,2) + t305 * t288 - t306 * t290 - t382;
t257 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t391 + t381;
t184 = m(4) * t257 - t325 * mrSges(4,1) - t324 * mrSges(4,3) - t330 * t391 - t332 * t390 + t188;
t300 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t361 + Ifges(3,6) * t365) * qJD(1);
t301 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t361 - Ifges(4,6) * t365) * qJD(1);
t245 = Ifges(6,5) * t284 + Ifges(6,6) * t283 + Ifges(6,3) * t347;
t175 = mrSges(6,2) * t219 - mrSges(6,3) * t211 + Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t346 - pkin(9) * t190 - t359 * t194 + t363 * t195 + t283 * t245 - t347 * t246;
t373 = mrSges(7,1) * t205 - mrSges(7,2) * t206 + Ifges(7,5) * t226 + Ifges(7,6) * t225 + Ifges(7,3) * t237 + t269 * t228 - t268 * t229;
t176 = -mrSges(6,1) * t219 + mrSges(6,3) * t212 + Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t346 - pkin(5) * t190 - t284 * t245 + t347 * t247 - t373;
t277 = Ifges(5,5) * t306 + Ifges(5,6) * t305 + Ifges(5,3) * t347;
t162 = -mrSges(5,1) * t242 + mrSges(5,3) * t224 + Ifges(5,4) * t272 + Ifges(5,2) * t271 + Ifges(5,6) * t346 - pkin(4) * t382 + qJ(5) * t385 + t356 * t175 + t357 * t176 - t306 * t277 + t347 * t279;
t167 = mrSges(5,2) * t242 - mrSges(5,3) * t223 + Ifges(5,1) * t272 + Ifges(5,4) * t271 + Ifges(5,5) * t346 - qJ(5) * t181 + t357 * t175 - t356 * t176 + t305 * t277 - t347 * t278;
t374 = -mrSges(4,1) * t257 + mrSges(4,2) * t262 - pkin(3) * t188 - pkin(8) * t174 - t364 * t162 - t360 * t167;
t159 = -mrSges(3,1) * t307 + mrSges(3,3) * t287 - pkin(2) * t184 + (Ifges(3,2) + Ifges(4,3)) * t325 + (Ifges(3,4) - Ifges(4,5)) * t324 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t392 * qJD(2) + (-t300 - t301) * t391 + t374;
t375 = mrSges(4,2) * t270 - mrSges(4,3) * t257 + Ifges(4,1) * t324 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t325 - pkin(8) * t173 + qJD(2) * t299 - t360 * t162 + t364 * t167 + t301 * t390;
t161 = mrSges(3,2) * t307 - mrSges(3,3) * t286 + Ifges(3,1) * t324 + Ifges(3,4) * t325 + Ifges(3,5) * qJDD(2) - qJ(3) * t184 - qJD(2) * t302 + t300 * t390 + t375;
t370 = -m(3) * t307 + t325 * mrSges(3,1) - t324 * mrSges(3,2) - t329 * t391 + t331 * t390 - t184;
t379 = mrSges(2,1) * t334 - mrSges(2,2) * t335 + Ifges(2,3) * qJDD(1) + pkin(1) * t370 + pkin(7) * t386 + t365 * t159 + t361 * t161;
t182 = m(2) * t334 + qJDD(1) * mrSges(2,1) - t368 * mrSges(2,2) + t370;
t165 = t361 * t169 + t365 * t170;
t163 = m(2) * t335 - t368 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t386;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t335 + t368 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t399;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t334 + Ifges(2,5) * qJDD(1) - t368 * Ifges(2,6) - pkin(7) * t165 - t361 * t159 + t365 * t161;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t366 * t156 - t362 * t157 - pkin(6) * (t362 * t163 + t366 * t182), t156, t161, t375, t167, t175, t195; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t362 * t156 + t366 * t157 + pkin(6) * (t366 * t163 - t362 * t182), t157, t159, t371 + (-t361 * t299 - t365 * t303) * qJD(1), t162, t176, t194; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t379, t379, t399, Ifges(4,5) * t324 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t325 - qJD(2) * t303 + t301 * t391 - t374, t372, t376, t373;];
m_new  = t1;
