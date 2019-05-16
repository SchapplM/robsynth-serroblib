% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-06 22:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:16:28
% EndTime: 2019-05-06 22:16:58
% DurationCPUTime: 14.40s
% Computational Cost: add. (258606->382), mult. (524755->470), div. (0->0), fcn. (343430->10), ass. (0->145)
t363 = sin(qJ(1));
t368 = cos(qJ(1));
t336 = -t368 * g(1) - t363 * g(2);
t370 = qJD(1) ^ 2;
t308 = -t370 * pkin(1) + qJDD(1) * pkin(7) + t336;
t362 = sin(qJ(2));
t367 = cos(qJ(2));
t285 = -t367 * g(3) - t362 * t308;
t286 = -t362 * g(3) + t367 * t308;
t299 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t362 - Ifges(4,3) * t367) * qJD(1);
t302 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t362 + Ifges(3,2) * t367) * qJD(1);
t323 = (-mrSges(4,1) * t367 - mrSges(4,3) * t362) * qJD(1);
t390 = qJD(1) * qJD(2);
t388 = t367 * t390;
t325 = t362 * qJDD(1) + t388;
t389 = t362 * t390;
t326 = t367 * qJDD(1) - t389;
t322 = (-pkin(2) * t367 - qJ(3) * t362) * qJD(1);
t369 = qJD(2) ^ 2;
t391 = qJD(1) * t367;
t397 = 2 * qJD(3);
t262 = -t369 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t397 + t322 * t391 + t286;
t392 = qJD(1) * t362;
t334 = -qJD(2) * pkin(3) - pkin(8) * t392;
t395 = t367 ^ 2 * t370;
t245 = -pkin(3) * t395 - t326 * pkin(8) + qJD(2) * t334 + t262;
t269 = -qJDD(2) * pkin(2) - t369 * qJ(3) + t322 * t392 + qJDD(3) - t285;
t246 = (-t325 + t388) * pkin(8) + (-t362 * t367 * t370 - qJDD(2)) * pkin(3) + t269;
t361 = sin(qJ(4));
t366 = cos(qJ(4));
t229 = t366 * t245 + t361 * t246;
t306 = (-t361 * t367 + t362 * t366) * qJD(1);
t270 = -t306 * qJD(4) - t361 * t325 - t366 * t326;
t305 = -t361 * t392 - t366 * t391;
t280 = -t305 * mrSges(5,1) + t306 * mrSges(5,2);
t349 = -qJD(2) + qJD(4);
t288 = t349 * mrSges(5,1) - t306 * mrSges(5,3);
t348 = -qJDD(2) + qJDD(4);
t335 = t363 * g(1) - t368 * g(2);
t307 = -qJDD(1) * pkin(1) - t370 * pkin(7) - t335;
t384 = -t326 * pkin(2) + t307 + (-t325 - t388) * qJ(3);
t236 = -pkin(2) * t389 + t326 * pkin(3) - pkin(8) * t395 - t384 + (t334 + t397) * t392;
t271 = t305 * qJD(4) + t366 * t325 - t361 * t326;
t221 = t236 + (t306 * t349 - t270) * pkin(4) + (-t305 * t349 - t271) * pkin(9);
t281 = -t305 * pkin(4) - t306 * pkin(9);
t347 = t349 ^ 2;
t224 = -t347 * pkin(4) + t348 * pkin(9) + t305 * t281 + t229;
t360 = sin(qJ(5));
t365 = cos(qJ(5));
t211 = t365 * t221 - t360 * t224;
t283 = -t360 * t306 + t365 * t349;
t239 = t283 * qJD(5) + t365 * t271 + t360 * t348;
t268 = qJDD(5) - t270;
t284 = t365 * t306 + t360 * t349;
t298 = qJD(5) - t305;
t209 = (t283 * t298 - t239) * pkin(10) + (t283 * t284 + t268) * pkin(5) + t211;
t212 = t360 * t221 + t365 * t224;
t238 = -t284 * qJD(5) - t360 * t271 + t365 * t348;
t274 = t298 * pkin(5) - t284 * pkin(10);
t282 = t283 ^ 2;
t210 = -t282 * pkin(5) + t238 * pkin(10) - t298 * t274 + t212;
t359 = sin(qJ(6));
t364 = cos(qJ(6));
t207 = t364 * t209 - t359 * t210;
t252 = t364 * t283 - t359 * t284;
t218 = t252 * qJD(6) + t359 * t238 + t364 * t239;
t253 = t359 * t283 + t364 * t284;
t234 = -t252 * mrSges(7,1) + t253 * mrSges(7,2);
t292 = qJD(6) + t298;
t243 = -t292 * mrSges(7,2) + t252 * mrSges(7,3);
t260 = qJDD(6) + t268;
t202 = m(7) * t207 + t260 * mrSges(7,1) - t218 * mrSges(7,3) - t253 * t234 + t292 * t243;
t208 = t359 * t209 + t364 * t210;
t217 = -t253 * qJD(6) + t364 * t238 - t359 * t239;
t244 = t292 * mrSges(7,1) - t253 * mrSges(7,3);
t203 = m(7) * t208 - t260 * mrSges(7,2) + t217 * mrSges(7,3) + t252 * t234 - t292 * t244;
t195 = t364 * t202 + t359 * t203;
t254 = -t283 * mrSges(6,1) + t284 * mrSges(6,2);
t272 = -t298 * mrSges(6,2) + t283 * mrSges(6,3);
t193 = m(6) * t211 + t268 * mrSges(6,1) - t239 * mrSges(6,3) - t284 * t254 + t298 * t272 + t195;
t273 = t298 * mrSges(6,1) - t284 * mrSges(6,3);
t385 = -t359 * t202 + t364 * t203;
t194 = m(6) * t212 - t268 * mrSges(6,2) + t238 * mrSges(6,3) + t283 * t254 - t298 * t273 + t385;
t386 = -t360 * t193 + t365 * t194;
t185 = m(5) * t229 - t348 * mrSges(5,2) + t270 * mrSges(5,3) + t305 * t280 - t349 * t288 + t386;
t228 = -t361 * t245 + t366 * t246;
t287 = -t349 * mrSges(5,2) + t305 * mrSges(5,3);
t223 = -t348 * pkin(4) - t347 * pkin(9) + t306 * t281 - t228;
t213 = -t238 * pkin(5) - t282 * pkin(10) + t284 * t274 + t223;
t381 = m(7) * t213 - t217 * mrSges(7,1) + t218 * mrSges(7,2) - t252 * t243 + t253 * t244;
t375 = -m(6) * t223 + t238 * mrSges(6,1) - t239 * mrSges(6,2) + t283 * t272 - t284 * t273 - t381;
t198 = m(5) * t228 + t348 * mrSges(5,1) - t271 * mrSges(5,3) - t306 * t280 + t349 * t287 + t375;
t175 = t361 * t185 + t366 * t198;
t230 = Ifges(7,5) * t253 + Ifges(7,6) * t252 + Ifges(7,3) * t292;
t232 = Ifges(7,1) * t253 + Ifges(7,4) * t252 + Ifges(7,5) * t292;
t196 = -mrSges(7,1) * t213 + mrSges(7,3) * t208 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t260 - t253 * t230 + t292 * t232;
t231 = Ifges(7,4) * t253 + Ifges(7,2) * t252 + Ifges(7,6) * t292;
t197 = mrSges(7,2) * t213 - mrSges(7,3) * t207 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t260 + t252 * t230 - t292 * t231;
t247 = Ifges(6,5) * t284 + Ifges(6,6) * t283 + Ifges(6,3) * t298;
t249 = Ifges(6,1) * t284 + Ifges(6,4) * t283 + Ifges(6,5) * t298;
t178 = -mrSges(6,1) * t223 + mrSges(6,3) * t212 + Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t268 - pkin(5) * t381 + pkin(10) * t385 + t364 * t196 + t359 * t197 - t284 * t247 + t298 * t249;
t248 = Ifges(6,4) * t284 + Ifges(6,2) * t283 + Ifges(6,6) * t298;
t180 = mrSges(6,2) * t223 - mrSges(6,3) * t211 + Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t268 - pkin(10) * t195 - t359 * t196 + t364 * t197 + t283 * t247 - t298 * t248;
t276 = Ifges(5,4) * t306 + Ifges(5,2) * t305 + Ifges(5,6) * t349;
t277 = Ifges(5,1) * t306 + Ifges(5,4) * t305 + Ifges(5,5) * t349;
t378 = mrSges(5,1) * t228 - mrSges(5,2) * t229 + Ifges(5,5) * t271 + Ifges(5,6) * t270 + Ifges(5,3) * t348 + pkin(4) * t375 + pkin(9) * t386 + t365 * t178 + t360 * t180 + t306 * t276 - t305 * t277;
t374 = -mrSges(4,1) * t269 + mrSges(4,3) * t262 + Ifges(4,4) * t325 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t326 - pkin(3) * t175 - t378;
t333 = mrSges(4,2) * t391 + qJD(2) * mrSges(4,3);
t380 = -m(4) * t269 + qJDD(2) * mrSges(4,1) + qJD(2) * t333 - t175;
t176 = t366 * t185 - t361 * t198;
t331 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t392;
t383 = m(4) * t262 + qJDD(2) * mrSges(4,3) + qJD(2) * t331 + t323 * t391 + t176;
t303 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t362 - Ifges(4,5) * t367) * qJD(1);
t393 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t362 + Ifges(3,4) * t367) * qJD(1) + t303;
t399 = -((t299 - t302) * t362 + t393 * t367) * qJD(1) + mrSges(3,1) * t285 - mrSges(3,2) * t286 + Ifges(3,5) * t325 + Ifges(3,6) * t326 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t325 * mrSges(4,2) - t323 * t392 + t380) + qJ(3) * (t326 * mrSges(4,2) + t383) + t374;
t396 = mrSges(3,3) + mrSges(4,2);
t188 = t365 * t193 + t360 * t194;
t324 = (-mrSges(3,1) * t367 + mrSges(3,2) * t362) * qJD(1);
t330 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t392;
t171 = m(3) * t286 - qJDD(2) * mrSges(3,2) - qJD(2) * t330 + t324 * t391 + t396 * t326 + t383;
t332 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t391;
t172 = m(3) * t285 + qJDD(2) * mrSges(3,1) + qJD(2) * t332 - t396 * t325 + (-t323 - t324) * t392 + t380;
t387 = t367 * t171 - t362 * t172;
t186 = -m(5) * t236 + t270 * mrSges(5,1) - t271 * mrSges(5,2) + t305 * t287 - t306 * t288 - t188;
t251 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t392 + t384;
t183 = m(4) * t251 - t326 * mrSges(4,1) - t325 * mrSges(4,3) - t331 * t392 - t333 * t391 + t186;
t300 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t362 + Ifges(3,6) * t367) * qJD(1);
t301 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t362 - Ifges(4,6) * t367) * qJD(1);
t275 = Ifges(5,5) * t306 + Ifges(5,6) * t305 + Ifges(5,3) * t349;
t165 = mrSges(5,2) * t236 - mrSges(5,3) * t228 + Ifges(5,1) * t271 + Ifges(5,4) * t270 + Ifges(5,5) * t348 - pkin(9) * t188 - t360 * t178 + t365 * t180 + t305 * t275 - t349 * t276;
t379 = -mrSges(7,1) * t207 + mrSges(7,2) * t208 - Ifges(7,5) * t218 - Ifges(7,6) * t217 - Ifges(7,3) * t260 - t253 * t231 + t252 * t232;
t372 = mrSges(6,1) * t211 - mrSges(6,2) * t212 + Ifges(6,5) * t239 + Ifges(6,6) * t238 + Ifges(6,3) * t268 + pkin(5) * t195 + t284 * t248 - t283 * t249 - t379;
t169 = -mrSges(5,1) * t236 + mrSges(5,3) * t229 + Ifges(5,4) * t271 + Ifges(5,2) * t270 + Ifges(5,6) * t348 - pkin(4) * t188 - t306 * t275 + t349 * t277 - t372;
t376 = -mrSges(4,1) * t251 + mrSges(4,2) * t262 - pkin(3) * t186 - pkin(8) * t176 - t361 * t165 - t366 * t169;
t161 = -mrSges(3,1) * t307 + mrSges(3,3) * t286 - pkin(2) * t183 + (Ifges(3,2) + Ifges(4,3)) * t326 + (Ifges(3,4) - Ifges(4,5)) * t325 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t393 * qJD(2) + (-t300 - t301) * t392 + t376;
t377 = mrSges(4,2) * t269 - mrSges(4,3) * t251 + Ifges(4,1) * t325 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t326 - pkin(8) * t175 + qJD(2) * t299 + t366 * t165 - t361 * t169 + t301 * t391;
t163 = mrSges(3,2) * t307 - mrSges(3,3) * t285 + Ifges(3,1) * t325 + Ifges(3,4) * t326 + Ifges(3,5) * qJDD(2) - qJ(3) * t183 - qJD(2) * t302 + t300 * t391 + t377;
t373 = -m(3) * t307 + t326 * mrSges(3,1) - t325 * mrSges(3,2) - t330 * t392 + t332 * t391 - t183;
t382 = mrSges(2,1) * t335 - mrSges(2,2) * t336 + Ifges(2,3) * qJDD(1) + pkin(1) * t373 + pkin(7) * t387 + t367 * t161 + t362 * t163;
t181 = m(2) * t335 + qJDD(1) * mrSges(2,1) - t370 * mrSges(2,2) + t373;
t168 = t362 * t171 + t367 * t172;
t166 = m(2) * t336 - t370 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t387;
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t336 + t370 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t168 - t399;
t158 = -mrSges(2,2) * g(3) - mrSges(2,3) * t335 + Ifges(2,5) * qJDD(1) - t370 * Ifges(2,6) - pkin(7) * t168 - t362 * t161 + t367 * t163;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t368 * t158 - t363 * t159 - pkin(6) * (t363 * t166 + t368 * t181), t158, t163, t377, t165, t180, t197; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t363 * t158 + t368 * t159 + pkin(6) * (t368 * t166 - t363 * t181), t159, t161, t374 + (-t362 * t299 - t367 * t303) * qJD(1), t169, t178, t196; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t382, t382, t399, Ifges(4,5) * t325 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t326 - qJD(2) * t303 + t301 * t392 - t376, t378, t372, -t379;];
m_new  = t1;
