% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 10:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:52:33
% EndTime: 2019-05-06 10:52:59
% DurationCPUTime: 13.66s
% Computational Cost: add. (223853->380), mult. (511562->473), div. (0->0), fcn. (338684->10), ass. (0->143)
t363 = sin(qJ(1));
t367 = cos(qJ(1));
t336 = -g(1) * t367 - g(2) * t363;
t369 = qJD(1) ^ 2;
t309 = -pkin(1) * t369 + qJDD(1) * pkin(7) + t336;
t362 = sin(qJ(2));
t366 = cos(qJ(2));
t285 = -t366 * g(3) - t362 * t309;
t286 = -g(3) * t362 + t366 * t309;
t302 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t362 - Ifges(4,3) * t366) * qJD(1);
t305 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t362 + Ifges(3,2) * t366) * qJD(1);
t323 = (-mrSges(4,1) * t366 - mrSges(4,3) * t362) * qJD(1);
t389 = qJD(1) * qJD(2);
t388 = t366 * t389;
t325 = qJDD(1) * t362 + t388;
t387 = t362 * t389;
t326 = qJDD(1) * t366 - t387;
t322 = (-pkin(2) * t366 - qJ(3) * t362) * qJD(1);
t368 = qJD(2) ^ 2;
t390 = qJD(1) * t366;
t396 = 2 * qJD(3);
t262 = -pkin(2) * t368 + qJDD(2) * qJ(3) + qJD(2) * t396 + t322 * t390 + t286;
t391 = qJD(1) * t362;
t330 = -qJD(2) * pkin(3) - qJ(4) * t391;
t394 = t366 ^ 2 * t369;
t255 = -pkin(3) * t394 - qJ(4) * t326 + qJD(2) * t330 + t262;
t266 = -qJDD(2) * pkin(2) - t368 * qJ(3) + t322 * t391 + qJDD(3) - t285;
t256 = (-t325 + t388) * qJ(4) + (-t362 * t366 * t369 - qJDD(2)) * pkin(3) + t266;
t357 = sin(pkin(10));
t358 = cos(pkin(10));
t301 = (-t357 * t366 + t358 * t362) * qJD(1);
t220 = -0.2e1 * qJD(4) * t301 - t357 * t255 + t358 * t256;
t284 = t325 * t358 - t326 * t357;
t300 = (-t357 * t362 - t358 * t366) * qJD(1);
t215 = (-qJD(2) * t300 - t284) * pkin(8) + (t300 * t301 - qJDD(2)) * pkin(4) + t220;
t221 = 0.2e1 * qJD(4) * t300 + t358 * t255 + t357 * t256;
t283 = -t325 * t357 - t326 * t358;
t289 = -qJD(2) * pkin(4) - pkin(8) * t301;
t299 = t300 ^ 2;
t217 = -pkin(4) * t299 + pkin(8) * t283 + qJD(2) * t289 + t221;
t361 = sin(qJ(5));
t365 = cos(qJ(5));
t212 = t361 * t215 + t365 * t217;
t274 = t300 * t361 + t301 * t365;
t238 = -qJD(5) * t274 + t283 * t365 - t284 * t361;
t273 = t300 * t365 - t301 * t361;
t253 = -mrSges(6,1) * t273 + mrSges(6,2) * t274;
t347 = -qJD(2) + qJD(5);
t268 = mrSges(6,1) * t347 - mrSges(6,3) * t274;
t346 = -qJDD(2) + qJDD(5);
t254 = -pkin(5) * t273 - pkin(9) * t274;
t345 = t347 ^ 2;
t208 = -pkin(5) * t345 + pkin(9) * t346 + t254 * t273 + t212;
t335 = t363 * g(1) - t367 * g(2);
t308 = -qJDD(1) * pkin(1) - t369 * pkin(7) - t335;
t382 = -t326 * pkin(2) + t308 + (-t325 - t388) * qJ(3);
t242 = -pkin(2) * t387 + t326 * pkin(3) - qJ(4) * t394 + qJDD(4) - t382 + (t330 + t396) * t391;
t223 = -t283 * pkin(4) - t299 * pkin(8) + t301 * t289 + t242;
t239 = qJD(5) * t273 + t283 * t361 + t284 * t365;
t213 = t223 + (-t273 * t347 - t239) * pkin(9) + (t274 * t347 - t238) * pkin(5);
t360 = sin(qJ(6));
t364 = cos(qJ(6));
t205 = -t208 * t360 + t213 * t364;
t263 = -t274 * t360 + t347 * t364;
t226 = qJD(6) * t263 + t239 * t364 + t346 * t360;
t237 = qJDD(6) - t238;
t264 = t274 * t364 + t347 * t360;
t240 = -mrSges(7,1) * t263 + mrSges(7,2) * t264;
t269 = qJD(6) - t273;
t243 = -mrSges(7,2) * t269 + mrSges(7,3) * t263;
t201 = m(7) * t205 + mrSges(7,1) * t237 - mrSges(7,3) * t226 - t240 * t264 + t243 * t269;
t206 = t208 * t364 + t213 * t360;
t225 = -qJD(6) * t264 - t239 * t360 + t346 * t364;
t244 = mrSges(7,1) * t269 - mrSges(7,3) * t264;
t202 = m(7) * t206 - mrSges(7,2) * t237 + mrSges(7,3) * t225 + t240 * t263 - t244 * t269;
t384 = -t201 * t360 + t364 * t202;
t187 = m(6) * t212 - mrSges(6,2) * t346 + mrSges(6,3) * t238 + t253 * t273 - t268 * t347 + t384;
t211 = t215 * t365 - t217 * t361;
t267 = -mrSges(6,2) * t347 + mrSges(6,3) * t273;
t207 = -pkin(5) * t346 - pkin(9) * t345 + t254 * t274 - t211;
t379 = -m(7) * t207 + t225 * mrSges(7,1) - mrSges(7,2) * t226 + t263 * t243 - t244 * t264;
t197 = m(6) * t211 + mrSges(6,1) * t346 - mrSges(6,3) * t239 - t253 * t274 + t267 * t347 + t379;
t181 = t361 * t187 + t365 * t197;
t278 = -mrSges(5,1) * t300 + mrSges(5,2) * t301;
t287 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t300;
t178 = m(5) * t220 - qJDD(2) * mrSges(5,1) - mrSges(5,3) * t284 - qJD(2) * t287 - t278 * t301 + t181;
t288 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t301;
t385 = t365 * t187 - t197 * t361;
t179 = m(5) * t221 + qJDD(2) * mrSges(5,2) + mrSges(5,3) * t283 + qJD(2) * t288 + t278 * t300 + t385;
t173 = t358 * t178 + t357 * t179;
t271 = Ifges(5,4) * t301 + Ifges(5,2) * t300 - Ifges(5,6) * qJD(2);
t272 = Ifges(5,1) * t301 + Ifges(5,4) * t300 - Ifges(5,5) * qJD(2);
t227 = Ifges(7,5) * t264 + Ifges(7,6) * t263 + Ifges(7,3) * t269;
t229 = Ifges(7,1) * t264 + Ifges(7,4) * t263 + Ifges(7,5) * t269;
t194 = -mrSges(7,1) * t207 + mrSges(7,3) * t206 + Ifges(7,4) * t226 + Ifges(7,2) * t225 + Ifges(7,6) * t237 - t227 * t264 + t229 * t269;
t228 = Ifges(7,4) * t264 + Ifges(7,2) * t263 + Ifges(7,6) * t269;
t195 = mrSges(7,2) * t207 - mrSges(7,3) * t205 + Ifges(7,1) * t226 + Ifges(7,4) * t225 + Ifges(7,5) * t237 + t227 * t263 - t228 * t269;
t246 = Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * t347;
t247 = Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * t347;
t377 = mrSges(6,1) * t211 - mrSges(6,2) * t212 + Ifges(6,5) * t239 + Ifges(6,6) * t238 + Ifges(6,3) * t346 + pkin(5) * t379 + pkin(9) * t384 + t364 * t194 + t360 * t195 + t274 * t246 - t273 * t247;
t373 = mrSges(5,1) * t220 - mrSges(5,2) * t221 + Ifges(5,5) * t284 + Ifges(5,6) * t283 - Ifges(5,3) * qJDD(2) + pkin(4) * t181 + t301 * t271 - t300 * t272 + t377;
t372 = -mrSges(4,1) * t266 + mrSges(4,3) * t262 + Ifges(4,4) * t325 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t326 - pkin(3) * t173 - t373;
t334 = mrSges(4,2) * t390 + qJD(2) * mrSges(4,3);
t378 = -m(4) * t266 + qJDD(2) * mrSges(4,1) + qJD(2) * t334 - t173;
t174 = -t357 * t178 + t358 * t179;
t332 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t391;
t381 = m(4) * t262 + qJDD(2) * mrSges(4,3) + qJD(2) * t332 + t323 * t390 + t174;
t306 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t362 - Ifges(4,5) * t366) * qJD(1);
t392 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t362 + Ifges(3,4) * t366) * qJD(1) + t306;
t398 = -((t302 - t305) * t362 + t366 * t392) * qJD(1) + mrSges(3,1) * t285 - mrSges(3,2) * t286 + Ifges(3,5) * t325 + Ifges(3,6) * t326 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t325 * mrSges(4,2) - t323 * t391 + t378) + qJ(3) * (mrSges(4,2) * t326 + t381) + t372;
t395 = mrSges(3,3) + mrSges(4,2);
t190 = t364 * t201 + t360 * t202;
t324 = (-mrSges(3,1) * t366 + mrSges(3,2) * t362) * qJD(1);
t331 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t391;
t169 = m(3) * t286 - qJDD(2) * mrSges(3,2) - qJD(2) * t331 + t324 * t390 + t326 * t395 + t381;
t333 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t390;
t170 = m(3) * t285 + qJDD(2) * mrSges(3,1) + qJD(2) * t333 - t395 * t325 + (-t323 - t324) * t391 + t378;
t386 = t366 * t169 - t170 * t362;
t383 = m(6) * t223 - t238 * mrSges(6,1) + t239 * mrSges(6,2) - t273 * t267 + t274 * t268 + t190;
t188 = -m(5) * t242 + t283 * mrSges(5,1) - t284 * mrSges(5,2) + t300 * t287 - t301 * t288 - t383;
t257 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t391 + t382;
t184 = m(4) * t257 - t326 * mrSges(4,1) - t325 * mrSges(4,3) - t332 * t391 - t334 * t390 + t188;
t303 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t362 + Ifges(3,6) * t366) * qJD(1);
t304 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t362 - Ifges(4,6) * t366) * qJD(1);
t245 = Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * t347;
t175 = mrSges(6,2) * t223 - mrSges(6,3) * t211 + Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t346 - pkin(9) * t190 - t194 * t360 + t195 * t364 + t245 * t273 - t246 * t347;
t374 = mrSges(7,1) * t205 - mrSges(7,2) * t206 + Ifges(7,5) * t226 + Ifges(7,6) * t225 + Ifges(7,3) * t237 + t228 * t264 - t229 * t263;
t176 = -mrSges(6,1) * t223 + mrSges(6,3) * t212 + Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t346 - pkin(5) * t190 - t245 * t274 + t247 * t347 - t374;
t270 = Ifges(5,5) * t301 + Ifges(5,6) * t300 - Ifges(5,3) * qJD(2);
t162 = -mrSges(5,1) * t242 + mrSges(5,3) * t221 + Ifges(5,4) * t284 + Ifges(5,2) * t283 - Ifges(5,6) * qJDD(2) - pkin(4) * t383 + pkin(8) * t385 - qJD(2) * t272 + t361 * t175 + t365 * t176 - t301 * t270;
t167 = mrSges(5,2) * t242 - mrSges(5,3) * t220 + Ifges(5,1) * t284 + Ifges(5,4) * t283 - Ifges(5,5) * qJDD(2) - pkin(8) * t181 + qJD(2) * t271 + t175 * t365 - t176 * t361 + t270 * t300;
t375 = -mrSges(4,1) * t257 + mrSges(4,2) * t262 - pkin(3) * t188 - qJ(4) * t174 - t358 * t162 - t357 * t167;
t159 = -mrSges(3,1) * t308 + mrSges(3,3) * t286 - pkin(2) * t184 + (Ifges(3,2) + Ifges(4,3)) * t326 + (Ifges(3,4) - Ifges(4,5)) * t325 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t392 * qJD(2) + (-t303 - t304) * t391 + t375;
t376 = mrSges(4,2) * t266 - mrSges(4,3) * t257 + Ifges(4,1) * t325 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t326 - qJ(4) * t173 + qJD(2) * t302 - t162 * t357 + t358 * t167 + t304 * t390;
t161 = mrSges(3,2) * t308 - mrSges(3,3) * t285 + Ifges(3,1) * t325 + Ifges(3,4) * t326 + Ifges(3,5) * qJDD(2) - qJ(3) * t184 - qJD(2) * t305 + t303 * t390 + t376;
t371 = -m(3) * t308 + t326 * mrSges(3,1) - t325 * mrSges(3,2) - t331 * t391 + t333 * t390 - t184;
t380 = mrSges(2,1) * t335 - mrSges(2,2) * t336 + Ifges(2,3) * qJDD(1) + pkin(1) * t371 + pkin(7) * t386 + t366 * t159 + t362 * t161;
t182 = m(2) * t335 + qJDD(1) * mrSges(2,1) - t369 * mrSges(2,2) + t371;
t165 = t169 * t362 + t170 * t366;
t163 = m(2) * t336 - mrSges(2,1) * t369 - qJDD(1) * mrSges(2,2) + t386;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t336 + t369 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t398;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t335 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t369 - pkin(7) * t165 - t159 * t362 + t161 * t366;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t367 * t156 - t363 * t157 - pkin(6) * (t163 * t363 + t182 * t367), t156, t161, t376, t167, t175, t195; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t363 * t156 + t367 * t157 + pkin(6) * (t163 * t367 - t182 * t363), t157, t159, t372 + (-t362 * t302 - t366 * t306) * qJD(1), t162, t176, t194; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t380, t380, t398, Ifges(4,5) * t325 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t326 - qJD(2) * t306 + t304 * t391 - t375, t373, t377, t374;];
m_new  = t1;
