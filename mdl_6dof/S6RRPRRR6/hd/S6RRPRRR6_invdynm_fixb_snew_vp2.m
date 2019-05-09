% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR6
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
% Datum: 2019-05-06 22:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:58:23
% EndTime: 2019-05-06 21:58:56
% DurationCPUTime: 14.86s
% Computational Cost: add. (252531->381), mult. (545632->471), div. (0->0), fcn. (367916->10), ass. (0->145)
t362 = sin(qJ(1));
t367 = cos(qJ(1));
t334 = -g(1) * t367 - g(2) * t362;
t369 = qJD(1) ^ 2;
t307 = -pkin(1) * t369 + qJDD(1) * pkin(7) + t334;
t361 = sin(qJ(2));
t366 = cos(qJ(2));
t285 = -t366 * g(3) - t361 * t307;
t286 = -g(3) * t361 + t366 * t307;
t298 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t361 - Ifges(4,3) * t366) * qJD(1);
t301 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t361 + Ifges(3,2) * t366) * qJD(1);
t321 = (-mrSges(4,1) * t366 - mrSges(4,3) * t361) * qJD(1);
t389 = qJD(1) * qJD(2);
t388 = t366 * t389;
t323 = qJDD(1) * t361 + t388;
t387 = t361 * t389;
t324 = qJDD(1) * t366 - t387;
t320 = (-pkin(2) * t366 - qJ(3) * t361) * qJD(1);
t368 = qJD(2) ^ 2;
t390 = qJD(1) * t366;
t396 = 2 * qJD(3);
t262 = -pkin(2) * t368 + qJDD(2) * qJ(3) + qJD(2) * t396 + t320 * t390 + t286;
t391 = qJD(1) * t361;
t332 = -qJD(2) * pkin(3) - pkin(8) * t391;
t394 = t366 ^ 2 * t369;
t253 = -pkin(3) * t394 - pkin(8) * t324 + qJD(2) * t332 + t262;
t270 = -qJDD(2) * pkin(2) - t368 * qJ(3) + t320 * t391 + qJDD(3) - t285;
t255 = (-t323 + t388) * pkin(8) + (-t361 * t366 * t369 - qJDD(2)) * pkin(3) + t270;
t360 = sin(qJ(4));
t365 = cos(qJ(4));
t225 = -t360 * t253 + t365 * t255;
t304 = (-t360 * t361 - t365 * t366) * qJD(1);
t272 = qJD(4) * t304 + t323 * t365 - t324 * t360;
t305 = (-t360 * t366 + t361 * t365) * qJD(1);
t347 = -qJDD(2) + qJDD(4);
t348 = -qJD(2) + qJD(4);
t215 = (t304 * t348 - t272) * pkin(9) + (t304 * t305 + t347) * pkin(4) + t225;
t226 = t365 * t253 + t360 * t255;
t271 = -qJD(4) * t305 - t323 * t360 - t324 * t365;
t289 = pkin(4) * t348 - pkin(9) * t305;
t297 = t304 ^ 2;
t217 = -pkin(4) * t297 + pkin(9) * t271 - t289 * t348 + t226;
t359 = sin(qJ(5));
t364 = cos(qJ(5));
t213 = t359 * t215 + t364 * t217;
t283 = t304 * t359 + t305 * t364;
t236 = -qJD(5) * t283 + t271 * t364 - t272 * t359;
t282 = t304 * t364 - t305 * t359;
t254 = -mrSges(6,1) * t282 + mrSges(6,2) * t283;
t342 = qJD(5) + t348;
t274 = mrSges(6,1) * t342 - mrSges(6,3) * t283;
t341 = qJDD(5) + t347;
t256 = -pkin(5) * t282 - pkin(10) * t283;
t340 = t342 ^ 2;
t208 = -pkin(5) * t340 + pkin(10) * t341 + t256 * t282 + t213;
t333 = t362 * g(1) - t367 * g(2);
t306 = -qJDD(1) * pkin(1) - t369 * pkin(7) - t333;
t382 = -t324 * pkin(2) + t306 + (-t323 - t388) * qJ(3);
t242 = -pkin(2) * t387 + t324 * pkin(3) - pkin(8) * t394 - t382 + (t332 + t396) * t391;
t222 = -t271 * pkin(4) - t297 * pkin(9) + t305 * t289 + t242;
t237 = qJD(5) * t282 + t271 * t359 + t272 * t364;
t209 = (t283 * t342 - t236) * pkin(5) + (-t282 * t342 - t237) * pkin(10) + t222;
t358 = sin(qJ(6));
t363 = cos(qJ(6));
t205 = -t208 * t358 + t209 * t363;
t263 = -t283 * t358 + t342 * t363;
t220 = qJD(6) * t263 + t237 * t363 + t341 * t358;
t234 = qJDD(6) - t236;
t264 = t283 * t363 + t342 * t358;
t240 = -mrSges(7,1) * t263 + mrSges(7,2) * t264;
t278 = qJD(6) - t282;
t243 = -mrSges(7,2) * t278 + mrSges(7,3) * t263;
t201 = m(7) * t205 + mrSges(7,1) * t234 - mrSges(7,3) * t220 - t240 * t264 + t243 * t278;
t206 = t208 * t363 + t209 * t358;
t219 = -qJD(6) * t264 - t237 * t358 + t341 * t363;
t244 = mrSges(7,1) * t278 - mrSges(7,3) * t264;
t202 = m(7) * t206 - mrSges(7,2) * t234 + mrSges(7,3) * t219 + t240 * t263 - t244 * t278;
t384 = -t201 * t358 + t363 * t202;
t188 = m(6) * t213 - mrSges(6,2) * t341 + mrSges(6,3) * t236 + t254 * t282 - t274 * t342 + t384;
t212 = t215 * t364 - t217 * t359;
t273 = -mrSges(6,2) * t342 + mrSges(6,3) * t282;
t207 = -pkin(5) * t341 - pkin(10) * t340 + t256 * t283 - t212;
t379 = -m(7) * t207 + t219 * mrSges(7,1) - mrSges(7,2) * t220 + t263 * t243 - t244 * t264;
t197 = m(6) * t212 + mrSges(6,1) * t341 - mrSges(6,3) * t237 - t254 * t283 + t273 * t342 + t379;
t181 = t359 * t188 + t364 * t197;
t284 = -mrSges(5,1) * t304 + mrSges(5,2) * t305;
t287 = -mrSges(5,2) * t348 + mrSges(5,3) * t304;
t178 = m(5) * t225 + mrSges(5,1) * t347 - mrSges(5,3) * t272 - t284 * t305 + t287 * t348 + t181;
t288 = mrSges(5,1) * t348 - mrSges(5,3) * t305;
t385 = t364 * t188 - t197 * t359;
t179 = m(5) * t226 - mrSges(5,2) * t347 + mrSges(5,3) * t271 + t284 * t304 - t288 * t348 + t385;
t173 = t365 * t178 + t360 * t179;
t276 = Ifges(5,4) * t305 + Ifges(5,2) * t304 + Ifges(5,6) * t348;
t277 = Ifges(5,1) * t305 + Ifges(5,4) * t304 + Ifges(5,5) * t348;
t227 = Ifges(7,5) * t264 + Ifges(7,6) * t263 + Ifges(7,3) * t278;
t229 = Ifges(7,1) * t264 + Ifges(7,4) * t263 + Ifges(7,5) * t278;
t194 = -mrSges(7,1) * t207 + mrSges(7,3) * t206 + Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t234 - t227 * t264 + t229 * t278;
t228 = Ifges(7,4) * t264 + Ifges(7,2) * t263 + Ifges(7,6) * t278;
t195 = mrSges(7,2) * t207 - mrSges(7,3) * t205 + Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t234 + t227 * t263 - t228 * t278;
t246 = Ifges(6,4) * t283 + Ifges(6,2) * t282 + Ifges(6,6) * t342;
t247 = Ifges(6,1) * t283 + Ifges(6,4) * t282 + Ifges(6,5) * t342;
t377 = mrSges(6,1) * t212 - mrSges(6,2) * t213 + Ifges(6,5) * t237 + Ifges(6,6) * t236 + Ifges(6,3) * t341 + pkin(5) * t379 + pkin(10) * t384 + t363 * t194 + t358 * t195 + t283 * t246 - t282 * t247;
t373 = mrSges(5,1) * t225 - mrSges(5,2) * t226 + Ifges(5,5) * t272 + Ifges(5,6) * t271 + Ifges(5,3) * t347 + pkin(4) * t181 + t305 * t276 - t304 * t277 + t377;
t372 = -mrSges(4,1) * t270 + mrSges(4,3) * t262 + Ifges(4,4) * t323 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t324 - pkin(3) * t173 - t373;
t331 = mrSges(4,2) * t390 + qJD(2) * mrSges(4,3);
t378 = -m(4) * t270 + qJDD(2) * mrSges(4,1) + qJD(2) * t331 - t173;
t174 = -t360 * t178 + t365 * t179;
t329 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t391;
t381 = m(4) * t262 + qJDD(2) * mrSges(4,3) + qJD(2) * t329 + t321 * t390 + t174;
t302 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t361 - Ifges(4,5) * t366) * qJD(1);
t392 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t361 + Ifges(3,4) * t366) * qJD(1) + t302;
t398 = -(t392 * t366 + (t298 - t301) * t361) * qJD(1) + mrSges(3,1) * t285 - mrSges(3,2) * t286 + Ifges(3,5) * t323 + Ifges(3,6) * t324 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t323 * mrSges(4,2) - t321 * t391 + t378) + qJ(3) * (mrSges(4,2) * t324 + t381) + t372;
t395 = mrSges(3,3) + mrSges(4,2);
t190 = t363 * t201 + t358 * t202;
t322 = (-mrSges(3,1) * t366 + mrSges(3,2) * t361) * qJD(1);
t328 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t391;
t169 = m(3) * t286 - qJDD(2) * mrSges(3,2) - qJD(2) * t328 + t322 * t390 + t395 * t324 + t381;
t330 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t390;
t170 = m(3) * t285 + qJDD(2) * mrSges(3,1) + qJD(2) * t330 - t395 * t323 + (-t321 - t322) * t391 + t378;
t386 = t366 * t169 - t170 * t361;
t383 = m(6) * t222 - t236 * mrSges(6,1) + t237 * mrSges(6,2) - t282 * t273 + t283 * t274 + t190;
t185 = -m(5) * t242 + t271 * mrSges(5,1) - t272 * mrSges(5,2) + t304 * t287 - t305 * t288 - t383;
t257 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t391 + t382;
t184 = m(4) * t257 - t324 * mrSges(4,1) - t323 * mrSges(4,3) - t329 * t391 - t331 * t390 + t185;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t361 + Ifges(3,6) * t366) * qJD(1);
t300 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t361 - Ifges(4,6) * t366) * qJD(1);
t245 = Ifges(6,5) * t283 + Ifges(6,6) * t282 + Ifges(6,3) * t342;
t175 = mrSges(6,2) * t222 - mrSges(6,3) * t212 + Ifges(6,1) * t237 + Ifges(6,4) * t236 + Ifges(6,5) * t341 - pkin(10) * t190 - t194 * t358 + t195 * t363 + t245 * t282 - t246 * t342;
t374 = mrSges(7,1) * t205 - mrSges(7,2) * t206 + Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t234 + t228 * t264 - t229 * t263;
t176 = -mrSges(6,1) * t222 + mrSges(6,3) * t213 + Ifges(6,4) * t237 + Ifges(6,2) * t236 + Ifges(6,6) * t341 - pkin(5) * t190 - t245 * t283 + t247 * t342 - t374;
t275 = Ifges(5,5) * t305 + Ifges(5,6) * t304 + Ifges(5,3) * t348;
t162 = -mrSges(5,1) * t242 + mrSges(5,3) * t226 + Ifges(5,4) * t272 + Ifges(5,2) * t271 + Ifges(5,6) * t347 - pkin(4) * t383 + pkin(9) * t385 + t359 * t175 + t364 * t176 - t305 * t275 + t348 * t277;
t167 = mrSges(5,2) * t242 - mrSges(5,3) * t225 + Ifges(5,1) * t272 + Ifges(5,4) * t271 + Ifges(5,5) * t347 - pkin(9) * t181 + t175 * t364 - t176 * t359 + t275 * t304 - t276 * t348;
t375 = -mrSges(4,1) * t257 + mrSges(4,2) * t262 - pkin(3) * t185 - pkin(8) * t174 - t365 * t162 - t360 * t167;
t159 = -mrSges(3,1) * t306 + mrSges(3,3) * t286 - pkin(2) * t184 + (Ifges(3,2) + Ifges(4,3)) * t324 + (Ifges(3,4) - Ifges(4,5)) * t323 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t392 * qJD(2) + (-t299 - t300) * t391 + t375;
t376 = mrSges(4,2) * t270 - mrSges(4,3) * t257 + Ifges(4,1) * t323 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t324 - pkin(8) * t173 + qJD(2) * t298 - t162 * t360 + t365 * t167 + t300 * t390;
t161 = mrSges(3,2) * t306 - mrSges(3,3) * t285 + Ifges(3,1) * t323 + Ifges(3,4) * t324 + Ifges(3,5) * qJDD(2) - qJ(3) * t184 - qJD(2) * t301 + t299 * t390 + t376;
t371 = -m(3) * t306 + t324 * mrSges(3,1) - t323 * mrSges(3,2) - t328 * t391 + t330 * t390 - t184;
t380 = mrSges(2,1) * t333 - mrSges(2,2) * t334 + Ifges(2,3) * qJDD(1) + pkin(1) * t371 + pkin(7) * t386 + t366 * t159 + t361 * t161;
t182 = m(2) * t333 + qJDD(1) * mrSges(2,1) - t369 * mrSges(2,2) + t371;
t165 = t169 * t361 + t170 * t366;
t163 = m(2) * t334 - mrSges(2,1) * t369 - qJDD(1) * mrSges(2,2) + t386;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t334 + t369 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t398;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t333 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t369 - pkin(7) * t165 - t159 * t361 + t161 * t366;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t367 * t156 - t362 * t157 - pkin(6) * (t163 * t362 + t182 * t367), t156, t161, t376, t167, t175, t195; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t362 * t156 + t367 * t157 + pkin(6) * (t163 * t367 - t182 * t362), t157, t159, t372 + (-t361 * t298 - t366 * t302) * qJD(1), t162, t176, t194; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t380, t380, t398, Ifges(4,5) * t323 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t324 - qJD(2) * t302 + t300 * t391 - t375, t373, t377, t374;];
m_new  = t1;
