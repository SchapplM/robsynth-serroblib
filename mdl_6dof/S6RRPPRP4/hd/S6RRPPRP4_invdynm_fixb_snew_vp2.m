% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-05-06 09:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:18:47
% EndTime: 2019-05-06 09:19:03
% DurationCPUTime: 6.88s
% Computational Cost: add. (99829->377), mult. (217824->443), div. (0->0), fcn. (140690->8), ass. (0->133)
t393 = -2 * qJD(3);
t392 = -2 * qJD(4);
t349 = sin(qJ(1));
t351 = cos(qJ(1));
t333 = -t351 * g(1) - t349 * g(2);
t353 = qJD(1) ^ 2;
t308 = -t353 * pkin(1) + qJDD(1) * pkin(7) + t333;
t348 = sin(qJ(2));
t350 = cos(qJ(2));
t279 = -t350 * g(3) - t348 * t308;
t325 = (-pkin(2) * t350 - qJ(3) * t348) * qJD(1);
t352 = qJD(2) ^ 2;
t379 = qJD(1) * t348;
t256 = -qJDD(2) * pkin(2) - t352 * qJ(3) + t325 * t379 + qJDD(3) - t279;
t375 = qJD(1) * qJD(2);
t372 = t350 * t375;
t327 = t348 * qJDD(1) + t372;
t346 = sin(pkin(9));
t386 = cos(pkin(9));
t293 = -t386 * qJDD(2) + t346 * t327;
t294 = t346 * qJDD(2) + t386 * t327;
t318 = t346 * qJD(2) + t386 * t379;
t317 = -t386 * qJD(2) + t346 * t379;
t376 = t350 * qJD(1);
t373 = t317 * t376;
t216 = t256 - (t294 + t373) * qJ(4) - (t318 * t376 - t293) * pkin(3) + t318 * t392;
t332 = t349 * g(1) - t351 * g(2);
t307 = -qJDD(1) * pkin(1) - t353 * pkin(7) - t332;
t337 = t348 * t375;
t328 = t350 * qJDD(1) - t337;
t250 = (-t327 - t372) * qJ(3) + (-t328 + t337) * pkin(2) + t307;
t280 = -t348 * g(3) + t350 * t308;
t257 = -t352 * pkin(2) + qJDD(2) * qJ(3) + t325 * t376 + t280;
t217 = t386 * t250 - t346 * t257 + t318 * t393;
t218 = t346 * t250 + t386 * t257 + t317 * t393;
t269 = Ifges(4,1) * t318 - Ifges(4,4) * t317 - Ifges(4,5) * t376;
t276 = t317 * mrSges(5,1) - t318 * mrSges(5,3);
t289 = -t317 * mrSges(5,2) - mrSges(5,3) * t376;
t292 = mrSges(5,1) * t376 + t318 * mrSges(5,2);
t275 = t317 * pkin(3) - t318 * qJ(4);
t385 = t350 ^ 2 * t353;
t215 = t328 * pkin(3) - qJ(4) * t385 + t318 * t275 + qJDD(4) - t217;
t206 = (-t294 + t373) * pkin(8) + (t317 * t318 + t328) * pkin(4) + t215;
t213 = -pkin(3) * t385 - t328 * qJ(4) - t317 * t275 + t376 * t392 + t218;
t295 = pkin(4) * t376 - t318 * pkin(8);
t315 = t317 ^ 2;
t208 = -t315 * pkin(4) + t293 * pkin(8) - t295 * t376 + t213;
t347 = sin(qJ(5));
t389 = cos(qJ(5));
t202 = t347 * t206 + t389 * t208;
t274 = t347 * t317 + t389 * t318;
t234 = t274 * qJD(5) - t389 * t293 + t347 * t294;
t335 = qJD(5) + t376;
t259 = t335 * mrSges(6,1) - t274 * mrSges(6,3);
t273 = -t389 * t317 + t347 * t318;
t324 = qJDD(5) + t328;
t246 = t273 * pkin(5) - t274 * qJ(6);
t334 = t335 ^ 2;
t195 = -t334 * pkin(5) + t324 * qJ(6) + 0.2e1 * qJD(6) * t335 - t273 * t246 + t202;
t260 = -t335 * mrSges(7,1) + t274 * mrSges(7,2);
t374 = m(7) * t195 + t324 * mrSges(7,3) + t335 * t260;
t247 = t273 * mrSges(7,1) - t274 * mrSges(7,3);
t383 = -t273 * mrSges(6,1) - t274 * mrSges(6,2) - t247;
t387 = -mrSges(6,3) - mrSges(7,2);
t184 = m(6) * t202 - t324 * mrSges(6,2) + t387 * t234 - t335 * t259 + t383 * t273 + t374;
t201 = t389 * t206 - t347 * t208;
t235 = -t273 * qJD(5) + t347 * t293 + t389 * t294;
t258 = -t335 * mrSges(6,2) - t273 * mrSges(6,3);
t198 = -t324 * pkin(5) - t334 * qJ(6) + t274 * t246 + qJDD(6) - t201;
t261 = -t273 * mrSges(7,2) + t335 * mrSges(7,3);
t370 = -m(7) * t198 + t324 * mrSges(7,1) + t335 * t261;
t185 = m(6) * t201 + t324 * mrSges(6,1) + t387 * t235 + t335 * t258 + t383 * t274 + t370;
t179 = t347 * t184 + t389 * t185;
t268 = Ifges(5,1) * t318 - Ifges(5,4) * t376 + Ifges(5,5) * t317;
t239 = Ifges(6,4) * t274 - Ifges(6,2) * t273 + Ifges(6,6) * t335;
t241 = Ifges(6,1) * t274 - Ifges(6,4) * t273 + Ifges(6,5) * t335;
t236 = Ifges(7,5) * t274 + Ifges(7,6) * t335 + Ifges(7,3) * t273;
t240 = Ifges(7,1) * t274 + Ifges(7,4) * t335 + Ifges(7,5) * t273;
t365 = -mrSges(7,1) * t198 + mrSges(7,3) * t195 + Ifges(7,4) * t235 + Ifges(7,2) * t324 + Ifges(7,6) * t234 - t274 * t236 + t273 * t240;
t358 = Ifges(6,5) * t235 - Ifges(6,6) * t234 + Ifges(6,3) * t324 + t274 * t239 + t273 * t241 + mrSges(6,1) * t201 - mrSges(6,2) * t202 + pkin(5) * (-t235 * mrSges(7,2) - t274 * t247 + t370) + qJ(6) * (-t234 * mrSges(7,2) - t273 * t247 + t374) + t365;
t355 = mrSges(5,1) * t215 - mrSges(5,3) * t213 - Ifges(5,4) * t294 + Ifges(5,2) * t328 - Ifges(5,6) * t293 + pkin(4) * t179 - t317 * t268 + t358;
t362 = -m(5) * t215 - t328 * mrSges(5,1) - t179;
t180 = t389 * t184 - t347 * t185;
t366 = m(5) * t213 - t328 * mrSges(5,3) + t180;
t264 = Ifges(5,5) * t318 - Ifges(5,6) * t376 + Ifges(5,3) * t317;
t382 = -Ifges(4,4) * t318 + Ifges(4,2) * t317 + Ifges(4,6) * t376 + t264;
t391 = t382 * t318 - mrSges(4,1) * t217 + mrSges(4,2) * t218 - Ifges(4,5) * t294 + Ifges(4,6) * t293 - pkin(3) * (-t294 * mrSges(5,2) - t318 * t276 - t289 * t376 + t362) - qJ(4) * (-t293 * mrSges(5,2) - t317 * t276 - t292 * t376 + t366) - t317 * t269 + t355;
t210 = -t293 * pkin(4) - t315 * pkin(8) + t318 * t295 - t216;
t204 = t210 - 0.2e1 * qJD(6) * t274 + (t273 * t335 - t235) * qJ(6) + (t274 * t335 + t234) * pkin(5);
t192 = m(7) * t204 + t234 * mrSges(7,1) - t235 * mrSges(7,3) - t274 * t260 + t273 * t261;
t187 = -m(6) * t210 - t234 * mrSges(6,1) - t235 * mrSges(6,2) - t273 * t258 - t274 * t259 - t192;
t186 = m(5) * t216 + t293 * mrSges(5,1) - t294 * mrSges(5,3) + t317 * t289 - t318 * t292 + t187;
t369 = -mrSges(7,1) * t204 + mrSges(7,2) * t195;
t238 = Ifges(7,4) * t274 + Ifges(7,2) * t335 + Ifges(7,6) * t273;
t384 = -Ifges(6,5) * t274 + Ifges(6,6) * t273 - Ifges(6,3) * t335 - t238;
t176 = -mrSges(6,1) * t210 + mrSges(6,3) * t202 - pkin(5) * t192 + (t240 + t241) * t335 + (Ifges(6,6) - Ifges(7,6)) * t324 + t384 * t274 + (Ifges(6,4) - Ifges(7,5)) * t235 + (-Ifges(6,2) - Ifges(7,3)) * t234 + t369;
t364 = mrSges(7,2) * t198 - mrSges(7,3) * t204 + Ifges(7,1) * t235 + Ifges(7,4) * t324 + Ifges(7,5) * t234 + t335 * t236;
t178 = mrSges(6,2) * t210 - mrSges(6,3) * t201 + Ifges(6,1) * t235 - Ifges(6,4) * t234 + Ifges(6,5) * t324 - qJ(6) * t192 - t335 * t239 + t384 * t273 + t364;
t359 = mrSges(5,1) * t216 - mrSges(5,2) * t213 + pkin(4) * t187 + pkin(8) * t180 + t389 * t176 + t347 * t178;
t266 = Ifges(5,4) * t318 - Ifges(5,2) * t376 + Ifges(5,6) * t317;
t381 = -Ifges(4,5) * t318 + Ifges(4,6) * t317 + Ifges(4,3) * t376 - t266;
t160 = -mrSges(4,1) * t256 + mrSges(4,3) * t218 - pkin(3) * t186 + (-Ifges(4,6) + Ifges(5,6)) * t328 + t381 * t318 + (Ifges(4,4) - Ifges(5,5)) * t294 + (-Ifges(4,2) - Ifges(5,3)) * t293 + (-t268 - t269) * t376 - t359;
t361 = mrSges(5,2) * t215 - mrSges(5,3) * t216 + Ifges(5,1) * t294 - Ifges(5,4) * t328 + Ifges(5,5) * t293 - pkin(8) * t179 - t347 * t176 + t389 * t178;
t161 = mrSges(4,2) * t256 - mrSges(4,3) * t217 + Ifges(4,1) * t294 - Ifges(4,4) * t293 - Ifges(4,5) * t328 - qJ(4) * t186 + t381 * t317 - t382 * t376 + t361;
t291 = -mrSges(4,1) * t376 - t318 * mrSges(4,3);
t380 = -t317 * mrSges(4,1) - t318 * mrSges(4,2) - t276;
t388 = -mrSges(4,3) - mrSges(5,2);
t172 = m(4) * t218 + t328 * mrSges(4,2) + t380 * t317 + t388 * t293 + (t291 - t292) * t376 + t366;
t290 = mrSges(4,2) * t376 - t317 * mrSges(4,3);
t173 = m(4) * t217 - t328 * mrSges(4,1) + t380 * t318 + t388 * t294 + (-t289 - t290) * t376 + t362;
t170 = t386 * t172 - t346 * t173;
t182 = -m(4) * t256 - t293 * mrSges(4,1) - t294 * mrSges(4,2) - t317 * t290 - t318 * t291 - t186;
t305 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t348 + Ifges(3,2) * t350) * qJD(1);
t306 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t348 + Ifges(3,4) * t350) * qJD(1);
t390 = mrSges(3,1) * t279 - mrSges(3,2) * t280 + Ifges(3,5) * t327 + Ifges(3,6) * t328 + Ifges(3,3) * qJDD(2) + pkin(2) * t182 + qJ(3) * t170 + (t305 * t348 - t306 * t350) * qJD(1) + t386 * t160 + t346 * t161;
t326 = (-mrSges(3,1) * t350 + mrSges(3,2) * t348) * qJD(1);
t330 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t379;
t168 = m(3) * t280 - qJDD(2) * mrSges(3,2) + t328 * mrSges(3,3) - qJD(2) * t330 + t326 * t376 + t170;
t331 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t376;
t181 = m(3) * t279 + qJDD(2) * mrSges(3,1) - t327 * mrSges(3,3) + qJD(2) * t331 - t326 * t379 + t182;
t371 = t350 * t168 - t348 * t181;
t169 = t346 * t172 + t386 * t173;
t304 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t348 + Ifges(3,6) * t350) * qJD(1);
t157 = mrSges(3,2) * t307 - mrSges(3,3) * t279 + Ifges(3,1) * t327 + Ifges(3,4) * t328 + Ifges(3,5) * qJDD(2) - qJ(3) * t169 - qJD(2) * t305 - t346 * t160 + t386 * t161 + t304 * t376;
t159 = Ifges(3,6) * qJDD(2) - t304 * t379 + Ifges(3,4) * t327 + qJD(2) * t306 - mrSges(3,1) * t307 + mrSges(3,3) * t280 - pkin(2) * t169 + (Ifges(4,3) + Ifges(3,2)) * t328 + t391;
t357 = -m(3) * t307 + t328 * mrSges(3,1) - t327 * mrSges(3,2) - t330 * t379 + t331 * t376 - t169;
t363 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t357 + pkin(7) * t371 + t348 * t157 + t350 * t159;
t165 = m(2) * t332 + qJDD(1) * mrSges(2,1) - t353 * mrSges(2,2) + t357;
t164 = t348 * t168 + t350 * t181;
t162 = m(2) * t333 - t353 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t371;
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + t353 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t164 - t390;
t154 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - t353 * Ifges(2,6) - pkin(7) * t164 + t350 * t157 - t348 * t159;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t351 * t154 - t349 * t155 - pkin(6) * (t349 * t162 + t351 * t165), t154, t157, t161, -t264 * t376 - t317 * t266 + t361, t178, -t273 * t238 + t364; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t349 * t154 + t351 * t155 + pkin(6) * (t351 * t162 - t349 * t165), t155, t159, t160, -t318 * t264 - t355, t176, t365; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t363, t363, t390, -Ifges(4,3) * t328 - t391, Ifges(5,5) * t294 - Ifges(5,6) * t328 + Ifges(5,3) * t293 + t318 * t266 + t268 * t376 + t359, t358, Ifges(7,5) * t235 + Ifges(7,6) * t324 + Ifges(7,3) * t234 + t274 * t238 - t335 * t240 - t369;];
m_new  = t1;
