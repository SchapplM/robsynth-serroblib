% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:52:55
% EndTime: 2019-05-06 18:53:17
% DurationCPUTime: 8.71s
% Computational Cost: add. (132108->386), mult. (271370->452), div. (0->0), fcn. (167398->8), ass. (0->137)
t340 = sin(qJ(4));
t344 = cos(qJ(4));
t345 = cos(qJ(2));
t380 = qJD(1) * t345;
t309 = t344 * qJD(2) - t340 * t380;
t341 = sin(qJ(2));
t379 = qJD(1) * qJD(2);
t373 = t341 * t379;
t314 = t345 * qJDD(1) - t373;
t267 = -t309 * qJD(4) - t340 * qJDD(2) - t344 * t314;
t308 = -t340 * qJD(2) - t344 * t380;
t268 = t308 * qJD(4) + t344 * qJDD(2) - t340 * t314;
t339 = sin(qJ(5));
t343 = cos(qJ(5));
t270 = t343 * t308 - t339 * t309;
t228 = t270 * qJD(5) + t339 * t267 + t343 * t268;
t271 = t339 * t308 + t343 * t309;
t246 = -mrSges(7,1) * t270 + mrSges(7,2) * t271;
t330 = t341 * qJD(1);
t322 = pkin(3) * t330 - qJD(2) * pkin(8);
t338 = t345 ^ 2;
t348 = qJD(1) ^ 2;
t372 = t345 * t379;
t313 = t341 * qJDD(1) + t372;
t342 = sin(qJ(1));
t346 = cos(qJ(1));
t323 = t342 * g(1) - t346 * g(2);
t367 = -qJDD(1) * pkin(1) - t323;
t393 = -2 * qJD(3);
t355 = pkin(2) * t373 + t330 * t393 + (-t313 - t372) * qJ(3) + t367;
t231 = -t322 * t330 + (-pkin(3) * t338 - pkin(7)) * t348 + (-pkin(2) - pkin(8)) * t314 + t355;
t324 = -t346 * g(1) - t342 * g(2);
t296 = -t348 * pkin(1) + qJDD(1) * pkin(7) + t324;
t275 = -t345 * g(3) - t341 * t296;
t310 = (-pkin(2) * t345 - qJ(3) * t341) * qJD(1);
t347 = qJD(2) ^ 2;
t253 = -qJDD(2) * pkin(2) - t347 * qJ(3) + t310 * t330 + qJDD(3) - t275;
t237 = (-t341 * t345 * t348 - qJDD(2)) * pkin(8) + (t313 - t372) * pkin(3) + t253;
t209 = -t340 * t231 + t344 * t237;
t307 = qJDD(4) + t313;
t327 = t330 + qJD(4);
t205 = (t308 * t327 - t268) * pkin(9) + (t308 * t309 + t307) * pkin(4) + t209;
t210 = t344 * t231 + t340 * t237;
t277 = t327 * pkin(4) - t309 * pkin(9);
t306 = t308 ^ 2;
t207 = -t306 * pkin(4) + t267 * pkin(9) - t327 * t277 + t210;
t198 = t343 * t205 - t339 * t207;
t298 = qJDD(5) + t307;
t325 = qJD(5) + t327;
t193 = -0.2e1 * qJD(6) * t271 + (t270 * t325 - t228) * qJ(6) + (t270 * t271 + t298) * pkin(5) + t198;
t254 = -mrSges(7,2) * t325 + t270 * mrSges(7,3);
t377 = m(7) * t193 + t298 * mrSges(7,1) + t325 * t254;
t190 = -t228 * mrSges(7,3) - t271 * t246 + t377;
t199 = t339 * t205 + t343 * t207;
t227 = -t271 * qJD(5) + t343 * t267 - t339 * t268;
t241 = Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t325;
t242 = Ifges(7,1) * t271 + Ifges(7,4) * t270 + Ifges(7,5) * t325;
t243 = Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t325;
t256 = pkin(5) * t325 - t271 * qJ(6);
t269 = t270 ^ 2;
t196 = -t269 * pkin(5) + t227 * qJ(6) + 0.2e1 * qJD(6) * t270 - t256 * t325 + t199;
t240 = Ifges(7,4) * t271 + Ifges(7,2) * t270 + Ifges(7,6) * t325;
t363 = -mrSges(7,1) * t193 + mrSges(7,2) * t196 - Ifges(7,5) * t228 - Ifges(7,6) * t227 - Ifges(7,3) * t298 - t271 * t240;
t397 = mrSges(6,1) * t198 - mrSges(6,2) * t199 + Ifges(6,5) * t228 + Ifges(6,6) * t227 + Ifges(6,3) * t298 + pkin(5) * t190 + t271 * t241 - t363 - (t242 + t243) * t270;
t276 = -t341 * g(3) + t345 * t296;
t251 = t347 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t393 - t310 * t380 - t276;
t236 = -t338 * t348 * pkin(8) + t314 * pkin(3) + qJD(2) * t322 - t251;
t212 = -t267 * pkin(4) - t306 * pkin(9) + t309 * t277 + t236;
t255 = -mrSges(6,2) * t325 + t270 * mrSges(6,3);
t258 = mrSges(6,1) * t325 - t271 * mrSges(6,3);
t202 = -t227 * pkin(5) - t269 * qJ(6) + t271 * t256 + qJDD(6) + t212;
t257 = mrSges(7,1) * t325 - t271 * mrSges(7,3);
t375 = m(7) * t202 + t228 * mrSges(7,2) + t271 * t257;
t396 = m(6) * t212 + t228 * mrSges(6,2) + t271 * t258 + t375 - (t254 + t255) * t270 - (mrSges(6,1) + mrSges(7,1)) * t227;
t273 = -t327 * mrSges(5,2) + t308 * mrSges(5,3);
t274 = t327 * mrSges(5,1) - t309 * mrSges(5,3);
t395 = -m(5) * t236 + t267 * mrSges(5,1) - t268 * mrSges(5,2) + t308 * t273 - t309 * t274 - t396;
t247 = -mrSges(6,1) * t270 + mrSges(6,2) * t271;
t181 = m(6) * t198 + t298 * mrSges(6,1) + t255 * t325 + (-t246 - t247) * t271 + (-mrSges(6,3) - mrSges(7,3)) * t228 + t377;
t376 = m(7) * t196 + t227 * mrSges(7,3) + t270 * t246;
t184 = m(6) * t199 + t227 * mrSges(6,3) + t270 * t247 + (-t257 - t258) * t325 + (-mrSges(6,2) - mrSges(7,2)) * t298 + t376;
t179 = t343 * t181 + t339 * t184;
t260 = Ifges(5,4) * t309 + Ifges(5,2) * t308 + Ifges(5,6) * t327;
t261 = Ifges(5,1) * t309 + Ifges(5,4) * t308 + Ifges(5,5) * t327;
t394 = mrSges(5,1) * t209 - mrSges(5,2) * t210 + Ifges(5,5) * t268 + Ifges(5,6) * t267 + Ifges(5,3) * t307 + pkin(4) * t179 + t309 * t260 - t308 * t261 + t397;
t272 = -mrSges(5,1) * t308 + mrSges(5,2) * t309;
t175 = m(5) * t209 + t307 * mrSges(5,1) - t268 * mrSges(5,3) - t309 * t272 + t327 * t273 + t179;
t368 = -t339 * t181 + t343 * t184;
t176 = m(5) * t210 - t307 * mrSges(5,2) + t267 * mrSges(5,3) + t308 * t272 - t327 * t274 + t368;
t170 = t344 * t175 + t340 * t176;
t388 = t348 * pkin(7);
t248 = -t314 * pkin(2) + t355 - t388;
t392 = mrSges(4,1) * t253 - mrSges(4,3) * t248 + pkin(3) * t170 + t394;
t291 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t341 + Ifges(3,4) * t345) * qJD(1);
t311 = (mrSges(4,2) * t345 - mrSges(4,3) * t341) * qJD(1);
t320 = -mrSges(4,1) * t380 - qJD(2) * mrSges(4,3);
t321 = mrSges(4,1) * t330 + qJD(2) * mrSges(4,2);
t351 = -m(4) * t251 + qJDD(2) * mrSges(4,3) + qJD(2) * t321 + t311 * t380 - t395;
t238 = Ifges(7,5) * t271 + Ifges(7,6) * t270 + Ifges(7,3) * t325;
t239 = Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t325;
t364 = -mrSges(7,1) * t202 + mrSges(7,3) * t196 + Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t298 + t325 * t242;
t172 = Ifges(6,4) * t228 + Ifges(6,2) * t227 + Ifges(6,6) * t298 + t325 * t243 - mrSges(6,1) * t212 + mrSges(6,3) * t199 - pkin(5) * (-t227 * mrSges(7,1) - t270 * t254 + t375) + qJ(6) * (-t298 * mrSges(7,2) - t257 * t325 + t376) + (-t239 - t238) * t271 + t364;
t362 = mrSges(7,2) * t202 - mrSges(7,3) * t193 + Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t298 + t270 * t238;
t177 = mrSges(6,2) * t212 - mrSges(6,3) * t198 + Ifges(6,1) * t228 + Ifges(6,4) * t227 + Ifges(6,5) * t298 - qJ(6) * t190 + t270 * t239 + (-t240 - t241) * t325 + t362;
t259 = Ifges(5,5) * t309 + Ifges(5,6) * t308 + Ifges(5,3) * t327;
t159 = -mrSges(5,1) * t236 + mrSges(5,3) * t210 + Ifges(5,4) * t268 + Ifges(5,2) * t267 + Ifges(5,6) * t307 - pkin(4) * t396 + pkin(9) * t368 + t343 * t172 + t339 * t177 - t309 * t259 + t327 * t261;
t161 = mrSges(5,2) * t236 - mrSges(5,3) * t209 + Ifges(5,1) * t268 + Ifges(5,4) * t267 + Ifges(5,5) * t307 - pkin(9) * t179 - t339 * t172 + t343 * t177 + t308 * t259 - t327 * t260;
t293 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t341 - Ifges(4,6) * t345) * qJD(1);
t358 = -mrSges(4,2) * t253 + mrSges(4,3) * t251 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t313 + Ifges(4,5) * t314 + pkin(8) * t170 + t340 * t159 - t344 * t161 - t293 * t380;
t360 = -m(4) * t253 - t313 * mrSges(4,1) - t170;
t292 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t341 - Ifges(4,3) * t345) * qJD(1);
t381 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t341 + Ifges(3,2) * t345) * qJD(1) - t292;
t391 = (-t345 * t291 + t341 * t381) * qJD(1) + mrSges(3,1) * t275 - mrSges(3,2) * t276 + Ifges(3,5) * t313 + Ifges(3,6) * t314 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t320 - t311 * t330 + t360) + qJ(3) * (t314 * mrSges(4,1) + t351) - t358;
t386 = Ifges(3,4) + Ifges(4,6);
t171 = -t340 * t175 + t344 * t176;
t294 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t341 - Ifges(4,5) * t345) * qJD(1);
t382 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t341 + Ifges(3,6) * t345) * qJD(1) + t294;
t312 = (-mrSges(3,1) * t345 + mrSges(3,2) * t341) * qJD(1);
t319 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t380;
t167 = m(3) * t275 - t313 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t319 - t320) * qJD(2) + (-t311 - t312) * t330 + t360;
t318 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t330;
t186 = -qJDD(2) * mrSges(3,2) + t351 - qJD(2) * t318 + m(3) * t276 + t312 * t380 + (mrSges(4,1) + mrSges(3,3)) * t314;
t369 = -t341 * t167 + t345 * t186;
t365 = -m(4) * t248 - t314 * mrSges(4,2) + t321 * t330 - t171;
t168 = -t313 * mrSges(4,3) + t320 * t380 - t365;
t295 = t367 - t388;
t356 = -mrSges(4,1) * t251 + mrSges(4,2) * t248 - pkin(3) * t395 - pkin(8) * t171 - t344 * t159 - t340 * t161;
t156 = -mrSges(3,1) * t295 + mrSges(3,3) * t276 - pkin(2) * t168 + (Ifges(3,2) + Ifges(4,3)) * t314 + t386 * t313 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t291 - t293) * qJD(2) - t382 * t330 + t356;
t158 = t382 * t380 + mrSges(3,2) * t295 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - mrSges(3,3) * t275 + (Ifges(3,1) + Ifges(4,2)) * t313 + t386 * t314 - t381 * qJD(2) - qJ(3) * t168 + t392;
t353 = -m(3) * t295 + t319 * t380 + t314 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t313 + (-t318 * t341 - t320 * t345) * qJD(1) + t365;
t359 = mrSges(2,1) * t323 - mrSges(2,2) * t324 + Ifges(2,3) * qJDD(1) + pkin(1) * t353 + pkin(7) * t369 + t345 * t156 + t341 * t158;
t165 = m(2) * t323 + qJDD(1) * mrSges(2,1) - t348 * mrSges(2,2) + t353;
t164 = t345 * t167 + t341 * t186;
t162 = m(2) * t324 - t348 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t369;
t154 = mrSges(2,1) * g(3) + mrSges(2,3) * t324 + t348 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t164 - t391;
t153 = -mrSges(2,2) * g(3) - mrSges(2,3) * t323 + Ifges(2,5) * qJDD(1) - t348 * Ifges(2,6) - pkin(7) * t164 - t341 * t156 + t345 * t158;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t346 * t153 - t342 * t154 - pkin(6) * (t342 * t162 + t346 * t165), t153, t158, -t292 * t330 - t358, t161, t177, -t240 * t325 + t362; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t153 + t346 * t154 + pkin(6) * (t346 * t162 - t342 * t165), t154, t156, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t313 - Ifges(4,6) * t314 - qJD(2) * t292 - t294 * t380 - t392, t159, t172, -t271 * t238 + t364; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, t391, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t313 - Ifges(4,3) * t314 + qJD(2) * t293 + t294 * t330 - t356, t394, t397, -t270 * t242 - t363;];
m_new  = t1;
