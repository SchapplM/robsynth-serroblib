% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 09:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:31:08
% EndTime: 2019-05-06 09:31:35
% DurationCPUTime: 13.49s
% Computational Cost: add. (217080->381), mult. (507075->466), div. (0->0), fcn. (354640->10), ass. (0->145)
t397 = -2 * qJD(3);
t356 = sin(qJ(2));
t360 = cos(qJ(2));
t383 = qJD(1) * qJD(2);
t326 = t356 * qJDD(1) + t360 * t383;
t363 = qJD(1) ^ 2;
t357 = sin(qJ(1));
t361 = cos(qJ(1));
t333 = -t361 * g(1) - t357 * g(2);
t320 = -t363 * pkin(1) + qJDD(1) * pkin(7) + t333;
t390 = t356 * t320;
t259 = qJDD(2) * pkin(2) - t326 * qJ(3) - t390 + (pkin(2) * t356 * t363 + qJ(3) * t383 - g(3)) * t360;
t301 = -t356 * g(3) + t360 * t320;
t327 = t360 * qJDD(1) - t356 * t383;
t387 = qJD(1) * t356;
t329 = qJD(2) * pkin(2) - qJ(3) * t387;
t391 = t360 ^ 2 * t363;
t260 = -pkin(2) * t391 + t327 * qJ(3) - qJD(2) * t329 + t301;
t352 = sin(pkin(10));
t392 = cos(pkin(10));
t314 = (t352 * t360 + t356 * t392) * qJD(1);
t233 = t259 * t392 - t352 * t260 + t314 * t397;
t386 = qJD(1) * t360;
t313 = t352 * t387 - t386 * t392;
t234 = t352 * t259 + t392 * t260 + t313 * t397;
t296 = t352 * t326 - t327 * t392;
t303 = qJD(2) * mrSges(4,1) - t314 * mrSges(4,3);
t282 = t313 * pkin(3) - t314 * qJ(4);
t362 = qJD(2) ^ 2;
t222 = -qJDD(2) * pkin(3) - t362 * qJ(4) + t314 * t282 + qJDD(4) - t233;
t297 = t326 * t392 + t352 * t327;
t385 = qJD(2) * t313;
t213 = (-t297 - t385) * pkin(8) + (t313 * t314 - qJDD(2)) * pkin(4) + t222;
t395 = 2 * qJD(4);
t220 = -t362 * pkin(3) + qJDD(2) * qJ(4) + qJD(2) * t395 - t313 * t282 + t234;
t306 = -qJD(2) * pkin(4) - t314 * pkin(8);
t312 = t313 ^ 2;
t215 = -t312 * pkin(4) + t296 * pkin(8) + qJD(2) * t306 + t220;
t355 = sin(qJ(5));
t359 = cos(qJ(5));
t211 = t355 * t213 + t359 * t215;
t277 = t355 * t313 + t359 * t314;
t242 = -t277 * qJD(5) + t359 * t296 - t355 * t297;
t276 = t359 * t313 - t355 * t314;
t252 = -t276 * mrSges(6,1) + t277 * mrSges(6,2);
t342 = -qJD(2) + qJD(5);
t267 = t342 * mrSges(6,1) - t277 * mrSges(6,3);
t341 = -qJDD(2) + qJDD(5);
t253 = -t276 * pkin(5) - t277 * pkin(9);
t340 = t342 ^ 2;
t206 = -t340 * pkin(5) + t341 * pkin(9) + t276 * t253 + t211;
t332 = t357 * g(1) - t361 * g(2);
t319 = -qJDD(1) * pkin(1) - t363 * pkin(7) - t332;
t263 = -t327 * pkin(2) - qJ(3) * t391 + t329 * t387 + qJDD(3) + t319;
t372 = t296 * pkin(3) + t263 + (-t297 + t385) * qJ(4);
t393 = pkin(3) * qJD(2);
t217 = -t296 * pkin(4) - t312 * pkin(8) - t372 + (t306 - t393 + t395) * t314;
t243 = t276 * qJD(5) + t355 * t296 + t359 * t297;
t207 = (-t276 * t342 - t243) * pkin(9) + (t277 * t342 - t242) * pkin(5) + t217;
t354 = sin(qJ(6));
t358 = cos(qJ(6));
t203 = -t354 * t206 + t358 * t207;
t264 = -t354 * t277 + t358 * t342;
t225 = t264 * qJD(6) + t358 * t243 + t354 * t341;
t241 = qJDD(6) - t242;
t265 = t358 * t277 + t354 * t342;
t244 = -t264 * mrSges(7,1) + t265 * mrSges(7,2);
t269 = qJD(6) - t276;
t245 = -t269 * mrSges(7,2) + t264 * mrSges(7,3);
t199 = m(7) * t203 + t241 * mrSges(7,1) - t225 * mrSges(7,3) - t265 * t244 + t269 * t245;
t204 = t358 * t206 + t354 * t207;
t224 = -t265 * qJD(6) - t354 * t243 + t358 * t341;
t246 = t269 * mrSges(7,1) - t265 * mrSges(7,3);
t200 = m(7) * t204 - t241 * mrSges(7,2) + t224 * mrSges(7,3) + t264 * t244 - t269 * t246;
t380 = -t354 * t199 + t358 * t200;
t186 = m(6) * t211 - t341 * mrSges(6,2) + t242 * mrSges(6,3) + t276 * t252 - t342 * t267 + t380;
t210 = t359 * t213 - t355 * t215;
t266 = -t342 * mrSges(6,2) + t276 * mrSges(6,3);
t205 = -t341 * pkin(5) - t340 * pkin(9) + t277 * t253 - t210;
t375 = -m(7) * t205 + t224 * mrSges(7,1) - t225 * mrSges(7,2) + t264 * t245 - t265 * t246;
t195 = m(6) * t210 + t341 * mrSges(6,1) - t243 * mrSges(6,3) - t277 * t252 + t342 * t266 + t375;
t181 = t359 * t186 - t355 * t195;
t304 = -qJD(2) * mrSges(5,1) + t314 * mrSges(5,2);
t377 = m(5) * t220 + qJDD(2) * mrSges(5,3) + qJD(2) * t304 + t181;
t283 = t313 * mrSges(5,1) - t314 * mrSges(5,3);
t388 = -t313 * mrSges(4,1) - t314 * mrSges(4,2) - t283;
t394 = -mrSges(4,3) - mrSges(5,2);
t174 = m(4) * t234 - qJDD(2) * mrSges(4,2) - qJD(2) * t303 + t296 * t394 + t313 * t388 + t377;
t302 = -qJD(2) * mrSges(4,2) - t313 * mrSges(4,3);
t180 = t355 * t186 + t359 * t195;
t305 = -t313 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t374 = -m(5) * t222 + qJDD(2) * mrSges(5,1) + qJD(2) * t305 - t180;
t175 = m(4) * t233 + qJDD(2) * mrSges(4,1) + qJD(2) * t302 + t297 * t394 + t314 * t388 + t374;
t167 = t352 * t174 + t392 * t175;
t300 = -t360 * g(3) - t390;
t316 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t356 + Ifges(3,2) * t360) * qJD(1);
t317 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t356 + Ifges(3,4) * t360) * qJD(1);
t273 = Ifges(4,4) * t314 - Ifges(4,2) * t313 + Ifges(4,6) * qJD(2);
t275 = Ifges(4,1) * t314 - Ifges(4,4) * t313 + Ifges(4,5) * qJD(2);
t270 = Ifges(5,5) * t314 + Ifges(5,6) * qJD(2) + Ifges(5,3) * t313;
t274 = Ifges(5,1) * t314 + Ifges(5,4) * qJD(2) + Ifges(5,5) * t313;
t228 = Ifges(7,5) * t265 + Ifges(7,6) * t264 + Ifges(7,3) * t269;
t230 = Ifges(7,1) * t265 + Ifges(7,4) * t264 + Ifges(7,5) * t269;
t193 = -mrSges(7,1) * t205 + mrSges(7,3) * t204 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t241 - t265 * t228 + t269 * t230;
t229 = Ifges(7,4) * t265 + Ifges(7,2) * t264 + Ifges(7,6) * t269;
t194 = mrSges(7,2) * t205 - mrSges(7,3) * t203 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t241 + t264 * t228 - t269 * t229;
t248 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t342;
t249 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t342;
t373 = mrSges(6,1) * t210 - mrSges(6,2) * t211 + Ifges(6,5) * t243 + Ifges(6,6) * t242 + Ifges(6,3) * t341 + pkin(5) * t375 + pkin(9) * t380 + t358 * t193 + t354 * t194 + t277 * t248 - t276 * t249;
t367 = mrSges(5,1) * t222 - mrSges(5,3) * t220 - Ifges(5,4) * t297 - Ifges(5,2) * qJDD(2) - Ifges(5,6) * t296 + pkin(4) * t180 + t314 * t270 - t313 * t274 + t373;
t365 = -mrSges(4,2) * t234 + t313 * t275 + qJ(4) * (-t296 * mrSges(5,2) - t313 * t283 + t377) + pkin(3) * (-t297 * mrSges(5,2) - t314 * t283 + t374) + mrSges(4,1) * t233 + t314 * t273 - Ifges(4,6) * t296 + Ifges(4,5) * t297 + Ifges(4,3) * qJDD(2) - t367;
t396 = mrSges(3,1) * t300 - mrSges(3,2) * t301 + Ifges(3,5) * t326 + Ifges(3,6) * t327 + Ifges(3,3) * qJDD(2) + pkin(2) * t167 + (t356 * t316 - t360 * t317) * qJD(1) + t365;
t189 = t358 * t199 + t354 * t200;
t272 = Ifges(5,4) * t314 + Ifges(5,2) * qJD(2) + Ifges(5,6) * t313;
t389 = -Ifges(4,5) * t314 + Ifges(4,6) * t313 - Ifges(4,3) * qJD(2) - t272;
t325 = (-mrSges(3,1) * t360 + mrSges(3,2) * t356) * qJD(1);
t331 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t386;
t165 = m(3) * t300 + qJDD(2) * mrSges(3,1) - t326 * mrSges(3,3) + qJD(2) * t331 - t325 * t387 + t167;
t330 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t387;
t381 = t392 * t174 - t352 * t175;
t166 = m(3) * t301 - qJDD(2) * mrSges(3,2) + t327 * mrSges(3,3) - qJD(2) * t330 + t325 * t386 + t381;
t382 = -t356 * t165 + t360 * t166;
t187 = -m(6) * t217 + t242 * mrSges(6,1) - t243 * mrSges(6,2) + t276 * t266 - t277 * t267 - t189;
t227 = (-(2 * qJD(4)) + t393) * t314 + t372;
t184 = m(5) * t227 + t296 * mrSges(5,1) - t297 * mrSges(5,3) - t314 * t304 + t313 * t305 + t187;
t247 = Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * t342;
t169 = mrSges(6,2) * t217 - mrSges(6,3) * t210 + Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t341 - pkin(9) * t189 - t354 * t193 + t358 * t194 + t276 * t247 - t342 * t248;
t369 = mrSges(7,1) * t203 - mrSges(7,2) * t204 + Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t241 + t265 * t229 - t264 * t230;
t170 = -mrSges(6,1) * t217 + mrSges(6,3) * t211 + Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t341 - pkin(5) * t189 - t277 * t247 + t342 * t249 - t369;
t370 = -mrSges(5,1) * t227 + mrSges(5,2) * t220 - pkin(4) * t187 - pkin(8) * t181 - t355 * t169 - t359 * t170;
t159 = -mrSges(4,1) * t263 + mrSges(4,3) * t234 - pkin(3) * t184 + t389 * t314 + (Ifges(4,4) - Ifges(5,5)) * t297 + (-Ifges(4,2) - Ifges(5,3)) * t296 + (Ifges(4,6) - Ifges(5,6)) * qJDD(2) + (t275 + t274) * qJD(2) + t370;
t371 = mrSges(5,2) * t222 - mrSges(5,3) * t227 + Ifges(5,1) * t297 + Ifges(5,4) * qJDD(2) + Ifges(5,5) * t296 - pkin(8) * t180 + qJD(2) * t270 + t359 * t169 - t355 * t170;
t160 = mrSges(4,2) * t263 - mrSges(4,3) * t233 + Ifges(4,1) * t297 - Ifges(4,4) * t296 + Ifges(4,5) * qJDD(2) - qJ(4) * t184 - qJD(2) * t273 + t313 * t389 + t371;
t315 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t356 + Ifges(3,6) * t360) * qJD(1);
t368 = m(4) * t263 + t296 * mrSges(4,1) + t297 * mrSges(4,2) + t313 * t302 + t314 * t303 + t184;
t155 = -mrSges(3,1) * t319 + mrSges(3,3) * t301 + Ifges(3,4) * t326 + Ifges(3,2) * t327 + Ifges(3,6) * qJDD(2) - pkin(2) * t368 + qJ(3) * t381 + qJD(2) * t317 + t159 * t392 + t352 * t160 - t315 * t387;
t157 = mrSges(3,2) * t319 - mrSges(3,3) * t300 + Ifges(3,1) * t326 + Ifges(3,4) * t327 + Ifges(3,5) * qJDD(2) - qJ(3) * t167 - qJD(2) * t316 - t352 * t159 + t160 * t392 + t315 * t386;
t366 = -m(3) * t319 + t327 * mrSges(3,1) - t326 * mrSges(3,2) - t330 * t387 + t331 * t386 - t368;
t376 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t366 + pkin(7) * t382 + t360 * t155 + t356 * t157;
t182 = m(2) * t332 + qJDD(1) * mrSges(2,1) - t363 * mrSges(2,2) + t366;
t163 = t360 * t165 + t356 * t166;
t161 = m(2) * t333 - t363 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t382;
t158 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t163 - t396;
t153 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - pkin(7) * t163 - t356 * t155 + t360 * t157;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t361 * t153 - t357 * t158 - pkin(6) * (t357 * t161 + t361 * t182), t153, t157, t160, -t313 * t272 + t371, t169, t194; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t357 * t153 + t361 * t158 + pkin(6) * (t361 * t161 - t357 * t182), t158, t155, t159, -t367, t170, t193; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t376, t376, t396, t365, Ifges(5,5) * t297 + Ifges(5,6) * qJDD(2) + Ifges(5,3) * t296 - qJD(2) * t274 + t314 * t272 - t370, t373, t369;];
m_new  = t1;
