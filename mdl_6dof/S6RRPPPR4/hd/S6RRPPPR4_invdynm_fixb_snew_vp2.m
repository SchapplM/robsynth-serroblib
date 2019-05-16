% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-05-06 08:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:41:42
% EndTime: 2019-05-06 08:41:55
% DurationCPUTime: 7.02s
% Computational Cost: add. (79361->391), mult. (175422->455), div. (0->0), fcn. (99821->8), ass. (0->142)
t400 = -2 * qJD(4);
t350 = sin(qJ(1));
t353 = cos(qJ(1));
t327 = -g(1) * t353 - g(2) * t350;
t355 = qJD(1) ^ 2;
t294 = -pkin(1) * t355 + qJDD(1) * pkin(7) + t327;
t349 = sin(qJ(2));
t352 = cos(qJ(2));
t261 = -t352 * g(3) - t349 * t294;
t262 = -t349 * g(3) + t352 * t294;
t288 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t349 + Ifges(3,4) * t352) * qJD(1);
t314 = (mrSges(4,2) * t352 - mrSges(4,3) * t349) * qJD(1);
t383 = qJD(1) * qJD(2);
t380 = t352 * t383;
t316 = qJDD(1) * t349 + t380;
t379 = t349 * t383;
t317 = qJDD(1) * t352 - t379;
t385 = qJD(1) * t352;
t324 = -mrSges(4,1) * t385 - qJD(2) * mrSges(4,3);
t345 = sin(pkin(9));
t346 = cos(pkin(9));
t275 = qJDD(2) * t345 + t346 * t317;
t301 = qJD(2) * t346 - t345 * t385;
t386 = qJD(1) * t349;
t277 = -pkin(5) * t386 - pkin(8) * t301;
t300 = qJD(2) * t345 + t346 * t385;
t298 = t300 ^ 2;
t382 = qJD(3) * qJD(2);
t335 = -0.2e1 * t382;
t323 = pkin(3) * t386 - qJD(2) * qJ(4);
t344 = t352 ^ 2;
t313 = (-pkin(2) * t352 - qJ(3) * t349) * qJD(1);
t354 = qJD(2) ^ 2;
t375 = t354 * pkin(2) - qJDD(2) * qJ(3) - t313 * t385 - t262;
t367 = t344 * t355 * qJ(4) - t317 * pkin(3) - qJD(2) * t323 - qJDD(4) + t375;
t396 = 2 * qJD(5);
t276 = qJDD(2) * t346 - t317 * t345;
t381 = t300 * t386;
t398 = (-t276 + t381) * qJ(5);
t210 = -t298 * pkin(8) + t335 + (-pkin(4) - pkin(5)) * t275 - t398 + (-pkin(4) * t386 + t277 + t396) * t301 + t367;
t348 = sin(qJ(6));
t351 = cos(qJ(6));
t255 = t300 * t348 + t301 * t351;
t224 = -qJD(6) * t255 + t275 * t351 - t276 * t348;
t254 = t300 * t351 - t301 * t348;
t225 = qJD(6) * t254 + t275 * t348 + t276 * t351;
t330 = qJD(6) - t386;
t241 = -mrSges(7,2) * t330 + mrSges(7,3) * t254;
t242 = mrSges(7,1) * t330 - mrSges(7,3) * t255;
t202 = -m(7) * t210 + t224 * mrSges(7,1) - t225 * mrSges(7,2) + t254 * t241 - t255 * t242;
t228 = 0.2e1 * t382 - t367;
t215 = -0.2e1 * qJD(5) * t301 + t398 + (t301 * t386 + t275) * pkin(4) + t228;
t271 = -mrSges(6,2) * t300 + mrSges(6,3) * t386;
t274 = -mrSges(6,1) * t386 + mrSges(6,2) * t301;
t197 = m(6) * t215 + t275 * mrSges(6,1) - t276 * mrSges(6,3) + t300 * t271 - t301 * t274 + t202;
t272 = -mrSges(5,2) * t386 - mrSges(5,3) * t300;
t273 = mrSges(5,1) * t386 - mrSges(5,3) * t301;
t193 = -m(5) * t228 - t275 * mrSges(5,1) - t276 * mrSges(5,2) - t300 * t272 - t301 * t273 - t197;
t238 = t335 + t375;
t325 = mrSges(4,1) * t386 + qJD(2) * mrSges(4,2);
t359 = -m(4) * t238 + qJDD(2) * mrSges(4,3) + qJD(2) * t325 + t314 * t385 - t193;
t326 = t350 * g(1) - t353 * g(2);
t374 = -qJDD(1) * pkin(1) - t326;
t363 = pkin(2) * t379 - 0.2e1 * qJD(3) * t386 + (-t316 - t380) * qJ(3) + t374;
t220 = -t323 * t386 + (-pkin(3) * t344 - pkin(7)) * t355 + (-pkin(2) - qJ(4)) * t317 + t363;
t240 = -qJDD(2) * pkin(2) - t354 * qJ(3) + t313 * t386 + qJDD(3) - t261;
t232 = (-t349 * t352 * t355 - qJDD(2)) * qJ(4) + (t316 - t380) * pkin(3) + t240;
t213 = t346 * t220 + t345 * t232 + t300 * t400;
t250 = Ifges(6,1) * t301 + Ifges(6,4) * t386 + Ifges(6,5) * t300;
t251 = Ifges(5,1) * t301 - Ifges(5,4) * t300 + Ifges(5,5) * t386;
t212 = -t345 * t220 + t346 * t232 + t301 * t400;
t256 = pkin(4) * t300 - qJ(5) * t301;
t392 = t349 ^ 2 * t355;
t209 = -t316 * pkin(4) - qJ(5) * t392 + t301 * t256 + qJDD(5) - t212;
t203 = (-t276 - t381) * pkin(8) + (t300 * t301 - t316) * pkin(5) + t209;
t207 = -pkin(4) * t392 + t316 * qJ(5) - t300 * t256 + t386 * t396 + t213;
t204 = -pkin(5) * t298 + pkin(8) * t275 + t277 * t386 + t207;
t200 = t203 * t351 - t204 * t348;
t234 = -mrSges(7,1) * t254 + mrSges(7,2) * t255;
t312 = qJDD(6) - t316;
t195 = m(7) * t200 + mrSges(7,1) * t312 - mrSges(7,3) * t225 - t234 * t255 + t241 * t330;
t201 = t203 * t348 + t204 * t351;
t196 = m(7) * t201 - mrSges(7,2) * t312 + mrSges(7,3) * t224 + t234 * t254 - t242 * t330;
t186 = -t348 * t195 + t351 * t196;
t229 = Ifges(7,5) * t255 + Ifges(7,6) * t254 + Ifges(7,3) * t330;
t231 = Ifges(7,1) * t255 + Ifges(7,4) * t254 + Ifges(7,5) * t330;
t188 = -mrSges(7,1) * t210 + mrSges(7,3) * t201 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t312 - t229 * t255 + t231 * t330;
t230 = Ifges(7,4) * t255 + Ifges(7,2) * t254 + Ifges(7,6) * t330;
t189 = mrSges(7,2) * t210 - mrSges(7,3) * t200 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t312 + t229 * t254 - t230 * t330;
t364 = -mrSges(6,1) * t215 + mrSges(6,2) * t207 - pkin(5) * t202 - pkin(8) * t186 - t351 * t188 - t348 * t189;
t248 = Ifges(6,4) * t301 + Ifges(6,2) * t386 + Ifges(6,6) * t300;
t390 = -Ifges(5,5) * t301 + Ifges(5,6) * t300 - Ifges(5,3) * t386 - t248;
t164 = -mrSges(5,1) * t228 + mrSges(5,3) * t213 - pkin(4) * t197 + (Ifges(5,6) - Ifges(6,6)) * t316 + t390 * t301 + (Ifges(5,4) - Ifges(6,5)) * t276 + (-Ifges(5,2) - Ifges(6,3)) * t275 + (t250 + t251) * t386 + t364;
t249 = Ifges(5,4) * t301 - Ifges(5,2) * t300 + Ifges(5,6) * t386;
t185 = t351 * t195 + t348 * t196;
t246 = Ifges(6,5) * t301 + Ifges(6,6) * t386 + Ifges(6,3) * t300;
t366 = mrSges(6,2) * t209 - mrSges(6,3) * t215 + Ifges(6,1) * t276 + Ifges(6,4) * t316 + Ifges(6,5) * t275 - pkin(8) * t185 - t348 * t188 + t351 * t189 + t246 * t386;
t169 = mrSges(5,2) * t228 - mrSges(5,3) * t212 + Ifges(5,1) * t276 - Ifges(5,4) * t275 + Ifges(5,5) * t316 - qJ(5) * t197 - t249 * t386 + t390 * t300 + t366;
t372 = m(6) * t207 + t316 * mrSges(6,3) + t274 * t386 + t186;
t257 = mrSges(6,1) * t300 - mrSges(6,3) * t301;
t389 = -mrSges(5,1) * t300 - mrSges(5,2) * t301 - t257;
t394 = -mrSges(5,3) - mrSges(6,2);
t178 = m(5) * t213 - t316 * mrSges(5,2) - t273 * t386 + t394 * t275 + t389 * t300 + t372;
t368 = -m(6) * t209 + t316 * mrSges(6,1) + t271 * t386 - t185;
t180 = m(5) * t212 + t316 * mrSges(5,1) + t272 * t386 + t394 * t276 + t389 * t301 + t368;
t175 = t345 * t178 + t346 * t180;
t290 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t349 - Ifges(4,6) * t352) * qJD(1);
t365 = -mrSges(4,2) * t240 + mrSges(4,3) * t238 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t316 + Ifges(4,5) * t317 + qJ(4) * t175 + t345 * t164 - t346 * t169 - t290 * t385;
t370 = -m(4) * t240 - t316 * mrSges(4,1) - t175;
t289 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t349 - Ifges(4,3) * t352) * qJD(1);
t387 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t349 + Ifges(3,2) * t352) * qJD(1) - t289;
t399 = qJD(1) * (-t352 * t288 + t387 * t349) + mrSges(3,1) * t261 - mrSges(3,2) * t262 + Ifges(3,5) * t316 + Ifges(3,6) * t317 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t324 - t314 * t386 + t370) + qJ(3) * (t317 * mrSges(4,1) + t359) - t365;
t395 = t355 * pkin(7);
t393 = Ifges(3,4) + Ifges(4,6);
t176 = t346 * t178 - t345 * t180;
t291 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t349 - Ifges(4,5) * t352) * qJD(1);
t388 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t349 + Ifges(3,6) * t352) * qJD(1) + t291;
t315 = (-mrSges(3,1) * t352 + mrSges(3,2) * t349) * qJD(1);
t322 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t385;
t172 = m(3) * t261 - t316 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t322 - t324) * qJD(2) + (-t314 - t315) * t386 + t370;
t321 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t386;
t191 = t359 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t317 + t315 * t385 - qJD(2) * t321 + m(3) * t262;
t378 = -t172 * t349 + t352 * t191;
t235 = -t317 * pkin(2) + t363 - t395;
t373 = -m(4) * t235 - t317 * mrSges(4,2) + t325 * t386 - t176;
t371 = mrSges(7,1) * t200 - mrSges(7,2) * t201 + Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t312 + t255 * t230 - t254 * t231;
t173 = -t316 * mrSges(4,3) + t324 * t385 - t373;
t293 = t374 - t395;
t362 = -mrSges(4,1) * t238 + mrSges(4,2) * t235 - pkin(3) * t193 - qJ(4) * t176 - t346 * t164 - t345 * t169;
t161 = -mrSges(3,1) * t293 + mrSges(3,3) * t262 - pkin(2) * t173 + (Ifges(3,2) + Ifges(4,3)) * t317 + t393 * t316 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t288 - t290) * qJD(2) - t388 * t386 + t362;
t361 = mrSges(6,1) * t209 - mrSges(6,3) * t207 - Ifges(6,4) * t276 - Ifges(6,2) * t316 - Ifges(6,6) * t275 + pkin(5) * t185 + t301 * t246 - t300 * t250 + t371;
t357 = mrSges(5,2) * t213 - t300 * t251 - qJ(5) * (-t275 * mrSges(6,2) - t300 * t257 + t372) - pkin(4) * (-t276 * mrSges(6,2) - t301 * t257 + t368) - mrSges(5,1) * t212 - t301 * t249 + Ifges(5,6) * t275 - Ifges(5,5) * t276 - Ifges(5,3) * t316 + t361;
t356 = -mrSges(4,1) * t240 + mrSges(4,3) * t235 - pkin(3) * t175 + t357;
t163 = -t356 + t393 * t317 + (Ifges(4,2) + Ifges(3,1)) * t316 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t387 * qJD(2) + mrSges(3,2) * t293 - mrSges(3,3) * t261 - qJ(3) * t173 + t388 * t385;
t360 = -m(3) * t293 + t322 * t385 + t317 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t316 + (-t321 * t349 - t324 * t352) * qJD(1) + t373;
t369 = mrSges(2,1) * t326 - mrSges(2,2) * t327 + Ifges(2,3) * qJDD(1) + pkin(1) * t360 + pkin(7) * t378 + t352 * t161 + t349 * t163;
t170 = m(2) * t326 + qJDD(1) * mrSges(2,1) - t355 * mrSges(2,2) + t360;
t167 = t172 * t352 + t191 * t349;
t165 = m(2) * t327 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t378;
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t327 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t399;
t158 = -mrSges(2,2) * g(3) - mrSges(2,3) * t326 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t355 - pkin(7) * t167 - t161 * t349 + t163 * t352;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t353 * t158 - t350 * t159 - pkin(6) * (t165 * t350 + t170 * t353), t158, t163, -t289 * t386 - t365, t169, -t248 * t300 + t366, t189; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t350 * t158 + t353 * t159 + pkin(6) * (t165 * t353 - t170 * t350), t159, t161, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t316 - Ifges(4,6) * t317 - qJD(2) * t289 - t291 * t385 + t356, t164, -t361, t188; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t369, t369, t399, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t316 - Ifges(4,3) * t317 + qJD(2) * t290 + t291 * t386 - t362, -t357, Ifges(6,5) * t276 + Ifges(6,6) * t316 + Ifges(6,3) * t275 + t301 * t248 - t250 * t386 - t364, t371;];
m_new  = t1;
