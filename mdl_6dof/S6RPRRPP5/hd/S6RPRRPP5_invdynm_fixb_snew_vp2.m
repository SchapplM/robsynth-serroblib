% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:36:19
% EndTime: 2019-05-05 21:36:34
% DurationCPUTime: 6.91s
% Computational Cost: add. (104405->363), mult. (244494->421), div. (0->0), fcn. (175623->8), ass. (0->142)
t345 = qJD(1) ^ 2;
t338 = sin(pkin(9));
t380 = qJD(1) * t338;
t339 = cos(pkin(9));
t400 = qJD(1) * t339;
t342 = sin(qJ(1));
t343 = cos(qJ(1));
t323 = -t343 * g(1) - t342 * g(2);
t316 = -t345 * pkin(1) + qJDD(1) * qJ(2) + t323;
t377 = qJD(1) * qJD(2);
t372 = -t339 * g(3) - 0.2e1 * t338 * t377;
t389 = pkin(2) * t339;
t279 = (-pkin(7) * qJDD(1) + t345 * t389 - t316) * t338 + t372;
t303 = -t338 * g(3) + (t316 + 0.2e1 * t377) * t339;
t375 = qJDD(1) * t339;
t335 = t339 ^ 2;
t384 = t335 * t345;
t280 = -pkin(2) * t384 + pkin(7) * t375 + t303;
t341 = sin(qJ(3));
t392 = cos(qJ(3));
t222 = t341 * t279 + t392 * t280;
t373 = t339 * t392;
t314 = qJD(1) * t373 - t341 * t380;
t360 = t392 * t338 + t339 * t341;
t315 = t360 * qJD(1);
t298 = -t314 * pkin(3) - t315 * pkin(8);
t344 = qJD(3) ^ 2;
t213 = -t344 * pkin(3) + qJDD(3) * pkin(8) + t314 * t298 + t222;
t334 = t338 ^ 2;
t322 = t342 * g(1) - t343 * g(2);
t366 = qJDD(2) - t322;
t299 = (-pkin(1) - t389) * qJDD(1) + (-qJ(2) + (-t334 - t335) * pkin(7)) * t345 + t366;
t376 = qJDD(1) * t338;
t378 = t315 * qJD(3);
t300 = qJDD(1) * t373 - t341 * t376 - t378;
t379 = t314 * qJD(3);
t301 = qJDD(1) * t360 + t379;
t215 = (-t301 - t379) * pkin(8) + (-t300 + t378) * pkin(3) + t299;
t340 = sin(qJ(4));
t391 = cos(qJ(4));
t209 = t391 * t213 + t340 * t215;
t305 = -t391 * qJD(3) + t340 * t315;
t306 = t340 * qJD(3) + t391 * t315;
t263 = t305 * pkin(4) - t306 * qJ(5);
t297 = qJDD(4) - t300;
t312 = qJD(4) - t314;
t311 = t312 ^ 2;
t393 = 2 * qJD(5);
t204 = -t311 * pkin(4) + t297 * qJ(5) - t305 * t263 + t312 * t393 + t209;
t276 = -t312 * mrSges(6,1) + t306 * mrSges(6,2);
t399 = m(6) * t204 + t297 * mrSges(6,3) + t312 * t276;
t208 = -t340 * t213 + t391 * t215;
t206 = -t297 * pkin(4) - t311 * qJ(5) + t306 * t263 + qJDD(5) - t208;
t270 = -t305 * mrSges(6,2) + t312 * mrSges(6,3);
t398 = -m(6) * t206 + t297 * mrSges(6,1) + t312 * t270;
t259 = -t305 * qJD(4) + t340 * qJDD(3) + t391 * t301;
t221 = t392 * t279 - t341 * t280;
t356 = qJDD(3) * pkin(3) + t344 * pkin(8) - t315 * t298 + t221;
t386 = t305 * t312;
t397 = (-t259 + t386) * qJ(5) - t356;
t271 = t312 * mrSges(7,2) + t305 * mrSges(7,3);
t394 = -0.2e1 * t306;
t194 = qJD(6) * t394 + (-t259 - t386) * qJ(6) + (t305 * t306 - t297) * pkin(5) + t206;
t265 = -t305 * mrSges(7,1) + t306 * mrSges(7,2);
t367 = -m(7) * t194 + t259 * mrSges(7,3) + t306 * t265;
t190 = -t297 * mrSges(7,1) - t312 * t271 - t367;
t274 = -t312 * mrSges(7,1) - t306 * mrSges(7,3);
t258 = t306 * qJD(4) - t391 * qJDD(3) + t340 * t301;
t273 = -t312 * pkin(5) - t306 * qJ(6);
t304 = t305 ^ 2;
t198 = -t304 * pkin(5) + t258 * qJ(6) + 0.2e1 * qJD(6) * t305 + t312 * t273 + t204;
t374 = m(7) * t198 + t258 * mrSges(7,3) + t305 * t265;
t191 = t297 * mrSges(7,2) + t312 * t274 + t374;
t230 = Ifges(6,5) * t306 + Ifges(6,6) * t312 + Ifges(6,3) * t305;
t234 = Ifges(5,4) * t306 - Ifges(5,2) * t305 + Ifges(5,6) * t312;
t237 = Ifges(5,1) * t306 - Ifges(5,4) * t305 + Ifges(5,5) * t312;
t264 = t305 * mrSges(6,1) - t306 * mrSges(6,3);
t236 = Ifges(6,1) * t306 + Ifges(6,4) * t312 + Ifges(6,5) * t305;
t232 = Ifges(7,4) * t306 + Ifges(7,2) * t305 - Ifges(7,6) * t312;
t235 = Ifges(7,1) * t306 + Ifges(7,4) * t305 - Ifges(7,5) * t312;
t359 = mrSges(7,1) * t194 - mrSges(7,2) * t198 + Ifges(7,5) * t259 + Ifges(7,6) * t258 - Ifges(7,3) * t297 + t306 * t232 - t305 * t235;
t351 = mrSges(6,1) * t206 - mrSges(6,3) * t204 - Ifges(6,4) * t259 - Ifges(6,2) * t297 - Ifges(6,6) * t258 + pkin(5) * t190 - t305 * t236 + t359;
t396 = (t234 - t230) * t306 + mrSges(5,1) * t208 - mrSges(5,2) * t209 + Ifges(5,5) * t259 - Ifges(5,6) * t258 + Ifges(5,3) * t297 + pkin(4) * (-t259 * mrSges(6,2) - t306 * t264 - t190 + t398) + qJ(5) * (-t258 * mrSges(6,2) - t305 * t264 + t191 + t399) + t305 * t237 - t351;
t293 = -t314 * mrSges(4,1) + t315 * mrSges(4,2);
t308 = qJD(3) * mrSges(4,1) - t315 * mrSges(4,3);
t272 = -t312 * mrSges(5,2) - t305 * mrSges(5,3);
t381 = -t305 * mrSges(5,1) - t306 * mrSges(5,2) - t264;
t388 = -mrSges(5,3) - mrSges(6,2);
t181 = m(5) * t208 + (t271 + t272) * t312 + t381 * t306 + (mrSges(5,1) + mrSges(7,1)) * t297 + t388 * t259 + t367 + t398;
t275 = t312 * mrSges(5,1) - t306 * mrSges(5,3);
t184 = m(5) * t209 + (t274 - t275) * t312 + t381 * t305 + (-mrSges(5,2) + mrSges(7,2)) * t297 + t388 * t258 + t374 + t399;
t369 = -t340 * t181 + t391 * t184;
t174 = m(4) * t222 - qJDD(3) * mrSges(4,2) + t300 * mrSges(4,3) - qJD(3) * t308 + t314 * t293 + t369;
t307 = -qJD(3) * mrSges(4,2) + t314 * mrSges(4,3);
t201 = -t304 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t258 + (-pkin(4) * t312 + t273 + t393) * t306 - t397;
t192 = -m(7) * t201 + t258 * mrSges(7,1) - t259 * mrSges(7,2) + t305 * t271 - t306 * t274;
t207 = qJD(5) * t394 + (t306 * t312 + t258) * pkin(4) + t397;
t189 = m(6) * t207 + t258 * mrSges(6,1) - t259 * mrSges(6,3) + t305 * t270 - t306 * t276 + t192;
t347 = m(5) * t356 - t258 * mrSges(5,1) - t259 * mrSges(5,2) - t305 * t272 - t306 * t275 - t189;
t179 = m(4) * t221 + qJDD(3) * mrSges(4,1) - t301 * mrSges(4,3) + qJD(3) * t307 - t315 * t293 + t347;
t167 = t341 * t174 + t392 * t179;
t302 = -t338 * t316 + t372;
t229 = Ifges(7,5) * t306 + Ifges(7,6) * t305 - Ifges(7,3) * t312;
t358 = mrSges(7,1) * t201 - mrSges(7,3) * t198 - Ifges(7,4) * t259 - Ifges(7,2) * t258 + Ifges(7,6) * t297 + t312 * t235;
t350 = mrSges(6,1) * t207 - mrSges(6,2) * t204 + pkin(5) * t192 + qJ(6) * t191 - t358;
t233 = Ifges(6,4) * t306 + Ifges(6,2) * t312 + Ifges(6,6) * t305;
t383 = -Ifges(5,5) * t306 + Ifges(5,6) * t305 - Ifges(5,3) * t312 - t233;
t163 = mrSges(5,3) * t209 + mrSges(5,1) * t356 - pkin(4) * t189 - t350 + (t237 + t236) * t312 + (t229 + t383) * t306 + (Ifges(5,6) - Ifges(6,6)) * t297 + (Ifges(5,4) - Ifges(6,5)) * t259 + (-Ifges(5,2) - Ifges(6,3)) * t258;
t357 = mrSges(7,2) * t201 - mrSges(7,3) * t194 + Ifges(7,1) * t259 + Ifges(7,4) * t258 - Ifges(7,5) * t297 + t305 * t229;
t349 = mrSges(6,2) * t206 - mrSges(6,3) * t207 + Ifges(6,1) * t259 + Ifges(6,4) * t297 + Ifges(6,5) * t258 - qJ(6) * t190 + t312 * t230 + t357;
t169 = Ifges(5,5) * t297 - Ifges(5,4) * t258 + Ifges(5,1) * t259 - mrSges(5,2) * t356 - mrSges(5,3) * t208 - qJ(5) * t189 + t349 + (-t234 + t232) * t312 + t383 * t305;
t282 = Ifges(4,4) * t315 + Ifges(4,2) * t314 + Ifges(4,6) * qJD(3);
t283 = Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * qJD(3);
t353 = -mrSges(4,1) * t221 + mrSges(4,2) * t222 - Ifges(4,5) * t301 - Ifges(4,6) * t300 - Ifges(4,3) * qJDD(3) - pkin(3) * t347 - pkin(8) * t369 - t391 * t163 - t340 * t169 - t315 * t282 + t314 * t283;
t364 = Ifges(3,4) * t338 + Ifges(3,2) * t339;
t365 = Ifges(3,1) * t338 + Ifges(3,4) * t339;
t395 = -mrSges(3,1) * t302 + mrSges(3,2) * t303 - pkin(2) * t167 - (t364 * t380 - t365 * t400) * qJD(1) + t353;
t387 = mrSges(3,2) * t338;
t385 = t312 * t232;
t176 = t391 * t181 + t340 * t184;
t361 = mrSges(3,3) * qJDD(1) + t345 * (-mrSges(3,1) * t339 + t387);
t165 = m(3) * t302 - t338 * t361 + t167;
t370 = t392 * t174 - t341 * t179;
t166 = m(3) * t303 + t339 * t361 + t370;
t371 = -t338 * t165 + t339 * t166;
t363 = Ifges(3,5) * t338 + Ifges(3,6) * t339;
t281 = Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * qJD(3);
t157 = mrSges(4,2) * t299 - mrSges(4,3) * t221 + Ifges(4,1) * t301 + Ifges(4,4) * t300 + Ifges(4,5) * qJDD(3) - pkin(8) * t176 - qJD(3) * t282 - t340 * t163 + t391 * t169 + t314 * t281;
t161 = -mrSges(4,1) * t299 + mrSges(4,3) * t222 + Ifges(4,4) * t301 + Ifges(4,2) * t300 + Ifges(4,6) * qJDD(3) - pkin(3) * t176 + qJD(3) * t283 - t315 * t281 - t396;
t313 = -qJDD(1) * pkin(1) - t345 * qJ(2) + t366;
t318 = t363 * qJD(1);
t354 = m(4) * t299 - t300 * mrSges(4,1) + t301 * mrSges(4,2) - t314 * t307 + t315 * t308 + t176;
t153 = -mrSges(3,1) * t313 + mrSges(3,3) * t303 - pkin(2) * t354 + pkin(7) * t370 + t364 * qJDD(1) + t341 * t157 + t392 * t161 - t318 * t380;
t156 = mrSges(3,2) * t313 - mrSges(3,3) * t302 - pkin(7) * t167 + t365 * qJDD(1) + t392 * t157 - t341 * t161 + t318 * t400;
t352 = -m(3) * t313 + mrSges(3,1) * t375 - t354 + (t334 * t345 + t384) * mrSges(3,3);
t355 = -mrSges(2,2) * t323 + qJ(2) * t371 + t339 * t153 + t338 * t156 + pkin(1) * (-mrSges(3,2) * t376 + t352) + mrSges(2,1) * t322 + Ifges(2,3) * qJDD(1);
t170 = -t345 * mrSges(2,2) + m(2) * t322 + (mrSges(2,1) - t387) * qJDD(1) + t352;
t160 = t339 * t165 + t338 * t166;
t158 = m(2) * t323 - t345 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t371;
t154 = t345 * Ifges(2,5) + mrSges(2,3) * t323 + mrSges(2,1) * g(3) - pkin(1) * t160 + (Ifges(2,6) - t363) * qJDD(1) + t395;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t322 + Ifges(2,5) * qJDD(1) - t345 * Ifges(2,6) - qJ(2) * t160 - t338 * t153 + t339 * t156;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t343 * t151 - t342 * t154 - pkin(6) * (t342 * t158 + t343 * t170), t151, t156, t157, t169, -t305 * t233 + t349 + t385, t357 + t385; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t151 + t343 * t154 + pkin(6) * (t343 * t158 - t342 * t170), t154, t153, t161, t163, -t306 * t230 - t351, -t306 * t229 - t358; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t355, t355, qJDD(1) * t363 - t395, -t353, t396, -t312 * t236 + Ifges(6,6) * t297 + Ifges(6,3) * t258 + Ifges(6,5) * t259 + t350 + (-t229 + t233) * t306, t359;];
m_new  = t1;
