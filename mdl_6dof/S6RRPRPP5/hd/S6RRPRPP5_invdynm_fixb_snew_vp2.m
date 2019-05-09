% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-05-06 12:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:41:48
% EndTime: 2019-05-06 12:42:00
% DurationCPUTime: 4.99s
% Computational Cost: add. (44215->393), mult. (89030->435), div. (0->0), fcn. (48189->6), ass. (0->135)
t347 = cos(qJ(2));
t344 = sin(qJ(2));
t383 = qJD(1) * qJD(2);
t378 = t344 * t383;
t314 = qJDD(1) * t347 - t378;
t385 = qJD(1) * t344;
t323 = pkin(3) * t385 - qJD(2) * pkin(8);
t341 = t347 ^ 2;
t350 = qJD(1) ^ 2;
t345 = sin(qJ(1));
t348 = cos(qJ(1));
t325 = -g(1) * t348 - g(2) * t345;
t288 = -pkin(1) * t350 + qJDD(1) * pkin(7) + t325;
t274 = -t344 * g(3) + t347 * t288;
t310 = (-pkin(2) * t347 - qJ(3) * t344) * qJD(1);
t349 = qJD(2) ^ 2;
t384 = qJD(1) * t347;
t372 = t349 * pkin(2) - qJDD(2) * qJ(3) - t310 * t384 - t274;
t363 = pkin(8) * t341 * t350 - pkin(3) * t314 - qJD(2) * t323 + t372;
t382 = qJD(3) * qJD(2);
t212 = 0.2e1 * t382 - t363;
t343 = sin(qJ(4));
t346 = cos(qJ(4));
t309 = qJD(2) * t346 - t343 * t384;
t258 = qJD(4) * t309 + qJDD(2) * t343 + t346 * t314;
t329 = qJD(4) + t385;
t400 = -0.2e1 * t309;
t308 = qJD(2) * t343 + t346 * t384;
t259 = -qJD(4) * t308 + qJDD(2) * t346 - t314 * t343;
t393 = t308 * t329;
t405 = (-t259 + t393) * qJ(5);
t203 = qJD(5) * t400 + t405 + (t309 * t329 + t258) * pkin(4) + t212;
t268 = -mrSges(6,2) * t308 + mrSges(6,3) * t329;
t270 = -mrSges(7,1) * t329 - mrSges(7,3) * t309;
t272 = -mrSges(6,1) * t329 + mrSges(6,2) * t309;
t269 = -pkin(5) * t329 - qJ(6) * t309;
t306 = t308 ^ 2;
t333 = -0.2e1 * t382;
t399 = 0.2e1 * qJD(5);
t196 = -qJ(6) * t306 + qJDD(6) + t333 + (-pkin(4) - pkin(5)) * t258 - t405 + (-pkin(4) * t329 + t269 + t399) * t309 + t363;
t266 = mrSges(7,2) * t329 + mrSges(7,3) * t308;
t373 = m(7) * t196 - t258 * mrSges(7,1) - t308 * t266;
t185 = m(6) * t203 + t258 * mrSges(6,1) + t308 * t268 - t373 - (t270 + t272) * t309 - (mrSges(7,2) + mrSges(6,3)) * t259;
t267 = -mrSges(5,2) * t329 - mrSges(5,3) * t308;
t271 = mrSges(5,1) * t329 - mrSges(5,3) * t309;
t407 = -m(5) * t212 - t258 * mrSges(5,1) - t259 * mrSges(5,2) - t308 * t267 - t309 * t271 - t185;
t273 = -t347 * g(3) - t344 * t288;
t283 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t344 + Ifges(3,4) * t347) * qJD(1);
t311 = (mrSges(4,2) * t347 - mrSges(4,3) * t344) * qJD(1);
t379 = t347 * t383;
t313 = qJDD(1) * t344 + t379;
t320 = -mrSges(4,1) * t384 - qJD(2) * mrSges(4,3);
t220 = t333 + t372;
t321 = mrSges(4,1) * t385 + qJD(2) * mrSges(4,2);
t354 = -m(4) * t220 + qJDD(2) * mrSges(4,3) + qJD(2) * t321 + t311 * t384 - t407;
t324 = g(1) * t345 - t348 * g(2);
t371 = -qJDD(1) * pkin(1) - t324;
t360 = pkin(2) * t378 - 0.2e1 * qJD(3) * t385 + (-t313 - t379) * qJ(3) + t371;
t209 = -t323 * t385 + (-pkin(3) * t341 - pkin(7)) * t350 + (-pkin(2) - pkin(8)) * t314 + t360;
t222 = -qJDD(2) * pkin(2) - qJ(3) * t349 + t310 * t385 + qJDD(3) - t273;
t213 = (-t344 * t347 * t350 - qJDD(2)) * pkin(8) + (t313 - t379) * pkin(3) + t222;
t206 = t346 * t209 + t343 * t213;
t229 = Ifges(7,5) * t309 + Ifges(7,6) * t308 - Ifges(7,3) * t329;
t236 = Ifges(6,1) * t309 + Ifges(6,4) * t329 + Ifges(6,5) * t308;
t237 = Ifges(5,1) * t309 - Ifges(5,4) * t308 + Ifges(5,5) * t329;
t307 = qJDD(4) + t313;
t262 = pkin(4) * t308 - qJ(5) * t309;
t326 = t329 ^ 2;
t199 = -pkin(4) * t326 + t307 * qJ(5) - t308 * t262 + t329 * t399 + t206;
t194 = -pkin(5) * t306 + qJ(6) * t258 + 0.2e1 * qJD(6) * t308 + t269 * t329 + t199;
t264 = -mrSges(7,1) * t308 + mrSges(7,2) * t309;
t381 = m(7) * t194 + t258 * mrSges(7,3) + t308 * t264;
t187 = mrSges(7,2) * t307 + t270 * t329 + t381;
t235 = Ifges(7,1) * t309 + Ifges(7,4) * t308 - Ifges(7,5) * t329;
t367 = mrSges(7,1) * t196 - mrSges(7,3) * t194 - Ifges(7,4) * t259 - Ifges(7,2) * t258 + Ifges(7,6) * t307 + t329 * t235;
t358 = mrSges(6,1) * t203 - mrSges(6,2) * t199 + pkin(5) * (-t259 * mrSges(7,2) - t309 * t270 - t373) + qJ(6) * t187 - t367;
t233 = Ifges(6,4) * t309 + Ifges(6,2) * t329 + Ifges(6,6) * t308;
t390 = -Ifges(5,5) * t309 + Ifges(5,6) * t308 - Ifges(5,3) * t329 - t233;
t163 = -t358 + (t237 + t236) * t329 + (t229 + t390) * t309 + (Ifges(5,6) - Ifges(6,6)) * t307 + (Ifges(5,4) - Ifges(6,5)) * t259 + (-Ifges(5,2) - Ifges(6,3)) * t258 - pkin(4) * t185 + mrSges(5,3) * t206 - mrSges(5,1) * t212;
t205 = -t343 * t209 + t213 * t346;
t232 = Ifges(7,4) * t309 + Ifges(7,2) * t308 - Ifges(7,6) * t329;
t234 = Ifges(5,4) * t309 - Ifges(5,2) * t308 + Ifges(5,6) * t329;
t201 = -pkin(4) * t307 - qJ(5) * t326 + t309 * t262 + qJDD(5) - t205;
t190 = qJD(6) * t400 + (-t259 - t393) * qJ(6) + (t308 * t309 - t307) * pkin(5) + t201;
t374 = -m(7) * t190 + t259 * mrSges(7,3) + t309 * t264;
t186 = -mrSges(7,1) * t307 - t266 * t329 - t374;
t230 = Ifges(6,5) * t309 + Ifges(6,6) * t329 + Ifges(6,3) * t308;
t366 = mrSges(7,2) * t196 - mrSges(7,3) * t190 + Ifges(7,1) * t259 + Ifges(7,4) * t258 - Ifges(7,5) * t307 + t308 * t229;
t357 = mrSges(6,2) * t201 - mrSges(6,3) * t203 + Ifges(6,1) * t259 + Ifges(6,4) * t307 + Ifges(6,5) * t258 - qJ(6) * t186 + t329 * t230 + t366;
t165 = t357 + t390 * t308 + (-t234 + t232) * t329 + Ifges(5,5) * t307 - qJ(5) * t185 - mrSges(5,3) * t205 + mrSges(5,2) * t212 - Ifges(5,4) * t258 + Ifges(5,1) * t259;
t263 = mrSges(6,1) * t308 - mrSges(6,3) * t309;
t389 = -mrSges(5,1) * t308 - mrSges(5,2) * t309 - t263;
t395 = -mrSges(5,3) - mrSges(6,2);
t402 = -m(6) * t201 + t307 * mrSges(6,1) + t329 * t268;
t177 = m(5) * t205 + (t266 + t267) * t329 + t389 * t309 + (mrSges(5,1) + mrSges(7,1)) * t307 + t395 * t259 + t374 + t402;
t403 = m(6) * t199 + t307 * mrSges(6,3) + t329 * t272;
t178 = m(5) * t206 + (t270 - t271) * t329 + t389 * t308 + (-mrSges(5,2) + mrSges(7,2)) * t307 + t395 * t258 + t381 + t403;
t171 = t177 * t346 + t178 * t343;
t285 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t344 - Ifges(4,6) * t347) * qJD(1);
t362 = -mrSges(4,2) * t222 + mrSges(4,3) * t220 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t313 + Ifges(4,5) * t314 + pkin(8) * t171 + t343 * t163 - t346 * t165 - t285 * t384;
t365 = -m(4) * t222 - t313 * mrSges(4,1) - t171;
t284 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t344 - Ifges(4,3) * t347) * qJD(1);
t386 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t344 + Ifges(3,2) * t347) * qJD(1) - t284;
t406 = (-t347 * t283 + t386 * t344) * qJD(1) + mrSges(3,1) * t273 - mrSges(3,2) * t274 + Ifges(3,5) * t313 + Ifges(3,6) * t314 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t320 - t311 * t385 + t365) + qJ(3) * (t314 * mrSges(4,1) + t354) - t362;
t397 = pkin(7) * t350;
t394 = Ifges(3,4) + Ifges(4,6);
t392 = t329 * t232;
t172 = -t343 * t177 + t346 * t178;
t286 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t344 - Ifges(4,5) * t347) * qJD(1);
t387 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t344 + Ifges(3,6) * t347) * qJD(1) + t286;
t312 = (-mrSges(3,1) * t347 + mrSges(3,2) * t344) * qJD(1);
t319 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t384;
t168 = m(3) * t273 - mrSges(3,3) * t313 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t319 - t320) * qJD(2) + (-t311 - t312) * t385 + t365;
t318 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t385;
t174 = t354 + t312 * t384 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t314 - qJD(2) * t318 + m(3) * t274;
t376 = -t168 * t344 + t347 * t174;
t214 = -pkin(2) * t314 + t360 - t397;
t370 = -m(4) * t214 - t314 * mrSges(4,2) + t321 * t385 - t172;
t369 = mrSges(7,1) * t190 - mrSges(7,2) * t194 + Ifges(7,5) * t259 + Ifges(7,6) * t258 - Ifges(7,3) * t307 + t309 * t232 - t308 * t235;
t169 = -mrSges(4,3) * t313 + t320 * t384 - t370;
t287 = t371 - t397;
t361 = -mrSges(4,1) * t220 + mrSges(4,2) * t214 - pkin(3) * t407 - pkin(8) * t172 - t346 * t163 - t343 * t165;
t157 = -mrSges(3,1) * t287 + mrSges(3,3) * t274 - pkin(2) * t169 + (Ifges(3,2) + Ifges(4,3)) * t314 + t394 * t313 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t283 - t285) * qJD(2) - t387 * t385 + t361;
t356 = mrSges(6,1) * t201 - mrSges(6,3) * t199 - Ifges(6,4) * t259 - Ifges(6,2) * t307 - Ifges(6,6) * t258 + pkin(5) * t186 + t309 * t230 - t308 * t236 + t369;
t352 = mrSges(5,2) * t206 - t308 * t237 - pkin(4) * (-mrSges(6,2) * t259 - t263 * t309 - t186 + t402) - qJ(5) * (-mrSges(6,2) * t258 - t263 * t308 + t187 + t403) - mrSges(5,1) * t205 - t309 * t234 + Ifges(5,6) * t258 - Ifges(5,5) * t259 - Ifges(5,3) * t307 + t356;
t351 = -mrSges(4,1) * t222 + mrSges(4,3) * t214 - pkin(3) * t171 + t352;
t159 = -t351 + t387 * t384 + t394 * t314 + (Ifges(3,1) + Ifges(4,2)) * t313 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t386 * qJD(2) - qJ(3) * t169 - mrSges(3,3) * t273 + mrSges(3,2) * t287;
t355 = -m(3) * t287 + t319 * t384 + t314 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t313 + (-t318 * t344 - t320 * t347) * qJD(1) + t370;
t364 = mrSges(2,1) * t324 - mrSges(2,2) * t325 + Ifges(2,3) * qJDD(1) + pkin(1) * t355 + pkin(7) * t376 + t347 * t157 + t344 * t159;
t166 = m(2) * t324 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t350 + t355;
t162 = t168 * t347 + t174 * t344;
t160 = m(2) * t325 - mrSges(2,1) * t350 - qJDD(1) * mrSges(2,2) + t376;
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t325 + t350 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t406;
t154 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t350 - pkin(7) * t162 - t157 * t344 + t159 * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t348 * t154 - t345 * t155 - pkin(6) * (t160 * t345 + t166 * t348), t154, t159, -t284 * t385 - t362, t165, -t308 * t233 + t357 + t392, t366 + t392; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t345 * t154 + t348 * t155 + pkin(6) * (t160 * t348 - t166 * t345), t155, t157, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t313 - Ifges(4,6) * t314 - qJD(2) * t284 - t286 * t384 + t351, t163, -t356, -t309 * t229 - t367; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t364, t364, t406, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t313 - Ifges(4,3) * t314 + qJD(2) * t285 + t286 * t385 - t361, -t352, t358 - t329 * t236 + Ifges(6,6) * t307 + (-t229 + t233) * t309 + Ifges(6,3) * t258 + Ifges(6,5) * t259, t369;];
m_new  = t1;
