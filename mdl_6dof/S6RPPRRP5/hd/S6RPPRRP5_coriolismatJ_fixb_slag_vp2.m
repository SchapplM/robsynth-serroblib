% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:43
% EndTime: 2019-03-09 02:07:51
% DurationCPUTime: 4.21s
% Computational Cost: add. (4609->456), mult. (9289->621), div. (0->0), fcn. (6933->4), ass. (0->227)
t391 = -m(7) / 0.2e1;
t256 = cos(qJ(5));
t411 = Ifges(6,6) + Ifges(7,6);
t418 = t256 * t411;
t254 = sin(qJ(5));
t248 = t254 ^ 2;
t249 = t256 ^ 2;
t331 = t248 + t249;
t366 = -qJ(6) - pkin(8);
t192 = t366 * t254;
t195 = t366 * t256;
t257 = cos(qJ(4));
t255 = sin(qJ(4));
t338 = t256 * t257;
t185 = t255 * mrSges(7,1) - mrSges(7,3) * t338;
t344 = t254 * t257;
t181 = -t255 * mrSges(7,2) - mrSges(7,3) * t344;
t340 = t256 * t181;
t381 = t254 / 0.2e1;
t335 = t185 * t381 - t340 / 0.2e1;
t251 = pkin(1) + qJ(3);
t371 = pkin(4) * t255;
t187 = -pkin(8) * t257 + t251 + t371;
t156 = t256 * t187;
t250 = qJ(2) - pkin(7);
t301 = -t250 * t254 + pkin(5);
t321 = qJ(6) * t338;
t47 = t255 * t301 + t156 - t321;
t342 = t255 * t256;
t62 = t187 * t254 + t250 * t342;
t52 = -qJ(6) * t344 + t62;
t417 = t335 + ((-t192 * t257 + t52) * t256 + (t195 * t257 - t47) * t254) * t391;
t343 = t255 * t254;
t61 = -t250 * t343 + t156;
t51 = t61 - t321;
t365 = t47 - t51;
t413 = m(7) * t365;
t416 = -t413 / 0.2e1;
t415 = mrSges(7,3) + mrSges(6,3);
t369 = t254 * pkin(5);
t300 = -t250 + t369;
t172 = t300 * t257;
t375 = m(7) * t172;
t414 = m(7) * (t254 * t52 + t256 * t47);
t412 = Ifges(6,5) + Ifges(7,5);
t410 = Ifges(6,3) + Ifges(7,3);
t389 = m(7) * pkin(5);
t409 = -mrSges(7,1) - t389;
t238 = t254 * mrSges(7,1);
t241 = t256 * mrSges(7,2);
t404 = t241 + t238;
t149 = t404 * t257;
t403 = t254 * mrSges(6,1) + t256 * mrSges(6,2);
t150 = t403 * t257;
t408 = t149 + t150;
t179 = -t257 * mrSges(7,2) + mrSges(7,3) * t343;
t180 = -t257 * mrSges(6,2) + mrSges(6,3) * t343;
t407 = t179 + t180;
t182 = -t255 * mrSges(6,2) - mrSges(6,3) * t344;
t406 = t181 + t182;
t183 = t257 * mrSges(7,1) + mrSges(7,3) * t342;
t184 = t257 * mrSges(6,1) + mrSges(6,3) * t342;
t405 = t183 + t184;
t186 = t255 * mrSges(6,1) - mrSges(6,3) * t338;
t333 = t185 + t186;
t247 = pkin(4) * t257;
t203 = pkin(8) * t255 + t247;
t175 = t256 * t203;
t77 = -t250 * t344 + t175;
t78 = t254 * t203 + t250 * t338;
t282 = -t254 * t77 + t256 * t78;
t48 = qJ(6) * t342 + t257 * t301 + t175;
t56 = qJ(6) * t343 + t78;
t401 = -t254 * t48 + t256 * t56;
t245 = Ifges(7,4) * t256;
t400 = -t254 * Ifges(7,1) - t245;
t246 = Ifges(6,4) * t256;
t399 = -t254 * Ifges(6,1) - t246;
t398 = -t241 / 0.2e1 - t238 / 0.2e1;
t286 = -Ifges(7,2) * t254 + t245;
t105 = t255 * Ifges(7,6) + t257 * t286;
t287 = -Ifges(6,2) * t254 + t246;
t107 = t255 * Ifges(6,6) + t257 * t287;
t324 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t397 = t324 * t255 + t105 / 0.2e1 + t107 / 0.2e1;
t396 = t255 ^ 2;
t395 = t257 ^ 2;
t394 = -m(5) / 0.4e1;
t393 = -m(6) / 0.2e1;
t392 = m(6) / 0.2e1;
t390 = m(7) / 0.2e1;
t388 = -t48 / 0.2e1;
t145 = mrSges(7,1) * t338 - mrSges(7,2) * t344;
t387 = t145 / 0.2e1;
t386 = t182 / 0.2e1;
t385 = t183 / 0.2e1;
t384 = -t186 / 0.2e1;
t383 = t192 / 0.2e1;
t380 = t255 / 0.2e1;
t378 = t256 / 0.2e1;
t377 = -t257 / 0.2e1;
t376 = t257 / 0.2e1;
t368 = t256 * pkin(5);
t228 = -pkin(4) - t368;
t374 = m(7) * t228;
t373 = m(7) * t255;
t372 = m(7) * t257;
t370 = pkin(5) * t395;
t364 = Ifges(6,4) * t254;
t363 = Ifges(7,4) * t254;
t240 = t255 * mrSges(5,2);
t304 = t249 / 0.2e1 + t248 / 0.2e1;
t268 = t415 * t304;
t280 = -t192 * t254 - t195 * t256;
t302 = pkin(8) * t331;
t261 = t268 * t255 + (t255 * t302 + t247) * t392 + (-t228 * t257 + t255 * t280) * t390;
t263 = (-t254 * t78 - t256 * t77) * t393 + (-t254 * t56 - t256 * t48) * t391;
t193 = -t256 * mrSges(7,1) + t254 * mrSges(7,2);
t357 = t256 * mrSges(6,1);
t292 = -t254 * mrSges(6,2) + t357;
t308 = t292 / 0.2e1 - t193 / 0.2e1;
t310 = t385 + t184 / 0.2e1;
t312 = t179 / 0.2e1 + t180 / 0.2e1;
t9 = -t240 + t310 * t256 + t312 * t254 + (mrSges(5,1) + t308) * t257 + t261 + t263;
t362 = qJD(1) * t9;
t104 = Ifges(7,6) * t257 - t255 * t286;
t106 = Ifges(6,6) * t257 - t255 * t287;
t288 = Ifges(7,1) * t256 - t363;
t108 = Ifges(7,5) * t257 - t255 * t288;
t289 = Ifges(6,1) * t256 - t364;
t110 = Ifges(6,5) * t257 - t255 * t289;
t147 = t404 * t255;
t148 = t403 * t255;
t171 = t300 * t255;
t109 = t255 * Ifges(7,5) + t257 * t288;
t111 = t255 * Ifges(6,5) + t257 * t289;
t325 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t269 = -t109 / 0.2e1 - t111 / 0.2e1 - t325 * t255;
t3 = -t172 * t147 - t171 * t149 + t52 * t179 + t62 * t180 + t56 * t181 + t78 * t182 + t47 * t183 + t61 * t184 + t48 * t185 + t77 * t186 - t251 * t240 + m(7) * (-t171 * t172 + t47 * t48 + t52 * t56) + m(6) * (t61 * t77 + t62 * t78) + (t251 * mrSges(5,1) - Ifges(5,4) * t257 + t250 * t148 + (t108 / 0.2e1 + t110 / 0.2e1 + t325 * t257) * t256 + (-t104 / 0.2e1 - t106 / 0.2e1 - t324 * t257) * t254) * t257 + (Ifges(5,4) * t255 + t250 * t150 + t269 * t256 + t397 * t254 + (-m(6) * t250 ^ 2 - Ifges(5,1) + Ifges(5,2) + t410) * t257) * t255;
t354 = t3 * qJD(1);
t146 = t292 * t257;
t197 = t256 * Ifges(7,2) + t363;
t151 = t257 * t197;
t198 = t256 * Ifges(6,2) + t364;
t152 = t257 * t198;
t153 = t400 * t257;
t154 = t399 * t257;
t4 = t172 * t145 + t51 * t181 + t61 * t182 - t62 * t186 + (-t413 - t185) * t52 + (-t250 * t146 + (t47 * mrSges(7,3) + t61 * mrSges(6,3) + t152 / 0.2e1 + t151 / 0.2e1 + t269) * t254 + (-t52 * mrSges(7,3) - t62 * mrSges(6,3) + t153 / 0.2e1 + t154 / 0.2e1 + (t149 + t375) * pkin(5) - t397) * t256) * t257;
t353 = t4 * qJD(1);
t326 = -mrSges(6,2) / 0.2e1 - mrSges(7,2) / 0.2e1;
t295 = t326 * t254;
t309 = t185 / 0.2e1 + t186 / 0.2e1;
t311 = t181 / 0.2e1 + t386;
t327 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t7 = (t146 / 0.2e1 + t387) * t257 + t295 + (t311 * t254 + t268 * t257) * t255 + (t309 * t255 + (pkin(5) / 0.2e1 + t370 / 0.2e1 + t47 * t380 - t255 * t51 / 0.2e1) * m(7) + t327) * t256;
t352 = t7 * qJD(1);
t351 = -t292 - mrSges(5,1);
t283 = -t254 * t61 + t256 * t62;
t284 = -t254 * t47 + t256 * t52;
t296 = m(5) * (t395 + t396);
t339 = t256 * t182;
t11 = -mrSges(4,2) - mrSges(3,3) + (mrSges(5,3) * t257 + t375 + t408) * t257 + (-m(6) * t395 - t296) * t250 + (-m(4) - m(3)) * qJ(2) + (-m(6) * t283 - m(7) * t284 + mrSges(5,3) * t255 + t254 * t333 - t339 - t340) * t255;
t350 = qJD(1) * t11;
t15 = t255 * mrSges(5,1) + t257 * mrSges(5,2) + mrSges(4,3) + t333 * t256 + t406 * t254 + t414 + m(6) * (t254 * t62 + t256 * t61) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t251;
t349 = qJD(1) * t15;
t23 = (-t254 * t181 - t256 * t185 - t414) * t257;
t348 = qJD(1) * t23;
t329 = t389 / 0.2e1;
t278 = t329 + t327;
t12 = (t255 * t326 + t386) * t256 + (-t278 * t255 + t384 + t416) * t254 - t335;
t347 = t12 * qJD(1);
t346 = t250 * t257;
t341 = t255 * t257;
t328 = m(6) / 0.4e1 + m(7) / 0.4e1;
t274 = -m(5) / 0.2e1 - 0.2e1 * t328 * t331;
t303 = m(7) * t331;
t28 = -m(4) + 0.2e1 * (t394 - t328) * t395 + 0.2e1 * (-m(6) * t331 / 0.4e1 - t303 / 0.4e1 + t394) * t396 + t274;
t337 = t28 * qJD(1);
t75 = (0.1e1 / 0.2e1 + t304) * t372;
t336 = t75 * qJD(1);
t330 = qJD(4) * t255;
t155 = -m(7) * t369 - t404;
t320 = -t344 / 0.2e1;
t319 = -t343 / 0.2e1;
t317 = -t342 / 0.2e1;
t307 = t197 / 0.2e1 + t198 / 0.2e1;
t306 = t400 / 0.2e1 + t399 / 0.2e1;
t243 = Ifges(7,5) * t256;
t244 = Ifges(6,5) * t256;
t305 = t244 / 0.4e1 + t243 / 0.4e1;
t299 = t331 * mrSges(7,3);
t294 = t326 * t256;
t293 = t303 / 0.2e1;
t35 = 0.4e1 * t328 * (-0.1e1 + t331) * t341;
t272 = (t192 * t256 - t195 * t254) * t390;
t273 = t282 - 0.2e1 * t346;
t275 = t172 + t401;
t276 = t171 + t284;
t5 = t272 + (-t148 / 0.2e1 - t147 / 0.2e1 - t311 * t256 + t309 * t254 + t276 * t391 + t283 * t393) * t257 + (-t150 / 0.2e1 - t149 / 0.2e1 - t312 * t256 + t310 * t254 + t275 * t391 + t273 * t393) * t255;
t281 = -t5 * qJD(1) + t35 * qJD(3);
t72 = -t338 * t389 - t145;
t279 = qJD(1) * t72 + qJD(4) * t155;
t14 = pkin(4) * t403 - t228 * t404 - t193 * t369 + (-t245 / 0.2e1 - t246 / 0.2e1 + t306) * t256 + (-pkin(5) * t374 + (Ifges(6,4) / 0.2e1 + Ifges(7,4) / 0.2e1) * t254 + (Ifges(7,2) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1) * t256 + t307) * t254;
t259 = -t250 * t403 / 0.2e1 - t304 * pkin(8) * mrSges(6,3) + (mrSges(7,3) * t383 + t399 / 0.4e1 + t400 / 0.4e1 - t246 / 0.4e1 - t245 / 0.4e1 + (Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.4e1) * t254) * t254 + (t195 * mrSges(7,3) / 0.2e1 - t198 / 0.4e1 - t197 / 0.4e1 + (Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1) * t256 + (-Ifges(6,4) / 0.4e1 - Ifges(7,4) / 0.4e1) * t254 + (t193 / 0.2e1 + t374 / 0.2e1) * pkin(5)) * t256;
t260 = -(t416 - t185 / 0.2e1) * t195 + (t111 / 0.4e1 + t109 / 0.4e1 - t152 / 0.4e1 - t151 / 0.4e1 + pkin(8) * t384 + (t51 / 0.2e1 - t47 / 0.2e1) * mrSges(7,3)) * t256 - pkin(4) * t146 / 0.2e1 + t172 * t404 / 0.2e1 + t181 * t383 + t228 * t387;
t262 = t154 / 0.4e1 + t153 / 0.4e1 - t107 / 0.4e1 - t105 / 0.4e1 - pkin(8) * t182 / 0.2e1 + (t375 / 0.2e1 + t149 / 0.2e1) * pkin(5);
t266 = mrSges(7,1) * t388 + t56 * mrSges(7,2) / 0.2e1 - t77 * mrSges(6,1) / 0.2e1 + t78 * mrSges(6,2) / 0.2e1;
t2 = (-t183 / 0.2e1 + m(7) * t388) * pkin(5) + t262 * t254 + (t325 * t256 + (-0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t254 + t305) * t255 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t259) * t257 + t260 + t266;
t25 = (t403 / 0.2e1 + t294 - t327 * t254 - t398) * t257;
t271 = t2 * qJD(1) - t25 * qJD(3) - t14 * qJD(4);
t267 = -t171 * t390 + t255 * t398;
t18 = t267 + t417;
t41 = m(7) * t280 + t299;
t74 = (-0.1e1 / 0.2e1 + t304) * t373;
t270 = -qJD(1) * t18 + qJD(3) * t74 + qJD(4) * t41;
t264 = -t254 * t278 + t294;
t76 = -t372 / 0.2e1 + t257 * t293;
t73 = t373 / 0.2e1 + t255 * t293;
t29 = t296 / 0.2e1 + t274 + (t392 + t390) * (t331 * t396 + t395);
t26 = t264 * t257 + t320 * t389 + (t403 + t404) * t377;
t19 = t267 - t417;
t13 = -t339 / 0.2e1 + t264 * t255 + t335 + (t413 + t186) * t381;
t10 = t308 * t257 + t261 - t263 - t407 * t254 / 0.2e1 - t405 * t256 / 0.2e1;
t8 = m(7) * (-t255 * t365 - t370) * t378 + t295 + t278 * t256 + (t146 + t145) * t377 + t406 * t319 + t333 * t317 - t415 * t331 * t341 / 0.2e1;
t6 = (t255 * t275 + t257 * t276) * t390 + (t255 * t273 + t257 * t283) * t392 + t272 + t408 * t380 + (-t148 - t147) * t377 + t333 * t320 + t405 * t319 + t407 * t342 / 0.2e1 + t406 * t338 / 0.2e1;
t1 = t305 * t255 + t260 + ((-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t255 + t262) * t254 + t259 * t257 + pkin(5) * t385 + t48 * t329 - t266 + t410 * t376 + t411 * t343 / 0.2e1 + t412 * t317;
t16 = [-qJD(2) * t11 + qJD(3) * t15 + qJD(4) * t3 + qJD(5) * t4 + qJD(6) * t23, qJD(3) * t29 + qJD(4) * t10 + qJD(5) * t13 + qJD(6) * t76 - t350, qJD(2) * t29 + qJD(4) * t6 + qJD(5) * t8 + t349, t354 + t10 * qJD(2) + t6 * qJD(3) + t1 * qJD(5) + t19 * qJD(6) + (-Ifges(5,5) + t306 * t256 + t307 * t254 + (-m(6) * pkin(4) + t351) * t250) * t330 + (-Ifges(5,6) * t257 - t228 * t147 + t192 * t183 - t171 * t193 - t195 * t179 + pkin(4) * t148 + m(7) * (-t171 * t228 + t192 * t48 - t195 * t56) - mrSges(5,2) * t346 + (m(6) * t282 + t256 * t180 - t254 * t184) * pkin(8) + (t108 + t110) * t381 + (t104 + t106) * t378 + (t254 * t412 + t418) * t376 + t401 * mrSges(7,3) + t282 * mrSges(6,3)) * qJD(4), t13 * qJD(2) + t8 * qJD(3) + t1 * qJD(4) + t353 + (-mrSges(6,1) * t62 - mrSges(6,2) * t61 - mrSges(7,2) * t51 + (-t418 + (mrSges(7,3) * pkin(5) - t412) * t254) * t257 + t409 * t52) * qJD(5), qJD(2) * t76 + qJD(4) * t19 + t348; qJD(3) * t28 - qJD(4) * t9 - qJD(5) * t12 + qJD(6) * t75 + t350, 0, t337, -t362, -t347 + (-t155 + t403) * qJD(5), t336; -qJD(2) * t28 - qJD(4) * t5 - qJD(5) * t7 - t349, -t337, t35 * qJD(4), t26 * qJD(5) + t73 * qJD(6) + (t193 + t351) * t330 + (mrSges(6,3) * t331 - mrSges(5,2) + t299) * qJD(4) * t257 + 0.2e1 * ((t228 * t255 + t257 * t280) * t390 + (t257 * t302 - t371) * t392) * qJD(4) + t281, -t352 + t26 * qJD(4) + ((mrSges(6,2) + mrSges(7,2)) * t254 + (-mrSges(6,1) + t409) * t256) * qJD(5) * t255, t73 * qJD(4); qJD(2) * t9 + qJD(3) * t5 + qJD(5) * t2 - qJD(6) * t18 - t354, t362, -qJD(5) * t25 + qJD(6) * t74 - t281, -qJD(5) * t14 + qJD(6) * t41, t271 + (-mrSges(7,2) * t192 - mrSges(7,3) * t368 - pkin(8) * t357 + t243 + t244 + (mrSges(6,2) * pkin(8) - t411) * t254 - t409 * t195) * qJD(5), t270; qJD(2) * t12 + qJD(3) * t7 - qJD(4) * t2 + qJD(6) * t72 - t353, t347, qJD(4) * t25 + t352, t155 * qJD(6) - t271, 0, t279; -qJD(2) * t75 + qJD(4) * t18 - qJD(5) * t72 - t348, -t336, -t74 * qJD(4), -qJD(5) * t155 - t270, -t279, 0;];
Cq  = t16;
