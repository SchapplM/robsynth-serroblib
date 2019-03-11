% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:16
% EndTime: 2019-03-09 02:00:23
% DurationCPUTime: 4.43s
% Computational Cost: add. (9866->392), mult. (19207->529), div. (0->0), fcn. (19965->8), ass. (0->217)
t420 = mrSges(7,2) + mrSges(6,3);
t234 = sin(qJ(5));
t229 = t234 ^ 2;
t236 = cos(qJ(5));
t230 = t236 ^ 2;
t424 = t229 + t230;
t426 = t420 * t424;
t419 = Ifges(6,5) + Ifges(7,4);
t425 = -Ifges(6,6) + Ifges(7,6);
t231 = sin(pkin(10));
t232 = cos(pkin(10));
t235 = sin(qJ(4));
t381 = cos(qJ(4));
t210 = t231 * t235 - t381 * t232;
t212 = t381 * t231 + t235 * t232;
t343 = t212 * t234;
t328 = mrSges(6,3) * t343;
t264 = -t210 * mrSges(6,2) - t328;
t207 = t210 * mrSges(7,3);
t325 = mrSges(7,2) * t343;
t297 = -t207 + t325;
t423 = t264 - t297;
t354 = t236 * mrSges(7,3);
t359 = t234 * mrSges(7,1);
t292 = -t354 + t359;
t151 = t292 * t212;
t422 = t151 / 0.2e1;
t421 = -mrSges(6,1) - mrSges(7,1);
t418 = Ifges(7,2) + Ifges(6,3);
t375 = t210 * pkin(5);
t269 = -cos(pkin(9)) * pkin(1) - pkin(3) * t232 - pkin(2);
t376 = pkin(8) * t212;
t129 = pkin(4) * t210 + t269 - t376;
t296 = sin(pkin(9)) * pkin(1) + qJ(3);
t268 = pkin(7) + t296;
t206 = t268 * t232;
t257 = t231 * t268;
t136 = t381 * t206 - t235 * t257;
t63 = t129 * t236 - t234 * t136;
t49 = -t63 - t375;
t369 = t49 + t63;
t295 = t236 * mrSges(6,1) - t234 * mrSges(6,2);
t417 = -mrSges(5,1) - t295;
t344 = t210 * t234;
t154 = -mrSges(6,2) * t212 + mrSges(6,3) * t344;
t157 = mrSges(7,2) * t344 + t212 * mrSges(7,3);
t416 = t154 + t157;
t225 = Ifges(6,4) * t236;
t221 = Ifges(6,1) * t234 + t225;
t413 = -Ifges(6,2) * t234 + t225;
t415 = t221 + t413;
t224 = Ifges(7,5) * t234;
t414 = Ifges(7,1) * t236 + t224;
t365 = Ifges(6,4) * t234;
t291 = Ifges(6,1) * t236 - t365;
t260 = t291 * t212;
t411 = t210 * t419 + t212 * t414 + t260;
t410 = t425 * t234 + t236 * t419;
t409 = -t424 * t376 / 0.2e1;
t282 = Ifges(7,3) * t236 - t224;
t408 = t414 - t282;
t219 = Ifges(6,2) * t236 + t365;
t364 = Ifges(7,5) * t236;
t220 = Ifges(7,1) * t234 - t364;
t383 = t236 / 0.2e1;
t387 = -t234 / 0.2e1;
t407 = t219 * t387 + t220 * t383;
t135 = t235 * t206 + t381 * t257;
t377 = pkin(8) * t210;
t378 = pkin(4) * t212;
t158 = t377 + t378;
t68 = t234 * t135 + t158 * t236;
t69 = -t236 * t135 + t234 * t158;
t274 = -t234 * t68 + t236 * t69;
t59 = qJ(6) * t212 + t69;
t60 = -t212 * pkin(5) - t68;
t276 = t234 * t60 + t236 * t59;
t321 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t342 = t212 * t236;
t327 = mrSges(6,3) * t342;
t266 = t210 * mrSges(6,1) - t327;
t405 = t234 * t423 + t236 * t266;
t223 = m(7) * qJ(6) + mrSges(7,3);
t349 = qJ(6) * t236;
t374 = t234 * pkin(5);
t280 = -t349 + t374;
t71 = t280 * t212 + t135;
t404 = -m(7) * t71 - t151;
t403 = t212 ^ 2;
t402 = 2 * qJD(4);
t401 = m(6) / 0.2e1;
t400 = m(6) / 0.4e1;
t399 = m(7) / 0.2e1;
t120 = t234 * t129;
t122 = t236 * t136;
t64 = t122 + t120;
t395 = t64 / 0.2e1;
t394 = m(7) * t60;
t393 = -t210 / 0.2e1;
t390 = t212 / 0.2e1;
t293 = t236 * mrSges(7,1) + t234 * mrSges(7,3);
t389 = -t293 / 0.2e1;
t388 = -t220 / 0.2e1;
t386 = t234 / 0.2e1;
t384 = -t236 / 0.2e1;
t380 = m(7) * t212;
t379 = m(7) * t280;
t373 = t60 * mrSges(7,1);
t372 = t68 * mrSges(6,1);
t371 = t69 * mrSges(6,2);
t346 = t210 * qJ(6);
t48 = t64 + t346;
t370 = t48 - t64;
t368 = m(7) * qJD(6);
t361 = t212 * mrSges(7,1);
t360 = t234 * mrSges(6,1);
t355 = t236 * mrSges(6,2);
t262 = t295 * t212;
t252 = -t262 / 0.2e1;
t324 = mrSges(7,2) * t342;
t265 = -t210 * mrSges(7,1) + t324;
t253 = t236 * t265;
t261 = t293 * t212;
t7 = -t212 * t253 / 0.2e1 + t210 * t252 + t261 * t393 - ((t369 + t375) * t236 + (t346 - t370) * t234) * t380 / 0.2e1 + t405 * t390 + t403 * t426 / 0.2e1;
t351 = t7 * qJD(1);
t244 = t280 * t399 - t354 / 0.2e1;
t250 = m(7) * (t369 * t234 + t370 * t236);
t294 = t355 + t360;
t256 = t294 + t359;
t340 = t236 * t207;
t312 = -t340 / 0.2e1;
t9 = t312 - t250 / 0.2e1 + (t244 + t256) * t210;
t350 = t9 * qJD(1);
t22 = t151 * t342 + t210 * t297 - m(7) * (t48 * t210 - t71 * t342);
t348 = qJD(1) * t22;
t208 = t210 * mrSges(5,2);
t337 = t424 * t377;
t281 = pkin(5) * t236 + t234 * qJ(6);
t213 = -pkin(4) - t281;
t341 = t213 * t212;
t267 = (-t337 - t378) * t401 + (-t337 + t341) * t399;
t241 = t267 + t420 * t210 * (-t229 / 0.2e1 - t230 / 0.2e1);
t243 = -m(6) * (t234 * t69 + t236 * t68) / 0.2e1 - m(7) * (t234 * t59 - t236 * t60) / 0.2e1;
t306 = -t295 / 0.2e1 + t389;
t307 = -t157 / 0.2e1 - t154 / 0.2e1;
t339 = t236 * t210;
t155 = t212 * mrSges(6,1) + mrSges(6,3) * t339;
t326 = mrSges(7,2) * t339;
t156 = -t326 - t361;
t308 = -t155 / 0.2e1 + t156 / 0.2e1;
t11 = t208 + t308 * t236 + t307 * t234 + (-mrSges(5,1) + t306) * t212 + t241 + t243;
t347 = t11 * qJD(1);
t345 = t210 * t212;
t299 = (t400 + m(7) / 0.4e1) * (0.1e1 - t424) * t345;
t24 = 0.2e1 * t299;
t338 = t24 * qJD(1);
t335 = qJD(4) * t210;
t334 = qJD(5) * t234;
t333 = qJD(5) * t236;
t146 = m(7) * t344;
t332 = t146 * qJD(1);
t331 = t146 * qJD(6);
t322 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t317 = -t342 / 0.2e1;
t310 = t293 * t384;
t300 = pkin(5) * mrSges(7,2) - t419;
t298 = -qJ(6) * mrSges(7,2) + t425;
t286 = Ifges(7,4) * t234 - Ifges(7,6) * t236;
t284 = Ifges(6,5) * t234 + Ifges(6,6) * t236;
t283 = Ifges(7,3) * t234 + t364;
t110 = Ifges(7,6) * t212 - t283 * t210;
t111 = Ifges(6,6) * t212 - t210 * t413;
t112 = Ifges(7,4) * t212 - t210 * t414;
t113 = Ifges(6,5) * t212 - t291 * t210;
t149 = t210 * t292;
t150 = t210 * t294;
t152 = t294 * t212;
t72 = -pkin(5) * t344 + qJ(6) * t339 + t136;
t1 = -t71 * t149 - t135 * t150 + t72 * t151 + t136 * t152 + t64 * t154 + t63 * t155 + t49 * t156 + t48 * t157 - t269 * t208 + t59 * t207 + m(6) * (t135 * t136 + t63 * t68 + t64 * t69) + m(7) * (t48 * t59 + t49 * t60 + t71 * t72) + (t269 * mrSges(5,1) - Ifges(5,4) * t212 + (t112 / 0.2e1 + t113 / 0.2e1 + t60 * mrSges(7,2) - t68 * mrSges(6,3) + t322 * t212) * t236 + (t110 / 0.2e1 - t111 / 0.2e1 - t59 * mrSges(7,2) - t69 * mrSges(6,3) - t321 * t212) * t234) * t212 + (-t371 - t373 + t372 + (-Ifges(5,1) + Ifges(5,2) + (-Ifges(6,1) / 0.2e1 - Ifges(7,1) / 0.2e1) * t230 + ((-Ifges(6,2) / 0.2e1 - Ifges(7,3) / 0.2e1) * t234 + (Ifges(6,4) - Ifges(7,5)) * t236) * t234 + t418) * t212 + (Ifges(5,4) - t410) * t210) * t210;
t249 = t355 / 0.2e1 + t360 / 0.2e1 + t359 / 0.2e1;
t275 = t234 * t63 - t236 * t64;
t277 = -t234 * t49 - t236 * t48;
t6 = (t312 - t150 / 0.2e1 - t149 / 0.2e1 + t249 * t210 + (t136 + t275) * t401 + (t277 + t72) * t399) * t210 + (t152 / 0.2e1 + t422 - t307 * t236 + t308 * t234 + (t135 + t274) * t401 + (t276 + t71) * t399) * t212;
t279 = t1 * qJD(1) + t6 * qJD(2);
t259 = t283 * t212;
t245 = Ifges(7,6) * t210 + t259;
t246 = Ifges(6,6) * t210 + t212 * t413;
t258 = t281 * t212;
t5 = t246 * t342 / 0.2e1 + t245 * t317 - t71 * t261 + t48 * t324 + t49 * t325 - t135 * t262 + t407 * t403 + (t284 / 0.2e1 + t286 / 0.2e1) * t345 - (-t236 * t221 + t234 * t282) * t403 / 0.2e1 + t411 * t343 / 0.2e1 + t404 * t258 + (-m(7) * t49 - t265 + t266 + t327) * t64 + (-m(7) * t48 - t328 - t423) * t63;
t278 = -t5 * qJD(1) - t7 * qJD(2);
t23 = 0.4e1 * t299;
t273 = t6 * qJD(1) + t23 * qJD(2);
t8 = (mrSges(5,3) * t212 + t152 + 0.4e1 * (t400 + m(5) / 0.4e1) * t135 - t404) * t212 + (-t340 + (mrSges(5,3) + t256) * t210 + m(7) * t277 + m(6) * t275 - m(5) * t136) * t210 + (m(4) * t296 + mrSges(4,3)) * (t231 ^ 2 + t232 ^ 2);
t272 = qJD(1) * t8 + qJD(2) * t24;
t159 = (m(7) * t213 - t293) * t234;
t242 = (-t234 * t71 + (-t341 + t377) * t236) * t399 + t151 * t387;
t21 = t326 + (-t310 + mrSges(7,1) / 0.2e1) * t212 - t394 / 0.2e1 + t242;
t271 = -qJD(1) * t21 + qJD(4) * t159;
t27 = -t207 + 0.2e1 * (t64 / 0.4e1 - t346 / 0.2e1 - t120 / 0.4e1 - t122 / 0.4e1) * m(7);
t270 = qJD(1) * t27 - qJD(5) * t223;
t237 = t71 * t292 / 0.2e1 + t135 * t294 / 0.2e1 + t213 * t261 / 0.2e1 + pkin(4) * t252 + t219 * t317 + t258 * t389 + (t213 * t258 + t280 * t71) * t399 - t234 * t246 / 0.4e1 + t280 * t422 + t410 * t210 / 0.4e1 + (t259 + t245) * t234 / 0.4e1 + t409 * mrSges(6,3) + (t260 + t411) * t236 / 0.4e1 + (-t221 / 0.4e1 + t388 - t415 / 0.4e1) * t343 + (-t282 / 0.4e1 + t408 / 0.4e1) * t342 + (t253 / 0.2e1 + (-t370 * t234 + t369 * t236) * t399 - t405 / 0.2e1) * pkin(8) + ((-t48 / 0.2e1 + t395) * t234 + t369 * t383 + t409) * mrSges(7,2);
t238 = (-pkin(5) * t60 + qJ(6) * t59) * t399 - pkin(5) * t156 / 0.2e1 + qJ(6) * t157 / 0.2e1 + t59 * mrSges(7,3) / 0.2e1 - t373 / 0.2e1 + t372 / 0.2e1 - t371 / 0.2e1;
t3 = -t237 + t238 + t418 * t390 + t321 * t344 - t419 * t339 / 0.2e1;
t32 = (t374 / 0.2e1 - t349 / 0.2e1 - t280 / 0.2e1) * m(7) * t210;
t36 = -t280 * t293 + t283 * t384 - pkin(4) * t294 + (t292 + t379) * t213 + t415 * t383 + (t291 + t408) * t386 + t407;
t251 = t3 * qJD(1) + t32 * qJD(2) - t36 * qJD(4);
t240 = (-m(7) * t281 - t293 - t295) * qJD(5);
t239 = (t244 + t249) * t210;
t214 = (m(7) * pkin(8) + mrSges(7,2)) * t236;
t33 = -t379 * t393 + t239 + (t292 + t294) * t210 / 0.2e1;
t25 = (t64 + 0.2e1 * t346) * t399 + m(7) * t395 - t297;
t20 = -t212 * t310 + t394 / 0.2e1 - t361 / 0.2e1 + t242;
t12 = t155 * t383 + t156 * t384 + t306 * t212 + t416 * t386 + t241 - t243;
t10 = t297 * t384 + t250 / 0.2e1 + t266 * t387 + t265 * t386 + t264 * t383 + t239;
t4 = (t321 * t234 - t322 * t236) * t210 + t237 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t212 + t238;
t2 = qJD(3) * t24 + qJD(4) * t6 - qJD(5) * t7;
t13 = [qJD(3) * t8 + qJD(4) * t1 - qJD(5) * t5 - qJD(6) * t22, t2, qJD(4) * t12 + qJD(5) * t10 + t272, t12 * qJD(3) + t4 * qJD(5) + t20 * qJD(6) + ((pkin(8) * t276 + t213 * t72) * t399 + (-pkin(4) * t136 + pkin(8) * t274) * t401) * t402 + (-Ifges(5,5) + (t388 - t221 / 0.2e1) * t236 + (t282 / 0.2e1 + t219 / 0.2e1) * t234) * t335 + t279 + (t135 * mrSges(5,2) - Ifges(5,6) * t212 + pkin(4) * t150 + t110 * t384 + t111 * t383 - t213 * t149 - t72 * t293 + (t416 * t236 + (-t155 + t156) * t234) * pkin(8) + (t284 + t286) * t390 + (t112 + t113) * t386 + t417 * t136 + t274 * mrSges(6,3) + t276 * mrSges(7,2)) * qJD(4), t10 * qJD(3) + t4 * qJD(4) + t25 * qJD(6) + t278 + ((-m(7) * pkin(5) + t421) * t64 + (-mrSges(6,2) + t223) * t63 + (t234 * t300 + t236 * t298) * t212) * qJD(5), qJD(4) * t20 + qJD(5) * t25 - t348; t2, t23 * qJD(4), t338, t33 * qJD(5) + t267 * t402 + t273 - t331 - t335 * t426 + (t208 + (-t293 + t417) * t212) * qJD(4), -t351 + t33 * qJD(4) + (t236 * t368 + t240) * t212, -t146 * qJD(4) + t333 * t380; -qJD(4) * t11 - qJD(5) * t9 - t272 + t331, -t338, 0, -t347, -t350 + (t354 - t355 - t379) * qJD(5) + (t421 * qJD(5) + t368) * t234, m(7) * t334 + t332; qJD(3) * t11 - qJD(5) * t3 + qJD(6) * t21 - t279, -qJD(5) * t32 - t273, t347, qJD(5) * t36 - qJD(6) * t159, pkin(8) * t240 + t214 * qJD(6) + t298 * t334 - t300 * t333 - t251, qJD(5) * t214 - t271; qJD(3) * t9 + qJD(4) * t3 - qJD(6) * t27 - t278, qJD(4) * t32 + t351, t350, t251, t223 * qJD(6), -t270; -qJD(3) * t146 - qJD(4) * t21 + qJD(5) * t27 + t348, 0, -t332, t271, t270, 0;];
Cq  = t13;
