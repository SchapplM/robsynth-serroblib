% Calculate vector of inverse dynamics joint torques for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:08
% EndTime: 2019-12-05 16:54:38
% DurationCPUTime: 13.52s
% Computational Cost: add. (3016->495), mult. (7211->678), div. (0->0), fcn. (5103->10), ass. (0->236)
t385 = Ifges(5,4) + Ifges(6,4);
t386 = Ifges(5,1) + Ifges(6,1);
t366 = -Ifges(6,5) - Ifges(5,5);
t384 = Ifges(5,2) + Ifges(6,2);
t383 = Ifges(5,6) + Ifges(6,6);
t168 = cos(qJ(4));
t390 = t385 * t168;
t165 = sin(qJ(4));
t389 = t385 * t165;
t270 = qJD(3) * t168;
t166 = sin(qJ(3));
t274 = qJD(2) * t166;
t128 = -t165 * t274 + t270;
t169 = cos(qJ(3));
t264 = qJD(2) * qJD(3);
t134 = qJDD(2) * t166 + t169 * t264;
t59 = qJD(4) * t128 + qJDD(3) * t165 + t134 * t168;
t328 = t59 / 0.2e1;
t129 = qJD(3) * t165 + t168 * t274;
t60 = -qJD(4) * t129 + qJDD(3) * t168 - t134 * t165;
t327 = t60 / 0.2e1;
t388 = -t128 / 0.2e1;
t322 = t129 / 0.2e1;
t272 = qJD(2) * t169;
t375 = t272 - qJD(4);
t320 = -t375 / 0.2e1;
t387 = qJD(3) / 0.2e1;
t382 = -Ifges(6,3) - Ifges(5,3);
t347 = -t165 * t383 - t168 * t366;
t346 = -t165 * t384 + t390;
t343 = t168 * t386 - t389;
t133 = qJDD(2) * t169 - t166 * t264;
t125 = qJDD(4) - t133;
t381 = t385 * t327 + t386 * t328 - t366 * t125 / 0.2e1;
t170 = cos(qJ(2));
t162 = sin(pkin(5));
t275 = qJD(2) * t162;
t241 = qJD(1) * t275;
t148 = t170 * t241;
t167 = sin(qJ(2));
t263 = qJDD(1) * t162;
t101 = t167 * t263 + t148;
t163 = cos(pkin(5));
t276 = qJD(1) * t163;
t380 = qJDD(2) * pkin(7) + qJD(3) * t276 + t101;
t161 = Ifges(4,4) * t272;
t378 = t385 * t128;
t362 = t386 * t129 + t366 * t375 + t378;
t379 = Ifges(4,1) * t274 + Ifges(4,5) * qJD(3) + t168 * t362 + t161;
t160 = pkin(4) * t168 + pkin(3);
t164 = -qJ(5) - pkin(8);
t209 = mrSges(4,1) * t169 - mrSges(4,2) * t166;
t377 = -m(6) * (t160 * t169 - t164 * t166) - t166 * mrSges(6,3) - mrSges(3,1) - t209;
t376 = t385 * t129;
t305 = Ifges(4,4) * t166;
t200 = Ifges(4,2) * t169 + t305;
t374 = t382 * t320 + t366 * t322 + t383 * t388 + Ifges(4,6) * t387 + qJD(2) * t200 / 0.2e1;
t329 = m(6) * pkin(4);
t373 = t383 * t125 + t384 * t60 + t385 * t59;
t371 = t133 / 0.2e1;
t370 = t134 / 0.2e1;
t369 = -mrSges(5,1) - mrSges(6,1);
t368 = mrSges(3,2) - mrSges(4,3);
t367 = mrSges(5,2) + mrSges(6,2);
t363 = t384 * t128 - t375 * t383 + t376;
t307 = mrSges(6,3) * t128;
t84 = mrSges(6,2) * t375 + t307;
t309 = mrSges(5,3) * t128;
t85 = mrSges(5,2) * t375 + t309;
t361 = t85 + t84;
t306 = mrSges(6,3) * t129;
t86 = -mrSges(6,1) * t375 - t306;
t308 = mrSges(5,3) * t129;
t87 = -mrSges(5,1) * t375 - t308;
t360 = t86 + t87;
t358 = mrSges(6,1) + t329;
t226 = qJD(4) * t164;
t248 = t165 * t272;
t265 = qJD(5) * t168;
t210 = pkin(3) * t166 - pkin(8) * t169;
t131 = t210 * qJD(2);
t277 = qJD(1) * t162;
t253 = t167 * t277;
t136 = qJD(2) * pkin(7) + t253;
t88 = -t166 * t136 + t169 * t276;
t46 = t165 * t131 + t168 * t88;
t357 = qJ(5) * t248 + t165 * t226 + t265 - t46;
t282 = t168 * t169;
t189 = pkin(4) * t166 - qJ(5) * t282;
t45 = t168 * t131 - t165 * t88;
t356 = -qJD(2) * t189 - qJD(5) * t165 + t168 * t226 - t45;
t13 = -mrSges(5,1) * t60 + mrSges(5,2) * t59;
t355 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t134 + t13;
t251 = t166 * t276;
t268 = qJD(4) * t165;
t315 = pkin(4) * t165;
t354 = -t251 - (qJD(2) * t315 + t136) * t169 + pkin(4) * t268;
t255 = mrSges(4,3) * t274;
t353 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t128 + mrSges(5,2) * t129 + t255;
t211 = pkin(3) * t169 + pkin(8) * t166;
t138 = -pkin(2) - t211;
t158 = pkin(7) * t282;
t97 = t165 * t138 + t158;
t352 = -t382 * t166 + t347 * t169;
t351 = t383 * t166 + t346 * t169;
t350 = -t366 * t166 + t343 * t169;
t348 = -t366 * t165 + t383 * t168;
t345 = t384 * t168 + t389;
t344 = t386 * t165 + t390;
t342 = -t382 * t125 - t366 * t59 + t383 * t60;
t262 = qJDD(1) * t163;
t271 = qJD(3) * t166;
t22 = -t136 * t271 + t166 * t262 + t169 * t380;
t269 = qJD(3) * t169;
t23 = -t136 * t269 - t166 * t380 + t169 * t262;
t341 = -t166 * t23 + t169 * t22;
t20 = qJDD(3) * pkin(8) + t22;
t266 = qJD(4) * t168;
t147 = t167 * t241;
t100 = t170 * t263 - t147;
t94 = -qJDD(2) * pkin(2) - t100;
t38 = -pkin(3) * t133 - pkin(8) * t134 + t94;
t89 = t136 * t169 + t251;
t78 = qJD(3) * pkin(8) + t89;
t252 = t170 * t277;
t91 = qJD(2) * t138 - t252;
t3 = t165 * t38 + t168 * t20 + t91 * t266 - t268 * t78;
t26 = t165 * t91 + t168 * t78;
t4 = -qJD(4) * t26 - t165 * t20 + t168 * t38;
t340 = -t165 * t4 + t168 * t3;
t339 = -m(4) - m(5) - m(6);
t338 = -mrSges(5,1) - t358;
t206 = -mrSges(6,1) * t168 + mrSges(6,2) * t165;
t208 = -mrSges(5,1) * t168 + mrSges(5,2) * t165;
t337 = m(5) * pkin(3) + m(6) * t160 + mrSges(4,1) - t206 - t208;
t336 = -m(5) * pkin(8) + m(6) * t164 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t77 = -qJD(3) * pkin(3) - t88;
t47 = -pkin(4) * t128 + qJD(5) + t77;
t63 = -mrSges(6,1) * t128 + mrSges(6,2) * t129;
t335 = -m(6) * t47 - t63;
t333 = m(5) * t211 + t166 * mrSges(5,3);
t332 = t333 - t377;
t331 = -m(5) * t77 - t353;
t2 = qJ(5) * t60 + qJD(5) * t128 + t3;
t330 = -t4 * mrSges(5,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t171 = qJD(2) ^ 2;
t326 = t125 / 0.2e1;
t324 = t128 / 0.2e1;
t304 = Ifges(4,4) * t169;
t21 = -qJDD(3) * pkin(3) - t23;
t295 = t166 * t21;
t292 = cos(pkin(9));
t291 = sin(pkin(9));
t220 = t291 * t170;
t223 = t292 * t167;
t107 = t163 * t223 + t220;
t290 = t107 * t165;
t221 = t291 * t167;
t222 = t292 * t170;
t109 = -t163 * t221 + t222;
t289 = t109 * t165;
t288 = t162 * t167;
t287 = t162 * t170;
t286 = t165 * t166;
t285 = t165 * t167;
t284 = t165 * t169;
t283 = t166 * t168;
t281 = t169 * t170;
t132 = t210 * qJD(3);
t280 = t165 * t132 + t138 * t266;
t258 = pkin(7) * t271;
t279 = t168 * t132 + t165 * t258;
t273 = qJD(2) * t167;
t267 = qJD(4) * t166;
t257 = pkin(7) * t269;
t254 = mrSges(4,3) * t272;
t250 = t162 * t273;
t249 = t170 * t275;
t247 = t165 * t269;
t12 = -t60 * mrSges(6,1) + t59 * mrSges(6,2);
t25 = -t165 * t78 + t168 * t91;
t225 = t162 * t292;
t224 = t162 * t291;
t219 = t264 / 0.2e1;
t207 = mrSges(5,1) * t165 + mrSges(5,2) * t168;
t205 = mrSges(6,1) * t165 + mrSges(6,2) * t168;
t195 = Ifges(4,5) * t169 - Ifges(4,6) * t166;
t17 = -qJ(5) * t129 + t25;
t112 = t163 * t166 + t169 * t288;
t71 = -t112 * t165 - t168 * t287;
t188 = -t112 * t168 + t165 * t287;
t111 = -t163 * t169 + t166 * t288;
t137 = -qJD(2) * pkin(2) - t252;
t185 = t137 * (mrSges(4,1) * t166 + mrSges(4,2) * t169);
t184 = t166 * (Ifges(4,1) * t169 - t305);
t181 = -t165 * t267 + t168 * t269;
t180 = t166 * t266 + t247;
t93 = (t168 * t281 + t285) * t162;
t92 = (-t165 * t281 + t167 * t168) * t162;
t142 = t164 * t168;
t141 = t164 * t165;
t140 = -qJD(3) * mrSges(4,2) + t254;
t135 = (pkin(7) + t315) * t166;
t130 = t209 * qJD(2);
t127 = t168 * t138;
t108 = t163 * t220 + t223;
t106 = -t163 * t222 + t221;
t102 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t133;
t96 = -pkin(7) * t284 + t127;
t90 = pkin(4) * t180 + t257;
t83 = qJD(1) * t93;
t82 = qJD(1) * t92;
t74 = -mrSges(4,1) * t133 + mrSges(4,2) * t134;
t73 = -qJ(5) * t286 + t97;
t70 = -qJD(3) * t111 + t169 * t249;
t69 = qJD(3) * t112 + t166 * t249;
t68 = t109 * t169 + t166 * t224;
t67 = t109 * t166 - t169 * t224;
t66 = t107 * t169 - t166 * t225;
t65 = t107 * t166 + t169 * t225;
t61 = -qJ(5) * t283 + t127 + (-pkin(7) * t165 - pkin(4)) * t169;
t40 = -qJD(4) * t97 + t279;
t39 = (-t166 * t270 - t169 * t268) * pkin(7) + t280;
t34 = -mrSges(5,2) * t125 + mrSges(5,3) * t60;
t33 = -mrSges(6,2) * t125 + mrSges(6,3) * t60;
t32 = mrSges(5,1) * t125 - mrSges(5,3) * t59;
t31 = mrSges(6,1) * t125 - mrSges(6,3) * t59;
t18 = qJ(5) * t128 + t26;
t16 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t283 + (-qJD(5) * t166 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t169) * t165 + t280;
t15 = qJD(4) * t71 + t165 * t250 + t168 * t70;
t14 = qJD(4) * t188 - t165 * t70 + t168 * t250;
t11 = -t166 * t265 + t189 * qJD(3) + (-t158 + (qJ(5) * t166 - t138) * t165) * qJD(4) + t279;
t10 = -pkin(4) * t375 + t17;
t5 = -pkin(4) * t60 + qJDD(5) + t21;
t1 = pkin(4) * t125 - qJ(5) * t59 - qJD(5) * t129 + t4;
t6 = [m(2) * qJDD(1) + t112 * t102 + t70 * t140 - (t33 + t34) * t188 + (t31 + t32) * t71 + t361 * t15 + t360 * t14 + (t63 + t353) * t69 + (t12 + t355) * t111 + (-m(2) - m(3) + t339) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t171 - t74) * t170 + (-mrSges(3,1) * t171 - mrSges(3,2) * qJDD(2) - qJD(2) * t130) * t167) * t162 + m(3) * (qJDD(1) * t163 ^ 2 + (t100 * t170 + t101 * t167) * t162) + m(4) * (-t111 * t23 + t112 * t22 - t69 * t88 + t70 * t89 + (t137 * t273 - t170 * t94) * t162) + m(5) * (t111 * t21 + t14 * t25 + t15 * t26 - t188 * t3 + t4 * t71 + t69 * t77) + m(6) * (t1 * t71 + t10 * t14 + t111 * t5 + t15 * t18 - t188 * t2 + t47 * t69); -t360 * t82 - t361 * t83 + qJDD(3) * (Ifges(4,5) * t166 + Ifges(4,6) * t169) - (t165 * t362 + t168 * t363) * t267 / 0.2e1 + (-t101 + t148) * mrSges(3,2) + t135 * t12 + t97 * t34 + t39 * t85 + t11 * t86 + t40 * t87 + t90 * t63 + t96 * t32 + t73 * t33 - pkin(2) * t74 + t16 * t84 + t61 * t31 + (t331 + t335) * t166 * t252 + (t100 + t147) * mrSges(3,1) + (-(t137 * t167 + (-t166 * t88 + t169 * t89) * t170) * t277 - pkin(2) * t94 + ((-t166 * t89 - t169 * t88) * qJD(3) + t341) * pkin(7)) * m(4) - t342 * t169 / 0.2e1 + t130 * t253 + (-t333 * t287 + t369 * t93 - t367 * t92 + t339 * (pkin(2) * t287 + pkin(7) * t288) + (t368 * t167 + t170 * t377 - t285 * t329) * t162) * g(3) - t94 * t209 + (t25 * mrSges(5,1) + t10 * mrSges(6,1) - t26 * mrSges(5,2) - t18 * mrSges(6,2) - t89 * mrSges(4,3) - t374) * t271 + (-t1 * t283 - t10 * t181 - t18 * t180 - t2 * t286) * mrSges(6,3) + (-t269 * t88 + t341) * mrSges(4,3) + (-t320 * t348 - t322 * t344 - t324 * t345) * t267 + t379 * t269 / 0.2e1 + t304 * t370 + t200 * t371 + t353 * t257 + t47 * (mrSges(6,1) * t180 + mrSges(6,2) * t181) + t77 * (mrSges(5,1) * t180 + mrSges(5,2) * t181) + t355 * pkin(7) * t166 + (-t180 * t26 - t181 * t25 - t283 * t4 - t286 * t3) * mrSges(5,3) + t283 * t381 - t363 * t247 / 0.2e1 + (t195 * t387 + t352 * t320 + t350 * t322 + t351 * t324 + t185) * qJD(3) - t140 * t258 + (Ifges(4,1) * t134 + Ifges(4,4) * t371 + t5 * t205 + t326 * t347 + t327 * t346 + t328 * t343) * t166 - t373 * t286 / 0.2e1 + (-t290 * t329 + t369 * (-t106 * t282 + t290) - t367 * (t106 * t284 + t107 * t168) + t339 * (-t106 * pkin(2) + pkin(7) * t107) + t368 * t107 + t332 * t106) * g(2) + (-t289 * t329 + t369 * (-t108 * t282 + t289) - t367 * (t108 * t284 + t109 * t168) + t339 * (-t108 * pkin(2) + pkin(7) * t109) + t368 * t109 + t332 * t108) * g(1) + (t3 * t97 + t4 * t96 + (t269 * t77 + t295) * pkin(7) + (t39 - t83) * t26 + (-t82 + t40) * t25) * m(5) + t184 * t219 + t207 * t295 + ((-Ifges(4,2) * t166 + t304) * t219 + pkin(7) * t102 - t140 * t252 - t1 * mrSges(6,1) + Ifges(4,4) * t370 + Ifges(4,2) * t371 + t366 * t328 - t383 * t327 + t382 * t326 + t330) * t169 + Ifges(3,3) * qJDD(2) + (t1 * t61 + t135 * t5 + t2 * t73 + t47 * t90 + (-t83 + t16) * t18 + (-t82 + t11) * t10) * m(6); t362 * t266 / 0.2e1 - t160 * t12 + (t248 / 0.2e1 - t268 / 0.2e1) * t363 + (t336 * t66 + t337 * t65) * g(2) + (t111 * t337 + t112 * t336) * g(3) + (t336 * t68 + t337 * t67) * g(1) + t141 * t31 - t142 * t33 + Ifges(4,6) * t133 + Ifges(4,5) * t134 - t45 * t87 - t46 * t85 - t22 * mrSges(4,2) + t23 * mrSges(4,1) + (-t140 + t254) * t88 + t165 * t381 + (-t85 * t268 + t168 * t34 - t165 * t32 - t87 * t266 + m(5) * ((-t165 * t26 - t168 * t25) * qJD(4) + t340)) * pkin(8) + (-t25 * t266 - t26 * t268 + t340) * mrSges(5,3) - pkin(3) * t13 + t375 * (-t47 * t205 - t77 * t207) + (-t1 * t165 - t10 * t266 + t168 * t2 - t18 * t268) * mrSges(6,3) + t5 * t206 + t21 * t208 + t374 * t274 + t344 * t328 + t345 * t327 + t348 * t326 - (-Ifges(4,2) * t274 + t161 + t379) * t272 / 0.2e1 + (t128 * t346 + t129 * t343 - t347 * t375) * qJD(4) / 0.2e1 - (t128 * t351 + t129 * t350 - t352 * t375) * qJD(2) / 0.2e1 + t354 * t63 + t356 * t86 + t357 * t84 + (t1 * t141 + t10 * t356 - t142 * t2 - t160 * t5 + t18 * t357 + t354 * t47) * m(6) + (t255 + t331) * t89 + (-t10 * (mrSges(6,1) * t166 - mrSges(6,3) * t282) - t25 * (mrSges(5,1) * t166 - mrSges(5,3) * t282) - t18 * (-mrSges(6,2) * t166 - mrSges(6,3) * t284) - t26 * (-mrSges(5,2) * t166 - mrSges(5,3) * t284) - t185) * qJD(2) + (-pkin(3) * t21 - t25 * t45 - t26 * t46) * m(5) + t373 * t168 / 0.2e1 - t195 * t264 / 0.2e1 - t171 * t184 / 0.2e1 + Ifges(4,3) * qJDD(3); -t330 + (t87 + t308) * t26 - t77 * (mrSges(5,1) * t129 + mrSges(5,2) * t128) - t47 * (mrSges(6,1) * t129 + mrSges(6,2) * t128) - t17 * t84 + t363 * t322 + (-t85 + t309) * t25 + (t86 + t306 - m(6) * (-t10 + t17)) * t18 + (-t367 * (-t106 * t165 - t168 * t66) + t338 * (t106 * t168 - t165 * t66)) * g(2) + (-t367 * (-t108 * t165 - t168 * t68) + t338 * (t108 * t168 - t165 * t68)) * g(1) + (-t128 * t366 - t129 * t383) * t375 / 0.2e1 - (t386 * t128 - t376) * t129 / 0.2e1 + (-t384 * t129 + t362 + t378) * t388 + (-t188 * t367 + t338 * t71) * g(3) + t358 * t1 + t10 * t307 + t342 + (t129 * t335 + t31) * pkin(4); -t128 * t84 + t129 * t86 + (-g(1) * t67 - g(2) * t65 - g(3) * t111 + t10 * t129 - t18 * t128 + t5) * m(6) + t12;];
tau = t6;
