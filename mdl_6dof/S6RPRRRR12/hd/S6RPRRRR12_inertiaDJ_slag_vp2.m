% Calculate time derivative of joint inertia matrix for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_inertiaDJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:35
% EndTime: 2019-03-09 07:51:02
% DurationCPUTime: 11.26s
% Computational Cost: add. (37748->850), mult. (114122->1255), div. (0->0), fcn. (124553->16), ass. (0->346)
t280 = sin(pkin(7));
t284 = cos(pkin(7));
t285 = cos(pkin(6));
t278 = sin(pkin(14));
t281 = sin(pkin(6));
t282 = cos(pkin(14));
t363 = t281 * t282;
t394 = pkin(1) * t285;
t351 = qJ(2) * t363 + t278 * t394;
t185 = (t280 * t285 + t284 * t363) * pkin(10) + t351;
t270 = t282 * t394;
t369 = t278 * t281;
t188 = pkin(2) * t285 + t270 + (-pkin(10) * t284 - qJ(2)) * t369;
t207 = (-pkin(10) * t278 * t280 - pkin(2) * t282 - pkin(1)) * t281;
t289 = sin(qJ(3));
t293 = cos(qJ(3));
t359 = t284 * t293;
t364 = t280 * t293;
t119 = -t185 * t289 + t188 * t359 + t207 * t364;
t350 = qJD(2) * t281;
t332 = t278 * t350;
t316 = t284 * t332;
t331 = t282 * t350;
t110 = t119 * qJD(3) - t316 * t289 + t293 * t331;
t360 = t284 * t289;
t365 = t280 * t289;
t187 = t281 * (t278 * t293 + t282 * t360) + t285 * t365;
t183 = t187 * qJD(3);
t231 = -t280 * t363 + t284 * t285;
t283 = cos(pkin(8));
t392 = pkin(11) * t283;
t102 = pkin(3) * t231 - t187 * t392 + t119;
t377 = t102 * t283;
t434 = qJD(4) * t377 - t183 * t392 + t110;
t286 = sin(qJ(6));
t290 = cos(qJ(6));
t287 = sin(qJ(5));
t342 = qJD(6) * t287;
t291 = cos(qJ(5));
t344 = qJD(5) * t291;
t296 = -t286 * t342 + t290 * t344;
t433 = -m(6) * pkin(4) - mrSges(6,1) * t291 + mrSges(6,2) * t287;
t288 = sin(qJ(4));
t362 = t283 * t288;
t279 = sin(pkin(8));
t292 = cos(qJ(4));
t366 = t279 * t292;
t237 = pkin(3) * t362 + pkin(11) * t366;
t219 = pkin(12) * t283 + t237;
t220 = (-pkin(4) * t292 - pkin(12) * t288 - pkin(3)) * t279;
t428 = t291 * t219 + t287 * t220;
t432 = qJD(5) * t428;
t186 = t285 * t364 + (-t278 * t289 + t282 * t359) * t281;
t431 = 0.2e1 * t186;
t430 = 0.2e1 * t187;
t397 = t286 / 0.2e1;
t395 = t290 / 0.2e1;
t372 = t183 * t279;
t429 = Ifges(5,6) * t372;
t368 = t279 * t288;
t271 = pkin(11) * t368;
t361 = t283 * t292;
t236 = pkin(3) * t361 - t271;
t355 = t292 * t293;
t358 = t288 * t289;
t427 = t283 * t355 - t358;
t426 = -mrSges(5,1) + t433;
t150 = -t186 * t279 + t231 * t283;
t144 = -t188 * t280 + t284 * t207;
t393 = pkin(11) * t279;
t114 = -pkin(3) * t186 - t187 * t393 + t144;
t178 = t293 * t185;
t120 = t188 * t360 + t207 * t365 + t178;
t305 = t186 * t283 + t231 * t279;
t99 = pkin(11) * t305 + t120;
t49 = t102 * t362 + t114 * t368 + t292 * t99;
t44 = pkin(12) * t150 + t49;
t371 = t187 * t288;
t123 = -t186 * t361 - t231 * t366 + t371;
t124 = t187 * t292 + t288 * t305;
t67 = -t102 * t279 + t283 * t114;
t47 = pkin(4) * t123 - pkin(12) * t124 + t67;
t389 = t287 * t47 + t291 * t44;
t425 = qJD(5) * t389;
t48 = -t288 * t99 + t292 * (t114 * t279 + t377);
t182 = t186 * qJD(3);
t317 = t280 * t332;
t131 = pkin(3) * t183 - t182 * t393 + t317;
t347 = qJD(4) * t288;
t328 = t279 * t347;
t346 = qJD(4) * t292;
t111 = (-t278 * t359 - t282 * t289) * t350 + (-t178 + (-t188 * t284 - t207 * t280) * t289) * qJD(3);
t92 = -t182 * t392 + t111;
t27 = (t131 * t279 + t283 * t92) * t292 - t114 * t328 - t99 * t346 - t434 * t288;
t348 = qJD(4) * t279;
t228 = (pkin(4) * t288 - pkin(12) * t292) * t348;
t229 = t236 * qJD(4);
t130 = t228 * t291 - t229 * t287 - t432;
t327 = t279 * t346;
t26 = t114 * t327 + t131 * t368 + t434 * t292 - t347 * t99 + t92 * t362;
t24 = pkin(12) * t372 + t26;
t70 = t283 * t131 - t279 * t92;
t86 = qJD(4) * t124 + t182 * t288 + t183 * t361;
t87 = -t183 * t362 + t182 * t292 + (t292 * t305 - t371) * qJD(4);
t38 = pkin(4) * t86 - pkin(12) * t87 + t70;
t6 = -t24 * t287 + t291 * t38 - t425;
t256 = -mrSges(7,1) * t290 + mrSges(7,2) * t286;
t424 = -m(7) * pkin(5) - mrSges(6,1) + t256;
t423 = 2 * m(5);
t422 = 0.2e1 * m(6);
t421 = 0.2e1 * m(7);
t420 = 0.2e1 * pkin(12);
t419 = -2 * mrSges(5,3);
t101 = t124 * t291 + t150 * t287;
t367 = t279 * t291;
t61 = qJD(5) * t101 - t183 * t367 + t287 * t87;
t100 = t124 * t287 - t150 * t291;
t62 = -qJD(5) * t100 + t287 * t372 + t291 * t87;
t23 = Ifges(6,1) * t62 - Ifges(6,4) * t61 + Ifges(6,5) * t86;
t417 = t23 / 0.2e1;
t54 = Ifges(6,1) * t101 - Ifges(6,4) * t100 + Ifges(6,5) * t123;
t416 = t54 / 0.2e1;
t233 = -t291 * t283 + t287 * t368;
t197 = -qJD(5) * t233 + t291 * t327;
t234 = t283 * t287 + t288 * t367;
t198 = qJD(5) * t234 + t287 * t327;
t143 = Ifges(6,1) * t197 - Ifges(6,4) * t198 + Ifges(6,5) * t328;
t415 = t143 / 0.2e1;
t166 = Ifges(6,1) * t234 - Ifges(6,4) * t233 - Ifges(6,5) * t366;
t414 = t166 / 0.2e1;
t382 = Ifges(7,4) * t286;
t260 = Ifges(7,2) * t290 + t382;
t381 = Ifges(7,4) * t290;
t311 = -Ifges(7,2) * t286 + t381;
t172 = -t260 * t342 + (Ifges(7,6) * t287 + t291 * t311) * qJD(5);
t413 = t172 / 0.2e1;
t262 = Ifges(7,1) * t286 + t381;
t312 = Ifges(7,1) * t290 - t382;
t173 = -t262 * t342 + (Ifges(7,5) * t287 + t291 * t312) * qJD(5);
t412 = t173 / 0.2e1;
t199 = -t234 * t286 - t290 * t366;
t411 = t199 / 0.2e1;
t301 = -t234 * t290 + t286 * t366;
t410 = -t301 / 0.2e1;
t222 = -Ifges(7,6) * t291 + t287 * t311;
t409 = t222 / 0.2e1;
t223 = -Ifges(7,5) * t291 + t287 * t312;
t408 = t223 / 0.2e1;
t341 = qJD(6) * t290;
t276 = Ifges(7,5) * t341;
t343 = qJD(6) * t286;
t244 = -Ifges(7,6) * t343 + t276;
t407 = t244 / 0.2e1;
t246 = t311 * qJD(6);
t406 = t246 / 0.2e1;
t248 = t312 * qJD(6);
t405 = t248 / 0.2e1;
t384 = Ifges(6,4) * t287;
t249 = (Ifges(6,1) * t291 - t384) * qJD(5);
t404 = t249 / 0.2e1;
t403 = Ifges(7,5) * t397 + Ifges(7,6) * t395;
t402 = Ifges(6,5) * t287 / 0.2e1 + Ifges(6,6) * t291 / 0.2e1;
t401 = t260 / 0.2e1;
t400 = t262 / 0.2e1;
t383 = Ifges(6,4) * t291;
t263 = Ifges(6,1) * t287 + t383;
t399 = t263 / 0.2e1;
t398 = -t286 / 0.2e1;
t396 = -t290 / 0.2e1;
t391 = pkin(12) * t291;
t71 = -t101 * t286 + t123 * t290;
t72 = t101 * t290 + t123 * t286;
t41 = -mrSges(7,1) * t71 + mrSges(7,2) * t72;
t74 = mrSges(6,1) * t123 - mrSges(6,3) * t101;
t390 = t41 - t74;
t388 = m(7) * qJD(5);
t387 = mrSges(7,3) * t287;
t386 = Ifges(5,4) * t288;
t385 = Ifges(5,4) * t292;
t380 = Ifges(7,6) * t286;
t104 = mrSges(5,1) * t150 - mrSges(5,3) * t124;
t64 = mrSges(6,1) * t100 + mrSges(6,2) * t101;
t378 = t64 - t104;
t148 = t284 * t327 + (t427 * qJD(4) + (-t283 * t358 + t355) * qJD(3)) * t280;
t356 = t289 * t292;
t357 = t288 * t293;
t299 = t283 * t357 + t356;
t191 = t280 * t299 + t284 * t368;
t232 = -t279 * t364 + t283 * t284;
t151 = t191 * t287 - t291 * t232;
t349 = qJD(3) * t280;
t330 = t289 * t349;
t315 = t279 * t330;
t112 = -qJD(5) * t151 + t148 * t291 + t287 * t315;
t376 = t112 * t291;
t152 = t191 * t291 + t232 * t287;
t113 = qJD(5) * t152 + t148 * t287 - t291 * t315;
t375 = t113 * t151;
t374 = t113 * t287;
t149 = t284 * t328 + (t299 * qJD(4) + (t283 * t356 + t357) * qJD(3)) * t280;
t190 = -t280 * t427 - t284 * t366;
t373 = t149 * t190;
t230 = t237 * qJD(4);
t370 = t190 * t230;
t146 = -mrSges(7,1) * t199 - mrSges(7,2) * t301;
t203 = -mrSges(6,1) * t366 - mrSges(6,3) * t234;
t354 = t146 - t203;
t353 = Ifges(4,5) * t182 - Ifges(4,6) * t183;
t181 = mrSges(6,1) * t233 + mrSges(6,2) * t234;
t240 = mrSges(5,1) * t283 - mrSges(5,3) * t368;
t352 = t181 - t240;
t345 = qJD(5) * t287;
t340 = qJD(6) * t291;
t31 = -qJD(6) * t72 - t286 * t62 + t290 * t86;
t32 = qJD(6) * t71 + t286 * t86 + t290 * t62;
t7 = Ifges(7,5) * t32 + Ifges(7,6) * t31 + Ifges(7,3) * t61;
t21 = Ifges(6,5) * t62 - Ifges(6,6) * t61 + Ifges(6,3) * t86;
t22 = Ifges(6,4) * t62 - Ifges(6,2) * t61 + Ifges(6,6) * t86;
t338 = t7 / 0.2e1 - t22 / 0.2e1;
t34 = Ifges(7,5) * t72 + Ifges(7,6) * t71 + Ifges(7,3) * t100;
t53 = Ifges(6,4) * t101 - Ifges(6,2) * t100 + Ifges(6,6) * t123;
t337 = t34 / 0.2e1 - t53 / 0.2e1;
t13 = -mrSges(7,1) * t31 + mrSges(7,2) * t32;
t4 = -pkin(5) * t86 - t6;
t336 = -m(7) * t4 - t13;
t142 = Ifges(6,4) * t197 - Ifges(6,2) * t198 + Ifges(6,6) * t328;
t139 = qJD(6) * t199 + t197 * t290 + t286 * t328;
t140 = qJD(6) * t301 - t197 * t286 + t290 * t328;
t80 = Ifges(7,5) * t139 + Ifges(7,6) * t140 + Ifges(7,3) * t198;
t333 = -t142 / 0.2e1 + t80 / 0.2e1;
t141 = Ifges(6,5) * t197 - Ifges(6,6) * t198 + Ifges(6,3) * t328;
t132 = -Ifges(7,5) * t301 + Ifges(7,6) * t199 + Ifges(7,3) * t233;
t165 = Ifges(6,4) * t234 - Ifges(6,2) * t233 - Ifges(6,6) * t366;
t324 = -t165 / 0.2e1 + t132 / 0.2e1;
t295 = t286 * t344 + t287 * t341;
t171 = t296 * Ifges(7,5) - Ifges(7,6) * t295 + Ifges(7,3) * t345;
t247 = (-Ifges(6,2) * t287 + t383) * qJD(5);
t323 = -t247 / 0.2e1 + t171 / 0.2e1;
t221 = -Ifges(7,3) * t291 + (Ifges(7,5) * t290 - t380) * t287;
t261 = Ifges(6,2) * t291 + t384;
t322 = -t261 / 0.2e1 + t221 / 0.2e1;
t105 = -mrSges(7,1) * t140 + mrSges(7,2) * t139;
t122 = -pkin(5) * t328 - t130;
t321 = -m(7) * t122 - t105;
t252 = (pkin(5) * t287 - pkin(13) * t291) * qJD(5);
t255 = -pkin(5) * t291 - pkin(13) * t287 - pkin(4);
t160 = t255 * t341 + t252 * t286 + (-t286 * t340 - t290 * t345) * pkin(12);
t213 = t255 * t290 - t286 * t391;
t319 = -qJD(6) * t213 + t160;
t161 = -t255 * t343 + t252 * t290 + (t286 * t345 - t290 * t340) * pkin(12);
t214 = t255 * t286 + t290 * t391;
t318 = -qJD(6) * t214 - t161;
t17 = pkin(13) * t123 + t389;
t43 = -pkin(4) * t150 - t48;
t28 = pkin(5) * t100 - pkin(13) * t101 + t43;
t10 = -t17 * t286 + t28 * t290;
t25 = -pkin(4) * t372 - t27;
t12 = pkin(5) * t61 - pkin(13) * t62 + t25;
t5 = t291 * t24 + t287 * t38 + t47 * t344 - t345 * t44;
t3 = pkin(13) * t86 + t5;
t1 = qJD(6) * t10 + t12 * t286 + t290 * t3;
t11 = t17 * t290 + t28 * t286;
t2 = -qJD(6) * t11 + t12 * t290 - t286 * t3;
t314 = t1 * t290 - t2 * t286;
t313 = mrSges(7,1) * t286 + mrSges(7,2) * t290;
t218 = t271 + (-pkin(3) * t292 - pkin(4)) * t283;
t153 = pkin(5) * t233 - pkin(13) * t234 + t218;
t155 = -pkin(13) * t366 + t428;
t115 = t153 * t290 - t155 * t286;
t129 = -t219 * t345 + t220 * t344 + t287 * t228 + t291 * t229;
t121 = pkin(13) * t328 + t129;
t135 = pkin(5) * t198 - pkin(13) * t197 + t230;
t65 = qJD(6) * t115 + t121 * t290 + t135 * t286;
t116 = t153 * t286 + t155 * t290;
t66 = -qJD(6) * t116 - t121 * t286 + t135 * t290;
t309 = -t66 * t286 + t65 * t290;
t18 = -t287 * t44 + t291 * t47;
t127 = t152 * t290 + t190 * t286;
t126 = -t152 * t286 + t190 * t290;
t162 = -t219 * t287 + t220 * t291;
t55 = Ifges(5,5) * t87 - Ifges(5,6) * t86 + Ifges(5,3) * t372;
t35 = Ifges(7,4) * t72 + Ifges(7,2) * t71 + Ifges(7,6) * t100;
t36 = Ifges(7,1) * t72 + Ifges(7,4) * t71 + Ifges(7,5) * t100;
t302 = t35 * t398 + t36 * t395;
t133 = -Ifges(7,4) * t301 + Ifges(7,2) * t199 + Ifges(7,6) * t233;
t134 = -Ifges(7,1) * t301 + Ifges(7,4) * t199 + Ifges(7,5) * t233;
t298 = t133 * t398 + t134 * t395;
t297 = t151 * t344 + t374;
t277 = Ifges(6,5) * t344;
t266 = Ifges(5,5) * t327;
t251 = -mrSges(7,1) * t291 - t290 * t387;
t250 = mrSges(7,2) * t291 - t286 * t387;
t245 = -Ifges(6,6) * t345 + t277;
t243 = (mrSges(6,1) * t287 + mrSges(6,2) * t291) * qJD(5);
t242 = t313 * qJD(6);
t241 = -mrSges(5,2) * t283 + mrSges(5,3) * t366;
t238 = t313 * t287;
t235 = (-mrSges(5,1) * t292 + mrSges(5,2) * t288) * t279;
t227 = (Ifges(5,1) * t292 - t386) * t348;
t226 = (-Ifges(5,2) * t288 + t385) * t348;
t225 = -Ifges(5,6) * t328 + t266;
t224 = (mrSges(5,1) * t288 + mrSges(5,2) * t292) * t348;
t216 = Ifges(5,5) * t283 + (Ifges(5,1) * t288 + t385) * t279;
t215 = Ifges(5,6) * t283 + (Ifges(5,2) * t292 + t386) * t279;
t206 = -mrSges(7,2) * t345 - mrSges(7,3) * t295;
t205 = mrSges(7,1) * t345 - mrSges(7,3) * t296;
t202 = mrSges(6,2) * t366 - mrSges(6,3) * t233;
t184 = mrSges(7,1) * t295 + mrSges(7,2) * t296;
t169 = -mrSges(6,2) * t328 - mrSges(6,3) * t198;
t168 = mrSges(6,1) * t328 - mrSges(6,3) * t197;
t164 = Ifges(6,5) * t234 - Ifges(6,6) * t233 - Ifges(6,3) * t366;
t159 = mrSges(7,1) * t233 + mrSges(7,3) * t301;
t158 = -mrSges(7,2) * t233 + mrSges(7,3) * t199;
t157 = mrSges(4,1) * t231 - mrSges(4,3) * t187;
t156 = -mrSges(4,2) * t231 + mrSges(4,3) * t186;
t154 = pkin(5) * t366 - t162;
t145 = mrSges(6,1) * t198 + mrSges(6,2) * t197;
t138 = mrSges(4,1) * t183 + mrSges(4,2) * t182;
t118 = -mrSges(7,2) * t198 + mrSges(7,3) * t140;
t117 = mrSges(7,1) * t198 - mrSges(7,3) * t139;
t103 = -mrSges(5,2) * t150 - mrSges(5,3) * t123;
t82 = Ifges(7,1) * t139 + Ifges(7,4) * t140 + Ifges(7,5) * t198;
t81 = Ifges(7,4) * t139 + Ifges(7,2) * t140 + Ifges(7,6) * t198;
t79 = mrSges(5,1) * t123 + mrSges(5,2) * t124;
t78 = mrSges(5,1) * t372 - mrSges(5,3) * t87;
t77 = -mrSges(5,2) * t372 - mrSges(5,3) * t86;
t76 = Ifges(5,1) * t124 - Ifges(5,4) * t123 + Ifges(5,5) * t150;
t75 = Ifges(5,4) * t124 - Ifges(5,2) * t123 + Ifges(5,6) * t150;
t73 = -mrSges(6,2) * t123 - mrSges(6,3) * t100;
t69 = -qJD(6) * t127 - t112 * t286 + t149 * t290;
t68 = qJD(6) * t126 + t112 * t290 + t149 * t286;
t63 = mrSges(5,1) * t86 + mrSges(5,2) * t87;
t57 = Ifges(5,1) * t87 - Ifges(5,4) * t86 + Ifges(5,5) * t372;
t56 = Ifges(5,4) * t87 - Ifges(5,2) * t86 + t429;
t52 = Ifges(6,5) * t101 - Ifges(6,6) * t100 + Ifges(6,3) * t123;
t51 = mrSges(7,1) * t100 - mrSges(7,3) * t72;
t50 = -mrSges(7,2) * t100 + mrSges(7,3) * t71;
t40 = mrSges(6,1) * t86 - mrSges(6,3) * t62;
t39 = -mrSges(6,2) * t86 - mrSges(6,3) * t61;
t33 = mrSges(6,1) * t61 + mrSges(6,2) * t62;
t16 = -pkin(5) * t123 - t18;
t15 = mrSges(7,1) * t61 - mrSges(7,3) * t32;
t14 = -mrSges(7,2) * t61 + mrSges(7,3) * t31;
t9 = Ifges(7,1) * t32 + Ifges(7,4) * t31 + Ifges(7,5) * t61;
t8 = Ifges(7,4) * t32 + Ifges(7,2) * t31 + Ifges(7,6) * t61;
t19 = [-0.2e1 * (mrSges(3,1) * t285 - mrSges(3,3) * t369) * t332 + 0.2e1 * (-mrSges(3,2) * t285 + mrSges(3,3) * t363) * t331 + t231 * t353 + (t1 * t11 + t10 * t2 + t16 * t4) * t421 + (t26 * t49 + t27 * t48 + t67 * t70) * t423 + 0.2e1 * (-mrSges(4,1) * t186 + mrSges(4,2) * t187) * t317 + 0.2e1 * m(4) * (t110 * t120 + t111 * t119 + t144 * t317) + (Ifges(5,5) * t124 + Ifges(5,3) * t150) * t372 + (t7 - t22) * t100 + 0.2e1 * t110 * t156 + 0.2e1 * t111 * t157 + t150 * t55 + 0.2e1 * t144 * t138 + t124 * t57 + 0.2e1 * t26 * t103 + 0.2e1 * t27 * t104 + t101 * t23 + t87 * t76 + 0.2e1 * t49 * t77 + 0.2e1 * t48 * t78 + 0.2e1 * t70 * t79 + t71 * t8 + t72 * t9 + 0.2e1 * t5 * t73 + 0.2e1 * t6 * t74 + 0.2e1 * t25 * t64 + 0.2e1 * t67 * t63 + t62 * t54 + 0.2e1 * t1 * t50 + 0.2e1 * t2 * t51 + 0.2e1 * t43 * t33 + 0.2e1 * t18 * t40 + 0.2e1 * t4 * t41 + t32 * t36 + t31 * t35 + (t21 - t56 - t429) * t123 + 0.2e1 * t16 * t13 + 0.2e1 * t10 * t15 + 0.2e1 * t11 * t14 + 0.2e1 * m(3) * (t351 * t282 + (qJ(2) * t369 - t270) * t278) * t350 + (t52 - t75) * t86 + (-t53 + t34) * t61 - (0.2e1 * mrSges(4,3) * t120 + Ifges(4,4) * t430 + Ifges(4,2) * t431 + Ifges(4,6) * t231) * t183 + (-0.2e1 * mrSges(4,3) * t119 + Ifges(4,1) * t430 + Ifges(4,4) * t431 + Ifges(4,5) * t231) * t182 + (t18 * t6 + t25 * t43 + t389 * t5) * t422 + 0.2e1 * t389 * t39; t148 * t103 + t112 * t73 + t126 * t15 + t127 * t14 + t284 * t138 + t152 * t39 + t191 * t77 + t232 * t63 + t68 * t50 + t69 * t51 + (t33 - t78) * t190 + (t13 - t40) * t151 + t378 * t149 + t390 * t113 + m(5) * (t148 * t49 - t149 * t48 - t190 * t27 + t191 * t26 + t232 * t70) + m(7) * (t1 * t127 + t10 * t69 + t11 * t68 + t113 * t16 + t126 * t2 + t151 * t4) + m(6) * (t112 * t389 - t113 * t18 + t149 * t43 - t151 * t6 + t152 * t5 + t190 * t25) + (m(4) * (t110 * t289 + t111 * t293 + t316) + (-t182 * t293 - t183 * t289) * mrSges(4,3) + ((m(4) * t120 + t156) * t293 + (-m(4) * t119 - t157 + (m(5) * t67 + t79) * t279) * t289) * qJD(3)) * t280; 0.2e1 * m(7) * (t126 * t69 + t127 * t68 + t375) + 0.2e1 * m(6) * (t112 * t152 + t373 + t375) + 0.2e1 * m(5) * (t148 * t191 + t232 * t315 + t373); t353 + (-t226 / 0.2e1 + t141 / 0.2e1) * t123 + t378 * t230 + t9 * t410 + t8 * t411 + t62 * t414 + t101 * t415 + t197 * t416 + t234 * t417 + t338 * t233 + t337 * t198 + t333 * t100 + t324 * t61 + t283 * t55 / 0.2e1 + t236 * t78 + t237 * t77 + t27 * t240 + t26 * t241 + t70 * t235 + t67 * t224 + t150 * t225 / 0.2e1 + t124 * t227 / 0.2e1 + t229 * t103 + t87 * t216 / 0.2e1 + t218 * t33 + t5 * t202 + t6 * t203 + t25 * t181 + t18 * t168 + t154 * t13 + t1 * t158 + t2 * t159 + t162 * t40 + t43 * t145 + t4 * t146 + t31 * t133 / 0.2e1 + t32 * t134 / 0.2e1 + t139 * t36 / 0.2e1 + t140 * t35 / 0.2e1 + t129 * t73 + t130 * t74 + t122 * t41 + t115 * t15 + t116 * t14 + t10 * t117 + t11 * t118 - t110 * mrSges(4,2) + t111 * mrSges(4,1) + t16 * t105 + t71 * t81 / 0.2e1 + t72 * t82 / 0.2e1 + t65 * t50 + t66 * t51 + t428 * t39 + m(6) * (t129 * t389 + t130 * t18 + t162 * t6 + t218 * t25 + t230 * t43 + t428 * t5) + (t288 * t57 / 0.2e1 + (-t21 / 0.2e1 + t56 / 0.2e1) * t292 - (-Ifges(5,3) * t283 / 0.2e1 - (Ifges(5,5) * t288 + Ifges(5,6) * t292) * t279 / 0.2e1) * t183 + (-m(5) * t70 - t63) * pkin(3) + ((-t48 * mrSges(5,3) + t76 / 0.2e1) * t292 + (-t49 * mrSges(5,3) + t52 / 0.2e1 - t75 / 0.2e1) * t288) * qJD(4)) * t279 + m(5) * (t229 * t49 - t230 * t48 + t236 * t27 + t237 * t26) + t389 * t169 + (-t215 / 0.2e1 + t164 / 0.2e1) * t86 + m(7) * (t1 * t116 + t10 * t66 + t11 * t65 + t115 * t2 + t122 * t16 + t154 * t4); t112 * t202 + t126 * t117 + t127 * t118 + t190 * t145 + t148 * t241 + t152 * t169 + t68 * t158 + t69 * t159 + t232 * t224 + (t105 - t168) * t151 + t352 * t149 + t354 * t113 + (t190 * t292 - t191 * t288) * mrSges(5,3) * t348 + (-mrSges(4,2) * t293 + (t235 * t279 - mrSges(4,1)) * t289) * t349 + m(7) * (t113 * t154 + t115 * t69 + t116 * t68 + t122 * t151 + t126 * t66 + t127 * t65) + m(6) * (t112 * t428 - t113 * t162 + t129 * t152 - t130 * t151 + t149 * t218 + t370) + m(5) * (-pkin(3) * t279 ^ 2 * t330 + t148 * t237 - t149 * t236 + t191 * t229 + t370); (-0.2e1 * pkin(3) * t224 + t288 * t227 + (-t141 + t226) * t292 + ((t236 * t419 + t216) * t292 + (t237 * t419 + t164 - t215) * t288) * qJD(4)) * t279 + (t115 * t66 + t116 * t65 + t122 * t154) * t421 + (t229 * t237 - t230 * t236) * t423 + t283 * t225 + 0.2e1 * t229 * t241 + t234 * t143 + 0.2e1 * t218 * t145 + t197 * t166 + t199 * t81 + 0.2e1 * t129 * t202 + 0.2e1 * t130 * t203 + 0.2e1 * t162 * t168 + 0.2e1 * t154 * t105 + 0.2e1 * t65 * t158 + 0.2e1 * t66 * t159 + 0.2e1 * t122 * t146 + t139 * t134 + t140 * t133 + 0.2e1 * t115 * t117 + 0.2e1 * t116 * t118 + (t129 * t428 + t130 * t162 + t218 * t230) * t422 + 0.2e1 * t428 * t169 - t301 * t82 + (t132 - t165) * t198 + (t80 - t142) * t233 + 0.2e1 * t352 * t230; t433 * t25 + (t417 + t9 * t395 + t8 * t398 - t6 * mrSges(6,3) + (t35 * t396 + t36 * t398) * qJD(6) + (-mrSges(6,3) * t389 + t337) * qJD(5) + (-qJD(5) * t73 - t40 + m(6) * (-t6 - t425) - t336) * pkin(12)) * t287 + (t5 * mrSges(6,3) + (-t18 * mrSges(6,3) + t302 + t416) * qJD(5) + (t39 + t390 * qJD(5) + m(6) * (-qJD(5) * t18 + t5) + t16 * t388) * pkin(12) - t338) * t291 + m(7) * (t1 * t214 + t10 * t161 + t11 * t160 + t2 * t213) + t32 * t408 + t31 * t409 + t72 * t412 + t71 * t413 + t62 * t399 + t86 * t402 + t101 * t404 + t322 * t61 + t323 * t100 + t2 * t251 + t4 * t238 + t43 * t243 + t123 * t245 / 0.2e1 + t1 * t250 + t10 * t205 + t11 * t206 + t213 * t15 + t214 * t14 + t16 * t184 + t160 * t50 + t161 * t51 - pkin(4) * t33 - t26 * mrSges(5,2) + t27 * mrSges(5,1) + t55; -t148 * mrSges(5,2) + t113 * t238 + t126 * t205 + t127 * t206 + t151 * t184 + t190 * t243 + t68 * t250 + t69 * t251 + m(7) * (t126 * t161 + t127 * t160 + t213 * t69 + t214 * t68) + (m(7) * t297 / 0.2e1 + m(6) * (-t152 * t345 + t297 + t376) / 0.2e1) * t420 + (t376 + t374 + (t151 * t291 - t152 * t287) * qJD(5)) * mrSges(6,3) + t426 * t149; (t415 + t82 * t395 + t81 * t398 - t130 * mrSges(6,3) + (t133 * t396 + t134 * t398) * qJD(6) + (-mrSges(6,3) * t428 + t324) * qJD(5) + (-qJD(5) * t202 - t168 + m(6) * (-t130 - t432) - t321) * pkin(12)) * t287 + t426 * t230 + m(7) * (t115 * t161 + t116 * t160 + t213 * t66 + t214 * t65) + (t129 * mrSges(6,3) + (-t162 * mrSges(6,3) + t298 + t414) * qJD(5) + (t169 + t354 * qJD(5) + m(6) * (-qJD(5) * t162 + t129) + t154 * t388) * pkin(12) - t333) * t291 + (-t292 * t245 / 0.2e1 + (-Ifges(5,6) + t402) * t347) * t279 + t140 * t409 + t173 * t410 + t172 * t411 + t197 * t399 + t234 * t404 + t139 * t408 + t322 * t198 + t323 * t233 + t66 * t251 + t122 * t238 + t218 * t243 + t65 * t250 - t229 * mrSges(5,2) + t115 * t205 + t116 * t206 + t213 * t117 + t214 * t118 + t154 * t184 + t160 * t158 + t161 * t159 - pkin(4) * t145 + t266; (t160 * t214 + t161 * t213) * t421 + 0.2e1 * t160 * t250 + 0.2e1 * t214 * t206 + 0.2e1 * t161 * t251 + 0.2e1 * t213 * t205 - 0.2e1 * pkin(4) * t243 + (-t171 + t247 + (-t222 * t286 + t223 * t290 + t238 * t420 + t263) * qJD(5)) * t291 + (t184 * t420 - t286 * t172 + t290 * t173 + t249 + (-t222 * t290 - t223 * t286) * qJD(6) + (pkin(12) ^ 2 * t291 * t421 + t221 - t261) * qJD(5)) * t287; t8 * t395 + t9 * t397 + t4 * t256 + t61 * t403 + t31 * t401 + t32 * t400 + t16 * t242 + t100 * t407 + t71 * t406 + t72 * t405 + t6 * mrSges(6,1) - t5 * mrSges(6,2) + t302 * qJD(6) + t336 * pkin(5) + ((-t10 * t290 - t11 * t286) * qJD(6) + t314) * mrSges(7,3) + (m(7) * (-t10 * t341 - t11 * t343 + t314) + t290 * t14 - t286 * t15 - t50 * t343 - t51 * t341) * pkin(13) + t21; -t112 * mrSges(6,2) + t151 * t242 + (m(7) * pkin(13) + mrSges(7,3)) * (-t69 * t286 + t68 * t290 + (-t126 * t290 - t127 * t286) * qJD(6)) + t424 * t113; t81 * t395 + t82 * t397 + t122 * t256 + t198 * t403 + t140 * t401 + t139 * t400 + t154 * t242 + t233 * t407 + t199 * t406 - t301 * t405 - t129 * mrSges(6,2) + t130 * mrSges(6,1) + t298 * qJD(6) + t321 * pkin(5) + ((-t115 * t290 - t116 * t286) * qJD(6) + t309) * mrSges(7,3) + (m(7) * (-t115 * t341 - t116 * t343 + t309) + t290 * t118 - t286 * t117 - t158 * t343 - t159 * t341) * pkin(13) + t141; -pkin(5) * t184 + t277 + (-t244 / 0.2e1 + t424 * qJD(5) * pkin(12)) * t291 + (t413 + qJD(6) * t408 + t344 * t400 + t319 * mrSges(7,3) + (m(7) * t319 - qJD(6) * t251 + t206) * pkin(13)) * t290 + (t412 - qJD(6) * t222 / 0.2e1 - t260 * t344 / 0.2e1 + t318 * mrSges(7,3) + (m(7) * t318 - qJD(6) * t250 - t205) * pkin(13)) * t286 + (t248 * t395 + t246 * t398 + pkin(12) * t242 + (t260 * t396 + t262 * t398) * qJD(6) + (pkin(12) * mrSges(6,2) - Ifges(6,6) + t403) * qJD(5)) * t287; -0.2e1 * pkin(5) * t242 + t246 * t290 + t248 * t286 + (-t260 * t286 + t262 * t290) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t69 - mrSges(7,2) * t68; mrSges(7,1) * t66 - mrSges(7,2) * t65 + t80; mrSges(7,1) * t161 - mrSges(7,2) * t160 + t171; t276 + (pkin(13) * t256 - t380) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
