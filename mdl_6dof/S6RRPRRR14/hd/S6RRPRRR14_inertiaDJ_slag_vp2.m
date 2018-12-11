% Calculate time derivative of joint inertia matrix for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaDJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:10:08
% EndTime: 2018-12-10 18:10:57
% DurationCPUTime: 15.91s
% Computational Cost: add. (42948->916), mult. (133028->1354), div. (0->0), fcn. (143005->16), ass. (0->370)
t314 = sin(pkin(6));
t467 = 0.2e1 * t314;
t326 = cos(qJ(2));
t318 = cos(pkin(6));
t428 = pkin(1) * t318;
t306 = t326 * t428;
t299 = qJD(2) * t306;
t322 = sin(qJ(2));
t317 = cos(pkin(7));
t357 = -qJ(3) * t317 - pkin(10);
t350 = t357 * t322;
t313 = sin(pkin(7));
t384 = qJD(3) * t313;
t388 = t317 * t326;
t199 = t318 * t384 + t299 + (qJD(2) * t350 + qJD(3) * t388) * t314;
t305 = t322 * t428;
t383 = qJD(3) * t322;
t393 = t314 * t326;
t207 = -t314 * t317 * t383 + (t357 * t393 - t305) * qJD(2);
t408 = qJ(3) * t313;
t231 = (-t313 * t383 + (pkin(2) * t322 - t326 * t408) * qJD(2)) * t314;
t315 = cos(pkin(14));
t311 = sin(pkin(14));
t400 = t311 * t317;
t401 = t311 * t313;
t140 = t315 * t199 + t207 * t400 + t231 * t401;
t385 = qJD(2) * t314;
t389 = t317 * t322;
t249 = (-t311 * t326 - t315 * t389) * t385;
t316 = cos(pkin(8));
t312 = sin(pkin(8));
t366 = t322 * t385;
t353 = t313 * t366;
t343 = t312 * t353;
t328 = t249 * t316 + t343;
t272 = pkin(10) * t393 + t305;
t395 = t313 * t318;
t333 = t314 * t388 + t395;
t230 = qJ(3) * t333 + t272;
t238 = pkin(2) * t318 + t314 * t350 + t306;
t254 = (-pkin(2) * t326 - t322 * t408 - pkin(1)) * t314;
t392 = t315 * t317;
t396 = t313 * t315;
t159 = -t230 * t311 + t238 * t392 + t254 * t396;
t394 = t314 * t322;
t233 = t311 * t333 + t315 * t394;
t265 = -t313 * t393 + t318 * t317;
t426 = pkin(11) * t316;
t127 = pkin(3) * t265 - t233 * t426 + t159;
t407 = t127 * t316;
t466 = pkin(11) * t328 + qJD(4) * t407 + t140;
t268 = pkin(2) * t400 + qJ(3) * t396;
t229 = (t312 * t317 + t316 * t396) * pkin(11) + t268;
t321 = sin(qJ(4));
t325 = cos(qJ(4));
t303 = pkin(2) * t392;
t237 = pkin(3) * t317 + t303 + (-qJ(3) - t426) * t401;
t427 = pkin(11) * t312;
t253 = (-pkin(3) * t315 - t311 * t427 - pkin(2)) * t313;
t337 = t237 * t316 + t253 * t312;
t157 = -t321 * t229 + t337 * t325;
t319 = sin(qJ(6));
t323 = cos(qJ(6));
t320 = sin(qJ(5));
t377 = qJD(6) * t320;
t324 = cos(qJ(5));
t379 = qJD(5) * t324;
t330 = -t319 * t377 + t323 * t379;
t187 = -t237 * t312 + t316 * t253;
t391 = t316 * t321;
t398 = t312 * t321;
t235 = t317 * t398 + (t311 * t325 + t315 * t391) * t313;
t390 = t316 * t325;
t397 = t312 * t325;
t464 = t313 * (-t311 * t321 + t315 * t390) + t317 * t397;
t152 = -pkin(4) * t464 - pkin(12) * t235 + t187;
t218 = t325 * t229;
t158 = t237 * t391 + t253 * t398 + t218;
t264 = -t312 * t396 + t316 * t317;
t155 = pkin(12) * t264 + t158;
t463 = t320 * t152 + t324 * t155;
t465 = qJD(5) * t463;
t431 = t319 / 0.2e1;
t429 = t323 / 0.2e1;
t232 = t315 * t395 + (-t311 * t322 + t315 * t388) * t314;
t190 = -t232 * t312 + t265 * t316;
t160 = t315 * t230 + t238 * t400 + t254 * t401;
t338 = t232 * t316 + t265 * t312;
t124 = pkin(11) * t338 + t160;
t188 = -t238 * t313 + t317 * t254;
t151 = -pkin(3) * t232 - t233 * t427 + t188;
t61 = t325 * t124 + t127 * t391 + t151 * t398;
t53 = pkin(12) * t190 + t61;
t405 = t233 * t321;
t166 = -t232 * t390 - t265 * t397 + t405;
t167 = t233 * t325 + t321 * t338;
t81 = -t127 * t312 + t316 * t151;
t56 = pkin(4) * t166 - pkin(12) * t167 + t81;
t423 = t320 * t56 + t324 * t53;
t461 = qJD(5) * t423;
t60 = -t321 * t124 + t325 * (t151 * t312 + t407);
t139 = -t199 * t311 + t207 * t392 + t231 * t396;
t250 = (-t311 * t389 + t315 * t326) * t385;
t112 = pkin(3) * t353 - t250 * t426 + t139;
t176 = -t207 * t313 + t317 * t231;
t145 = -pkin(3) * t249 - t250 * t427 + t176;
t382 = qJD(4) * t321;
t363 = t312 * t382;
t381 = qJD(4) * t325;
t26 = t325 * (t112 * t316 + t145 * t312) - t124 * t381 - t151 * t363 - t466 * t321;
t147 = (-t311 * t391 + t315 * t325) * t384 + t157 * qJD(4);
t225 = t464 * qJD(4);
t226 = t235 * qJD(4);
t365 = t311 * t384;
t352 = t312 * t365;
t173 = pkin(4) * t226 - pkin(12) * t225 + t352;
t58 = -t147 * t320 + t173 * t324 - t465;
t206 = -t249 * t312 + t316 * t353;
t362 = t312 * t381;
t25 = t112 * t391 - t124 * t382 + t145 * t398 + t151 * t362 + t325 * t466;
t21 = pkin(12) * t206 + t25;
t116 = t250 * t325 + t328 * t321 + (t325 * t338 - t405) * qJD(4);
t117 = qJD(4) * t167 - t249 * t390 + t250 * t321 - t325 * t343;
t80 = -t112 * t312 + t316 * t145;
t40 = pkin(4) * t117 - pkin(12) * t116 + t80;
t6 = -t21 * t320 + t324 * t40 - t461;
t289 = -mrSges(7,1) * t323 + mrSges(7,2) * t319;
t460 = -m(7) * pkin(5) - mrSges(6,1) + t289;
t459 = 2 * m(4);
t458 = 2 * m(5);
t457 = 2 * m(6);
t456 = 0.2e1 * m(7);
t455 = 0.2e1 * pkin(12);
t454 = -2 * mrSges(3,3);
t453 = -2 * mrSges(5,3);
t452 = m(6) / 0.2e1;
t451 = m(6) * pkin(4);
t125 = t167 * t320 - t324 * t190;
t73 = -qJD(5) * t125 + t116 * t324 + t206 * t320;
t126 = t167 * t324 + t190 * t320;
t74 = qJD(5) * t126 + t116 * t320 - t324 * t206;
t29 = Ifges(6,1) * t73 - Ifges(6,4) * t74 + Ifges(6,5) * t117;
t450 = t29 / 0.2e1;
t66 = Ifges(6,1) * t126 - Ifges(6,4) * t125 + Ifges(6,5) * t166;
t449 = t66 / 0.2e1;
t191 = t235 * t320 - t324 * t264;
t164 = -qJD(5) * t191 + t225 * t324;
t192 = t235 * t324 + t264 * t320;
t165 = qJD(5) * t192 + t225 * t320;
t102 = Ifges(6,1) * t164 - Ifges(6,4) * t165 + Ifges(6,5) * t226;
t448 = t102 / 0.2e1;
t134 = Ifges(6,1) * t192 - Ifges(6,4) * t191 - Ifges(6,5) * t464;
t447 = t134 / 0.2e1;
t418 = Ifges(7,4) * t319;
t293 = Ifges(7,2) * t323 + t418;
t417 = Ifges(7,4) * t323;
t347 = -Ifges(7,2) * t319 + t417;
t209 = -t293 * t377 + (Ifges(7,6) * t320 + t324 * t347) * qJD(5);
t446 = t209 / 0.2e1;
t295 = Ifges(7,1) * t319 + t417;
t348 = Ifges(7,1) * t323 - t418;
t210 = -t295 * t377 + (Ifges(7,5) * t320 + t324 * t348) * qJD(5);
t445 = t210 / 0.2e1;
t260 = -Ifges(7,6) * t324 + t320 * t347;
t444 = t260 / 0.2e1;
t261 = -Ifges(7,5) * t324 + t320 * t348;
t443 = t261 / 0.2e1;
t376 = qJD(6) * t323;
t309 = Ifges(7,5) * t376;
t378 = qJD(6) * t319;
t278 = -Ifges(7,6) * t378 + t309;
t442 = t278 / 0.2e1;
t310 = Ifges(6,5) * t379;
t380 = qJD(5) * t320;
t441 = -Ifges(6,6) * t380 / 0.2e1 + t310 / 0.2e1;
t280 = t347 * qJD(6);
t440 = t280 / 0.2e1;
t282 = t348 * qJD(6);
t439 = t282 / 0.2e1;
t420 = Ifges(6,4) * t320;
t283 = (Ifges(6,1) * t324 - t420) * qJD(5);
t438 = t283 / 0.2e1;
t437 = Ifges(7,5) * t431 + Ifges(7,6) * t429;
t436 = Ifges(6,5) * t320 / 0.2e1 + Ifges(6,6) * t324 / 0.2e1;
t435 = t293 / 0.2e1;
t434 = t295 / 0.2e1;
t419 = Ifges(6,4) * t324;
t296 = Ifges(6,1) * t320 + t419;
t433 = t296 / 0.2e1;
t432 = -t319 / 0.2e1;
t430 = -t323 / 0.2e1;
t425 = pkin(12) * t324;
t84 = -t126 * t319 + t166 * t323;
t85 = t126 * t323 + t166 * t319;
t48 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t87 = mrSges(6,1) * t166 - mrSges(6,3) * t126;
t424 = t48 - t87;
t422 = m(7) * qJD(5);
t421 = mrSges(7,3) * t320;
t416 = Ifges(7,6) * t319;
t415 = t206 * Ifges(5,6);
t414 = t249 * Ifges(4,6);
t413 = t250 * Ifges(4,5);
t262 = -pkin(10) * t366 + t299;
t412 = t262 * mrSges(3,2);
t263 = t272 * qJD(2);
t411 = t263 * mrSges(3,1);
t290 = -mrSges(6,1) * t324 + mrSges(6,2) * t320;
t410 = -mrSges(5,1) + t290;
t129 = mrSges(5,1) * t190 - mrSges(5,3) * t167;
t79 = mrSges(6,1) * t125 + mrSges(6,2) * t126;
t409 = -t129 + t79;
t148 = (t311 * t390 + t315 * t321) * t384 + (t321 * t337 + t218) * qJD(4);
t406 = t148 * t325;
t269 = -t324 * t316 + t320 * t398;
t242 = -qJD(5) * t269 + t324 * t362;
t404 = t242 * t324;
t270 = t316 * t320 + t324 * t398;
t243 = qJD(5) * t270 + t320 * t362;
t403 = t243 * t269;
t402 = t243 * t320;
t168 = -t192 * t319 - t323 * t464;
t169 = t192 * t323 - t319 * t464;
t107 = -mrSges(7,1) * t168 + mrSges(7,2) * t169;
t171 = -mrSges(6,1) * t464 - mrSges(6,3) * t192;
t387 = t107 - t171;
t156 = mrSges(6,1) * t191 + mrSges(6,2) * t192;
t197 = mrSges(5,1) * t264 - mrSges(5,3) * t235;
t386 = t156 - t197;
t178 = Ifges(5,5) * t225 - Ifges(5,6) * t226;
t375 = qJD(6) * t324;
t33 = qJD(6) * t84 + t117 * t319 + t323 * t73;
t34 = -qJD(6) * t85 + t117 * t323 - t319 * t73;
t7 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t74;
t28 = Ifges(6,4) * t73 - Ifges(6,2) * t74 + Ifges(6,6) * t117;
t374 = t7 / 0.2e1 - t28 / 0.2e1;
t27 = Ifges(6,5) * t73 - Ifges(6,6) * t74 + Ifges(6,3) * t117;
t98 = qJD(6) * t168 + t164 * t323 + t226 * t319;
t99 = -qJD(6) * t169 - t164 * t319 + t226 * t323;
t41 = Ifges(7,5) * t98 + Ifges(7,6) * t99 + Ifges(7,3) * t165;
t36 = Ifges(7,5) * t85 + Ifges(7,6) * t84 + Ifges(7,3) * t125;
t65 = Ifges(6,4) * t126 - Ifges(6,2) * t125 + Ifges(6,6) * t166;
t373 = t36 / 0.2e1 - t65 / 0.2e1;
t13 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t4 = -pkin(5) * t117 - t6;
t372 = -m(7) * t4 - t13;
t101 = Ifges(6,4) * t164 - Ifges(6,2) * t165 + Ifges(6,6) * t226;
t369 = -t101 / 0.2e1 + t41 / 0.2e1;
t133 = Ifges(6,4) * t192 - Ifges(6,2) * t191 - Ifges(6,6) * t464;
t92 = Ifges(7,5) * t169 + Ifges(7,6) * t168 + Ifges(7,3) * t191;
t368 = -t133 / 0.2e1 + t92 / 0.2e1;
t67 = Ifges(5,5) * t116 - Ifges(5,6) * t117 + Ifges(5,3) * t206;
t100 = Ifges(6,5) * t164 - Ifges(6,6) * t165 + Ifges(6,3) * t226;
t51 = -pkin(5) * t226 - t58;
t59 = -mrSges(7,1) * t99 + mrSges(7,2) * t98;
t367 = -m(7) * t51 - t59;
t329 = t319 * t379 + t320 * t376;
t208 = Ifges(7,5) * t330 - Ifges(7,6) * t329 + Ifges(7,3) * t380;
t281 = (-Ifges(6,2) * t320 + t419) * qJD(5);
t359 = -t281 / 0.2e1 + t208 / 0.2e1;
t259 = -Ifges(7,3) * t324 + (Ifges(7,5) * t323 - t416) * t320;
t294 = Ifges(6,2) * t324 + t420;
t358 = -t294 / 0.2e1 + t259 / 0.2e1;
t193 = -t249 * mrSges(4,1) + t250 * mrSges(4,2);
t287 = (pkin(5) * t320 - pkin(13) * t324) * qJD(5);
t288 = -pkin(5) * t324 - pkin(13) * t320 - pkin(4);
t200 = t288 * t376 + t287 * t319 + (-t319 * t375 - t323 * t380) * pkin(12);
t256 = t288 * t323 - t319 * t425;
t356 = -qJD(6) * t256 + t200;
t201 = -t288 * t378 + t287 * t323 + (t319 * t380 - t323 * t375) * pkin(12);
t257 = t288 * t319 + t323 * t425;
t355 = -qJD(6) * t257 - t201;
t19 = pkin(13) * t166 + t423;
t52 = -pkin(4) * t190 - t60;
t30 = pkin(5) * t125 - pkin(13) * t126 + t52;
t10 = -t19 * t319 + t30 * t323;
t22 = -pkin(4) * t206 - t26;
t12 = pkin(5) * t74 - pkin(13) * t73 + t22;
t5 = t324 * t21 + t320 * t40 + t56 * t379 - t380 * t53;
t3 = pkin(13) * t117 + t5;
t1 = qJD(6) * t10 + t12 * t319 + t3 * t323;
t11 = t19 * t323 + t30 * t319;
t2 = -qJD(6) * t11 + t12 * t323 - t3 * t319;
t351 = t1 * t323 - t2 * t319;
t349 = mrSges(7,1) * t319 + mrSges(7,2) * t323;
t83 = -pkin(13) * t464 + t463;
t154 = -pkin(4) * t264 - t157;
t97 = pkin(5) * t191 - pkin(13) * t192 + t154;
t46 = -t319 * t83 + t323 * t97;
t57 = t324 * t147 + t152 * t379 - t155 * t380 + t320 * t173;
t50 = pkin(13) * t226 + t57;
t78 = pkin(5) * t165 - pkin(13) * t164 + t148;
t14 = qJD(6) * t46 + t319 * t78 + t323 * t50;
t47 = t319 * t97 + t323 * t83;
t15 = -qJD(6) * t47 - t319 * t50 + t323 * t78;
t346 = t14 * t323 - t15 * t319;
t23 = -t320 * t53 + t324 * t56;
t88 = t152 * t324 - t155 * t320;
t37 = Ifges(7,4) * t85 + Ifges(7,2) * t84 + Ifges(7,6) * t125;
t38 = Ifges(7,1) * t85 + Ifges(7,4) * t84 + Ifges(7,5) * t125;
t336 = t37 * t432 + t38 * t429;
t93 = Ifges(7,4) * t169 + Ifges(7,2) * t168 + Ifges(7,6) * t191;
t94 = Ifges(7,1) * t169 + Ifges(7,4) * t168 + Ifges(7,5) * t191;
t335 = t429 * t94 + t432 * t93;
t244 = -t270 * t319 - t323 * t397;
t334 = -t270 * t323 + t319 * t397;
t332 = t269 * t379 + t402;
t298 = Ifges(3,5) * t326 * t385;
t286 = -mrSges(7,1) * t324 - t323 * t421;
t285 = mrSges(7,2) * t324 - t319 * t421;
t277 = (mrSges(6,1) * t320 + mrSges(6,2) * t324) * qJD(5);
t276 = t349 * qJD(6);
t275 = -mrSges(4,2) * t317 + mrSges(4,3) * t396;
t274 = mrSges(4,1) * t317 - mrSges(4,3) * t401;
t273 = t349 * t320;
t271 = -pkin(10) * t394 + t306;
t267 = -qJ(3) * t401 + t303;
t252 = -mrSges(7,2) * t380 - mrSges(7,3) * t329;
t251 = mrSges(7,1) * t380 - mrSges(7,3) * t330;
t228 = mrSges(7,1) * t329 + mrSges(7,2) * t330;
t212 = mrSges(4,1) * t353 - mrSges(4,3) * t250;
t211 = -mrSges(4,2) * t353 + mrSges(4,3) * t249;
t196 = mrSges(4,1) * t265 - mrSges(4,3) * t233;
t195 = -mrSges(5,2) * t264 + mrSges(5,3) * t464;
t194 = -mrSges(4,2) * t265 + mrSges(4,3) * t232;
t186 = t250 * Ifges(4,1) + t249 * Ifges(4,4) + Ifges(4,5) * t353;
t185 = t250 * Ifges(4,4) + t249 * Ifges(4,2) + Ifges(4,6) * t353;
t184 = Ifges(4,3) * t353 + t413 + t414;
t183 = -mrSges(5,1) * t464 + mrSges(5,2) * t235;
t182 = qJD(6) * t334 - t242 * t319 + t323 * t363;
t181 = qJD(6) * t244 + t242 * t323 + t319 * t363;
t180 = Ifges(5,1) * t225 - Ifges(5,4) * t226;
t179 = Ifges(5,4) * t225 - Ifges(5,2) * t226;
t177 = mrSges(5,1) * t226 + mrSges(5,2) * t225;
t175 = Ifges(5,1) * t235 + Ifges(5,4) * t464 + Ifges(5,5) * t264;
t174 = Ifges(5,4) * t235 + Ifges(5,2) * t464 + Ifges(5,6) * t264;
t170 = mrSges(6,2) * t464 - mrSges(6,3) * t191;
t142 = -mrSges(6,2) * t226 - mrSges(6,3) * t165;
t141 = mrSges(6,1) * t226 - mrSges(6,3) * t164;
t132 = Ifges(6,5) * t192 - Ifges(6,6) * t191 - Ifges(6,3) * t464;
t131 = mrSges(7,1) * t191 - mrSges(7,3) * t169;
t130 = -mrSges(7,2) * t191 + mrSges(7,3) * t168;
t128 = -mrSges(5,2) * t190 - mrSges(5,3) * t166;
t106 = mrSges(5,1) * t166 + mrSges(5,2) * t167;
t105 = mrSges(6,1) * t165 + mrSges(6,2) * t164;
t104 = -mrSges(5,2) * t206 - mrSges(5,3) * t117;
t103 = mrSges(5,1) * t206 - mrSges(5,3) * t116;
t91 = Ifges(5,1) * t167 - Ifges(5,4) * t166 + Ifges(5,5) * t190;
t90 = Ifges(5,4) * t167 - Ifges(5,2) * t166 + Ifges(5,6) * t190;
t86 = -mrSges(6,2) * t166 - mrSges(6,3) * t125;
t82 = pkin(5) * t464 - t88;
t77 = -mrSges(7,2) * t165 + mrSges(7,3) * t99;
t76 = mrSges(7,1) * t165 - mrSges(7,3) * t98;
t75 = mrSges(5,1) * t117 + mrSges(5,2) * t116;
t69 = Ifges(5,1) * t116 - Ifges(5,4) * t117 + t206 * Ifges(5,5);
t68 = Ifges(5,4) * t116 - Ifges(5,2) * t117 + t415;
t64 = Ifges(6,5) * t126 - Ifges(6,6) * t125 + Ifges(6,3) * t166;
t63 = mrSges(7,1) * t125 - mrSges(7,3) * t85;
t62 = -mrSges(7,2) * t125 + mrSges(7,3) * t84;
t45 = -mrSges(6,2) * t117 - mrSges(6,3) * t74;
t44 = mrSges(6,1) * t117 - mrSges(6,3) * t73;
t43 = Ifges(7,1) * t98 + Ifges(7,4) * t99 + Ifges(7,5) * t165;
t42 = Ifges(7,4) * t98 + Ifges(7,2) * t99 + Ifges(7,6) * t165;
t35 = mrSges(6,1) * t74 + mrSges(6,2) * t73;
t18 = -pkin(5) * t166 - t23;
t17 = -mrSges(7,2) * t74 + mrSges(7,3) * t34;
t16 = mrSges(7,1) * t74 - mrSges(7,3) * t33;
t9 = Ifges(7,1) * t33 + Ifges(7,4) * t34 + Ifges(7,5) * t74;
t8 = Ifges(7,4) * t33 + Ifges(7,2) * t34 + Ifges(7,6) * t74;
t20 = [t249 * (Ifges(4,4) * t233 + Ifges(4,2) * t232 + Ifges(4,6) * t265) + t250 * (Ifges(4,1) * t233 + Ifges(4,4) * t232 + Ifges(4,5) * t265) + t265 * t184 + t232 * t185 + 0.2e1 * t176 * (-mrSges(4,1) * t232 + mrSges(4,2) * t233) + t233 * t186 + 0.2e1 * t160 * t211 + 0.2e1 * t159 * t212 + 0.2e1 * t139 * t196 + t190 * t67 + 0.2e1 * t188 * t193 + 0.2e1 * t140 * t194 + t167 * t69 + t126 * t29 + 0.2e1 * t25 * t128 + 0.2e1 * t26 * t129 + t116 * t91 + (t27 - t68 - t415) * t166 + 0.2e1 * t60 * t103 + 0.2e1 * t61 * t104 + 0.2e1 * t80 * t106 + 0.2e1 * t6 * t87 + t84 * t8 + t85 * t9 + 0.2e1 * t5 * t86 + 0.2e1 * t81 * t75 + 0.2e1 * t22 * t79 + t73 * t66 + 0.2e1 * t2 * t63 + 0.2e1 * t1 * t62 + 0.2e1 * t52 * t35 + 0.2e1 * t4 * t48 + 0.2e1 * t23 * t44 + t33 * t38 + t34 * t37 + 0.2e1 * t18 * t13 + 0.2e1 * t11 * t17 + 0.2e1 * t10 * t16 + 0.2e1 * m(3) * (t262 * t272 - t263 * t271) + (t298 - 0.2e1 * t411 - 0.2e1 * t412) * t318 + 0.2e1 * t423 * t45 + (t22 * t52 + t23 * t6 + t423 * t5) * t457 + (-t65 + t36) * t74 + (0.2e1 * (t262 * t326 + t263 * t322) * mrSges(3,3) + ((t271 * t454 + Ifges(3,5) * t318 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t326) * t467) * t326 + (-0.2e1 * Ifges(3,6) * t318 + t313 * (Ifges(4,5) * t233 + Ifges(4,6) * t232 + Ifges(4,3) * t265) + t272 * t454 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t322 + (Ifges(3,1) - Ifges(3,2)) * t326) * t467) * t322) * qJD(2)) * t314 + t206 * (Ifges(5,5) * t167 + Ifges(5,3) * t190) + (t1 * t11 + t10 * t2 + t18 * t4) * t456 + (t25 * t61 + t26 * t60 + t80 * t81) * t458 + (t139 * t159 + t140 * t160 + t176 * t188) * t459 + (t64 - t90) * t117 + (-t28 + t7) * t125; (-t174 / 0.2e1 + t132 / 0.2e1) * t117 + t140 * t275 + t267 * t212 + t268 * t211 + t139 * t274 + t264 * t67 / 0.2e1 + t235 * t69 / 0.2e1 + t25 * t195 + t26 * t197 + t190 * t178 / 0.2e1 + t80 * t183 + t187 * t75 + t5 * t170 + t6 * t171 + t116 * t175 / 0.2e1 + t81 * t177 + t167 * t180 / 0.2e1 + t168 * t8 / 0.2e1 + t169 * t9 / 0.2e1 + t22 * t156 + t157 * t103 + t158 * t104 + t154 * t35 + t23 * t141 + t147 * t128 + t2 * t131 + t1 * t130 + t464 * t68 / 0.2e1 - t464 * t27 / 0.2e1 + t206 * (Ifges(5,5) * t235 + Ifges(5,6) * t464 + Ifges(5,3) * t264) / 0.2e1 + (t184 / 0.2e1 + t414 / 0.2e1 + t413 / 0.2e1) * t317 + t409 * t148 + t52 * t105 + t4 * t107 + t98 * t38 / 0.2e1 + t99 * t37 / 0.2e1 + t34 * t93 / 0.2e1 + t33 * t94 / 0.2e1 + t88 * t44 + t84 * t42 / 0.2e1 + t85 * t43 / 0.2e1 + t57 * t86 + t58 * t87 + t82 * t13 + t11 * t77 + t10 * t76 + t15 * t63 + t14 * t62 + t18 * t59 + t373 * t165 + t374 * t191 + t51 * t48 + t46 * t16 + t47 * t17 + t368 * t74 + t369 * t125 + (t176 * (-mrSges(4,1) * t315 + mrSges(4,2) * t311) + t249 * (Ifges(4,4) * t311 + Ifges(4,2) * t315) / 0.2e1 + t250 * (Ifges(4,1) * t311 + Ifges(4,4) * t315) / 0.2e1 + t311 * t186 / 0.2e1 - pkin(2) * t193 + t315 * t185 / 0.2e1 + (Ifges(4,3) * t317 / 0.2e1 + (Ifges(4,5) * t311 + Ifges(4,6) * t315) * t313 / 0.2e1) * t366 + (t315 * t194 + (t106 * t312 - t196) * t311) * qJD(3)) * t313 - Ifges(3,6) * t366 + m(5) * (t147 * t61 - t148 * t60 + t157 * t26 + t158 * t25 + t187 * t80 + t352 * t81) + t298 + m(4) * (t139 * t267 + t140 * t268 + (-pkin(2) * t176 + (-t159 * t311 + t160 * t315) * qJD(3)) * t313) + t423 * t142 + m(7) * (t1 * t47 + t10 * t15 + t11 * t14 + t18 * t51 + t2 * t46 + t4 * t82) + t463 * t45 + m(6) * (t148 * t52 + t154 * t22 + t23 * t58 + t423 * t57 + t463 * t5 + t6 * t88) - t411 - t412 + t73 * t447 + t126 * t448 + t164 * t449 + t192 * t450 + (-t179 / 0.2e1 + t100 / 0.2e1) * t166 + (t91 / 0.2e1 - t60 * mrSges(5,3)) * t225 + (t64 / 0.2e1 - t90 / 0.2e1 - t61 * mrSges(5,3)) * t226; t264 * t178 + t235 * t180 + (t158 * t453 + t132 - t174) * t226 + (t157 * t453 + t175) * t225 + 0.2e1 * t147 * t195 + t192 * t102 + 0.2e1 * t187 * t177 + 0.2e1 * t57 * t170 + 0.2e1 * t58 * t171 + t168 * t42 + t169 * t43 + t164 * t134 + 0.2e1 * t154 * t105 + 0.2e1 * t88 * t141 + 0.2e1 * t14 * t130 + 0.2e1 * t15 * t131 - (t100 - t179) * t464 + 0.2e1 * t51 * t107 + t98 * t94 + t99 * t93 + 0.2e1 * t82 * t59 + 0.2e1 * t46 * t76 + 0.2e1 * t47 * t77 + (t41 - t101) * t191 + (t92 - t133) * t165 + 0.2e1 * t386 * t148 + ((-t267 * t311 + t268 * t315) * t459 + 0.2e1 * t275 * t315 + 0.2e1 * (t183 * t312 - t274) * t311) * t384 + 0.2e1 * t463 * t142 + (t148 * t154 + t463 * t57 + t58 * t88) * t457 + (t14 * t47 + t15 * t46 + t51 * t82) * t456 + (t147 * t158 - t148 * t157 + t187 * t352) * t458; t244 * t16 - t334 * t17 + t181 * t62 + t182 * t63 + t242 * t86 + t270 * t45 + (t13 - t44) * t269 + t424 * t243 + m(7) * (-t1 * t334 + t10 * t182 + t11 * t181 + t18 * t243 + t2 * t244 + t269 * t4) + m(6) * (-t23 * t243 + t242 * t423 - t269 * t6 + t270 * t5) + m(4) * t176 + t193 + (m(5) * t80 + t75) * t316 + (t321 * t104 + (t103 - t35) * t325 + (t128 * t325 + t321 * t409) * qJD(4) + 0.2e1 * (-t22 * t325 + t382 * t52) * t452 + m(5) * (t25 * t321 + t26 * t325 + t381 * t61 - t382 * t60)) * t312; t181 * t130 + t182 * t131 + t270 * t142 + t242 * t170 + t316 * t177 + t244 * t76 - t334 * t77 + (-t141 + t59) * t269 + t387 * t243 + m(7) * (-t14 * t334 + t15 * t244 + t181 * t47 + t182 * t46 + t243 * t82 + t269 * t51) + m(6) * (t242 * t463 - t243 * t88 - t269 * t58 + t270 * t57) + (-t325 * t105 + (-t225 * t325 - t226 * t321) * mrSges(5,3) + (t325 * t195 + t321 * t386) * qJD(4) + m(6) * (t154 * t382 - t406) + m(5) * (t147 * t321 - t157 * t382 + t158 * t381 + t316 * t365 - t406)) * t312; 0.2e1 * m(7) * (-t181 * t334 + t182 * t244 + t403) + 0.2e1 * m(6) * (-t312 ^ 2 * t321 * t381 + t242 * t270 + t403); t67 + m(7) * (t1 * t257 + t10 * t201 + t11 * t200 + t2 * t256) + (t450 + t9 * t429 + t8 * t432 - t6 * mrSges(6,3) + (t37 * t430 + t38 * t432) * qJD(6) + (-mrSges(6,3) * t423 + t373) * qJD(5) + (-qJD(5) * t86 - t44 + m(6) * (-t6 - t461) - t372) * pkin(12)) * t320 + t52 * t277 + t1 * t285 + t2 * t286 + t4 * t273 + t256 * t16 + t257 * t17 + t10 * t251 + t11 * t252 + (t5 * mrSges(6,3) + (-t23 * mrSges(6,3) + t336 + t449) * qJD(5) + (t45 + t424 * qJD(5) + t18 * t422 + m(6) * (-qJD(5) * t23 + t5)) * pkin(12) - t374) * t324 + t18 * t228 + t200 * t62 + t201 * t63 - pkin(4) * t35 - t25 * mrSges(5,2) + t26 * mrSges(5,1) + t358 * t74 + t359 * t125 + (t290 - t451) * t22 + t73 * t433 + t117 * t436 + t126 * t438 + t166 * t441 + t33 * t443 + t34 * t444 + t85 * t445 + t84 * t446; t154 * t277 + t14 * t285 + t15 * t286 + t51 * t273 + t256 * t76 + t257 * t77 + t46 * t251 + t47 * t252 + (t410 - t451) * t148 + t82 * t228 + (t57 * mrSges(6,3) + (-t88 * mrSges(6,3) + t335 + t447) * qJD(5) + (t142 + t387 * qJD(5) + t82 * t422 + m(6) * (-qJD(5) * t88 + t57)) * pkin(12) - t369) * t324 + t200 * t130 + t201 * t131 - t147 * mrSges(5,2) - t464 * t441 - pkin(4) * t105 + t358 * t165 + t359 * t191 + m(7) * (t14 * t257 + t15 * t256 + t200 * t47 + t201 * t46) + t178 + (t448 + t43 * t429 + t42 * t432 - t58 * mrSges(6,3) + (t430 * t93 + t432 * t94) * qJD(6) + (-mrSges(6,3) * t463 + t368) * qJD(5) + (-qJD(5) * t170 - t141 + m(6) * (-t58 - t465) - t367) * pkin(12)) * t320 + t164 * t433 + t226 * t436 + t192 * t438 + t98 * t443 + t99 * t444 + t169 * t445 + t168 * t446; t181 * t285 + t182 * t286 + t269 * t228 + t243 * t273 + t244 * t251 - t334 * t252 + (-t325 * t277 + (-mrSges(5,2) * t325 + t321 * t410) * qJD(4)) * t312 + m(7) * (t181 * t257 + t182 * t256 - t200 * t334 + t201 * t244) - t363 * t451 + (m(7) * t332 / 0.2e1 + (-t270 * t380 + t332 + t404) * t452) * t455 + (t404 + t402 + (t269 * t324 - t270 * t320) * qJD(5)) * mrSges(6,3); (t200 * t257 + t201 * t256) * t456 + 0.2e1 * t200 * t285 + 0.2e1 * t257 * t252 + 0.2e1 * t201 * t286 + 0.2e1 * t256 * t251 - 0.2e1 * pkin(4) * t277 + (-t208 + t281 + (-t260 * t319 + t261 * t323 + t273 * t455 + t296) * qJD(5)) * t324 + (t228 * t455 - t319 * t209 + t323 * t210 + t283 + (-t260 * t323 - t261 * t319) * qJD(6) + (pkin(12) ^ 2 * t324 * t456 + t259 - t294) * qJD(5)) * t320; t8 * t429 + t9 * t431 + t4 * t289 + t74 * t437 + t34 * t435 + t33 * t434 + t18 * t276 + t125 * t442 + t84 * t440 + t85 * t439 + t6 * mrSges(6,1) - t5 * mrSges(6,2) + t336 * qJD(6) + t372 * pkin(5) + ((-t10 * t323 - t11 * t319) * qJD(6) + t351) * mrSges(7,3) + (m(7) * (-t10 * t376 - t11 * t378 + t351) + t323 * t17 - t319 * t16 - t62 * t378 - t63 * t376) * pkin(13) + t27; t42 * t429 + t43 * t431 + t51 * t289 + t165 * t437 + t99 * t435 + t98 * t434 + t82 * t276 + t191 * t442 + t168 * t440 + t169 * t439 + t58 * mrSges(6,1) - t57 * mrSges(6,2) + t335 * qJD(6) + t367 * pkin(5) + ((-t319 * t47 - t323 * t46) * qJD(6) + t346) * mrSges(7,3) + (m(7) * (-t376 * t46 - t378 * t47 + t346) + t323 * t77 - t319 * t76 - t130 * t378 - t131 * t376) * pkin(13) + t100; -t242 * mrSges(6,2) + t269 * t276 + (m(7) * pkin(13) + mrSges(7,3)) * (t181 * t323 - t182 * t319 + (-t244 * t323 + t319 * t334) * qJD(6)) + t460 * t243; -pkin(5) * t228 + t310 + (-t278 / 0.2e1 + t460 * qJD(5) * pkin(12)) * t324 + (t446 + qJD(6) * t443 + t379 * t434 + t356 * mrSges(7,3) + (m(7) * t356 - qJD(6) * t286 + t252) * pkin(13)) * t323 + (t445 - qJD(6) * t260 / 0.2e1 - t293 * t379 / 0.2e1 + t355 * mrSges(7,3) + (m(7) * t355 - qJD(6) * t285 - t251) * pkin(13)) * t319 + (t282 * t429 + t280 * t432 + pkin(12) * t276 + (t293 * t430 + t295 * t432) * qJD(6) + (pkin(12) * mrSges(6,2) - Ifges(6,6) + t437) * qJD(5)) * t320; -0.2e1 * pkin(5) * t276 + t280 * t323 + t282 * t319 + (-t293 * t319 + t295 * t323) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t15 - mrSges(7,2) * t14 + t41; mrSges(7,1) * t182 - mrSges(7,2) * t181; mrSges(7,1) * t201 - mrSges(7,2) * t200 + t208; t309 + (pkin(13) * t289 - t416) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
