% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:08:10
% EndTime: 2018-11-23 16:08:28
% DurationCPUTime: 17.76s
% Computational Cost: add. (29538->721), mult. (103071->1038), div. (0->0), fcn. (88925->14), ass. (0->328)
t260 = sin(pkin(12));
t262 = sin(pkin(6));
t263 = cos(pkin(12));
t268 = sin(qJ(3));
t264 = cos(pkin(7));
t271 = cos(qJ(3));
t345 = t264 * t271;
t276 = (-t260 * t345 - t263 * t268) * t262;
t228 = qJD(1) * t276;
t333 = qJD(1) * t262;
t346 = t264 * t268;
t229 = (-t260 * t346 + t263 * t271) * t333;
t259 = sin(pkin(13));
t261 = sin(pkin(7));
t353 = cos(pkin(13));
t311 = t353 * t271;
t450 = t259 * t228 + t229 * t353 - (-t259 * t268 + t311) * t261 * qJD(3);
t265 = cos(pkin(6));
t348 = t261 * t271;
t215 = t265 * t348 + (-t260 * t268 + t263 * t345) * t262;
t207 = t215 * qJD(1);
t349 = t261 * t268;
t216 = t262 * (t260 * t271 + t263 * t346) + t265 * t349;
t208 = t216 * qJD(1);
t310 = t353 * t207 - t208 * t259;
t163 = qJD(5) - t310;
t384 = pkin(3) * t259;
t257 = pkin(10) + t384;
t267 = sin(qJ(5));
t330 = qJD(5) * t267;
t280 = t259 * t207 + t208 * t353;
t385 = pkin(3) * t208;
t104 = pkin(4) * t280 - pkin(10) * t310 + t385;
t270 = cos(qJ(5));
t347 = t262 * t263;
t253 = qJ(2) * t347;
t386 = pkin(1) * t265;
t327 = qJD(1) * t386;
t238 = qJD(1) * t253 + t260 * t327;
t279 = (t261 * t265 + t264 * t347) * pkin(9);
t196 = qJD(1) * t279 + t238;
t252 = t263 * t327;
t350 = t260 * t262;
t275 = pkin(2) * t265 + (-pkin(9) * t264 - qJ(2)) * t350;
t205 = qJD(1) * t275 + t252;
t230 = (-pkin(9) * t260 * t261 - pkin(2) * t263 - pkin(1)) * t262;
t225 = qJD(1) * t230 + qJD(2);
t141 = t268 * (t205 * t264 + t225 * t261) + t271 * t196;
t127 = t207 * qJ(4) + t141;
t121 = t259 * t127;
t140 = -t196 * t268 + t205 * t345 + t225 * t348;
t126 = -qJ(4) * t208 + t140;
t69 = t126 * t353 - t121;
t41 = t267 * t104 + t270 * t69;
t452 = -t257 * t330 - t41;
t236 = (t259 * t271 + t268 * t353) * t261;
t221 = t236 * t267 - t270 * t264;
t323 = t260 * t333;
t309 = t261 * t323;
t451 = qJD(5) * t221 + t267 * t309 + t450 * t270;
t442 = -qJD(3) * t236 - t228 * t353 + t229 * t259;
t312 = t353 * t127;
t68 = t126 * t259 + t312;
t449 = -t68 + t163 * (pkin(5) * t267 - pkin(11) * t270);
t448 = pkin(11) * t280 - t452;
t239 = -t261 * t347 + t264 * t265;
t233 = qJD(1) * t239 + qJD(3);
t146 = t233 * t267 + t270 * t280;
t209 = t215 * qJD(3);
t194 = qJD(1) * t209;
t210 = t216 * qJD(3);
t195 = qJD(1) * t210;
t161 = t194 * t353 - t259 * t195;
t103 = qJD(5) * t146 + t161 * t267;
t400 = t103 / 0.2e1;
t401 = -t103 / 0.2e1;
t447 = Ifges(6,4) * t401;
t109 = mrSges(6,1) * t163 - mrSges(6,3) * t146;
t266 = sin(qJ(6));
t269 = cos(qJ(6));
t106 = -t146 * t266 + t163 * t269;
t107 = t146 * t269 + t163 * t266;
t62 = -mrSges(7,1) * t106 + mrSges(7,2) * t107;
t354 = t109 - t62;
t145 = t233 * t270 - t267 * t280;
t429 = mrSges(5,1) * t233 + mrSges(6,1) * t145 - mrSges(6,2) * t146 - mrSges(5,3) * t280;
t222 = t236 * t270 + t264 * t267;
t235 = t259 * t349 - t261 * t311;
t185 = -t222 * t266 + t235 * t269;
t446 = qJD(6) * t185 - t266 * t442 - t269 * t451;
t186 = t222 * t269 + t235 * t266;
t445 = -qJD(6) * t186 + t266 * t451 - t269 * t442;
t444 = -qJD(5) * t222 + t267 * t450 - t270 * t309;
t343 = t266 * t270;
t119 = t269 * t280 - t310 * t343;
t329 = qJD(5) * t270;
t440 = t266 * t329 + t119;
t102 = qJD(5) * t145 + t161 * t270;
t160 = t194 * t259 + t195 * t353;
t50 = qJD(6) * t106 + t102 * t269 + t160 * t266;
t407 = t50 / 0.2e1;
t51 = -qJD(6) * t107 - t102 * t266 + t160 * t269;
t406 = t51 / 0.2e1;
t113 = pkin(3) * t233 + t126;
t66 = t259 * t113 + t312;
t64 = pkin(10) * t233 + t66;
t173 = -t205 * t261 + t264 * t225;
t142 = -pkin(3) * t207 + qJD(4) + t173;
t83 = -pkin(4) * t310 - pkin(10) * t280 + t142;
t32 = t267 * t83 + t270 * t64;
t28 = pkin(11) * t163 + t32;
t65 = t113 * t353 - t121;
t63 = -t233 * pkin(4) - t65;
t37 = -t145 * pkin(5) - t146 * pkin(11) + t63;
t9 = -t266 * t28 + t269 * t37;
t439 = t9 * mrSges(7,1);
t402 = t102 / 0.2e1;
t393 = t160 / 0.2e1;
t438 = t239 / 0.2e1;
t436 = t280 / 0.2e1;
t10 = t266 * t37 + t269 * t28;
t434 = t10 * mrSges(7,2);
t22 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t80 = mrSges(6,1) * t160 - mrSges(6,3) * t102;
t374 = t22 - t80;
t433 = Ifges(5,4) * t280;
t432 = Ifges(5,4) * t310;
t431 = Ifges(5,5) * t161;
t430 = Ifges(5,6) * t160;
t324 = t353 * pkin(3);
t258 = -t324 - pkin(4);
t242 = -t270 * pkin(5) - t267 * pkin(11) + t258;
t340 = t269 * t270;
t220 = t242 * t266 + t257 * t340;
t428 = -qJD(6) * t220 + t448 * t266 + t269 * t449;
t219 = t242 * t269 - t257 * t343;
t427 = qJD(6) * t219 + t266 * t449 - t448 * t269;
t328 = qJD(6) * t267;
t426 = t269 * t328 + t440;
t120 = t266 * t280 + t310 * t340;
t425 = t266 * t328 - t269 * t329 + t120;
t424 = Ifges(4,5) * t194 - Ifges(4,6) * t195 - t430 + t431;
t273 = t141 * qJD(3);
t283 = -t194 * qJ(4) - t208 * qJD(4);
t332 = qJD(2) * t262;
t315 = qJD(1) * t332;
t308 = t260 * t315;
t284 = t264 * t308;
t307 = t263 * t315;
t331 = qJD(3) * t271;
t319 = t264 * t331;
t320 = t261 * t331;
t129 = t205 * t319 + t225 * t320 + t271 * t307 + (-qJD(3) * t196 - t284) * t268;
t90 = -qJ(4) * t195 + qJD(4) * t207 + t129;
t55 = t353 * t90 + (-t268 * t307 - t271 * t284 - t273 + t283) * t259;
t244 = t261 * t308;
t187 = pkin(3) * t195 + t244;
t89 = pkin(4) * t160 - pkin(10) * t161 + t187;
t7 = t267 * t89 + t270 * t55 + t83 * t329 - t330 * t64;
t8 = -qJD(5) * t32 - t267 * t55 + t270 * t89;
t423 = -t267 * t8 + t270 * t7;
t274 = qJD(2) * t276;
t130 = qJD(1) * t274 - t273;
t54 = t259 * t90 - t353 * (t130 + t283);
t25 = pkin(5) * t103 - pkin(11) * t102 + t54;
t5 = pkin(11) * t160 + t7;
t1 = qJD(6) * t9 + t25 * t266 + t269 * t5;
t2 = -qJD(6) * t10 + t25 * t269 - t266 * t5;
t304 = t1 * t269 - t2 * t266;
t31 = -t267 * t64 + t270 * t83;
t27 = -pkin(5) * t163 - t31;
t302 = mrSges(7,1) * t266 + mrSges(7,2) * t269;
t303 = t10 * t266 + t9 * t269;
t105 = Ifges(7,4) * t106;
t144 = qJD(6) - t145;
t45 = Ifges(7,1) * t107 + Ifges(7,5) * t144 + t105;
t355 = t269 * t45;
t367 = Ifges(7,4) * t107;
t44 = Ifges(7,2) * t106 + Ifges(7,6) * t144 + t367;
t409 = -t44 / 0.2e1;
t294 = Ifges(7,5) * t269 - Ifges(7,6) * t266;
t365 = Ifges(7,4) * t269;
t297 = -Ifges(7,2) * t266 + t365;
t366 = Ifges(7,4) * t266;
t300 = Ifges(7,1) * t269 - t366;
t394 = t144 / 0.2e1;
t396 = t107 / 0.2e1;
t398 = t106 / 0.2e1;
t417 = t294 * t394 + t297 * t398 + t300 * t396;
t422 = -t303 * mrSges(7,3) + t27 * t302 + t266 * t409 + t355 / 0.2e1 + t417;
t334 = t260 * t386 + t253;
t213 = t279 + t334;
t255 = t263 * t386;
t217 = t255 + t275;
t147 = -t213 * t268 + t217 * t345 + t230 * t348;
t128 = pkin(3) * t239 - qJ(4) * t216 + t147;
t203 = t271 * t213;
t148 = t217 * t346 + t230 * t349 + t203;
t137 = qJ(4) * t215 + t148;
t76 = t259 * t128 + t353 * t137;
t72 = pkin(10) * t239 + t76;
t176 = -t217 * t261 + t264 * t230;
t153 = -pkin(3) * t215 + t176;
t174 = -t215 * t353 + t216 * t259;
t175 = t259 * t215 + t216 * t353;
t86 = pkin(4) * t174 - pkin(10) * t175 + t153;
t373 = t267 * t86 + t270 * t72;
t321 = t263 * t332;
t322 = t260 * t332;
t138 = t217 * t319 + t230 * t320 + t271 * t321 + (-qJD(3) * t213 - t264 * t322) * t268;
t96 = -qJ(4) * t210 + qJD(4) * t215 + t138;
t139 = t274 + (-t203 + (-t217 * t264 - t230 * t261) * t268) * qJD(3);
t97 = -t209 * qJ(4) - t216 * qJD(4) + t139;
t60 = t259 * t97 + t353 * t96;
t166 = t209 * t259 + t210 * t353;
t169 = t209 * t353 - t259 * t210;
t247 = t261 * t322;
t188 = pkin(3) * t210 + t247;
t94 = pkin(4) * t166 - pkin(10) * t169 + t188;
t14 = -qJD(5) * t373 - t267 * t60 + t270 * t94;
t421 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t420 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t419 = -t434 + t439;
t360 = t144 * Ifges(7,3);
t363 = t107 * Ifges(7,5);
t364 = t106 * Ifges(7,6);
t43 = t360 + t363 + t364;
t358 = t163 * Ifges(6,6);
t370 = Ifges(6,4) * t146;
t78 = t145 * Ifges(6,2) + t358 + t370;
t418 = t78 / 0.2e1 - t43 / 0.2e1 - t419;
t416 = mrSges(6,2) * t54 + 0.2e1 * Ifges(6,1) * t402 + 0.2e1 * Ifges(6,5) * t393 + t447;
t415 = t130 * mrSges(4,1) - t54 * mrSges(5,1) - t129 * mrSges(4,2) - t55 * mrSges(5,2);
t15 = Ifges(7,5) * t50 + Ifges(7,6) * t51 + Ifges(7,3) * t103;
t414 = Ifges(7,5) * t407 + Ifges(7,6) * t406 + Ifges(7,3) * t400 + t54 * mrSges(6,1) + t15 / 0.2e1 - t102 * Ifges(6,4) / 0.2e1 + t420 + (-t160 / 0.2e1 - t393) * Ifges(6,6) + (t400 - t401) * Ifges(6,2);
t16 = Ifges(7,4) * t50 + Ifges(7,2) * t51 + Ifges(7,6) * t103;
t413 = t16 / 0.2e1;
t412 = Ifges(7,1) * t407 + Ifges(7,4) * t406 + Ifges(7,5) * t400;
t410 = t43 / 0.2e1;
t405 = -t78 / 0.2e1;
t143 = Ifges(6,4) * t145;
t359 = t163 * Ifges(6,5);
t79 = t146 * Ifges(6,1) + t143 + t359;
t403 = -t79 / 0.2e1;
t399 = -t106 / 0.2e1;
t397 = -t107 / 0.2e1;
t395 = -t144 / 0.2e1;
t391 = t208 / 0.2e1;
t387 = mrSges(5,3) * t65;
t6 = -pkin(5) * t160 - t8;
t381 = t267 * t6;
t378 = t31 * mrSges(6,3);
t377 = t32 * mrSges(6,3);
t376 = t66 * mrSges(5,3);
t372 = mrSges(5,3) * t160;
t371 = mrSges(5,3) * t161;
t369 = Ifges(6,4) * t267;
t368 = Ifges(6,4) * t270;
t362 = t140 * mrSges(4,3);
t361 = t141 * mrSges(4,3);
t357 = t208 * Ifges(4,4);
t356 = t235 * t54;
t352 = t310 * t270;
t344 = t266 * t267;
t342 = t267 * t310;
t341 = t267 * t269;
t326 = Ifges(6,5) * t102 - Ifges(6,6) * t103 + Ifges(6,3) * t160;
t317 = t257 * t329;
t59 = t259 * t96 - t353 * t97;
t116 = t160 * mrSges(5,1) + t161 * mrSges(5,2);
t301 = Ifges(6,1) * t270 - t369;
t299 = Ifges(7,1) * t266 + t365;
t298 = -Ifges(6,2) * t267 + t368;
t296 = Ifges(7,2) * t269 + t366;
t295 = Ifges(6,5) * t270 - Ifges(6,6) * t267;
t293 = Ifges(7,5) * t266 + Ifges(7,6) * t269;
t29 = mrSges(7,1) * t103 - mrSges(7,3) * t50;
t30 = -mrSges(7,2) * t103 + mrSges(7,3) * t51;
t292 = -t266 * t29 + t269 * t30;
t34 = pkin(11) * t174 + t373;
t154 = t175 * t267 - t270 * t239;
t155 = t175 * t270 + t239 * t267;
t75 = t128 * t353 - t259 * t137;
t71 = -t239 * pkin(4) - t75;
t42 = t154 * pkin(5) - t155 * pkin(11) + t71;
t19 = t266 * t42 + t269 * t34;
t18 = -t266 * t34 + t269 * t42;
t73 = -mrSges(7,2) * t144 + mrSges(7,3) * t106;
t74 = mrSges(7,1) * t144 - mrSges(7,3) * t107;
t291 = -t266 * t73 - t269 * t74;
t38 = -t267 * t72 + t270 * t86;
t40 = t104 * t270 - t267 * t69;
t118 = t155 * t269 + t174 * t266;
t117 = -t155 * t266 + t174 * t269;
t285 = -(-qJ(2) * t323 + t252) * t260 + t238 * t263;
t13 = t267 * t94 + t270 * t60 + t86 * t329 - t330 * t72;
t240 = (mrSges(3,1) * t265 - mrSges(3,3) * t350) * qJD(1);
t241 = (-mrSges(3,2) * t265 + mrSges(3,3) * t347) * qJD(1);
t204 = Ifges(4,4) * t207;
t180 = mrSges(4,1) * t233 - mrSges(4,3) * t208;
t179 = -mrSges(4,2) * t233 + mrSges(4,3) * t207;
t172 = -mrSges(4,1) * t207 + mrSges(4,2) * t208;
t162 = mrSges(4,1) * t195 + mrSges(4,2) * t194;
t152 = t208 * Ifges(4,1) + t233 * Ifges(4,5) + t204;
t151 = t207 * Ifges(4,2) + t233 * Ifges(4,6) + t357;
t149 = -mrSges(5,2) * t233 + mrSges(5,3) * t310;
t124 = -mrSges(5,1) * t310 + mrSges(5,2) * t280;
t115 = qJD(5) * t155 + t169 * t267;
t114 = -qJD(5) * t154 + t169 * t270;
t111 = Ifges(5,1) * t280 + t233 * Ifges(5,5) + t432;
t110 = Ifges(5,2) * t310 + t233 * Ifges(5,6) + t433;
t108 = -mrSges(6,2) * t163 + mrSges(6,3) * t145;
t92 = pkin(5) * t146 - pkin(11) * t145;
t81 = -mrSges(6,2) * t160 - mrSges(6,3) * t103;
t77 = t146 * Ifges(6,5) + t145 * Ifges(6,6) + t163 * Ifges(6,3);
t61 = mrSges(6,1) * t103 + mrSges(6,2) * t102;
t57 = -qJD(6) * t118 - t114 * t266 + t166 * t269;
t56 = qJD(6) * t117 + t114 * t269 + t166 * t266;
t35 = -pkin(5) * t280 - t40;
t33 = -pkin(5) * t174 - t38;
t26 = pkin(5) * t115 - pkin(11) * t114 + t59;
t24 = t266 * t92 + t269 * t31;
t23 = -t266 * t31 + t269 * t92;
t12 = -pkin(5) * t166 - t14;
t11 = pkin(11) * t166 + t13;
t4 = -qJD(6) * t19 - t11 * t266 + t26 * t269;
t3 = qJD(6) * t18 + t11 * t269 + t26 * t266;
t17 = [(Ifges(7,5) * t118 + Ifges(7,6) * t117) * t400 + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t406 + m(6) * (t13 * t32 + t14 * t31 + t373 * t7 + t38 * t8 + t54 * t71) + t373 * t81 + t310 * (Ifges(5,4) * t169 - Ifges(5,2) * t166) / 0.2e1 + (Ifges(4,5) * t209 + Ifges(5,5) * t169 - Ifges(4,6) * t210 - Ifges(5,6) * t166) * t233 / 0.2e1 + (Ifges(4,1) * t209 - Ifges(4,4) * t210) * t391 + t173 * (mrSges(4,1) * t210 + mrSges(4,2) * t209) + t207 * (Ifges(4,4) * t209 - Ifges(4,2) * t210) / 0.2e1 + (Ifges(7,5) * t56 + Ifges(7,6) * t57 + Ifges(7,3) * t115) * t394 + (Ifges(7,1) * t56 + Ifges(7,4) * t57 + Ifges(7,5) * t115) * t396 + (Ifges(7,4) * t56 + Ifges(7,2) * t57 + Ifges(7,6) * t115) * t398 + t115 * t405 + (Ifges(7,1) * t118 + Ifges(7,4) * t117) * t407 + m(7) * (t1 * t19 + t10 * t3 + t12 * t27 + t18 * t2 + t33 * t6 + t4 * t9) + (mrSges(5,2) * t187 + mrSges(5,3) * t54 + Ifges(5,1) * t161 - Ifges(5,4) * t160) * t175 + t115 * t439 - t210 * t361 + (Ifges(5,1) * t169 - Ifges(5,4) * t166) * t436 + t424 * t438 + (-mrSges(4,1) * t215 + mrSges(4,2) * t216) * t244 + (t1 * t117 + t10 * t57 - t118 * t2 - t56 * t9) * mrSges(7,3) + m(3) * ((t263 * t334 + (qJ(2) * t350 - t255) * t260) * qJD(1) + t285) * t332 + m(5) * (t142 * t188 + t153 * t187 - t54 * t75 + t55 * t76 + t60 * t66) + (t129 * t215 - t130 * t216) * mrSges(4,3) - t169 * t387 - t166 * t376 - t75 * t371 - t76 * t372 - t209 * t362 + t172 * t247 + t18 * t29 + (-mrSges(6,3) * t8 + t416 + t447) * t155 + t19 * t30 + t33 * t22 + t56 * t45 / 0.2e1 + m(4) * (t129 * t148 + t130 * t147 + t138 * t141 + t139 * t140 + (qJD(1) * t176 + t173) * t247) + t57 * t44 / 0.2e1 + t27 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t12 * t62 + t71 * t61 + t3 * t73 + t4 * t74 + t38 * t80 + (-t7 * mrSges(6,3) - Ifges(6,4) * t402 + t414) * t154 + t13 * t108 + t14 * t109 + t114 * t79 / 0.2e1 + t63 * (mrSges(6,1) * t115 + mrSges(6,2) * t114) + t6 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) + (Ifges(6,3) * t393 + Ifges(6,6) * t401 + Ifges(6,5) * t402 - Ifges(5,4) * t161 + Ifges(5,2) * t160 + t187 * mrSges(5,1) - t55 * mrSges(5,3) + t326 / 0.2e1 + t421) * t174 + t60 * t149 + t153 * t116 + t166 * t77 / 0.2e1 + t31 * (mrSges(6,1) * t166 - mrSges(6,3) * t114) + t32 * (-mrSges(6,2) * t166 - mrSges(6,3) * t115) + t145 * (Ifges(6,4) * t114 - Ifges(6,2) * t115 + Ifges(6,6) * t166) / 0.2e1 + t146 * (Ifges(6,1) * t114 - Ifges(6,4) * t115 + Ifges(6,5) * t166) / 0.2e1 + t163 * (Ifges(6,5) * t114 - Ifges(6,6) * t115 + Ifges(6,3) * t166) / 0.2e1 - t166 * t110 / 0.2e1 + t169 * t111 / 0.2e1 + t142 * (mrSges(5,1) * t166 + mrSges(5,2) * t169) + t176 * t162 + (-m(5) * t65 + m(6) * t63 - t429) * t59 + t138 * t179 + t139 * t180 + t188 * t124 + 0.2e1 * t241 * t321 - 0.2e1 * t240 * t322 + (t431 / 0.2e1 - t430 / 0.2e1 + t415) * t239 + t209 * t152 / 0.2e1 - t210 * t151 / 0.2e1 - t115 * t434 - (t148 * mrSges(4,3) + Ifges(4,4) * t216 + Ifges(4,2) * t215 + Ifges(4,6) * t438) * t195 + (-t147 * mrSges(4,3) + Ifges(4,1) * t216 + Ifges(4,4) * t215 + Ifges(4,5) * t438) * t194 + t115 * t410 + t118 * t412 + t117 * t413; -t450 * t149 + (-t160 * t236 + t161 * t235) * mrSges(5,3) + t446 * t73 + (t116 + t162) * t264 + t445 * t74 - t451 * t108 + t374 * t221 + t185 * t29 + t186 * t30 + t222 * t81 - t228 * t180 - t229 * t179 + t235 * t61 - m(4) * (t140 * t228 + t141 * t229) + t354 * t444 + t429 * t442 + (-m(3) * t285 + t260 * t240 - t263 * t241) * t333 + (t1 * t186 + t10 * t446 + t185 * t2 + t221 * t6 - t27 * t444 + t445 * t9) * m(7) + (-t221 * t8 + t222 * t7 + t31 * t444 - t32 * t451 - t442 * t63 + t356) * m(6) + (t187 * t264 + t236 * t55 + t442 * t65 - t450 * t66 + t356) * m(5) + (m(4) * (t130 * t271 + t141 * t331 + t284 + (-qJD(3) * t140 + t129) * t268) + (t179 * t271 - t180 * t268) * qJD(3) + (-t194 * t271 - t195 * t268) * mrSges(4,3) + (-m(4) * t173 - m(5) * t142 - t124 - t172) * t323) * t261; -t31 * mrSges(6,1) * t280 + t32 * mrSges(6,2) * t280 + t369 * t401 + (t317 - t35) * t62 - (-Ifges(5,2) * t280 + t111 + t432) * t310 / 0.2e1 - (Ifges(5,1) * t310 - t433 + t77) * t280 / 0.2e1 + t452 * t108 + (-t293 * t394 - t296 * t398 - t299 * t396) * t328 + (-t40 - t317) * t109 - (t266 * t45 + t269 * t44) * t328 / 0.2e1 + (t355 + t79) * t329 / 0.2e1 + (t145 * t298 + t146 * t301 + t163 * t295) * qJD(5) / 0.2e1 - (-Ifges(4,2) * t208 + t152 + t204) * t207 / 0.2e1 + t280 * t376 + t151 * t391 + t352 * t403 + t440 * t409 + t424 + t368 * t402 - (Ifges(4,5) * t207 + Ifges(5,5) * t310 - Ifges(4,6) * t208 - Ifges(5,6) * t280) * t233 / 0.2e1 - t146 * (Ifges(6,5) * t280 + t301 * t310) / 0.2e1 - t163 * (Ifges(6,3) * t280 + t295 * t310) / 0.2e1 - t145 * (Ifges(6,6) * t280 + t298 * t310) / 0.2e1 - t142 * (mrSges(5,1) * t280 + mrSges(5,2) * t310) + t310 * t387 + t302 * t381 + t208 * t361 + t207 * t362 + t110 * t436 + t415 + (Ifges(7,1) * t120 + Ifges(7,4) * t119) * t397 + (Ifges(7,5) * t120 + Ifges(7,6) * t119) * t395 + ((t259 * t55 - t353 * t54) * pkin(3) - t142 * t385 + t65 * t68 - t66 * t69) * m(5) - t124 * t385 - t372 * t384 - t329 * t378 - t324 * t371 - t208 * (Ifges(4,1) * t207 - t357) / 0.2e1 - t16 * t344 / 0.2e1 + (Ifges(7,4) * t120 + Ifges(7,2) * t119) * t399 + t163 * t63 * (mrSges(6,1) * t267 + mrSges(6,2) * t270) + t374 * t257 * t267 + (t294 * t400 + t297 * t406 + t300 * t407 + (Ifges(7,5) * t396 + Ifges(7,6) * t398 + Ifges(7,3) * t394) * qJD(5) + t416) * t267 + (qJD(5) * t417 + t257 * t81 - t414) * t270 + (Ifges(7,5) * t397 + Ifges(7,6) * t399 + Ifges(7,3) * t395 + t418) * t342 - t120 * t45 / 0.2e1 + (-t377 + t405 + t410 + t419) * t330 - t69 * t149 + (-t31 * t40 - t32 * t41 - t63 * t68 + t258 * t54 + ((-t267 * t32 - t270 * t31) * qJD(5) + t423) * t257) * m(6) + (t31 * t352 + t32 * t342 + t423) * mrSges(6,3) + (t426 * mrSges(7,1) - t425 * mrSges(7,2)) * t27 + (-t1 * t344 - t10 * t426 - t2 * t341 + t425 * t9) * mrSges(7,3) + t427 * t73 + (-t27 * t35 + t1 * t220 + t2 * t219 + (t27 * t329 + t381) * t257 + t428 * t9 + t427 * t10) * m(7) + t428 * t74 + t429 * t68 - t140 * t179 + t141 * t180 - t173 * (mrSges(4,1) * t208 + mrSges(4,2) * t207) + t219 * t29 + t220 * t30 + t258 * t61 + t341 * t412; -t119 * t74 - t120 * t73 - t310 * t149 + t429 * t280 + (-t310 * t108 + (-t266 * t74 + t269 * t73 + t108) * qJD(5) - t374) * t270 + (qJD(6) * t291 - t163 * t354 + t292 + t81) * t267 + t116 + (-t10 * t120 - t119 * t9 - t27 * t342 + (-t6 + (t10 * t269 - t9 * t266) * qJD(5)) * t270 + (qJD(5) * t27 - qJD(6) * t303 + t304) * t267) * m(7) + (-t280 * t63 + t267 * t7 + t270 * t8 + t163 * (-t267 * t31 + t270 * t32)) * m(6) + (t280 * t65 - t310 * t66 + t187) * m(5); t354 * t32 + (t378 + t403 - t63 * mrSges(6,2) - t143 / 0.2e1 - t359 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t146 - t422) * t145 + t422 * qJD(6) + t304 * mrSges(7,3) - pkin(5) * t22 + t326 + (-t364 / 0.2e1 - t363 / 0.2e1 - t360 / 0.2e1 + t358 / 0.2e1 + t377 - t63 * mrSges(6,1) + t370 / 0.2e1 + t418) * t146 + t293 * t400 + t296 * t406 + t299 * t407 + t269 * t413 + t6 * (-mrSges(7,1) * t269 + mrSges(7,2) * t266) - t24 * t73 - t23 * t74 - t31 * t108 + t266 * t412 + (-pkin(5) * t6 - t10 * t24 - t23 * t9 - t27 * t32) * m(7) + (t292 + m(7) * t304 + (-m(7) * t303 + t291) * qJD(6)) * pkin(11) + t421; (Ifges(7,1) * t106 - t367) * t397 + t44 * t396 + (Ifges(7,5) * t106 - Ifges(7,6) * t107) * t395 - t9 * t73 + t10 * t74 - t27 * (mrSges(7,1) * t107 + mrSges(7,2) * t106) + (t10 * t107 + t106 * t9) * mrSges(7,3) + t15 + (-Ifges(7,2) * t107 + t105 + t45) * t399 + t420;];
tauc  = t17(:);
