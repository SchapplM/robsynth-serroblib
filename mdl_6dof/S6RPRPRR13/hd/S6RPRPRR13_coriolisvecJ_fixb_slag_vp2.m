% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:10:37
% EndTime: 2018-11-23 16:10:54
% DurationCPUTime: 17.25s
% Computational Cost: add. (17768->711), mult. (58945->976), div. (0->0), fcn. (49596->12), ass. (0->322)
t241 = sin(pkin(12));
t345 = cos(pkin(6));
t319 = pkin(1) * t345;
t239 = t241 * t319;
t243 = sin(pkin(6));
t343 = cos(pkin(12));
t312 = t243 * t343;
t300 = qJD(1) * t312;
t210 = qJ(2) * t300 + qJD(1) * t239;
t344 = cos(pkin(7));
t293 = t344 * t343;
t242 = sin(pkin(7));
t310 = t345 * t242;
t444 = t243 * t293 + t310;
t254 = t444 * pkin(9);
t168 = qJD(1) * t254 + t210;
t240 = t343 * t319;
t236 = qJD(1) * t240;
t339 = t241 * t243;
t256 = t345 * pkin(2) + (-pkin(9) * t344 - qJ(2)) * t339;
t177 = qJD(1) * t256 + t236;
t202 = (-pkin(9) * t241 * t242 - pkin(2) * t343 - pkin(1)) * t243;
t195 = qJD(1) * t202 + qJD(2);
t246 = sin(qJ(3));
t372 = cos(qJ(3));
t301 = t344 * t372;
t318 = t242 * t372;
t102 = t168 * t246 - t177 * t301 - t195 * t318;
t432 = -qJD(4) - t102;
t189 = t246 * t310 + (t241 * t372 + t246 * t293) * t243;
t180 = t189 * qJD(1);
t245 = sin(qJ(5));
t248 = cos(qJ(5));
t299 = pkin(5) * t248 + pkin(11) * t245;
t446 = qJD(5) * t299 - (-pkin(4) - t299) * t180 - t432;
t330 = qJD(1) * t243;
t317 = t241 * t330;
t414 = t444 * t372;
t415 = t414 * qJD(1);
t179 = t246 * t317 - t415;
t327 = qJD(5) * t248;
t395 = pkin(3) + pkin(10);
t342 = qJ(4) * t179;
t106 = t180 * t395 + t342;
t158 = t372 * t168;
t309 = t344 * t177;
t338 = t242 * t246;
t103 = t195 * t338 + t246 * t309 + t158;
t83 = -pkin(4) * t179 + t103;
t44 = t248 * t106 + t245 * t83;
t445 = -pkin(11) * t179 + t395 * t327 + t44;
t255 = t243 * (t241 * t301 + t246 * t343);
t200 = qJD(1) * t255;
t440 = qJD(3) * t338 - t200;
t426 = mrSges(5,2) - mrSges(4,1);
t435 = Ifges(5,4) - Ifges(4,5);
t434 = Ifges(5,5) - Ifges(4,6);
t230 = pkin(5) * t245 - pkin(11) * t248 + qJ(4);
t244 = sin(qJ(6));
t247 = cos(qJ(6));
t336 = t245 * t395;
t205 = t247 * t230 + t244 * t336;
t443 = qJD(6) * t205 + t446 * t244 - t445 * t247;
t206 = t244 * t230 - t247 * t336;
t442 = -qJD(6) * t206 + t445 * t244 + t446 * t247;
t260 = -t245 * t344 - t248 * t318;
t305 = t242 * t317;
t417 = qJD(5) * t260 + t440 * t245 - t248 * t305;
t276 = t372 * t312;
t302 = t344 * t339;
t201 = (-t246 * t302 + t276) * qJD(1);
t313 = qJD(3) * t372;
t303 = t242 * t313;
t441 = -t201 + t303;
t176 = qJD(5) + t180;
t294 = t345 * t344;
t207 = -qJD(1) * t294 + t242 * t300 - qJD(3);
t262 = pkin(4) * t180 + t102;
t59 = t207 * t395 + qJD(4) + t262;
t138 = -t177 * t242 + t344 * t195;
t270 = -qJ(4) * t180 + t138;
t63 = t179 * t395 + t270;
t28 = -t245 * t63 + t248 * t59;
t26 = -pkin(5) * t176 - t28;
t144 = t179 * t245 - t207 * t248;
t116 = mrSges(6,1) * t176 - mrSges(6,3) * t144;
t113 = -t144 * t244 + t176 * t247;
t114 = t144 * t247 + t176 * t244;
t56 = -mrSges(7,1) * t113 + mrSges(7,2) * t114;
t346 = -t56 + t116;
t439 = m(7) * t26 - t346;
t29 = t245 * t59 + t248 * t63;
t27 = pkin(11) * t176 + t29;
t143 = t179 * t248 + t207 * t245;
t204 = t207 * qJ(4);
t64 = -t204 + t83;
t39 = -pkin(5) * t143 - pkin(11) * t144 + t64;
t10 = t244 * t39 + t247 * t27;
t429 = t10 * mrSges(7,2);
t9 = -t244 * t27 + t247 * t39;
t431 = t9 * mrSges(7,1);
t438 = t429 - t431;
t436 = mrSges(5,1) + mrSges(4,3);
t175 = Ifges(4,4) * t179;
t422 = t176 * Ifges(6,3);
t423 = t144 * Ifges(6,5);
t424 = t143 * Ifges(6,6);
t433 = t180 * Ifges(4,1) - t207 * Ifges(4,5) - t175 + t422 + t423 + t424;
t229 = qJD(2) * t276;
t275 = qJD(2) * t302;
t263 = qJD(1) * t275;
t274 = qJD(3) * t301;
t88 = -t246 * (qJD(3) * t168 + t263) + qJD(1) * t229 + t177 * t274 + t195 * t303;
t70 = t207 * qJD(4) - t88;
t212 = -t242 * t312 + t294;
t331 = qJ(2) * t312 + t239;
t186 = t254 + t331;
t190 = t240 + t256;
t95 = -t246 * (qJD(3) * t186 + t275) + t190 * t274 + t202 * t303 + t229;
t90 = -t212 * qJD(4) - t95;
t350 = t180 * mrSges(4,3);
t418 = -mrSges(5,1) * t180 + t207 * t426 - t350;
t182 = t189 * qJD(3);
t167 = qJD(1) * t182;
t111 = qJD(5) * t143 + t167 * t245;
t324 = t246 * t339;
t304 = qJD(3) * t324;
t166 = qJD(1) * t304 - qJD(3) * t415;
t48 = qJD(6) * t113 + t111 * t247 - t166 * t244;
t401 = t48 / 0.2e1;
t49 = -qJD(6) * t114 - t111 * t244 - t166 * t247;
t400 = t49 / 0.2e1;
t394 = t111 / 0.2e1;
t112 = qJD(5) * t144 - t248 * t167;
t393 = -t112 / 0.2e1;
t392 = t112 / 0.2e1;
t384 = -t166 / 0.2e1;
t430 = t212 / 0.2e1;
t428 = t28 * mrSges(6,1);
t427 = t29 * mrSges(6,2);
t22 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t86 = -mrSges(6,1) * t166 - mrSges(6,3) * t111;
t425 = t22 - t86;
t295 = t10 * t244 + t247 * t9;
t271 = t295 * mrSges(7,3);
t216 = -t245 * t318 + t248 * t344;
t197 = t216 * t247 + t244 * t338;
t421 = -qJD(6) * t197 - t244 * t417 + t247 * t441;
t196 = -t216 * t244 + t247 * t338;
t420 = qJD(6) * t196 + t244 * t441 + t247 * t417;
t351 = t179 * mrSges(4,3);
t146 = mrSges(4,2) * t207 - t351;
t147 = mrSges(5,1) * t179 + mrSges(5,3) * t207;
t419 = t146 - t147;
t282 = t245 * t28 - t248 * t29;
t328 = qJD(5) * t245;
t252 = qJD(2) * t255;
t89 = qJD(1) * t252 + (t158 + (t195 * t242 + t309) * t246) * qJD(3);
t61 = -t166 * pkin(4) + t89;
t329 = qJD(2) * t243;
t316 = t241 * t329;
t226 = t242 * t316;
t222 = qJD(1) * t226;
t268 = qJ(4) * t166 - qJD(4) * t180 + t222;
t71 = t167 * t395 + t268;
t7 = t245 * t61 + t248 * t71 + t59 * t327 - t328 * t63;
t8 = -qJD(5) * t29 - t245 * t71 + t248 * t61;
t413 = t282 * t180 - t245 * t7 - t248 * t8;
t410 = t166 * t435 + t434 * t167;
t57 = -t167 * pkin(4) - t70;
t25 = t112 * pkin(5) - t111 * pkin(11) + t57;
t5 = -pkin(11) * t166 + t7;
t1 = qJD(6) * t9 + t244 * t25 + t247 * t5;
t2 = -qJD(6) * t10 - t244 * t5 + t247 * t25;
t298 = t1 * t247 - t2 * t244;
t97 = -mrSges(6,1) * t143 + mrSges(6,2) * t144;
t409 = t97 + t419;
t291 = mrSges(7,1) * t244 + mrSges(7,2) * t247;
t373 = t247 / 0.2e1;
t374 = -t244 / 0.2e1;
t142 = qJD(6) - t143;
t362 = Ifges(7,4) * t114;
t41 = Ifges(7,2) * t113 + Ifges(7,6) * t142 + t362;
t110 = Ifges(7,4) * t113;
t42 = Ifges(7,1) * t114 + Ifges(7,5) * t142 + t110;
t408 = -t26 * t291 - t373 * t42 - t374 * t41;
t253 = t246 * t186 - t190 * t301 - t202 * t318;
t72 = t189 * pkin(4) - t212 * t395 + t253;
t188 = t324 - t414;
t145 = -t190 * t242 + t344 * t202;
t269 = -qJ(4) * t189 + t145;
t80 = t188 * t395 + t269;
t366 = t245 * t72 + t248 * t80;
t181 = -qJD(3) * t414 + t304;
t173 = t372 * t186;
t308 = t344 * t190;
t96 = t252 + (t173 + (t202 * t242 + t308) * t246) * qJD(3);
t77 = -t181 * pkin(4) + t96;
t267 = qJ(4) * t181 - qJD(4) * t189 + t226;
t84 = t182 * t395 + t267;
t14 = -qJD(5) * t366 - t245 * t84 + t248 * t77;
t407 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t141 = Ifges(6,4) * t143;
t353 = t176 * Ifges(6,5);
t355 = t144 * Ifges(6,1);
t76 = t141 + t353 + t355;
t396 = -t76 / 0.2e1;
t406 = t271 - t353 / 0.2e1 + t28 * mrSges(6,3) - t141 / 0.2e1 + t396 - t64 * mrSges(6,2);
t17 = Ifges(7,5) * t48 + Ifges(7,6) * t49 + Ifges(7,3) * t112;
t405 = t17 / 0.2e1;
t404 = Ifges(7,1) * t401 + Ifges(7,4) * t400 + Ifges(7,5) * t392;
t357 = t142 * Ifges(7,3);
t358 = t114 * Ifges(7,5);
t359 = t113 * Ifges(7,6);
t40 = t357 + t358 + t359;
t403 = t40 / 0.2e1;
t402 = -t41 / 0.2e1;
t399 = -t111 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t392 + t166 * Ifges(6,6) / 0.2e1;
t398 = Ifges(6,1) * t394 + Ifges(6,4) * t393 + Ifges(6,5) * t384;
t352 = t176 * Ifges(6,6);
t356 = t143 * Ifges(6,2);
t365 = Ifges(6,4) * t144;
t75 = t352 + t356 + t365;
t397 = -t75 / 0.2e1;
t391 = -t113 / 0.2e1;
t390 = t113 / 0.2e1;
t389 = -t114 / 0.2e1;
t388 = t114 / 0.2e1;
t386 = -t142 / 0.2e1;
t385 = t142 / 0.2e1;
t383 = -t179 / 0.2e1;
t382 = t179 / 0.2e1;
t381 = -t180 / 0.2e1;
t380 = t180 / 0.2e1;
t378 = t189 / 0.2e1;
t377 = -t207 / 0.2e1;
t364 = Ifges(6,4) * t245;
t363 = Ifges(6,4) * t248;
t361 = Ifges(7,4) * t244;
t360 = Ifges(7,4) * t247;
t354 = t166 * mrSges(5,1);
t349 = t180 * Ifges(4,4);
t348 = t180 * Ifges(5,6);
t347 = -t147 + t97;
t341 = t180 * t245;
t337 = t244 * t248;
t335 = t247 * t248;
t334 = t248 * t395;
t326 = t372 * t89;
t322 = Ifges(6,5) * t111 - Ifges(6,6) * t112 - Ifges(6,3) * t166;
t118 = t202 * t338 + t246 * t308 + t173;
t100 = -t212 * qJ(4) - t118;
t296 = t10 * t247 - t244 * t9;
t292 = mrSges(7,1) * t247 - mrSges(7,2) * t244;
t290 = Ifges(7,1) * t247 - t361;
t289 = Ifges(7,1) * t244 + t360;
t288 = -Ifges(7,2) * t244 + t360;
t287 = Ifges(7,2) * t247 + t361;
t286 = Ifges(7,5) * t247 - Ifges(7,6) * t244;
t285 = Ifges(7,5) * t244 + Ifges(7,6) * t247;
t31 = mrSges(7,1) * t112 - mrSges(7,3) * t48;
t32 = -mrSges(7,2) * t112 + mrSges(7,3) * t49;
t284 = -t244 * t31 + t247 * t32;
t34 = pkin(11) * t189 + t366;
t151 = t188 * t245 + t212 * t248;
t278 = t248 * t188 - t212 * t245;
t85 = -pkin(4) * t188 - t100;
t45 = -pkin(5) * t278 - pkin(11) * t151 + t85;
t16 = t244 * t45 + t247 * t34;
t15 = -t244 * t34 + t247 * t45;
t68 = -mrSges(7,2) * t142 + mrSges(7,3) * t113;
t69 = mrSges(7,1) * t142 - mrSges(7,3) * t114;
t283 = -t244 * t68 - t247 * t69;
t35 = -t245 * t80 + t248 * t72;
t43 = -t106 * t245 + t248 * t83;
t126 = t151 * t247 + t189 * t244;
t125 = -t151 * t244 + t189 * t247;
t13 = t245 * t77 + t248 * t84 + t72 * t327 - t328 * t80;
t265 = -(-qJ(2) * t317 + t236) * t241 + t210 * t343;
t213 = (mrSges(3,1) * t345 - mrSges(3,3) * t339) * qJD(1);
t214 = (-mrSges(3,2) * t345 + mrSges(3,3) * t312) * qJD(1);
t62 = -t182 * pkin(4) - t90;
t251 = -t357 / 0.2e1 - t358 / 0.2e1 - t359 / 0.2e1 + t352 / 0.2e1 - t40 / 0.2e1 + t75 / 0.2e1 - t64 * mrSges(6,1) + t29 * mrSges(6,3) + t365 / 0.2e1 + t438;
t250 = t286 * t385 + t288 * t390 + t290 * t388 - t408;
t174 = Ifges(5,6) * t179;
t140 = -t180 * t335 - t207 * t244;
t139 = t180 * t337 - t207 * t247;
t133 = -mrSges(5,2) * t179 - mrSges(5,3) * t180;
t132 = mrSges(4,1) * t179 + mrSges(4,2) * t180;
t131 = pkin(3) * t180 + t342;
t128 = mrSges(4,1) * t167 - mrSges(4,2) * t166;
t127 = -mrSges(5,2) * t167 + mrSges(5,3) * t166;
t124 = qJD(5) * t151 - t248 * t182;
t123 = qJD(5) * t278 + t182 * t245;
t121 = -t179 * Ifges(4,2) - t207 * Ifges(4,6) + t349;
t120 = -t207 * Ifges(5,4) - t180 * Ifges(5,2) + t174;
t119 = -t207 * Ifges(5,5) + t179 * Ifges(5,3) - t348;
t115 = -mrSges(6,2) * t176 + mrSges(6,3) * t143;
t104 = pkin(3) * t182 + t267;
t101 = -t212 * pkin(3) + t253;
t99 = pkin(3) * t188 + t269;
t98 = pkin(5) * t144 - pkin(11) * t143;
t94 = pkin(3) * t167 + t268;
t93 = t204 - t103;
t92 = pkin(3) * t207 - t432;
t91 = pkin(3) * t179 + t270;
t87 = mrSges(6,2) * t166 - mrSges(6,3) * t112;
t55 = mrSges(6,1) * t112 + mrSges(6,2) * t111;
t54 = -qJD(6) * t126 - t123 * t244 - t181 * t247;
t53 = qJD(6) * t125 + t123 * t247 - t181 * t244;
t37 = pkin(5) * t179 - t43;
t33 = -pkin(5) * t189 - t35;
t30 = t124 * pkin(5) - t123 * pkin(11) + t62;
t24 = t244 * t98 + t247 * t28;
t23 = -t244 * t28 + t247 * t98;
t18 = t48 * Ifges(7,4) + t49 * Ifges(7,2) + t112 * Ifges(7,6);
t12 = pkin(5) * t181 - t14;
t11 = -pkin(11) * t181 + t13;
t6 = pkin(5) * t166 - t8;
t4 = -qJD(6) * t16 - t11 * t244 + t247 * t30;
t3 = qJD(6) * t15 + t11 * t247 + t244 * t30;
t19 = [(t124 * t64 + t189 * t8) * mrSges(6,1) + (t1 * t125 + t10 * t54 - t126 * t2 - t53 * t9) * mrSges(7,3) + (Ifges(7,5) * t126 + Ifges(7,6) * t125) * t392 + (Ifges(7,4) * t126 + Ifges(7,2) * t125) * t400 + (Ifges(6,4) * t151 + Ifges(6,6) * t189) * t393 + (-Ifges(4,1) * t378 + Ifges(5,4) * t430 - Ifges(5,2) * t189 - t253 * mrSges(4,3) + (Ifges(5,6) + Ifges(4,4) / 0.2e1) * t188) * t166 - (mrSges(6,1) * t57 - Ifges(6,4) * t394 + Ifges(7,5) * t401 - Ifges(6,2) * t393 - Ifges(6,6) * t384 + Ifges(7,6) * t400 + Ifges(7,3) * t392 + t399 + t405 + t407) * t278 + (-t123 * t28 - t124 * t29 - t151 * t8 + t278 * t7) * mrSges(6,3) + m(7) * (t1 * t16 + t10 * t3 + t12 * t26 + t15 * t2 + t33 * t6 + t4 * t9) + t322 * t378 + (t100 * mrSges(5,1) - t118 * mrSges(4,3) - Ifges(5,6) * t189 + (-Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1) * t212 + (Ifges(4,2) + Ifges(5,3)) * t188 + (-t378 - t189 / 0.2e1) * Ifges(4,4)) * t167 + m(6) * (t13 * t29 + t14 * t28 + t35 * t8 + t366 * t7 + t57 * t85 + t62 * t64) + t366 * t87 + (Ifges(6,1) * t151 + Ifges(6,5) * t189) * t394 + (mrSges(4,1) * t188 + mrSges(4,2) * t189) * t222 + (-Ifges(4,4) * t188 + Ifges(4,5) * t212 + Ifges(6,5) * t151 + (Ifges(4,1) + Ifges(6,3)) * t189) * t384 + (t123 * t64 + t151 * t57 - t189 * t7) * mrSges(6,2) + t151 * t398 + t124 * t403 + t126 * t404 + t132 * t226 + t143 * (Ifges(6,4) * t123 - Ifges(6,2) * t124) / 0.2e1 + m(4) * (t103 * t95 + t118 * t88 + (qJD(1) * t145 + t138) * t226) + m(5) * (t100 * t70 + t104 * t91 + t90 * t93 + t94 * t99) + (Ifges(7,5) * t53 + Ifges(7,6) * t54 + Ifges(7,3) * t124) * t385 + (Ifges(7,1) * t53 + Ifges(7,4) * t54 + Ifges(7,5) * t124) * t388 + (Ifges(7,4) * t53 + Ifges(7,2) * t54 + Ifges(7,6) * t124) * t390 + t124 * t397 + (-t433 / 0.2e1 - Ifges(4,1) * t380 + Ifges(5,2) * t381 + Ifges(5,6) * t382 - Ifges(4,4) * t383 + t91 * mrSges(5,3) - t138 * mrSges(4,2) + t120 / 0.2e1 - t422 / 0.2e1 - t423 / 0.2e1 - t424 / 0.2e1 + t427 - t428 - t92 * mrSges(5,1) - t102 * mrSges(4,3) + t435 * t377) * t181 - t101 * t354 - 0.2e1 * t213 * t316 + t144 * (Ifges(6,1) * t123 - Ifges(6,4) * t124) / 0.2e1 + t70 * (mrSges(5,1) * t188 - mrSges(5,3) * t212) + t88 * (-mrSges(4,2) * t212 - mrSges(4,3) * t188) + t94 * (-mrSges(5,2) * t188 - mrSges(5,3) * t189) + t145 * t128 + t95 * t146 + t90 * t147 + t104 * t133 + t99 * t127 + t125 * t18 / 0.2e1 + t6 * (-mrSges(7,1) * t125 + mrSges(7,2) * t126) + t123 * t76 / 0.2e1 + t13 * t115 + t14 * t116 + t62 * t97 + t85 * t55 + t35 * t86 + t3 * t68 + t4 * t69 + t12 * t56 + t410 * t430 + t124 * t431 + t176 * (Ifges(6,5) * t123 - Ifges(6,6) * t124) / 0.2e1 + (m(4) * t102 + m(5) * t92 - t418) * t96 + (Ifges(7,1) * t126 + Ifges(7,4) * t125) * t401 - t124 * t429 + 0.2e1 * qJD(2) * t214 * t312 + (-Ifges(4,4) * t380 + Ifges(5,6) * t381 + Ifges(5,3) * t382 - Ifges(4,2) * t383 - t91 * mrSges(5,2) + t119 / 0.2e1 + t138 * mrSges(4,1) - t121 / 0.2e1 + t93 * mrSges(5,1) - t103 * mrSges(4,3) + t434 * t377) * t182 + (m(4) * t253 + m(5) * t101 + t189 * t436 + t426 * t212) * t89 + m(3) * ((t343 * t331 + (qJ(2) * t339 - t240) * t241) * qJD(1) + t265) * t329 + t15 * t31 + t16 * t32 + t33 * t22 + t53 * t42 / 0.2e1 + t54 * t41 / 0.2e1 + t26 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53); (-m(4) * t103 + m(5) * t93 - m(6) * t64 - t409) * t201 + t409 * t303 + (t127 + t128) * t344 + t213 * t317 + t55 * t338 + (-t133 - t132) * t305 + (t94 * t344 + (-t326 - t246 * t70 + (t246 * t92 - t372 * t93) * qJD(3)) * t242 - t200 * t92 - t91 * t305) * m(5) + (-t102 * t200 - t138 * t305 + (t263 - t326 + t246 * t88 + (t102 * t246 + t103 * t372) * qJD(3)) * t242) * m(4) + t216 * t87 + t196 * t31 + t197 * t32 - t214 * t300 + t417 * t115 + (t8 * t260 + t7 * t216 + (t246 * t57 + t313 * t64) * t242 + t417 * t29) * m(6) + t420 * t68 + (t1 * t197 + t10 * t420 + t196 * t2 - t260 * t6 + t421 * t9) * m(7) + t421 * t69 - t425 * t260 - m(3) * t265 * t330 + (-m(6) * t28 + t439) * (qJD(5) * t216 - t245 * t305 - t248 * t440) + t436 * (t166 * t318 - t167 * t338) - t418 * t440; (-Ifges(4,1) * t179 + t119 - t349) * t381 + ((-t356 / 0.2e1 - t251) * t248 - (-m(6) * t282 + t248 * t115) * t395 + (t286 * t386 + t290 * t389 + t288 * t391 - t355 / 0.2e1 - t439 * t395 + t406 + t408) * t245) * qJD(5) + (-Ifges(7,5) * t389 - Ifges(7,6) * t391 - Ifges(7,3) * t386 + t397 + t403 - t438) * t180 * t248 - (t245 * t87 - t248 * t425) * t395 + (t26 * mrSges(7,1) - t10 * mrSges(7,3) + Ifges(7,4) * t389 + Ifges(7,2) * t391 + Ifges(7,6) * t386 + t402) * (-t179 * t247 - t244 * t341) + (Ifges(7,5) * t386 + Ifges(7,1) * t389 + Ifges(7,4) * t391 + t9 * mrSges(7,3) - t26 * mrSges(7,2) - t42 / 0.2e1) * (-t179 * t244 + t247 * t341) - m(6) * (-t262 * t64 + t28 * t43 + t29 * t44) + t262 * t97 + (Ifges(5,2) * t179 + t121 + t348) * t380 + t410 + t245 * t399 + (Ifges(7,6) * t245 + t248 * t288) * t400 + (Ifges(7,5) * t245 + t248 * t290) * t401 + t335 * t404 + t245 * t405 - t176 * (-Ifges(6,3) * t179 + (Ifges(6,5) * t245 + Ifges(6,6) * t248) * t180) / 0.2e1 + t442 * t69 + (t1 * t206 + t10 * t443 + t2 * t205 - t26 * t37 + t334 * t6 + t442 * t9) * m(7) + t443 * t68 + (Ifges(6,5) * t248 - Ifges(6,6) * t245) * t384 + (Ifges(7,3) * t245 + t248 * t286) * t392 + (-Ifges(6,2) * t245 + t363) * t393 + (Ifges(6,1) * t248 - t364) * t394 + t341 * t396 + t248 * t398 + (pkin(3) * t166 - qJ(4) * t167 + t179 * t92 - t180 * t93) * mrSges(5,1) + t57 * (mrSges(6,1) * t245 + mrSges(6,2) * t248) + m(6) * (t57 * qJ(4) + t64 * qJD(4) - t334 * t8 - t336 * t7) - t144 * (-Ifges(6,5) * t179 + (Ifges(6,1) * t245 + t363) * t180) / 0.2e1 - t143 * (-Ifges(6,6) * t179 + (Ifges(6,2) * t248 + t364) * t180) / 0.2e1 + t347 * qJD(4) - t18 * t337 / 0.2e1 + t1 * (-mrSges(7,2) * t245 - mrSges(7,3) * t337) + t2 * (mrSges(7,1) * t245 - mrSges(7,3) * t335) + t205 * t31 + t206 * t32 - t91 * (-mrSges(5,2) * t180 + mrSges(5,3) * t179) - t138 * (mrSges(4,1) * t180 - mrSges(4,2) * t179) - t131 * t133 + (-mrSges(7,3) * t296 + t247 * t402 + t26 * t292 + t285 * t386 + t287 * t391 + t289 * t389 + t374 * t42) * t248 * qJD(6) - t44 * t115 - t43 * t116 - t88 * mrSges(4,2) - t70 * mrSges(5,3) + qJ(4) * t55 - t37 * t56 + t179 * t428 - t64 * (-mrSges(6,1) * t248 + mrSges(6,2) * t245) * t180 + t413 * mrSges(6,3) + (t350 + t418) * t103 + (t351 + t419) * t102 + t6 * t291 * t248 + (-pkin(3) * t89 - qJ(4) * t70 - t103 * t92 - t131 * t91 + t432 * t93) * m(5) + t426 * t89 + (Ifges(5,3) * t180 + t120 + t174) * t383 - t179 * t427 + (-Ifges(4,2) * t180 - t175 + t433) * t382 + (t179 * t435 + t180 * t434) * t207 / 0.2e1; -t354 + t180 * t133 - t139 * t69 - t140 * t68 + t347 * t207 + (t180 * t115 + (-t244 * t69 + t247 * t68 + t115) * qJD(5) - t425) * t248 + (qJD(6) * t283 - t176 * t346 + t284 + t87) * t245 + ((qJD(5) * t296 - t6) * t248 + (qJD(5) * t26 - qJD(6) * t295 + t298) * t245 - t10 * t140 - t139 * t9 + t26 * t341) * m(7) + (-qJD(5) * t282 + t207 * t64 - t413) * m(6) + (t180 * t91 - t207 * t93 + t89) * m(5); t322 + t18 * t373 + t251 * t144 + t346 * t29 + t244 * t404 - t28 * t115 + t298 * mrSges(7,3) + (t250 - t271) * qJD(6) - t6 * t292 + t289 * t401 + t285 * t392 + t287 * t400 - t24 * t68 - t23 * t69 + ((Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t144 - t250 + t406) * t143 - t7 * mrSges(6,2) + t8 * mrSges(6,1) - pkin(5) * t22 + (-pkin(5) * t6 - t10 * t24 - t23 * t9 - t26 * t29) * m(7) + (m(7) * t298 + (-m(7) * t295 + t283) * qJD(6) + t284) * pkin(11); -t26 * (mrSges(7,1) * t114 + mrSges(7,2) * t113) + (Ifges(7,1) * t113 - t362) * t389 + t41 * t388 + (Ifges(7,5) * t113 - Ifges(7,6) * t114) * t386 - t9 * t68 + t10 * t69 + (t10 * t114 + t113 * t9) * mrSges(7,3) + t17 + (-Ifges(7,2) * t114 + t110 + t42) * t391 + t407;];
tauc  = t19(:);
